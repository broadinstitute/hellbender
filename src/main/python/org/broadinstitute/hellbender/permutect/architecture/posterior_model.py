from collections import defaultdict
from itertools import chain
from math import ceil

import torch
from matplotlib import pyplot as plt
from torch.utils.tensorboard import SummaryWriter
from tqdm.autonotebook import trange, tqdm

from permutect import utils
from permutect.architecture.artifact_spectra import ArtifactSpectra
from permutect.architecture.overdispersed_binomial_mixture import OverdispersedBinomialMixture
from permutect.architecture.normal_seq_error_spectrum import NormalSeqErrorSpectrum
from permutect.architecture.somatic_spectrum import SomaticSpectrum
from permutect.data.base_datum import DEFAULT_GPU_FLOAT, DEFAULT_CPU_FLOAT
from permutect.data.posterior import PosteriorBatch
from permutect.metrics import plotting
from permutect.utils import Variation, Call
from permutect.metrics.evaluation_metrics import MAX_COUNT, NUM_COUNT_BINS, multiple_of_three_bin_index, multiple_of_three_bin_index_to_count


# TODO: write unit test asserting that this comes out to zero when counts are zero
# given germline, the probability of these particular reads being alt
def germline_log_likelihood(afs, mafs, alt_counts, depths, het_beta=None):
    hom_alpha, hom_beta = torch.tensor([98.0], device=depths.device), torch.tensor([2.0], device=depths.device)
    het_alpha, het_beta_to_use = (None, None) if het_beta is None else (torch.tensor([het_beta], device=depths.device), torch.tensor([het_beta], device=depths.device))
    het_probs = 2 * afs * (1 - afs)
    hom_probs = afs * afs
    het_proportion = het_probs / (het_probs + hom_probs)
    hom_proportion = 1 - het_proportion

    log_mafs = torch.log(mafs)
    log_1m_mafs = torch.log(1 - mafs)
    log_half_het_prop = torch.log(het_proportion / 2)

    ref_counts = depths - alt_counts

    combinatorial_term = torch.lgamma(depths + 1) - torch.lgamma(alt_counts + 1) - torch.lgamma(ref_counts + 1)
    # the following should both be 1D tensors of length batch size
    alt_minor_binomial = combinatorial_term + alt_counts * log_mafs + ref_counts * log_1m_mafs
    alt_major_binomial = combinatorial_term + ref_counts * log_mafs + alt_counts * log_1m_mafs
    alt_minor_ll = log_half_het_prop + (alt_minor_binomial if het_beta is None else utils.beta_binomial(depths, alt_counts, het_alpha, het_beta_to_use))
    alt_major_ll = log_half_het_prop + (alt_major_binomial if het_beta is None else utils.beta_binomial(depths, alt_counts, het_alpha, het_beta_to_use))
    hom_ll = torch.log(hom_proportion) + utils.beta_binomial(depths, alt_counts, hom_alpha, hom_beta)

    return torch.logsumexp(torch.vstack((alt_minor_ll, alt_major_ll, hom_ll)), dim=0)


# TODO: max_mean is hard-coded magic constant!!
def initialize_normal_artifact_spectra():
    return OverdispersedBinomialMixture(num_components=1, max_mean=0.1, mode='beta')


# this works for ArtifactSpectra and OverdispersedBinomialMixture
def plot_artifact_spectra(artifact_spectra, depth: int = None):
    # plot AF spectra in two-column grid with as many rows as needed
    art_spectra_fig, art_spectra_axs = plt.subplots(ceil(len(Variation) / 2), 2, sharex='all', sharey='all')
    for variant_type in Variation:
        n = variant_type
        row, col = int(n / 2), n % 2
        frac, dens = artifact_spectra.spectrum_density_vs_fraction(variant_type, depth)
        art_spectra_axs[row, col].plot(frac.detach().numpy(), dens.detach().numpy(), label=variant_type.name)
        art_spectra_axs[row, col].set_title(variant_type.name + " artifact AF spectrum")
    for ax in art_spectra_fig.get_axes():
        ax.label_outer()
    return art_spectra_fig, art_spectra_axs


class PosteriorModel(torch.nn.Module):
    """

    """
    def __init__(self, variant_log_prior: float, artifact_log_prior: float, num_base_features: int, no_germline_mode: bool = False, device=utils.gpu_if_available(), het_beta: float = None):
        super(PosteriorModel, self).__init__()

        self._device = device
        self._dtype = DEFAULT_GPU_FLOAT if device != torch.device("cpu") else DEFAULT_CPU_FLOAT
        self.no_germline_mode = no_germline_mode
        self.num_base_features = num_base_features
        self.het_beta = het_beta

        # TODO introduce parameters class so that num_components is not hard-coded
        self.somatic_spectrum = SomaticSpectrum(num_components=5)

        # artifact spectra for each variant type.  Variant type encoded as one-hot input vector.
        self.artifact_spectra = ArtifactSpectra(num_components=2)

        # normal sequencing error spectra for each variant type.
        self.normal_seq_error_spectra = torch.nn.ModuleList([NormalSeqErrorSpectrum(num_samples=50, max_mean=0.001) for _ in Variation])

        self.normal_artifact_spectra = initialize_normal_artifact_spectra()

        # pre-softmax priors of different call types [log P(variant), log P(artifact), log P(seq error)] for each variant type
        self._unnormalized_priors_vc = torch.nn.Parameter(torch.ones(len(Variation), len(Call)))
        with torch.no_grad():
            self._unnormalized_priors_vc[:, Call.SOMATIC] = variant_log_prior
            self._unnormalized_priors_vc[:, Call.ARTIFACT] = artifact_log_prior
            self._unnormalized_priors_vc[:, Call.SEQ_ERROR] = 0
            self._unnormalized_priors_vc[:, Call.GERMLINE] = -9999 if self.no_germline_mode else 0
            self._unnormalized_priors_vc[:, Call.NORMAL_ARTIFACT] = artifact_log_prior

        self.to(device=self._device, dtype=self._dtype)

    def make_unnormalized_priors(self, variant_types_b: torch.IntTensor, allele_frequencies_1d: torch.Tensor) -> torch.Tensor:
        result_bc = self._unnormalized_priors_vc[variant_types_b.long(), :].to(device=self._device, dtype=self._dtype)
        result_bc[:, Call.SEQ_ERROR] = 0
        result_bc[:, Call.GERMLINE] = -9999 if self.no_germline_mode else torch.log(1 - torch.square(1 - allele_frequencies_1d))     # 1 minus hom ref probability
        return result_bc   # batch size x len(CallType)

    def posterior_probabilities(self, batch: PosteriorBatch) -> torch.Tensor:
        """
        :param batch:
        :return: non-log probabilities as a 2D tensor, 1st index is batch, 2nd is variant/artifact/seq error
        """
        return torch.nn.functional.softmax(self.log_relative_posteriors(batch), dim=1)

    def error_probabilities(self, batch: PosteriorBatch, germline_mode: bool = False) -> torch.Tensor:
        """
        :param germline_mode: if True, germline classification is not considered an error mode
        :param batch:
        :return: non-log error probabilities as a 1D tensor with length batch size
        """
        assert not (germline_mode and self.no_germline_mode), "germline mode and no-germline mode are incompatible"
        return 1 - self.posterior_probabilities(batch)[:, Call.GERMLINE if germline_mode else Call.SOMATIC]     # 0th column is variant

    def log_posterior_and_ingredients(self, batch: PosteriorBatch) -> torch.Tensor:
        """
        :param batch:
        :batch.seq_error_log_likelihoods() is the probability that these *particular* reads exhibit the alt allele given a
        sequencing error ie an error explainable in terms of base qualities.  For example if we have two alt reads with error
        probability of 0.1 and 0.2, and two ref reads with error probabilities 0.05 and 0.06 this quantity would be
        log(0.1*0.2*0.95*0.94).  This is an annotation emitted by the GATK and by the time it reaches here is a 1D tensor
        of length batch_size.
        :return:
        """
        variant_types = batch.get_variant_types().to(device=self._device, dtype=self._dtype)

        # All log likelihood/relative posterior tensors below have shape batch.size() x len(CallType)
        # spectra tensors contain the likelihood that these *particular* reads (that is, not just the read count) are alt
        # normal log likelihoods contain everything going on in the matched normal sample
        # note that the call to make_unnormalized_priors ensures that no_germline_mode works
        log_priors = torch.nn.functional.log_softmax(self.make_unnormalized_priors(variant_types, batch.get_allele_frequencies()), dim=1)

        # defined as log [ int_0^1 Binom(alt count | depth, f) df ], including the combinatorial N choose N_alt factor
        depths, alt_counts = batch.get_depths(), batch.get_alt_counts()
        normal_depths, normal_alt_counts = batch.get_normal_depths(), batch.get_normal_alt_counts()
        flat_prior_spectra_log_likelihoods = -torch.log(depths + 1)
        somatic_spectrum_log_likelihoods = self.somatic_spectrum.forward(depths, alt_counts)
        tumor_artifact_spectrum_log_likelihood = self.artifact_spectra.forward(batch.get_variant_types(), depths, alt_counts)
        spectra_log_likelihoods = torch.zeros_like(log_priors, device=self._device, dtype=self._dtype)

        # essentially, this corrects the TLOD from M2, computed with a flat prior, to account for the precises somatic spectrum
        spectra_log_likelihoods[:, Call.SOMATIC] = somatic_spectrum_log_likelihoods - flat_prior_spectra_log_likelihoods
        spectra_log_likelihoods[:, Call.ARTIFACT] = tumor_artifact_spectrum_log_likelihood - flat_prior_spectra_log_likelihoods
        spectra_log_likelihoods[:, Call.NORMAL_ARTIFACT] = tumor_artifact_spectrum_log_likelihood - flat_prior_spectra_log_likelihoods # yup, it's the same spectrum
        spectra_log_likelihoods[:, Call.SEQ_ERROR] = -batch.get_tlods_from_m2()
        # spectra_log_likelihoods[:, Call.GERMLINE] is computed below

        normal_log_likelihoods = torch.zeros_like(log_priors)
        normal_seq_error_log_likelihoods = torch.zeros_like(alt_counts).float()

        for var_index, _ in enumerate(Variation):
            mask = (variant_types == var_index)
            log_likelihoods_for_this_type = self.normal_seq_error_spectra[var_index].forward(normal_alt_counts, batch.get_normal_ref_counts())
            normal_seq_error_log_likelihoods += mask * log_likelihoods_for_this_type

        normal_log_likelihoods[:, Call.SOMATIC] = normal_seq_error_log_likelihoods
        normal_log_likelihoods[:, Call.ARTIFACT] = normal_seq_error_log_likelihoods
        normal_log_likelihoods[:, Call.SEQ_ERROR] = normal_seq_error_log_likelihoods

        no_alt_in_normal_mask = normal_alt_counts < 1
        normal_log_likelihoods[:, Call.NORMAL_ARTIFACT] = -9999 * no_alt_in_normal_mask + \
            torch.logical_not(no_alt_in_normal_mask) * self.normal_artifact_spectra.forward(variant_types, normal_depths, normal_alt_counts)

        afs = batch.get_allele_frequencies()
        spectra_log_likelihoods[:, Call.GERMLINE] = germline_log_likelihood(afs, batch.get_mafs(), alt_counts, depths, self.het_beta) - flat_prior_spectra_log_likelihoods

        # it is correct not to subtract the flat prior likelihood from the normal term because this is an absolute likelihood, not
        # relative to seq error as the M2 TLOD is defined
        normal_log_likelihoods[:, Call.GERMLINE] = germline_log_likelihood(afs, batch.get_normal_mafs(), normal_alt_counts, normal_depths, self.het_beta)

        log_posteriors = log_priors + spectra_log_likelihoods + normal_log_likelihoods
        log_posteriors[:, Call.ARTIFACT] += batch.get_artifact_logits()

        log_posteriors[:, Call.NORMAL_ARTIFACT] += batch.get_artifact_logits()

        return log_priors, spectra_log_likelihoods, normal_log_likelihoods, log_posteriors

    def log_relative_posteriors(self, batch: PosteriorBatch) -> torch.Tensor:
        _, _, _, log_posteriors = self.log_posterior_and_ingredients(batch)
        return log_posteriors

    def learn_priors_and_spectra(self, posterior_loader, num_iterations, ignored_to_non_ignored_ratio: float,
                                 summary_writer: SummaryWriter = None, learning_rate: float = 0.001):
        """
        :param summary_writer:
        :param num_iterations:
        :param posterior_loader:
        :param ignored_to_non_ignored_ratio: ratio of sites in which no evidence of variation was found to sites in which
        sufficient evidence was found to emit test data.  Without this parameter (i.e. if it were set to zero) we would
        underestimate the frequency of sequencing error, hence overestimate the prior probability of variation.
        :param artifact_spectra_state_dict: (possibly None) if given, pretrained parameters of self.artifact_spectra
        from train_model.py.  In this case we make sure to freeze this part of the model
        :param artifact_log_priors: (possibly None) 1D tensor with length len(utils.Variation) containing log prior probabilities
        of artifacts for each variation type, from train_model.py.  If given, freeze these parameters.
        :return:
        """
        spectra_and_prior_params = chain(self.somatic_spectrum.parameters(), self.artifact_spectra.parameters(),
                                         [self._unnormalized_priors_vc], self.normal_seq_error_spectra.parameters(),
                                         self.normal_artifact_spectra.parameters())
        optimizer = torch.optim.Adam(spectra_and_prior_params, lr=learning_rate)

        for epoch in trange(1, num_iterations + 1, desc="AF spectra epoch"):
            epoch_loss = utils.StreamingAverage()

            # store posteriors as a list (to be stacked at the end of the epoch) for an M step
            # 'l' for loader, 'b' for batch, 'c' for call type
            posteriors_lbc = []
            alt_counts_lb = []
            depths_lb = []
            types_lb = []

            pbar = tqdm(enumerate(posterior_loader), mininterval=10)
            batch_cpu: PosteriorBatch
            for n, batch_cpu in pbar:
                batch = batch_cpu.copy_to(self._device, self._dtype, non_blocking=self._device.type == 'cuda')
                relative_posteriors = self.log_relative_posteriors(batch)
                log_evidence = torch.logsumexp(relative_posteriors, dim=1)

                posteriors_lbc.append(torch.softmax(relative_posteriors, dim=-1).detach())
                alt_counts_lb.append(batch.get_alt_counts().detach())
                depths_lb.append(batch.get_depths().detach())
                types_lb.append(batch.get_variant_types().detach())

                confidence_mask = torch.abs(batch.get_artifact_logits()) > 3.0
                #loss = -torch.mean(confidence_mask * log_evidence)
                loss = - torch.sum(confidence_mask * log_evidence) / (torch.sum(confidence_mask) + 0.000001)

                # note that we don't multiply by batch size because we take the mean of log evidence above
                # however, we must sum over variant types since each ignored site is simultaneously a missing non-SNV,
                # a missing non-INSERTION etc
                # we use a germline allele frequency of 0.001 for the missing sites but it doesn't really matter
                for var_type_idx, variant_type in enumerate(Variation):
                    log_priors = torch.nn.functional.log_softmax(self.make_unnormalized_priors(torch.LongTensor([var_type_idx]).to(device=self._device, dtype=self._dtype), torch.tensor([0.001], device=self._device)), dim=1)
                    log_seq_error_prior = log_priors.squeeze()[Call.SEQ_ERROR]
                    missing_loss = -ignored_to_non_ignored_ratio * log_seq_error_prior  
                    loss += missing_loss

                utils.backpropagate(optimizer, loss)

                epoch_loss.record_sum(batch.size() * loss.detach().item(), batch.size())
            # iteration over posterior dataloader finished

            # 'n' denotes index of data within entire Posterior Dataset
            posteriors_nc = torch.vstack(posteriors_lbc)
            alt_counts_n = torch.hstack(alt_counts_lb)
            depths_n = torch.hstack(depths_lb)
            types_n = torch.hstack(types_lb)

            self.update_priors_m_step(posteriors_nc, types_n, ignored_to_non_ignored_ratio)
            self.somatic_spectrum.update_m_step(posteriors_nc[:, Call.SOMATIC], alt_counts_n, depths_n)

            if summary_writer is not None:
                summary_writer.add_scalar("spectrum negative log evidence", epoch_loss.get(), epoch)

                for variant_index, variant_type in enumerate(Variation):
                    mean = self.normal_seq_error_spectra[variant_index].get_mean()
                    summary_writer.add_scalar("normal seq error mean fraction for " + variant_type.name, mean, epoch)

                for depth in [10, 20, 30, 50, 100]:
                    art_spectra_fig, art_spectra_axs = plot_artifact_spectra(self.artifact_spectra, depth)
                    summary_writer.add_figure("Artifact AF Spectra at depth = " + str(depth), art_spectra_fig, epoch)

                #normal_seq_error_spectra_fig, normal_seq_error_spectra_axs = plot_artifact_spectra(self.normal_seq_error_spectra)
                #summary_writer.add_figure("Normal Seq Error AF Spectra", normal_seq_error_spectra_fig, epoch)

                normal_artifact_spectra_fig, normal_artifact_spectra_axs = plot_artifact_spectra(self.normal_artifact_spectra)
                summary_writer.add_figure("Normal Artifact AF Spectra", normal_artifact_spectra_fig, epoch)

                var_spectra_fig, var_spectra_axs = plt.subplots()
                frac, dens = self.somatic_spectrum.spectrum_density_vs_fraction()
                var_spectra_axs.plot(frac.detach().numpy(), dens.detach().numpy(), label="spectrum")
                var_spectra_axs.set_title("Variant AF Spectrum")
                summary_writer.add_figure("Variant AF Spectra", var_spectra_fig, epoch)

                # bar plot of log priors -- data is indexed by call type name, and x ticks are variant types
                log_prior_bar_plot_data = defaultdict(list)
                for var_type_idx, variant_type in enumerate(Variation):
                    log_priors = torch.nn.functional.log_softmax(self.make_unnormalized_priors(torch.LongTensor([var_type_idx]).to(device=self._device, dtype=self._dtype), torch.Tensor([0.001])), dim=-1)
                    log_priors_cpu = log_priors.squeeze().detach().cpu()
                    for call_type in (Call.SOMATIC, Call.ARTIFACT, Call.NORMAL_ARTIFACT):
                        log_prior_bar_plot_data[call_type.name].append(log_priors_cpu[call_type])

                prior_fig, prior_ax = plotting.grouped_bar_plot(log_prior_bar_plot_data, [v_type.name for v_type in Variation], "log priors")
                summary_writer.add_figure("log priors", prior_fig, epoch)

                # normal artifact joint tumor-normal spectra
                # na_fig, na_axes = plt.subplots(1, len(Variation), sharex='all', sharey='all', squeeze=False)
                # for variant_index, variant_type in enumerate(Variation):
                #    self.normal_artifact_spectra[variant_index].density_plot_on_axis(na_axes[0, variant_index])
                # plotting.tidy_subplots(na_fig, na_axes, x_label="tumor fraction", y_label="normal fraction",
                #                       row_labels=[""], column_labels=[var_type.name for var_type in Variation])
                # summary_writer.add_figure("normal artifact spectra", na_fig, epoch)

    # map of Variant type to probability threshold that maximizes F1 score
    # loader is a Dataloader whose collate_fn is the PosteriorBatch constructor
    def calculate_probability_thresholds(self, loader, summary_writer: SummaryWriter = None, germline_mode: bool = False):
        self.train(False)
        error_probs_by_type = {var_type: [] for var_type in Variation}   # includes both artifact and seq errors

        error_probs_by_type_by_cnt = {var_type: [[] for _ in range(NUM_COUNT_BINS)] for var_type in Variation}

        pbar = tqdm(enumerate(loader), mininterval=10)
        for n, batch_cpu in pbar:
            batch = batch_cpu.copy_to(self._device, self._dtype, non_blocking=self._device.type == 'cuda')
            alt_counts = batch_cpu.get_alt_counts().tolist()
            # 0th column is true variant, subtract it from 1 to get error prob
            error_probs = self.error_probabilities(batch, germline_mode).cpu().tolist()

            for var_type, alt_count, error_prob in zip(batch_cpu.get_variant_types().tolist(), alt_counts, error_probs):
                error_probs_by_type[var_type].append(error_prob)
                error_probs_by_type_by_cnt[var_type][multiple_of_three_bin_index(min(alt_count, MAX_COUNT))].append(error_prob)

        thresholds_by_type = {}
        roc_fig, roc_axes = plt.subplots(1, len(Variation), sharex='all', sharey='all', squeeze=False)
        roc_by_cnt_fig, roc_by_cnt_axes = plt.subplots(1, len(Variation), sharex='all', sharey='all', squeeze=False, figsize=(10, 6), dpi=100)
        for var_type in Variation:
            # plot all count ROC curves for this variant type
            count_bin_labels = [str(multiple_of_three_bin_index_to_count(count_bin)) for count_bin in range(NUM_COUNT_BINS)]
            _ = plotting.plot_theoretical_roc_on_axis(error_probs_by_type_by_cnt[var_type], count_bin_labels, roc_by_cnt_axes[0, var_type])
            best_threshold = plotting.plot_theoretical_roc_on_axis([error_probs_by_type[var_type]], [""], roc_axes[0, var_type])[0][0]

            # TODO: the theoretical ROC might need to return the best threshold for this
            thresholds_by_type[var_type] = best_threshold

        variation_types = [var_type.name for var_type in Variation]
        plotting.tidy_subplots(roc_by_cnt_fig, roc_by_cnt_axes, x_label="sensitivity", y_label="precision",
                               row_labels=[""], column_labels=variation_types)
        plotting.tidy_subplots(roc_fig, roc_axes, x_label="sensitivity", y_label="precision",
                               row_labels=[""], column_labels=variation_types)
        if summary_writer is not None:
            summary_writer.add_figure("theoretical ROC by variant type ", roc_fig)
            summary_writer.add_figure("theoretical ROC by variant type and alt count ", roc_by_cnt_fig)

        return thresholds_by_type

    def update_priors_m_step(self, posteriors_nc, types_n, ignored_to_non_ignored_ratio):
        # update the priors in an EM-style M step.  We'll need the counts of each call type vs variant type
        total_nonignored = torch.sum(posteriors_nc).item()
        total_ignored = ignored_to_non_ignored_ratio * total_nonignored
        overall_total = total_ignored + total_nonignored

        with torch.no_grad():
            for c, call_type in enumerate(Call):
                if call_type == Call.SEQ_ERROR or call_type == Call.GERMLINE:
                    continue
                posteriors_n = posteriors_nc[:, c]

                for t, var_type in enumerate(Variation):
                    var_type_mask = (types_n == t)
                    total_for_this_call_and_var_type = torch.sum(posteriors_n * var_type_mask)

                    self._unnormalized_priors_vc[t, c] = torch.log(
                        total_for_this_call_and_var_type / (total_for_this_call_and_var_type + overall_total)).item()

            self._unnormalized_priors_vc[:, Call.SEQ_ERROR] = 0
            self._unnormalized_priors_vc[:, Call.GERMLINE] = -9999 if self.no_germline_mode else 0
