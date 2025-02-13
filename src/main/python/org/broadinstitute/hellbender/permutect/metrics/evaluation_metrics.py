import math
from collections import defaultdict
from typing import List

import numpy as np
import torch
from matplotlib import pyplot as plt
from torch.utils.tensorboard import SummaryWriter

from permutect.data.base_datum import BaseBatch, ArtifactBatch
from permutect.metrics import plotting
from permutect.utils import Variation, Call, Epoch, StreamingAverage

MAX_COUNT = 18  # counts above this will be truncated
MAX_LOGIT = 15
NUM_DATA_FOR_TENSORBOARD_PROJECTION = 10000


def round_up_to_nearest_three(x: int):
    return math.ceil(x / 3) * 3


def multiple_of_three_bin_index(x: int):
    return (round_up_to_nearest_three(x)//3) - 1    # -1 because zero is not a bin


MAX_BIN = multiple_of_three_bin_index(MAX_COUNT)

def multiple_of_three_bin_indices(counts: torch.Tensor):
    return (torch.ceil(counts/3) - 1).int()


def multiple_of_three_bin_index_to_count(idx: int):
    return 3 * (idx + 1)


# round logit to nearest int, truncate to range, ending up with bins 0. . . 2*max_logit
def logit_to_bin(logit):
    return min(max(round(logit), -MAX_LOGIT), MAX_LOGIT) + MAX_LOGIT


def bin_center(bin_idx):
    return bin_idx - MAX_LOGIT


NUM_COUNT_BINS = round_up_to_nearest_three(MAX_COUNT) // 3    # zero is not a bin


def make_count_bin_mask(bin_index: int, counts: torch.Tensor):
    assert bin_index < NUM_COUNT_BINS
    count_bin_bottom = 3*bin_index + 1
    count_bin_top = 3*bin_index + 3
    return (count_bin_bottom <= counts) * (counts <= count_bin_top)


# simple container class for holding results of the posterior model and other things that get output to the VCF and
# tensorboard analysis
class PosteriorResult:
    def __init__(self, artifact_logit: float, posterior_probabilities, log_priors, spectra_lls, normal_lls, label, alt_count, depth, var_type, embedding):
        self.artifact_logit = artifact_logit
        self.posterior_probabilities = posterior_probabilities
        self.log_priors = log_priors
        self.spectra_lls = spectra_lls
        self.normal_lls = normal_lls
        self.label = label
        self.alt_count = alt_count
        self.depth = depth
        self.variant_type = var_type
        self.embedding = embedding


# keep track of losses during training of artifact model
class LossMetrics:
    def __init__(self):
        self.labeled_loss = StreamingAverage()
        self.unlabeled_loss = StreamingAverage()

        self.labeled_loss_by_type = {variant_type: StreamingAverage() for variant_type in Variation}
        self.labeled_loss_by_count = {bin_idx: StreamingAverage() for bin_idx in range(NUM_COUNT_BINS)}

    def get_labeled_loss(self) -> float:
        return self.labeled_loss.get()

    def get_unlabeled_loss(self) -> float:
        return self.unlabeled_loss.get()

    def write_to_summary_writer(self, epoch_type: Epoch, epoch: int, summary_writer: SummaryWriter, prefix: str = ""):
        if not self.labeled_loss.is_empty():
            summary_writer.add_scalar(prefix + epoch_type.name + "/Labeled Loss", self.labeled_loss.get(), epoch)

        if not self.unlabeled_loss.is_empty():
            summary_writer.add_scalar(prefix + epoch_type.name + "/Unlabeled Loss", self.unlabeled_loss.get(), epoch)

        for bin_idx, loss in self.labeled_loss_by_count.items():
            if not loss.is_empty():
                summary_writer.add_scalar(
                    prefix + epoch_type.name + "/Labeled Loss/By Count/" + str(multiple_of_three_bin_index_to_count(bin_idx)), loss.get(), epoch)

        for var_type, loss in self.labeled_loss_by_type.items():
            if not loss.is_empty():
                summary_writer.add_scalar(prefix + epoch_type.name + "/Labeled Loss/By Type/" + var_type.name, loss.get(), epoch)

    # record the losses (indexed by batch dimension) by type and count, as well as the total loss not stratified by type and count
    # input losses are NOT weighted, but when recorded they are multiplied by weights if given
    # losses are divided into labeled and unlabeled
    # TODO: put type hint batch: BaseBatch | ArtifactBatch once docker update
    def record_losses(self, losses: torch.Tensor, batch, weights: torch.Tensor):
        # handle total loss
        labeled_weights, unlabeled_weights = batch.get_is_labeled_mask() * weights, (1 - batch.get_is_labeled_mask()) * weights

        self.labeled_loss.record_with_weights(losses, labeled_weights)
        self.unlabeled_loss.record_with_weights(losses, unlabeled_weights)

        # Note that we currently do not track unlabeled loss by type or by count
        # by type
        variant_types = batch.get_variant_types()

        # weight for losses is product of 1) the weights 2) the is_labeled mask, 3) the variant type mask
        for var_type_idx, var_type in enumerate(Variation):
            variant_type_mask = (variant_types == var_type_idx)
            self.labeled_loss_by_type[var_type].record_with_weights(losses, labeled_weights * variant_type_mask)

        # by count
        if isinstance(batch, BaseBatch):
            # rather than individually record each count, and therefore send lots of stuff off the GPU, we
            # send everything with the same bin simultaneously
            bins = multiple_of_three_bin_indices(batch.get_alt_counts())
            for count_bin in range(MAX_BIN + 1):
                indices = (bins == count_bin)
                self.labeled_loss_by_count[count_bin].record_with_weights(losses[indices], labeled_weights[indices])

        elif isinstance(batch, ArtifactBatch):
            for count_bin_index in range(NUM_COUNT_BINS):
                count_bin_mask = make_count_bin_mask(count_bin_index, batch.get_alt_counts())
                self.labeled_loss_by_count[count_bin_index].record_with_weights(losses, labeled_weights * count_bin_mask)


# predictions_and_labels is list of (predicted logit, actual label) tuples
# adjustment is the logit threshold that maximizes accuracy -- basically we're trying to find the shift such that
# a logit of 0 expresses 50/50 confidence
# output is the amount to be SUBTRACTED from logit to get a final adjusted logit
def calculate_logit_adjustment(predictions_and_labels, use_harmonic_mean: bool = False):
    _, adjustment = plotting.get_roc_data(predictions_and_labels, given_threshold=None, sens_prec=False, use_harmonic_mean=use_harmonic_mean)
    return adjustment


class EvaluationMetricsForOneEpochType:
    def __init__(self):
        # indexed by variant type, then count bin, then logit bin
        self.acc_vs_logit = {
            var_type: [[StreamingAverage() for _ in range(2 * MAX_LOGIT + 1)] for _ in range(NUM_COUNT_BINS)] for
            var_type in Variation}

        self.acc_vs_logit_all_counts = {
            var_type: [StreamingAverage() for _ in range(2 * MAX_LOGIT + 1)] for var_type in Variation}

        # indexed by variant type, then call type (artifact vs variant), then count bin
        self.acc_vs_cnt = {var_type: defaultdict(lambda: [StreamingAverage() for _ in range(NUM_COUNT_BINS)]) for
                      var_type in Variation}

        # variant type -> (predicted logit, actual label)
        self.roc_data = {var_type: [] for var_type in Variation}

        # variant type, count -> (predicted logit, actual label)
        self.roc_data_by_cnt = {var_type: [[] for _ in range(NUM_COUNT_BINS)] for var_type in Variation}

    # Variant is an IntEnum, so variant_type can also be integer
    # label is 1 for artifact / error; 0 for non-artifact / true variant
    # correct_call is boolean -- was the prediction correct?
    # the predicted logit is the logit corresponding to the predicted probability that call in question is an artifact / error
    def record_call(self, variant_type: Variation, predicted_logit: float, label: float, correct_call, alt_count: int, weight: float = 1.0):
        count_bin_index = multiple_of_three_bin_index(min(MAX_COUNT, alt_count))
        self.acc_vs_cnt[variant_type][Call.SOMATIC if label < 0.5 else Call.ARTIFACT][count_bin_index].record(correct_call, weight)
        self.acc_vs_logit[variant_type][count_bin_index][logit_to_bin(predicted_logit)].record(correct_call, weight)
        self.acc_vs_logit_all_counts[variant_type][logit_to_bin(predicted_logit)].record(correct_call, weight)

        self.roc_data[variant_type].append((predicted_logit, label))
        self.roc_data_by_cnt[variant_type][count_bin_index].append((predicted_logit, label))

    # return a list of tuples.  This outer list is over the two labels, Call.SOMATIC and Call.ARTIFACT.  Each tuple consists of
    # (list of alt counts (x axis), list of accuracies (y axis), the label)
    def make_data_for_accuracy_plot(self, var_type: Variation):
        non_empty_count_bins_by_label = {
            label: [idx for idx in range(NUM_COUNT_BINS) if not self.acc_vs_cnt[var_type][label][idx].is_empty()]
            for label in self.acc_vs_cnt[var_type].keys()}

        return [([multiple_of_three_bin_index_to_count(idx) for idx in non_empty_count_bins_by_label[label]],
                    [self.acc_vs_cnt[var_type][label][idx].get() for idx in non_empty_count_bins_by_label[label]],
                    label.name) for label in self.acc_vs_cnt[var_type].keys()]

    # similar tuple format but now it's (list of logits, list of accuracies, count)
    def make_data_for_calibration_plot(self, var_type: Variation):
        non_empty_logit_bins = [
            [idx for idx in range(2 * MAX_LOGIT + 1) if not self.acc_vs_logit[var_type][count_idx][idx].is_empty()]
            for count_idx in range(NUM_COUNT_BINS)]
        return [([bin_center(idx) for idx in non_empty_logit_bins[count_idx]],
                                        [self.acc_vs_logit[var_type][count_idx][idx].get() for idx in
                                         non_empty_logit_bins[count_idx]],
                                        str(multiple_of_three_bin_index_to_count(count_idx))) for count_idx in
                                       range(NUM_COUNT_BINS)]

    # now it's (list of logits, list of accuracies)
    def make_data_for_calibration_plot_all_counts(self, var_type: Variation):
        non_empty_logit_bins = [idx for idx in range(2 * MAX_LOGIT + 1) if not self.acc_vs_logit_all_counts[var_type][idx].is_empty()]
        return ([bin_center(idx) for idx in non_empty_logit_bins],
                    [self.acc_vs_logit_all_counts[var_type][idx].get() for idx in non_empty_logit_bins])

    def plot_accuracy(self, var_type: Variation, axis):
        acc_vs_cnt_x_y_lab_tuples = self.make_data_for_accuracy_plot(var_type)
        plotting.simple_plot_on_axis(axis, acc_vs_cnt_x_y_lab_tuples, None, None)

    def plot_calibration(self, var_type: Variation, axis):
        acc_vs_logit_x_y_lab_tuples = self.make_data_for_calibration_plot(var_type)
        plotting.simple_plot_on_axis(axis, acc_vs_logit_x_y_lab_tuples, None, None)

    def plot_calibration_all_counts(self, var_type: Variation, axis):
        logits_list, accuracies_list = self.make_data_for_calibration_plot_all_counts(var_type)
        plotting.simple_plot_on_axis(axis, [(logits_list, accuracies_list, "calibration")], None, None)

    def plot_roc_curve(self, var_type: Variation, axis, given_threshold: float = None, sens_prec: bool = False):
        plotting.plot_accuracy_vs_accuracy_roc_on_axis([self.roc_data[var_type]], [None], axis, given_threshold, sens_prec)

    def plot_roc_curves_by_count(self, var_type: Variation, axis, given_threshold: float = None, sens_prec: bool = False):
        plotting.plot_accuracy_vs_accuracy_roc_on_axis(self.roc_data_by_cnt[var_type],
                                                       [str(multiple_of_three_bin_index_to_count(idx)) for idx in
                                                        range(NUM_COUNT_BINS)], axis, given_threshold, sens_prec)

    # return variant type, count bin -> logit adjustment to be subtracted (so that maximum accuracy is at threshold of logit = 0)
    def calculate_logit_adjustments(self, use_harmonic_mean: bool = False):
        result = {var_type: [0.0 for _ in range(NUM_COUNT_BINS)] for var_type in Variation}
        for var_type in Variation:
            for cbin in range(NUM_COUNT_BINS):
                data = self.roc_data_by_cnt[var_type][cbin]
                if data:    # leave adjustment at 0 if no data
                    result[var_type][cbin] = calculate_logit_adjustment(data, use_harmonic_mean)

        return result


class EvaluationMetrics:
    def __init__(self):
        # we will have a map from epoch type to EvaluationMetricsForOneEpochType
        self.metrics = defaultdict(EvaluationMetricsForOneEpochType)

        # list of (PosteriorResult, Call) tuples
        self.mistakes = []

    # Variant is an IntEnum, so variant_type can also be integer
    # label is 1 for artifact / error; 0 for non-artifact / true variant
    # correct_call is boolean -- was the prediction correct?
    # the predicted logit is the logit corresponding to the predicted probability that call in question is an artifact / error
    def record_call(self, epoch_type: Epoch, variant_type: Variation, predicted_logit: float, label: float, correct_call, alt_count: int, weight: float = 1.0):
        self.metrics[epoch_type].record_call(variant_type, predicted_logit, label, correct_call, alt_count, weight)

    # track bad calls when filtering is given an optional evaluation truth VCF
    def record_mistake(self, posterior_result: PosteriorResult, call: Call):
        self.mistakes.append((posterior_result, call))

    def make_mistake_histograms(self, summary_writer: SummaryWriter):
        # indexed by call then var_type, inner is a list of posterior results with that call and var type
        posterior_result_mistakes_by_call_and_var_type = defaultdict(lambda: defaultdict(list))
        for posterior_result, call in self.mistakes:
            posterior_result_mistakes_by_call_and_var_type[call][posterior_result.variant_type].append(posterior_result)

        mistake_calls = posterior_result_mistakes_by_call_and_var_type.keys()
        num_rows = len(mistake_calls)

        af_fig, af_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='none', squeeze=False)
        logit_fig, logit_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='none', squeeze=False)
        ac_fig, ac_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='none', squeeze=False)
        prob_fig, prob_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='none', squeeze=False)

        for row_idx, mistake_call in enumerate(mistake_calls):
            for var_type in Variation:
                posterior_results = posterior_result_mistakes_by_call_and_var_type[mistake_call][var_type]

                af_data = [pr.alt_count / pr.depth for pr in posterior_results]
                plotting.simple_histograms_on_axis(af_axes[row_idx, var_type], [af_data], [""], 20)

                ac_data = [pr.alt_count for pr in posterior_results]
                plotting.simple_histograms_on_axis(ac_axes[row_idx, var_type], [ac_data], [""], 20)

                logit_data = [pr.artifact_logit for pr in posterior_results]
                plotting.simple_histograms_on_axis(logit_axes[row_idx, var_type], [logit_data], [""], 20)

                # posterior probability assigned to this incorrect call
                prob_data = [pr.posterior_probabilities[mistake_call] for pr in posterior_results]
                plotting.simple_histograms_on_axis(prob_axes[row_idx, var_type], [prob_data], [""], 20)

        variation_types = [var_type.name for var_type in Variation]
        row_names = [mistake.name for mistake in mistake_calls]

        plotting.tidy_subplots(af_fig, af_axes, x_label="alt allele fraction", y_label="", row_labels=row_names, column_labels=variation_types)
        plotting.tidy_subplots(ac_fig, ac_axes, x_label="alt count", y_label="", row_labels=row_names,
                               column_labels=variation_types)
        plotting.tidy_subplots(logit_fig, logit_axes, x_label="artifact logit", y_label="", row_labels=row_names,
                               column_labels=variation_types)
        plotting.tidy_subplots(prob_fig, prob_axes, x_label="mistake call probability", y_label="", row_labels=row_names,
                               column_labels=variation_types)

        summary_writer.add_figure("mistake allele fractions", af_fig)
        summary_writer.add_figure("mistake alt counts", ac_fig)
        summary_writer.add_figure("mistake artifact logits", logit_fig)
        summary_writer.add_figure("probability assigned to mistake calls", prob_fig)

    def make_plots(self, summary_writer: SummaryWriter, given_thresholds=None, sens_prec: bool = False, epoch: int = None):
        # given_thresholds is a dict from Variation to float (logit-scaled) used in the ROC curves
        keys = self.metrics.keys()
        num_rows = len(keys)
        # grid of figures -- rows are epoch types, columns are variant types
        # each subplot has two line graphs of accuracy vs alt count, one each for artifact, non-artifact
        acc_vs_cnt_fig, acc_vs_cnt_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='all', squeeze=False)
        roc_fig, roc_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='all', squeeze=False, figsize=(4 * len(Variation), 4 * len(keys)), dpi=200)
        cal_fig, cal_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='all', squeeze=False)
        cal_fig_all_counts, cal_axes_all_counts = plt.subplots(num_rows, len(Variation), sharex='all', sharey='all', squeeze=False)
        roc_by_cnt_fig, roc_by_cnt_axes = plt.subplots(num_rows, len(Variation), sharex='all', sharey='all', squeeze=False, figsize=(4 * len(Variation), 4 * len(keys)), dpi=200)

        for row_idx, key in enumerate(keys):
            metric = self.metrics[key]
            for var_type in Variation:
                given_threshold = None if given_thresholds is None else given_thresholds[var_type]
                metric.plot_accuracy(var_type, acc_vs_cnt_axes[row_idx, var_type])
                metric.plot_calibration(var_type, cal_axes[row_idx, var_type])
                metric.plot_calibration_all_counts(var_type, cal_axes_all_counts[row_idx, var_type])
                metric.plot_roc_curve(var_type, roc_axes[row_idx, var_type], given_threshold, sens_prec)
                metric.plot_roc_curves_by_count(var_type, roc_by_cnt_axes[row_idx, var_type], given_threshold, sens_prec)
        # done collecting stats for all loaders and filling in subplots

        nonart_label = "sensitivity" if sens_prec else "non-artifact accuracy"
        art_label = "precision" if sens_prec else "artifact accuracy"

        variation_types = [var_type.name for var_type in Variation]
        row_names = [epoch_type.name for epoch_type in self.metrics.keys()]
        plotting.tidy_subplots(acc_vs_cnt_fig, acc_vs_cnt_axes, x_label="alt count", y_label="accuracy", row_labels=row_names, column_labels=variation_types)
        plotting.tidy_subplots(roc_fig, roc_axes, x_label=nonart_label, y_label=art_label, row_labels=row_names, column_labels=variation_types)
        plotting.tidy_subplots(roc_by_cnt_fig, roc_by_cnt_axes, x_label=nonart_label, y_label=art_label, row_labels=row_names, column_labels=variation_types)
        plotting.tidy_subplots(cal_fig, cal_axes, x_label="predicted logit", y_label="accuracy", row_labels=row_names, column_labels=variation_types)
        plotting.tidy_subplots(cal_fig_all_counts, cal_axes_all_counts, x_label="predicted logit", y_label="accuracy", row_labels=row_names, column_labels=variation_types)

        summary_writer.add_figure("accuracy by alt count", acc_vs_cnt_fig, global_step=epoch)
        summary_writer.add_figure(" accuracy by logit output by count", cal_fig, global_step=epoch)
        summary_writer.add_figure(" accuracy by logit output", cal_fig_all_counts, global_step=epoch)
        summary_writer.add_figure("sensitivity vs precision" if sens_prec else "variant accuracy vs artifact accuracy", roc_fig, global_step=epoch)
        summary_writer.add_figure("sensitivity vs precision by alt count" if sens_prec else "variant accuracy vs artifact accuracy by alt count", roc_by_cnt_fig, global_step=epoch)


def sample_indices_for_tensorboard(indices: List[int]):
    indices_np = np.array(indices)

    if len(indices_np) <= NUM_DATA_FOR_TENSORBOARD_PROJECTION:
        return indices_np

    idx = np.random.choice(len(indices_np), size=NUM_DATA_FOR_TENSORBOARD_PROJECTION, replace=False)
    return indices_np[idx]


class EmbeddingMetrics:
    TRUE_POSITIVE = "true-positive"
    FALSE_POSITIVE = "false-positive"
    TRUE_NEGATIVE_ARTIFACT = "true-negative-artifact"   # distinguish these because artifact and eg germline should embed differently
    TRUE_NEGATIVE_NONARTIFACT = "true-negative-nonartifact"
    TRUE_NEGATIVE = "true-negative"
    FALSE_NEGATIVE_ARTIFACT = "false-negative-artifact"
    TRUE_NEGATIVE_SEQ_ERROR = "true-negative-seq-error"

    def __init__(self):
        # things we will collect for the projections
        self.label_metadata = []  # list (extended by each batch) 1 if artifact, 0 if not
        self.correct_metadata = []  # list (extended by each batch), 1 if correct prediction, 0 if not
        self.type_metadata = []  # list of lists, strings of variant type
        self.truncated_count_metadata = []  # list of lists
        self.representations = []  # list of 2D tensors (to be stacked into a single 2D tensor), representations over batches

    def output_to_summary_writer(self, summary_writer: SummaryWriter, prefix: str = "", is_filter_variants: bool = False, epoch: int = None):
        # downsample to a reasonable amount of UMAP data
        all_metadata = list(zip(self.label_metadata, self.correct_metadata, self.type_metadata, self.truncated_count_metadata))

        indices_by_correct_status = defaultdict(list)

        for n, correct_status in enumerate(self.correct_metadata):
            indices_by_correct_status[correct_status].append(n)

        # note that if we don't have labeled truth, everything is boring
        all_indices = set(range(len(all_metadata)))
        interesting_indices = set(indices_by_correct_status[EmbeddingMetrics.TRUE_POSITIVE] +
                                       indices_by_correct_status[EmbeddingMetrics.FALSE_POSITIVE] +
                                       indices_by_correct_status[EmbeddingMetrics.FALSE_NEGATIVE_ARTIFACT])
        boring_indices = all_indices - interesting_indices

        '''if is_filter_variants:
            boring_indices = np.array(indices_by_correct_status["unknown"] + indices_by_correct_status[EmbeddingMetrics.TRUE_NEGATIVE_ARTIFACT])

            # if we have labeled truth, keep a few "boring" true negatives around; otherwise we only have "unknown"s
            boring_count = len(interesting_indices) // 3 if len(interesting_indices) > 0 else len(boring_indices)
            boring_to_keep = boring_indices[np.random.choice(len(boring_indices), size=boring_count, replace=False)]
            idx = np.hstack((boring_to_keep, interesting_indices))

        idx = np.random.choice(len(all_metadata), size=min(NUM_DATA_FOR_TENSORBOARD_PROJECTION, len(all_metadata)), replace=False)
'''

        stacked_representations = torch.vstack(self.representations)

        # read average embeddings stratified by variant type
        for variant_type in Variation:
            variant_name = variant_type.name
            indices = set([n for n, type_name in enumerate(self.type_metadata) if type_name == variant_name])

            interesting = interesting_indices & indices
            boring = boring_indices & indices
            boring_count = max(len(interesting) // 3, 100) if is_filter_variants else len(boring)
            boring_to_keep = np.array([int(n) for n in boring])[np.random.choice(len(boring), size=boring_count, replace=False)]
            idx = sample_indices_for_tensorboard(np.hstack((boring_to_keep, np.array([int(n) for n in interesting]))))

            summary_writer.add_embedding(stacked_representations[idx],
                                         metadata=[all_metadata[round(n)] for n in idx.tolist()],
                                         metadata_header=["Labels", "Correctness", "Types", "Counts"],
                                         tag=prefix+"embedding for variant type " + variant_name, global_step=epoch)

        # read average embeddings stratified by alt count
        for count_bin in range(NUM_COUNT_BINS):
            count = multiple_of_three_bin_index_to_count(count_bin)
            indices = set([n for n, alt_count in enumerate(self.truncated_count_metadata) if alt_count == str(count)])
            interesting = interesting_indices & indices
            boring = boring_indices & indices
            boring_count = max(len(interesting) // 3, 100) if is_filter_variants else len(boring)
            boring_to_keep = np.array([int(n) for n in boring])[np.random.choice(len(boring), size=boring_count, replace=False)]
            idx = sample_indices_for_tensorboard(np.hstack((boring_to_keep, np.array([int(n) for n in interesting]))))

            if len(idx) > 0:
                summary_writer.add_embedding(stacked_representations[idx],
                                        metadata=[all_metadata[round(n)] for n in idx.tolist()],
                                        metadata_header=["Labels", "Correctness", "Types", "Counts"],
                                        tag=prefix+"embedding for alt count " + str(count), global_step=epoch)



