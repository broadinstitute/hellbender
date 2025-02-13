# bug before PyTorch 1.7.1 that warns when constructing ParameterList
import math
import time
import warnings
from collections import defaultdict
from typing import List

import psutil
import torch
from torch import nn, Tensor
import numpy as np
from torch.utils.tensorboard import SummaryWriter
from queue import PriorityQueue


from tqdm.autonotebook import trange, tqdm
from itertools import chain
from matplotlib import pyplot as plt

from permutect.architecture.base_model import calculate_batch_weights, BaseModel, base_model_from_saved_dict, calculate_batch_source_weights
from permutect.architecture.gradient_reversal.module import GradientReversal
from permutect.architecture.mlp import MLP
from permutect.architecture.monotonic import MonoDense
from permutect.data.base_datum import ArtifactBatch, DEFAULT_GPU_FLOAT, DEFAULT_CPU_FLOAT
from permutect.data.artifact_dataset import ArtifactDataset
from permutect import utils, constants
from permutect.metrics.evaluation_metrics import LossMetrics, EvaluationMetrics, MAX_COUNT, round_up_to_nearest_three, \
    EmbeddingMetrics, multiple_of_three_bin_index_to_count, multiple_of_three_bin_index
from permutect.parameters import TrainingParameters, ArtifactModelParameters
from permutect.utils import Variation, Epoch, Label
from permutect.metrics import plotting

warnings.filterwarnings("ignore", message="Setting attributes on ParameterList is not supported.")


WORST_OFFENDERS_QUEUE_SIZE = 100


def effective_count(weights: Tensor):
    return (torch.square(torch.sum(weights)) / torch.sum(torch.square(weights))).item()


# group rows into consecutive chunks to yield a 3D tensor, average over dim=1 to get
# 2D tensor of sums within each chunk
def sums_over_chunks(tensor2d: Tensor, chunk_size: int):
    assert len(tensor2d) % chunk_size == 0
    return torch.sum(tensor2d.reshape([len(tensor2d) // chunk_size, chunk_size, -1]), dim=1)


class Calibration(nn.Module):

    def __init__(self, hidden_layer_sizes: List[int]):
        super(Calibration, self).__init__()

        # calibration takes [logit, ref count, alt count] as input and maps it to [calibrated logit]
        # it is monotonically increasing in logit, unconstrained in ref and alt count
        # we initialize it to calibrated logit = input logit

        # likewise, we cap the effective alt and ref counts and input logits to avoid arbitrarily large confidence
        self.max_alt = nn.Parameter(torch.tensor(20.0))
        self.max_ref = nn.Parameter(torch.tensor(20.0))
        self.max_input_logit = nn.Parameter(torch.tensor(20.0))

        center_spacing = 1
        ref_center_spacing = 5

        # centers of Gaussian comb featurizations
        # note: even though they aren't learned and requires_grad is False, we still wrap them in nn.Parameter
        # so that they can be sent to GPU recursively when the grandparent ArtifactModel is
        self.alt_centers = nn.Parameter(torch.arange(start=1, end=20, step=center_spacing), requires_grad=False)
        self.ref_centers = nn.Parameter(torch.arange(start=1, end=20, step=ref_center_spacing), requires_grad=False)

        # increasing in the 1st feature, logits
        # logit is one feature, then the Gaussian comb for alt and ref counts is the other
        self.monotonic = MonoDense(1 + len(self.ref_centers) + len(self.alt_centers), hidden_layer_sizes + [1], 1, 0)

        self.is_turned_on = True

        self.max_alt_count_for_adjustment = 20
        # after training we compute one final calibration adjustment, which depends on alt count
        # the nth element is the adjustment for alt count n
        # note that this is NOT a learnable parameter!!!! It is *set* but not learned!!
        self.final_adjustments = nn.Parameter(torch.zeros(self.max_alt_count_for_adjustment + 1), requires_grad=False)

    def set_adjustments(self, adjustments: torch.Tensor):
        current_device, current_dtype = self.final_adjustments.device, self.final_adjustments.dtype
        clipped_adjustments = adjustments[:len(self.final_adjustments)]
        padding_needed = len(self.final_adjustments) - len(clipped_adjustments)
        padded_adjustments = torch.hstack((clipped_adjustments, clipped_adjustments[-1] * torch.ones(padding_needed))) if padding_needed else clipped_adjustments
        self.final_adjustments = nn.Parameter(padded_adjustments.to(device=current_device, dtype=current_dtype), requires_grad=False)
        self.max_alt_count_for_adjustment = len(adjustments) - 1

    def calibrated_logits(self, logits_b: Tensor, ref_counts_b: Tensor, alt_counts_b: Tensor):
        if self.is_turned_on:
            logits_bc = torch.tanh(logits_b / self.max_input_logit)[:, None]

            ref_comb_bc = torch.softmax(-torch.square(ref_counts_b[:, None] - self.ref_centers[None, :]).float(), dim=1)
            alt_comb_bc = torch.softmax(-torch.square(alt_counts_b[:, None] - self.alt_centers[None, :]).float(), dim=1)
            input_2d = torch.hstack([logits_bc, ref_comb_bc, alt_comb_bc])
            calibrated_b = self.monotonic.forward(input_2d).squeeze()

            counts_for_adjustment = torch.clamp(alt_counts_b, max=self.max_alt_count_for_adjustment).long()
            adjustments = self.final_adjustments[counts_for_adjustment]

            return calibrated_b + adjustments
        else:   # should never happen
            return logits_b

    def forward(self, logits, ref_counts: Tensor, alt_counts: Tensor):
        return self.calibrated_logits(logits, ref_counts, alt_counts)

    def plot_calibration(self):
        device, dtype = self.final_adjustments.device, self.final_adjustments.dtype
        alt_counts = [1, 3, 5, 10, 15, 20]
        ref_counts = [1, 3, 5, 10, 15, 20]
        logits = torch.arange(-10, 10, 0.1, device=device, dtype=dtype)
        cal_fig,cal_axes = plt.subplots(len(alt_counts), len(ref_counts), sharex='all', sharey='all',
                                        squeeze=False, figsize=(10, 6), dpi=100)

        for row_idx, alt_count in enumerate(alt_counts):
            for col_idx, ref_count in enumerate(ref_counts):
                calibrated = self.forward(logits, ref_count * torch.ones_like(logits, device=device, dtype=dtype), alt_count * torch.ones_like(logits, device=device, dtype=dtype))
                plotting.simple_plot_on_axis(cal_axes[row_idx, col_idx], [(logits.detach().cpu(), calibrated.detach().cpu(), "")], None, None)

        plotting.tidy_subplots(cal_fig, cal_axes, x_label="alt count", y_label="ref count",
                               row_labels=[str(n) for n in ref_counts], column_labels=[str(n) for n in alt_counts])

        return cal_fig, cal_axes


class ArtifactModel(nn.Module):
    """
    aggregation_layers: dimensions of layers for aggregation, excluding its input which is determined by the
    representation model.

    output_layers: dimensions of layers after aggregation, excluding the output dimension,
    which is 1 for a single logit representing artifact/non-artifact.  This is not part of the aggregation layers
    because we have different output layers for each variant type.
    """

    def __init__(self, params: ArtifactModelParameters, num_base_features: int, num_ref_alt_features: int, device=utils.gpu_if_available()):
        super(ArtifactModel, self).__init__()

        self._device = device
        self._dtype = DEFAULT_GPU_FLOAT if device != torch.device("cpu") else DEFAULT_CPU_FLOAT
        self.num_base_features = num_base_features
        self.num_ref_alt_features = num_ref_alt_features
        self.params = params

        # feature layers before the domain adaptation source classifier splits from the artifact classifier
        self.feature_layers = MLP([num_base_features] + params.aggregation_layers, batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)

        # TODO: artifact classifier hidden layers are hard-coded!!!
        # The [1] is for the output logit
        self.artifact_classifier = MLP([self.feature_layers.output_dimension()] + [-1, -1, 1], batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)

        # one Calibration module for each variant type; that is, calibration depends on both count and type
        self.calibration = nn.ModuleList([Calibration(params.calibration_layers) for variant_type in Variation])

        self.to(device=self._device, dtype=self._dtype)

    def training_parameters(self):
        return chain(self.feature_layers.parameters(), self.artifact_classifier.parameters(), self.calibration.parameters())

    def calibration_parameters(self):
        return self.calibration.parameters()

    def freeze_all(self):
        utils.freeze(self.parameters())

    def set_epoch_type(self, epoch_type: utils.Epoch):
        if epoch_type == utils.Epoch.TRAIN:
            self.train(True)
            utils.freeze(self.parameters())
            utils.unfreeze(self.training_parameters())
        else:
            self.freeze_all()

    # returns 1D tensor of length batch_size of log odds ratio (logits) between artifact and non-artifact
    def forward(self, batch: ArtifactBatch):
        # batch has already gotten copy_to(self._device, self._dtype)
        features = self.feature_layers.forward(batch.get_representations_2d())
        uncalibrated_logits = self.artifact_classifier.forward(features).reshape(batch.size())
        calibrated_logits = torch.zeros_like(uncalibrated_logits, device=self._device)
        variant_types = batch.get_variant_types()
        for n, _ in enumerate(Variation):
            mask = (variant_types == n)
            calibrated_logits += mask * self.calibration[n].forward(uncalibrated_logits, batch.get_ref_counts(), batch.get_alt_counts())
        return calibrated_logits, uncalibrated_logits, features

    def learn(self, dataset: ArtifactDataset, training_params: TrainingParameters, summary_writer: SummaryWriter, validation_fold: int = None, epochs_per_evaluation: int = None):
        bce = nn.BCEWithLogitsLoss(reduction='none')  # no reduction because we may want to first multiply by weights for unbalanced data
        # cross entropy (with logit inputs) loss for adversarial source classification task
        ce = nn.CrossEntropyLoss(reduction='none')

        num_sources = len(dataset.counts_by_source.keys())
        if num_sources == 1:
            print("Training data come from a single source (this could be multiple files with the same source annotation applied in preprocessing)")
        else:
            sources_list = list(dataset.counts_by_source.keys())
            sources_list.sort()
            assert sources_list[0] == 0, "There is no source 0"
            assert sources_list[-1] == num_sources - 1, f"sources should be 0, 1, 2. . . without gaps, but sources are {sources_list}."

            print(f"Training data come from multiple sources, with counts {dataset.counts_by_source}.")
        source_classifier = MLP([self.feature_layers.output_dimension()] + [-1, -1, num_sources],
                                    batch_normalize=self.params.batch_normalize, dropout_p=self.params.dropout_p)
        source_classifier.to(device=self._device, dtype=self._dtype)
        source_gradient_reversal = GradientReversal(alpha=0.01)  # initialize as barely active
        source_gradient_reversal.to(device=self._device, dtype=self._dtype)

        # TODO: fused = is_cuda?
        train_optimizer = torch.optim.AdamW(chain(self.training_parameters(), source_classifier.parameters()), lr=training_params.learning_rate,
                                            weight_decay=training_params.weight_decay)
        train_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(train_optimizer, factor=0.2, patience=5,
            threshold=0.001, min_lr=(training_params.learning_rate / 100), verbose=True)

        for idx, variation_type in enumerate(utils.Variation):
            print(f"For variation type {variation_type.name}, there are {int(dataset.totals[-1][Label.ARTIFACT][idx].item())} \
                artifacts, {int(dataset.totals[-1][Label.VARIANT][idx].item())} \
                non-artifacts, and {int(dataset.totals[-1][Label.UNLABELED][idx].item())} unlabeled data.")

        is_cuda = self._device.type == 'cuda'
        print(f"Is CUDA available? {is_cuda}")

        validation_fold_to_use = (dataset.num_folds - 1) if validation_fold is None else validation_fold
        train_loader = dataset.make_data_loader(dataset.all_but_one_fold(validation_fold_to_use), training_params.batch_size, is_cuda, training_params.num_workers)
        print(f"Train loader created, memory usage percent: {psutil.virtual_memory().percent:.1f}")
        valid_loader = dataset.make_data_loader([validation_fold_to_use], training_params.inference_batch_size, is_cuda, training_params.num_workers)
        print(f"Validation loader created, memory usage percent: {psutil.virtual_memory().percent:.1f}")

        first_epoch, last_epoch = 1, training_params.num_epochs + training_params.num_calibration_epochs
        for epoch in trange(1, last_epoch + 1, desc="Epoch"):
            start_of_epoch = time.time()
            print(f"Epoch {epoch}, memory usage percent: {psutil.virtual_memory().percent:.1f}")
            is_calibration_epoch = epoch > training_params.num_epochs

            p = epoch - 1
            new_alpha = (2 / (1 + math.exp(-0.1 * p))) - 1
            source_gradient_reversal.set_alpha(new_alpha)

            for epoch_type in [utils.Epoch.TRAIN, utils.Epoch.VALID]:
                self.set_epoch_type(epoch_type)
                # in calibration epoch, freeze the model except for calibration
                if is_calibration_epoch and epoch_type == utils.Epoch.TRAIN:
                    utils.freeze(self.parameters())
                    utils.unfreeze(self.calibration_parameters())  # unfreeze calibration but everything else stays frozen

                loss_metrics = LossMetrics()    # based on calibrated logits
                source_prediction_loss_metrics = LossMetrics()  # based on calibrated logits
                uncalibrated_loss_metrics = LossMetrics()  # based on uncalibrated logits

                loader = train_loader if epoch_type == utils.Epoch.TRAIN else valid_loader
                loader_iter = iter(loader)

                next_batch_cpu = next(loader_iter)
                next_batch = next_batch_cpu.copy_to(self._device, self._dtype, non_blocking=is_cuda)

                pbar = tqdm(range(len(loader)), mininterval=60)
                for n in pbar:
                    # forward and backward pass on batch, which is the last iteration's prefetched "next_batch"
                    batch_cpu = next_batch_cpu
                    batch = next_batch

                    # Optimization: Asynchronously send the next batch to the device while the model does work
                    next_batch_cpu = next(loader_iter)
                    next_batch = next_batch_cpu.copy_to(self._device, self._dtype, non_blocking=is_cuda)

                    logits, precalibrated_logits, features = self.forward(batch)

                    # one-hot prediction of sources
                    if num_sources > 1:
                        # gradient reversal means parameters before the features try to maximize source prediction loss, i.e. features
                        # try to forget the source, while parameters after the features try to minimize it, i.e. they try
                        # to achieve the adversarial task of distinguishing sources
                        source_prediction_logits = source_classifier.forward(source_gradient_reversal(features))
                        source_prediction_probs = torch.nn.functional.softmax(source_prediction_logits, dim=-1)
                        source_prediction_targets = torch.nn.functional.one_hot(batch.get_sources().long(), num_sources)
                        source_prediction_losses = torch.sum(torch.square(source_prediction_probs - source_prediction_targets), dim=-1)

                        # TODO: always by count?
                        source_prediction_weights = calculate_batch_source_weights(batch_cpu, dataset, by_count=is_calibration_epoch)
                        source_prediction_weights = source_prediction_weights.to(device=self._device, dtype=self._dtype, non_blocking=True)
                    else:
                        source_prediction_losses = torch.zeros_like(logits, device=self._device)
                        source_prediction_weights = torch.zeros_like(logits, device=self._device)

                    # TODO: we need a parameter to control the relative weight of unlabeled loss to labeled loss
                    weights = calculate_batch_weights(batch_cpu, dataset, by_count=True)
                    weights = weights.to(device=self._device, dtype=self._dtype, non_blocking=True)

                    labels = batch.get_training_labels()
                    uncalibrated_cross_entropies = bce(precalibrated_logits, labels)
                    calibrated_cross_entropies = bce(logits, labels)
                    labeled_losses = batch.get_is_labeled_mask() * (uncalibrated_cross_entropies + calibrated_cross_entropies) / 2

                    # unlabeled loss: entropy regularization. We use the uncalibrated logits because otherwise entropy
                    # regularization simply biases calibration to be overconfident.
                    probabilities = torch.sigmoid(precalibrated_logits)
                    entropies = torch.nn.functional.binary_cross_entropy_with_logits(precalibrated_logits, probabilities, reduction='none')
                    unlabeled_losses = (1 - batch.get_is_labeled_mask()) * entropies

                    # these losses include weights and take labeled vs unlabeled into account
                    losses = (labeled_losses + unlabeled_losses) * weights + (source_prediction_losses * source_prediction_weights)
                    loss = torch.sum(losses)

                    # at this point, losses, weights are on GPU (if available), while metrics are on CPU
                    # if we have done things right, this is okay and record_losses handles GPU <--> CPU efficiently
                    loss_metrics.record_losses(calibrated_cross_entropies.detach(), batch, weights * batch.get_is_labeled_mask())
                    uncalibrated_loss_metrics.record_losses(uncalibrated_cross_entropies.detach(), batch, weights * batch.get_is_labeled_mask())
                    uncalibrated_loss_metrics.record_losses(entropies.detach(), batch, weights * (1 - batch.get_is_labeled_mask()))
                    source_prediction_loss_metrics.record_losses(source_prediction_losses.detach(), batch, source_prediction_weights)

                    # calibration epochs freeze the model up to calibration, so I wonder if a purely unlabeled batch
                    # would cause lack of gradient problems. . .
                    if epoch_type == utils.Epoch.TRAIN:
                        utils.backpropagate(train_optimizer, loss)

                # done with one epoch type -- training or validation -- for this epoch
                loss_metrics.write_to_summary_writer(epoch_type, epoch, summary_writer)
                source_prediction_loss_metrics.write_to_summary_writer(epoch_type, epoch, summary_writer, prefix="source prediction")
                uncalibrated_loss_metrics.write_to_summary_writer(epoch_type, epoch, summary_writer, prefix="uncalibrated")
                if epoch_type == utils.Epoch.TRAIN:
                    train_scheduler.step(loss_metrics.get_labeled_loss())

                print(f"Labeled loss for {epoch_type.name} epoch {epoch}: {loss_metrics.get_labeled_loss():.3f}")
                print(f"Unlabeled loss for {epoch_type.name} epoch {epoch}: {uncalibrated_loss_metrics.get_unlabeled_loss():.3f}")
                if num_sources > 1:
                    print(f"Adversarial source prediction loss on labeled data for {epoch_type.name} epoch {epoch}: {source_prediction_loss_metrics.get_labeled_loss():.3f}")
                    print(f"Adversarial source prediction loss on unlabeled data for {epoch_type.name} epoch {epoch}: {source_prediction_loss_metrics.get_unlabeled_loss():.3f}")
            # done with training and validation for this epoch
            print(f"End of epoch {epoch}, memory usage percent: {psutil.virtual_memory().percent:.1f}, time elapsed(s): {time.time() - start_of_epoch:.2f}")
            is_last = (epoch == last_epoch)
            if (epochs_per_evaluation is not None and epoch % epochs_per_evaluation == 0) or is_last:
                print(f"performing evaluation on epoch {epoch}")
                self.evaluate_model(epoch, dataset, train_loader, valid_loader, summary_writer, collect_embeddings=False, report_worst=False)
            if is_last:
                # collect data in order to do final calibration
                print("collecting data for final calibration")
                evaluation_metrics, _ = self.collect_evaluation_data(dataset, train_loader, valid_loader, report_worst=False)

                logit_adjustments_by_var_type_and_count_bin = evaluation_metrics.metrics[Epoch.VALID].calculate_logit_adjustments(use_harmonic_mean=False)
                print("here are the logit adjustments:")
                for var_type_idx, var_type in enumerate(Variation):
                    adjustments_by_count_bin = logit_adjustments_by_var_type_and_count_bin[var_type]
                    max_bin_idx = len(adjustments_by_count_bin) - 1
                    max_count = multiple_of_three_bin_index_to_count(max_bin_idx)
                    adjustments_by_count = torch.zeros(max_count + 1)
                    for count in range(max_count + 1):
                        bin_idx = multiple_of_three_bin_index(count)
                        # negative sign because these are subtractive adjustments
                        adjustments_by_count[count] = -adjustments_by_count_bin[bin_idx]
                    print(f"for variant type {var_type.name} the adjustments are ")
                    print(adjustments_by_count.tolist())
                    self.calibration[var_type_idx].set_adjustments(adjustments_by_count)

                # consider this an extra post-postprocessing/final calibration epoch, hence epoch+1
                print("doing one final evaluation after the last logit adjustment")
                self.evaluate_model(epoch + 1, dataset, train_loader, valid_loader, summary_writer, collect_embeddings=True, report_worst=True)

            # note that we have not learned the AF spectrum yet
        # done with training

    def evaluate_model_after_training(self, dataset: ArtifactDataset, batch_size, num_workers, summary_writer: SummaryWriter):
        train_loader = dataset.make_data_loader(dataset.all_but_the_last_fold(), batch_size, self._device.type == 'cuda', num_workers)
        valid_loader = dataset.make_data_loader(dataset.last_fold_only(), batch_size, self._device.type == 'cuda', num_workers)
        self.evaluate_model(None, dataset, train_loader, valid_loader, summary_writer, collect_embeddings=True, report_worst=True)

    @torch.inference_mode()
    def collect_evaluation_data(self, dataset: ArtifactDataset, train_loader, valid_loader, report_worst: bool):
        # the keys are tuples of (true label -- 1 for variant, 0 for artifact; rounded alt count)
        worst_offenders_by_truth_and_alt_count = defaultdict(lambda: PriorityQueue(WORST_OFFENDERS_QUEUE_SIZE))

        evaluation_metrics = EvaluationMetrics()
        epoch_types = [Epoch.TRAIN, Epoch.VALID]
        for epoch_type in epoch_types:
            assert epoch_type == Epoch.TRAIN or epoch_type == Epoch.VALID  # not doing TEST here
            loader = train_loader if epoch_type == Epoch.TRAIN else valid_loader
            pbar = tqdm(enumerate(loader), mininterval=60)
            for n, batch_cpu in pbar:
                batch = batch_cpu.copy_to(self._device, self._dtype, non_blocking=self._device.type == 'cuda')

                # these are the same weights used in training
                # TODO: we need a parameter to control the relative weight of unlabeled loss to labeled loss
                weights = calculate_batch_weights(batch_cpu, dataset, by_count=True)
                weights = weights.to(dtype=self._dtype)     # not sent to GPU!

                logits, _, _ = self.forward(batch)
                # logits are calculated on the GPU (when available), so we must detach AND send back to CPU (if applicable)
                pred = logits.detach().cpu()

                # note that for metrics we use batch_cpu
                labels = batch_cpu.get_training_labels()
                correct = ((pred > 0) == (labels > 0.5)).tolist()

                for variant_type, predicted_logit, label, is_labeled, correct_call, alt_count, variant, weight in zip(
                        batch_cpu.get_variant_types().tolist(), pred.tolist(), labels.tolist(), batch_cpu.get_is_labeled_mask().tolist(), correct,
                        batch_cpu.get_alt_counts().tolist(), batch_cpu.get_variants(), weights.tolist()):
                    if is_labeled < 0.5:    # we only evaluate labeled data
                        continue
                    evaluation_metrics.record_call(epoch_type, variant_type, predicted_logit, label, correct_call, alt_count, weight)
                    if report_worst and not correct_call:
                        rounded_count = round_up_to_nearest_three(alt_count)
                        label_name = Label.ARTIFACT.name if label > 0.5 else Label.VARIANT.name
                        confidence = abs(predicted_logit)

                        # the 0th aka highest priority element in the queue is the one with the lowest confidence
                        pqueue = worst_offenders_by_truth_and_alt_count[(label_name, rounded_count)]

                        # clear space if this confidence is more egregious
                        if pqueue.full() and pqueue.queue[0][0] < confidence:
                            pqueue.get()  # discards the least confident bad call

                        if not pqueue.full():  # if space was cleared or if it wasn't full already
                            pqueue.put((confidence, str(variant.contig) + ":" + str(
                                variant.position) + ':' + variant.ref + "->" + variant.alt))
            # done with this epoch type
        # done collecting data
        return evaluation_metrics, worst_offenders_by_truth_and_alt_count

    @torch.inference_mode()
    def evaluate_model(self, epoch: int, dataset: ArtifactDataset, train_loader, valid_loader, summary_writer: SummaryWriter,
                                      collect_embeddings: bool = False, report_worst: bool = False):

        # self.freeze_all()
        evaluation_metrics, worst_offenders_by_truth_and_alt_count = self.collect_evaluation_data(dataset, train_loader, valid_loader, report_worst)
        evaluation_metrics.make_plots(summary_writer, epoch=epoch)

        if report_worst:
            for (true_label, rounded_count), pqueue in worst_offenders_by_truth_and_alt_count.items():
                tag = "True label: " + true_label + ", rounded alt count: " + str(rounded_count)

                lines = []
                while not pqueue.empty():   # this goes from least to most egregious, FYI
                    confidence, var_string = pqueue.get()
                    lines.append(f"{var_string} ({confidence:.2f})")

                summary_writer.add_text(tag, "\n".join(lines), global_step=epoch)

        if collect_embeddings:
            embedding_metrics = EmbeddingMetrics()

            # now go over just the validation data and generate feature vectors / metadata for tensorboard projectors (UMAP)
            pbar = tqdm(enumerate(valid_loader), mininterval=60)

            for n, batch_cpu in pbar:
                batch = batch_cpu.copy_to(self._device, self._dtype, non_blocking=self._device.type == 'cuda')
                logits, _, _ = self.forward(batch)
                pred = logits.detach().cpu()
                labels = batch_cpu.get_training_labels()
                correct = ((pred > 0) == (labels > 0.5)).tolist()

                label_strings = [("artifact" if label > 0.5 else "non-artifact") if is_labeled > 0.5 else "unlabeled"
                                 for (label, is_labeled) in zip(labels.tolist(), batch_cpu.get_is_labeled_mask().tolist())]

                correct_strings = [str(correctness) if is_labeled > 0.5 else "-1"
                                 for (correctness, is_labeled) in zip(correct, batch_cpu.get_is_labeled_mask().tolist())]

                for (metrics, embedding) in [(embedding_metrics, batch_cpu.get_representations_2d().detach())]:
                    metrics.label_metadata.extend(label_strings)
                    metrics.correct_metadata.extend(correct_strings)
                    metrics.type_metadata.extend([Variation(idx).name for idx in batch_cpu.get_variant_types().tolist()])
                    metrics.truncated_count_metadata.extend([str(round_up_to_nearest_three(min(MAX_COUNT, alt_count))) for alt_count in batch_cpu.get_alt_counts().tolist()])
                    metrics.representations.append(embedding)
            embedding_metrics.output_to_summary_writer(summary_writer, epoch=epoch)
        # done collecting data

    def make_dict_for_saving(self, artifact_log_priors, artifact_spectra, prefix: str = "artifact"):
        return {(prefix + constants.STATE_DICT_NAME): self.state_dict(),
                (prefix + constants.NUM_BASE_FEATURES_NAME): self.num_base_features,
                (prefix + constants.NUM_REF_ALT_FEATURES_NAME): self.num_ref_alt_features,
                (prefix + constants.HYPERPARAMS_NAME): self.params,
                (prefix + constants.ARTIFACT_LOG_PRIORS_NAME): artifact_log_priors,
                (prefix + constants.ARTIFACT_SPECTRA_STATE_DICT_NAME): artifact_spectra.state_dict()}

    def save(self, path, artifact_log_priors, artifact_spectra, prefix: str = "artifact"):
        torch.save(self.make_dict_for_saving(artifact_log_priors, artifact_spectra, prefix), path)

    def save_with_base_model(self, base_model: BaseModel, path, artifact_log_priors, artifact_spectra):
        artifact_dict = self.make_dict_for_saving(artifact_log_priors, artifact_spectra, prefix="artifact")
        base_dict = base_model.make_dict_for_saving(prefix="base")
        torch.save({**artifact_dict, **base_dict}, path)


def artifact_model_from_saved_dict(saved, prefix: str = "artifact"):
    model_params = saved[prefix + constants.HYPERPARAMS_NAME]
    num_base_features = saved[prefix + constants.NUM_BASE_FEATURES_NAME]
    num_ref_alt_features = saved[prefix + constants.NUM_REF_ALT_FEATURES_NAME]
    model = ArtifactModel(model_params, num_base_features, num_ref_alt_features)
    model.load_state_dict(saved[prefix + constants.STATE_DICT_NAME])

    artifact_log_priors = saved[prefix + constants.ARTIFACT_LOG_PRIORS_NAME]  # possibly None
    artifact_spectra_state_dict = saved[prefix + constants.ARTIFACT_SPECTRA_STATE_DICT_NAME]  # possibly None
    return model, artifact_log_priors, artifact_spectra_state_dict


# log artifact priors and artifact spectra may be None
def load_artifact_model(path,  device, prefix: str = "artifact") -> ArtifactModel:
    saved = torch.load(path, map_location=device)
    return artifact_model_from_saved_dict(saved, prefix)


def load_base_model_and_artifact_model(path, device) -> ArtifactModel:
    saved = torch.load(path, map_location=device)
    base_model = base_model_from_saved_dict(saved, prefix="base")
    artifact_model, artifact_log_priors, artifact_spectra = artifact_model_from_saved_dict(saved, prefix="artifact")
    return base_model, artifact_model, artifact_log_priors, artifact_spectra

