import math
from abc import ABC, abstractmethod
from enum import Enum
from itertools import chain
import time
from typing import List

import psutil
import torch
import numpy as np
from torch.utils.tensorboard import SummaryWriter
from torch.nn.parameter import Parameter
from tqdm.autonotebook import trange, tqdm

from permutect import utils, constants
from permutect.architecture.dna_sequence_convolution import DNASequenceConvolution
from permutect.architecture.gated_mlp import GatedMLP, GatedRefAltMLP
from permutect.architecture.gradient_reversal.module import GradientReversal
from permutect.architecture.mlp import MLP
from permutect.data.base_datum import BaseBatch, DEFAULT_GPU_FLOAT, DEFAULT_CPU_FLOAT
from permutect.data.base_dataset import BaseDataset, ALL_COUNTS_SENTINEL
from permutect.metrics.evaluation_metrics import LossMetrics, EmbeddingMetrics, round_up_to_nearest_three, MAX_COUNT
from permutect.parameters import BaseModelParameters, TrainingParameters


# group rows into consecutive chunks to yield a 3D tensor, average over dim=1 to get
# 2D tensor of sums within each chunk
from permutect.utils import Variation, Label


def sums_over_chunks(tensor2d: torch.Tensor, chunk_size: int):
    assert len(tensor2d) % chunk_size == 0
    return torch.sum(tensor2d.reshape([len(tensor2d) // chunk_size, chunk_size, -1]), dim=1)


# note: this works for both BaseBatch/BaseDataset AND ArtifactBatch/ArtifactDataset
# if by_count is True, each count is weighted separately for balanced loss within that count
def calculate_batch_weights(batch, dataset, by_count: bool):
    # TODO: we need a parameter to control the relative weight of unlabeled loss to labeled loss
    # For batch index n, we want weight[n] = dataset.weights[alt_counts[n], labels[n], variant_types[n]]
    counts = batch.get_alt_counts()
    labels = batch.get_labels()
    variant_types = batch.get_variant_types()

    return utils.index_3d_array(dataset.weights, counts, labels, variant_types) if by_count else \
        utils.index_2d_array(dataset.weights[ALL_COUNTS_SENTINEL], labels, variant_types)


# note: this works for both BaseBatch/BaseDataset AND ArtifactBatch/ArtifactDataset
# if by_count is True, each count is weighted separately for balanced loss within that count
def calculate_batch_source_weights(batch, dataset, by_count: bool):
    # For batch index n, we want weight[n] = dataset.source_weights[alt_counts[n], sources[n], variant_types[n]]
    counts = batch.get_alt_counts()
    sources = batch.get_sources()
    variant_types = batch.get_variant_types()

    return utils.index_3d_array(dataset.source_weights, counts, sources, variant_types) if by_count else \
        utils.index_2d_array(dataset.source_weights[ALL_COUNTS_SENTINEL], sources, variant_types)


class LearningMethod(Enum):
    # train the embedding by minimizing cross-entropy loss of binary predictor on labeled data
    SUPERVISED = "SUPERVISED"

    # same but use entropy regularization loss on unlabeled data
    SEMISUPERVISED = "SEMISUPERVISED"

    # TODO: IMPLEMENT THIS
    # optimize a clustering model with center triplet loss
    SUPERVISED_CLUSTERING = "SUPERVISED_CLUSTERING"

    # TODO: IMPLEMENT THIS
    # modify data via a finite set of affine transformations and train the embedding to recognize which was applied
    AFFINE_TRANSFORMATION = "AFFINE"

    # modify data via a finite set of affine transformations and train the embedding to recognize which was applied
    MASK_PREDICTION = "MASK_PREDICTION"

    AUTOENCODER = "AUTOENCODER"

    DEEPSAD = "DEEPSAD"

    MARS = "MARS"


def make_gated_ref_alt_mlp_encoder(input_dimension: int, params: BaseModelParameters):
    return GatedRefAltMLP(d_model=input_dimension, d_ffn=params.self_attention_hidden_dimension, num_blocks=params.num_self_attention_layers)


class BaseModel(torch.nn.Module):
    """
    DeepSets framework for reads and variant info.  We embed each read and concatenate the mean ref read
    embedding, mean alt read embedding, and variant info embedding, then apply an aggregation function to
    this concatenation to obtain an embedding / representation of the read set for downstream use such as
    variant filtering and clustering.

    hidden_read_layers: dimensions of layers for embedding reads, excluding input dimension, which is the
    size of each read's 1D tensor

    hidden_info_layers: dimensions of layers for embedding variant info, excluding input dimension, which is the
    size of variant info 1D tensor

    aggregation_layers: dimensions of layers for aggregation, excluding its input which is determined by the
    read and info embeddings.

    output_layers: dimensions of layers after aggregation, excluding the output dimension,
    which is 1 for a single logit representing artifact/non-artifact.  This is not part of the aggregation layers
    because we have different output layers for each variant type.
    """

    def __init__(self, params: BaseModelParameters, num_read_features: int, num_info_features: int, ref_sequence_length: int, device=utils.gpu_if_available()):
        super(BaseModel, self).__init__()

        self._device = device
        self._dtype = DEFAULT_GPU_FLOAT if device != torch.device("cpu") else DEFAULT_CPU_FLOAT
        self._ref_sequence_length = ref_sequence_length
        self._params = params

        # embeddings of reads, info, and reference sequence prior to the transformer layers
        self.read_embedding = MLP([num_read_features] + params.read_layers, batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)
        self.info_embedding = MLP([num_info_features] + params.info_layers, batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)
        self.ref_seq_cnn = DNASequenceConvolution(params.ref_seq_layer_strings, ref_sequence_length)

        embedding_dim = self.read_embedding.output_dimension() + self.info_embedding.output_dimension() + self.ref_seq_cnn.output_dimension()

        self.ref_alt_reads_encoder = make_gated_ref_alt_mlp_encoder(embedding_dim, params)

        # after encoding alt reads (along with info and ref seq embeddings and with self-attention to ref reads)
        # pass through another MLP
        self.aggregation = MLP([embedding_dim] + params.aggregation_layers, batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)

        self.to(device=self._device, dtype=self._dtype)

    def output_dimension(self) -> int:
        return self.aggregation.output_dimension()

    def ref_alt_seq_embedding_dimension(self) -> int:
        return self.ref_seq_cnn.output_dimension()

    def ref_sequence_length(self) -> int:
        return self._ref_sequence_length

    def set_epoch_type(self, epoch_type: utils.Epoch):
        if epoch_type == utils.Epoch.TRAIN:
            self.train(True)
            utils.unfreeze(self.parameters())
        else:
            self.train(False)
            utils.freeze(self.parameters())

    # I really don't like the forward method of torch.nn.Module with its implicit calling that PyCharm doesn't recognize
    def forward(self, batch: BaseBatch):
        pass

    # here 'v' means "variant index within a batch", 'r' means "read index within a variant or the batch", 'e' means "index within an embedding"
    # so, for example, "re" means a 2D tensor with all reads in the batch stacked and "vre" means a 3D tensor indexed
    # first by variant within the batch, then the read
    def calculate_representations(self, batch: BaseBatch, weight_range: float = 0) -> torch.Tensor:
        ref_counts, alt_counts = batch.get_ref_counts(), batch.get_alt_counts()
        total_ref, total_alt = torch.sum(ref_counts).item(), torch.sum(alt_counts).item()

        read_embeddings_re = self.read_embedding.forward(batch.get_reads_2d().to(dtype=self._dtype))
        info_embeddings_ve = self.info_embedding.forward(batch.get_info_2d().to(dtype=self._dtype))
        ref_seq_embeddings_ve = self.ref_seq_cnn(batch.get_ref_sequences_2d().to(dtype=self._dtype))
        info_and_seq_ve = torch.hstack((info_embeddings_ve, ref_seq_embeddings_ve))
        info_and_seq_re = torch.vstack((torch.repeat_interleave(info_and_seq_ve, repeats=ref_counts, dim=0),
                                       torch.repeat_interleave(info_and_seq_ve, repeats=alt_counts, dim=0)))
        reads_info_seq_re = torch.hstack((read_embeddings_re, info_and_seq_re))

        # TODO: might be a bug if every datum in batch has zero ref reads?
        ref_reads_info_seq_re = reads_info_seq_re[:total_ref]
        alt_reads_info_seq_re = reads_info_seq_re[total_ref:]

        # TODO: make sure it handles ref count = 0 case
        transformed_ref_re, transformed_alt_re = self.ref_alt_reads_encoder.forward(ref_reads_info_seq_re, alt_reads_info_seq_re, ref_counts, alt_counts)

        alt_weights_r = 1 + weight_range * (1 - 2 * torch.rand(total_alt, device=self._device, dtype=self._dtype))

        # normalize so read weights within each variant sum to 1
        alt_wt_sums_v = utils.sums_over_rows(alt_weights_r, alt_counts)
        normalized_alt_weights_r = alt_weights_r / torch.repeat_interleave(alt_wt_sums_v, repeats=alt_counts, dim=0)

        alt_means_ve = utils.sums_over_rows(transformed_alt_re * normalized_alt_weights_r[:,None], alt_counts)

        result_ve = self.aggregation.forward(alt_means_ve)

        return result_ve, ref_seq_embeddings_ve # ref seq embeddings are useful later

    def make_dict_for_saving(self, prefix: str = ""):
        return {(prefix + constants.STATE_DICT_NAME): self.state_dict(),
                (prefix + constants.HYPERPARAMS_NAME): self._params,
                (prefix + constants.NUM_READ_FEATURES_NAME): self.read_embedding.input_dimension(),
                (prefix + constants.NUM_INFO_FEATURES_NAME): self.info_embedding.input_dimension(),
                (prefix + constants.REF_SEQUENCE_LENGTH_NAME): self.ref_sequence_length()}

    def save(self, path):
        torch.save(self.make_dict_for_saving(), path)


def base_model_from_saved_dict(saved, prefix: str = "", device: torch.device = utils.gpu_if_available()):
    hyperparams = saved[prefix + constants.HYPERPARAMS_NAME]
    num_read_features = saved[prefix + constants.NUM_READ_FEATURES_NAME]
    num_info_features = saved[prefix + constants.NUM_INFO_FEATURES_NAME]
    ref_sequence_length = saved[prefix + constants.REF_SEQUENCE_LENGTH_NAME]

    model = BaseModel(hyperparams, num_read_features=num_read_features, num_info_features=num_info_features,
                      ref_sequence_length=ref_sequence_length, device=device)
    model.load_state_dict(saved[prefix + constants.STATE_DICT_NAME])

    # in case the state dict had the wrong dtype for the device we're on now eg base model was pretrained on GPU
    # and we're now on CPU
    model.to(model._dtype)

    return model


def load_base_model(path, prefix: str = "", device: torch.device = utils.gpu_if_available()) -> BaseModel:
    saved = torch.load(path, map_location=device)
    return base_model_from_saved_dict(saved, prefix, device)


# outputs a 1D tensor of losses over the batch.  We assume it needs the representations of the batch data from the base
# model.  We nonetheless also use the model as an input because there are some learning strategies that involve
# computing representations of a modified batch.
class BaseModelLearningStrategy(ABC):
    @abstractmethod
    def loss_function(self, base_model: BaseModel, base_batch: BaseBatch, base_model_representations: torch.Tensor):
        pass


class BaseModelSemiSupervisedLoss(torch.nn.Module, BaseModelLearningStrategy):
    def __init__(self, input_dim: int, hidden_top_layers: List[int], params: BaseModelParameters):
        super(BaseModelSemiSupervisedLoss, self).__init__()

        self.bce = torch.nn.BCEWithLogitsLoss(reduction='none')  # no reduction because we may want to first multiply by weights for unbalanced data

        # go from base model output representation to artifact logit for supervised loss
        self.logit_predictor = MLP([input_dim] + hidden_top_layers + [1], batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)

    def loss_function(self, base_model: BaseModel, base_batch: BaseBatch, base_model_representations: torch.Tensor):
        logits = self.logit_predictor.forward(base_model_representations).reshape((base_batch.size()))
        labels = base_batch.get_training_labels()

        # base batch always has labels, but for unlabeled elements these labels are meaningless and is_labeled_mask is zero
        cross_entropies = self.bce(logits, labels)
        probabilities = torch.sigmoid(logits)
        entropies = self.bce(logits, probabilities)

        return base_batch.get_is_labeled_mask() * cross_entropies + (1 - base_batch.get_is_labeled_mask()) * entropies

    # I don't like implicit forward!!
    def forward(self):
        pass


def permute_columns_independently(mat: torch.Tensor):
    assert mat.dim() == 2
    num_rows, num_cols = mat.size()
    weights = torch.ones(num_rows)

    result = torch.clone(mat)
    for col in range(num_cols):
        idx = torch.multinomial(weights, num_rows, replacement=True)
        result[:, col] = result[:, col][idx]
    return result


# randomly choose read features to "mask" -- where a masked feature is permuted randomly over all the reads in the batch.
# this essentially means drawing masked features from the empirical marginal distribution
# the pretext self-supervision task is, for each datum, to predict which features were masked
# note that this basically means destroy correlations for a random selection of features
class BaseModelMaskPredictionLoss(torch.nn.Module, BaseModelLearningStrategy):
    def __init__(self, num_read_features: int, base_model_output_dim: int, hidden_top_layers: List[int], params: BaseModelParameters):
        super(BaseModelMaskPredictionLoss, self).__init__()

        self.num_read_features = num_read_features

        self.bce = torch.nn.BCEWithLogitsLoss(reduction='none')  # no reduction because we may want to first multiply by weights for unbalanced data

        # go from base model output representation to artifact logit for supervised loss
        self.mask_predictor = MLP([base_model_output_dim] + hidden_top_layers + [num_read_features], batch_normalize=params.batch_normalize, dropout_p=params.dropout_p)

    def loss_function(self, base_model: BaseModel, base_batch: BaseBatch, base_model_representations):
        # TODO: this is broken now that batches have mixed counts
        '''ref_count, alt_count = base_batch.ref_count, base_batch.alt_count
        total_ref, total_alt = ref_count * base_batch.size(), alt_count * base_batch.size()

        alt_reads_2d = base_batch.get_reads_2d()[total_ref:]
        permuted_reads = permute_columns_independently(base_batch.get_reads_2d())
        permuted_alt_reads = permuted_reads[:total_alt]

        datum_mask = torch.bernoulli(0.1 * torch.ones(base_batch.size(), self.num_read_features))

        # each read within a datum gets the same mask
        reads_mask = torch.repeat_interleave(datum_mask, repeats=alt_count, dim=0)

        original_reads_2d = base_batch.reads_2d
        modified_alt_reads = alt_reads_2d * (1 - reads_mask) + permuted_alt_reads * reads_mask
        base_batch.reads_2d = torch.vstack((original_reads_2d[:total_ref], modified_alt_reads))

        # TODO: is there any reason to fix the batch with base_batch.reads_2d = original_reads_2d?

        # shape is batch size x num read features, each entry being a logit for "was this feature masked in this datum?"
        mask_prediction_logits = self.mask_predictor.forward(base_model_representations)

        # by batch index and feature
        losses_bf = self.bce(mask_prediction_logits, datum_mask)
        return torch.mean(losses_bf, dim=1)   # average over read features
        '''
        pass

    # I don't like implicit forward!!
    def forward(self):
        pass


# chamfer distance between two 3D tensors B x N1 x E and B x N2 x E, where B is the batch size, N1/2 are the number
# of items in the two sets, and E is the dimensionality of each item
# returns a 1D tensor of length B
def chamfer_distance(set1_bne, set2_bne):
    diffs_bnne = torch.unsqueeze(set1_bne, dim=2) - torch.unsqueeze(set2_bne, dim=1)
    l1_dists_bnn = torch.mean(torch.abs(diffs_bnne), dim=-1)

    chamfer_dists12_bn = torch.min(l1_dists_bnn, dim=-2).values
    chamfer_dists21_bn = torch.min(l1_dists_bnn, dim=-1).values
    symmetric_chamfer_b = torch.mean(chamfer_dists12_bn, dim=-1) + torch.mean(chamfer_dists21_bn, dim=-1)
    return symmetric_chamfer_b


# self-supervision approach where we use the base model embedding to regenerate the set and use Chamfer distance as the
# reconstruction error.  We regenerate the set via the Transformer Set Prediction Network approach of Kosiorek et al -- seed a set
# of N reads by concatenated the embedding with N random vectors, then map it so the final reconstructed set with transformers.
class BaseModelAutoencoderLoss(torch.nn.Module, BaseModelLearningStrategy):
    def __init__(self, read_dim: int, hidden_top_layers: List[int], params: BaseModelParameters):
        super(BaseModelAutoencoderLoss, self).__init__()
        self.base_model_output_dimension = params.output_dimension()

        # TODO: explore making random seed dimension different from the base model embedding dimension
        self.random_seed_dimension = self.base_model_output_dimension
        self.transformer_dimension = self.base_model_output_dimension + self.random_seed_dimension

        # TODO: maybe also a parameter to scale the random vectors?

        # TODO: should these decoder params be the same as the base model encoder params?  It seems reasonable.

        # TODO: this is broken -- use the ref_alt_encoder
        #self.alt_decoder = make_gated_mlp_encoder(self.transformer_dimension, params)
        #self.ref_decoder = make_gated_mlp_encoder(self.transformer_dimension, params)

        self.mapping_back_to_reads = MLP([self.transformer_dimension] + hidden_top_layers + [read_dim])

    def loss_function(self, base_model: BaseModel, base_batch: BaseBatch, base_model_representations):
        # TODO: this is broken now that batches have mixed counts
        '''var_count, alt_count, ref_count = base_batch.size(), base_batch.alt_count, base_batch.ref_count

        total_ref, total_alt = ref_count * var_count, alt_count * var_count

        representations_ve = base_model_representations
        random_alt_seeds_vre = torch.randn(var_count, alt_count, self.random_seed_dimension)
        random_ref_seeds_vre = torch.randn(var_count, ref_count, self.random_seed_dimension) if ref_count > 0 else None
        alt_representations_vre = torch.unsqueeze(representations_ve, dim=1).expand(-1, alt_count, -1) # repeat over the dummy read index
        ref_representations_vre = torch.unsqueeze(representations_ve, dim=1).expand(-1, ref_count, -1)

        alt_vre = torch.cat((alt_representations_vre, random_alt_seeds_vre), dim=-1)
        ref_vre = torch.cat((ref_representations_vre, random_ref_seeds_vre), dim=-1) if ref_count > 0 else None

        # TODO: update these to reflect mixed-count batches.  Gated MLPs now take inputs flattened over batch dimension
        # TODO: and have an extra input of ref and alt read counts
        decoded_alt_vre = self.alt_decoder.forward(alt_vre)
        decoded_ref_vre = self.ref_decoder.forward(ref_vre) if ref_count > 0 else None

        decoded_alt_re = torch.reshape(decoded_alt_vre, (var_count * alt_count, -1))
        decoded_ref_re = torch.reshape(decoded_ref_vre, (var_count * ref_count, -1)) if ref_count > 0 else None

        # the raw read tensors are quantile normalized with Gaussian output
        reconstructed_alt_vre = torch.reshape(self.mapping_back_to_reads(decoded_alt_re),(var_count, alt_count, -1))
        reconstructed_ref_vre = torch.reshape(self.mapping_back_to_reads(decoded_ref_re), (var_count, ref_count, -1)) if ref_count > 0 else None

        original_alt_vre = base_batch.get_reads_2d()[total_ref:].reshape(var_count, alt_count, -1)
        original_ref_vre = base_batch.get_reads_2d()[:total_ref].reshape(var_count, ref_count, -1) if ref_count > 0 else None

        alt_chamfer_dist = chamfer_distance(original_alt_vre, reconstructed_alt_vre)
        ref_chamfer_dist = chamfer_distance(original_ref_vre, reconstructed_ref_vre) if ref_count > 0 else 0
        return alt_chamfer_dist + ref_chamfer_dist
        '''
        pass

    # I don't like implicit forward!!
    def forward(self):
        pass


class BaseModelDeepSADLoss(torch.nn.Module, BaseModelLearningStrategy):
    def __init__(self, embedding_dim: int):
        super(BaseModelDeepSADLoss, self).__init__()
        self.embedding_dim = embedding_dim
        self.normal_centroid = Parameter(torch.zeros(embedding_dim))

    # normal embeddings should cluster near the origin and artifact embeddings should be far
    def loss_function(self, base_model: BaseModel, base_batch: BaseBatch, base_model_representations):
        dist_squared = torch.square(torch.norm(base_model_representations - self.normal_centroid, dim=1))

        # labels are 1 for artifact, 0 otherwise.  We convert to +1 if normal, -1 if artifact
        # DeepSAD assumes most unlabeled data are normal and so the unlabeled loss is identical to the normal loss, that is,
        # squared Euclidean distance from the centroid
        signs = (1 - 2 * base_batch.get_training_labels()) * base_batch.get_is_labeled_mask() + 1 * (1 - base_batch.get_is_labeled_mask())

        # distance squared for normal and unlabeled, inverse distance squared for artifact
        return dist_squared ** signs

    # I don't like implicit forward!!
    def forward(self):
        pass


class BaseModelMARSLoss(torch.nn.Module, BaseModelLearningStrategy):
    def __init__(self, embedding_dim: int):
        super(BaseModelMARSLoss, self).__init__()
        self.embedding_dim = embedding_dim
        # TODO: magic constants!!!!!
        self.num_normal_clusters = 3
        self.num_artifact_clusters = 5

        # weight of centroid-centroid loss vs embedding-centroid loss
        self.tau = 0.2

        # ce denotes indexing by cluster, then embedding dimension
        self.centroids_ce = Parameter(torch.zeros((self.num_normal_clusters + self.num_artifact_clusters, embedding_dim)))

    # normal embeddings should cluster near the origin and artifact embeddings should be far
    def loss_function(self, base_model: BaseModel, base_batch: BaseBatch, base_model_representations):
        embeddings_be = base_model_representations
        seps_bce = torch.unsqueeze(embeddings_be, dim=1) - torch.unsqueeze(self.centroids_ce, dim=0)
        dist_squared_bc = torch.square(torch.norm(seps_bce, dim=-1))

        normal_dist_squared_b = torch.min(dist_squared_bc[:, :self.num_normal_clusters], dim=-1).values
        artifact_dist_squared_b = torch.min(dist_squared_bc[:, self.num_normal_clusters:], dim=-1).values
        min_dist_squared_b = torch.min(dist_squared_bc, dim=-1).values

        # closest centroid with correct label is labeled, otherwise just the closest centroid
        labeled_losses_b = (base_batch.get_training_labels() * artifact_dist_squared_b + (1 - base_batch.get_training_labels()) * normal_dist_squared_b)
        unlabeled_losses_b = min_dist_squared_b
        embedding_centroid_losses_b =  base_batch.get_is_labeled_mask() * labeled_losses_b + (1 - base_batch.get_is_labeled_mask()) * unlabeled_losses_b

        # average distance between centroids
        centroid_seps_cce = torch.unsqueeze(self.centroids_ce, dim=0) - torch.unsqueeze(self.centroids_ce, dim=1)
        centroid_dist_squared_cc = torch.square(torch.norm(centroid_seps_cce, dim=-1))
        centroid_centroid_loss = torch.mean(centroid_dist_squared_cc)

        # TODO: need to control arbitrarily negative loss achieved by making an unused centroid go far away from the others
        # note that broadcasting the subtraction means the centroid-centroid loss is repeated once for each datum in the batch
        return embedding_centroid_losses_b - self.tau * centroid_centroid_loss

    # I don't like implicit forward!!
    def forward(self):
        pass


# artifact model parameters are for simultaneously training an artifact model on top of the base model
# to measure quality, especially in unsupervised training when the loss metric isn't directly related to accuracy or cross-entropy
def learn_base_model(base_model: BaseModel, dataset: BaseDataset, learning_method: LearningMethod, training_params: TrainingParameters,
                     summary_writer: SummaryWriter, validation_fold: int = None):
    print(f"Memory usage percent: {psutil.virtual_memory().percent:.1f}")
    is_cuda = base_model._device.type == 'cuda'
    print(f"Is CUDA available? {is_cuda}")

    for idx, variation_type in enumerate(utils.Variation):
        print(f"For variation type {variation_type.name}, there are {int(dataset.totals[ALL_COUNTS_SENTINEL][Label.ARTIFACT][idx].item())} \
            artifacts, {int(dataset.totals[ALL_COUNTS_SENTINEL][Label.VARIANT][idx].item())} \
            non-artifacts, and {int(dataset.totals[ALL_COUNTS_SENTINEL][Label.UNLABELED][idx].item())} unlabeled data.")

    # TODO: use Python's match syntax, but this requires updating Python version in the docker
    # TODO: hidden_top_layers are hard-coded!
    if learning_method == LearningMethod.SUPERVISED or learning_method == LearningMethod.SEMISUPERVISED:
        learning_strategy = BaseModelSemiSupervisedLoss(input_dim=base_model.output_dimension(), hidden_top_layers=[30,-1,-1,-1,10], params=base_model._params)
    elif learning_method == LearningMethod.MASK_PREDICTION:
        learning_strategy = BaseModelMaskPredictionLoss(num_read_features=dataset.num_read_features,
                                                        base_model_output_dim=base_model.output_dimension(), hidden_top_layers=[10,10,10], params=base_model._params)
    elif learning_method == LearningMethod.AUTOENCODER:
        learning_strategy = BaseModelAutoencoderLoss(read_dim=dataset.num_read_features, hidden_top_layers=[20,20,20], params=base_model._params)
    elif learning_method == LearningMethod.DEEPSAD:
        learning_strategy = BaseModelDeepSADLoss(embedding_dim=base_model.output_dimension())
    elif learning_method == LearningMethod.MARS:
        learning_strategy = BaseModelMARSLoss(embedding_dim=base_model.output_dimension())
    else:
        raise Exception("not implemented yet")
    learning_strategy.to(device=base_model._device, dtype=base_model._dtype)

    # adversarial loss to learn features that forget the alt count
    alt_count_gradient_reversal = GradientReversal(alpha=0.01)  #initialize as barely active
    alt_count_predictor = MLP([base_model.output_dimension()] + [30, -1, -1, -1, 1]).to(device=base_model._device, dtype=base_model._dtype)
    alt_count_loss_func = torch.nn.MSELoss(reduction='none')
    alt_count_adversarial_metrics = LossMetrics()

    # TODO: fused = is_cuda?
    train_optimizer = torch.optim.AdamW(chain(base_model.parameters(), learning_strategy.parameters(), alt_count_predictor.parameters()),
                                        lr=training_params.learning_rate, weight_decay=training_params.weight_decay)
    # train scheduler needs to be given the thing that's supposed to decrease at the end of each epoch
    train_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        train_optimizer, factor=0.2, patience=5, threshold=0.001, min_lr=(training_params.learning_rate/100), verbose=True)

    classifier_on_top = MLP([base_model.output_dimension()] + [30, -1, -1, -1, 10] + [1])\
        .to(device=base_model._device, dtype=base_model._dtype)
    classifier_bce = torch.nn.BCEWithLogitsLoss(reduction='none')

    classifier_optimizer = torch.optim.AdamW(classifier_on_top.parameters(),
                                             lr=training_params.learning_rate,
                                             weight_decay=training_params.weight_decay,
                                             fused=is_cuda)
    classifier_metrics = LossMetrics()

    validation_fold_to_use = (dataset.num_folds - 1) if validation_fold is None else validation_fold
    train_loader = dataset.make_data_loader(dataset.all_but_one_fold(validation_fold_to_use), training_params.batch_size, is_cuda, training_params.num_workers)
    valid_loader = dataset.make_data_loader([validation_fold_to_use], training_params.batch_size, is_cuda, training_params.num_workers)

    for epoch in trange(1, training_params.num_epochs + 1, desc="Epoch"):
        p = epoch - 1
        new_alpha = (2/(1 + math.exp(-0.1*p))) - 1
        alt_count_gradient_reversal.set_alpha(new_alpha) # alpha increases linearly
        start_epoch = time.time()
        print(f"Start of epoch {epoch}, memory usage percent: {psutil.virtual_memory().percent:.1f}")
        for epoch_type in (utils.Epoch.TRAIN, utils.Epoch.VALID):
            base_model.set_epoch_type(epoch_type)
            loss_metrics = LossMetrics()

            loader = train_loader if epoch_type == utils.Epoch.TRAIN else valid_loader
            loader_iter = iter(loader)

            next_batch_cpu = next(loader_iter)
            next_batch = next_batch_cpu.copy_to(base_model._device, non_blocking=is_cuda)

            pbar = tqdm(range(len(loader)), mininterval=60)
            for n in pbar:
                batch_cpu = next_batch_cpu
                batch = next_batch

                # Optimization: Asynchronously send the next batch to the device while the model does work
                next_batch_cpu = next(loader_iter)
                next_batch = next_batch_cpu.copy_to(base_model._device, non_blocking=is_cuda)

                # TODO: we need a parameter to control the relative weight of unlabeled loss to labeled loss
                weights = calculate_batch_weights(batch_cpu, dataset, by_count=True)
                weights = weights.to(device=base_model._device, dtype=base_model._dtype, non_blocking=True)

                # unused output is the embedding of ref and alt alleles with context
                representations, _ = base_model.calculate_representations(batch, weight_range=base_model._params.reweighting_range)
                losses = learning_strategy.loss_function(base_model, batch, representations)

                if losses is None:
                    continue

                loss_metrics.record_losses(losses.detach(), batch, weights)

                # gradient reversal means parameters before the representation try to maximize alt count prediction loss, i.e. features
                # try to forget alt count, while parameters after the representation try to minimize it, i.e. they try
                # to achieve the adversarial task
                alt_count_pred = torch.sigmoid(alt_count_predictor.forward(alt_count_gradient_reversal(representations)).squeeze())
                alt_count_target = batch.get_alt_counts().to(dtype=alt_count_pred.dtype)/20
                alt_count_losses = alt_count_loss_func(alt_count_pred, alt_count_target)

                alt_count_adversarial_metrics.record_losses(alt_count_losses.detach(), batch, weights=torch.ones_like(alt_count_losses))

                loss = torch.sum((weights * losses) + alt_count_losses)

                classification_logits = classifier_on_top.forward(representations.detach()).reshape(batch.size())
                classification_losses = classifier_bce(classification_logits, batch.get_training_labels())
                classification_loss = torch.sum(batch.get_is_labeled_mask() * weights * classification_losses)
                classifier_metrics.record_losses(classification_losses.detach(), batch, batch.get_is_labeled_mask() * weights)

                if epoch_type == utils.Epoch.TRAIN:
                    utils.backpropagate(train_optimizer, loss)
                    utils.backpropagate(classifier_optimizer, classification_loss)

            # done with one epoch type -- training or validation -- for this epoch
            loss_metrics.write_to_summary_writer(epoch_type, epoch, summary_writer)
            classifier_metrics.write_to_summary_writer(epoch_type, epoch, summary_writer, prefix="auxiliary-classifier-")
            alt_count_adversarial_metrics.write_to_summary_writer(epoch_type, epoch, summary_writer, prefix="alt-count-adversarial-predictor")

            if epoch_type == utils.Epoch.TRAIN:
                train_scheduler.step(loss_metrics.get_labeled_loss())

            print(f"Labeled base model loss for {epoch_type.name} epoch {epoch}: {loss_metrics.get_labeled_loss():.3f}")
            print(f"Labeled auxiliary classifier loss for {epoch_type.name} epoch {epoch}: {classifier_metrics.get_labeled_loss():.3f}")
            print(f"Alt count adversarial loss for {epoch_type.name} epoch {epoch}: {alt_count_adversarial_metrics.get_labeled_loss():.3f}")
        print(f"End of epoch {epoch}, memory usage percent: {psutil.virtual_memory().percent:.1f}, time elapsed(s): {time.time() - start_epoch:.2f}")
        # done with training and validation for this epoch
        # note that we have not learned the AF spectrum yet
    # done with training

    record_embeddings(base_model, train_loader, summary_writer)


# after training for visualizing clustering etc of base model embeddings
def record_embeddings(base_model: BaseModel, loader, summary_writer: SummaryWriter):
    # base_model.freeze_all() whoops -- it doesn't have freeze_all
    embedding_metrics = EmbeddingMetrics()
    ref_alt_seq_metrics = EmbeddingMetrics()

    pbar = tqdm(enumerate(loader), mininterval=60)
    for n, batch_cpu in pbar:
        batch = batch_cpu.copy_to(base_model._device, non_blocking=base_model._device.type=='cuda')
        representations, ref_alt_seq_embeddings = base_model.calculate_representations(batch, weight_range=base_model._params.reweighting_range)

        representations = representations.cpu()
        ref_alt_seq_embeddings = ref_alt_seq_embeddings.cpu()

        labels = [("artifact" if label > 0.5 else "non-artifact") if is_labeled > 0.5 else "unlabeled" for (label, is_labeled) in
                  zip(batch.get_training_labels().tolist(), batch.get_is_labeled_mask().tolist())]
        for (metrics, embeddings) in [(embedding_metrics, representations), (ref_alt_seq_metrics, ref_alt_seq_embeddings)]:
            metrics.label_metadata.extend(labels)
            metrics.correct_metadata.extend(["unknown"] * batch.size())
            metrics.type_metadata.extend([Variation(idx).name for idx in batch.get_variant_types().tolist()])
            alt_count_strings = [str(round_up_to_nearest_three(min(MAX_COUNT, ac))) for ac in batch.get_alt_counts().tolist()]
            metrics.truncated_count_metadata.extend(alt_count_strings)
            metrics.representations.append(embeddings)
    embedding_metrics.output_to_summary_writer(summary_writer)
    ref_alt_seq_metrics.output_to_summary_writer(summary_writer, prefix="ref and alt allele context")

