from permutect.test.test_utils import artificial_data
from permutect.data.base_dataset import BaseDataset, make_test_data_loader
from permutect.data.base_datum import BaseDatum
from typing import Iterable
from permutect.architecture.artifact_model import ArtifactModel
from permutect.parameters import ArtifactModelParameters
from permutect import utils
from permutect.tools.train_model import TrainingParameters

from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
import tempfile
from torch.utils.tensorboard import SummaryWriter

BATCH_SIZE = 64
CHUNK_SIZE = 100000
NUM_EPOCHS = 50
NUM_CALIBRATION_EPOCHS=10
NUM_SPECTRUM_ITERATIONS = 100
TRAINING_PARAMS = TrainingParameters(batch_size=BATCH_SIZE, num_epochs=NUM_EPOCHS, num_calibration_epochs=NUM_CALIBRATION_EPOCHS, reweighting_range=0.3)

REF_SEQ_LAYER_STRINGS = ['convolution/kernel_size=3/out_channels=64',
                     'pool/kernel_size=3',
                     'leaky_relu',
                     'convolution/kernel_size=3/dilation=1/out_channels=5',
                     'leaky_relu',
                     'flatten',
                     'linear/out_features=10']
SMALL_MODEL_PARAMS = ArtifactModelParameters(read_embedding_dimension=12,
    num_transformer_heads=3, transformer_hidden_dimension=10, num_transformer_layers=2,
    info_layers=[5, 5], aggregation_layers=[5, 5, 5, 5], calibration_layers=[6],
    ref_seq_layers_strings=REF_SEQ_LAYER_STRINGS,
    dropout_p=0.2, batch_normalize=False, learning_rate=0.001, weight_decay=0.01, alt_downsample=20)


# Note that the test methods in this class also cover batching, samplers, datasets, and data loaders
def train_model_and_write_summary(hyperparams: ArtifactModelParameters, training_params: TrainingParameters,
                                  data: Iterable[BaseDatum], summary_writer: SummaryWriter = None):
    dataset = BaseDataset(data=data)
    big_dataset = BigReadSetDataset(batch_size=training_params.batch_size, dataset=dataset, num_workers=2)
    model = ArtifactModel(params=hyperparams, num_read_features=dataset.num_read_features(), num_info_features=dataset.num_info_features(), ref_sequence_length=dataset.ref_sequence_length()).float()

    model.learn(big_dataset, training_params.num_epochs, training_params.num_calibration_epochs, summary_writer=summary_writer,
                reweighting_range=training_params.reweighting_range, hyperparams=hyperparams)
    model.evaluate_model_after_training({"training": big_dataset.generate_batches(utils.Epoch.TRAIN)}, summary_writer)
    return model


def test_big_data():
    training_dataset_file = "/Users/davidben/mutect3/just-dream-1/dream1-normal-medium-training.dataset"
    big_dataset = BigReadSetDataset(batch_size=64, max_bytes_per_chunk=int(100*1e6), dataset_files=[training_dataset_file], num_workers=2)
    params = SMALL_MODEL_PARAMS
    training_params = TRAINING_PARAMS

    with tempfile.TemporaryDirectory() as tensorboard_dir:
        summary_writer = SummaryWriter(tensorboard_dir)
        model = ArtifactModel(params=params, num_read_features=big_dataset.num_read_features, num_info_features=big_dataset.num_info_features, ref_sequence_length=big_dataset.ref_sequence_length).float()
        model.learn(big_dataset, training_params.num_epochs, training_params.num_calibration_epochs, summary_writer=summary_writer,
                    reweighting_range=training_params.reweighting_range, hyperparams=params)
        model.evaluate_model_after_training({"training": big_dataset.generate_batches(utils.Epoch.TRAIN)}, summary_writer)

        events = EventAccumulator(tensorboard_dir)
        events.Reload()


def test_separate_gaussian_data():
    # in the test for alt count agnostic, we make training data where variant alt counts are much larger than artifact
    # alt counts and test data with a low alt allele fraction
    for test_alt_fraction_agnostic in (False, True):
        data = artificial_data.make_two_gaussian_data(1000) if not test_alt_fraction_agnostic else \
            artificial_data.make_two_gaussian_data(1000, vaf=0.5, downsample_variants_to_match_artifacts=False, alt_downsampling=20)
        params = SMALL_MODEL_PARAMS
        training_params = TRAINING_PARAMS

        with tempfile.TemporaryDirectory() as tensorboard_dir:
            summary_writer = SummaryWriter(tensorboard_dir)
            model = train_model_and_write_summary(hyperparams=params, training_params=training_params, data=data, summary_writer=summary_writer)

            # TODO: migrate this old stuff to test for PosteriorModel
            # test_vaf = 0.05 if test_alt_fraction_agnostic else 0.5
            # test_data = artificial_data.make_two_gaussian_data(1000, is_training_data=False, vaf=test_vaf, unlabeled_fraction=0.0)
            # test_dataset = ReadSetDataset(data=test_data)
            # test_loader = make_test_data_loader(test_dataset, BATCH_SIZE)
            # model.learn_spectra(test_loader, NUM_SPECTRUM_ITERATIONS, summary_writer=summary_writer)

            events = EventAccumulator(tensorboard_dir)
            events.Reload()

            # TODO: these have been replaced with images, so it's not so simple to check the output quality from the tensorboard
            # TODO: for now I can put in a breakpoint and manually run tensorboard --logdir <tensorboard temp directory>
            # TODO: to spot check the figures
            # assert events.Scalars('Variant Sensitivity')[0].value > 0.98
            # assert events.Scalars('Artifact Sensitivity')[0].value > 0.98


def test_wide_and_narrow_gaussian_data():
    data = artificial_data.make_wide_and_narrow_gaussian_data(10000)
    params = SMALL_MODEL_PARAMS
    training_params = TRAINING_PARAMS

    with tempfile.TemporaryDirectory() as tensorboard_dir:
        summary_writer = SummaryWriter(tensorboard_dir)
        model = train_model_and_write_summary(hyperparams=params, training_params=training_params, data=data, summary_writer=summary_writer)

        events = EventAccumulator(tensorboard_dir)
        events.Reload()


# TODO: this test currently fails -- almost everything is considered an artifact
# TODO: I must investigate
def test_strand_bias_data():
    data = artificial_data.make_random_strand_bias_data(1000, is_training_data=True)
    params = SMALL_MODEL_PARAMS # TODO: change!!!!!!!
    training_params = TRAINING_PARAMS

    with tempfile.TemporaryDirectory() as tensorboard_dir:
        summary_writer = SummaryWriter(tensorboard_dir)
        model = train_model_and_write_summary(hyperparams=params, training_params=training_params, data=data, summary_writer=summary_writer)

        test_data = artificial_data.make_random_strand_bias_data(1000, is_training_data=False, vaf=0.25, unlabeled_fraction=0.0)
        test_dataset = BaseDataset(data=test_data)
        test_loader = make_test_data_loader(test_dataset, BATCH_SIZE)
        model.learn_spectra(test_loader, NUM_SPECTRUM_ITERATIONS, summary_writer=summary_writer)

        events = EventAccumulator(tensorboard_dir)
        events.Reload()

        assert events.Scalars('Variant Sensitivity')[0].value > 0.90
        assert events.Scalars('Artifact Sensitivity')[0].value > 0.90
