import tempfile
from argparse import Namespace

from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
from permutect import constants
from permutect.data.base_dataset import BaseDataset
from permutect.tools import prune_dataset


def test_prune_dataset():
    # Inputs
    training_data_tarfile = '/Users/davidben/mutect3/permutect/integration-tests/singular-10-Mb/preprocessed-dataset.tar'
    base_model = '/Users/davidben/mutect3/permutect/integration-tests/singular-10-Mb/base-model.pt'

    # Outputs
    pruned_dataset = tempfile.NamedTemporaryFile()
    training_tensorboard_dir = tempfile.TemporaryDirectory()

    # STEP 2: train a model
    prune_dataset_args = Namespace()
    setattr(prune_dataset_args, constants.AGGREGATION_LAYERS_NAME, [20, 20, 20])
    setattr(prune_dataset_args, constants.CALIBRATION_LAYERS_NAME, [6,6])
    setattr(prune_dataset_args, constants.DROPOUT_P_NAME, 0.0)
    setattr(prune_dataset_args, constants.BATCH_NORMALIZE_NAME, False)
    setattr(prune_dataset_args, constants.LEARN_ARTIFACT_SPECTRA_NAME, True)  # could go either way
    setattr(prune_dataset_args, constants.GENOMIC_SPAN_NAME, 100000)

    # Training data inputs
    setattr(prune_dataset_args, constants.TRAIN_TAR_NAME, training_data_tarfile)
    setattr(prune_dataset_args, constants.BASE_MODEL_NAME, base_model)

    setattr(prune_dataset_args, constants.CHUNK_SIZE_NAME, 2e9)

    # training hyperparameters
    setattr(prune_dataset_args, constants.BATCH_SIZE_NAME, 64)
    setattr(prune_dataset_args, constants.NUM_WORKERS_NAME, 2)
    setattr(prune_dataset_args, constants.NUM_EPOCHS_NAME, 2)
    setattr(prune_dataset_args, constants.NUM_CALIBRATION_EPOCHS_NAME, 1)
    setattr(prune_dataset_args, constants.LEARNING_RATE_NAME, 0.001)
    setattr(prune_dataset_args, constants.WEIGHT_DECAY_NAME, 0.01)

    # path to saved model
    setattr(prune_dataset_args, constants.OUTPUT_NAME, pruned_dataset.name)
    setattr(prune_dataset_args, constants.TENSORBOARD_DIR_NAME, training_tensorboard_dir.name)

    prune_dataset.main_without_parsing(prune_dataset_args)

    events = EventAccumulator(training_tensorboard_dir.name)
    events.Reload()

    pruned_base_dataset = BaseDataset(data_tarfile=pruned_dataset, num_folds=10)
    h = 99

