import tempfile
from argparse import Namespace

from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
from permutect import constants, utils
from permutect.tools import train_model
from permutect.architecture.artifact_model import load_artifact_model


def test_train_model():
    # Inputs
    training_data_tarfile = '/Users/davidben/mutect3/permutect/integration-tests/singular-10-Mb/preprocessed-dataset.tar'
    base_model = '/Users/davidben/mutect3/permutect/integration-tests/singular-10-Mb/base-model.pt'

    # Outputs
    saved_artifact_model = tempfile.NamedTemporaryFile()
    training_tensorboard_dir = tempfile.TemporaryDirectory()

    # STEP 2: train a model
    train_model_args = Namespace()
    setattr(train_model_args, constants.AGGREGATION_LAYERS_NAME, [20, 20, 20])
    setattr(train_model_args, constants.CALIBRATION_LAYERS_NAME, [6,6])
    setattr(train_model_args, constants.DROPOUT_P_NAME, 0.0)
    setattr(train_model_args, constants.BATCH_NORMALIZE_NAME, False)
    setattr(train_model_args, constants.LEARN_ARTIFACT_SPECTRA_NAME, True)  # could go either way
    setattr(train_model_args, constants.GENOMIC_SPAN_NAME, 100000)

    # Training data inputs
    setattr(train_model_args, constants.TRAIN_TAR_NAME, training_data_tarfile)
    setattr(train_model_args, constants.BASE_MODEL_NAME, base_model)

    # training hyperparameters
    setattr(train_model_args, constants.BATCH_SIZE_NAME, 64)
    setattr(train_model_args, constants.INFERENCE_BATCH_SIZE_NAME, 64)
    setattr(train_model_args, constants.NUM_WORKERS_NAME, 0)
    setattr(train_model_args, constants.NUM_EPOCHS_NAME, 2)
    setattr(train_model_args, constants.NUM_CALIBRATION_EPOCHS_NAME, 1)
    setattr(train_model_args, constants.LEARNING_RATE_NAME, 0.001)
    setattr(train_model_args, constants.WEIGHT_DECAY_NAME, 0.01)

    # path to saved model
    setattr(train_model_args, constants.OUTPUT_NAME, saved_artifact_model.name)
    setattr(train_model_args, constants.TENSORBOARD_DIR_NAME, training_tensorboard_dir.name)

    train_model.main_without_parsing(train_model_args)

    events = EventAccumulator(training_tensorboard_dir.name)
    events.Reload()

    device = utils.gpu_if_available()
    loaded_artifact_model, artifact_log_priors, artifact_spectra_state_dict = load_artifact_model(saved_artifact_model, device=device)
    assert artifact_log_priors is not None
    assert artifact_spectra_state_dict is not None

    print(artifact_log_priors)
    h = 99

