from argparse import Namespace
import tempfile

from permutect.tools import preprocess_dataset, train_model, filter_variants
from permutect import constants


def test_on_dream1():
    # Input Files
    training_datasets = ["/Users/davidben/permutect/just-dream-1/dream1-normal-small-training.dataset"]
    #mutect2_vcf = "/Users/davidben/permutect/dream-vcfs/dream1-50000.vcf"
    mutect2_vcf = "/Users/davidben/permutect/integration-test/dream1-mutect2-small.vcf"
    #filtering_dataset = "/Users/davidben/permutect/just-dream-1/dream1-test.dataset"
    filtering_dataset = "/Users/davidben/permutect/integration-test/dream1-small-test.dataset"

    # Intermediate and Output Files
    training_data_tarfile = tempfile.NamedTemporaryFile()
    saved_artifact_model = tempfile.NamedTemporaryFile()
    training_tensorboard_dir = tempfile.TemporaryDirectory()
    filtering_tensorboard_dir = tempfile.TemporaryDirectory()
    filtered_mutect3_vcf = tempfile.NamedTemporaryFile()

    # STEP 1: preprocess the plain text training dataset yielding a training tarfile
    preprocess_args = Namespace()
    setattr(preprocess_args, constants.CHUNK_SIZE_NAME, 1e6)
    setattr(preprocess_args, constants.TRAINING_DATASETS_NAME, training_datasets)
    setattr(preprocess_args, constants.OUTPUT_NAME, training_data_tarfile.name)
    setattr(preprocess_args, constants.SOURCES_NAME, [0])
    preprocess_dataset.main_without_parsing(preprocess_args)

    # STEP 2: train a model
    train_model_args = Namespace()
    setattr(train_model_args, constants.READ_LAYERS_NAME, [10, 10, 10])
    setattr(train_model_args, constants.SELF_ATTENTION_HIDDEN_DIMENSION_NAME, 20)
    setattr(train_model_args, constants.NUM_SELF_ATTENTION_LAYERS_NAME, 3)
    setattr(train_model_args, constants.INFO_LAYERS_NAME, [30, 30, 30])
    setattr(train_model_args, constants.AGGREGATION_LAYERS_NAME, [30, 30, 30, 30])
    setattr(train_model_args, constants.CALIBRATION_LAYERS_NAME, [6,6])
    cnn_layer_strings = ['convolution/kernel_size=3/out_channels=64',
                     'pool/kernel_size=2',
                     'leaky_relu',
                     'convolution/kernel_size=3/dilation=2/out_channels=5',
                     'leaky_relu',
                     'flatten',
                     'linear/out_features=10']
    setattr(train_model_args, constants.REF_SEQ_LAYER_STRINGS_NAME, cnn_layer_strings)
    setattr(train_model_args, constants.DROPOUT_P_NAME, 0.0)
    setattr(train_model_args, constants.LEARNING_RATE_NAME, 0.001)
    setattr(train_model_args, constants.WEIGHT_DECAY_NAME, 0.01)
    setattr(train_model_args, constants.BATCH_NORMALIZE_NAME, False)
    setattr(train_model_args, constants.LEARN_ARTIFACT_SPECTRA_NAME, True)  # could go either way
    setattr(train_model_args, constants.GENOMIC_SPAN_NAME, 100000)

    # Training data inputs
    setattr(train_model_args, constants.TRAIN_TAR_NAME, training_data_tarfile.name)

    # training hyperparameters
    setattr(train_model_args, constants.REWEIGHTING_RANGE_NAME, 0.3)
    setattr(train_model_args, constants.BATCH_SIZE_NAME, 64)
    setattr(train_model_args, constants.NUM_WORKERS_NAME, 2)
    setattr(train_model_args, constants.NUM_EPOCHS_NAME, 2)
    setattr(train_model_args, constants.NUM_CALIBRATION_EPOCHS_NAME, 1)
    setattr(train_model_args, constants.NUM_REFLESS_EPOCHS_NAME, 2)

    # path to saved model
    setattr(train_model_args, constants.OUTPUT_NAME, saved_artifact_model.name)
    setattr(train_model_args, constants.TENSORBOARD_DIR_NAME, training_tensorboard_dir.name)

    train_model.main_without_parsing(train_model_args)

    # STEP 3: call variants
    filtering_args = Namespace()
    setattr(filtering_args, constants.INPUT_NAME, mutect2_vcf)
    setattr(filtering_args, constants.TEST_DATASET_NAME, filtering_dataset)
    setattr(filtering_args, constants.M3_MODEL_NAME, saved_artifact_model.name)
    setattr(filtering_args, constants.OUTPUT_NAME, filtered_mutect3_vcf.name)
    setattr(filtering_args, constants.TENSORBOARD_DIR_NAME, filtering_tensorboard_dir.name)
    setattr(filtering_args, constants.BATCH_SIZE_NAME, 64)
    setattr(filtering_args, constants.CHUNK_SIZE_NAME, 100000)
    setattr(filtering_args, constants.NUM_SPECTRUM_ITERATIONS_NAME, 10)
    setattr(filtering_args, constants.INITIAL_LOG_VARIANT_PRIOR_NAME, -10.0)
    setattr(filtering_args, constants.INITIAL_LOG_ARTIFACT_PRIOR_NAME, -10.0)
    setattr(filtering_args, constants.GENOMIC_SPAN_NAME, 100000)
    setattr(filtering_args, constants.MAF_SEGMENTS_NAME, None)
    setattr(filtering_args, constants.NORMAL_MAF_SEGMENTS_NAME, None)
    setattr(filtering_args, constants.GERMLINE_MODE_NAME, False)
    setattr(filtering_args, constants.NO_GERMLINE_MODE_NAME, False)

    filter_variants.main_without_parsing(filtering_args)
    h = 9
