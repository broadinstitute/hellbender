from argparse import Namespace
import tempfile

from permutect.data import base_datum
from permutect.data.base_dataset import BaseDataset
from permutect.tools import preprocess_dataset
from permutect import constants, utils
from permutect.utils import extract_to_temp_dir


def test_on_10_megabases_singular():
    training_datasets = ["/Users/davidben/mutect3/permutect/integration-tests/singular-10-Mb/training-dataset.txt"]
    training_data_tarfile = tempfile.NamedTemporaryFile()

    preprocess_args = Namespace()
    setattr(preprocess_args, constants.CHUNK_SIZE_NAME, 1e6)
    setattr(preprocess_args, constants.TRAINING_DATASETS_NAME, training_datasets)
    setattr(preprocess_args, constants.OUTPUT_NAME, training_data_tarfile.name)
    setattr(preprocess_args, constants.SOURCES_NAME, [0])
    preprocess_dataset.main_without_parsing(preprocess_args)

    with tempfile.TemporaryDirectory() as train_temp_dir:
        training_files = extract_to_temp_dir(training_data_tarfile.name, train_temp_dir)
        for training_file in training_files:
            base_data_list = base_datum.load_list_of_base_data(training_file)

    dataset = BaseDataset(data_tarfile=training_data_tarfile.name, num_folds=10)