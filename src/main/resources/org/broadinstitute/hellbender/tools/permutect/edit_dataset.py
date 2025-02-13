import argparse
import os
import tarfile
import tempfile
from enum import Enum

import psutil
import torch.utils.data

from permutect.data import base_datum
from tqdm.autonotebook import tqdm

from permutect import constants
from permutect.data.base_dataset import BaseDataset
from permutect.utils import Label


class EditType(Enum):
    UNLABEL_ARTIFACTS = "unlabel_artifacts"
    UNLABEL_VARIANTS = "unlabel_variants"
    UNLABEL_EVERYTHING = "unlabel_everything"
    REMOVE_ARTIFACTS = "remove_artifacts"
    REMOVE_VARIANTS = "remove_variants"
    KEEP_EVERYTHING = "keep_everything"


# generates BaseDatum(s) from the original dataset that *pass* the pruning thresholds
def generate_edited_data(base_datasets, edit_type: str, source: int):
    pbar = tqdm(enumerate(torch.utils.data.ConcatDataset(base_datasets)), mininterval=60)

    for n, base_datum in pbar:
        if source is not None:
            base_datum.set_source(source)

        if edit_type == EditType.UNLABEL_ARTIFACTS.value:
            if base_datum.label == Label.ARTIFACT:
                base_datum.set_label(Label.UNLABELED)
            yield base_datum
        elif edit_type == EditType.UNLABEL_VARIANTS.value:
            if base_datum.label == Label.VARIANT:
                base_datum.set_label(Label.UNLABELED)
            yield base_datum
        elif edit_type == EditType.UNLABEL_EVERYTHING.value:
            base_datum.set_label(Label.UNLABELED)
            yield base_datum
        elif edit_type == EditType.REMOVE_ARTIFACTS.value:
            if base_datum.label != Label.ARTIFACT:
                yield base_datum
        elif edit_type == EditType.REMOVE_VARIANTS.value:
            if base_datum.label != Label.VARIANT:
                yield base_datum
        elif edit_type == EditType.KEEP_EVERYTHING.value:
            yield base_datum
        else:
            raise Exception(f"edit type {edit_type} not implemented yet")


# takes a ReadSet generator and organizes into buffers.
# TODO: probably code duplication since the generator is already pruned
def generate_output_data_buffers(output_data_generator, max_bytes_per_chunk: int):
    buffer, bytes_in_buffer = [], 0
    for datum in output_data_generator:

        buffer.append(datum)
        bytes_in_buffer += datum.size_in_bytes()
        if bytes_in_buffer > max_bytes_per_chunk:
            print(f"Memory usage percent: {psutil.virtual_memory().percent:.1f}")
            print(f"{bytes_in_buffer} bytes in chunk")
            yield buffer
            buffer, bytes_in_buffer = [], 0

    # There will be some data left over, in general.
    if buffer:
        yield buffer


def make_output_training_dataset(pruned_data_buffer_generator, output_tarfile):
    pruned_data_files = []
    for base_data_list in pruned_data_buffer_generator:
        with tempfile.NamedTemporaryFile(delete=False) as train_data_file:
            base_datum.save_list_base_data(base_data_list, train_data_file)
            pruned_data_files.append(train_data_file.name)

    # bundle them in a tarfile
    with tarfile.open(output_tarfile, "w") as train_tar:
        for train_file in pruned_data_files:
            train_tar.add(train_file, arcname=os.path.basename(train_file))


def parse_arguments():
    parser = argparse.ArgumentParser(description='train the Mutect3 artifact model')
    parser.add_argument('--' + constants.CHUNK_SIZE_NAME, type=int, default=int(2e9), required=False,
                        help='size in bytes of output binary data files')
    parser.add_argument('--' + constants.DATASET_EDIT_TYPE_NAME, type=str, required=True,
                        help='how to modify the dataset')
    parser.add_argument('--' + constants.SOURCE_NAME, type=int, required=False, help='new source integer to apply')

    # input / output
    parser.add_argument('--' + constants.TRAIN_TAR_NAME, nargs='+', type=str, required=True,
                        help='tarfile(s) of training/validation datasets produced by preprocess_dataset.py')
    parser.add_argument('--' + constants.OUTPUT_NAME, type=str, required=True, help='path to pruned dataset file')

    return parser.parse_args()


def main_without_parsing(args):
    original_tarfiles = getattr(args, constants.TRAIN_TAR_NAME) # list of files
    output_tarfile = getattr(args, constants.OUTPUT_NAME)
    chunk_size = getattr(args, constants.CHUNK_SIZE_NAME)
    edit_type = getattr(args, constants.DATASET_EDIT_TYPE_NAME)
    new_source = getattr(args, constants.SOURCE_NAME)
    base_datasets = map(lambda original_tarfile: BaseDataset(data_tarfile=original_tarfile), original_tarfiles)

    # generate ReadSets
    output_data_generator = generate_edited_data(base_datasets, edit_type, new_source)

    # generate List[ReadSet]s
    output_data_buffer_generator = generate_output_data_buffers(output_data_generator, chunk_size)

    make_output_training_dataset(output_data_buffer_generator, output_tarfile=output_tarfile)


def main():
    args = parse_arguments()
    main_without_parsing(args)


if __name__ == '__main__':
    main()
