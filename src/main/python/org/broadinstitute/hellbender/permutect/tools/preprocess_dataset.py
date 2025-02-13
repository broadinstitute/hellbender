import argparse
import os
import tarfile
import tempfile
from typing import List

from permutect import constants
from permutect.data import base_datum
from permutect.data.plain_text_data import generate_normalized_data
from permutect.utils import ConsistentValue

"""
This tool takes as input a list of text file Mutect3 training datasets, reads them in chunks that fit in memory,
normalizes each chunk, outputs each chunk as a binary PyTorch file, and bundles the output as a tarfile.
"""


def parse_arguments():
    parser = argparse.ArgumentParser(description='preprocess plain text training dataset into tarfile of nprmalized binary data')
    parser.add_argument('--' + constants.TRAINING_DATASETS_NAME, nargs='+', type=str, required=True,
                        help='list of plain text data files')
    parser.add_argument('--' + constants.CHUNK_SIZE_NAME, type=int, default=int(2e9), required=False,
                        help='size in bytes of output binary data files')
    parser.add_argument('--' + constants.SOURCES_NAME, nargs='+', type=int, required=False,
                        help='integer sources corresponding to plain text data files for distinguishing different sequencing conditions')
    parser.add_argument('--' + constants.OUTPUT_NAME, type=str, default=None, required=True,
                        help='path to output tarfile of training data')
    return parser.parse_args()


def do_work(training_datasets, training_output_file, chunk_size, sources: List[int]):
    data_files = []
    num_read_features, num_info_features, ref_sequence_length = ConsistentValue(), ConsistentValue(), ConsistentValue()

    # save all the lists of read sets to tempfiles. . .
    # TODO: left off here.  Need to give it sources, which will need to be command line argument
    for base_data_list in generate_normalized_data(training_datasets, max_bytes_per_chunk=chunk_size, sources=sources):
        num_read_features.check(base_data_list[0].get_reads_2d().shape[1])
        num_info_features.check(base_data_list[0].get_info_tensor_1d().shape[0])
        ref_sequence_length.check(base_data_list[0].get_ref_sequence_1d().shape[0])

        with tempfile.NamedTemporaryFile(delete=False) as train_data_file:
            base_datum.save_list_base_data(base_data_list, train_data_file)
            data_files.append(train_data_file.name)

    # . . . and bundle them in a tarfile
    with tarfile.open(training_output_file, "w") as train_tar:
        for train_file in data_files:
            train_tar.add(train_file, arcname=os.path.basename(train_file))


def main_without_parsing(args):
    chunk_size = getattr(args, constants.CHUNK_SIZE_NAME)
    training_datasets = getattr(args, constants.TRAINING_DATASETS_NAME)
    output_file = getattr(args, constants.OUTPUT_NAME)
    sources = getattr(args, constants.SOURCES_NAME)

    do_work(training_datasets, output_file, chunk_size, sources)


def main():
    args = parse_arguments()
    main_without_parsing(args)


if __name__ == '__main__':
    main()
