"""
Functions for reading from plain text dataset files of the format

UNLABELED                                       # label
1:13118,A->G                                    # locus and mutation
GAGGAAAGTGAGGTTGCCTGC                           # reference context
0.12 0.09 2.00 0.57 1.00 1.00 1.00 1.00 1.00    # variant-level info vector
5 4 0 0                                         # ref count, alt count, matched normal ref count, matched normal alt count
27 30 0 1 11 29 333 321 12 0 0                  # one ref read vector per line
50 22 1 0 30 19 342 70 272 0 0
27 31 0 1 32 17 236 203 33 0 0
27 20 0 1 32 17 141 72 69 1 0
21 28 1 0 49 0 232 49 183 1 0
23 29 1 1 40 9 335 294 41 0 0                   # one alt read vector per line
24 29 0 1 38 11 354 315 39 0 0
24 30 0 1 36 13 351 314 37 0 0
23 30 1 1 42 7 341 298 43 0 0
51 13 0 0                                       # original ref, alt, normal ref, normal alt counts before downsampling
-108.131                                        # sequencing error log likelihood
-0.000                                          # matched normal sequencing error log likelihood
VARIANT
1:13302,C->T
GTCCTGGACACGCTGTTGGCC
0.00 0.00 0.00 1.00 1.00 2.00 1.00 1.00 1.00
2 1 0 0
24 29 0 0 11 21 338 11 327 0 0
50 25 0 1 49 -8 355 303 52 0 0
23 33 1 0 13 21 312 87 225 0 0
69 4 0 0
-11.327
-0.000
"""
from typing import List

import numpy as np
import psutil
from sklearn.preprocessing import QuantileTransformer

from permutect import utils
from permutect.data.base_datum import BaseDatum, Variant, CountsAndSeqLks, DEFAULT_NUMPY_FLOAT

from permutect.utils import Label, Variation

MAX_VALUE = 10000
EPSILON = 0.00001
QUANTILE_DATA_COUNT = 10000


def read_data(dataset_file, only_artifacts: bool = False, source: int=0):
    """
    generator that yields data from a plain text dataset file.
    """
    with open(dataset_file) as file:
        n = 0
        while label_str := file.readline().strip():
            label = Label.get_label(label_str)
            passes_label_filter = (label == Label.ARTIFACT or not only_artifacts)
            n += 1

            # contig:position,ref->alt
            variant_line = file.readline().strip()
            locus, mutation = variant_line.split(",")
            contig, position = map(int, locus.split(":"))   # contig is an integer *index* from a sequence dictionary
            # TODO: replace with tqdm progress bar by counting file in initial pass.  It can't be that expensive.
            if n % 100000 == 0:
                print(f"{contig}:{position}")
            ref_allele, alt_allele = mutation.strip().split("->")

            ref_sequence_string = file.readline().strip()
            gatk_info_tensor = line_to_tensor(file.readline())
            ref_tensor_size, alt_tensor_size, normal_ref_tensor_size, normal_alt_tensor_size = map(int, file.readline().strip().split())

            # the first column is read group index, which we currently discard
            # later we're going to want to use this
            ref_tensor = read_2d_tensor(file, ref_tensor_size)[:,1:] if ref_tensor_size > 0 else None
            alt_tensor = read_2d_tensor(file, alt_tensor_size)[:,1:] if alt_tensor_size > 0 else None

            # normal_ref_tensor = read_2d_tensor(file, normal_ref_tensor_size)  # not currently used
            # normal_alt_tensor = read_2d_tensor(file, normal_alt_tensor_size)  # not currently used
            # round down normal tensors as well

            depth, alt_count, normal_depth, normal_alt_count = read_integers(file.readline())
            seq_error_log_lk = read_float(file.readline())
            normal_seq_error_log_lk = read_float(file.readline())

            if alt_tensor_size > 0 and passes_label_filter:
                variant = Variant(contig, position, ref_allele, alt_allele)
                counts_and_seq_lks = CountsAndSeqLks(depth, alt_count, normal_depth, normal_alt_count, seq_error_log_lk, normal_seq_error_log_lk)
                yield BaseDatum.from_gatk(ref_sequence_string, Variation.get_type(ref_allele, alt_allele), ref_tensor,
                                          alt_tensor, gatk_info_tensor, label, source, variant, counts_and_seq_lks)


# if sources is None, source is set to zero
# if List is length-1, that's the source for all files
# otherwise each file has its own source int
def generate_normalized_data(dataset_files, max_bytes_per_chunk: int, sources: List[int]=None):
    """
    given text dataset files, generate normalized lists of read sets that fit in memory

    In addition to quantile-normalizing read tensors it also enlarges the info tensors
    :param dataset_files:
    :param max_bytes_per_chunk:
    :return:
    """
    for n, dataset_file in enumerate(dataset_files):
        buffer, bytes_in_buffer = [], 0
        read_quantile_transform = QuantileTransformer(n_quantiles=100, output_distribution='normal')
        info_quantile_transform = QuantileTransformer(n_quantiles=100, output_distribution='normal')

        num_buffers_filled = 0
        source = 0 if sources is None else (sources[0] if len(sources) == 1 else sources[n])
        for base_datum in read_data(dataset_file, source=source):
            buffer.append(base_datum)
            bytes_in_buffer += base_datum.size_in_bytes()
            if bytes_in_buffer > max_bytes_per_chunk:
                print(f"Memory usage percent: {psutil.virtual_memory().percent:.1f}")
                print(f"{bytes_in_buffer} bytes in chunk")

                normalize_buffer(buffer, read_quantile_transform, info_quantile_transform)
                yield buffer
                num_buffers_filled += 1
                buffer, bytes_in_buffer = [], 0
        # There will be some data left over, in general.  Since it's small, use the last buffer's
        # quantile transforms for better statistical power if it's from the same text file
        if buffer:
            normalize_buffer(buffer, read_quantile_transform, info_quantile_transform, refit_transforms=(num_buffers_filled==0))
            yield buffer


# this normalizes the buffer and also prepends new features to the info tensor
def normalize_buffer(buffer, read_quantile_transform, info_quantile_transform, refit_transforms=True):
    EPSILON = 0.00001   # tiny quantity for jitter

    # 2D array.  Rows are ref/alt reads, columns are read features
    all_ref = np.vstack([datum.get_ref_reads_2d() for datum in buffer])
    all_reads = np.vstack([datum.reads_2d for datum in buffer])

    # 2D array.  Rows are read sets, columns are info features
    all_info = np.vstack([datum.get_info_tensor_1d() for datum in buffer])

    all_ref_jittered = all_ref + EPSILON * np.random.randn(*all_ref.shape)
    all_reads_jittered = all_reads + EPSILON * np.random.randn(*all_reads.shape)
    all_info_jittered = all_info + EPSILON * np.random.randn(*all_info.shape)

    num_read_features = all_ref.shape[1]
    binary_read_columns = binary_column_indices(all_ref)    # make sure not to use jittered arrays here!

    # 1 if is binary, 0 if not binary
    binary_read_column_mask = np.zeros(num_read_features)
    binary_read_column_mask[binary_read_columns] = 1

    binary_info_columns = binary_column_indices(all_info)   # make sure not to use jittered arrays here!

    if refit_transforms:    # fit quantiles column by column (aka feature by feature)
        read_quantile_transform.fit(all_ref_jittered)
        info_quantile_transform.fit(all_info_jittered)

    # it's more efficient to apply the quantile transform to all reads at once, then split it back into read sets
    all_reads_transformed = transform_except_for_binary_columns(all_reads_jittered, read_quantile_transform, binary_read_columns)
    all_info_transformed = transform_except_for_binary_columns(all_info_jittered, info_quantile_transform, binary_info_columns)

    read_counts = np.array([len(datum.reads_2d) for datum in buffer])
    read_index_ranges = np.cumsum(read_counts)

    for n, datum in enumerate(buffer):
        datum.reads_2d = all_reads_transformed[0 if n == 0 else read_index_ranges[n-1]:read_index_ranges[n]]

        # medians are an appropriate outlier-tolerant summary, except for binary columns where the mean makes more sense
        alt_medians = np.median(datum.get_alt_reads_2d(), axis=0)
        alt_means = np.mean(datum.get_alt_reads_2d(), axis=0)

        extra_info = binary_read_column_mask * alt_means + (1 - binary_read_column_mask) * alt_medians
        datum.set_info_tensor_1d(np.hstack([extra_info, all_info_transformed[n]]))


def line_to_tensor(line: str) -> np.ndarray:
    tokens = line.strip().split()
    floats = [float(token) for token in tokens]
    return np.clip(np.array(floats, dtype=DEFAULT_NUMPY_FLOAT), -MAX_VALUE, MAX_VALUE)


def read_2d_tensor(file, num_lines: int) -> np.ndarray:
    if num_lines == 0:
        return None
    lines = [file.readline() for _ in range(num_lines)]
    return np.vstack([line_to_tensor(line) for line in lines])


def read_integers(line: str):
    return map(int, line.strip().split())


def read_float(line: str):
    return float(line.strip().split()[0])


def is_binary(column_tensor_1d: np.ndarray):
    assert len(column_tensor_1d.shape) == 1
    return all(el.item() == 0 or el.item() == 1 for el in column_tensor_1d)


def binary_column_indices(tensor_2d: np.ndarray):
    assert len(tensor_2d.shape) == 2
    return [n for n in range(tensor_2d.shape[1]) if is_binary(tensor_2d[:, n])]


def non_binary_column_indices(tensor_2d: np.ndarray):
    assert len(tensor_2d.shape) == 2
    return [n for n in range(tensor_2d.shape[1]) if not is_binary(tensor_2d[:, n])]


# copy the unnormalized values of binary features (columns)
# we modify the normalized values in-place
def restore_binary_columns(normalized, original, binary_columns):
    result = normalized
    if len(normalized.shape) == 2:
        result[:, binary_columns] = original[:, binary_columns]
    elif len(normalized.shape) == 1:
        result[binary_columns] = original[binary_columns]
    else:
        raise Exception("This is only for 1D or 2D tensors")
    return result


def transform_except_for_binary_columns(tensor_2d, quantile_transform: QuantileTransformer, binary_column_indices):
    return restore_binary_columns(quantile_transform.transform(tensor_2d), tensor_2d, binary_column_indices)