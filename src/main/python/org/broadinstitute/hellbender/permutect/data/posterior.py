import copy
import random
import math
from typing import List, Iterable

import torch
from torch import IntTensor
from torch.utils.data import Dataset, DataLoader
from permutect.data.base_datum import Variant, CountsAndSeqLks, bases5_as_base_string

from permutect import utils
from permutect.utils import Label, Variation


def variant_from_int_array(subarray) -> Variant:
    contig = subarray[0].item()
    position = subarray[1].item()
    ref = bases5_as_base_string(subarray[2].item())  # ref and alt are the base-5 encoding as integers
    alt = bases5_as_base_string(subarray[3].item())
    return Variant(contig, position, ref, alt)


class PosteriorDatum:
    CONTIG = 0
    POSITION = 1
    REF = 2
    ALT = 3
    VAR_TYPE = 4
    DEPTH = 5
    ALT_COUNT = 6
    NORMAL_DEPTH = 7
    NORMAL_ALT_COUNT = 8
    LABEL = 9

    SEQ_ERROR_LOG_LK = 0
    TLOD_FROM_M2 = 1
    NORMAL_SEQ_ERROR_LOG_LK = 2
    ALLELE_FREQUENCY = 3
    ARTIFACT_LOGIT = 4
    MAF = 5
    NORMAL_MAF = 6

    def __init__(self, variant: Variant, counts_and_seq_lks: CountsAndSeqLks, allele_frequency: float,
                 artifact_logit: float, embedding: torch.Tensor, label: Label, maf: float, normal_maf: float):
        self.embedding = embedding

        this_class = self.__class__
        self.int_array = torch.zeros(10, dtype=int)
        self.int_array[this_class.CONTIG] = variant.contig
        self.int_array[this_class.POSITION] = variant.position
        self.int_array[this_class.REF] = variant.get_ref_as_int()    # ref and alt are the base-5 encoding as integers
        self.int_array[this_class.ALT] = variant.get_alt_as_int()
        self.int_array[this_class.VAR_TYPE] = utils.Variation.get_type(variant.ref, variant.alt)  # Variation is IntEnum so this is int
        self.int_array[this_class.DEPTH] = counts_and_seq_lks.depth
        self.int_array[this_class.ALT_COUNT] = counts_and_seq_lks.alt_count
        self.int_array[this_class.NORMAL_DEPTH] = counts_and_seq_lks.normal_depth
        self.int_array[this_class.NORMAL_ALT_COUNT] = counts_and_seq_lks.normal_alt_count
        self.int_array[this_class.LABEL] = label

        self.float_array = torch.zeros(7, dtype=torch.float16)
        self.float_array[this_class.SEQ_ERROR_LOG_LK] = counts_and_seq_lks.seq_error_log_lk
        self.float_array[this_class.TLOD_FROM_M2] = -counts_and_seq_lks.seq_error_log_lk - math.log(counts_and_seq_lks.depth + 1)
        self.float_array[this_class.NORMAL_SEQ_ERROR_LOG_LK] = counts_and_seq_lks.normal_seq_error_log_lk
        self.float_array[this_class.ALLELE_FREQUENCY] = allele_frequency
        self.float_array[this_class.ARTIFACT_LOGIT] = artifact_logit
        self.float_array[this_class.MAF] = maf
        self.float_array[this_class.NORMAL_MAF] = normal_maf

    def get_variant(self) -> Variant:
        this_class = self.__class__
        subarray = self.int_array[this_class.CONTIG:this_class.ALT + 1]
        return variant_from_int_array(subarray)

    def get_artifact_logit(self) -> float:
        return self.float_array[self.__class__.ARTIFACT_LOGIT]


class PosteriorBatch:

    def __init__(self, data: List[PosteriorDatum]):
        self.embeddings = torch.vstack([item.embedding for item in data]).float()
        self.int_tensor = torch.vstack([item.int_array for item in data])
        self.float_tensor = torch.vstack([item.float_array for item in data]).float()

        self._size = len(data)

    def pin_memory(self):
        self.embeddings = self.embeddings.pin_memory()
        self.int_tensor = self.int_tensor.pin_memory()
        self.float_tensor = self.float_tensor.pin_memory()
        return self

    # dtype is just for floats!!! Better not convert the int tensor to a float accidentally!
    def copy_to(self, device, dtype, non_blocking):
        # For all non-tensor attributes, shallow copy is sufficient
        new_batch = copy.copy(self)

        new_batch.embeddings = self.embeddings.to(device=device, dtype=dtype, non_blocking=non_blocking)
        new_batch.int_tensor = self.int_tensor.to(device=device, non_blocking=non_blocking)
        new_batch.float_tensor = self.float_tensor.to(device=device, dtype=dtype, non_blocking=non_blocking)

        return new_batch

    def get_variants(self) -> List[Variant]:
        subarray_2d = self.int_tensor[:, PosteriorDatum.CONTIG:PosteriorDatum.ALT + 1]
        return [variant_from_int_array(subarray) for subarray in subarray_2d]

    def get_variant_types(self) -> torch.IntTensor:
        return self.int_tensor[:, PosteriorDatum.VAR_TYPE]

    def get_alt_counts(self) -> torch.Tensor:
        return self.int_tensor[:, PosteriorDatum.ALT_COUNT]

    def get_depths(self) -> torch.Tensor:
        return self.int_tensor[:, PosteriorDatum.DEPTH]

    def get_labels(self) -> torch.Tensor:
        return self.int_tensor[:, PosteriorDatum.LABEL]

    def get_normal_alt_counts(self) -> torch.Tensor:
        return self.int_tensor[:, PosteriorDatum.NORMAL_ALT_COUNT]

    def get_normal_depths(self) -> torch.Tensor:
        return self.int_tensor[:, PosteriorDatum.NORMAL_DEPTH]

    def get_tlods_from_m2(self) -> torch.Tensor:
        return self.float_tensor[:, PosteriorDatum.TLOD_FROM_M2]

    def get_allele_frequencies(self) -> torch.Tensor:
        return self.float_tensor[:, PosteriorDatum.ALLELE_FREQUENCY]

    def get_artifact_logits(self) -> torch.Tensor:
        return self.float_tensor[:, PosteriorDatum.ARTIFACT_LOGIT]

    def get_mafs(self) -> torch.Tensor:
        return self.float_tensor[:, PosteriorDatum.MAF]

    def get_normal_mafs(self) -> torch.Tensor:
        return self.float_tensor[:, PosteriorDatum.NORMAL_MAF]

    def size(self) -> int:
        return self._size

    def get_normal_ref_counts(self) -> IntTensor:
        return self.get_normal_depths() - self.get_normal_alt_counts()


class PosteriorDataset(Dataset):
    def __init__(self, data: Iterable[PosteriorDatum], shuffle: bool = True):
        self.data = data

        if shuffle:
            random.shuffle(self.data)

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, index) -> PosteriorDatum:
        return self.data[index]

    def make_data_loader(self, batch_size: int, pin_memory: bool = False, num_workers: int = 0):
        return DataLoader(dataset=self, batch_size=batch_size, pin_memory=pin_memory, num_workers=num_workers, collate_fn=PosteriorBatch)
