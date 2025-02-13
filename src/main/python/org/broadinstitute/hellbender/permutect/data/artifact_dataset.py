import math
import random
from typing import List

import torch
from tqdm.autonotebook import tqdm
from torch.utils.data import Dataset, DataLoader, Sampler

from permutect.architecture.base_model import BaseModel
from permutect.data.base_datum import ArtifactDatum, ArtifactBatch
from permutect.data.base_dataset import BaseDataset, chunk


# given a ReadSetDataset, apply a BaseModel to get an ArtifactDataset (in RAM, maybe implement memory map later)
# of RepresentationReadSets
class ArtifactDataset(Dataset):
    def __init__(self, base_dataset: BaseDataset,
                 base_model: BaseModel,
                 folds_to_use: List[int] = None,
                 base_loader_num_workers=0,
                 base_loader_batch_size=8192):
        self.counts_by_source = base_dataset.counts_by_source
        self.totals = base_dataset.totals
        self.source_totals = base_dataset.source_totals
        self.weights = base_dataset.weights
        self.source_weights = base_dataset.source_weights

        self.artifact_data = []
        self.num_folds = base_dataset.num_folds
        self.labeled_indices = [[] for _ in range(self.num_folds)]  # one list for each fold
        self.unlabeled_indices = [[] for _ in range(self.num_folds)]    # ditto
        self.num_base_features = base_model.output_dimension()
        self.num_ref_alt_features = base_model.ref_alt_seq_embedding_dimension()

        index = 0

        loader = base_dataset.make_data_loader(base_dataset.all_folds() if folds_to_use is None else folds_to_use,
                                               batch_size=base_loader_batch_size,
                                               num_workers=base_loader_num_workers)
        print("making artifact dataset from base dataset")

        is_cuda = base_model._device.type == 'cuda'
        print(f"Is base model using CUDA? {is_cuda}")

        pbar = tqdm(enumerate(loader), mininterval=60)
        for n, base_batch_cpu in pbar:
            base_batch = base_batch_cpu.copy_to(base_model._device, non_blocking=is_cuda)
            with torch.inference_mode():
                representations, _ = base_model.calculate_representations(base_batch)

            for representation, base_datum in zip(representations.detach().cpu(), base_batch_cpu.original_list()):
                artifact_datum = ArtifactDatum(base_datum, representation.detach())
                self.artifact_data.append(artifact_datum)
                fold = index % self.num_folds
                if artifact_datum.is_labeled():
                    self.labeled_indices[fold].append(index)
                else:
                    self.unlabeled_indices[fold].append(index)
                index += 1

    def __len__(self):
        return len(self.artifact_data)

    def __getitem__(self, index):
        return self.artifact_data[index]

    # it is often convenient to arbitrarily use the last fold for validation
    def last_fold_only(self):
        return [self.num_folds - 1]  # use the last fold for validation

    def all_but_the_last_fold(self):
        return list(range(self.num_folds - 1))

    def all_but_one_fold(self, fold_to_exclude: int):
        return list(range(fold_to_exclude)) + list(range(fold_to_exclude + 1, self.num_folds))

    def all_folds(self):
        return list(range(self.num_folds))

    def make_data_loader(self, folds_to_use: List[int], batch_size: int, pin_memory=False, num_workers: int = 0, labeled_only: bool = False):
        sampler = SemiSupervisedArtifactBatchSampler(self, batch_size, folds_to_use, labeled_only)
        return DataLoader(dataset=self, batch_sampler=sampler, collate_fn=ArtifactBatch, pin_memory=pin_memory, num_workers=num_workers)


# make ArtifactBatches that mix different ref, alt counts, labeled, unlabeled
# with an option to emit only labeled data
class SemiSupervisedArtifactBatchSampler(Sampler):
    def __init__(self, dataset: ArtifactDataset, batch_size, folds_to_use: List[int], labeled_only: bool = False):
        # combine the index lists of all relevant folds
        self.indices_to_use = []

        for fold in folds_to_use:
            self.indices_to_use.extend(dataset.labeled_indices[fold])
            if not labeled_only:
                self.indices_to_use.extend(dataset.unlabeled_indices[fold])

        self.batch_size = batch_size
        self.num_batches = math.ceil(len(self.indices_to_use) // self.batch_size)

    def __iter__(self):
        random.shuffle(self.indices_to_use)
        batches = chunk(self.indices_to_use, self.batch_size)   # list of lists of indices -- each sublist is a batch
        random.shuffle(batches)

        return iter(batches)

    def __len__(self):
        return self.num_batches

