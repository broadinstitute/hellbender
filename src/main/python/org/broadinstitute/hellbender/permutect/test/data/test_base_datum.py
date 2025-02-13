import numpy as np
import torch

from permutect.data import base_datum
from permutect.utils import Variation, Label


def test_base_datum():
    num_ref_reads = 6
    num_alt_reads = 8
    num_read_features = 11
    num_info_features = 9

    ref_tensor = torch.rand(num_ref_reads, num_read_features)
    alt_tensor = torch.rand(num_alt_reads, num_read_features)
    gatk_info_tensor = torch.rand(num_info_features)
    label = Label.ARTIFACT
    source = 0

    snv_datum = base_datum.BaseDatum.from_gatk("AC", Variation.SNV, ref_tensor, alt_tensor, gatk_info_tensor, label, source)

    assert torch.equal(snv_datum.get_ref_sequence_1d(), torch.Tensor([0,1]))
    assert torch.equal(snv_datum.reads_2d, np.vstack([ref_tensor, alt_tensor]))
    assert snv_datum.get_variant_type() == Variation.SNV
    assert snv_datum.label == label

    insertion_datum = base_datum.BaseDatum.from_gatk("GT", Variation.INSERTION, ref_tensor, alt_tensor, gatk_info_tensor, label, source)
    deletion_datum = base_datum.BaseDatum.from_gatk("TT", Variation.DELETION, ref_tensor, alt_tensor, gatk_info_tensor, label, source)
    assert insertion_datum.get_variant_type() == Variation.INSERTION
    assert deletion_datum.get_variant_type() == Variation.DELETION

