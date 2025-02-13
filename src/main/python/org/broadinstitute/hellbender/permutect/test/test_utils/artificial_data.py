import torch
import random
from permutect.data.base_datum import BaseDatum
from permutect.utils import Variation, Label
from numpy.random import binomial


NUM_READ_FEATURES = 5


# random isotropic Gaussian tensor, dilated by different amounts in each dimension
def make_random_tensor(mean: torch.Tensor, std: torch.Tensor) -> torch.Tensor:
    assert mean.size() == std.size()

    # TODO: random normal needs same length as mean
    return mean + std * torch.randn(len(mean))


class RandomGATKInfoGenerator:
    def __init__(self, mean: torch.Tensor, std: torch.Tensor):
        assert len(mean) == len(std)
        self.mean = mean
        self.std = std

    def generate(self) -> torch.Tensor:
        return make_random_tensor(self.mean, self.std)


class RandomReadGenerator:
    def __init__(self, mean: torch.Tensor, std: torch.Tensor):
        assert len(mean) == len(std)
        self.mean = mean
        self.std = std

    def generate(self, num_reads: int) -> torch.Tensor:
        return torch.vstack([make_random_tensor(self.mean, self.std) for _ in range(num_reads)])


def make_random_data(art_gatk_info_gen: RandomGATKInfoGenerator, var_gatk_info_gen: RandomGATKInfoGenerator, art_read_gen: RandomReadGenerator,
                     var_read_gen: RandomReadGenerator, num_data: int, artifact_fraction=0.5, unlabeled_fraction=0.1,
                     indel_fraction=0.2, ref_downsampling=10, alt_downsampling=10, is_training_data=True, vaf=0.5,
                     downsample_variants_to_match_artifacts=True):
    data = []
    for _ in range(0, num_data):
        position = random.randint(1, 1000000)

        # generate label
        artifact = random.uniform(0,1) < artifact_fraction
        unlabeled = random.uniform(0,1) < unlabeled_fraction
        label = Label.UNLABELED if unlabeled else (Label.ARTIFACT if artifact else Label.VARIANT)

        # generate variant type
        indel = random.uniform(0,1) < indel_fraction
        variant_type = (Variation.DELETION if random.uniform(0, 1) < 0.5 else Variation.INSERTION) if indel else Variation.SNV

        ref_count = ref_downsampling
        # if it's test data and a variant, we draw the alt count from the AF spectrum
        # the pd_alt_count is only relevant for variant test data
        # we assume artifact used the original alts but variant was downsampled from a het
        pd_tumor_depth = 100
        pd_alt_count = random.randint(3, 10) if artifact else binomial(pd_tumor_depth, vaf)
        alt_count = random.randint(3, 10) if (artifact or (is_training_data and downsample_variants_to_match_artifacts)) \
            else min(alt_downsampling, pd_alt_count)

        if alt_count == 0:
            continue

        gatk_info_tensor = (art_gatk_info_gen if artifact else var_gatk_info_gen).generate()

        ref_tensor = var_read_gen.generate(ref_count)
        alt_tensor = (art_read_gen if artifact else var_read_gen).generate(alt_count)

        # TODO: vary the reference sequence string?
        data.append(BaseDatum.from_gatk("GTAAAGT", variant_type, ref_tensor, alt_tensor, gatk_info_tensor, label))

    return data


# artifacts and variants are identical except 0th component of artifact read tensors all have the same sign, whereas
# each non-artifact read is randomly + or -
def make_random_strand_bias_data(num_data: int, artifact_fraction=0.5, unlabeled_fraction=0.1,
                                 ref_downsampling=10, alt_downsampling=10, is_training_data=True, vaf=0.5, num_gatk_info_features=5):
    data = []
    for _ in range(0, num_data):
        # generate label
        artifact = random.uniform(0,1) < artifact_fraction
        unlabeled = random.uniform(0, 1) < unlabeled_fraction
        label = Label.UNLABELED if unlabeled else (Label.ARTIFACT if artifact else Label.VARIANT)

        # if it's test data and a variant, we draw the alt count from the AF spectrum
        # the alt_count is only relevant for variant test data
        # we assume artifact used the original alts but variant was downsampled from a het
        depth = 100
        alt_count = random.randint(3, 10) if artifact else binomial(depth, vaf)
        alt_tensor_size = random.randint(3, 10) if (artifact or is_training_data) else min(alt_downsampling, alt_count)

        if alt_count == 0:
            continue

        gatk_info_tensor = torch.zeros(num_gatk_info_features)

        # before modifying the 0th element, it's all uniform Gaussian data
        ref_tensor = torch.randn(ref_downsampling, NUM_READ_FEATURES)
        alt_tensor = torch.randn(alt_tensor_size, NUM_READ_FEATURES)

        if artifact:
            sign = 1 if random.uniform(0,1) < 0.5 else -1
            alt_tensor[:, 0] = sign * torch.abs(alt_tensor[:, 0])

        # TODO: vary the reference sequence string?
        data.append(BaseDatum.from_gatk("TGGGAATG", Variation.SNV, ref_tensor, alt_tensor, gatk_info_tensor, label))

    return data


# good and bad data are generated by distinct gaussians
def make_two_gaussian_data(num_data, is_training_data=True, vaf=0.5, artifact_fraction=0.5, unlabeled_fraction=0.1,
                     indel_fraction=0.2, ref_downsampling=10, alt_downsampling=10, downsample_variants_to_match_artifacts=True):
    var_gatk_info_mean = torch.tensor([-1]*9)
    var_gatk_info_std = torch.tensor([1]*9)
    art_gatk_info_mean = torch.tensor([1] * 9)
    art_gatk_info_std = torch.tensor([1] * 9)

    var_gatk_info_gen = RandomGATKInfoGenerator(var_gatk_info_mean, var_gatk_info_std)
    art_gatk_info_gen = RandomGATKInfoGenerator(art_gatk_info_mean, art_gatk_info_std)

    var_read_mean = torch.tensor([-1] * 11)
    var_read_std = torch.tensor([1] * 11)
    art_read_mean = torch.tensor([1] * 11)
    art_read_std = torch.tensor([1] * 11)

    var_read_gen = RandomReadGenerator(var_read_mean, var_read_std)
    art_read_gen = RandomReadGenerator(art_read_mean, art_read_std)

    return make_random_data(art_gatk_info_gen=art_gatk_info_gen, var_gatk_info_gen=var_gatk_info_gen, art_read_gen=art_read_gen,
                            var_read_gen=var_read_gen, num_data=num_data, is_training_data=is_training_data, vaf=vaf,
                            artifact_fraction=artifact_fraction, unlabeled_fraction=unlabeled_fraction,
                            indel_fraction=indel_fraction, ref_downsampling=ref_downsampling, alt_downsampling=alt_downsampling,
                            downsample_variants_to_match_artifacts=downsample_variants_to_match_artifacts)


# good and bad data are generated by gaussians with same mean (0) but artifacts are much more spread out
def make_wide_and_narrow_gaussian_data(num_data, is_training_data=True, vaf=0.5, artifact_fraction=0.5, unlabeled_fraction=0.1,
                     indel_fraction=0.2, ref_downsampling=10, alt_downsampling=10, downsample_variants_to_match_artifacts=True):
    var_gatk_info_mean = torch.tensor([0]*9)
    var_gatk_info_std = torch.tensor([1]*9)
    art_gatk_info_mean = torch.tensor([0] * 9)
    art_gatk_info_std = torch.tensor([2] * 9)

    var_gatk_info_gen = RandomGATKInfoGenerator(var_gatk_info_mean, var_gatk_info_std)
    art_gatk_info_gen = RandomGATKInfoGenerator(art_gatk_info_mean, art_gatk_info_std)

    var_read_mean = torch.tensor([0] * 11)
    var_read_std = torch.tensor([1] * 11)
    art_read_mean = torch.tensor([0] * 11)
    art_read_std = torch.tensor([2] * 11)

    var_read_gen = RandomReadGenerator(var_read_mean, var_read_std)
    art_read_gen = RandomReadGenerator(art_read_mean, art_read_std)

    return make_random_data(art_gatk_info_gen=art_gatk_info_gen, var_gatk_info_gen=var_gatk_info_gen, art_read_gen=art_read_gen,
                            var_read_gen=var_read_gen, num_data=num_data, is_training_data=is_training_data, vaf=vaf,
                            artifact_fraction=artifact_fraction, unlabeled_fraction=unlabeled_fraction,
                            indel_fraction=indel_fraction, ref_downsampling=ref_downsampling,
                            alt_downsampling=alt_downsampling, downsample_variants_to_match_artifacts=downsample_variants_to_match_artifacts)


