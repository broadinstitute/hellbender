import enum
import numpy as np
import cyvcf2
import tarfile
import os
import torch


class ConsistentValue:
    """
    Tracks a value that once initialized, is consistent among eg all members of a dataset.  For example, all tensors
    must have the same number of columns.
    """
    def __init__(self, value=None):
        self.value = value

    def check(self, value):
        if self.value is None:
            self.value = value
        else:
            assert self.value == value, "inconsistent values"


class MutableInt:
    def __init__(self, value:int = 0):
        self.value = value

    def __str__(self):
        return str(self.value)

    def increment(self, amount: int = 1):
        self.value += amount

    def decrement(self, amount: int = 1):
        self.value -= amount

    def get_and_then_increment(self):
        self.value += 1
        return self.value - 1

    def get(self):
        return self.value

    def set(self, value: int):
        self.value = value


def gpu_if_available(exploit_mps=False) -> torch.device:
    if torch.cuda.is_available():
        d = 'cuda'
    elif exploit_mps and torch.mps.is_available():
        d = 'mps'
    else:
        d = 'cpu'
    return torch.device(d)


def downsample_tensor(tensor2d: np.ndarray, new_length: int):
    if tensor2d is None or new_length >= len(tensor2d):
        return tensor2d
    perm = np.random.permutation(len(tensor2d))
    return tensor2d[perm[:new_length]]


def get_variant_type(alt_allele, ref_allele):
    variant_size = len(alt_allele) - len(ref_allele)
    if variant_size == 0:
        return Variation.SNV
    else:
        return Variation.INSERTION if variant_size > 0 else Variation.DELETION


class Variation(enum.IntEnum):
    SNV = 0
    INSERTION = 1
    DELETION = 2
    BIG_INSERTION = 3
    BIG_DELETION = 4

    @staticmethod
    def get_type(ref_allele: str, alt_allele: str):
        diff = len(alt_allele) - len(ref_allele)
        if diff == 0:
            return Variation.SNV
        elif diff > 0:
            return Variation.BIG_INSERTION if diff > 1 else Variation.INSERTION
        else:
            return Variation.BIG_DELETION if diff < -1 else Variation.DELETION


class Call(enum.IntEnum):
    SOMATIC = 0
    ARTIFACT = 1
    SEQ_ERROR = 2
    GERMLINE = 3
    NORMAL_ARTIFACT = 4


class Epoch(enum.IntEnum):
    TRAIN = 0
    VALID = 1
    TEST = 2


class Label(enum.IntEnum):
    ARTIFACT = 0
    VARIANT = 1
    UNLABELED = 2

    @staticmethod
    def get_label(label_str: str):
        for label in Label:
            if label_str == label.name:
                return label

        raise ValueError('label is invalid: %s' % label_str)

    @staticmethod
    def is_label(label_str: str):
        for label in Label:
            if label_str == label.name:
                return True

        return False


def freeze(parameters):
    for parameter in parameters:
        parameter.requires_grad = False


def unfreeze(parameters):
    for parameter in parameters:
        if parameter.dtype.is_floating_point:   # an integer parameter isn't trainable by gradient descent
            parameter.requires_grad = True


def f_score(tp, fp, total_true):
    fn = total_true - tp
    return tp / (tp + (fp + fn) / 2)


# note: this function works for n, k, alpha, beta tensors of the same shape
# the result is computed element-wise ie result[i,j. . .] = beta_binomial(n[i,j..], k[i,j..], alpha[i,j..], beta[i,j..)
# often n, k will correspond to a batch dimension and alpha, beta correspond to a model, in which case
# unsqueezing is necessary
# NOTE: this includes the nCk factor
def beta_binomial(n, k, alpha, beta):
    combinatorial_term = torch.lgamma(n + 1) - torch.lgamma(n - k + 1) - torch.lgamma(k + 1)
    return combinatorial_term + torch.lgamma(k + alpha) + torch.lgamma(n - k + beta) + torch.lgamma(alpha + beta) \
           - torch.lgamma(n + alpha + beta) - torch.lgamma(alpha) - torch.lgamma(beta)


# note: this function works for n, k, p tensors of the same shape
# the result is computed element-wise ie result[i,j. . .] = binomial(n[i,j..], k[i,j..], p[i,j..])
# often n, k will correspond to a batch dimension and p correspond to a model, in which case
# unsqueezing is necessary
# NOTE: this includes the nCk factor
def binomial(n, k, p):
    combinatorial_term = torch.lgamma(n + 1) - torch.lgamma(n - k + 1) - torch.lgamma(k + 1)
    return combinatorial_term + k * torch.log(p) + (n - k) * torch.log(1 - p)


# note: this function works for n, k, alpha, beta tensors of the same shape
# the result is computed element-wise ie result[i,j. . .] = gamma_binomial(n[i,j..], k[i,j..], alpha[i,j..], beta[i,j..)
# often n, k will correspond to a batch dimension and alpha, beta correspond to a model, in which case
# unsqueezing is necessary
# NOTE: this includes the nCk factor
# WARNING: the approximations here only work if Gamma(f|alpha, beta) has very little density for f > 1
# see pp 2 - 4 of my notebook
def gamma_binomial(n, k, alpha, beta):
    alpha_tilde = (k + 1) * (n + 2) / (n - k + 1)
    beta_tilde = (n + 1) * (n + 2) / (n - k + 1)

    exponent_term = alpha_tilde * torch.log(beta_tilde) + alpha * torch.log(beta) -\
                    (alpha + alpha_tilde - 1) * torch.log(beta + beta_tilde)
    gamma_term = torch.lgamma(alpha + alpha_tilde - 1) - torch.lgamma(alpha) - torch.lgamma(alpha_tilde)
    return exponent_term + gamma_term - torch.log(n + 1)


# for tensor of shape (R, C...) and row counts n1, n2. . nK, return a tensor of shape (K, C...) whose 1st row is the sum of the
# first n1 rows of the input, 2nd row is the sum of the next n2 rows etc
# note that this works for arbitrary C, including empty.  That is, it works for 1D, 2D, 3D etc input.
def sums_over_rows(input_tensor: torch.Tensor, counts: torch.IntTensor):
    range_ends = torch.cumsum(counts, dim=0)
    assert range_ends[-1] == len(input_tensor)   # the counts need to add up!

    row_cumsums = torch.cumsum(input_tensor, dim=0)

    # if counts are eg 1, 2, 3 then range ends are 1, 3, 6 and we are interested in cumsums[0, 2, 5]
    relevant_cumsums = row_cumsums[(range_ends - 1).long()]

    # if counts are eg 1, 2, 3 we now have, the sum of the first 1, 3, and 6 rows.  To get the sums of row 0, rows 1-2, rows 3-5
    # we need the consecutive differences, with a row of zeroes prepended
    row_of_zeroes = torch.zeros_like(relevant_cumsums[0])[None] # the [None] makes it (1xC)
    relevant_sums = torch.diff(relevant_cumsums, dim=0, prepend=row_of_zeroes)
    return relevant_sums


# same but divide by the counts to get means
def means_over_rows(input_tensor: torch.Tensor, counts: torch.IntTensor, keepdim: bool = False):
    extra_dims = (1,) * (input_tensor.dim() - 1)
    result = sums_over_rows(input_tensor, counts) / counts.view(-1, *extra_dims)

    return torch.repeat_interleave(result, dim=0, repeats=counts) if keepdim else result


# given 3d tensor T_ijk and 1D index tensors I, J, K, return the 1D tensor:
# result[n] = T[I[n], J[n], K[n]]
def index_3d_array(tens, idx0, idx1, idx2):
    dim0, dim1, dim2 = tens.shape
    flattened_indices = (idx0 * dim1 * dim2) + (idx1 * dim2) + idx2
    return tens.view(-1)[flattened_indices]


def index_2d_array(tens, idx0, idx1):
    dim0, dim1 = tens.shape
    flattened_indices = (idx0 * dim1) + idx1
    return tens.view(-1)[flattened_indices]


# same but include a regularizer in case of zeros in the counts vector
# regularizer has the dimension of one row of the input tensor
def means_over_rows_with_regularizer(input_tensor: torch.Tensor, counts: torch.IntTensor, regularizer, regularizer_weight, keepdim: bool = False):
    extra_dims = (1,) * (input_tensor.dim() - 1)

    regularized_sums = sums_over_rows(input_tensor, counts) + (regularizer_weight * regularizer)[None, :]
    regularized_counts = counts + regularizer_weight
    result = regularized_sums / regularized_counts.view(-1, *extra_dims)

    return torch.repeat_interleave(result, dim=0, repeats=counts) if keepdim else result


class StreamingAverage:
    def __init__(self):
        self._count = 0.0
        self._sum = 0.0

    def is_empty(self) -> bool:
        return self._count == 0.0

    def get(self) -> float:
        return self._sum / (self._count + 0.0001)

    def record(self, value: float, weight: float=1):
        self._count += weight
        self._sum += value * weight

    def record_sum(self, value_sum: float, count):
        self._count += count
        self._sum += value_sum

    # record only values masked as true
    def record_with_mask(self, values: torch.Tensor, mask: torch.Tensor):
        self._count += torch.sum(mask).item()
        self._sum += torch.sum(values*mask).item()

    # record values with different weights
    # values and mask should live on same device as self._sum
    def record_with_weights(self, values: torch.Tensor, weights: torch.Tensor):
        self._count += torch.sum(weights).item()
        self._sum += torch.sum(values * weights).item()


def log_binomial_coefficient(n: torch.Tensor, k: torch.Tensor):
    return (n + 1).lgamma() - (k + 1).lgamma() - ((n - k) + 1).lgamma()


def backpropagate(optimizer: torch.optim.Optimizer, loss: torch.Tensor):
    optimizer.zero_grad(set_to_none=True)
    loss.backward()
    optimizer.step()


def find_variant_type(v: cyvcf2.Variant):
    alt = v.ALT[0]  # TODO: we're assuming biallelic
    ref = v.REF
    return Variation.get_type(ref, alt)


def extract_to_temp_dir(tar_file, directory):
    tar = tarfile.open(tar_file)
    tar.extractall(directory)
    tar.close()
    return [os.path.abspath(os.path.join(directory, p)) for p in os.listdir(directory)]


def trim_alleles_on_right(ref: str, alt: str):
    trimmed_ref, trimmed_alt = ref, alt
    while len(trimmed_ref) > 1 and len(trimmed_alt) > 1 and trimmed_alt[-1] == trimmed_ref[-1]:
        trimmed_ref, trimmed_alt = trimmed_ref[:-1], trimmed_alt[:-1]
    return trimmed_ref, trimmed_alt