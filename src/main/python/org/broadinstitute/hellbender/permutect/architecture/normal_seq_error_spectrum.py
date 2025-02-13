import math

import torch
from torch import nn

import matplotlib.pyplot as plt
import numpy as np

EPSILON = 0.0001

# the mean of a half-normal distribution is related to the standard deviation sigma of its corresponding normal distribution by
# sigma = mean * sqrt(pi/2)
SQRT_PI_OVER_2 = math.sqrt(math.pi / 2)


# we can't use a beta binomial for normal seq error because betas have such long tails that even if we constrain the mean
# to be small there is too large a probability of a large allele fraction.  Here we assume an underlying half normal distribution on
# the allele fraction ie it is a half normal-binomial.  Since these are not conjugate we have to explicitly sample and
# essentially perform a brute force Monte Carlo integral.
class NormalSeqErrorSpectrum(nn.Module):
    def __init__(self, num_samples: int, max_mean: float):
        super(NormalSeqErrorSpectrum, self).__init__()

        self.num_samples = num_samples

        self.max_mean = max_mean

        # this is 1/lambda parameter
        # TODO: magic constant initialization!!!
        self.mean_pre_sigmoid = torch.nn.Parameter(torch.tensor(0.0))

    def forward(self, alt_counts_1d: torch.Tensor, ref_counts_1d: torch.Tensor):
        batch_size = len(alt_counts_1d)
        fractions_2d = self.get_fractions(batch_size, self.num_samples)

        log_likelihoods_2d = torch.reshape(alt_counts_1d, (batch_size, 1)) * torch.log(fractions_2d) \
            + torch.reshape(ref_counts_1d, (batch_size, 1)) * torch.log(1 - fractions_2d)

        # average over sample dimension
        log_likelihoods_1d = torch.logsumexp(log_likelihoods_2d, dim=1) - math.log(self.num_samples)

        combinatorial_term = torch.lgamma(alt_counts_1d + ref_counts_1d + 1) - torch.lgamma(alt_counts_1d + 1) - torch.lgamma(ref_counts_1d + 1)

        return combinatorial_term + log_likelihoods_1d

    def get_mean(self):
        return torch.sigmoid(self.mean_pre_sigmoid) * self.max_mean

    def get_fractions(self, batch_size, num_samples):
        actual_mean = torch.sigmoid(self.mean_pre_sigmoid) * self.max_mean
        actual_sigma = SQRT_PI_OVER_2 * actual_mean
        normal_samples = torch.randn(batch_size, num_samples, device=actual_sigma.device)
        half_normal_samples = torch.abs(normal_samples)
        fractions_2d_unbounded = actual_sigma * half_normal_samples
        # apply tanh to constrain fractions to [0, 1), and then to [EPSILON, 1 - EPSILON] for numerical stability
        fractions_2d = EPSILON + (1 - 2*EPSILON)*torch.tanh(fractions_2d_unbounded)
        return fractions_2d

    # TODO: move this method to plotting
    def density_plot_on_axis(self, ax):
        fractions = torch.squeeze(self.get_fractions(1, 100000)).detach().numpy()
        ax.hist(fractions, bins=1000, range=[0, 1])
