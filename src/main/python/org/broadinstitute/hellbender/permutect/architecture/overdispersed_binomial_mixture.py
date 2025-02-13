import math
from typing import List

import torch
from permutect import utils
from torch import nn, exp, unsqueeze, logsumexp
from torch.nn.functional import softmax, log_softmax

from permutect.architecture.mlp import MLP
from permutect.metrics.plotting import simple_plot
from permutect.utils import beta_binomial, gamma_binomial, binomial, Variation


class OverdispersedBinomialMixture(nn.Module):
    """
    This model takes in 1D tensor inputs (variant type indices by batch) and as a function of input has a Beta OR Gamma mixture model.  That is, it computes for each input
    vector 1) a vector of mixture component weights 2) a vector of the alpha shape parameters of each component and 3) a
    vector of beta shape parameters of each component.  Due to batching these are all represented as 2D tensors.

    Note that both the Beta and Gamma distributions' shape parameters are traditionally called "alpha" and "beta".

    The computed likelihoods take in a 1D batch of total counts n and 1D batch of "success" counts k.

    It optionally has a max mean that scales every mean to some amount less than or equal to 1, which is useful when we want
    to force the mixture to represent only small fractions.

    When using a BetaBinomial mixture, due to conjugacy the integral over the latent probability of success (in our uses,
    this is the allele fraction of variants or artifacts) is exact and we use a closed form analytic expression for the
    density of a BetaBinomial.  That is, the probability (k = alt count, n = depth, f = latent allele fraction)

    P(k|n, alpha, beta) = integral{Beta(f|alpha, beta) * Binomial(k|n, f)}

    is exact.

    When using a GammaBinomial mixture, i.e. one with a Gamma prior Gamma(f, alpha, beta) the cannot do the integral exactly.
    However, the *binomial* factor Binom(k|n,f), which as a function of f is a Beta distribution, is extremely well-approximated
    by a Gamma distribution, and the product of this Gamma approximation and the Gamma prior on f *is* exactly integrable.

    The approximation breaks down if the allele fractions are not small (since then the support of the Gamma for f > 1
    becaomes significant), so we should only use the Gamma prior version to model artifact allele fractions.

    In addition to 'beta' and 'gamma' modes, there is also the 'none' mode which has no overdispersion in the individual components.
    That is, each component is a plain binomial, though of course by virtue of being a mixture the distribution as a whole is overdispersed.
    """

    def __init__(self, num_components: int, max_mean: float = 1, mode: str = 'beta'):
        super(OverdispersedBinomialMixture, self).__init__()
        self.mode = mode
        self.K = num_components
        self.V = len(Variation)
        self.max_mean = max_mean

        # parameters for each component and variant type:
        self.weights_pre_softmax_vk = torch.nn.Parameter(torch.ones(self.V, self.K))
        self.mean_pre_sigmoid_vk = torch.nn.Parameter(torch.randn(self.V, self.K))
        self.concentration_pre_sigmoid_vk = torch.nn.Parameter(torch.randn(self.V, self.K))
        self.max_concentration = torch.nn.Parameter(torch.tensor(50.0))

    '''
    here x is a 2D tensor, 1st dimension batch, 2nd dimension being features that determine which Beta mixture to use
    n and k are 1D tensors, the only dimension being batch.
    '''
    def forward(self, types_b, n_b, k_b):
        types_idx = types_b.long()
        log_weights_bk = log_softmax(self.weights_pre_softmax_vk[types_idx, :], dim=-1)

        # we make them 2D, with 1st dim batch, to match alpha and beta.  A single column is OK because the single value of
        # n/k are broadcast over all mixture components
        n_bk = n_b[:, None]
        k_bk = k_b[:, None]

        # 2D tensors -- 1st dim batch, 2nd dim mixture component
        mean_bk = self.max_mean * torch.sigmoid(self.mean_pre_sigmoid_vk[types_idx, :])
        concentration_bk = self.get_concentration(types_b)

        if self.mode == 'beta':
            alpha_bk = mean_bk * concentration_bk
            beta_bk = (1 - mean_bk) * concentration_bk
            log_likelihoods_bk = beta_binomial(n_bk, k_bk, alpha_bk, beta_bk)
        elif self.mode == 'gamma':
            alpha_bk = mean_bk * concentration_bk
            beta_bk = concentration_bk
            log_likelihoods_bk = gamma_binomial(n_bk, k_bk, alpha_bk, beta_bk)
        elif self.mode == 'none':
            # each mean is the center of a binomial
            log_likelihoods_bk = binomial(n_bk, k_bk, mean_bk)
        else:
            raise Exception("we don't have that kind of mode!")

        log_weighted_likelihoods_bk = log_weights_bk + log_likelihoods_bk

        # yields one number per batch, squeezed into 1D output tensor
        return logsumexp(log_weighted_likelihoods_bk, dim=-1, keepdim=False)

    def get_concentration(self, types_b):
        return self.max_concentration * torch.sigmoid(self.concentration_pre_sigmoid_vk[types_b.long(),:])

    # given 1D input tensor, return 1D tensors of component alphas and betas
    def component_shapes(self, var_type: int):
        means_k = self.max_mean * torch.sigmoid(self.mean_pre_sigmoid_vk[var_type])
        concentrations_k = self.max_concentration * torch.sigmoid(self.concentration_pre_sigmoid_vk[var_type])
        alphas_k = means_k * concentrations_k
        betas_k = (1 - means_k) * concentrations_k if self.mode == 'beta' else concentrations_k
        return alphas_k, betas_k

    def component_weights(self, var_type: int):
        return softmax(self.weights_pre_softmax_vk[var_type], dim=-1)

    # given variant type, return the moments E[x], E[ln(x)], and E[x ln(x)] of the underlying beta mixture
    def moments_of_underlying_beta_mixture(self, var_type: int):
        assert self.mode == 'beta'
        alphas, betas = self.component_shapes(var_type)
        weights = self.component_weights(var_type)

        # E[x]
        component_means = alphas / (alphas + betas)
        mixture_mean = torch.sum(weights * component_means)

        # E[ln(x)]
        component_log_means = torch.digamma(alphas) - torch.digamma(
            alphas + betas)  # digamma broadcasts to make 1D tensor
        mixture_log_mean = torch.sum(weights * component_log_means)

        # E[x ln(x)]
        component_log_linear_means = component_means * (torch.digamma(alphas + 1) - torch.digamma(alphas + betas + 1))
        mixture_log_linear_mean = torch.sum(weights * component_log_linear_means)

        return mixture_mean, mixture_log_mean, mixture_log_linear_mean

    '''
    here x is a 2D tensor, 1st dimension batch, 2nd dimension being features that determine which Beta mixture to use
    n is a 1D tensor, the only dimension being batch, and we sample a 1D tensor of k's
    '''
    def sample(self, types_b, n):
        # compute weights and select one mixture component from the corresponding multinomial for each datum / row
        weights = softmax(self.weights_pre_softmax_vk[types_b, :], dim=-1)
        component_indices = torch.multinomial(weights, num_samples=1, replacement=True)  # 2D tensor with one column

        # get 1D tensors of one selected alpha and beta shape parameter per datum / row, then sample a fraction from each
        # It may be very wasteful computing everything and only using one component, but this is just for unit testing
        means = self.max_mean * torch.sigmoid(self.mean_pre_sigmoid_vk[types_b, :].detach()).gather(dim=1, index=component_indices).squeeze()
        concentrations = self.get_concentration(types_b).detach().gather(dim=1, index=component_indices).squeeze()
        alphas = means * concentrations
        betas = (1 - means) * concentrations if self.mode == 'beta' else concentrations
        dist = torch.distributions.beta.Beta(alphas, betas) if self.mode == 'beta' else torch.distributions.gamma.Gamma(alphas, betas)
        fractions = dist.sample()  # 1D tensor

        # recall, n and fractions are 1D tensors; result is also 1D tensor, one "success" count per datum
        return torch.distributions.binomial.Binomial(total_count=n, probs=fractions).sample()

    def fit(self, num_epochs, types_b, depths_1d_tensor, alt_counts_1d_tensor, batch_size=64):
        optimizer = torch.optim.Adam(self.parameters())
        num_batches = math.ceil(len(alt_counts_1d_tensor) / batch_size)

        for epoch in range(num_epochs):
            for batch in range(num_batches):
                batch_start = batch * batch_size
                batch_end = min(batch_start + batch_size, len(alt_counts_1d_tensor))
                batch_slice = slice(batch_start, batch_end)
                loss = -torch.mean(self.forward(types_b[batch_slice], depths_1d_tensor[batch_slice],
                                                alt_counts_1d_tensor[batch_slice]))
                utils.backpropagate(optimizer, loss)

    '''
    get raw data for a spectrum plot of probability density vs allele fraction.  
    here x is a 1D tensor, a single datum/row of the 2D tensors as above
    '''
    def spectrum_density_vs_fraction(self, variant_type: Variation, depth: int):
        # device = self.mean_pre_sigmoid_vk.device
        fractions = torch.arange(0.01, 0.99, 0.001)  # 1D tensor on CPU

        log_weights_k = log_softmax(self.weights_pre_softmax_vk[variant_type].detach(), dim=-1).cpu()
        means_k = self.max_mean * torch.sigmoid(self.mean_pre_sigmoid_vk[variant_type].detach()).cpu()

        # now we're on CPU
        if self.mode == 'none':
            # this is copied from the beta case below -- basically we smear each delta function / discrete binomial
            # into a narrow Gaussian
            dist = torch.distributions.normal.Normal(means_k, 0.01 * torch.ones_like(means_k))
            densities = exp(torch.logsumexp(log_weights_k + dist.log_prob(fractions.unsqueeze(dim=1)), dim=1,
                                            keepdim=False))  # 1D tensor
            return fractions, densities
        else:
            concentrations_k = self.max_concentration.cpu() * torch.sigmoid(self.concentration_pre_sigmoid_vk[variant_type]).detach().cpu()
            alphas_k = means_k * concentrations_k
            betas_k = (1 - means_k) * concentrations_k if self.mode == 'beta' else concentrations_k

            # since f.unsqueeze(dim=1) is 2D column vector, log_prob produces 2D tensor where row index is f and column index is mixture component
            # adding the single-row 2D tensor log_weights broadcasts to each row / value of f
            # then we apply log_sum_exp, dim= 1, to sum over components and get a log_density for each f
            dist = torch.distributions.beta.Beta(alphas_k, betas_k) if self.mode == 'beta' else torch.distributions.gamma.Gamma(alphas_k, betas_k)
            densities = exp(torch.logsumexp(log_weights_k + dist.log_prob(fractions.unsqueeze(dim=1)), dim=1, keepdim=False))  # 1D tensor

            return fractions, densities

    '''
    here x is a 1D tensor, a single datum/row of the 2D tensors as above
    '''
    def plot_spectrum(self, x, title, depth: int):
        fractions, densities = self.spectrum_density_vs_fraction(x, depth)
        return simple_plot([(fractions.numpy(), densities.numpy(), " ")], "AF", "density", title)


