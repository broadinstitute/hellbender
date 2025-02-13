from .functional import revgrad
import torch
from torch import nn


class GradientReversal(nn.Module):
    def __init__(self, alpha):
        super().__init__()
        self.alpha = torch.tensor(alpha, requires_grad=False)

    def forward(self, x):
        return revgrad(x, self.alpha)

    def set_alpha(self, alpha_new):
        self.alpha = torch.tensor(alpha_new, requires_grad=False)
