from torch import nn
import torch.nn.functional as F
import torch
import math
from typing import List


class MonoDenseLayer(nn.Module):
    """
    MonoDenseLayer from Constrained Monotonic Neural Networks, Runje and Shankaranarayana, https://arxiv.org/abs/2205.11775

    It is a modification of a plain old linear layer.

    1) The output is constrained to be monotonically increasing, decreasing, or unconstrained with respect to each input

    2) Input vectors are assumed ordered with increasing features, then decreasing, then unconstrained
    """

    def __init__(self, input_dimension: int, output_dimension: int, num_increasing: int, num_decreasing: int, omit_activation: bool = False):
        super(MonoDenseLayer, self).__init__()

        self.convex_activation = torch.relu

        self.omit_activation = omit_activation

        self.num_constrained = num_increasing + num_decreasing
        num_free = input_dimension - self.num_constrained
        assert self.num_constrained <= input_dimension
        assert self.num_constrained > 0

        self.input_dimension = input_dimension
        self.output_dimension = output_dimension

        # mask has -1's for decreasing features, otherwise 1's
        # in the forward pass we multiply by the mask for convenience so that monotonically increasing AND decreasing can both
        # be treated as increasing
        self.mask = nn.Parameter(torch.ones(input_dimension), requires_grad=False)
        self.mask[num_increasing: num_increasing + num_decreasing] = -1

        self.monotonic_W = nn.Parameter(torch.empty((output_dimension, self.num_constrained)))
        nn.init.kaiming_uniform_(self.monotonic_W, a=math.sqrt(5))

        self.free_W = nn.Parameter(torch.empty((output_dimension, input_dimension - self.num_constrained))) if num_free > 0 else None
        if self.free_W is not None:
            nn.init.kaiming_uniform_(self.free_W, a=math.sqrt(5))

        self.b = nn.Parameter(torch.empty(output_dimension))
        bound = 1 / math.sqrt(input_dimension)
        nn.init.uniform_(self.b, -bound, bound)

    def forward(self, x):
        flipped = x * self.mask

        # note that monotonicity is enforced by taking the absolute value of the monotonic weight matrix
        monotonic_contribution = F.linear(flipped[:, :self.num_constrained], torch.abs(self.monotonic_W))
        free_contribution = F.linear(flipped[:, self.num_constrained:], self.free_W) if self.free_W is not None else 0

        before_activation = monotonic_contribution + free_contribution + self.b

        if self.omit_activation:
            return before_activation

        # as in the paper, we apply three nonlinear activation functions: 1) an ordinary convex activation g(x) which
        # could be a ReLU, leaky ReLU, tanh etc; 2) the concave reflection -g(-x); 3) g(x+1)-g(1) (if x < 0) or g(1) - g(1-x) (if x > 0)

        features_per_activation = self.output_dimension // 3

        left = before_activation[:, :features_per_activation]
        middle = before_activation[:, features_per_activation:(2*features_per_activation)]
        right = before_activation[:, (2*features_per_activation):]

        output1 = self.convex_activation(left)
        output2 = -self.convex_activation(-middle)
        output3 = torch.sgn(right)*(self.convex_activation(torch.ones_like(right)) - self.convex_activation(1-torch.abs(right)))

        return torch.hstack([output1, output2, output3])


class MonotonicHighwayLayer(nn.Module):
    """
    This is a purely monotonic increasing layer, like all layers but the first in the MonoDense architecture below.

    It is a highway network layer of the form:
       output = (1 - gate) * input + gate * nonlinear(input)
    where the gate is the sigmoid of some linear transformation of the input, gating is element-by-element, and the
    nonlinear function of the input is one or more monotonic dense layers as above.
    """
    def __init__(self, dim: int, num_layers: int):
        super(MonotonicHighwayLayer, self).__init__()
        self.nonlinear = MonoDense(input_dimension=dim, output_dimensions=(num_layers * [dim]), num_increasing=dim, num_decreasing=0)

        # initialize with negative bias so behavior starts near identity with gates almost closed
        self.gate_pre_sigmoid = nn.Parameter(torch.tensor(-2.0))

    def forward(self, x):
        gate = torch.sigmoid(self.gate_pre_sigmoid)

        return (1 - gate) * x + gate * self.nonlinear(x)


class MonoDense(nn.Module):
    """

    """

    def __init__(self, input_dimension: int, output_dimensions: List[int], num_increasing: int, num_decreasing):
        super(MonoDense, self).__init__()

        self.input_dimension = input_dimension
        self.layers = torch.nn.Sequential()

        last_layer_dim = input_dimension
        for layer, output_dim in enumerate(output_dimensions):
            omit_activation = (layer == len(output_dimensions) - 1)

            if output_dim > 0:
                # layers after the first are purely monotonic increasing
                n_increasing = num_increasing if layer == 0 else last_layer_dim
                n_decreasing = num_decreasing if layer == 0 else 0
                self.layers.append(MonoDenseLayer(last_layer_dim, output_dim, n_increasing, n_decreasing, omit_activation=omit_activation))
                last_layer_dim = output_dim
            else:   # negative output dimension denotes monotonic highway layer
                if layer == 0:
                    assert num_increasing == input_dimension, "initial highway layer is only valid for purely increasing network"
                num_hidden_layers = -output_dim
                self.layers.append(MonotonicHighwayLayer(dim=last_layer_dim, num_layers=num_hidden_layers))

    def forward(self, x):
        return self.layers.forward(x)
