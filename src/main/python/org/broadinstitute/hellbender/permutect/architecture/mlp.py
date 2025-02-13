import torch
from torch import nn, Tensor
from typing import List


class DenseSkipBlock(nn.Module):
    """
    computes x + f(x) where f(x) has some given number of linear layers, each with input and output dimension equal
    to that of the input x.  As suggested in arxiv:1603.05027, Identity Maps in Deep Residual Networks, nonlinearities come before each linear transformation
    """
    def __init__(self, input_size: int, num_layers: int, batch_normalize: bool = False, dropout_p: float = 0):
        super(DenseSkipBlock, self).__init__()
        self.mlp = MLP((num_layers + 1) * [input_size], batch_normalize, dropout_p, prepend_activation=True)

        # scale the MLP and initially set it to a small amount so that the block is close to an identity map early in learning
        self.alpha = nn.Parameter(torch.tensor(0.1))

    def forward(self, x):
        return x + self.alpha * self.mlp.forward(x)


class MLP(nn.Module):
    """
    A fully-connected network (multi-layer perceptron) that we need frequently
    as a sub-network.  It is parameterized by the dimensions of its layers, starting with
    the input layer and ending with the output.  Output is logits and as such no non-linearity
    is applied after the last linear transformation.
    """

    def __init__(self, layer_sizes: List[int], batch_normalize: bool = False, dropout_p: float = 0, prepend_activation: bool = False):
        super(MLP, self).__init__()

        layers = [nn.SELU()] if prepend_activation else []
        self._input_dim = layer_sizes[0]
        input_dim = layer_sizes[0]
        for k, output_dim in enumerate(layer_sizes[1:]):
            # negative output dimension -d will denote a d-layer residual skip connection
            # the output dimension of which equals the current input dimension
            if output_dim < 0:
                layers.append(DenseSkipBlock(input_dim, -output_dim, batch_normalize, dropout_p))
                continue

            if batch_normalize:
                layers.append(nn.BatchNorm1d(num_features=input_dim))

            layers.append(nn.Linear(input_dim, output_dim))

            if dropout_p > 0:
                layers.append(nn.Dropout(p=dropout_p))

            # k runs from 0 to len(layer_sizes) - 2.  Omit the nonlinearity after the last layer.
            if k < len(layer_sizes) - 2:
                layers.append(nn.SELU())

            input_dim = output_dim  # note that this does not happen for a residual skip connection

        self._output_dim = input_dim
        self._model = nn.Sequential(*layers)

    def input_dimension(self) -> int:
        return self._input_dim

    def output_dimension(self) -> int:
        return self._output_dim

    def forward(self, x: Tensor) -> Tensor:
        return self._model.forward(x)