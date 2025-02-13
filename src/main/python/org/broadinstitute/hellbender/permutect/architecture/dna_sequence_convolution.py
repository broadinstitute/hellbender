from torch import nn
from math import floor


def conv_output_length(input_length, kernel_size=1, stride=1, pad=0, dilation=1, **kwargs):
    """
    Output length of 1D convolution given input length and various options.  Copied from PyTorch docs
    """
    return floor(((input_length + (2 * pad) - (dilation * (kernel_size - 1)) - 1) / stride) + 1)


def pool_output_length(input_length, kernel_size=1, stride=None, pad=0, dilation=1, **kwargs):
    """
    Output length of 1D pooling given input length and various options.  Copied from PyTorch docs.
    Differs from convolution in that stride equals kernel_size by default.
    """
    return floor(((input_length + (2 * pad) - (dilation * (kernel_size - 1)) - 1) / (kernel_size if stride is None else stride)) + 1)


INITIAL_NUM_CHANNELS = 10    # 2 x (4 + 1); the '2' is for ref and alt; the 4 is A/C/G/T; the +1 is for insertion/deletion
TOKEN_SEPARATOR = '/'
KEY_VALUE_SEPARATOR = '='


class DNASequenceConvolution(nn.Module):
    """
    A fully-connected network (multi-layer perceptron) that we need frequently
    as a sub-network.  It is parameterized by the dimensions of its layers, starting with
    the input layer and ending with the output.  Output is logits and as such no non-linearity
    is applied after the last linear transformation.

    Layer strings have the format, for example:
                    ['convolution/kernel_size=3/out_channels=64',
                     'pool/kernel_size=2',
                     'leaky_relu',
                     'convolution/kernel_size=3/dilation=2/out_channels=5',
                     'selu',
                     'flatten',
                     'linear/out_features=10']
    """

    def __init__(self, layer_strings, sequence_length):
        super(DNASequenceConvolution, self).__init__()

        # note: convention for Pytorch convolutional tensors is (batch, channel, sequence)
        last_layer_channels, last_layer_length = (INITIAL_NUM_CHANNELS, sequence_length)   # we exclude the batch dimension, which is the first

        layers = []
        for layer_string in layer_strings:
            tokens = layer_string.split(TOKEN_SEPARATOR)
            layer_type_token = tokens[0]

            kwargs = {}
            for key_value_token in tokens[1:]:
                key, value = tuple(key_value_token.split(KEY_VALUE_SEPARATOR))
                kwargs[key] = int(value)    # we're assuming all params are integers

            if layer_type_token == "convolution":
                kwargs["in_channels"] = last_layer_channels
                layers.append(nn.Conv1d(**kwargs))
                last_layer_channels, last_layer_length = kwargs["out_channels"], conv_output_length(last_layer_length, **kwargs)
            elif layer_type_token == "pool":
                assert last_layer_length > 1, "You are trying to pool a length-1 sequence, which, while defined, is silly"
                layers.append(nn.MaxPool1d(**kwargs))
                last_layer_length = pool_output_length(last_layer_length, **kwargs)
            elif layer_type_token == "leaky_relu":
                layers.append(nn.LeakyReLU())
            elif layer_type_token == "selu":
                layers.append(nn.SELU())
            elif layer_type_token == "batch_norm":
                layers.append(nn.BatchNorm1d(last_layer_channels))
            elif layer_type_token == "flatten":
                layers.append(nn.Flatten()) # by default, batch dimension is not flattened
                last_layer_channels, last_layer_length = last_layer_channels * last_layer_length, 1   # no position left, everything is a "channel"
            elif layer_type_token == "linear":
                assert last_layer_length == 1, "Trying to use fully-connected layer before data have been flattened"
                kwargs["in_features"] = last_layer_channels
                layers.append(nn.Linear(**kwargs))
                last_layer_channels = kwargs["out_features"]
            else:
                raise Exception("unsupported layer_type: " + layer_type_token)

        assert last_layer_length == 1, "data have not been flattened"
        self._output_dimension = last_layer_channels
        self._model = nn.Sequential(*layers)

    def output_dimension(self):
        return self._output_dimension

    def forward(self, x):
        """
        :param x: a batch of DNA sequences represented as a 3D tensor -- 1st index batch, 2nd index channel (A, C, G, T),
                    3rd index position in the sequence.
        :return:
        """
        return self._model.forward(x)