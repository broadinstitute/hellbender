package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;
import java.util.List;

public class PermutectBaseModelArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;
    @Argument(
            doc = "Optional pretrained base model to initialize training.",
            fullName = PermutectArgumentConstants.PRETRAINED_MODEL_NAME,
            optional = true
    )
    public String pretrainedModelName = null;

    @Argument(
            doc = "Dimensions of hidden layers in the read embedding subnetwork, including the dimension of the embedding itself. Negative values indicate residual skip connections.",
            fullName = PermutectArgumentConstants.READ_LAYERS_NAME,
            optional = false
    )
    public List<String> readLayers = null;

    @Argument(
            doc = "Hidden dimension of transformer keys and values in the self-attention layers.",
            fullName = PermutectArgumentConstants.SELF_ATTENTION_HIDDEN_DIMENSION_NAME,
            optional = false
    )
    public String selfAttentionHiddenDimension = null;

    @Argument(
            doc = "Number of symmetric gated MLP self-attention layers.",
            fullName = PermutectArgumentConstants.NUM_SELF_ATTENTION_LAYERS_NAME,
            optional = false
    )
    public String numSelfAttentionLayers = null;

    @Argument(
            doc = "Dimensions of hidden layers in the info embedding subnetwork, including the dimension of the embedding itself. Negative values indicate residual skip connections.",
            fullName = PermutectArgumentConstants.INFO_LAYERS_NAME,
            optional = false
    )
    public List<String> infoLayers = null;

    @Argument(
            doc = "Dimensions of hidden layers in the aggregation subnetwork, excluding the dimension of input from lower subnetworks and the dimension (1) of the output logit. Negative values indicate residual skip connections.",
            fullName = PermutectArgumentConstants.AGGREGATION_LAYERS_NAME,
            optional = false
    )
    public List<String> aggregationLayers = null;

    @Argument(
            doc = "List of strings specifying convolution layers of the reference sequence embedding. For example: convolution/kernel_size=3/out_channels=64 pool/kernel_size=2 leaky_relu convolution/kernel_size=3/dilation=2/out_channels=5 leaky_relu flatten linear/out_features=10.",
            fullName = PermutectArgumentConstants.REF_SEQ_LAYER_STRINGS_NAME,
            optional = false
    )
    public List<String> refSeqLayerStrings = null;

    @Argument(
            doc = "Dropout probability (default: 0.0).",
            fullName = PermutectArgumentConstants.DROPOUT_P_NAME,
            optional = true
    )
    public String dropoutP = "0.0";

    @Argument(
            doc = "Magnitude of data augmentation by randomly weighted average of read embeddings. A value of x yields random weights between 1 - x and 1 + x (default: 0.3).",
            fullName = PermutectArgumentConstants.REWEIGHTING_RANGE_NAME,
            optional = true
    )
    public String reweightingRange = "0.3";

    @Argument(
            doc = "Flag to turn on batch normalization.",
            fullName = PermutectArgumentConstants.BATCH_NORMALIZE_NAME,
            optional = true
    )
    public Boolean batchNormalize = false;
}
