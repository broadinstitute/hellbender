package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;
import java.util.List;

public class PermutectArtifactModelArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;
    @Argument(
            doc = "Dimensions of hidden layers in the aggregation subnetwork, excluding the dimension of input from lower subnetworks and the dimension (1) of the output logit. Negative values indicate residual skip connections.",
            fullName = PermutectArgumentConstants.AGGREGATION_LAYERS_NAME,
            optional = false
    )
    public List<String> aggregationLayers;

    @Argument(
            doc = "Dimensions of hidden layers in the calibration subnetwork, excluding the dimension (1) of input logit and the dimension (also 1) of the output logit.",
            fullName = PermutectArgumentConstants.CALIBRATION_LAYERS_NAME,
            optional = false
    )
    public List<String> calibrationLayers;

    @Argument(
            doc = "Dropout probability.",
            fullName = PermutectArgumentConstants.DROPOUT_P_NAME,
            optional = true
    )
    public String dropoutP = "0.0";

    @Argument(
            doc = "Flag to turn on batch normalization.",
            fullName = PermutectArgumentConstants.BATCH_NORMALIZE_NAME,
            optional = true
    )
    public boolean batchNormalize = false;

}
