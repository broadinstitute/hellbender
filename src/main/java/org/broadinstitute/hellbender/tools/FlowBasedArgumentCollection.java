package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;

public class FlowBasedArgumentCollection implements Serializable {
    private static final long serialVersionUID = 0;

    public static final String FLOW_USE_T0_TAG = "flow-use-t0-tag";
    public static final String REMOVE_LONGER_THAN_ONE_INDELS_LONG_NAME = "flow-remove-non-single-base-pair-indels";
    public static final String REMOVE_ONE_TO_ZERO_PROBS_LONG_NAME = "flow-remove-one-zero-probs";
    public static final String NUMBER_OF_POSSIBLE_PROBS_LONG_NAME = "flow-quantization-bins";
    public static final String FILLING_VALUE_LONG_NAME = "flow-fill-empty-bins-value";
    public static final String SYMMETRIC_INDELS_LONG_NAME = "flow-symmetric-indel-probs";
    public static final String REPORT_INS_OR_DEL_LONG_NAME = "flow-report-insertion-or-deletion";
    public static final String DISALLOW_LARGER_PROBS_LONG_NAME = "flow-disallow-probs-larger-than-call";
    public static final String LUMP_PROBS_LONG_NAME = "flow-lump-probs";
    public static final String PROB_SF_LONG_NAME = "flow-probability-scaling-factor";
    public static final String RETAIN_MAX_N_PROBS_BASE_LONG_NAME = "flow-retain-max-n-probs-base-format";
    public static final String FLOW_ORDER_CYCLE_LENGTH_LONG_NAME = "flow-order-cycle-length";
    public static final String FLOW_MATRIX_MODS_LONG_NAME = "flow-matrix-mods";
    public static final String FLOW_KEEP_BOUNDARY_FLOWS_LONG_NAME = "keep-boundary-flows";



    public static final double DEFAULT_FILLING_VALUE = 0.001;
    private static final boolean DEFAULT_REMOVE_LONGER_INDELS = false;
    private static final boolean DEFAULT_REMOVE_ONE_TO_ZERO = false;
    private static final boolean DEFAULT_SYMMETRIC_INDELS = false;
    private static final int DEFAULT_QUANTIZATION = 121;
    private static final boolean DEFAULT_ONLY_INS_OR_DEL = false;
    private static final boolean DEFAULT_DISALLOW_LARGER_PROBS = false;
    private static final boolean DEFAULT_LUMP_PROBS = false;
    private static final boolean DEFAULT_RETAIN_MAX_N_PROBS = false;
    private static final int DEFAULT_PROB_SCALING_FACTOR = 10;
    private static final int DEFAULT_FLOW_ORDER_CYCLE_LENGTH = 4;
    private static final boolean DEFAULT_FLOW_USE_T0_TAG = false;

    @Advanced
    @Argument(fullName = FLOW_USE_T0_TAG, doc = "Use t0 tag if exists in the read to create flow matrix", optional = true)
    public boolean useT0Tag = DEFAULT_FLOW_USE_T0_TAG;

    @Advanced
    @Argument(fullName = REMOVE_LONGER_THAN_ONE_INDELS_LONG_NAME, doc = "Should the probabilities of more then 1 indel be used", optional = true)
    public boolean removeLongerThanOneIndels = DEFAULT_REMOVE_LONGER_INDELS;

    @Advanced
    @Argument(fullName = REMOVE_ONE_TO_ZERO_PROBS_LONG_NAME, doc = "Remove probabilities of basecall of zero from non-zero genome", optional = true)
    public boolean removeOneToZeroProbs = DEFAULT_REMOVE_ONE_TO_ZERO;

    @Advanced
    @Argument(fullName = NUMBER_OF_POSSIBLE_PROBS_LONG_NAME, doc = "Number of bins for probability quantization", optional = true)
    public int probabilityQuantization = DEFAULT_QUANTIZATION;

    @Advanced
    @Argument(fullName = FILLING_VALUE_LONG_NAME, doc = "Value to fill the zeros of the matrix with", optional=true)
    public double fillingValue = DEFAULT_FILLING_VALUE;

    @Advanced
    @Argument(fullName = SYMMETRIC_INDELS_LONG_NAME, doc = "Should indel probabilities be symmetric in flow", optional=true)
    public boolean symmetricIndels = DEFAULT_SYMMETRIC_INDELS;

    @Advanced
    @Argument(fullName = REPORT_INS_OR_DEL_LONG_NAME, doc = "Report either insertion or deletion, probability, not both", optional=true)
    public boolean onlyInsOrDel = DEFAULT_ONLY_INS_OR_DEL;

    @Advanced
    @Argument(fullName = DISALLOW_LARGER_PROBS_LONG_NAME, doc = "Cap probabilities of error to 1 relative to base call", optional=true)
    public boolean disallowLargerProbs = DEFAULT_DISALLOW_LARGER_PROBS;

    @Advanced
    @Argument(fullName = LUMP_PROBS_LONG_NAME, doc = "Should all probabilities of insertion or deletion in the flow be combined together", optional=true)
    public boolean lumpProbs = DEFAULT_LUMP_PROBS;

    @Advanced
    @Argument(fullName = RETAIN_MAX_N_PROBS_BASE_LONG_NAME, doc = "Keep only hmer/2 probabilities (like in base format)", optional=true)
    public boolean retainMaxNProbs = DEFAULT_RETAIN_MAX_N_PROBS;

    @Advanced
    @Argument(fullName = PROB_SF_LONG_NAME, doc = "probability scaling factor for (phred=10) for probability quantization", optional=true)
    public int probabilityScalingFactor = DEFAULT_PROB_SCALING_FACTOR;

    @Advanced
    @Hidden
    @Argument(fullName = FLOW_ORDER_CYCLE_LENGTH_LONG_NAME, doc = "Length of flow order cycle", optional=true)
    public int flowOrderCycleLength = DEFAULT_FLOW_ORDER_CYCLE_LENGTH;

    @Advanced
    @Argument(fullName=FLOW_MATRIX_MODS_LONG_NAME, doc="Modifications instructions to the read flow matrix. " +
            "Format is src,dst{,src,dst}+. Example: 10,12,11,12 - these instructions will copy element 10 into 11 and 12", optional = true)
    public String flowMatrixMods = null;

    @Advanced
    @Argument(fullName=FLOW_KEEP_BOUNDARY_FLOWS_LONG_NAME, doc="prevent spreading of boundary flows.", optional = true)
    public boolean keepBoundaryFlows = false;

    public FlowBasedArgumentCollection() {}
}
