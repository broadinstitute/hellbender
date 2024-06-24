package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;
import java.util.List;

/**
 * Set of arguments for the {@link  AddFlowSNVQuality}
 */
public class AddFlowSNVQualityArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;
    public static final String MAX_PHRED_SCORE_FULL_NAME = "max-phred-score";
    public static final String KEEP_SUPPLEMENTARY_ALIGNMENTS_FULL_NAME = "keep-supplementary-alignments";
    public static final String INCLUDE_QC_FAILED_READ_FULL_NAME = "include-qc-failed-read";
    public static final String SNVQ_MODE_FULL_NAME = "snvq-mode";
    public static final String OUTPUT_QUALITY_ATTRIBUTE_FULL_NAME = "output-quality-attribute";
    public static final String DEBUG_READ_NAME_FULL_NAME = "debug-read-name";
    public static final String DEBUG_COLLECT_STATS_INTO_FULL_NAME = "debug-collect-stats-into";

    public enum SnvqModeEnum {
        Legacy,
        Optimistic,
        Pessimistic,
        Geometric
    };

    /**
     *  maximum value for
     *  delta in score
     **/
    @Argument(fullName = MAX_PHRED_SCORE_FULL_NAME, doc = "Limit value for phred scores", optional = true)
    public double maxPhredScore = Double.NaN;

    /**
     *  keep supplementary alignments?
     **/
    @Argument(fullName = KEEP_SUPPLEMENTARY_ALIGNMENTS_FULL_NAME, doc = "keep supplementary alignments ?", optional = true)
    public boolean keepSupplementaryAlignments = true;

    @Advanced
    @Argument(fullName= INCLUDE_QC_FAILED_READ_FULL_NAME, doc = "include reads with QC failed flag", optional = true)
    public boolean includeQcFailedReads = true;

    /**
     * snvq computation mode
     */
    @Argument(fullName = SNVQ_MODE_FULL_NAME, doc = "snvq calculation mode.", optional = true)
    public SnvqModeEnum snvMode = SnvqModeEnum.Geometric;

    /**
     * By default this tool overwrites the QUAL field with the new qualities. Setting this argument saves the original qualities in the specified SAM tag.
     */
    @Argument(fullName = OUTPUT_QUALITY_ATTRIBUTE_FULL_NAME, doc = "alternate SAM tag to put original quality scores instead of overwriting the QUAL field. If not used, QUAL will be overwritten.", optional = true)
    public String outputQualityAttribute = null;

    /**
     *  debug read names?
     **/
    @Hidden
    @Argument(fullName = DEBUG_READ_NAME_FULL_NAME, doc = "Read names of reads to output details of as part of debugging. ", optional = true)
    public List<String> debugReadName = null;

    @Advanced
    @Hidden
    @Argument(fullName= DEBUG_COLLECT_STATS_INTO_FULL_NAME, doc = "Statistics about the reads will be output to given filename.", optional = true)
    public String debugCollectStatsInto = null;
}
