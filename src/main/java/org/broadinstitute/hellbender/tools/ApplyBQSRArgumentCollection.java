package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.QualityUtils;

/**
 * The collection of all arguments needed for ApplyBQSR.
 */
public class ApplyBQSRArgumentCollection extends ApplyBQSRUniqueArgumentCollection {
    private static final long serialVersionUID = 1L;
    public static final String USE_ORIGINAL_QUALITIES_LONG_NAME = "use-original-qualities";

    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out
     * unmodified in the recalibrated BAM file. In general it's unsafe to change qualities scores below < 6, since
     * base callers use these values to indicate random or bad bases. For example, Illumina writes Q2 bases when the
     * machine has really gone wrong. This would be fine in and of itself, but when you select a subset of these reads
     * based on their ability to align to the reference and their dinucleotide effect, your Q2 bin can be elevated to
     * Q8 or Q10, leading to issues downstream.
     */
    @Argument(fullName = "preserve-qscores-less-than", doc = "Don't recalibrate bases with quality scores less than this threshold", optional = true,
            minValue = 0, minRecommendedValue =  QualityUtils.MIN_USABLE_Q_SCORE)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;

    /**
     * This flag tells GATK to use the original base qualities (that were in the data before BQSR/recalibration) which
     * are stored in the OQ tag, if they are present, rather than use the post-recalibration quality scores. If no OQ
     * tag is present for a read, the standard quality score will be used.
     */
    @Argument(fullName=USE_ORIGINAL_QUALITIES_LONG_NAME, shortName = "OQ", doc = "Use the base quality scores from the OQ tag", optional = true)
    public Boolean useOriginalBaseQualities = false;
}
