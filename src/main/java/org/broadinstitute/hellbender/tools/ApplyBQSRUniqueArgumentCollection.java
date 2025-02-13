package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * The collection of those arguments for ApplyBQSR that are not already defined in RecalibrationArgumentCollection.
 * This is needed for tools (like {@link org.broadinstitute.hellbender.tools.spark.pipelines.ReadsPipelineSpark}
 * that use both ApplyBSQR and Recalibration argument collections, and which would otherwise have duplicate arguments.
 */
public class ApplyBQSRUniqueArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String STATIC_QUANTIZED_QUALS_LONG_NAME = "static-quantized-quals";
    public static final String ALLOW_MISSING_READ_GROUPS_LONG_NAME = "allow-missing-read-group";

    /**
     * Turns on the base quantization module. It requires a recalibration report.
     *
     * A value of 0 here means "do not quantize".
     * Any value greater than zero will be used to recalculate the quantization using that many levels.
     * Negative values mean that we should quantize using the recalibration report's quantization level.
     */
    @Argument(fullName="quantize-quals", doc = "Quantize quality scores to a given number of levels", optional=true)
    public int quantizationLevels = 0;


    /**
     * Static quantized quals are entirely separate from the quantize_qual option which uses dynamic binning.
     * The two types of binning should not be used together.
     *
     * For example, the Warp germline pipeline uses the static bins { 10, 20, 30, 40 }
     */
    @Advanced
    @Argument(fullName=STATIC_QUANTIZED_QUALS_LONG_NAME, doc = "Use static quantized quality scores to a given number of levels (with -"+ StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME+ ")", optional=true, mutex = "quantize-quals")
    public List<Integer> staticQuantizationQuals = new ArrayList<>();

    /**
     * Round down quantized only works with the static_quantized_quals option, and should not be used with
     * the dynamic binning option provided by quantize_quals.  When roundDown = false, rounding is done in
     * probability space to the nearest bin.  When roundDown = true, the value is rounded to the nearest bin
     * that is smaller than the current bin.
     */
    @Advanced
    @Argument(fullName="round-down-quantized", doc = "Round quals down to nearest quantized qual", optional=true, mutex = "quantize-quals")
    public boolean roundDown = false;

    /**
     * The tool is capable of writing out the original quality scores of each read in the recalibrated output file
     * under the "OQ" tag. By default, this behavior is disabled because emitting original qualities results in a
     * significant increase of the file size. Use this flag to turn on emission of original qualities.
     */
    @Argument(fullName="emit-original-quals", doc = "Emit original base qualities under the OQ tag", optional=true)
    public boolean emitOriginalQuals = false;

    /**
     * If specified, the value of this argument will be used as a flat prior for all mismatching quality scores instead
     * of the reported quality score (assigned by the sequencer).
     */
    @Argument(fullName = "global-qscore-prior", doc = "Global Qscore Bayesian prior to use for BQSR", optional = true)
    public double globalQScorePrior = -1.0;

    /**
     * If set to true, do not throw an error upon encountering a read with a read group that's not in the recalibration table.
     * Instead, simply set the quantized original base qualities as the recalibrated base qualities.
     */
    @Argument(fullName = ALLOW_MISSING_READ_GROUPS_LONG_NAME, doc = "Do not throw an error when encountering a read group not in the recal table", optional = true)
    public boolean allowMissingReadGroups = false;

    /**
     * Combine the extra arguments in {@link ApplyBQSRArgumentCollection} that are not in this {@link ApplyBQSRUniqueArgumentCollection}
     * from the given {@link RecalibrationArgumentCollection} to create a {@link ApplyBQSRArgumentCollection}.
     * @param bqsrArgs the recalibration arguments
     * @return an argument collection with all the combined arguments
     */
    public ApplyBQSRArgumentCollection toApplyBQSRArgumentCollection(RecalibrationArgumentCollection bqsrArgs) {
        ApplyBQSRArgumentCollection ret = new ApplyBQSRArgumentCollection();
        // include all the fields from ApplyBQSRArgumentCollection
        ret.quantizationLevels = this.quantizationLevels;
        ret.staticQuantizationQuals = this.staticQuantizationQuals;
        ret.roundDown = this.roundDown;
        ret.emitOriginalQuals = this.emitOriginalQuals;
        ret.globalQScorePrior = this.globalQScorePrior;
        ret.allowMissingReadGroups = this.allowMissingReadGroups;

        // include all the fields from RecalibrationArgumentCollection
        ret.PRESERVE_QSCORES_LESS_THAN = bqsrArgs.PRESERVE_QSCORES_LESS_THAN;
        ret.useOriginalBaseQualities = bqsrArgs.useOriginalBaseQualities;
        return ret;
    }
}
