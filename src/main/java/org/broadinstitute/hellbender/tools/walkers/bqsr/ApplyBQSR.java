package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;

/**
 * Apply base quality score recalibration
 *
 * <p>This tool performs the second pass in a two-stage process called Base Quality Score Recalibration (BQSR).
 * Specifically, it recalibrates the base qualities of the input reads based on the recalibration table produced by
 * the BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.</p>
 *
 * <h4>Summary of the BQSR procedure</h4>
 * <p>The goal of this procedure is to correct for systematic bias that affect the assignment of base quality scores
 * by the sequencer. The first pass consists of calculating error empirically and finding patterns in how error varies
 * with basecall features over all bases. The relevant observations are written to a recalibration table. The second
 * pass consists of applying numerical corrections to each individual basecall based on the patterns identified in the
 * first step (recorded in the recalibration table) and write out the recalibrated data to a new BAM or CRAM file.</p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A BAM or CRAM file containing input read data</li>
 *     <li>The covariates table (= recalibration file) generated by BaseRecalibrator on the input BAM or CRAM file</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p> A BAM or CRAM file containing the recalibrated read data</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk ApplyBQSR \
 *   -R reference.fasta \
 *   -I input.bam \
 *   --bqsr-recal-file recalibration.table \
 *   -O output.bam
 * </pre>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>This tool replaces the use of PrintReads for the application of base quality score recalibration as practiced
 * in earlier versions of GATK (2.x and 3.x).</li>
 *     <li>You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).</li>
 *     <li>Original qualities can be retained in the output file under the "OQ" tag if desired. See the
 *     `--emit-original-quals` argument for details.</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = ApplyBQSR.USAGE_SUMMARY,
        oneLineSummary = ApplyBQSR.USAGE_ONE_LINE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public final class ApplyBQSR extends ReadWalker{
    static final String USAGE_ONE_LINE_SUMMARY = "Apply base quality score recalibration";
    static final String USAGE_SUMMARY = "Apply a linear base quality recalibration model trained with the BaseRecalibrator tool.";

    private static final Logger logger = LogManager.getLogger(ApplyBQSR.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    /**
     * This argument is required for recalibration of base qualities. The recalibration table is a file produced by
     * the BaseRecalibrator tool. Please be aware that you should only run recalibration with the recalibration file
     * created on the same input data.
     */
    @Argument(fullName=StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, shortName=StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME, doc="Input recalibration table for BQSR")
    public File BQSR_RECAL_FILE;

    /**
     * Command-line arguments to fine tune the recalibration.
     */
    @ArgumentCollection
    public ApplyBQSRArgumentCollection bqsrArgs = new ApplyBQSRArgumentCollection();
    
    private SAMFileGATKReadWriter outputWriter;

    /**
     * Returns the BQSR post-transformer.
     */
    @Override
    public ReadTransformer makePostReadFilterTransformer(){
        return new BQSRReadTransformer(getHeaderForReads(), BQSR_RECAL_FILE, bqsrArgs);
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT, true);
        Utils.warnOnNonIlluminaReadGroups(getHeaderForReads(), logger);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
