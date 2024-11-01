package org.broadinstitute.hellbender.tools.walkers.coverage;


import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Collect statistics on callable, uncallable, poorly mapped, and other parts of the genome
 *
 * <p>
 * A very common question about a NGS set of reads is what areas of the genome are considered callable. This tool
 * considers the coverage at each locus and emits either a per base state or a summary interval BED file that
 * partitions the genomic intervals into the following callable states:
 * <dl>
 * <dt>REF_N</dt>
 * <dd>The reference base was an N, which is not considered callable the GATK</dd>
 * <dt>PASS</dt>
 * <dd>The base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE</dd>
 * <dt>NO_COVERAGE</dt>
 * <dd>Absolutely no reads were seen at this locus, regardless of the filtering parameters</dd>
 * <dt>LOW_COVERAGE</dt>
 * <dd>There were fewer than min. depth bases at the locus, after applying filters</dd>
 * <dt>EXCESSIVE_COVERAGE</dt>
 * <dd>More than -maxDepth read at the locus, indicating some sort of mapping problem</dd>
 * <dt>POOR_MAPPING_QUALITY</dt>
 * <dd>More than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads</dd>
 * </dl>
 * </p>
 * <p/>
 * <h3>Input</h3>
 * <p>
 * A BAM file containing <b>exactly one sample</b>.
 * </p>
 * <p/>
 * <h3>Output</h3>
 * <p>
 *     A file with the callable status covering each base and a table of callable status x count of all examined bases
 * </p>
 * <h3>Usage example</h3>
 * <pre>
 *  java -jar GenomeAnalysisTK.jar \
 *     -T CallableLoci \
 *     -R reference.fasta \
 *     -I myreads.bam \
 *     -summary table.txt \
 *     -o callable_status.bed
 * </pre>
 * <p/>
 * would produce a BED file that looks like:
 * <p/>
 * <pre>
 *     20 10000000 10000864 PASS
 *     20 10000865 10000985 POOR_MAPPING_QUALITY
 *     20 10000986 10001138 PASS
 *     20 10001139 10001254 POOR_MAPPING_QUALITY
 *     20 10001255 10012255 PASS
 *     20 10012256 10012259 POOR_MAPPING_QUALITY
 *     20 10012260 10012263 PASS
 *     20 10012264 10012328 POOR_MAPPING_QUALITY
 *     20 10012329 10012550 PASS
 *     20 10012551 10012551 LOW_COVERAGE
 *     20 10012552 10012554 PASS
 *     20 10012555 10012557 LOW_COVERAGE
 *     20 10012558 10012558 PASS
 * </pre>
 * as well as a summary table that looks like:
 * <p/>
 * <pre>
 *                        state nBases
 *                        REF_N 0
 *                         PASS 996046
 *                  NO_COVERAGE 121
 *                 LOW_COVERAGE 928
 *           EXCESSIVE_COVERAGE 0
 *         POOR_MAPPING_QUALITY 2906
 * </pre>
 *
 * @author Mark DePristo / Jonn Smith
 * @since May 7, 2010 / Nov 1, 2024
 */
@DocumentedFeature(groupName = "Coverage Analysis")
@CommandLineProgramProperties(
        summary = "Collect statistics on callable, uncallable, poorly mapped, and other parts of the genome",
        oneLineSummary = "Determine callable status of loci",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class CallableLoci extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (BED or per-base format)")
    private File outputFile = null;

    @Argument(fullName = "summary", 
            doc = "Name of file for output summary", 
            optional = false)
    private File summaryFile;

    @Argument(fullName = "maxLowMAPQ", shortName = "mlmq", 
            doc = "Maximum value for MAPQ to be considered a problematic mapped read")
    private byte maxLowMAPQ = 1;

    @Argument(fullName = "minMappingQuality", shortName = "mmq", 
            doc = "Minimum mapping quality of reads to count towards depth")
    private byte minMappingQuality = 10;

    @Argument(fullName = "minBaseQuality", shortName = "mbq", 
            doc = "Minimum quality of bases to count towards depth")
    private byte minBaseQuality = 20;

    @Advanced
    @Argument(fullName = "minDepth", shortName = "minDepth", 
            doc = "Minimum QC+ read depth before a locus is considered callable")
    private int minDepth = 4;

    @Argument(fullName = "maxDepth", shortName = "maxDepth", 
            doc = "Maximum read depth before a locus is considered poorly mapped")
    private int maxDepth = -1;

    @Advanced
    @Argument(fullName = "minDepthForLowMAPQ", shortName = "mdflmq", 
            doc = "Minimum read depth before a locus is considered a potential candidate for poorly mapped")
    private int minDepthLowMAPQ = 10;

    @Argument(fullName = "maxFractionOfReadsWithLowMAPQ", shortName = "frlmq", 
            doc = "If the fraction of reads at a base with low mapping quality exceeds this value, the site may be poorly mapped")
    private double maxLowMAPQFraction = 0.1;

    @Advanced
    @Argument(fullName = "format", shortName = "format", 
            doc = "Output format")
    private OutputFormat outputFormat = OutputFormat.BED;

    private PrintStream outputStream = null;
    private PrintStream summaryStream = null;
    private Integrator integrator;
    
    public enum OutputFormat {
        BED,
        STATE_PER_BASE
    }

    public enum CalledState {
        REF_N,
        CALLABLE,
        NO_COVERAGE,
        LOW_COVERAGE,
        EXCESSIVE_COVERAGE,
        POOR_MAPPING_QUALITY
    }

    protected static class Integrator {
        final long[] counts = new long[CalledState.values().length];
        CallableBaseState state = null;
    }

    protected static class CallableBaseState {
        private SimpleInterval interval;
        private final CalledState state;

        public CallableBaseState(SimpleInterval interval, CalledState state) {
            this.interval = interval;
            this.state = state;
        }

        public SimpleInterval getLocation() {
            return interval;
        }

        public CalledState getState() {
            return state;
        }

        public boolean changingState(CalledState newState) {
            return state != newState;
        }

        public void update(SimpleInterval newInterval) {
            this.interval = new SimpleInterval(
                interval.getContig(),
                interval.getStart(),
                newInterval.getEnd()
            );
        }

        @Override
        public String toString() {
            // BED format is 0-based, so subtract 1 from start
            return String.format("%s\t%d\t%d\t%s", 
                interval.getContig(), 
                interval.getStart() - 1, 
                interval.getEnd(), 
                state);
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // Validate sample count
        if (getHeaderForReads().getReadGroups().stream()
                .map(rg -> rg.getSample())
                .distinct()
                .count() != 1) {
            throw new UserException.BadInput("CallableLoci only works for a single sample");
        }

        try {
            outputStream = new PrintStream(outputFile);
            summaryStream = new PrintStream(summaryFile);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not create output file", e);
        }

        integrator = new Integrator();
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        CalledState state;

        if (BaseUtils.isNBase(referenceContext.getBase())) {
            state = CalledState.REF_N;
        } else {
            int rawDepth = 0, QCDepth = 0, lowMAPQDepth = 0;
            
            for (PileupElement e : alignmentContext.getBasePileup()) {
                rawDepth++;

                if (e.getMappingQual() <= maxLowMAPQ) {
                    lowMAPQDepth++;
                }

                if (e.getMappingQual() >= minMappingQuality && 
                    (e.getQual() >= minBaseQuality || e.isDeletion())) {
                    QCDepth++;
                }
            }

            if (rawDepth == 0) {
                state = CalledState.NO_COVERAGE;
            } else if (rawDepth >= minDepthLowMAPQ && (double)lowMAPQDepth / rawDepth >= maxLowMAPQFraction) {
                state = CalledState.POOR_MAPPING_QUALITY;
            } else if (QCDepth < minDepth) {
                state = CalledState.LOW_COVERAGE;
            } else if (maxDepth != -1 && rawDepth >= maxDepth) {
                state = CalledState.EXCESSIVE_COVERAGE;
            } else {
                state = CalledState.CALLABLE;
            }
        }

        CallableBaseState callableState = new CallableBaseState(
            new SimpleInterval(alignmentContext.getLocation()), 
            state
        );

        // Update counts
        integrator.counts[state.ordinal()]++;

        if (outputFormat == OutputFormat.STATE_PER_BASE) {
            outputStream.println(callableState.toString());
        } else {
            // BED format - integrate adjacent regions with same state
            if (integrator.state == null) {
                integrator.state = callableState;
            } else if (callableState.getLocation().getStart() != integrator.state.getLocation().getEnd() + 1 ||
                       integrator.state.changingState(callableState.getState())) {
                outputStream.println(integrator.state.toString());
                integrator.state = callableState;
            } else {
                integrator.state.update(callableState.getLocation());
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        // Print final state for BED format
        if (outputFormat == OutputFormat.BED && integrator.state != null) {
            outputStream.println(integrator.state.toString());
        }

        // Write summary statistics
        summaryStream.printf("%30s %s%n", "state", "nBases");
        for (CalledState state : CalledState.values()) {
            summaryStream.printf("%30s %d%n", state, integrator.counts[state.ordinal()]);
        }

        // Close streams
        outputStream.close();
        summaryStream.close();
        
        return null;
    }
}
