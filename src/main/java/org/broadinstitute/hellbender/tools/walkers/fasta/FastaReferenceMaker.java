package org.broadinstitute.hellbender.tools.walkers.fasta;

import com.google.common.primitives.Bytes;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.IOException;
import java.nio.file.Path;

/**
 * Create a subset of a FASTA reference sequence
 *
 * <p>This tool creates a new reference in FASTA format consisting of only those positions or intervals
 * provided in the input data set. The output format can be partially controlled using the provided command-line
 * arguments. Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).</p>
 *
 * <h3>Input</h3>
 * <p>
 * The reference and requested intervals.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A fasta file representing the requested intervals. Each interval has a description line starting with a greater-than (">") symbol followed by sequence data.
 * The description begins with the contig name followed by the beginning position on the contig.
 * <pre>
 * For example, the fasta file for contig 1 and intervals 1:3-1:4 and 1:6-1:9
 * >1 1:3
 * AT
 * >1 1:6
 * GGGG
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FastaReferenceMaker \
 *   -R reference.fasta \
 *   -O output.fasta \
 *   -L input.intervals
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
    summary = "Create snippets of a fasta file by subsetting to an interval list",
    oneLineSummary = "Create snippets of a fasta file",
    programGroup = ReferenceProgramGroup.class
)
public class FastaReferenceMaker extends ReferenceWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path to write the output fasta to")
    protected String output;

    public static final String LINE_WIDTH_LONG_NAME = "line-width";
    @Argument(fullName= LINE_WIDTH_LONG_NAME, doc="Maximum length of sequence to write per line", optional=true)
    public int basesPerLine = FastaReferenceWriter.DEFAULT_BASES_PER_LINE;

    public static final String KEEP_CONTIG_NAMES = "keep-contig-names";
    @Argument(fullName = KEEP_CONTIG_NAMES, doc="Keep the original contig names intact after modification", optional=true)
    public boolean keepContigNames = false;

    protected FastaReferenceWriter writer;
    private int contigCount = 0;
    private int currentSequenceStartPosition = 0;
    private SimpleInterval lastPosition = null;
    private ByteArrayList sequence = new ByteArrayList(10000);

    @Override
    public void onTraversalStart() {
        final Path path = IOUtils.getPath(output);
        try {
            writer = new FastaReferenceWriterBuilder()
                    .setFastaFile(path)
                    .setBasesPerLine(basesPerLine)
                    .build();
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create " + output + ", encountered exception: " + e.getMessage(), e);
        }
    }

    @Override
    public void apply(ReferenceContext referenceContext, ReadsContext read, FeatureContext featureContext) {
        addToReference(referenceContext.getInterval(), referenceContext.getBase());
    }

    protected void addToReference(SimpleInterval interval, byte base) {
        advancePosition(interval);
        sequence.add(base);
    }

    protected void advancePosition(SimpleInterval interval) {
        if(lastPosition == null ){
            initializeNewSequence(interval);
        } else if ( !lastPosition.withinDistanceOf(interval, 1)) {
            finalizeSequence();
            initializeNewSequence(interval);
        }
        lastPosition = interval;
    }

    private void finalizeSequence() {
        final String description = lastPosition.getContig() + ":" + currentSequenceStartPosition + "-" + lastPosition.getEnd();
        try {
            if(keepContigNames) {
                writer.appendSequence(lastPosition.getContig(), description, basesPerLine, Bytes.toArray(sequence));
            }
            else {
                writer.appendSequence(String.valueOf(contigCount), description, basesPerLine, Bytes.toArray(sequence));
            }
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Failed while writing " + output + ".", e);
        }
    }

    private void initializeNewSequence(SimpleInterval interval) {
        lastPosition = interval;
        contigCount++;
        currentSequenceStartPosition = lastPosition.getStart();
        sequence.clear();
    }

    @Override
    public Object onTraversalSuccess(){
        finalizeSequence();
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        try{
           if( writer != null ) {
               writer.close();
           }
        } catch (IllegalStateException e){
            //sink this
        } catch (IOException e) {
            throw new UserException("Failed to write fasta due to " + e.getMessage(), e);
        }
    }
}