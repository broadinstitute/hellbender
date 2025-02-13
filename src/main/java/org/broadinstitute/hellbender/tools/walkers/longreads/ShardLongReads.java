package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

/**
 * Shards long read BAMs, keeping reads from the same ZMW in the same file if ZMW information is present",
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A collection of BAM files where care has been keep reads from the same ZMW in the same file</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Split reads in BAM file by sample name, read group and library name</h4>
 * <pre>
 *   gatk ShardLongReads \
 *     -I input.bam \
 *     -O outputDirectory \
 *     --num-reads-per-split 10000
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
    summary = "Shards long read BAMs and keep subreads from the same ZMW in the same file if ZMW information is present",
    oneLineSummary = "Shards long read BAMs and keep subreads from the same ZMW in the same file if ZMW information is present",
    programGroup = LongReadProgramGroup.class
)
public final class ShardLongReads extends ReadWalker {
    public static final String READS_PER_SPLIT_LONG_NAME = "num-reads-per-split";
    public static final String READS_PER_SPLIT_SHORT_NAME = "nr";

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The directory to output SAM/BAM/CRAM files."
    )
    public File OUTPUT_DIRECTORY = new File("");

    @Argument(
            fullName = READS_PER_SPLIT_LONG_NAME,
            shortName = READS_PER_SPLIT_SHORT_NAME,
            doc = "Approximate number of reads per split"
    )
    public int NUM_READS_PER_SPLIT = 10000;

    private SAMFileGATKReadWriter writer;
    private int shardIndex = 0;

    private GATKRead lastRead = null;
    private int numWrittenToShard = 0;

    @Override
    public void onTraversalStart() {
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIRECTORY);
        if ( readArguments.getReadPathSpecifiers().size() != 1 ) {
            throw new UserException("This tool only accepts a single SAM/BAM/CRAM as input");
        }
        writer = createWriter(shardIndex);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        if (numWrittenToShard < NUM_READS_PER_SPLIT) {
            writer.addRead(read);
        } else {
            int lastZmw = lastRead.hasAttribute("zm") ? lastRead.getAttributeAsInteger("zm") : 0;
            int thisZmw = read.hasAttribute("zm") ? read.getAttributeAsInteger("zm") : 1;

            if (lastZmw != thisZmw) {
                writer.close();

                shardIndex++;
                numWrittenToShard = 0;
                writer = createWriter(shardIndex);
            }

            writer.addRead(read);
        }

        lastRead = read;
        numWrittenToShard++;
    }

    @Override
    public void closeTool() {
        if ( writer != null ) {
            writer.close();
        }
    }

    /**
     * Creates SAMFileWriter instances for the read shards
     * @param shardIndex index of shard
     * @return A SAMFileWriter.
     */
    private SAMFileGATKReadWriter createWriter(final int shardIndex) {
        final String base = readArguments.getReadPathSpecifiers().get(0).getBaseName().get();
        final String extension = readArguments.getReadPathSpecifiers().get(0).getExtension().get();
        final String   key     = String.format("%06d", shardIndex);

        final GATKPath outFile = new GATKPath(IOUtil.toPath(new File(OUTPUT_DIRECTORY, base + "." + key + "." + extension)).toUri().toString());

        return createSAMWriter(outFile, true);
    }
}
