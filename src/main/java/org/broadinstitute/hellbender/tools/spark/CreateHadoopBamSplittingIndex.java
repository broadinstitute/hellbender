package org.broadinstitute.hellbender.tools.spark;

import com.google.common.io.Files;
import htsjdk.samtools.*;
import htsjdk.samtools.BAMSBIIndexer;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import htsjdk.samtools.util.FileExtensions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;

/**
 * Create a Hadoop BAM splitting index and optionally a BAM index from a BAM file.
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>A BAM file</li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>BAM splitting index file</li>
 *     <li>BAM bai index (optional)</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *     gatk CreateHadoopBamSplittingIndex \
 *         -I input_reads.bam \
 *         -O input_reads.bam.sbi
 * </pre>
 * or if one wants to generate bai as well
 * <pre>
 *     gatk CreateHadoopBamSplittingIndex \
 *         -I input_reads.bam \
 *         -O input_reads.bam.sbi \
 *         --create-bai
 * </pre>
 */
@CommandLineProgramProperties(summary = "Create a Hadoop BAM splitting index and optionally a BAM index from a BAM file",
        oneLineSummary = "Create a Hadoop BAM splitting index" ,
        programGroup = OtherProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class CreateHadoopBamSplittingIndex extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger(CreateHadoopBamSplittingIndex.class);

    public static final String SPLITTING_INDEX_GRANULARITY_LONG_NAME = "splitting-index-granularity";
    public static final String CREATE_BAI_LONG_NAME = "create-bai";

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
             shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
             doc = "BAM file to create a HadoopBam splitting index for",
             optional = false )
    public File inputBam;

    @Argument(fullName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME,
            doc = "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true,
            optional=true)
    public ValidationStringency readValidationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The splitting index (SBI) file. If this is unspecified an index will be created with the same name as " +
                    "the input file but with the additional extension " + FileExtensions.SBI,
            optional = true)
    public File output;

    @Argument(fullName = SPLITTING_INDEX_GRANULARITY_LONG_NAME,
            doc = "Splitting index granularity, an entry is created in the index every this many reads.",
            optional = true)
    public long granularity = SBIIndexWriter.DEFAULT_GRANULARITY;

    @Argument(fullName = CREATE_BAI_LONG_NAME,
            doc = "Set this to create a bai index at the same time as creating a splitting index",
            optional = true)
    public boolean createBai = false;


    @Override
    public Object doWork() {
        if( granularity <= 0) {
            throw new CommandLineException.BadArgumentValue(SPLITTING_INDEX_GRANULARITY_LONG_NAME, Long.toString(granularity), "Granularity must be > 0");
        }
        final File index = getOutputFile(output, inputBam);
        if(createBai){
            createBaiAndSplittingIndex(inputBam, index, granularity, readValidationStringency);
        } else {
            createOnlySplittingIndex(inputBam, index, granularity);
        }

        return 0;
    }

    private static void createOnlySplittingIndex(final File inputBam, final File index, final long granularity) {
        assertIsBam(inputBam);
        try(SeekableStream in = new SeekableFileStream(inputBam);
            BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(index))) {
                BAMSBIIndexer.createIndex(in, out, granularity);
        } catch (final IOException e) {
            throw new UserException("Couldn't create splitting index", e);
        }
    }

    private static void createBaiAndSplittingIndex(final File inputBam, final File index, final long granularity, final ValidationStringency readValidationStringency) {
        assertIsBam(inputBam);
        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(readValidationStringency)
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .open(inputBam);
            BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(index))) {
                final SAMFileHeader header = reader.getFileHeader();
                assertBamIsCoordinateSorted(header);
                final SBIIndexWriter indexer = new SBIIndexWriter(out, granularity);

                final BAMIndexer bamIndexer = new BAMIndexer(IOUtils.replaceExtension(index, FileExtensions.BAI_INDEX), header);
                BAMFileSpan lastFilePointer = null;
                for(final SAMRecord read : reader){
                    BAMFileSpan filePointer = (BAMFileSpan) read.getFileSource().getFilePointer();
                    indexer.processRecord(filePointer.getFirstOffset());
                    bamIndexer.processAlignment(read);
                    lastFilePointer = filePointer;
                }
                long nextStart = 0;
                if (lastFilePointer != null && !lastFilePointer.getChunks().isEmpty()) {
                    nextStart = lastFilePointer.getChunks().get(0).getChunkEnd();
                }
                if (nextStart == 0) {
                    nextStart = BlockCompressedFilePointerUtil.makeFilePointer(inputBam.length()); // default to file length (in case of no reads)
                }
                indexer.finish(nextStart, inputBam.length()); // nextStart is start of next record that would be added
                bamIndexer.finish();
        } catch (final IOException e) {
            throw new UserException("Couldn't create splitting index", e);
        }
    }

    private static void assertBamIsCoordinateSorted(final SAMFileHeader header) {
        if( header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UserException.BadInput("Cannot create a " + FileExtensions.BAI_INDEX + " index for a file " +
                    "that isn't coordinate sorted.");
        }
    }


    private static void assertIsBam(final File inputBam) {
        if(!BamFileIoUtils.isBamFile(inputBam)) {
            throw new UserException.BadInput("A splitting index is only relevant for a bam file, but a "
                    + "file with extension "+ Files.getFileExtension(inputBam.getName()) + " was specified.");
        }
    }

    private static File getOutputFile(final File suggestedOutput, final File input) {
        if(suggestedOutput == null){
            return new File(input.getPath() + FileExtensions.SBI);
        } else {
            if (!suggestedOutput.getAbsolutePath().endsWith("bam" + FileExtensions.SBI)){
                logger.warn("Creating a splitting index (SBI) with an extension that doesn't match "
                        + "bam"+ FileExtensions.SBI + ".  Output file: "+suggestedOutput);
            }
            return suggestedOutput;
        }
    }
}
