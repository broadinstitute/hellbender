package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;


/**
 * This tool combines together rows of variant calls from multiple VCFs, e.g. those produced by scattering calling
 * across genomic intervals, into a single VCF. This tool enables scattering operations, e.g. in the cloud, and is
 * preferred for such contexts over Picard MergeVcfs or Picard GatherVCfs. The tool also runs locally.
 *
 * <p>The input files need to have the same set of samples but completely different sets of loci.
 * These input files must be supplied in genomic order and must not have events at overlapping positions.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A set of VCF files, each specified in genomic order with the -I option, or a .list text file listing the set of VCFs
 * to be merged, one file per line.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A single VCF file containing the variant call records from the multiple VCFs.
 * </p>
 *
 * <h3>Usage examples</h3>
 * Specify each VCF file within the command.
 * <pre>
 * gatk GatherVcfsCloud \
 *     -I cohortA_chr1.vcf.gz \
 *     -I cohortA_chr2.vcf.gz \
 *     -O cohortA_chr1chr2.vcf.gz
 * </pre>
 *
 * Specify the VCF files using the following input.list:
 * <pre>
 *     cohortA_chr1.vcf.gz
 *     cohortA_chr2.vcf.gz
 * </pre>
 *
 * <pre>
 * gatk GatherVcfsCloud \
 *     -I input.list
 *     -O cohortA_chr1chr2.vcf.gz
 * </pre>
 *
 * @author Tim Fennell
 */

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Gathers multiple VCF files from a scatter operation into a single VCF file. Input files " +
                  "must be supplied in genomic order and must not have events at overlapping positions.",
        oneLineSummary = "Gathers multiple VCF files from a scatter operation into a single VCF file",
        programGroup = VariantManipulationProgramGroup.class
)
public final class GatherVcfsCloud extends CommandLineProgram {

    public static final String IGNORE_SAFETY_CHECKS_LONG_NAME = "ignore-safety-checks";
    public static final String GATHER_TYPE_LONG_NAME = "gather-type";

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,  doc = "Input VCF file(s).")
	public List<String> inputs;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output VCF file.")
	public File output;

    @Argument(fullName = StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, shortName = StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_SHORT_NAME, doc = "Size of the cloud-only prefetch buffer (in MB; 0 to disable).", optional=true)
    public int cloudPrefetchBuffer = 2;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_SHORT_NAME,
            doc = "If true, create a VCF index when writing a coordinate-sorted VCF file.", optional=true)
    public boolean createIndex = true;

    @Argument(fullName = GATHER_TYPE_LONG_NAME, doc ="Choose which method should be used to gather: BLOCK gathering is faster but only" +
            "works when you have both bgzipped inputs and outputs, while CONVENTIONAL gather is much slower but should work on all vcf files. " +
            "AUTOMATIC chooses BLOCK if possible and CONVENTIONAL otherwise.")
    public GatherType gatherType = GatherType.AUTOMATIC;

    @Argument(fullName = "progress-logger-frequency",
            doc = "The frequency with which the progress is logged to output (i.e. every N records)")
    public Integer progressLoggerFrequency = 10000;

    @Advanced
    @Argument(fullName = IGNORE_SAFETY_CHECKS_LONG_NAME, doc = "Disable sanity checks to improve performance, may result in silently creating corrupted outputs data")
    public boolean ignoreSafetyChecks = false;

    @Advanced
    @Argument(fullName = "disable-contig-ordering-check", doc = "Don't check relative ordering of contigs when doing a conventional gather")
    public boolean disableContigOrderingCheck = false;

    private static final Logger log = LogManager.getLogger();

    public enum GatherType { BLOCK, CONVENTIONAL, AUTOMATIC}

    @Override
    protected Object doWork() {
        log.info("Checking inputs.");
        final List<Path> inputPaths = inputs.stream().map(IOUtils::getPath).collect(Collectors.toList());

        if(!ignoreSafetyChecks) {
            for (final Path f : inputPaths) {
                IOUtil.assertFileIsReadable(f);
            }
        }

        IOUtil.assertFileIsWritable(output);

        final SAMSequenceDictionary sequenceDictionary = getHeader(inputPaths.get(0)).getSequenceDictionary();

        if (createIndex && sequenceDictionary == null) {
            throw new UserException("In order to index the resulting VCF, the input VCFs must contain ##contig lines.");
        }

        if( !ignoreSafetyChecks) {
            log.info("Checking file headers and first records to ensure compatibility.");
            assertSameSamplesAndValidOrdering(inputPaths, disableContigOrderingCheck);
        }

        if(gatherType == GatherType.AUTOMATIC) {
            if ( canBlockCopy(inputPaths, output) ) {
                gatherType = GatherType.BLOCK;
            } else {
                gatherType = GatherType.CONVENTIONAL;
            }
        }

        if (gatherType == GatherType.BLOCK && !canBlockCopy(inputPaths, output)) {
            throw new UserException.BadInput(
                    "Requested block copy but some files are not bgzipped, all inputs and the output must be bgzipped to block copy");
        }

        switch (gatherType) {
            case BLOCK:
                log.info("Gathering by copying gzip blocks. Will not be able to validate position non-overlap of files.");
                if (createIndex) {
                    log.warn("Index creation not currently supported when gathering block compressed VCFs.");
                }
                gatherWithBlockCopying(inputPaths, output, cloudPrefetchBuffer);
                break;
            case CONVENTIONAL:
                log.info("Gathering by conventional means.");
                gatherConventionally(sequenceDictionary, createIndex, inputPaths, output, cloudPrefetchBuffer, disableContigOrderingCheck, progressLoggerFrequency);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Invalid gather type: " + gatherType + ".  Please report this bug to the developers.");
        }
        return null;
    }

    private static boolean canBlockCopy(final List<Path> inputPaths, final File output) {
        return areAllBlockCompressed(inputPaths) && areAllBlockCompressed(CollectionUtil.makeList(output.toPath()));
    }

    private static VCFHeader getHeader(final Path path) {
        try (FeatureReader<VariantContext> reader =  getReaderFromVCFUri(path, 0)) {
            return ((VCFHeader) reader.getHeader());
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(path, e.getMessage(), e);
        }
    }

    /** Checks (via filename checking) that all files appear to be block compressed files. */
    @VisibleForTesting
    static boolean areAllBlockCompressed(final List<Path> input) {
        for (final Path path : input) {
            if (path == null){
                return false;
            }
            final String pathString = path.toUri().toString();
            if ( pathString.endsWith(".bcf") || !IOUtil.hasBlockCompressedExtension(pathString)){
                return false;
            }
        }
        return true;
    }

    private static FeatureReader<VariantContext> getReaderFromVCFUri(final Path variantPath, final int cloudPrefetchBuffer) {
        final String variantURI = variantPath.toUri().toString();
        return AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(), false,
                BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer),
                Function.identity());
    }

    /** Validates that all headers contain the same set of genotyped samples and that files are in order by position of first record. */
    private static void assertSameSamplesAndValidOrdering(final List<Path> inputFiles, final boolean disableContigOrderingCheck) {
        final VCFHeader firstHeader = getHeader(inputFiles.get(0));
        final SAMSequenceDictionary dict = firstHeader.getSequenceDictionary();
        if ( dict == null) {
            throw new UserException.BadInput("The first VCF specified is missing the required sequence dictionary. " +
                                                     "This is required to perform validation.  You can skip this validation " +
                                                     "using --"+IGNORE_SAFETY_CHECKS_LONG_NAME +" but ignoring safety checks " +
                                                     "can result in invalid output.");
        }
        final VariantContextComparator comparator = new VariantContextComparator(dict);
        final List<String> samples = firstHeader.getGenotypeSamples();

        Path lastFile = null;
        VariantContext lastContext = null;

        for (final Path f : inputFiles) {
            final FeatureReader<VariantContext> in = getReaderFromVCFUri(f, 0);
            final VCFHeader header = (VCFHeader)in.getHeader();
            dict.assertSameDictionary(header.getSequenceDictionary());
            final List<String> theseSamples = header.getGenotypeSamples();

            if (!samples.equals(theseSamples)) {
                final SortedSet<String> s1 = new TreeSet<>(samples);
                final SortedSet<String> s2 = new TreeSet<>(theseSamples);
                s1.removeAll(theseSamples);
                s2.removeAll(samples);

                throw new IllegalArgumentException("VCFs do not have identical sample lists." +
                        " Samples unique to first file: " + s1 + ". Samples unique to " + f.toUri().toString() + ": " + s2 + ".");
            }

            try(final CloseableIterator<VariantContext> variantIterator = in.iterator()) {
                if (variantIterator.hasNext()) {
                    final VariantContext currentContext = variantIterator.next();
                    if (lastContext != null) {
                        if ( disableContigOrderingCheck ) {
                            if ( lastContext.getContig().equals(currentContext.getContig()) && lastContext.getStart() >= currentContext.getStart() ) {
                                throw new IllegalArgumentException(
                                        "First record in file " + f.toUri().toString() + " is not after first record in " +
                                                "previous file " + lastFile.toUri().toString());
                            }
                        }
                        else {
                            if ( comparator.compare(lastContext, currentContext) >= 0 ) {
                                throw new IllegalArgumentException(
                                        "First record in file " + f.toUri().toString() + " is not after first record in " +
                                                "previous file " + lastFile.toUri().toString());
                            }
                        }
                    }

                    lastContext = currentContext;
                    lastFile = f;
                }
            } catch (final IOException e) {
                throw new UserException.CouldNotReadInputFile(f, e.getMessage(), e);
            }

            CloserUtil.close(in);
        }
    }

    /** Code for gathering multiple VCFs that works regardless of input format and output format, but can be slow. */
    private static void gatherConventionally(final SAMSequenceDictionary sequenceDictionary,
                                             final boolean createIndex,
                                             final List<Path> inputFiles,
                                             final File outputFile,
                                             final int cloudPrefetchBuffer,
                                             final boolean disableContigOrderingCheck,
                                             final int progressLoggerFrequency) {
        final EnumSet<Options> options = EnumSet.copyOf(VariantContextWriterBuilder.DEFAULT_OPTIONS);
        if (createIndex) options.add(Options.INDEX_ON_THE_FLY); else options.remove(Options.INDEX_ON_THE_FLY);
        try (final VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .setReferenceDictionary(sequenceDictionary).setOptions(options).build()) {

            final ProgressLogger progress = new ProgressLogger(log, progressLoggerFrequency);
            VariantContext lastContext = null;
            Path lastFile = null;
            VCFHeader firstHeader = null;
            VariantContextComparator comparator = null;

            for (final Path f : inputFiles) {
                try {
                    log.debug("Gathering from file: ", f.toUri().toString());
                    final FeatureReader<VariantContext> variantReader = getReaderFromVCFUri(f, cloudPrefetchBuffer);
                    final PeekableIterator<VariantContext> variantIterator;

                        variantIterator = new PeekableIterator<>(variantReader.iterator());

                    final VCFHeader header = (VCFHeader)variantReader.getHeader();

                    if (firstHeader == null) {
                        firstHeader = header;
                        out.writeHeader(firstHeader);
                        comparator = new VariantContextComparator(firstHeader.getContigLines());
                    }

                    if (lastContext != null && variantIterator.hasNext()) {
                        final VariantContext vc = variantIterator.peek();
                        if ( disableContigOrderingCheck ) {
                            // Just check start positions
                            if ( vc.getContig().equals(lastContext.getContig()) && vc.getStart() <= lastContext.getStart() ) {
                                throw new IllegalStateException("First variant in file " + f.toUri().toString() + " is at start position " + vc.getStart() +
                                        " but last variant in earlier file " + lastFile.toUri().toString() + " is at start position " + lastContext.getStart());
                            }
                        }
                        else {
                            // Check contig ordering and start positions
                            if ( comparator.compare(vc, lastContext) <= 0 ) {
                                throw new IllegalStateException("First variant in file " + f.toUri().toString() + " is at " + String.format("%s:%d", vc.getContig(), vc.getStart()) +
                                        " but last variant in earlier file " + lastFile.toUri().toString() + " is at " + String.format("%s:%d", lastContext.getContig(), lastContext.getStart()));

                            }
                        }
                    }

                    while (variantIterator.hasNext()) {
                        lastContext = variantIterator.next();
                        out.add(lastContext);
                        progress.record(lastContext.getContig(), lastContext.getStart());
                    }

                    lastFile = f;

                    CloserUtil.close(variantIterator);
                    CloserUtil.close(variantReader);
                } catch (final IOException e) {
                    throw new UserException.CouldNotReadInputFile(f, e.getMessage(), e);
                }
            }
        }
    }

    /**
     * Assumes that all inputs and outputs are block compressed VCF files and copies them without decompressing and parsing
     * most of the gzip blocks. Will decompress and parse blocks up to the one containing the end of the header in each file
     * (often the first block) and re-compress any data remaining in that block into a new block in the output file. Subsequent
     * blocks (excluding a terminator block if present) are copied directly from input to output.
     */
    private static void gatherWithBlockCopying(final List<Path> vcfs, final File output, final int cloudPrefetchBuffer) {
         try (final FileOutputStream out = new FileOutputStream(output)) {
            boolean isFirstFile = true;

            for (final Path f : vcfs) {
                log.info("Gathering " + f.toUri());
                final Function<SeekableByteChannel, SeekableByteChannel> prefetcher = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);
                try (final SeekableStream in = SeekableStreamFactory.getInstance().getStreamFor(f.toUri().toString(), prefetcher)) {
                    // a) It's good to check that the end of the file is valid and b) we need to know if there's a terminator block and not copy it
                    final BlockCompressedInputStream.FileTermination term = BlockCompressedInputStream.checkTermination(f);
                    if (term == BlockCompressedInputStream.FileTermination.DEFECTIVE) {
                        throw new UserException.MalformedFile(f.toUri() + " does not have a valid GZIP block at the end of the file.");
                    }

                    if (!isFirstFile) {
                        final BlockCompressedInputStream blockIn = new BlockCompressedInputStream(in, false);
                        boolean lastByteNewline = true;

                        int firstNonHeaderByteIndex = -1;
                        while (blockIn.available() > 0) {
                            // Read a block - blockIn.available() is guaranteed to return the bytes remaining in the block that has been
                            // read, and since we haven't consumed any yet, that is the block size.
                            final int blockLength = blockIn.available();
                            final byte[] blockContents = new byte[blockLength];
                            final int read = blockIn.read(blockContents);
                            Utils.validate(blockLength > 0 && read == blockLength, "Could not read available bytes from BlockCompressedInputStream.");

                            // Scan forward within the block to see if we can find the end of the header within this block
                            firstNonHeaderByteIndex = -1;
                            for (int i = 0; i < read; ++i) {
                                final byte b = blockContents[i];
                                final boolean thisByteNewline = (b == '\n' || b == '\r');

                                if (lastByteNewline && !thisByteNewline && b != '#') {
                                    // Aha!  Found first byte of non-header data in file!
                                    firstNonHeaderByteIndex = i;
                                    break;
                                }

                                lastByteNewline = thisByteNewline;
                            }

                            // If we found the end of the header then write the remainder of this block out as a
                            // new gzip block and then break out of the while loop
                            if (firstNonHeaderByteIndex >= 0) {
                                final BlockCompressedOutputStream blockOut = new BlockCompressedOutputStream(out, (Path)null);
                                blockOut.write(blockContents, firstNonHeaderByteIndex, blockContents.length - firstNonHeaderByteIndex);
                                blockOut.flush();
                                // Don't close blockOut because closing underlying stream would break everything
                                break;
                            }
                        }

                        if( firstNonHeaderByteIndex == -1 ){
                            log.warn("Scanned the entire file " + f.toUri().toString() + " and found no variants");
                        }
                    }

                    // Copy remainder of input stream into output stream
                    final long currentPos = in.position();
                    final long length = in.length();
                    final long skipLast = (term == BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK) ?
                            BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length : 0;
                    final long bytesToWrite = length - skipLast - currentPos;

                    IOUtil.transferByStream(in, out, bytesToWrite);
                    isFirstFile = false;
                }
            }

            // And lastly add the Terminator block and close up
            out.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }

    }

}
