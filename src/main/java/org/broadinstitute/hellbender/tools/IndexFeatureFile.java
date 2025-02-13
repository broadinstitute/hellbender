
package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.*;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.OtherProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * This tool creates an index file for the various kinds of feature-containing files supported by GATK (such as VCF
 * and BED files). An index allows querying features by a genomic interval.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk IndexFeatureFile \
 *     -I cohort.vcf.gz
 * </pre>
 * This produces the corresponding index, cohort.vcf.gz.tbi.
 */

@CommandLineProgramProperties(
        summary = "Creates an index for a feature file, e.g. VCF or BED file.",
        oneLineSummary = "Creates an index for a feature file, e.g. VCF or BED file.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
public final class IndexFeatureFile extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger(IndexFeatureFile.class);

    @Argument(shortName =StandardArgumentDefinitions.INPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
              doc = "Feature file (eg., VCF or BED file) to index. Must be in a tribble-supported format")
    public GATKPath featurePath;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              doc = "The output index file. If missing, the tool will create an index file in the same directory " +
                     "as the input file.",
              optional = true)
    public GATKPath outputPath;

    public static final int OPTIMAL_GVCF_INDEX_BIN_SIZE = 128000;
    public static final String GVCF_FILE_EXTENSION = ".g.vcf";

    @Override
    protected Object doWork() {
        if (!Files.isReadable(featurePath.toPath()) ) {
            throw new UserException.CouldNotReadInputFile(featurePath.toPath());
        }

        // Get the right codec for the file to be indexed. This call will throw an appropriate exception
        // if featureFile is not in a supported format or is unreadable.
        final FeatureCodec<? extends Feature, ?> codec = new ProgressReportingDelegatingCodec<>(
                FeatureManager.getCodecForFile(featurePath.toPath()), ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES);

        final Index index = createAppropriateIndexInMemory(codec);
        final Path indexPath = determineFileName(index);

        try {
            index.write(indexPath);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not write index to file " + indexPath, e);
        }

        logger.info("Successfully wrote index to " + indexPath);
        return indexPath.toString();
    }

    private Path determineFileName(final Index index) {
        if (outputPath != null) {
            return outputPath.toPath();
        } else if (index instanceof TabixIndex) {
            return Tribble.tabixIndexPath(featurePath.toPath());
        } else {
            return Tribble.indexPath(featurePath.toPath());
        }
    }

    private Index createAppropriateIndexInMemory(final FeatureCodec<? extends Feature, ?> codec) {
        try {
            // For block-compression files, write a Tabix index
            if (IOUtil.hasBlockCompressedExtension(featurePath.toPath())) {
                // Creating tabix indices with a non standard extensions can cause problems so we disable it
                if (outputPath != null && !outputPath.getURIString().endsWith(FileExtensions.TABIX_INDEX)) {
                    throw new UserException("The index for " + featurePath + " must be written to a file with a \"" + FileExtensions.TABIX_INDEX + "\" extension");
                }

                // TODO: this could benefit from provided sequence dictionary from reference
                // TODO: this can be an optional parameter for the tool
                return IndexFactory.createIndex(featurePath.toPath(), codec, IndexFactory.IndexType.TABIX, null);
            }
            // TODO: detection of GVCF files should not be file-extension-based. Need to come up with canonical
            // TODO: way of detecting GVCFs based on the contents (may require changes to the spec!)
            else if (featurePath.getURIString().endsWith(GVCF_FILE_EXTENSION)) {
                // Optimize GVCF indices for the use case of having a large number of GVCFs open simultaneously
                return IndexFactory.createLinearIndex(featurePath.toPath(), codec, OPTIMAL_GVCF_INDEX_BIN_SIZE);
            } else {
                // Optimize indices for other kinds of files for seek time / querying
                return IndexFactory.createDynamicIndex(featurePath.toPath(), codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
            }
        } catch (TribbleException e) {
            // Underlying cause here is usually a malformed file, but can also be things like
            // "codec does not support tabix"
            throw new UserException.CouldNotIndexFile(featurePath.toPath(), e);
        }
    }

    /**
     * This class is useful when we want to report progress when indexing. The indexing code is in htsjdk and not available to us directly.
     * The workaround is to make a special decorator 'codec' that gets called on every decoded Feature and this can use used to track progress.
     * The codec delegates all calls and reports progress as it goes.
     */
    public static final class ProgressReportingDelegatingCodec<A extends Feature, B> implements FeatureCodec<A, B> {
        private final FeatureCodec<A, B> delegatee;

        private final ProgressMeter pm;

        public ProgressReportingDelegatingCodec(final FeatureCodec<A, B> delegatee, final double secondsBetweenUpdates){
            if ( secondsBetweenUpdates <= 0.0 ) {
                throw new IllegalArgumentException("secondsBetweenUpdates must be > 0.0");
            }
            this.delegatee = delegatee;
            this.pm = new ProgressMeter(secondsBetweenUpdates);
        }

        @Override
        public Feature decodeLoc(final B b) throws IOException {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            if (!pm.started()) {
                pm.start();
            }
            final Feature f = delegatee.decodeLoc(b);
            pm.update(f);
            return f;
        }

        @Override
        public A decode(final B b) throws IOException {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            if (!pm.started()) {
                pm.start();
            }

            final A result = delegatee.decode(b);
            pm.update(result);
            return result;
        }

        @Override
        public FeatureCodecHeader readHeader(final B b) throws IOException {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            return delegatee.readHeader(b);
        }

        @Override
        public Class<A> getFeatureType() {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            return delegatee.getFeatureType();
        }

        @Override
        public B makeSourceFromStream(final InputStream bufferedInputStream) {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            return delegatee.makeSourceFromStream(bufferedInputStream);
        }

        @Override
        public LocationAware makeIndexableSourceFromStream( final InputStream bufferedInputStream) {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            return delegatee.makeIndexableSourceFromStream(bufferedInputStream);
        }

        @Override
        public boolean isDone(final B b) {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            final boolean done = delegatee.isDone(b);

            // Make sure the progress meter has been started before trying to stop it (it might not
            // have been started if there were 0 records in the file):
            if (done && pm.started()){
                pm.stop();
            }

            return done;
        }

        @Override
        public void close(final B b) {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            delegatee.close(b);
        }

        @Override
        public boolean canDecode(final String path) {
            //If there's no delegatee then we're going to say no to all questions here
            return delegatee != null && delegatee.canDecode(path);
        }

        public FeatureCodec<A, B> getDelegatee() {
            return delegatee;
        }

        @Override
        public TabixFormat getTabixFormat() {
            if (delegatee == null) {
                throw new IllegalStateException("this codec cannot be used without a delegatee.");
            }
            return delegatee.getTabixFormat();
        }
    }
}
