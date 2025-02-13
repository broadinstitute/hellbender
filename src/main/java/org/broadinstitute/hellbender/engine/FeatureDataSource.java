package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.*;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBConstants;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.utils.IndexUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.codecs.FeatureReaderFactory;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.reader.GenomicsDBFeatureReader;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import static org.broadinstitute.hellbender.tools.genomicsdb.GATKGenomicsDBUtils.createExportConfiguration;
import static org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.BCI_FILE_EXTENSION;

/**
 * Enables traversals and queries over sources of Features, which are metadata associated with a location
 * on the genome in a format supported by our file parsing framework, Tribble. Examples of Features are
 * VCF records and hapmap records.
 * <p>
 * Two basic operations are available on this data source:
 * <p>
 * -Iteration over all Features in this data source, optionally restricted to Features overlapping
 * a set of intervals if intervals are provided via {@link #setIntervalsForTraversal(List)}. Traversal
 * by a set of intervals requires the file to have been indexed using the bundled tool IndexFeatureFile.
 * The set of intervals provided MUST be non-overlapping and sorted in increasing order of start position.
 * <p>
 * -Targeted queries by one interval at a time. This also requires the file to have been indexed using
 * the bundled tool IndexFeatureFile. Targeted queries by one interval at a time are unaffected by
 * any intervals for full traversal set via {@link #setIntervalsForTraversal(List)}.
 * <p>
 * To improve performance in the case of targeted queries by one interval at a time, this class uses a caching
 * scheme that is optimized for the common access pattern of multiple separate queries over intervals with
 * gradually increasing start positions. It optimizes for this use case by pre-fetching records immediately
 * following each interval during a query and caching them. Performance will suffer if the access pattern is
 * random, involves queries over intervals with DECREASING start positions instead of INCREASING start positions,
 * or involves lots of very large jumps forward on the genome or lots of contig switches. Query caching
 * can be disabled, if desired.
 *
 * @param <T> The type of Feature returned by this data source
 */
public final class FeatureDataSource<T extends Feature> implements GATKDataSource<T>, AutoCloseable {
    private static final Logger logger = LogManager.getLogger(FeatureDataSource.class);

    /**
     * Feature reader used to retrieve records from our file
     */
    private final FeatureReader<T> featureReader;

    /**
     * Iterator representing an open traversal over this data source initiated via a call to {@link #iterator}
     * (null if there is no open traversal). We need this to ensure that each iterator is properly closed,
     * and to enforce the constraint (required by Tribble) that we never have more than one iterator open
     * over our feature reader.
     */
    private CloseableTribbleIterator<T> currentIterator;

    /**
     * Our intervals for traversal. If set, restricts full traversals initiated via {@link #iterator} to
     * return only Features overlapping this set of intervals. Does not affect individual queries
     * initiated via {@link #query(SimpleInterval)} and/or {@link #queryAndPrefetch(Locatable)}.
     */
    private List<SimpleInterval> intervalsForTraversal;

    /**
     * Cache containing Features from recent queries initiated via {@link #query(SimpleInterval)} and/or
     * {@link #queryAndPrefetch(Locatable)}. This is guaranteed to start at the start position of the
     * most recent query, but will typically end well after the end of the most recent query. Designed to
     * improve performance of the common access pattern involving multiple queries across nearby intervals
     * with gradually increasing start positions.
     */
    private final FeatureCache<T> queryCache;

    /**
     * When we experience a cache miss (ie., a query interval not fully contained within our cache) and need
     * to re-populate the Feature cache from disk to satisfy a query, this controls the number of extra bases
     * AFTER the end of our interval to fetch. Should be sufficiently large so that typically a significant number
     * of subsequent queries will be cache hits (ie., query intervals fully contained within our cache) before
     * we have another cache miss and need to go to disk again.
     */
    private final int queryLookaheadBases;

    /**
     * Holds information about the path this datasource reads from.
     */
    private final FeatureInput<T> featureInput;

    /**
     * True if this datasource is backed by a file that has an associated index file, false if it doesn't
     */
    private final boolean hasIndex;

    /**
     * True if this datasource supports efficient random access queries.
     * <p>
     * For a file, this is the same as {@link #hasIndex}, but there are non-file data sources (eg., GenomicsDB)
     * that don't have a separate index file but do support random access.
     */
    private final boolean supportsRandomAccess;

    /**
     * Default value for queryLookaheadBases, if none is specified. This is designed to be large enough
     * so that in typical usage (ie., query intervals with gradually increasing start locations) there will
     * be a substantial number of cache hits between cache misses, reducing the number of times we need to
     * repopulate the cache from disk.
     */
    public static final int DEFAULT_QUERY_LOOKAHEAD_BASES = 1000;

    /**
     * Creates a FeatureDataSource backed by the provided File. The data source will have an automatically
     * generated name, and will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES})
     * during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     */
    public FeatureDataSource(final File featureFile) {
        this(featureFile, null);
    }

    /**
     * Creates a FeatureDataSource backed by the provided path. The data source will have an automatically
     * generated name, and will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES})
     * during queries that produce cache misses.
     *
     * @param featurePath path or URI to source of Features
     */
    public FeatureDataSource(final String featurePath) {
        this(featurePath, null, DEFAULT_QUERY_LOOKAHEAD_BASES, Feature.class);
    }

    /**
     * Creates a FeatureDataSource backed by the provided File and assigns this data source the specified logical
     * name. We will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES}) during queries
     * that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param name        logical name for this data source (may be null)
     */
    public FeatureDataSource(final File featureFile, final String name) {
        this(featureFile, name, DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    /**
     * Creates a FeatureDataSource backed by the provided File and assigns this data source the specified logical
     * name. We will look ahead the specified number of bases during queries that produce cache misses.
     *
     * @param featureFile         file containing Features
     * @param name                logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     */
    public FeatureDataSource(final File featureFile, final String name, final int queryLookaheadBases) {
        this(Utils.nonNull(featureFile).getAbsolutePath(), name, queryLookaheadBases, Feature.class);
    }

    /**
     * Creates a FeatureDataSource backed by the resource at the provided path.
     *
     * @param featurePath         path to file or GenomicsDB url containing features
     * @param name                logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType   When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                            that produce this type of Feature. May be null, which results in an unrestricted search.
     */
    public FeatureDataSource(final String featurePath, final String name, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType) {
        this(new FeatureInput<>(featurePath, name != null ? name : featurePath), queryLookaheadBases, targetFeatureType);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput        a FeatureInput specifying a source of Features
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType   When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                            that produce this type of Feature. May be null, which results in an unrestricted search.
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType) {
        this(featureInput, queryLookaheadBases, targetFeatureType, 0, 0);
    }

    /**
     * Creates a FeatureDataSource backed by the resource at the provided path.
     *
     * @param featurePath              path to file or GenomicsDB url containing features
     * @param name                     logical name for this data source (may be null)
     * @param queryLookaheadBases      look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType        When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                                 that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer      MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     */
    public FeatureDataSource(final String featurePath, final String name, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer) {
        this(new FeatureInput<>(featurePath, name != null ? name : featurePath), queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput             a FeatureInput specifying a source of Features
     * @param queryLookaheadBases      look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType        When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                                 that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer      MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer) {
        this(featureInput, queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
             new GenomicsDBOptions(), false);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput             a FeatureInput specifying a source of Features
     * @param queryLookaheadBases      look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType        When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                                 that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer      MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param reference                 the reference genome corresponding to the data to be read
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final Path reference) {
        this(featureInput, queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                new GenomicsDBOptions(reference), false);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput             a FeatureInput specifying a source of Features
     * @param queryLookaheadBases      look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType        When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                                 that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer      MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param reference                 the reference genome corresponding to the data to be read
     * @param setNameOnCodec            If true, and if this FeatureDataSource uses a NameAwareCodec, the name of the FeatureInput will be used to set the codec's name. This exists as a mechanism to store the FeatureInput name in the source field of VariantContexts
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final Path reference, final boolean setNameOnCodec) {
        this(featureInput, queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                new GenomicsDBOptions(reference), setNameOnCodec);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput             a FeatureInput specifying a source of Features
     * @param queryLookaheadBases      look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType        When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                                 that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer      MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param genomicsDBOptions         options and info for reading from a GenomicsDB; may be null
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final GenomicsDBOptions genomicsDBOptions) {
        this(featureInput, queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                genomicsDBOptions, false);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput             a FeatureInput specifying a source of Features
     * @param queryLookaheadBases      look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType        When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                                 that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer      MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param genomicsDBOptions         options and info for reading from a GenomicsDB; may be null
     * @param setNameOnCodec            If true, and if this FeatureDataSource uses a NameAwareCodec, the name of the FeatureInput will be used to set the codec's name. This exists as a mechanism to store the FeatureInput name in the source field of VariantContexts
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final GenomicsDBOptions genomicsDBOptions,
                             final boolean setNameOnCodec) {
        Utils.validateArg(queryLookaheadBases >= 0, "Query lookahead bases must be >= 0");
        this.featureInput = Utils.nonNull(featureInput, "featureInput must not be null");
        if (IOUtils.isGenomicsDBPath(featureInput)) {
            Utils.nonNull(genomicsDBOptions, "GenomicsDBOptions must not be null. Calling tool may not read from a GenomicsDB data source.");
        }

        // Create a feature reader without requiring an index.  We will require one ourselves as soon as
        // a query by interval is attempted.
        this.featureReader = getFeatureReader(featureInput, targetFeatureType,
                cloudPrefetchBuffer, cloudIndexPrefetchBuffer, genomicsDBOptions, setNameOnCodec);

        if (IOUtils.isGenomicsDBPath(featureInput)) {
            //genomics db uri's have no associated index file to read from, but they do support random access
            this.hasIndex = false;
            this.supportsRandomAccess = true;
        } else if (featureReader instanceof AbstractFeatureReader) {
            this.hasIndex = ((AbstractFeatureReader<T, ?>)featureReader).hasIndex();
            this.supportsRandomAccess = hasIndex;
        } else {
            this.hasIndex = false;
            this.supportsRandomAccess = featureReader.isQueryable();
        }

        this.currentIterator = null;
        this.intervalsForTraversal = null;
        this.queryCache = new FeatureCache<>();
        this.queryLookaheadBases = queryLookaheadBases;
    }

    final void printCacheStats() {
        queryCache.printCacheStatistics( getName() );
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static <T extends Feature> FeatureReader<T> getFeatureReader(final FeatureInput<T> featureInput, final Class<? extends Feature> targetFeatureType,
                                                                         final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer,
                                                                         final GenomicsDBOptions genomicsDBOptions, final boolean setNameOnCodec) {
        if (IOUtils.isGenomicsDBPath(featureInput.getFeaturePath())) {
            Utils.nonNull(genomicsDBOptions);
            try {
                if (genomicsDBOptions.getReference() == null) {
                    throw new UserException.MissingReference("You must provide a reference if you want to load from GenomicsDB");
                }
                try {
                    final File referenceAsFile = genomicsDBOptions.getReference().toFile();
                    return (FeatureReader<T>)getGenomicsDBFeatureReader(featureInput, referenceAsFile, genomicsDBOptions);
                } catch (final UnsupportedOperationException e){
                    throw new UserException.BadInput("GenomicsDB requires that the reference be a local file.", e);
                }
            } catch (final ClassCastException e) {
                throw new UserException("GenomicsDB inputs can only be used to provide VariantContexts.", e);
            }
        }
        final FeatureCodec<T, ?> codec = getCodecForFeatureInput(featureInput, targetFeatureType, setNameOnCodec);
        if ( codec instanceof FeatureReaderFactory ) {
            return ((FeatureReaderFactory)codec).getReader(featureInput, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
        }
        return getTribbleFeatureReader(featureInput, codec, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
    }

    /**
     * Get a new FeatureCodec instance to use for a FeatureInput. Avoid re-discovering which codec class to
     * use by checking to see if the FeatureInput already has a cached codec class. It not, discover the codec class
     * and cache it for next time.
     *
     * @return A new FeatureCodec instance to use for the FeatureInput.
     */
    @SuppressWarnings("unchecked")
    private static <T extends Feature> FeatureCodec<T, ?> getCodecForFeatureInput(final FeatureInput<T> featureInput,
                                                                                  final Class<? extends Feature> targetFeatureType,
                                                                                  final boolean setNameOnCodec) {
        final FeatureCodec<T, ?> codec;
        final Class<FeatureCodec<T, ?>> codecClass = featureInput.getFeatureCodecClass();
        if (codecClass == null) {
            final Path featurePath = featureInput.toPath();
            IOUtils.assertFileIsReadable(featurePath);
            codec = (FeatureCodec<T, ?>) FeatureManager.getCodecForFile(featurePath, targetFeatureType);
            featureInput.setFeatureCodecClass((Class<FeatureCodec<T, ?>>) codec.getClass());
        } else {
            try {
                codec = codecClass.getDeclaredConstructor().newInstance();
                if ( !codec.canDecode(featureInput.toPath().toAbsolutePath().toUri().toString()) ) {
                    throw new GATKException(codec.getClass().getSimpleName() + " cannot decode " +
                            featureInput.toPath().toString());
                }
            } catch (final InstantiationException | IllegalAccessException | NoSuchMethodException | InvocationTargetException e) {
                throw new GATKException("Unable to automatically instantiate codec " + codecClass.getName());
            }
        }

        if (setNameOnCodec && codec instanceof NameAwareCodec) {
            final NameAwareCodec namedCodec = (NameAwareCodec) codec;
            namedCodec.setName(featureInput.getName());
        }

        return codec;
    }

    private static <T extends Feature> AbstractFeatureReader<T, ?> getTribbleFeatureReader(final FeatureInput<T> featureInput, final FeatureCodec<T, ?> codec, final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer) {
        Utils.nonNull(codec);
        try {
            // Must get the path to the data file from the codec here:
            final String absoluteRawPath = featureInput.getRawInputString();

            // Instruct the reader factory to not require an index. We will require one ourselves as soon as
            // a query by interval is attempted.
            final boolean requireIndex = false;

            // Only apply the wrappers if the feature input is in a remote location which will benefit from prefetching.
            final boolean prefetch = BucketUtils.isEligibleForPrefetching(featureInput);
            return AbstractFeatureReader.getFeatureReader(absoluteRawPath, null, codec, requireIndex,
                    BucketUtils.getPrefetchingWrapper(prefetch ? cloudPrefetchBuffer : 0),
                    BucketUtils.getPrefetchingWrapper(prefetch ? cloudIndexPrefetchBuffer : 0));
        } catch (final TribbleException e) {
            throw new GATKException("Error initializing feature reader for path " + featureInput.getFeaturePath(), e);
        }
    }

    protected static FeatureReader<VariantContext> getGenomicsDBFeatureReader(final GATKPath path, final File reference, final GenomicsDBOptions genomicsDBOptions) {
        final String workspace = IOUtils.getGenomicsDBAbsolutePath(path) ;
        if (workspace == null) {
            throw new IllegalArgumentException("Trying to create a GenomicsDBReader from  non-GenomicsDB input path " + path);
        } else if (Files.notExists(IOUtils.getPath(workspace.endsWith("/") ? workspace : workspace + "/"))) {
            throw new UserException("GenomicsDB workspace " + path + " does not exist");
        }

        final String callsetJson = IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME);
        final String vidmapJson = IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME);
        final String vcfHeader = IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME);

        IOUtils.assertPathsAreReadable(callsetJson, vidmapJson, vcfHeader);

        try {
            final GenomicsDBExportConfiguration.ExportConfiguration exportConfigurationBuilder =
                    createExportConfiguration(workspace, callsetJson, vidmapJson, vcfHeader, genomicsDBOptions);
            if (genomicsDBOptions.useBCFCodec()) {
                return new GenomicsDBFeatureReader<>(exportConfigurationBuilder, new BCF2Codec(), Optional.empty());
            } else {
                return new GenomicsDBFeatureReader<>(exportConfigurationBuilder, new VCFCodec(), Optional.empty());
            }
        } catch (final IOException e) {
            throw new UserException("Couldn't create GenomicsDBFeatureReader", e);
        }
    }

    /**
     * Returns the sequence dictionary for this source of Features.
     * Uses the dictionary from the VCF header (if present) for variant inputs,
     * otherwise attempts to create a sequence dictionary from the index file (if present).
     * Returns null if no dictionary could be created from either the header or the index.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        SAMSequenceDictionary dict = null;
        final Object header = getHeader();
        if ( header instanceof SVFeaturesHeader ) {
            dict = ((SVFeaturesHeader)header).getDictionary();
        } else if (header instanceof VCFHeader) {
            dict = ((VCFHeader) header).getSequenceDictionary();
        }
        if (dict != null && !dict.isEmpty()) {
            return dict;
        }
        if (hasIndex) {
            return IndexUtils.createSequenceDictionaryFromFeatureIndex(IOUtils.getPath(featureInput.getFeaturePath()));
        }
        return null;
    }

    /**
     * Restricts traversals of this data source via {@link #iterator} to only return Features that overlap the provided
     * intervals. Calls to {@link #query(SimpleInterval)} and/or {@link #queryAndPrefetch(Locatable)} are not
     * affected by these intervals.
     * <p>
     * Intervals MUST be non-overlapping and sorted in order of increasing start position, otherwise traversal
     * results will be incorrect.
     * <p>
     * Passing in a null or empty interval List clears the intervals for traversal, making future iterations
     * over this data source unrestricted by intervals.
     *
     * @param intervals Our next full traversal will return only Features overlapping these intervals
     */
    public void setIntervalsForTraversal(final List<SimpleInterval> intervals) {
        // Treat null and empty interval lists the same
        intervalsForTraversal = (intervals != null && !intervals.isEmpty()) ? intervals : null;

        if (intervalsForTraversal != null && !supportsRandomAccess) {
            throw new UserException("Input " + featureInput.getFeaturePath() + " must support random access to enable traversal by intervals. " +
                    "If it's a file, please index it using the bundled tool " + IndexFeatureFile.class.getSimpleName());
        }
    }


    /**
     * Gets an iterator over all Features in this data source, restricting traversal to Features
     * overlapping our intervals if intervals were provided via {@link #setIntervalsForTraversal(List)}
     * <p>
     * Calling this method invalidates (closes) any previous iterator obtained from this method.
     *
     * @return an iterator over all Features in this data source, limited to Features that overlap the intervals supplied via {@link #setIntervalsForTraversal(List)} (if intervals were provided)
     */
    @Override
    public Iterator<T> iterator() {
        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        try {
            // Save the iterator returned so that we can close it properly later
            currentIterator = intervalsForTraversal != null ? new FeatureIntervalIterator<>(intervalsForTraversal, featureReader, featureInput.getFeaturePath())
                    : featureReader.iterator();
            return currentIterator;
        } catch (final IOException e) {
            throw new GATKException("Error creating iterator over file " + featureInput.getFeaturePath(), e);
        }
    }

    /**
     * Gets an iterator over all Features in this data source that overlap the provided interval.
     * <p>
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     * <p>
     * Requires the backing file to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     * <p>
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     * <p>
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Features overlapping this interval
     * @return an iterator over all Features in this data source that overlap the provided interval
     */
    @Override
    public Iterator<T> query(final SimpleInterval interval) {
        return queryAndPrefetch(interval).iterator();
    }

    /**
     * Returns a List of all Features in this data source that overlap the provided interval.
     * <p>
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     * <p>
     * Requires the backing file to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     * <p>
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     * <p>
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Features overlapping this interval
     * @return a List of all Features in this data source that overlap the provided interval
     */
    public List<T> queryAndPrefetch(final Locatable interval) {
        if (!supportsRandomAccess) {
            throw new UserException("Input " + featureInput.getFeaturePath() + " must support random access to enable queries by interval. " +
                    "If it's a file, please index it using the bundled tool " + IndexFeatureFile.class.getSimpleName());
        }

        // If the query can be satisfied using existing cache contents, prepare for retrieval
        // by discarding all Features at the beginning of the cache that end before the start
        // of our query interval.
        if (queryCache.cacheHit(interval)) {
            queryCache.trimToNewStartPosition(interval.getStart());
        }
        // Otherwise, we have a cache miss, so go to disk to refill our cache.
        else {
            refillQueryCache(interval);
        }

        // Return the subset of our cache that overlaps our query interval
        return queryCache.getCachedFeaturesUpToStopPosition(interval.getEnd());
    }

    /**
     * Refill our cache from disk after a cache miss. Will prefetch Features overlapping an additional
     * queryLookaheadBases bases after the end of the provided interval, in addition to those overlapping
     * the interval itself.
     * <p>
     * Calling this has the side effect of invalidating (closing) any currently-open iteration over
     * this data source.
     *
     * @param interval the query interval that produced a cache miss
     */
    private void refillQueryCache(final Locatable interval) {
        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        // Expand the end of our query by the configured number of bases, in anticipation of probable future
        // queries with slightly larger start/stop positions.
        //
        // Note that it doesn't matter if we go off the end of the contig in the process, since
        // our reader's query operation is not aware of (and does not care about) contig boundaries.
        // Note: we use addExact to blow up on overflow rather than propagate negative results downstream
        final SimpleInterval queryInterval = new SimpleInterval(interval.getContig(), interval.getStart(), Math.addExact(interval.getEnd(), queryLookaheadBases));

        // Query iterator over our reader will be immediately closed after re-populating our cache
        try (final CloseableTribbleIterator<T> queryIter = featureReader.query(queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd())) {
            queryCache.fill(queryIter, queryInterval);
        } catch (final IOException e) {
            throw new GATKException("Error querying file " + featureInput + " over interval " + interval, e);
        }
    }

    /**
     * Get the logical name of this data source.
     *
     * @return the logical name of this data source
     */
    public String getName() {
        return featureInput.getName();
    }

    /**
     * Gets the header associated with this data source
     *
     * @return header associated with this data source as an Object
     */
    public Object getHeader() {
        return featureReader.getHeader();
    }

    /**
     * Permanently close this data source, invalidating any open iteration over it, and making it invalid for future
     * iterations and queries.
     */
    @Override
    public void close() {
        closeOpenIterationIfNecessary();

        logger.debug(String.format("Cache statistics for FeatureInput %s:", featureInput));
        queryCache.printCacheStatistics();

        try {
            if (featureReader != null) {
                featureReader.close();
            }
        } catch (final IOException e) {
            throw new GATKException("Error closing Feature reader for input " + featureInput);
        }
    }

    /**
     * Close the iterator currently open over this data source, if there is one.
     */
    private void closeOpenIterationIfNecessary() {
        if (currentIterator != null) {
            currentIterator.close();
            currentIterator = null;
        }
    }
}
