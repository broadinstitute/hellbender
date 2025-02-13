package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentDefinition;
import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.config.GATKConfig;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;


/**
 * Handles discovery of available codecs and Feature arguments, file format detection and codec selection,
 * and creation/management/querying of FeatureDataSources for each source of Features.
 *
 * At startup, walks the packages specified in {@link GATKConfig#codec_packages} in the config file to discover what codecs are available
 * to decode Feature-containing files.
 *
 * Then, given a tool instance, it discovers what FeatureInput argument fields are declared in the
 * tool's class hierarchy (and associated ArgumentCollections), and for each argument actually specified
 * by the user on the command line, determines the type of the file and the codec required to decode it,
 * creates a FeatureDataSource for that file, and adds it to a query-able resource pool.
 *
 * Clients can then call {@link #getFeatures(FeatureInput, Locatable)} to query the data source for
 * a particular FeatureInput over a specific interval.
 */
public final class FeatureManager implements AutoCloseable {
    private static final Logger logger = LogManager.getLogger(FeatureManager.class);

    /**
     * All codecs descend from this class
     */
    private static final Class<FeatureCodec> CODEC_BASE_CLASS = FeatureCodec.class;

    /**
     * The codec classes we locate when searching codec packages
     */
    private static final Set<Class<?>> DISCOVERED_CODECS;

    /**
     * Feature arguments in tools are of this type
     */
    private static final Class<FeatureInput> FEATURE_ARGUMENT_CLASS = FeatureInput.class;

    /*
     * At startup, walk through the packages in codec packages, and save any (concrete) FeatureCodecs discovered
     * in DISCOVERED_CODECS
     */
    static {

        // Get our configuration:
        final GATKConfig config = ConfigFactory.getInstance().getGATKConfig();

        final ClassFinder finder = new ClassFinder();
        for ( final String codecPackage : config.codec_packages() ) {
            finder.find(codecPackage, CODEC_BASE_CLASS);
        }
        // Exclude abstract classes and interfaces from the list of discovered codec classes
        DISCOVERED_CODECS = Collections.unmodifiableSet(finder.getConcreteClasses());
    }

    /**
     * The simple class name of the tool instance containing the FeatureInput argument values that will form the basis of our
     * pool of FeatureDataSources
     */
    private final String toolInstanceSimpleClassName;

    /**
     * Mapping from FeatureInput argument to query-able FeatureDataSource for that source of Features
     */
    private final Map<FeatureInput<? extends Feature>, FeatureDataSource<? extends Feature>> featureSources;

    /**
     * Create a FeatureManager given a CommandLineProgram tool instance, discovering all FeatureInput
     * arguments in the tool and creating query-able FeatureDataSources for them. Uses the default
     * caching behavior of {@link FeatureDataSource}.
     *
     * @param toolInstance Instance of the tool to be run (potentially containing one or more FeatureInput arguments)
     *                     Must have undergone command-line argument parsing and argument value injection already.
     */
    public FeatureManager( final CommandLineProgram toolInstance ) {
        this(toolInstance, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    /**
     * Create a FeatureManager given a CommandLineProgram tool instance, discovering all FeatureInput
     * arguments in the tool and creating query-able FeatureDataSources for them. Allows control over
     * how much caching is performed by each {@link FeatureDataSource}.
     *
     * @param toolInstance Instance of the tool to be run (potentially containing one or more FeatureInput arguments)
     *                     Must have undergone command-line argument parsing and argument value injection already.
     * @param featureQueryLookahead When querying FeatureDataSources, cache this many extra bases of context beyond
     *                              the end of query intervals in anticipation of future queries (>= 0).
     */
    public FeatureManager( final CommandLineProgram toolInstance, final int featureQueryLookahead ) {
        this(toolInstance, featureQueryLookahead, 0, 0);
    }


    /**
     * Create a FeatureManager given a CommandLineProgram tool instance, discovering all FeatureInput
     * arguments in the tool and creating query-able FeatureDataSources for them. Allows control over
     * how much caching is performed by each {@link FeatureDataSource}.
     *
     * @param toolInstance Instance of the tool to be run (potentially containing one or more FeatureInput arguments)
     *                     Must have undergone command-line argument parsing and argument value injection already.
     * @param featureQueryLookahead When querying FeatureDataSources, cache this many extra bases of context beyond
     *                              the end of query intervals in anticipation of future queries (>= 0).
     * @param cloudPrefetchBuffer MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     *
     */
    public FeatureManager(final CommandLineProgram toolInstance, final int featureQueryLookahead, final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer) {
        this(toolInstance, featureQueryLookahead, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, null);
    }

    /**
     * Create a FeatureManager given a CommandLineProgram tool instance, discovering all FeatureInput
     * arguments in the tool and creating query-able FeatureDataSources for them. Allows control over
     * how much caching is performed by each {@link FeatureDataSource}.
     *  @param toolInstance Instance of the tool to be run (potentially containing one or more FeatureInput arguments)
     *                     Must have undergone command-line argument parsing and argument value injection already.
     * @param featureQueryLookahead When querying FeatureDataSources, cache this many extra bases of context beyond
     *                              the end of query intervals in anticipation of future queries (>= 0).
     * @param cloudPrefetchBuffer MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param gdbOptions settings for GenomicsDB to use when reading from a GenomicsDB workspace
     *
     */
    public FeatureManager(final CommandLineProgram toolInstance, final int featureQueryLookahead, final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final GenomicsDBOptions gdbOptions) {
        this.toolInstanceSimpleClassName = toolInstance.getClass().getSimpleName();
        this.featureSources = new LinkedHashMap<>();

        initializeFeatureSources(featureQueryLookahead, toolInstance, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, gdbOptions);
    }

    /**
     * Same as {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, GenomicsDBOptions)}, except used when the
     *  FeatureInputs (and associated types) are known.
     *
     *  This constructor should only be used in test code.
     *
     * @param featureInputsToTypeMap {@link Map} of a {@link FeatureInput} to the output type that must extend {@link Feature}.  Never {@code null}
     * @param toolInstanceName See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, GenomicsDBOptions)}
     * @param featureQueryLookahead See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, GenomicsDBOptions)}
     * @param cloudPrefetchBuffer See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, GenomicsDBOptions)}
     * @param cloudIndexPrefetchBuffer See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, GenomicsDBOptions)}
     * @param reference See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, GenomicsDBOptions)}
     */
    @VisibleForTesting
    FeatureManager(final Map<FeatureInput<? extends Feature>, Class<? extends Feature>> featureInputsToTypeMap, final String toolInstanceName, final int featureQueryLookahead, final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final Path reference) {

        Utils.nonNull(featureInputsToTypeMap);

        this.toolInstanceSimpleClassName = toolInstanceName;
        this.featureSources = new LinkedHashMap<>();
        Utils.nonNull(featureInputsToTypeMap);
        featureInputsToTypeMap.forEach((k,v) -> addToFeatureSources(featureQueryLookahead, k, v, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, reference));
    }

    /**
     * Given our tool instance, discover all argument of type FeatureInput (or Collections thereof), determine
     * the type of each Feature-containing file, and add a FeatureDataSource for each file to our query pool.
     *
     * @param featureQueryLookahead Set up each FeatureDataSource to cache this many extra bases of context beyond
     *                              the end of query intervals in anticipation of future queries (>= 0).
     * @param toolInstance Instance of the tool to be run (potentially containing one or more FeatureInput arguments)
     *                     Must have undergone command-line argument parsing and argument value injection already.
     * @param cloudPrefetchBuffer MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    private void initializeFeatureSources( final int featureQueryLookahead, final CommandLineProgram toolInstance, final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final GenomicsDBOptions gdbOptions) {

        // Discover all arguments of type FeatureInput (or Collections thereof) in our tool's class hierarchy
        // (and associated ArgumentCollections). Arguments not specified by the user on the command line will
        // come back to us with a null FeatureInput.
        final List<Pair<ArgumentDefinition, FeatureInput>> featureArgumentValues =
                toolInstance.getCommandLineParser().gatherArgumentValuesOfType(FEATURE_ARGUMENT_CLASS);

        for ( final Pair<ArgumentDefinition, FeatureInput> featureArgument : featureArgumentValues ) {
            final FeatureInput<? extends Feature> featureInput = featureArgument.getValue();

            // Only create a data source for Feature arguments that were actually specified
            if ( featureInput != null ) {
                final Class<? extends Feature> featureType = getFeatureTypeForFeatureInputArgument(featureArgument.getKey());
                addToFeatureSources(featureQueryLookahead, featureInput, featureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                        gdbOptions);
            }
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    public void dumpAllFeatureCacheStats() {
        for ( final FeatureDataSource f : featureSources.values() ) {
            f.printCacheStats();
        }
    }

    void addToFeatureSources(final int featureQueryLookahead, final FeatureInput<? extends Feature> featureInput,
                             final Class<? extends Feature> featureType, final int cloudPrefetchBuffer,
                             final int cloudIndexPrefetchBuffer, final Path reference) {
        // Create a new FeatureDataSource for this file, and add it to our query pool
        featureSources.put(featureInput, new FeatureDataSource<>(featureInput, featureQueryLookahead, featureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, new GenomicsDBOptions(reference)));
    }

    /**
     * Add the feature data source to the given feature input.
     *
     * @param featureQueryLookahead look ahead this many bases during queries that produce cache misses
     * @param featureInput source of features
     * @param featureType class of features
     * @param cloudPrefetchBuffer MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param genomicsDBOptions options and info for reading from a GenomicsDB
     *
     * Note: package-visible to enable access from the core walker classes
     * (but not actual tools, so it's not protected).
     */
    void addToFeatureSources(final int featureQueryLookahead, final FeatureInput<? extends Feature> featureInput,
                             final Class<? extends Feature> featureType, final int cloudPrefetchBuffer,
                             final int cloudIndexPrefetchBuffer, final GenomicsDBOptions genomicsDBOptions) {
        // Create a new FeatureDataSource for this file, and add it to our query pool
        featureSources.put(featureInput, new FeatureDataSource<>(featureInput, featureQueryLookahead, featureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, genomicsDBOptions));
    }

    <F extends Feature> void addToFeatureSources (final FeatureInput<F> featureInput,
                                                  final FeatureDataSource<F> featureDataSource) {
        featureSources.put(featureInput, featureDataSource);
    }

    /**
     * Given a ArgumentDefinition for an argument known to be of type FeatureInput (or a Collection thereof), retrieves the type
     * parameter for the FeatureInput (eg., for FeatureInput<VariantContext> or List<FeatureInput<VariantContext>>
     * this would be VariantContext).
     *
     * @param argDef an {@code ArgumentDefinition} for an argument known to be of type FeatureInput whose type parameter to retrieve
     * @return type parameter of the FeatureInput declaration represented by the given ArgumentDefinition
     */
    @SuppressWarnings("unchecked")
    static Class<? extends Feature> getFeatureTypeForFeatureInputArgument( final ArgumentDefinition argDef ) {
        final Type featureInputType = argDef.isCollection() ?
                                getNextTypeParameter((ParameterizedType)(argDef.getUnderlyingField().getGenericType())) :
                                argDef.getUnderlyingField().getGenericType();

        if ( ! (featureInputType instanceof ParameterizedType) ) {
            throw new GATKException(String.format("FeatureInput declaration for argument --%s lacks an explicit type parameter for the Feature type",
                                argDef.getUnderlyingField().getAnnotation(Argument.class).fullName()));
        }

        return (Class<? extends Feature>)getNextTypeParameter((ParameterizedType)featureInputType);
    }

    /**
     * Helper method for {@link #getFeatureTypeForFeatureInputArgument(ArgumentDefinition)} that "unpacks" a
     * parameterized type by one level of parameterization. Eg., given List<FeatureInput<VariantContext>>
     * would return FeatureInput<VariantContext>.
     *
     * @param parameterizedType parameterized type to unpack
     * @return the type parameter of the given parameterized type
     */
    private static Type getNextTypeParameter( final ParameterizedType parameterizedType ) {
        final Type[] typeParameters = parameterizedType.getActualTypeArguments();
        if ( typeParameters.length != 1 ) {
            throw new GATKException("Found a FeatureInput declaration with multiple type parameters, which is not supported");
        }
        return typeParameters[0];
    }

    /**
     * Does this manager have no sources of Features to query?
     *
     * @return true if there are no Feature sources available to query, otherwise false
     */
    public boolean isEmpty() {
        return featureSources.isEmpty();
    }


    /**
     * This method finds and returns all of the variant headers from the feature sources.
     *
     * @return A list of all variant headers for features.
     */
    public List<VCFHeader> getAllVariantHeaders() {
        return featureSources.values().stream()
                .map(feature -> feature.getHeader())
                .filter(header -> header instanceof VCFHeader)
                .map(header -> (VCFHeader)header).collect(Collectors.toList());
    }

    /**
     * Returns the list of sequence dictionaries retrieved from the VCF headers of variant Feature inputs.
     * Note: this method returns an empty list if the variant inputs
     * happen not to have sequence dictionaries (since they are optional in the VCF format).
     */
    public List<SAMSequenceDictionary> getVariantSequenceDictionaries() {
        return getAllVariantHeaders()
                .stream().map(h -> h.getSequenceDictionary())
                .filter(dict -> dict != null)
                .collect(Collectors.toList());
    }

    /**
     * Returns the sequence dictionaries associated with all feature sources.
     * This method will return an empty List if none of the feature sources have dictionaries.
     */
    public List<SAMSequenceDictionary> getAllSequenceDictionaries() {
        return featureSources.values().stream().map(fs -> fs.getSequenceDictionary())
                .filter(dict -> dict != null)
                .collect(Collectors.toList());
    }

    /**
     * Get all the sources managed by this FeatureManager.
     */
    public Set<FeatureInput<? extends Feature>> getAllInputs() {
        return featureSources.keySet();
    }

    /**
     * Given a FeatureInput argument field from our tool, queries the data source for that FeatureInput
     * over the specified interval, and returns a List of the Features overlapping that interval from
     * that data source.
     *
     * Will throw an exception if the provided FeatureInput did not come from the tool that this
     * FeatureManager was initialized with, or was not an @Argument-annotated field in the tool
     * (or parent classes).
     *
     * @param featureDescriptor FeatureInput argument from our tool representing the Feature source to query
     * @param interval interval to query over (returned Features will overlap this interval)
     * @param <T> type of Feature in the source represented by featureDescriptor
     * @return A List of all Features in the backing data source for the provided FeatureInput that overlap
     *         the provided interval (may be empty if there are none, but never null)
     */
    public <T extends Feature> List<T> getFeatures( final FeatureInput<T> featureDescriptor, final Locatable interval ) {
        final FeatureDataSource<T> dataSource = lookupDataSource(featureDescriptor);

        // No danger of a ClassCastException here, since we verified that the FeatureDataSource for this
        // FeatureInput will return Features of the expected type T when we first created the data source
        // in initializeFeatureSources()
        return dataSource.queryAndPrefetch(interval);
    }

    /**
     * Given a FeatureInput argument field from our tool, returns an iterator to its features starting
     * from the first one.
     * <p><b>Warning!</b>: calling this method a second time on the same {@link FeatureInput}
     * on the same FeatureManager instance will invalidate (close) the iterator returned from
     * the first call.
     * </p>
     * <p>
     * An exception will be thrown if the {@link FeatureInput} provided did not come from the tool that this
     * manager was initialized with, or was not an &#64;Argument-annotated field in the tool
     * (or parent classes).
     * </p>
     *
     * @param featureDescriptor FeatureInput argument from our tool representing the Feature source to query
     * @param <T> type of Feature in the source represented by featureDescriptor
     * @return never {@code null}, a iterator to all the features in the backing data source.
     * @throws GATKException if the feature-descriptor is not found in the manager or is {@code null}.
     */
    public <T extends Feature> Iterator<T> getFeatureIterator(final FeatureInput<T> featureDescriptor) {
        final FeatureDataSource<T> dataSource = lookupDataSource(featureDescriptor);
        return dataSource.iterator();
    }

    /**
     * As above, but takes an optional list of intervals to examine.
     * @param featureDescriptor FeatureInput to scan
     * @param intervals The userIntervals to examine (may be null)
     * @param <T> Feature type
     * @return An iterator over the Features
     */
    public <T extends Feature> Iterator<T> getFeatureIterator( final FeatureInput<T> featureDescriptor,
                                                               final List<SimpleInterval> intervals ) {
        final FeatureDataSource<T> dataSource = lookupDataSource(featureDescriptor);
        dataSource.setIntervalsForTraversal(intervals);
        final Iterator<T> itr = dataSource.iterator();
        dataSource.setIntervalsForTraversal(null);
        return itr;
    }

    /**
     * Get the header associated with a particular FeatureInput
     *
     * @param featureDescriptor the FeatureInput whose header we want to retrieve
     * @param <T> type of Feature in our FeatureInput
     * @return header for the provided FeatureInput
     */
    public <T extends Feature> Object getHeader( final FeatureInput<T> featureDescriptor ) {
        final FeatureDataSource<T> dataSource = lookupDataSource(featureDescriptor);
        return dataSource.getHeader();
    }

    /**
     * Retrieve the data source for a particular FeatureInput. Throws an exception if the provided
     * FeatureInput is not among our discovered sources of Features.
     *
     * @param featureDescriptor FeatureInput whose data source to retrieve
     * @param <T> type of Feature in our FeatureInput
     * @return query-able data source for the provided FeatureInput, if it was found
     */
    private <T extends Feature> FeatureDataSource<T> lookupDataSource( final FeatureInput<T> featureDescriptor ) {
        @SuppressWarnings("unchecked") final FeatureDataSource<T> dataSource = (FeatureDataSource<T>)featureSources.get(featureDescriptor);

        // Make sure the provided FeatureInput actually came from our tool as an @Argument-annotated field
        if ( dataSource == null ) {
            throw new GATKException(String.format("FeatureInput %s not found in feature manager's database for tool %s. " +
                                                  "In order to be detected, FeatureInputs must be declared in the tool class " +
                                                  "itself, a superclass of the tool class, or an @ArgumentCollection declared " +
                                                  "in the tool class or a superclass. They must also be annotated as an @Argument.",
                                                  featureDescriptor.getName(), toolInstanceSimpleClassName));
        }

        return dataSource;
    }

    /**
     * Utility method that determines the correct codec to use to read Features from the provided file.
     *
     * Codecs MUST correctly implement the {@link FeatureCodec#canDecode(String)} method
     * in order to be considered as candidates for decoding the file.
     *
     * Throws an exception if no suitable codecs are found (this is a user error, since the file is of
     * an unsupported format), or if more than one codec claims to be able to decode the file (this is
     * a configuration error on the codec authors' part).
     *
     * @param featurePath path for which to find the right codec
     * @return the codec suitable for decoding the provided file
     */
    public static FeatureCodec<? extends Feature, ?> getCodecForFile( final Path featurePath ) {
        return getCodecForFile(featurePath, Feature.class);
    }

    /**
     * Utility method that determines the correct codec to use to read Features from the provided file,
     * optionally considering only codecs that produce a particular type of Feature.
     *
     * Codecs MUST correctly implement the {@link FeatureCodec#canDecode(String)} method
     * in order to be considered as candidates for decoding the file, and must produce
     * Features of the specified type if featureType is non-null.
     *
     * Throws an exception if no suitable codecs are found (this is a user error, since the file is of
     * an unsupported format), or if more than one codec claims to be able to decode the file (this is
     * a configuration error on the codec authors' part).
     *
     * @param featurePath Path for which to find the right codec
     * @param featureType If specified, consider only codecs that produce Features of this type. May be null,
     *                    in which case all codecs are considered.
     * @return the codec suitable for decoding the provided file
     */
    public static FeatureCodec<? extends Feature, ?> getCodecForFile( final Path featurePath, final Class<? extends Feature> featureType ) {
        // Make sure Path exists/is readable
        if ( ! Files.isReadable(featurePath) ) {
            throw new UserException.CouldNotReadInputFile(featurePath.toUri().toString());
        }

        // Gather all discovered codecs that produce the right feature subtype and claim to be able
        // to decode the given file according to their canDecode() methods
        final List<FeatureCodec<? extends Feature, ?>> candidateCodecs = getCandidateCodecsForFile(featurePath, featureType);

        // If no codecs can handle the file, it's a user error (the user provided a file in an unsupported format)
        if ( candidateCodecs.isEmpty() ) {
            throw new UserException.NoSuitableCodecs(featurePath);
        }

        // If we still have multiple candidate codecs, it's a configuration error on the part of the codec authors
        if ( candidateCodecs.size() > 1 ) {
            final StringBuilder multiCodecMatches = new StringBuilder();
            for ( final FeatureCodec<? extends Feature, ?> candidateCodec : candidateCodecs ) {
                multiCodecMatches.append(candidateCodec.getClass().getCanonicalName());
                multiCodecMatches.append(' ');
            }
            throw new GATKException("Multiple codecs found able to decode file " + featurePath.toAbsolutePath().toUri() +
                                    ". This indicates a misconfiguration on the part of the codec authors. " +
                                    "Matching codecs are: " + multiCodecMatches.toString());
        }

        final FeatureCodec<? extends Feature, ?> selectedCodec = candidateCodecs.get(0);
        logger.info("Using codec " + selectedCodec.getClass().getSimpleName() + " to read file " + featurePath.toAbsolutePath().toUri());
        return selectedCodec;
    }

    /**
     * Returns a List of all codecs in DISCOVERED_CODECS that return a feature type assignable to the
     * requested feature subtype according to the {@link FeatureCodec#getFeatureType()} method,
     * and that claim to be able to decode the specified file according to the {@link FeatureCodec#canDecode(String)} methods.
     *
     * @param featureFile file for which to find potential codecs
     * @param featureType subtype of Feature for which codecs are sought
     * @return A List of all suitable codecs in DISCOVERED_CODECS
     */
    private static List<FeatureCodec<? extends Feature, ?>> getCandidateCodecsForFile( final Path featureFile, final Class<? extends Feature> featureType )  {
        final List<FeatureCodec<? extends Feature, ?>> candidateCodecs = new ArrayList<>();

        final String featurePathname = featureFile.toAbsolutePath().toUri().toString();
        for ( final Class<?> codecClass : DISCOVERED_CODECS ) {
            try {
                final FeatureCodec<? extends Feature, ?> codec =
                        (FeatureCodec<? extends Feature, ?>)codecClass.getDeclaredConstructor().newInstance();
                if ( featureType.isAssignableFrom(codec.getFeatureType()) ) {
                    if ( codec.canDecode(featurePathname) ) {
                        candidateCodecs.add(codec);
                    }
                }
            }
            catch (InstantiationException | IllegalAccessException | NoSuchMethodException | InvocationTargetException e ) {
                throw new GATKException("Unable to automatically instantiate codec " + codecClass.getName(), e);
            }
        }

        return candidateCodecs;
    }

    /**
     * @param file file to check
     * @return True if the file exists and contains Features (ie., we have a FeatureCodec that can decode it), otherwise false
     */
    public static boolean isFeatureFile( final Path file ) {
        return Files.exists(file) && ! getCandidateCodecsForFile(file, Feature.class).isEmpty();
    }

    /**
     * Permanently closes this manager by closing all backing data sources
     */
    @Override
    public void close() {
        featureSources.values().forEach(ds -> ds.close());
    }

}
