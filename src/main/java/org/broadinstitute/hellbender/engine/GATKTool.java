package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.nio.file.Path;
import java.time.ZonedDateTime;
import java.util.*;
import java.util.stream.Stream;


import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.config.GATKConfig;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.ShardingVCFWriter;
import org.broadinstitute.hellbender.utils.variant.writers.IntervalFilteringVcfWriter;

/**
 * Base class for all GATK tools. Tool authors that want to write a "GATK" tool but not use one of
 * the pre-packaged Walker traversals should feel free to extend this class directly. All other
 * GATK tools should extend one of the Walker classes instead.
 */
public abstract class GATKTool extends CommandLineProgram {

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = requiresIntervals() ? new RequiredIntervalArgumentCollection() : new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    protected final ReadInputArgumentCollection readArguments = requiresReads() ? new RequiredReadInputArgumentCollection() : new OptionalReadInputArgumentCollection();

    @ArgumentCollection
    protected final ReferenceInputArgumentCollection referenceArguments = requiresReference() ? new RequiredReferenceInputArgumentCollection() :  new OptionalReferenceInputArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME,
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME,
            doc = "Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.", optional = true, common = true)
    private GATKPath masterSequenceDictionaryFilename = null;

    public static final String SECONDS_BETWEEN_PROGRESS_UPDATES_NAME = "seconds-between-progress-updates";

    @Argument(fullName = SECONDS_BETWEEN_PROGRESS_UPDATES_NAME, shortName = SECONDS_BETWEEN_PROGRESS_UPDATES_NAME, doc = "Output traversal statistics every time this many seconds elapse", optional = true, common = true)
    private double secondsBetweenProgressUpdates = ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES;

    @ArgumentCollection
    protected SequenceDictionaryValidationArgumentCollection seqValidationArguments = getSequenceDictionaryValidationArgumentCollection();

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_SHORT_NAME,
            doc = "If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.", optional=true, common = true)
    public boolean createOutputBamIndex = ConfigFactory.getInstance().getGATKConfig().createOutputBamIndex();

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_SHORT_NAME,
            doc = "If true, create a MD5 digest for any BAM/SAM/CRAM file created", optional=true, common = true)
    public boolean createOutputBamMD5 = false;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_SHORT_NAME,
            doc = "If true, create a VCF index when writing a coordinate-sorted VCF file.", optional=true, common = true)
    public boolean createOutputVariantIndex = true;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_SHORT_NAME,
            doc = "If true, create a a MD5 digest any VCF file created.", optional=true, common = true)
    public boolean createOutputVariantMD5 = false;

    @Argument(fullName = StandardArgumentDefinitions.MAX_VARIANTS_PER_SHARD_LONG_NAME, optional = true, minValue = 0, common = true,
            doc = "If non-zero, partitions VCF output into shards, each containing up to the given number of records.")
    private int maxVariantsPerShard = 0;

    @Argument(fullName= StandardArgumentDefinitions.LENIENT_LONG_NAME,
            shortName = StandardArgumentDefinitions.LENIENT_SHORT_NAME,
            doc = "Lenient processing of VCF files", common = true, optional = true)
    protected boolean lenientVCFProcessing = false;

    @Argument(fullName = StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD, shortName = StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD, doc = "If true, adds a PG tag to created SAM/BAM/CRAM files.", optional=true, common = true)
    public boolean addOutputSAMProgramRecord = true;

    @Argument(fullName = StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, shortName = StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, doc = "If true, adds a command line header line to created VCF files.", optional=true, common = true)
    public boolean addOutputVCFCommandLine = true;

    // default value of 40MB based on a test with CountReads (it's 5x faster than no prefetching)
    @Argument(fullName = StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, shortName = StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_SHORT_NAME, doc = "Size of the cloud-only prefetch buffer (in MB; 0 to disable).", optional=true)
    public int cloudPrefetchBuffer = getDefaultCloudPrefetchBufferSize();

    @Argument(fullName = StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, shortName = StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_SHORT_NAME, doc = "Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.", optional=true)
    public int cloudIndexPrefetchBuffer = getDefaultCloudIndexPrefetchBufferSize();

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_BAM_INDEX_CACHING_LONG_NAME,
            shortName = StandardArgumentDefinitions.DISABLE_BAM_INDEX_CACHING_SHORT_NAME,
            doc = "If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.",
            optional = true)
    public boolean disableBamIndexCaching = false;

    @Argument(fullName = StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
            doc = "If true, don't emit genotype fields when writing vcf file output.", optional = true)
    public boolean outputSitesOnlyVCFs = false;

    /**
     * Master sequence dictionary to be used instead of all other dictionaries (if provided).
     */
    private SAMSequenceDictionary masterSequenceDictionary = null;

    /*
     * TODO: Feature arguments for the current tool are currently discovered through reflection via FeatureManager.
     * TODO: Perhaps we should eventually do the same auto-discovery for all input arguments (reads, reference, etc.)
     */

    /*
     * Engine-wide data structures. Package-private so that the various *Walker classes in the engine can access
     * them, but concrete tool implementations that extend from us cannot.
     */

    /**
     * Our source of reference data (null if no reference was provided)
     */
    ReferenceDataSource reference;

    /**
     * Our source of reads data (null if no source of reads was provided)
     */
    ReadsDataSource reads;

    /**
     * Our source of Feature data (null if no source of Features was provided)
     */
    public FeatureManager features;

    /**
     * Intervals to be used for traversal (null if no intervals were provided).
     *
     * Walker base classes (ReadWalker, etc.) are responsible for hooking these intervals up to
     * their particular driving data source.
     */
    List<SimpleInterval> userIntervals;

    /**
     * Get the {@link ReferenceDataSource} for this {@link GATKTool}.
     * Will throw a {@link GATKException} if the reference is null.
     * Clients are expected to call the {@link #hasReference()} method prior to calling this.
     *
     * Should only be called by walker base classes in the engine (such as {@link ReadWalker}), or by "free-form" tools that
     * extend the {@link GATKTool} class directly rather than one of the built-in walker types.
     * Tools that extend a walker type should get their data via {@code apply()} rather than directly accessing
     * the engine datasources.
     *
     * @return the {@link ReferenceDataSource} for this {@link GATKTool}.  Never {@code null}.
     */
    protected ReferenceDataSource directlyAccessEngineReferenceDataSource() {
        if ( reference == null ) {
            throw new GATKException("Attempted to retrieve null reference!");
        }
        return reference;
    }

    /**
     * Get the {@link ReadsDataSource} for this {@link GATKTool}.
     * Will throw a {@link GATKException} if the reads are null.
     * Clients are expected to call the {@link #hasReads()} method prior to calling this.
     *
     * Should only be called by walker base classes in the engine (such as {@link ReadWalker}), or by "free-form" tools that
     * extend the {@link GATKTool} class directly rather than one of the built-in walker types.
     * Tools that extend a walker type should get their data via {@code apply()} rather than directly accessing
     * the engine datasources.
     *
     * @return the {@link ReadsDataSource} for this {@link GATKTool}.  Never {@code null}.
     */
    protected ReadsDataSource directlyAccessEngineReadsDataSource() {
        if ( reads == null ) {
            throw new GATKException("Attempted to retrieve null reads!");
        }
        return reads;
    }

    /**
     * Get the {@link FeatureManager} for this {@link GATKTool}.
     * Will throw a {@link GATKException} if the features are null.
     * Clients are expected to call the {@link #hasFeatures()} method prior to calling this.
     *
     * Should only be called by walker base classes in the engine (such as {@link ReadWalker}), or by "free-form" tools that
     * extend the {@link GATKTool} class directly rather than one of the built-in walker types.
     * Tools that extend a walker type should get their data via {@code apply()} rather than directly accessing
     * the engine datasources.
     *
     * @return the {@link FeatureManager} for this {@link GATKTool}.  Never {@code null}.
     */
    protected FeatureManager directlyAccessEngineFeatureManager() {
        if ( features == null ) {
            throw new GATKException("Attempted to retrieve null features!");
        }
        return features;
    }

    /**
     * Progress meter to print out traversal statistics. Subclasses must invoke
     * {@link ProgressMeter#update(Locatable)} after each record processed from
     * the primary input in their {@link #traverse} method.
     */
    protected ProgressMeter progressMeter;

    /**
     * Return the list of GATKCommandLinePluginDescriptors to be used for this tool.
     * Uses the read filter plugin.
     */
    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        GATKReadFilterPluginDescriptor readFilterDescriptor = new GATKReadFilterPluginDescriptor(getDefaultReadFilters());
        return useVariantAnnotations()?
                Arrays.asList(readFilterDescriptor, new GATKAnnotationPluginDescriptor(
                        getDefaultVariantAnnotations(), getDefaultVariantAnnotationGroups())):
                Collections.singletonList(readFilterDescriptor);
    }

    /**
     * Returns the default list of ReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation returns an empty list. Subclasses can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {@link #getHeaderForReads}. The actual SAMFileHeader is propagated to the read filters
     * by {@link #makeReadFilter} after the filters have been merged with command line arguments.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new WellformedReadFilter());
    }

    /**
     * Returns a read filter (simple or composite) that can be applied to reads. This implementation combines
     * the default read filters for this tool (returned by {@link #getDefaultReadFilters} along with any read filter
     * command line directives specified by the user (such as enabling other filters or disabling default filters);
     * wraps each filter in the resulting list with a CountingReadFilter; and returns a single composite filter
     * resulting from the list by and'ing them together.
     *
     * NOTE: Most tools will not need to override the method, and should only do so in order to provide custom
     * behavior or processing of the final merged read filter. To change the default read filters used by the tool,
     * override {@link #getDefaultReadFilters} instead.
     *
     * Implementations of {@link #traverse()} should call this method once before iterating over the reads, in order to
     * unnecessary avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
     public CountingReadFilter makeReadFilter(){
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return hasReads() ?
                readFilterPlugin.getMergedCountingReadFilter(getHeaderForReads()) :
                new CountingReadFilter(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    /**
     * Must be overridden in order to add annotation arguments to the engine. If this is set to true the engine will
     * dynamically discover all {@link Annotation}s in the packages defined by {@link GATKAnnotationPluginDescriptor#getPackageNames()} and automatically
     * generate and add command line arguments allowing the user to specify which annotations or groups of annotations to use.
     *
     * To specify default annotations for a tool simply specify them using {@link #getDefaultVariantAnnotationGroups()} or {@link #getDefaultVariantAnnotations()}
     *
     * To access instantiated annotation objects simply use {@link #makeVariantAnnotations()}.
     */
    public boolean useVariantAnnotations() {
        return false;
    }

    /**
     * Returns the default list of {@link Annotation}s that are used for this tool. The annotations returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation returns an empty list. Subclasses can override to provide alternative annotations.
     *
     * @return List of individual annotations to be applied for this tool.
     */
    public List<Annotation> getDefaultVariantAnnotations() {
        return Collections.emptyList();
    }

    /**
     * Returns the default list of annotation groups that are used for this tool. The annotations returned
     * by this method will have default arguments, which can be overridden with specific arguments using
     * {@link #getDefaultVariantAnnotations()}. Returned annotation groups are subject to selective enabling/disabling
     * by the user via the command line. The default implementation returns an empty list.
     *
     * @return List of annotation groups to be applied for this tool.
     */
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Collections.emptyList();
    }

    /**
     * Returns a list of annotations that can be applied to VariantContexts. This implementation combines
     * the default annotations for this tool (returned by {@link #getDefaultVariantAnnotations()} and {@link #getDefaultVariantAnnotationGroups()})
     * along with any annotations command line directives specified by the user (such as enabling other annotations/groups
     * or disabling default annotations) and returns a collection of all the annotation arguments instantiated.
     *
     * NOTE: Most tools will not need to override the method, and should only do so in order to provide custom
     * behavior or processing of the final annotations based on other command line input. To change the default
     * annotations used by the tool, override {@link #getDefaultVariantAnnotations()} instead.
     *
     * To apply returned annotations to a VariantContext, simply use a {@link org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine}
     * constructed with the discovered annotations.
     */
    public Collection<Annotation> makeVariantAnnotations(){
        if (!useVariantAnnotations()) {
            throw new GATKException("Tool requested variant annotations but has not overridden 'useVariantAnnotations()' to return true");
        }
        final GATKAnnotationPluginDescriptor annotationPlugin =
                getCommandLineParser().getPluginDescriptor(GATKAnnotationPluginDescriptor.class);
        return annotationPlugin.getResolvedInstances();
    }

    /**
     * Returns the pre-filter read transformer (simple or composite) that will be applied to the reads before filtering.
     * The default implementation uses the {@link ReadTransformer#identity()}.
     * Default implementation of {@link #traverse()} calls this method once before iterating over the reads and reuses
     * the transformer object to avoid object allocation.
     *
     * Subclasses can extend to provide own transformers (ie override and call super).
     * Multiple transformers can be composed by using {@link ReadTransformer} composition methods.
     */
    public ReadTransformer makePreReadFilterTransformer() {
        return ReadTransformer.identity();
    }

    /**
     * Returns the post-filter read transformer (simple or composite) that will be applied to the reads after filtering.
     * The default implementation uses the {@link ReadTransformer#identity()}.
     * Default implementation of {@link #traverse()} calls this method once before iterating over the reads and reuses
     * the transformer object to avoid object allocation.
     *
     * Subclasses can extend to provide own transformers (ie override and call super).
     * Multiple transformers can be composed by using {@link ReadTransformer} composition methods.
     */
    public ReadTransformer makePostReadFilterTransformer(){
        return ReadTransformer.identity();
    }


    /**
     * Returns a stream over the reads, which are:
     *
     * 1. Transformed with {@link #makePreReadFilterTransformer()}.
     * 2. Filtered with {@code filter}.
     * 3. Transformed with {@link #makePostReadFilterTransformer()}.
     *
     * Note: the filter is passed to keep the state of {@link CountingReadFilter}, obtained with {@link #makeReadFilter()}.
     */
    protected Stream<GATKRead> getTransformedReadStream(final ReadFilter filter) {
        // if has reads, return an transformed/filtered/transformed stream
        if (hasReads()) {
            final ReadTransformer preTransformer = makePreReadFilterTransformer();
            final ReadTransformer postTransformer = makePostReadFilterTransformer();
            return Utils.stream(reads)
                    .map(preTransformer)
                    .filter(filter)
                    .map(postTransformer);
        }
        // returns an empty Stream if there are no reads
        return Stream.empty();
    }

    /**
     * @return Default size in MB of the cloud prefetch buffer. May be overridden by individual tools.
     *         The default implementation returns a value (40 MB) that is suitable for tools with a small
     *         number of large cloud inputs. Tools with large numbers of cloud inputs will likely want to
     *         override to specify a smaller size.
     *         This value is maintained in the {@link GATKConfig} file.
     */
    public int getDefaultCloudPrefetchBufferSize() {
        return ConfigFactory.getInstance().getGATKConfig().cloudPrefetchBuffer();
    }

    /**
     * @return Default size in MB of the cloud index prefetch buffer. May be overridden by individual tools.
     *         A return value of -1 means to use the same value as returned by {@link #getDefaultCloudPrefetchBufferSize()}.
     *         The default implementation returns -1.
     *         This value is maintained in the {@link GATKConfig} file.
     */
    public int getDefaultCloudIndexPrefetchBufferSize() {
        return ConfigFactory.getInstance().getGATKConfig().cloudIndexPrefetchBuffer();
    }

    /**
     * @return String label to use for records in progress meter output. Defaults to {@link ProgressMeter#DEFAULT_RECORD_LABEL},
     *         but tools may override to provide a more appropriate label (like "reads" or "regions")
     */
    public String getProgressMeterRecordLabel() { return ProgressMeter.DEFAULT_RECORD_LABEL; }

    /**
     * @return default value does no filtering. Override to change how variants are filtered against the intervals for your tools.
     */
    public IntervalFilteringVcfWriter.Mode getVariantOutputFilteringMode(){
        return IntervalFilteringVcfWriter.Mode.ANYWHERE;
    }

    protected List<SimpleInterval> transformTraversalIntervals(final List<SimpleInterval> getIntervals, final SAMSequenceDictionary sequenceDictionary) {
        return getIntervals;
    }

    /**
     * Get the GenomicsDB read settings for the current tool
     * @return By default, just return the vanilla options
     */
    protected GenomicsDBOptions getGenomicsDBOptions() {
        return new GenomicsDBOptions(referenceArguments.getReferencePath());
    }

    /**
     * Initialize our source of reference data (or set it to null if no reference argument was provided).
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of the reference data source.
     */
    void initializeReference() {
        reference = referenceArguments.getReferencePath() != null ? ReferenceDataSource.of(referenceArguments.getReferencePath()) : null;
    }

    /**
     * Initialize our source of reads data (or set it to null if no reads argument(s) were provided).
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of the reads data source.
     */
    void initializeReads() {
        if (! readArguments.getReadPathSpecifiers().isEmpty()) {
            final SamReaderFactory factory = makeSamReaderFactory();

            reads = new ReadsPathDataSource(readArguments.getReadPaths(), readArguments.getReadIndexPaths(), factory, cloudPrefetchBuffer,
                (cloudIndexPrefetchBuffer < 0 ? cloudPrefetchBuffer : cloudIndexPrefetchBuffer));
        }
        else {
            reads = null;
        }
    }

    protected final SamReaderFactory makeSamReaderFactory() {
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());
        if (hasReference()) { // pass in reference if available, because CRAM files need it
            factory = factory.referenceSequence(referenceArguments.getReferencePath());
        }
        else if (hasCramInput()) {
            throw UserException.MISSING_REFERENCE_FOR_CRAM;
        }

        if(bamIndexCachingShouldBeEnabled()) {
            factory = factory.enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES);
        }
        return factory;
    }


    protected final boolean bamIndexCachingShouldBeEnabled() {
        return intervalArgumentCollection.intervalsSpecified() && !disableBamIndexCaching;
    }

    /**
     * Helper method that simply returns a boolean regarding whether the input has CRAM files or not.
     */
    private boolean hasCramInput() {
        return readArguments.getReadPathSpecifiers().stream().anyMatch(GATKPath::isCram);
    }

    /**
     * Initialize our source of Feature data (or set it to null if no Feature argument(s) were provided).
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of Feature data sources.
     *
     * By default, this method initializes the FeatureManager to use the lookahead cache of {@link FeatureDataSource#DEFAULT_QUERY_LOOKAHEAD_BASES} bases.
     */
    void initializeFeatures() {
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer, getGenomicsDBOptions());
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Initialize our intervals for traversal.
     *
     * If intervals were specified, requires that another data source (reads or reference) be present
     * to supply a sequence dictionary. This requirement may be removed in the future.
     *
     * Must be called after other data sources have been initialized because of the sequence dictionary
     * requirement.
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of intervals.
     */
    void initializeIntervals() {
        if ( intervalArgumentCollection.intervalsSpecified() ) {
            final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
            if ( sequenceDictionary == null ) {
                throw new UserException("We require a sequence dictionary from a reference, a source of reads, or a source of variants to process intervals.  " +
                        "Since reference and reads files generally contain sequence dictionaries, this error most commonly occurs " +
                        "for VariantWalkers that do not require a reference or reads.  You can fix the problem by passing a reference file with a sequence dictionary " +
                        "via the -R argument or you can run the tool UpdateVCFSequenceDictionary on your vcf.");
            }

            userIntervals = transformTraversalIntervals(intervalArgumentCollection.getIntervals(sequenceDictionary), sequenceDictionary);
        }
    }

    /**
     * Is a source of reference data available?
     *
     * @return true if a reference is available, otherwise false
     */
    public final boolean hasReference() {
        return reference != null;
    }

    /**
     * Are sources of reads available?
     *
     * @return true if reads are available, otherwise false
     */
    public final boolean hasReads() {
        return reads != null;
    }

    /**
     * Are sources of Features available?
     *
     * @return true if sources of Features are available, otherwise false
     */
    public final boolean hasFeatures() {
        return features != null;
    }

    /**
     * Are sources of intervals available?
     *
     * @return true if user-supplied intervals are available, otherwise false
     */
    public final boolean hasUserSuppliedIntervals() {
        return userIntervals != null;
    }

    /**
     * Does this tool require reference data? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires a reference, otherwise false
     */
    public boolean requiresReference() {
        return false;
    }

    /**
     * Does this tool require features? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires reads, otherwise false
     */
    public boolean requiresFeatures() {
        return false;
    }

    /**
     * Does this tool require reads? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires reads, otherwise false
     */
    public boolean requiresReads() {
        return false;
    }

    /**
     * Does this tool require intervals? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires intervals, otherwise false
     */
    public boolean requiresIntervals() {
        return false;
    }

    /**
     * Does this tool want to disable the progress meter? If so, override here to return true
     *
     * @return true if this tools wants to disable progress meter output, otherwise false
     */
    public boolean disableProgressMeter() {
        return false;
    }

    /**
     * Get the {@link SequenceDictionaryValidationArgumentCollection} for the tool.
     * Subclasses may override this method in order to customize validation options.
     *
     * @return a SequenceDictionaryValidationArgumentCollection
     */
    protected SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection() {
        return new SequenceDictionaryValidationArgumentCollection.StandardValidationCollection();
    }

    /**
     * Load the master sequence dictionary as specified in {@code masterSequenceDictionaryFilename}.
     * Will only load the master sequence dictionary if it has not already been loaded.
     */
    private void loadMasterSequenceDictionary() {
        if ( (masterSequenceDictionary == null) && (masterSequenceDictionaryFilename != null) ) {
            masterSequenceDictionary = ReferenceUtils.loadFastaDictionary(masterSequenceDictionaryFilename);
        }
    }

    /**
     * Returns the reference sequence dictionary if there is a reference (hasReference() == true), otherwise null.
     * @return reference sequence dictionary if any, or null
     */
    public final SAMSequenceDictionary getReferenceDictionary() {
        return reference != null ? reference.getSequenceDictionary() : null;
    }

    /**
     * Returns the master sequence dictionary if it has been set, otherwise null.
     * In practice, this should only be called after engine initialization in {@link #onStartup()}
     * @return Master sequence dictionary if specified, or null.
     */
    public final SAMSequenceDictionary getMasterSequenceDictionary() {
        return masterSequenceDictionary;
    }

    /**
     * Returns the "best available" sequence dictionary or {@code null} if there is no single best dictionary.
     *
     * The algorithm for selecting the best dictionary is as follows:
     * 1) If a master sequence dictionary was specified, use that dictionary
     * 2) if there is a reference, then the best dictionary is the reference sequence dictionary
     * 3) Otherwise, if there are reads, then the best dictionary is the sequence dictionary constructed from the reads.
     * 4) Otherwise, if there are features and the feature data source has only one dictionary, then that one is the best dictionary.
     * 5) Otherwise, the result is {@code null}.
     *
     * TODO: check interval file(s) as well for a sequence dictionary
     *
     * Subclasses may override if they prefer a different algorithm.
     *
     * @return best available sequence dictionary given our inputs or {@code null} if no one dictionary is the best one.
     */
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        if (masterSequenceDictionary != null) {
            return masterSequenceDictionary;
        } else if (hasReference()){
            return reference.getSequenceDictionary();
        } else if (hasReads()){
            return reads.getSequenceDictionary();
        } else if (hasFeatures()){
            final List<SAMSequenceDictionary> dictionaries = features.getAllSequenceDictionaries();
            //If there is just one, it clearly is the best. Otherwise, none is best.
            if (dictionaries.size() == 1){
                return dictionaries.get(0);
            }
        }
        return null;
    }

    /**
     * Returns the SAM header for the reads data source. Will be a merged header if there are multiple inputs for
     * the reads. If there is only a single input, returns its header directly. Null if there are no reads.
     *
     * @return SAM header for our source of reads (null if there are no reads)
     */
    public final SAMFileHeader getHeaderForReads() {
        return hasReads() ? reads.getHeader() : null;
    }

    /**
     * Returns the header for the specified source of Features
     * @param featureDescriptor FeatureInput whose header to retrieve
     * @param <T> type of Feature in our FeatureInput
     * @return header for the provided FeatureInput (null if we have no sources of Features)
     */
    public final <T extends Feature> Object getHeaderForFeatures( final FeatureInput<T> featureDescriptor ) {
        return hasFeatures() ? features.getHeader(featureDescriptor) : null;
    }

    /**
     * Initialize our data sources, make sure that all tool requirements for input data have been satisfied
     * and start the progress meter.
     *
     * The data sources are initialized in the following order:
     *   initializeReference
     *   initializeReads (must be initialized after reference for CRAM files)
     *   initializeFeatures
     *   initializeIntervals (must be initialized after all other sources, because they require a sequence dictionary)
     *
     * Then, data sources are checked by calls to:
     *    validateSequenceDictionaries (unless disabled by disableSequenceDictionaryValidation)
     *    checkToolRequirements
     *
     * Subclasses may extend (ie, they should call super and add functionality).
     */
    @Override
    protected void onStartup() {
        super.onStartup();

        loadMasterSequenceDictionary();

        initializeReference();

        initializeReads(); // Must be initialized after reference, in case we are dealing with CRAM and a reference is required

        initializeFeatures();

        initializeIntervals(); // Must be initialized after reference, reads and features, since intervals currently require a sequence dictionary from another data source

        if (seqValidationArguments.performSequenceDictionaryValidation()) {
            validateSequenceDictionaries();
        }

        checkToolRequirements();

        if ((getVariantOutputFilteringMode() != IntervalFilteringVcfWriter.Mode.ANYWHERE ) && userIntervals == null){
            throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required if --" + StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME + " was specified or if the tool uses interval filtering.");
        }

        initializeProgressMeter(getProgressMeterRecordLabel());
    }

    /**
     * Helper method to initialize the progress meter without exposing engine level arguements.
     */
    protected final void initializeProgressMeter(final String progressMeterRecordLabel) {

        if (hasReference()) {
            progressMeter = new ProgressMeter(secondsBetweenProgressUpdates, disableProgressMeter(), reference.getSequenceDictionary());
        }
        else {
            progressMeter = new ProgressMeter(secondsBetweenProgressUpdates, disableProgressMeter());
        }
        progressMeter.setRecordLabel(progressMeterRecordLabel);
    }

    /**
     * Validates all sequence dictionaries by checking them against each other.
     *
     * Currently, the reads dictionary is checked against the reference dictionary,
     * and both are checked against any variant dictionaries. No interval dictionaries are
     * checked at this time.
     *
     * The reference dictionary is required to be a superset of the reads dictionary if there is at least
     * one cram input; otherwise, all dictionaries are required to have a common subset of equivalent contigs,
     * without respect to relative or absolute ordering (GATKTools use contig names rather than indices,
     * and so are not sensitive to contig ordering issues).
     */
    private void validateSequenceDictionaries() {

        final SAMSequenceDictionary refDict = hasReference() ? reference.getSequenceDictionary() : null;
        final SAMSequenceDictionary readDict = hasReads() ? reads.getSequenceDictionary() : null;
        final List<SAMSequenceDictionary> featureDicts = hasFeatures() ? features.getAllSequenceDictionaries() : Collections.emptyList();

        // Check the master dictionary against the reference / reads / features
        if (masterSequenceDictionary != null) {

            final boolean requireMasterDictionaryIsSuperSet = hasCramInput();

            // Check against the reads
            if ( hasReads() ) {
                SequenceDictionaryUtils.validateDictionaries("master sequence dictionary", masterSequenceDictionary,
                        "reads", readDict, requireMasterDictionaryIsSuperSet, false);
            }

            // Check against the reference
            if ( hasReference() ) {
                SequenceDictionaryUtils.validateDictionaries("master sequence dictionary", masterSequenceDictionary,
                        "reference", refDict, requireMasterDictionaryIsSuperSet, false);
            }

            // Check against the features
            for (final SAMSequenceDictionary featureDict : featureDicts) {
                SequenceDictionaryUtils.validateDictionaries("master sequence dictionary", masterSequenceDictionary, "features", featureDict);
            }
        }

        // Check the reference dictionary against the reads dictionary
        if ( hasReference() && hasReads() ) {
            if ( hasCramInput() ) {
                // Use stricter validation for CRAM vs. the reference
                SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(refDict, readDict);
            }
            else {
                // Use standard validation settings for non-CRAM reads input vs. the reference
                SequenceDictionaryUtils.validateDictionaries("reference", refDict, "reads", readDict);
            }
        }

        // Check all Feature dictionaries against the reference and/or reads dictionaries
        // TODO: pass file names associated with each sequence dictionary into validateDictionaries(); issue #660
        final SAMSequenceDictionary bestDict = getBestAvailableSequenceDictionary();
        for ( final SAMSequenceDictionary featureDict : featureDicts ) {
            if (hasReference()){
                SequenceDictionaryUtils.validateDictionaries("reference", refDict, "features", featureDict);
            }
            if (hasReads()) {
                SequenceDictionaryUtils.validateDictionaries("reads", readDict, "features", featureDict);
            }
            if (bestDict != null) { //VariantWalkers will use DrivingVariants for best dictionary, then check all other FeatureInputs
                SequenceDictionaryUtils.validateDictionaries("best available", bestDict, "features", featureDict);
            }
        }

        // we don't currently have access to seqdicts from intervals
        //if (hasUserSuppliedIntervals()) {}
    }

    /**
     * Check that all tool requirements for input data have been satisfied.
     *
     * Must be called after data source initialization.
     */
    private void checkToolRequirements() {

        if ( requiresFeatures() && ! hasFeatures() ) {
            throw new UserException("Tool " + getClass().getSimpleName() + " requires features, but none were provided");
        }

    }

    /*
     * Create a common SAMFileWriter using the reference and read header for this tool.
     *
     * @param outputPath    - if this path has a .cram extension then a reference is required. Can not be null.
     * @param preSorted     - if true then the records must already be sorted to match the header sort order
     *
     * @throws UserException if outputFile ends with ".cram" and no reference is provided
     * @return SAMFileWriter
     */
    public final SAMFileGATKReadWriter createSAMWriter(final GATKPath outputPathSpecifier, final boolean preSorted) {
        if (!hasReference() && outputPathSpecifier.isCram()) {
            throw UserException.MISSING_REFERENCE_FOR_CRAM;
        }

        return new SAMFileGATKReadWriter(
            ReadUtils.createCommonSAMWriter(
                outputPathSpecifier.toPath(),
                referenceArguments.getReferencePath(),
                getHeaderForSAMWriter(),
                preSorted,
                createOutputBamIndex,
                createOutputBamMD5
            )
        );
    }

    /**
     * Creates a VariantContextWriter whose outputFile type is determined by
     * the vcfOutput's extension, using the best available sequence dictionary for
     * this tool, and default index, leniency and md5 generation settings.
     *
     * Deprecated, use {@link #createVCFWriter(Path)} instead.
     *
     * @param outFile output File for this writer. May not be null.
     * @returns VariantContextWriter must be closed by the caller
     */
    public VariantContextWriter createVCFWriter(final File outFile) {
        return createVCFWriter(outFile == null ? null : outFile.toPath());
    }

    /**
     * Creates a VariantContextWriter whose outputFile type is determined by
     * the vcfOutput's extension, using the best available sequence dictionary for
     * this tool, and default index, leniency and md5 generation settings.
     *
     * @param outFile output GATKPath for this writer. May not be null.
     * @returns VariantContextWriter must be closed by the caller
     */
    public VariantContextWriter createVCFWriter(final GATKPath outFile) {
        return createVCFWriter(outFile == null ? null : outFile.toPath());
    }

    /**
     * Creates a VariantContextWriter whose outputFile type is determined by
     * vcfOutput's extension, using the best available sequence dictionary for
     * this tool, and default index, leniency and md5 generation settings.
     *
     * @param outPath output Path for this writer. May not be null.
     * @returns VariantContextWriter must be closed by the caller
     */
    public VariantContextWriter createVCFWriter(final Path outPath) {
        Utils.nonNull(outPath);

        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        List<Options> options = new ArrayList<>();
        if (lenientVCFProcessing) {
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        }

        if (createOutputVariantIndex) {
            if (null == sequenceDictionary) {
                logger.warn("A variant index will not be created - a sequence dictionary is required to create an output index");
                // fall through and create without an index
            } else {
                options.add(Options.INDEX_ON_THE_FLY);
            }
        }

        if (outputSitesOnlyVCFs) {
            options.add(Options.DO_NOT_WRITE_GENOTYPES);
        }
        final VariantContextWriter unfilteredWriter;
        if (maxVariantsPerShard > 0) {
            unfilteredWriter = new ShardingVCFWriter(
                    outPath,
                    maxVariantsPerShard,
                    sequenceDictionary,
                    createOutputVariantMD5,
                    options.toArray(new Options[0]));
        } else {
            unfilteredWriter = GATKVariantContextUtils.createVCFWriter(
                    outPath,
                    sequenceDictionary,
                    createOutputVariantMD5,
                    options.toArray(new Options[0]));
        }

        return getVariantOutputFilteringMode() == IntervalFilteringVcfWriter.Mode.ANYWHERE ?
                unfilteredWriter :
                new IntervalFilteringVcfWriter(unfilteredWriter,
                        intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()),
                        getVariantOutputFilteringMode());
    }

    /**
     * Returns the SAM header suitable for writing SAM/BAM/CRAM files produced by this tool.
     *
     * The default implementation calls {@link #getHeaderForReads} (and makes an empty header if that call returns null)
     * and optionally adds program tag to the header with a program version {@link #getVersion()}, program name {@link #getToolName()}
     * and command line {@link #getCommandLine()}.
     *
     * Subclasses may override.
     *
     * @return SAM header for the SAM writer with (optionally, if {@link #addOutputSAMProgramRecord} is true) program record appropriately.
     */
    protected SAMFileHeader getHeaderForSAMWriter(){
        final SAMFileHeader header = getHeaderForReads() == null ? new SAMFileHeader(): getHeaderForReads();
        if (addOutputSAMProgramRecord) {
            final SAMProgramRecord programRecord = new SAMProgramRecord(createProgramGroupID(header));
            programRecord.setProgramVersion(getVersion());
            programRecord.setCommandLine(getCommandLine());
            programRecord.setProgramName(getToolName());
            header.addProgramRecord(programRecord);
        }
        return header;
    }

    /**
     * @return If addOutputVCFCommandLine is true, a set of VCF header lines containing the tool name, version,
     * date and command line, otherwise an empty set.
     */
    protected Set<VCFHeaderLine> getDefaultToolVCFHeaderLines() {
        if (addOutputVCFCommandLine) {
            return GATKVariantContextUtils
                    .getDefaultVCFHeaderLines(getToolkitShortName(), this.getClass().getSimpleName(),
                            getVersion(), Utils.getDateTimeForDisplay((ZonedDateTime.now())), getCommandLine());
        } else {
            return new HashSet<>();
        }
    }

    /**
     * Returns the program group ID that will be used in the SAM writer.
     * Starts with {@link #getToolName} and looks for the first available ID by appending consecutive integers.
     */
    private String createProgramGroupID(final SAMFileHeader header) {
        final String toolName = getToolName();

        String pgID = toolName;
        SAMProgramRecord record = header.getProgramRecord(pgID);
        int count = 1;
        while (record != null){
            pgID = toolName + "." + String.valueOf(count++);
            record = header.getProgramRecord(pgID);
        }
        return pgID;
    }

    /**
     * A method to allow a user to inject {@link FeatureInput}s after initialization that were not
     * specified as command-line arguments.
     *
     * @param filePath path to the Feature file to register
     * @param name what to call the Feature input
     * @param featureType class of features
     * @param featureQueryLookahead look ahead this many bases during queries that produce cache misses
     * @return The {@link FeatureInput} used as the key for this data source.
     */
    public FeatureInput<? extends Feature> addFeatureInputsAfterInitialization(final String filePath,
                                                                               final String name,
                                                                               final Class<? extends Feature> featureType,
                                                                               final int featureQueryLookahead) {

        final FeatureInput<? extends Feature> featureInput = new FeatureInput<>(filePath, name);

        // Add the FeatureInput to our FeatureManager so that it will be available for FeatureContext queries
        // from the tool
        features.addToFeatureSources(
                featureQueryLookahead,
                featureInput,
                featureType,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer,
                referenceArguments.getReferencePath()
        );

        return featureInput;
    }

    /**
     * Returns the name of this tool.
     * The default implementation returns the result of calling {@link #getToolkitShortName} followed by the simple
     * name of the class. Subclasses may override.
     */
    public String getToolName() {
        return String.format("%s %s", getToolkitShortName(), getClass().getSimpleName());
    }

    /**
     * Expose a read-only version of the raw user-supplied intervals. This can be used by tools that need to explicitly
     * traverse the intervals themselves (rather than, for example, walking the reads based on the intervals)
     *
     * @return - the raw user-supplied intervals, as an unmodifiable list
     */
    public List<SimpleInterval> getUserSuppliedIntervals() {
        return Collections.unmodifiableList(userIntervals);
    }

    /**
     * Returns the list of intervals to iterate, either limited to the user-supplied intervals or the entire reference genome if none were specified.
     * If no reference was supplied, null is returned
     */
    public List<SimpleInterval> getTraversalIntervals() {
        return hasUserSuppliedIntervals() ? userIntervals : hasReference() ? IntervalUtils.getAllIntervalsForReference(getReferenceDictionary()) : null;
    }

    /**
     * Close all data sources on shutdown.
     * Subclasses may extend (ie they must call super).
     */
    @Override
    protected void onShutdown() {
        super.onShutdown();

        if ( hasReference() ) {
            reference.close();
        }

        if ( hasReads() ) {
            reads.close();
        }

        if ( hasFeatures() ) {
            features.close();
        }
    }

    /**
     * Operations performed just prior to the start of traversal. Should be overridden by tool authors
     * who need to process arguments local to their tool or perform other kinds of local initialization.
     *
     * Default implementation does nothing.
     */
    public void onTraversalStart() {}

    /**
     * A complete traversal from start to finish. Tool authors who wish to "roll their own" traversal
     * from scratch can extend this class directly and implement this method. Walker authors should
     * instead extend a Walker class and implement the Walker-appropriate apply() method, since the
     * Walker base classes implement the various kinds of traversals for you.
     */
    public abstract void traverse();

    /**
     * Operations performed immediately after a successful traversal (ie when no uncaught exceptions were thrown during the traversal).
     * Should be overridden by tool authors who need to close local resources, etc., after traversal.
     * Also allows tools to return a value representing the traversal result, which is printed by the engine.
     *
     * Default implementation does nothing and returns null.
     *
     * @return Object representing the traversal result, or null if a tool does not return a value
     */
    public Object onTraversalSuccess() { return null; }

    @Override
    @SuppressWarnings("try")
    protected final Object doWork() {
        //this makes use of try-with-resource exception suppression to ensure that onShutdown()
        //doesn't hide casual exceptions thrown during onStartup() or doWork()
        try(final AutoCloseableNoCheckedExceptions thisTool = this::closeTool){
            onTraversalStart();
            progressMeter.start();
            traverse();
            if (!progressMeter.stopped()) {
                progressMeter.stop();
            }
            return onTraversalSuccess();
        }
    }

    /**
     * This method is called by the GATK framework at the end of the {@link #doWork} template method.
     * It is called regardless of whether the {@link #traverse} has succeeded or not.
     * It is called <em>after</em> the {@link #onTraversalSuccess} has completed (successfully or not)
     * but before the {@link #doWork} method returns.
     *
     * In other words, on successful runs both {@link #onTraversalSuccess} and {@link #closeTool} will be called (in this order) while
     * on failed runs (when {@link #traverse} causes an exception), only {@link #closeTool} will be called.
     *
     * The default implementation does nothing.
     * Subclasses should override this method to close any resources that must be closed regardless of the success of traversal.
     */
    public void closeTool(){
    }
}
