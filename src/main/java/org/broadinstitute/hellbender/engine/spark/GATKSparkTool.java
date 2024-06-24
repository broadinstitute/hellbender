package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SBIIndexWriter;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.GZIIndex;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.ZonedDateTime;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Base class for GATK spark tools that accept standard kinds of inputs (reads, reference, and/or intervals).
 * Centralizes handling of tool inputs to enforce consistency, reduce duplicated boilerplate code
 * in tools, and apply standardized validation (such as sequence dictionary validation).
 *
 * Spark tools that do not fit into this pattern should extend SparkCommandLineProgram directly instead of
 * this class.
 *
 * USAGE:
 *
 * -Tools must implement {@link #runTool}.
 *
 * -Tools should override {@link #requiresReference}, {@link #requiresReads}, and/or {@link #requiresIntervals}
 *  as appropriate to indicate required inputs.
 *
 * -Tools can query whether certain inputs are present via {@link #hasReference}, {@link #hasReads}, and
 *  {@link #hasUserSuppliedIntervals}.
 *
 * -Tools can load the reads via {@link #getReads}, access the reference via {@link #getReference}, and
 *  access the intervals via {@link #getIntervals}. Any intervals specified are automatically applied
 *  to the reads. Input metadata is available via {@link #getHeaderForReads}, {@link #getReferenceSequenceDictionary},
 *  and {@link #getBestAvailableSequenceDictionary}.
 *
 * -Tools that require a custom reference window function (extra bases of reference context around each read)
 *  may override {@link #getReferenceWindowFunction} to supply one. This function will be propagated to the
 *  reference source returned by {@link #getReference}.
 */
public abstract class GATKSparkTool extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    public static final String BAM_PARTITION_SIZE_LONG_NAME = "bam-partition-size";
    public static final String NUM_REDUCERS_LONG_NAME = "num-reducers";
    public static final String SHARDED_OUTPUT_LONG_NAME = "sharded-output";
    public static final String OUTPUT_SHARD_DIR_LONG_NAME = "output-shard-tmp-dir";
    public static final String CREATE_OUTPUT_BAM_SPLITTING_INDEX_LONG_NAME = "create-output-bam-splitting-index";
    public static final String USE_NIO = "use-nio";
    public static final String SPLITTING_INDEX_GRANULARITY = "splitting-index-granularity";

    @ArgumentCollection
    public final ReferenceInputArgumentCollection referenceArguments = requiresReference() ? new RequiredReferenceInputArgumentCollection() :  new OptionalReferenceInputArgumentCollection();

    @ArgumentCollection
    public final ReadInputArgumentCollection readArguments = requiresReads() ? new RequiredReadInputArgumentCollection() : new OptionalReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = requiresIntervals() ? new RequiredIntervalArgumentCollection() : new OptionalIntervalArgumentCollection();

    @Argument(doc = "maximum number of bytes to read from a file into each partition of reads. " +
            "Setting this higher will result in fewer partitions. Note that this will not be equal to the size of the partition in memory. " +
            "Defaults to 0, which uses the default split size (determined by the Hadoop input format, typically the size of one HDFS block).",
            fullName = BAM_PARTITION_SIZE_LONG_NAME,
            optional = true)
    protected long bamPartitionSplitSize = 0;

    @Argument(doc = "Whether to use NIO or the Hadoop filesystem (default) for reading files. " +
            "(Note that the Hadoop filesystem is always used for writing files.)",
            fullName = USE_NIO,
            optional = true)
    protected boolean useNio = false;

    @ArgumentCollection
    protected SequenceDictionaryValidationArgumentCollection sequenceDictionaryValidationArguments = getSequenceDictionaryValidationArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, shortName = StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, doc = "If true, adds a command line header line to created VCF files.", optional=true, common = true)
    public boolean addOutputVCFCommandLine = true;

    @Argument(doc = "For tools that write an output, write the output in multiple pieces (shards)",
            fullName = SHARDED_OUTPUT_LONG_NAME,
            optional = true,
            mutex = {OUTPUT_SHARD_DIR_LONG_NAME})
    protected boolean shardedOutput = false;

    @Argument(doc = "when writing a bam, in single sharded mode this directory to write the temporary intermediate output shards, if not specified .parts/ will be used",
            fullName = OUTPUT_SHARD_DIR_LONG_NAME,
            optional = true,
            mutex = {SHARDED_OUTPUT_LONG_NAME})
    protected String shardedPartsDir = null;

    @Argument(doc="For tools that shuffle data or write an output, sets the number of reducers. Defaults to 0, which gives one partition per 10MB of input.",
            fullName = NUM_REDUCERS_LONG_NAME,
            optional = true)
    protected int numReducers = 0;

    @Argument(fullName = StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME,
            shortName = StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_SHORT_NAME,
            doc = "If true, create a BAM index when writing a coordinate-sorted BAM file.", optional = true, common = true)
    public boolean createOutputBamIndex = ConfigFactory.getInstance().getGATKConfig().createOutputBamIndex();

    @Argument(fullName = CREATE_OUTPUT_BAM_SPLITTING_INDEX_LONG_NAME,
            doc = "If true, create a BAM splitting index (SBI) when writing a coordinate-sorted BAM file.", optional = true, common = true)
    public boolean createOutputBamSplittingIndex = ConfigFactory.getInstance().getGATKConfig().createOutputBamIndex();

    @Argument(fullName = SPLITTING_INDEX_GRANULARITY,
             doc = "Granularity to use when writing a splitting index, one entry will be put into the index every n reads where n is this granularity value. Smaller granularity results in a larger index with more available split points.",
             optional = true, common = true,
             minValue = 1)
    public long splittingIndexGranularity = SBIIndexWriter.DEFAULT_GRANULARITY;

    @Argument(fullName = StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME,
            shortName = StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_SHORT_NAME,
            doc = "If true, create a VCF index when writing a coordinate-sorted VCF file.", optional = true, common = true)
    public boolean createOutputVariantIndex = true;

    private ReadsSparkSource readsSource;
    private SAMFileHeader readsHeader;
    private LinkedHashMap<GATKPath, SAMFileHeader> readInputs;
    private ReferenceMultiSparkSource referenceSource;
    private SAMSequenceDictionary referenceDictionary;
    private List<SimpleInterval> userIntervals;
    protected FeatureManager features;

    /**
     * Return the list of GATKCommandLinePluginDescriptor objects to be used for this CLP.
     * Use the read filter plugin.
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
     * Does this tool require reference data? Tools that do should override to return true.
     *
     * @return true if this tool requires a reference, otherwise false
     */
    public boolean requiresReference() {
        return false;
    }

    /**
     * Does this tool require reads? Tools that do should override to return true.
     *
     * @return true if this tool requires reads, otherwise false
     */
    public boolean requiresReads() {
        return false;
    }

    /**
     * Does this tool support multiple inputs? Tools that do should override this method with the desired {@link ReadInputMergingPolicy}.
     *
     * @return doNotMerge by default
     */
    public ReadInputMergingPolicy getReadInputMergingPolicy() {
        return ReadInputMergingPolicy.doNotMerge;
    }

    public static enum ReadInputMergingPolicy {
        doNotMerge,
        concatMerge
    }

    /**
     * Does this tool require intervals? Tools that do should override to return true.
     *
     * @return true if this tool requires intervals, otherwise false
     */
    public boolean requiresIntervals() {
        return false;
    }

    /**
     * Is a source of reference data available?
     *
     * @return true if a reference is available, otherwise false
     */
    public final boolean hasReference() {
        return referenceSource != null;
    }

    /**
     * Are sources of reads available?
     *
     * @return true if reads are available, otherwise false
     */
    public final boolean hasReads() {
        return readsSource != null;
    }

    /**
     * Are sources of intervals available?
     *
     * @return true if intervals are available, otherwise false
     */
    public final boolean hasUserSuppliedIntervals() {
        return userIntervals != null;
    }

    /**
     * Window function that controls how much reference context to return for each read when
     * using the reference source returned by {@link #getReference}. Tools should override
     * as appropriate. The default function is the identity function (ie., return exactly
     * the reference bases that span each read).
     *
     * @return reference window function used to initialize the reference source
     */
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return ReferenceWindowFunctions.IDENTITY_FUNCTION;
    }

    /**
     * subclasses can override this to provide different default behavior for sequence dictionary validation
     * @return a SequenceDictionaryValidationArgumentCollection
     */
    protected SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection() {
        return new SequenceDictionaryValidationArgumentCollection.StandardValidationCollection();
    }

    /**
     * Returns the "best available" sequence dictionary. This will be the reference sequence dictionary if
     * there is a reference, otherwise it will be the sequence dictionary constructed from the reads if
     * there are reads, otherwise it will be null.
     *
     * TODO: check interval file(s) as well for a sequence dictionary
     *
     * @return best available sequence dictionary given our inputs
     */
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        return hasReference() ? referenceDictionary : (hasReads() ? readsHeader.getSequenceDictionary() : null);
    }

    /**
     * @return sequence dictionary for the reference, or null if there is no reference available
     */
    public SAMSequenceDictionary getReferenceSequenceDictionary() {
        return referenceDictionary;
    }

    /**
     * @return header for the reads, or null if there are no reads available
     */
    public SAMFileHeader getHeaderForReads() {
        return readsHeader;
    }

    /**
     * Loads the reads into a {@link JavaRDD} using the intervals specified, and filters them using
     * the filter returned by {@link #makeReadFilter}.
     *
     * If no intervals were specified, returns all the reads (both mapped and unmapped).
     *
     * @return all reads from our reads input(s) as a {@link JavaRDD}, bounded by intervals if specified,
     *         and filtered using the filter from {@link #makeReadFilter}.
     */
    public JavaRDD<GATKRead> getReads() {
        final ReadFilter filter = makeReadFilter();
        return getUnfilteredReads().filter(read -> filter.test(read));
    }

    /**
     * Loads the reads into a {@link JavaRDD} using the intervals specified, and returns them
     * without applying any filtering.
     *
     * If no intervals were specified, returns all the reads (both mapped and unmapped).
     *
     * @return all reads from our reads input(s) as a {@link JavaRDD}, bounded by intervals if specified, and unfiltered.
     */
    public JavaRDD<GATKRead> getUnfilteredReads() {
        final TraversalParameters traversalParameters;
        if ( hasUserSuppliedIntervals() ) { // intervals may have been supplied by editIntervals
            final boolean traverseUnmapped;
            if (intervalArgumentCollection.intervalsSpecified()) {
                traverseUnmapped = intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()).traverseUnmappedReads();
            } else {
                traverseUnmapped = false;
            }
            traversalParameters = new TraversalParameters(getIntervals(), traverseUnmapped);
        } else {
            traversalParameters = null;
        }

        JavaRDD<GATKRead> output = null;
        ReadsSparkSource source = readsSource;
        for (final GATKPath inputPathSpecifier : readInputs.keySet()) {
            if (output == null) {
                output = getGatkReadJavaRDD(traversalParameters, source, inputPathSpecifier);
            } else {
                output = output.union(getGatkReadJavaRDD(traversalParameters, source, inputPathSpecifier));
            }
        }
        return output;
    }

    protected JavaRDD<GATKRead> getGatkReadJavaRDD(TraversalParameters traversalParameters, ReadsSparkSource source, GATKPath inputSpecifier) {
        if (hasCramInput() && !hasReference()) {
            throw UserException.MISSING_REFERENCE_FOR_CRAM;
        }

        return source.getParallelReads(inputSpecifier, referenceArguments.getReferenceSpecifier(), traversalParameters, bamPartitionSplitSize, useNio);
    }

    /**
     * Writes the reads from a {@link JavaRDD} to an output file.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam/cram.
     * @param reads reads to write.
     */
    public void writeReads(final JavaSparkContext ctx, final String outputFile, JavaRDD<GATKRead> reads) {
        writeReads(ctx, outputFile, reads, readsHeader, true);
    }

    /**
     * Writes the reads from a {@link JavaRDD} to an output file.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam/cram.
     * @param reads reads to write.
     * @param header the header to write.
     */
    public void writeReads(final JavaSparkContext ctx, final String outputFile, JavaRDD<GATKRead> reads, SAMFileHeader header, final boolean sortReadsToHeader) {
        try {
            ReadsSparkSink.writeReads(ctx, outputFile,
                    hasReference() ? referenceArguments.getReferenceSpecifier() : null,
                    reads, header, shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                    getRecommendedNumReducers(), shardedPartsDir, createOutputBamIndex, createOutputBamSplittingIndex, sortReadsToHeader, splittingIndexGranularity);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile,"writing failed", e);
        }
    }

    /**
     * Return the recommended number of reducers for a pipeline processing the reads. The number is
     * calculated by finding the total size (in bytes) of all the files in the input path, then
     * dividing by the target split size (determined by {@link #getTargetPartitionSize()}.
     * Subclasses that want to control the recommended number of reducers should typically override
     * {@link #getTargetPartitionSize()} rather than this method.
     * @return the recommended number of reducers
     */
    public int getRecommendedNumReducers() {
        if (numReducers != 0) {
            return numReducers;
        }
        long size = readInputs.keySet().stream().mapToLong(k -> BucketUtils.dirSize(k)).sum();
        final int targetPartitionSize = getTargetPartitionSize();
        return 1 + MathUtils.toIntExactOrThrow(size / targetPartitionSize,
                                               () -> new GATKException("getRecommendedNumReducers overflowed, size=" + size + " targetPartitionSize=" + targetPartitionSize));
    }

    /**
     * Returns the size of each input partition (in bytes) that is used to determine the recommended number of reducers
     * for running a processing pipeline. The larger the number of reducers used, the smaller the amount of memory
     * each one needs.
     *
     * Defaults to 10MB, but subclasses can override to change the value. Memory intensive pipelines should decrease
     * the partition size, while pipelines with lighter memory requirements may increase the partition size.
     */
    public int getTargetPartitionSize() {
        return 10 * 1024 * 1024; // 10MB
    }

    /**
     * Helper method that simply returns a boolean regarding whether the input has CRAM files or not.
     */
    private boolean hasCramInput() {
        return readArguments.getReadPathSpecifiers().stream().anyMatch(GATKPath::isCram);
    }

    /**
     * Returns a read filter (simple or composite) that can be applied to the reads returned from {@link #getReads}.
     * This implementation combines the default read filters for this tool (returned by {@link #getDefaultReadFilters}
     * along with any read filter command line directives specified by the user (such as enabling other filters or
     * disabling default filters); and returns a single composite filter resulting from the list by and'ing them together.
     *
     * NOTE: Most tools will not need to override the method, and should only do so in order to provide custom
     * behavior or processing of the final merged read filter. To change the default read filters used by the tool,
     * override {@link #getDefaultReadFilters} instead.
     *
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
    public ReadFilter makeReadFilter() {
        return makeReadFilter(getHeaderForReads());
    }

    /**
     * Like {@link #makeReadFilter()} but with the ability to pass a different SAMFileHeader.
     */
    protected ReadFilter makeReadFilter(SAMFileHeader samFileHeader) {
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return readFilterPlugin.getMergedReadFilter(samFileHeader);
    }

    /**
     * Returns the default list of ReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} filter with all default options. Subclasses
     * can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {@link #getHeaderForReads}. The actual SAMFileHeader is propagated to the read filters
     * by {@link #makeReadFilter} after the filters have been merged with command line arguments.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(new WellformedReadFilter());
    }

    /**
     * @see GATKTool#useVariantAnnotations()
     */
    public boolean useVariantAnnotations() {
        return false;
    }

    /**
     * @see GATKTool#getDefaultVariantAnnotations()
     */
    public List<Annotation> getDefaultVariantAnnotations() {
        return Collections.emptyList();
    }

    /**
     * @see GATKTool#getDefaultVariantAnnotationGroups()
     */
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Collections.emptyList();
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
     * @see GATKTool#makeVariantAnnotations()
     */
    public Collection<Annotation> makeVariantAnnotations() {
        final GATKAnnotationPluginDescriptor annotationPlugin =
                getCommandLineParser().getPluginDescriptor(GATKAnnotationPluginDescriptor.class);
        return annotationPlugin.getResolvedInstances();
    }

    /**
     * Returns the name of the source of reads data. It can be a file name or URL. Throws if the tool has more
     * than one source.
     */
    protected String getReadSourceName(){
        if (readInputs.size() > 1) {
            throw new GATKException("Multiple ReadsDataSources specified but a single source requested by the tool");
        }
        return readInputs.keySet().stream().findFirst().get().toString();
    }

    /**
     * Returns the header for a given input.
     */
    protected SAMFileHeader getHeaderForReadsInput(final GATKPath inputPathSpecifier){
        final SAMFileHeader header = readInputs.get(inputPathSpecifier);
        if (header == null) {
            throw new GATKException(String.format("Input %s not present in tool inputs", inputPathSpecifier.getRawInputString()));
        }
        return header;
    }

    /**
     * @return our reference source, or null if no reference is present
     */
    public ReferenceMultiSparkSource getReference() {
        return referenceSource;
    }

    /**
     * @return our intervals, or null if no intervals were specified
     */
    public List<SimpleInterval> getIntervals() {
        return userIntervals;
    }

    @Override
    protected void runPipeline( JavaSparkContext sparkContext ) {
        initializeToolInputs(sparkContext);
        validateSequenceDictionaries();
        runTool(sparkContext);
    }

    /**
     * Initialize standard tool inputs.
     */
    private void initializeToolInputs(final JavaSparkContext sparkContext) {
        initializeReference();
        initializeReads(sparkContext); // reference must be initialized before reads
        initializeFeatures();
        initializeIntervals();
    }

    /**
     * Initializes our reads source (but does not yet load the reads into a {@link JavaRDD}).
     * Does nothing if no reads inputs are present.
     */
    private void initializeReads(final JavaSparkContext sparkContext) {
        if ( readArguments.getReadPathSpecifiers().isEmpty() ) {
            return;
        }

        if (getReadInputMergingPolicy() == ReadInputMergingPolicy.doNotMerge && readArguments.getReadPathSpecifiers().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for for this spark tool.");
        }

        readInputs = new LinkedHashMap<>();
        readsSource = new ReadsSparkSource(sparkContext, readArguments.getReadValidationStringency());
        for (final GATKPath input : readArguments.getReadPathSpecifiers()) {
            readInputs.put(input, readsSource.getHeader(input, referenceArguments.getReferenceSpecifier()));
        }
        readsHeader = createHeaderMerger().getMergedHeader();
    }

    /**
     * Create a header merger from the individual SAM/BAM headers in our readers
     *
     * @return a header merger containing all individual headers in this data source
     */
    private SamFileHeaderMerger createHeaderMerger() {
        return new SamFileHeaderMerger(identifySortOrder(readInputs.values()), readInputs.values(), true);
    }
    // If multiple bams have had their contents merged make no assumption about the underlying sort order
    static SAMFileHeader.SortOrder identifySortOrder(final Collection<SAMFileHeader> headers){
        if (headers.size() > 1){
            return SAMFileHeader.SortOrder.unsorted;
        } else {
            return headers.iterator().next().getSortOrder();
        }
    }

    /**
     * Initializes our reference source. Does nothing if no reference was specified.
     */
    private void initializeReference() {
        final GATKPath referencePathSpecifier = referenceArguments.getReferenceSpecifier();
        if ( referencePathSpecifier != null ) {
            referenceSource = new ReferenceMultiSparkSource(referencePathSpecifier, getReferenceWindowFunction());
            referenceDictionary = referenceSource.getReferenceSequenceDictionary(readsHeader != null ? readsHeader.getSequenceDictionary() : null);
            if (referenceDictionary == null) {
                throw new UserException.MissingReferenceDictFile(referencePathSpecifier.getRawInputString());
            }
        }
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
        features = new FeatureManager(this);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Loads our intervals using the best available sequence dictionary (as returned by {@link #getBestAvailableSequenceDictionary})
     * to parse/verify them. Does nothing if no intervals were specified.
     */
    private void initializeIntervals() {
        if ( intervalArgumentCollection.intervalsSpecified() ) {
            final SAMSequenceDictionary intervalDictionary = getBestAvailableSequenceDictionary();
            if ( intervalDictionary == null ) {
                throw new UserException("We require at least one input source that has a sequence dictionary (reference or reads) when intervals are specified");
            }

            userIntervals = intervalArgumentCollection.getIntervals(intervalDictionary);
        }
        userIntervals = editIntervals(userIntervals);
    }

    /**
     * Transform the intervals during loading.
     *
     * Developers can override this method to do custom interval handling during initialization of their GATKSparkTool
     *
     * @param rawIntervals Intervals specified on command line by user (-L).  Can be {@code null}
     * @return Transformed set of intervals.  Allowed to return non-null, if null was specified in the input.
     */
    protected List<SimpleInterval> editIntervals(final List<SimpleInterval> rawIntervals) {
        return rawIntervals;
    }

    /**
     * Validates standard tool inputs against each other.
     */
    protected void validateSequenceDictionaries() {
        if ( sequenceDictionaryValidationArguments.performSequenceDictionaryValidation() ) {
            // Validate the reference sequence dictionary against the reads sequence dictionary, if both are present,
            // using standard GATK validation settings (requiring a common subset of equivalent contigs without respect
            // to ordering).
            // Check the reference dictionary against the reads dictionary
            if ( hasReference() && hasReads() ) {
                if ( hasCramInput() ) {
                    // Use stricter validation for CRAM vs. the reference
                    SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(referenceDictionary, readsHeader.getSequenceDictionary());
                }
                else {
                    // Use standard validation settings for non-CRAM reads input vs. the reference
                    SequenceDictionaryUtils.validateDictionaries("reference", referenceDictionary, "reads", readsHeader.getSequenceDictionary());
                }
            }
        }
    }

    /**
     * Register the reference file (and associated dictionary and index) to be downloaded to every node using Spark's
     * copying mechanism ({@code SparkContext#addFile()}).
     * @param ctx the Spark context
     * @param referencePath the reference file, can be a local file or a remote path
     * @return the reference file name; the absolute path of the file can be found by a Spark task using {@code SparkFiles#get()}
     */
    protected static String addReferenceFilesForSpark(JavaSparkContext ctx, Path referencePath) {
        if (referencePath == null) {
            return null;
        }
        Path indexPath = ReferenceSequenceFileFactory.getFastaIndexFileName(referencePath);
        Path dictPath = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(referencePath);
        Path gziPath = GZIIndex.resolveIndexNameForBgzipFile(referencePath);

        ctx.addFile(referencePath.toUri().toString());
        if (Files.exists(indexPath)) {
            ctx.addFile(indexPath.toUri().toString());
        }
        if (Files.exists(dictPath)) {
            ctx.addFile(dictPath.toUri().toString());
        }
        if (Files.exists(gziPath)) {
            ctx.addFile(gziPath.toUri().toString());
        }

        return referencePath.getFileName().toString();
    }

    /**
     * Register the VCF file (and associated index) to be downloaded to every node using Spark's copying mechanism
     * ({@code SparkContext#addFile()}).
     * @param ctx the Spark context
     * @param vcfFileNames the VCF files, can be local files or remote paths
     * @return the reference file name; the absolute path of the file can be found by a Spark task using {@code SparkFiles#get()}
     */
    protected static List<String> addVCFsForSpark(JavaSparkContext ctx, List<String> vcfFileNames) {
        for (String vcfFileName : vcfFileNames) {
            String vcfIndexFileName;
            if (vcfFileName.endsWith(FileExtensions.VCF)) {
                vcfIndexFileName = vcfFileName + FileExtensions.VCF_INDEX;
            } else if (vcfFileName.endsWith(FileExtensions.COMPRESSED_VCF)) {
                vcfIndexFileName = vcfFileName + FileExtensions.COMPRESSED_VCF_INDEX;
            } else {
                throw new IllegalArgumentException("Unrecognized known sites file extension. Must be .vcf or .vcf.gz");
            }
            ctx.addFile(vcfFileName);
            if (Files.exists(IOUtils.getPath(vcfIndexFileName))) {
                ctx.addFile(vcfIndexFileName);
            }
        }
        return vcfFileNames.stream().map(name -> IOUtils.getPath(name).getFileName().toString()).collect(Collectors.toList());
    }

    /**
     * Runs the tool itself after initializing and validating inputs. Must be implemented by subclasses.
     *
     * @param ctx our Spark context
     */
    protected abstract void runTool( JavaSparkContext ctx );

}
