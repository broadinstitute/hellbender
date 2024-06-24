package org.broadinstitute.hellbender.tools.genomicsdb;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.genomicsdb.Constants;
import org.genomicsdb.GenomicsDBUtils;
import org.genomicsdb.importer.GenomicsDBImporter;
import org.genomicsdb.model.BatchCompletionCallbackFunctionArgument;
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBImportConfiguration;
import org.genomicsdb.model.GenomicsDBVidMapProto;
import org.genomicsdb.model.ImportConfig;

import java.io.File;
import java.io.IOException;
import java.io.Serial;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.CompletionException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

/**
 * Import single-sample GVCFs into GenomicsDB before joint genotyping.
 *
 * <p>The GATK4 Best Practice Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from multiple samples.
 * GenomicsDBImport offers the same functionality as CombineGVCFs and initially came from the <i>Intel-Broad Center for Genomics</i>.
 * The datastore transposes sample-centric variant information across genomic loci to make data more accessible to tools.
 * </p>
 *
 * <p>To query the contents of the GenomicsDB datastore, use
 * <a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php'>SelectVariants</a>.
 * See <a href='https://software.broadinstitute.org/gatk/documentation/article?id=11813'>Tutorial#11813</a> to get started. </p>
 *
 * <p>Details on GenomicsDB are at
 * <a href='https://github.com/GenomicsDB/GenomicsDB/wiki'>https://github.com/GenomicsDB/GenomicsDB/wiki</a>.
 * In brief, GenomicsDB utilises a data storage system optimized for storing/querying sparse arrays.
 * Genomics data is typically sparse in that each sample has few variants with respect to the entire reference genome.
 * GenomicsDB contains specialized code for genomics applications, such as VCF parsing and INFO field annotation
 * calculation.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more GVCFs produced by in HaplotypeCaller with the `-ERC GVCF` or `-ERC BP_RESOLUTION` settings, containing
 * the samples to joint-genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A GenomicsDB workspace
 * </p>
 *
 *  <h3>Usage examples</h3>
 *
 *  Provide each sample GVCF separately.
 *  <pre>
 *    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
 *      -V data/gvcfs/mother.g.vcf.gz \
 *      -V data/gvcfs/father.g.vcf.gz \
 *      -V data/gvcfs/son.g.vcf.gz \
 *      --genomicsdb-workspace-path my_database \
 *      --tmp-dir=/path/to/large/tmp \
 *      -L 20
 *  </pre>
 *
 *  Provide sample GVCFs in a map file.
 *
 *  <pre>
 *    gatk --java-options "-Xmx4g -Xms4g" \
 *       GenomicsDBImport \
 *       --genomicsdb-workspace-path my_database \
 *       --batch-size 50 \
 *       -L chr1:1000-10000 \
 *       --sample-name-map cohort.sample_map \
 *       --tmp-dir /path/to/large/tmp \
 *       --reader-threads 5
 *  </pre>
 *
 *  The sample map is a tab-delimited text file with sample_name--tab--path_to_sample_vcf per line. Using a sample map
 *  saves the tool from having to download the GVCF headers in order to determine the sample names. Sample names in
 *  the sample name map file may have non-tab whitespace, but may not begin or end with whitespace.
 *
 *  <pre>
 *  sample1      sample1.vcf.gz
 *  sample2      sample2.vcf.gz
 *  sample3      sample3.vcf.gz
 *  </pre>
 *
 *  The sample name map file may optionally contain a third column with an explicit index path/URI for each VCF:
 *
 *  <pre>
 *  sample1      sample1.vcf.gz      sample1.vcf.gz.tbi
 *  sample2      sample2.vcf.gz      sample2.vcf.gz.tbi
 *  sample3      sample3.vcf.gz      sample3.vcf.gz.tbi
 *  </pre>
 *
 *  It is also possible to specify an explicit index for only a subset of the samples:
 *
 *  <pre>
 *  sample1      sample1.vcf.gz
 *  sample2      sample2.vcf.gz      sample2.vcf.gz.tbi
 *  sample3      sample3.vcf.gz
 *  </pre>
 *
 *  Add new samples to an existing genomicsdb workspace.
 *  In the incremental import case, no intervals are specified in the command because the tool will use the same
 *  intervals used in the initial import. Sample map is also supported for incremental import.
 *  <pre>
 *    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
 *      -V data/gvcfs/mother.g.vcf.gz \
 *      -V data/gvcfs/father.g.vcf.gz \
 *      -V data/gvcfs/son.g.vcf.gz \
 *      --genomicsdb-update-workspace-path my_database \
 *      --tmp-dir /path/to/large/tmp \
 *  </pre>
 *
 *  Get Picard-style interval_list from existing workspace
 *  <pre>
 *    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
 *      --genomicsdb-update-workspace-path my_database \
 *      --output-interval-list-to-file /output/path/to/file
 *  </pre>
 *  The interval_list for the specified/existing workspace will be written to /output/path/to/file.
 *  This will output a Picard-style interval_list (with a sequence dictionary header)
 * <h3>Caveats</h3>
 * <ul>
 *     <li>IMPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!</li>
 *     <li>At least one interval must be provided, unless incrementally importing new samples in which case specified intervals are ignored in favor of intervals specified in the existing workspace</li>
 *     <li>Input GVCFs cannot contain multiple entries for a single genomic position</li>
 *     <li>The --genomicsdb-workspace-path must point to a non-existent or empty directory</li>
 *     <li>The --genomicsdb-update-workspace-path must point to a existing genomicsdb workspace</li>
 *     <li>Either --genomicsdb-workspace-path or genomicsdb-update-workspace-path must be specified</li>
 *     <li>GenomicsDBImport uses temporary disk storage during import. The amount of temporary disk storage required can exceed the space available, especially when specifying a large number of intervals. The command line argument `--tmp-dir` can be used to specify an alternate temporary storage location with sufficient space..</li>
 *     <li>It is recommended that users backup existing genomicsdb workspaces before adding new samples using --genomicsdb-update-workspace-path. If the tool fails during incremental import for any reason, the workspace may be in an inconsistent/corrupted state</li>
 * </ul>
 *
 * <h3>Developer Note</h3>
 * To read data from GenomicsDB, use the query interface {@link org.genomicsdb.reader.GenomicsDBFeatureReader}
 */
@DocumentedFeature
@CommandLineProgramProperties(
    summary = "Import VCFs to GenomicsDB",
    oneLineSummary = "Import VCFs to GenomicsDB",
    programGroup = ShortVariantDiscoveryProgramGroup.class
)
public final class GenomicsDBImport extends GATKTool {

    private static final long DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE = 16*1024L;
    private static final long DEFAULT_SEGMENT_SIZE = 1048576L;
    private static final int DEFAULT_ZERO_BATCH_SIZE = 0;

    public static final String WORKSPACE_ARG_LONG_NAME = "genomicsdb-workspace-path";
    public static final String INCREMENTAL_WORKSPACE_ARG_LONG_NAME = "genomicsdb-update-workspace-path";
    public static final String SEGMENT_SIZE_ARG_LONG_NAME = "genomicsdb-segment-size";
    public static final String OVERWRITE_WORKSPACE_LONG_NAME = "overwrite-existing-genomicsdb-workspace";

    public static final String INTERVAL_LIST_LONG_NAME = "output-interval-list-to-file";

    public static final String VCF_BUFFER_SIZE_ARG_NAME = "genomicsdb-vcf-buffer-size";

    public static final String BATCHSIZE_ARG_LONG_NAME = "batch-size";
    public static final String CONSOLIDATE_ARG_NAME = "consolidate";
    public static final String SAMPLE_NAME_MAP_LONG_NAME = "sample-name-map";
    public static final String VALIDATE_SAMPLE_MAP_LONG_NAME = "validate-sample-name-map";
    public static final String MERGE_INPUT_INTERVALS_LONG_NAME = "merge-input-intervals";
    public static final String VCF_INITIALIZER_THREADS_LONG_NAME = "reader-threads";
    public static final String MAX_NUM_INTERVALS_TO_IMPORT_IN_PARALLEL = "max-num-intervals-to-import-in-parallel";
    public static final String MERGE_CONTIGS_INTO_NUM_PARTITIONS = "merge-contigs-into-num-partitions";
    public static final String BYPASS_FEATURE_READER = "bypass-feature-reader";
    public static final String VCF_HEADER_OVERRIDE = "header";
    public static final int INTERVAL_LIST_SIZE_WARNING_THRESHOLD = 100;
    public static final int ARRAY_COLUMN_BOUNDS_START = 0;
    public static final int ARRAY_COLUMN_BOUNDS_END = 1;

    public static final String SHARED_POSIXFS_OPTIMIZATIONS = GenomicsDBArgumentCollection.SHARED_POSIXFS_OPTIMIZATIONS;
    public static final String USE_GCS_HDFS_CONNECTOR = GenomicsDBArgumentCollection.USE_GCS_HDFS_CONNECTOR;
    public static final String AVOID_NIO = "avoid-nio";

    @Argument(fullName = WORKSPACE_ARG_LONG_NAME,
              doc = "Workspace for GenomicsDB. Can be a POSIX file system absolute or relative path or a HDFS/GCS URL. " +
                    "Use this argument when creating a new GenomicsDB workspace. " +
                    "Either this or "+INCREMENTAL_WORKSPACE_ARG_LONG_NAME+" must be specified." +
                    " Must be an empty or non-existent directory.",
              mutex = {INCREMENTAL_WORKSPACE_ARG_LONG_NAME,INTERVAL_LIST_LONG_NAME})
    private String workspace;

    @Argument(fullName = INCREMENTAL_WORKSPACE_ARG_LONG_NAME,
              doc = "Workspace when updating GenomicsDB. Can be a POSIX file system absolute or relative path or a HDFS/GCS URL. " +
                    "Use this argument when adding new samples to an existing GenomicsDB workspace or "+
                    "when using the "+INTERVAL_LIST_LONG_NAME+" option. " +
                    "Either this or "+WORKSPACE_ARG_LONG_NAME+" must be specified. " +
                    "Must point to an existing workspace.",
              mutex = {WORKSPACE_ARG_LONG_NAME, VCF_HEADER_OVERRIDE})
    private String incrementalImportWorkspace;

    @Argument(fullName = SEGMENT_SIZE_ARG_LONG_NAME,
              doc = "Buffer size in bytes allocated for GenomicsDB attributes during " +
                    "import. Should be large enough to hold data from one site. ",
              optional = true)
    private long segmentSize = DEFAULT_SEGMENT_SIZE;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
              shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
              doc = "GVCF files to be imported to GenomicsDB. Each file must contain" +
                    " data for only a single sample. Either this or " + SAMPLE_NAME_MAP_LONG_NAME +
                    " must be specified.",
              optional = true,
              mutex = {SAMPLE_NAME_MAP_LONG_NAME, AVOID_NIO})
    private List<String> variantPaths;

    @Argument(fullName = VCF_BUFFER_SIZE_ARG_NAME,
              doc = "Buffer size in bytes to store variant contexts." +
                    " Larger values are better as smaller values cause frequent disk writes." +
                    " Defaults to " + DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE + " which was empirically determined to work" +
                    " well for many inputs.",
              optional = true,
              minValue = 1024L,
              minRecommendedValue = 10 * 1024)
    private long vcfBufferSizePerSample = DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE;

    @Argument(fullName = OVERWRITE_WORKSPACE_LONG_NAME,
              doc = "Will overwrite given workspace if it exists. " +
                    "Otherwise a new workspace is created. " +
                    "Cannot be set to true if "+INCREMENTAL_WORKSPACE_ARG_LONG_NAME+" is also set. " +
                    "Defaults to false",
              optional = true)
    private Boolean overwriteExistingWorkspace = false;

    @Argument(fullName = BATCHSIZE_ARG_LONG_NAME,
              doc = "Batch size controls the number of samples for which readers are open at once " +
                    "and therefore provides a way to minimize memory consumption. However, it can take longer to complete. " +
                    "Use the consolidate flag if more than a hundred batches were used. This will improve feature read time. " +
                    "batchSize=0 means no batching (i.e. readers for all samples will be opened at once) " +
                    "Defaults to " + DEFAULT_ZERO_BATCH_SIZE,
              optional = true)
    private int batchSize = DEFAULT_ZERO_BATCH_SIZE;

    @Argument(fullName = CONSOLIDATE_ARG_NAME,
              doc = "Boolean flag to enable consolidation. If importing data in batches, a new fragment is created for " +
                    "each batch. In case thousands of fragments are created, GenomicsDB feature readers will try " +
                    "to open ~20x as many files. Also, internally GenomicsDB would consume more memory to maintain " +
                    "bookkeeping data from all fragments. Use this flag to merge all fragments into one. " +
                    "Merging can potentially improve read performance, however overall benefit might not be noticeable " +
                    "as the top Java layers have significantly higher overheads. This flag has no effect if only one " +
                    "batch is used. Defaults to false",
              optional = true)
    private Boolean doConsolidation = false;

    @Advanced
    @Argument(fullName = SAMPLE_NAME_MAP_LONG_NAME,
            doc = "Path to file containing a mapping of sample name to file uri in tab delimited format.  If this is " +
                    "specified then the header from the first sample will be treated as the merged header rather than " +
                    "merging the headers, and the sample names will be taken from this file.  This may be used to rename " +
                    "input samples. This is a performance optimization that relaxes the normal checks for consistent " +
                    "headers.  Using vcfs with incompatible headers may result in silent data corruption.",
            optional = true,
            mutex = {StandardArgumentDefinitions.VARIANT_LONG_NAME})
    private String sampleNameMapFile;

    @Argument(fullName = VALIDATE_SAMPLE_MAP_LONG_NAME,
            doc = "Boolean flag to enable checks on the sampleNameMap file. If true, tool checks whether" +
                "feature readers are valid and shows a warning if sample names do not match with the headers. " +
                "Defaults to false",
            optional = true)
    private Boolean validateSampleToReaderMap = false;

    @Argument(fullName = MERGE_INPUT_INTERVALS_LONG_NAME,
            doc = "Boolean flag to import all data in between intervals.  Improves performance using large lists of " +
                "intervals, as in exome sequencing, especially if GVCF data only exists for specified intervals.")
    private boolean mergeInputIntervals = false;

    @Argument(fullName = INTERVAL_LIST_LONG_NAME,
            doc = "Path to output file where intervals from existing workspace should be written." +
                  "If this option is specified, the tools outputs the interval_list of the workspace pointed to " +
                  "by "+INCREMENTAL_WORKSPACE_ARG_LONG_NAME+" at the path specified here in a Picard-style interval_list " +
                  "with a sequence dictionary header",
            optional= true,
            mutex = {WORKSPACE_ARG_LONG_NAME})
    private String intervalListOutputPathString;

    @Advanced
    @Argument(fullName = VCF_INITIALIZER_THREADS_LONG_NAME,
            doc = "How many simultaneous threads to use when opening VCFs in batches; higher values may improve performance " +
                    "when network latency is an issue. Multiple reader threads are not supported when running with multiple intervals.",
            optional = true,
            minValue = 1)
    private int vcfInitializerThreads = 1;

    @Advanced
    @Argument(fullName = MAX_NUM_INTERVALS_TO_IMPORT_IN_PARALLEL,
            doc = "Max number of intervals to import in parallel; higher values may improve performance, but require more" +
                  " memory and a higher number of file descriptors open at the same time",
            optional = true,
            minValue = 1)
    private int maxNumIntervalsToImportInParallel = 1;

    @Advanced
    @Argument(fullName = MERGE_CONTIGS_INTO_NUM_PARTITIONS,
            shortName = MERGE_CONTIGS_INTO_NUM_PARTITIONS,
            doc = "Number of GenomicsDB arrays to merge input intervals into. Defaults to 0, which disables this merging. " +
                  "This option can only be used if entire contigs are specified as intervals. The tool will not split up " +
                  "a contig into multiple arrays, which means the actual number of partitions may be less than what is " +
                  "specified for this argument. This can improve performance in the " +
                  "case where the user is trying to import a very large number of contigs - larger than " +
                  INTERVAL_LIST_SIZE_WARNING_THRESHOLD,
            optional = true,
            minValue = 0)
    private int mergeContigsIntoNumPartitions = 0;

    @Argument(fullName = SHARED_POSIXFS_OPTIMIZATIONS,
            doc = "Allow for optimizations to improve the usability and performance for shared Posix Filesystems(e.g. NFS, Lustre). " +
                  "If set, file level locking is disabled and file system writes are minimized by keeping a higher number of " +
                  "file descriptors open for longer periods of time. Use with " + BATCHSIZE_ARG_LONG_NAME + " option if keeping a " +
                  "large number of file descriptors open is an issue",
            optional = true)
    private boolean sharedPosixFSOptimizations = false;

    @Argument(fullName = VCF_HEADER_OVERRIDE,
        doc = "Specify a vcf file to use instead of reading and combining headers from the input vcfs",
        optional = true,
        mutex ={INCREMENTAL_WORKSPACE_ARG_LONG_NAME}
    )
    private FeatureInput<VariantContext> headerOverride = null;

    @Argument(fullName = BYPASS_FEATURE_READER,
            doc = "Use htslib to read input VCFs instead of GATK's FeatureReader. This will reduce memory usage and potentially speed up " +
                  "the import. Lower memory requirements may also enable parallelism through " + MAX_NUM_INTERVALS_TO_IMPORT_IN_PARALLEL +
                  ". To enable this option, VCFs must be normalized, block-compressed and indexed.",
            optional = true)
    private boolean bypassFeatureReader = false;

    @Argument(fullName = AVOID_NIO,
                doc = "Do not attempt to open the input vcf file paths in java.  This can only be used with " + BYPASS_FEATURE_READER
            + ".  It allows operating on file systems which GenomicsDB understands how to open but GATK does not.  This will disable "
            + "many of the sanity checks.",
            mutex = {StandardArgumentDefinitions.VARIANT_LONG_NAME}
    )
    @Advanced
    private boolean avoidNio = false;

    @Argument(fullName = USE_GCS_HDFS_CONNECTOR,
            doc = "Use the GCS HDFS Connector instead of the native GCS SDK client with gs:// URLs.",
            optional = true)
    public boolean useGcsHdfsConnector = false;

    //executor service used when vcfInitializerThreads > 1
    private ExecutorService inputPreloadExecutorService;

    /**
     * Get the largest interval per contig that contains the intervals specified on the command line.
     * @param getIntervals intervals to be transformed
     * @param sequenceDictionary used to validate intervals
     * @return a list of one interval per contig spanning the input intervals after processing and validation
     */
    @Override
    protected List<SimpleInterval> transformTraversalIntervals(final List<SimpleInterval> getIntervals, final SAMSequenceDictionary sequenceDictionary) {
        if (mergeInputIntervals) {
            return IntervalUtils.getSpanningIntervals(getIntervals, sequenceDictionary);
        } else {
            return getIntervals;
        }
    }

    @Override
    public boolean disableProgressMeter() {
        return true;
    }

    @Override
    public int getDefaultCloudPrefetchBufferSize() {
        // Empirical testing has shown that this tool performs best at scale with cloud buffering
        // disabled. With cloud buffering on and thousands of concurrent GenomicsDBImport tasks,
        // we do too many simultaneous GCS accesses (since the prefetcher spawns a new thread for each
        // reader upon a query) and start seeing intermittent failures, even with aggressive retries.
        return 0;
    }

    @Override
    public int getDefaultCloudIndexPrefetchBufferSize() {
        // Empirical testing has shown that this tool performs best at scale with cloud buffering
        // disabled. With cloud buffering on and thousands of concurrent GenomicsDBImport tasks,
        // we do too many simultaneous GCS accesses (since the prefetcher spawns a new thread for each
        // reader upon a query) and start seeing intermittent failures, even with aggressive retries.
        return 0;
    }

    @Override
    public String getProgressMeterRecordLabel() { return "batches"; }

    // Intervals from command line (merged if specified)
    private List<SimpleInterval> intervals;

    private SampleNameMap sampleNameMap;

    // Needed as smartMergeHeaders() returns a set of VCF header lines
    private Set<VCFHeaderLine> mergedHeaderLines = null;

    // sequence dictionary created from the merged header
    private SAMSequenceDictionary mergedHeaderSequenceDictionary;

    // Path to vidmap file to be written by GenomicsDBImporter
    private String vidMapJSONFile;

    // Path to callsetmap file to be written by GenomicsDBImporter
    private String callsetMapJSONFile;

    // Path to combined VCF header file to be written by GenomicsDBImporter
    private String vcfHeaderFile;

    //in-progress batchCount
    private int batchCount = 1;

    // true if genomicsdb-update-workspace-path is specified
    private Boolean doIncrementalImport = false;

    // true if --output-interval-list-file is specified
    private Boolean getIntervalsFromExistingWorkspace = false;

    /**
     * Before traversal starts, create the feature readers
     * for all the input GVCFs, create the merged header and
     * initialize the interval
     */
    @Override
    public void onStartup() {
        initializeWorkspaceAndToolMode();
        assertVariantPathsOrSampleNameFileWasSpecified();
        assertOverwriteWorkspaceAndIncrementalImportMutuallyExclusive();
        assertAvoidNioConditionsAreValid();
        initializeHeaderAndSampleMappings();
        initializeIntervals();
        super.onStartup();
    }



    private void initializeWorkspaceAndToolMode() {
        if (incrementalImportWorkspace != null && !incrementalImportWorkspace.isEmpty()) {
            doIncrementalImport = true;
            workspace = incrementalImportWorkspace;
            if (intervalListOutputPathString != null && !intervalListOutputPathString.isEmpty()) {
                getIntervalsFromExistingWorkspace = true;
            }
        }
        if (useGcsHdfsConnector) {
            GenomicsDBUtils.useGcsHdfsConnector(true);
        }
    }

    private void assertOverwriteWorkspaceAndIncrementalImportMutuallyExclusive() {
        if (overwriteExistingWorkspace && doIncrementalImport) {
            throw new CommandLineException(OVERWRITE_WORKSPACE_LONG_NAME + " cannot be set to true when " +
                    INCREMENTAL_WORKSPACE_ARG_LONG_NAME + " is set");
        }
    }

    private void assertVariantPathsOrSampleNameFileWasSpecified(){
        if ( (variantPaths == null || variantPaths.isEmpty()) && sampleNameMapFile == null && !getIntervalsFromExistingWorkspace) {
            throw new CommandLineException.MissingArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME,
                                                       "One of --" + StandardArgumentDefinitions.VARIANT_LONG_NAME + " or --" + SAMPLE_NAME_MAP_LONG_NAME + " must be specified" );
        }
    }

    private void assertAvoidNioConditionsAreValid() {
        if (avoidNio && (!bypassFeatureReader || headerOverride == null) ){
            final List<String> missing = new ArrayList<>();
            if(!bypassFeatureReader){
                missing.add(BYPASS_FEATURE_READER);
            }
            if(headerOverride == null){
                missing.add(VCF_HEADER_OVERRIDE);
            }
            final String missingArgs = String.join(" and ", missing);

            // this potentially produces and exception with bad grammar but that's probably ok
            throw new CommandLineException.MissingArgument(missingArgs, "If --" +AVOID_NIO + " is set then --" + BYPASS_FEATURE_READER
                    + " and --" + VCF_HEADER_OVERRIDE + " must also be specified.");

        }
    }

    private static void assertIntervalsCoverEntireContigs(GenomicsDBImporter importer,
                                                          List<SimpleInterval> intervals) {
        GenomicsDBVidMapProto.VidMappingPB vidMapPB = importer.getProtobufVidMapping();
        if (vidMapPB == null) {
            throw new UserException("Could not get protobuf vid mappping object from GenomicsDBImporter");
        }
        Map<String,GenomicsDBVidMapProto.Chromosome> vidContigs =
                vidMapPB.getContigsList().stream().collect(Collectors.toMap(item->item.getName(), item->item));
        for (SimpleInterval interval: intervals) {
            GenomicsDBVidMapProto.Chromosome vidContig = vidContigs.get(interval.getContig());
            long contigLength = vidContig.getLength();
            if (interval.getStart() != 1 || interval.getEnd() < contigLength) {
                String inputInterval = String.format("Contig:%s, Start:%d, End:%d",
                        interval.getContig(), interval.getStart(), interval.getEnd());
                String vidInterval = String.format("Contig:%s, Start:%d, End:%d",
                        vidContig.getName(), 1, vidContig.getLength());
                throw new UserException("--" + MERGE_CONTIGS_INTO_NUM_PARTITIONS + " requires that entire contigs " +
                        "be specified for input intervals. Input interval contained: " + inputInterval +
                        " while reference contig was: "+ vidInterval);
            }
        }
    }

    /**
     * sets the values of mergedHeaderLines, mergedHeaderSequenceDictionary, and sampleNameToVcfPath
     */
    private void initializeHeaderAndSampleMappings() {
        // Only one of -V and --sampleNameMapFile may be specified
        if (variantPaths != null && !variantPaths.isEmpty()) {
            // -V was specified
            final List<VCFHeader> headers = new ArrayList<>(variantPaths.size());
            sampleNameMap = new SampleNameMap();

            if(headerOverride == null) {
                for (final String variantPathString : variantPaths) {
                    final Path variantPath = IOUtils.getPath(variantPathString);
                    if (bypassFeatureReader) { // avoid-nio can't be set here because it requires headerOverride
                        GATKGenomicsDBUtils.assertVariantFileIsCompressedAndIndexed(variantPath);
                    }
                    final VCFHeader header = getHeaderFromPath(variantPath);
                    Utils.validate(header != null, "Null header was found in " + variantPath + ".");
                    assertGVCFHasOnlyOneSample(variantPathString, header);
                    headers.add(header);

                    final String sampleName = header.getGenotypeSamples().get(0);
                    try {
                        sampleNameMap.addSample(sampleName, new URI(variantPathString));
                    } catch (final URISyntaxException e) {
                        throw new UserException("Malformed URI " + e.getMessage(), e);
                    }
                }
                mergedHeaderLines = VCFUtils.smartMergeHeaders(headers, true);
                mergedHeaderSequenceDictionary = new VCFHeader(mergedHeaderLines).getSequenceDictionary();
            } else {
                final VCFHeader header = getHeaderFromPath(headerOverride.toPath());
                mergedHeaderLines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
                mergedHeaderSequenceDictionary = header.getSequenceDictionary();
            }
            mergedHeaderLines.addAll(getDefaultToolVCFHeaderLines());
        } else if (sampleNameMapFile != null) {
            // --sampleNameMap was specified

            //it's VERY IMPORTANT that this map is Sorted according to String's natural ordering, if it is not
            //the resulting database will have incorrect sample names
            //see https://github.com/broadinstitute/gatk/issues/3682 for more information
            // The SampleNameMap class guarantees that the samples will be sorted correctly.
            sampleNameMap = new SampleNameMap(IOUtils.getPath(sampleNameMapFile),
                    bypassFeatureReader && !avoidNio);

            final String firstSample = sampleNameMap.getSampleNameToVcfPath().entrySet().iterator().next().getKey();

            final VCFHeader header;
            if(headerOverride == null){
                final Path firstVCFPath = sampleNameMap.getVCFForSampleAsPath(firstSample);
                header = getHeaderFromPath(firstVCFPath);
            } else {
                header = getHeaderFromPath(headerOverride.toPath());
            }
            //getMetaDataInInputOrder() returns an ImmutableSet - LinkedHashSet is mutable and preserves ordering
            mergedHeaderLines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
            mergedHeaderSequenceDictionary = header.getSequenceDictionary();
            mergedHeaderLines.addAll(getDefaultToolVCFHeaderLines());
        } else if (getIntervalsFromExistingWorkspace){
            final String vcfHeader = IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME);
            IOUtils.assertPathsAreReadable(vcfHeader);
            final String header = GenomicsDBUtils.readEntireFile(vcfHeader);
            try {
                File tempHeader = IOUtils.createTempFile("tempheader", ".vcf");
                Files.writeString(tempHeader.toPath(), header);
                mergedHeaderSequenceDictionary = VCFFileReader.getSequenceDictionary(tempHeader);
            } catch (final IOException e) {
                throw new UserException("Unable to create temporary header file to get sequence dictionary", e);
            }
        } else {
            throw new UserException(StandardArgumentDefinitions.VARIANT_LONG_NAME+" or "+
                    SAMPLE_NAME_MAP_LONG_NAME+" must be specified unless "+
                    INTERVAL_LIST_LONG_NAME+" is specified");
        }


        if ( mergedHeaderSequenceDictionary == null) {
            throw new UserException("The merged vcf header has no sequence dictionary. Please provide a header that contains a sequence dictionary.");
        }

        // If any indices were specified in the sample name map file, make sure
        // that --bypass-feature-reader wasn't also specified:
        if ( sampleNameMap != null && sampleNameMap.indicesSpecified() && bypassFeatureReader ) {
            throw new UserException("Indices were specified for some VCFs in the sample name map file, but --" + BYPASS_FEATURE_READER +
                    " was also specified. Specifying explicit indices is not supported when running with --" + BYPASS_FEATURE_READER);
        }
    }

    private VCFHeader getHeaderFromPath(final Path variantPath) {
        //TODO make this mangling unecessary
        final String variantURI = variantPath.toAbsolutePath().toUri().toString();
        try(final FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(), false,
                BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer),
                BucketUtils.getPrefetchingWrapper(cloudIndexPrefetchBuffer))) {
            return (VCFHeader) reader.getHeader();
        } catch (final IOException e) {
            throw new UserException("Error while reading vcf header from " + variantPath.toUri(), e);
        }
    }

    private static void assertGVCFHasOnlyOneSample(final String variantPath, final VCFHeader header) {
        // A GVCF file must contain only one sample, throw an exception otherwise
        final int numberOfSamples = header.getNGenotypeSamples();
        if (numberOfSamples != 1) {
            throw new UserException("Input GVCF: " + variantPath + " was expected to contain a single sample but actually contained " + numberOfSamples + " samples.");
        }
    }



    /**
     * write out interval list to file
     */
    private void writeIntervalListToFile() {
        final IntervalList outputList = new IntervalList(getBestAvailableSequenceDictionary());
        intervals.forEach(i -> outputList.add(new Interval(i)));
        outputList.write(IOUtils.getPath(intervalListOutputPathString));
    }

    /**
     * Before traversal, fix configuration parameters and initialize
     * GenomicsDB. Hard-coded to handle only VCF files and headers
     */
    @Override
    public void onTraversalStart() {
        String workspaceDir = overwriteCreateOrCheckWorkspace();
        vidMapJSONFile = GATKGenomicsDBUtils.genomicsDBApppendPaths(workspaceDir, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME);
        callsetMapJSONFile = GATKGenomicsDBUtils.genomicsDBApppendPaths(workspaceDir, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME);
        vcfHeaderFile = GATKGenomicsDBUtils.genomicsDBApppendPaths(workspaceDir, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME);
        if (getIntervalsFromExistingWorkspace) {
            // intervals may be null if merge-contigs-into-num-partitions was used to create the workspace
            // if so, we need to wait for vid to be generated before writing out the interval list
            logger.info("Interval file list will be written out to " + intervalListOutputPathString);
        }
        else if (doIncrementalImport) {
            logger.info("Callset Map JSON file will be re-written to " + callsetMapJSONFile);
            logger.info("Incrementally importing to workspace - " + workspaceDir);
        } else {
            logger.info("Vid Map JSON file will be written to " + vidMapJSONFile);
            logger.info("Callset Map JSON file will be written to " + callsetMapJSONFile);
            logger.info("Complete VCF Header will be written to " + vcfHeaderFile);
            logger.info("Importing to workspace - " + workspaceDir);
        }
        initializeInputPreloadExecutorService();
    }

    private void initializeInputPreloadExecutorService() {
        if( vcfInitializerThreads > 1) {
            if( intervals != null && intervals.size() == 1) {
                final ThreadFactory threadFactory = new ThreadFactoryBuilder()
                    .setNameFormat("readerInitializer-thread-%d")
                    .setDaemon(true)
                    .build();
                this.inputPreloadExecutorService = Executors.newFixedThreadPool(vcfInitializerThreads, threadFactory);
            }
            else {
                logger.warn("GenomicsDBImport cannot use multiple VCF reader threads for initialization when the "
                    + "number of intervals is greater than 1. Falling back to serial VCF reader initialization.");
                inputPreloadExecutorService = null;
            }
        } else {
            inputPreloadExecutorService = null;
        }
    }

    private Map<String, FeatureReader<VariantContext>> createSampleToReaderMap(
            final Map<String, URI> sampleNameToVcfPath, final int batchSize, final int index) {
        // TODO: fix casting since it's really ugly
        return inputPreloadExecutorService != null ?
                getFeatureReadersInParallel((SortedMap<String, URI>) sampleNameToVcfPath, batchSize, index)
                : getFeatureReadersSerially(sampleNameToVcfPath, batchSize, index);
    }

    private Void logMessageOnBatchCompletion(final BatchCompletionCallbackFunctionArgument arg) {
        logger.info("Done importing batch " + arg.batchCount + "/" + arg.totalBatchCount);
        logger.debug("List of samples imported in batch " + arg.batchCount + ":");
        int index = 0;
        final int sampleCount = sampleNameMap.getNumSamples();
        final int updatedBatchSize = (batchSize == DEFAULT_ZERO_BATCH_SIZE) ? sampleCount : batchSize;
        final int startBatch = (arg.batchCount - 1) * updatedBatchSize;
        final int stopBatch = arg.batchCount * updatedBatchSize;
        for(String key : sampleNameMap.getSampleNamesInSortedOrder()) {
            index++;
            if (index <= startBatch || index > stopBatch) {
                continue;
            }
            logger.debug("\t"+key);
        }
        this.batchCount = arg.batchCount + 1;
        return null;
    }

    private GenomicsDBImportConfiguration.Partition createPartitionWithBeginAndEnd(
            Coordinates.GenomicsDBColumn begin, Coordinates.GenomicsDBColumn end) {
        GenomicsDBImportConfiguration.Partition.Builder partitionBuilder = GenomicsDBImportConfiguration.Partition.newBuilder();
        partitionBuilder.setBegin(begin);
        partitionBuilder.setEnd(end);
        partitionBuilder.setWorkspace(workspace);
        partitionBuilder.setGenerateArrayNameFromPartitionBounds(true);
        return partitionBuilder.build();
    }

    private List<GenomicsDBImportConfiguration.Partition> generatePartitionListFromWorkspace() {
        String[] partitions = GenomicsDBUtils.listGenomicsDBArrays(workspace);
        // sort here to keep partition order
        Arrays.sort(partitions);
        List<GenomicsDBImportConfiguration.Partition> configPartitions = new ArrayList<>();
        for (String partition : partitions) {
            long[] bounds = GenomicsDBUtils.getArrayColumnBounds(workspace, partition);
            Coordinates.GenomicsDBColumn.Builder beginBuilder = Coordinates.GenomicsDBColumn.newBuilder();
            Coordinates.GenomicsDBColumn.Builder endBuilder = Coordinates.GenomicsDBColumn.newBuilder();
            beginBuilder.setTiledbColumn(bounds[ARRAY_COLUMN_BOUNDS_START]);
            endBuilder.setTiledbColumn(bounds[ARRAY_COLUMN_BOUNDS_END]);
            configPartitions.add(createPartitionWithBeginAndEnd(beginBuilder.build(), endBuilder.build()));
        }
        return configPartitions;
    }

    private List<GenomicsDBImportConfiguration.Partition> generatePartitionListFromIntervals() {
        return intervals.stream().map(interval -> {
            Coordinates.ContigPosition.Builder contigPositionBuilder = Coordinates.ContigPosition.newBuilder();
            contigPositionBuilder.setContig(interval.getContig());
            Coordinates.GenomicsDBColumn.Builder beginBuilder = Coordinates.GenomicsDBColumn.newBuilder();
            Coordinates.GenomicsDBColumn.Builder endBuilder = Coordinates.GenomicsDBColumn.newBuilder();
            //begin
            contigPositionBuilder.setPosition(interval.getStart());
            beginBuilder.setContigPosition(contigPositionBuilder.build());
            //end
            contigPositionBuilder.setPosition(interval.getEnd());
            endBuilder.setContigPosition(contigPositionBuilder.build());
            return createPartitionWithBeginAndEnd(beginBuilder.build(), endBuilder.build());
        }).collect(Collectors.toList());
    }

    private List<SimpleInterval> generateIntervalListFromVidMap() {
        try {
            GenomicsDBVidMapProto.VidMappingPB vidMapPB = 
                GATKGenomicsDBUtils.getProtobufVidMappingFromJsonFile(vidMapJSONFile);
    
            List<String> partitions = Arrays.asList(GenomicsDBUtils.listGenomicsDBArrays(workspace));
            return partitions.stream().flatMap(partition -> {
                long[] bounds = GenomicsDBUtils.getArrayColumnBounds(workspace, partition);
                // merge-contigs-into-num-partitions ensures entire contigs are within a given partition
                // so we just check here that contig starts within the given bounds
                return vidMapPB.getContigsList().stream()
                        .filter(x -> x.getTiledbColumnOffset() >= bounds[ARRAY_COLUMN_BOUNDS_START] &&  
                        x.getTiledbColumnOffset() <= bounds[ARRAY_COLUMN_BOUNDS_END])
                        .map(x -> new SimpleInterval(x.getName(), 1, Math.toIntExact(x.getLength())));
            }).collect(Collectors.toList());
        } catch (final IOException e) {
            throw new UserException("Could not get vid map protobuf from file:" + vidMapJSONFile + 
                    ". Is the workspace corrupted?", e);
        }
    }

    private List<SimpleInterval> generateIntervalListFromWorkspace() {
        List<String> partitions = Arrays.asList(GenomicsDBUtils.listGenomicsDBArrays(workspace));
        return partitions.stream().map(partition -> {
            String[] partitionInfo = partition.split(Constants.CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX);
            if (partitionInfo.length != 3) {
                // if  merge-contigs-into-num-partitions was used we may need to use
                // some metadata instead of array names for partition bounds
                return null;
            }
            final String contig = partitionInfo[0];
            final int start = Integer.parseInt(partitionInfo[1]);
            final int end = Integer.parseInt(partitionInfo[2]);
            return new SimpleInterval(contig, start, end);
        }).filter(Objects::nonNull).collect(Collectors.toList());
    }

    private ImportConfig createImportConfig(final int batchSize) {
        final List<GenomicsDBImportConfiguration.Partition> partitions = (intervals == null || intervals.isEmpty()) ? 
                generatePartitionListFromWorkspace() : generatePartitionListFromIntervals();
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder =
                GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();
        importConfigurationBuilder.addAllColumnPartitions(partitions);
        importConfigurationBuilder.setSizePerColumnPartition(vcfBufferSizePerSample);
        importConfigurationBuilder.setFailIfUpdating(!doIncrementalImport);
        importConfigurationBuilder.setSegmentSize(segmentSize);
        importConfigurationBuilder.setConsolidateTiledbArrayAfterLoad(doConsolidation);
        importConfigurationBuilder.setEnableSharedPosixfsOptimizations(sharedPosixFSOptimizations);
        ImportConfig importConfig = new ImportConfig(importConfigurationBuilder.build(), validateSampleToReaderMap, true,
                batchSize, mergedHeaderLines, sampleNameMap.getSampleNameToVcfPath(), bypassFeatureReader ? null : this::createSampleToReaderMap,
                doIncrementalImport);
        importConfig.setOutputCallsetmapJsonFile(callsetMapJSONFile);
        importConfig.setOutputVidmapJsonFile(vidMapJSONFile);
        importConfig.setOutputVcfHeaderFile(vcfHeaderFile);
        importConfig.setUseSamplesInOrder(true);
        importConfig.setFunctionToCallOnBatchCompletion(this::logMessageOnBatchCompletion);
        return importConfig;
    }

    /**
     * A complete traversal from start to finish. This method will import all samples
     * specified in the input GVCF files.
     */
    @Override
    public void traverse() {
        // if we're just trying to get interval list, might need to generate intervals from vid
        // which only gets created with GenomicsDBImporter
        if (getIntervalsFromExistingWorkspace) {
            if (intervals == null || intervals.isEmpty()) {
                intervals = generateIntervalListFromVidMap();
            }
            writeIntervalListToFile();
            return;
        }
        // Force the progress meter to update after every batch
        progressMeter.setRecordsBetweenTimeChecks(1L);

        final int sampleCount = sampleNameMap.getNumSamples();
        final int updatedBatchSize = (batchSize == DEFAULT_ZERO_BATCH_SIZE) ? sampleCount : batchSize;
        final ImportConfig importConfig = createImportConfig(updatedBatchSize);

        GenomicsDBImporter importer;
        try {
            importer = new GenomicsDBImporter(importConfig);
            // Modify importer directly from updateImportProtobufVidMapping.
            GATKGenomicsDBUtils.updateImportProtobufVidMapping(importer);
            if (mergeContigsIntoNumPartitions != 0) {
                if (!doIncrementalImport) {
                    assertIntervalsCoverEntireContigs(importer, intervals);
                    importer.coalesceContigsIntoNumPartitions(mergeContigsIntoNumPartitions);
                }
                else {
                    logger.warn(INCREMENTAL_WORKSPACE_ARG_LONG_NAME+" was set, so ignoring " +
                        MERGE_CONTIGS_INTO_NUM_PARTITIONS + ". When updating workspaces, " +
                        "GenomicsDBImport must use the same partition boundaries/intervals as the original import"); 
                }
            }
            importer.executeImport(maxNumIntervalsToImportInParallel);
        } catch (final IOException e) {
            throw new UserException("Error initializing GenomicsDBImporter", e);
        } catch (final IllegalArgumentException iae) {
            throw new GATKException("Null feature reader found in sampleNameMap file: " + sampleNameMapFile, iae);
        } catch (final CompletionException ce) {
            throw (ce.getCause() instanceof RuntimeException ? (RuntimeException) ce.getCause() : ce);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (getIntervalsFromExistingWorkspace) {
            logger.info("Interval list generated!");
        } else if (batchSize == DEFAULT_ZERO_BATCH_SIZE) {
            logger.info("Import completed!");
        } else {
            logger.info("Import of all batches to GenomicsDB completed!");
        }
        return true;
    }

    /**
     * Method to create feature readers for input files or GCS URLs
     * in the current batch
     *
     * @param sampleNametoPath  Sample name to file name mapping
     * @param batchSize  Current batch size
     * @param lowerSampleIndex  0-based Lower bound of sample index -- inclusive
     * @return  Feature readers to be imported in the current batch, sorted by sample name
     */
    private SortedMap<String, FeatureReader<VariantContext>> getFeatureReadersInParallel(
            final SortedMap<String, URI> sampleNametoPath, final int batchSize, final int lowerSampleIndex) {
        final SortedMap<String, FeatureReader<VariantContext>> sampleToReaderMap = new TreeMap<>();
        logger.info("Starting batch input file preload");
        final Map<String, Future<FeatureReader<VariantContext>>> futures = new LinkedHashMap<>();
        final List<String> sampleNames = new ArrayList<>(sampleNametoPath.keySet());
        for(int i = lowerSampleIndex; i < sampleNametoPath.size() && i < lowerSampleIndex+batchSize; ++i) {
            final String sampleName = sampleNames.get(i);
            futures.put(sampleName, inputPreloadExecutorService.submit(() -> {
                final Path variantPath = IOUtils.getPath(sampleNametoPath.get(sampleName).toString());
                final Path variantIndexPath = sampleNameMap.getVCFIndexForSampleAsPath(sampleName);
                try {
                    return new InitializedQueryWrapper(getReaderFromPath(variantPath, variantIndexPath), intervals.get(0));
                } catch (final IOException e) {
                    throw new UserException.CouldNotReadInputFile("Couldn't read file: " + variantPath.toUri(), e);
                }
            }));
        }

        futures.forEach((sampleName, future) -> {
            try {
                final FeatureReader<VariantContext> reader = future.get();
                sampleToReaderMap.put(sampleName, reader);
            } catch (InterruptedException | ExecutionException e) {
                throw new UserException.CouldNotReadInputFile("Failure while waiting for FeatureReader to initialize ",
                                                              e);
            }
        });
        logger.info("Finished batch preload");
        logger.info("Importing batch " + this.batchCount + " with " + sampleToReaderMap.size() + " samples");
        return sampleToReaderMap;
    }

    private SortedMap<String, FeatureReader<VariantContext>> getFeatureReadersSerially(final Map<String, URI> sampleNameToPath,
                                                                                 final int batchSize, final int lowerSampleIndex){
        final SortedMap<String, FeatureReader<VariantContext>> sampleToReaderMap = new TreeMap<>();
        final List<String> sampleNames = new ArrayList<>(sampleNameToPath.keySet());
        for(int i = lowerSampleIndex; i < sampleNameToPath.size() && i < lowerSampleIndex+batchSize; ++i) {
            final String sampleName = sampleNames.get(i);
            final Path variantPath = IOUtils.getPath(sampleNameToPath.get(sampleName).toString());
            final Path variantIndexPath = sampleNameMap.getVCFIndexForSampleAsPath(sampleName);
            final FeatureReader<VariantContext> reader = getReaderFromPath(variantPath, variantIndexPath);
            sampleToReaderMap.put(sampleName, reader);
        }
        logger.info("Importing batch " + this.batchCount + " with " + sampleToReaderMap.size() + " samples");
        return sampleToReaderMap;
    }

    /**
     * Creates a feature reader object from a given VCF URI (can also be
     * a local file path) and returns it
     * @return  Feature reader
     * @param variantPath
     */
    private FeatureReader<VariantContext> getReaderFromPath(final Path variantPath, final Path variantIndexPath) {
        // TODO: we repeatedly convert between URI, Path, and String in this tool. Is this necessary?
        final String variantURI = variantPath.toAbsolutePath().toUri().toString();
        final String variantIndexURI = variantIndexPath == null ? null : variantIndexPath.toAbsolutePath().toUri().toString();

        try {
            final FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(variantURI, variantIndexURI, new VCFCodec(), true,
                    BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer),
                    BucketUtils.getPrefetchingWrapper(cloudIndexPrefetchBuffer));

            /* Anonymous FeatureReader subclass that wraps returned iterators to ensure that the GVCFs do not
             * contain MNPs.
             */
            return new FeatureReader<>() {
                /** Iterator that asserts that variants are not MNPs. */
                class NoMnpIterator implements CloseableTribbleIterator<VariantContext> {
                    private final CloseableTribbleIterator<VariantContext> inner;
                    NoMnpIterator(CloseableTribbleIterator<VariantContext> inner) { this.inner = inner; }
                    @Override public void close() { inner.close(); }
                    @Override public Iterator<VariantContext> iterator() { return this; }
                    @Override public boolean hasNext() { return inner.hasNext(); }
                    @Override
                    public VariantContext next() {
                        if (!hasNext()) throw new NoSuchElementException();
                        final VariantContext vc = inner.next();
                        if (GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(vc)) {
                            throw new UserException.BadInput(String.format(
                                    "GenomicsDBImport does not support GVCFs with MNPs. MNP found at %1s:%2d in VCF %3s",
                                    vc.getContig(), vc.getStart(), variantPath.toAbsolutePath()
                            ));
                        }

                        return vc;
                    }
                }

                @Override public void close() throws IOException { reader.close(); }
                @Override public List<String> getSequenceNames() { return reader.getSequenceNames(); }
                @Override public Object getHeader() { return reader.getHeader(); }
                @Override public boolean isQueryable() { return reader.isQueryable(); }

                @Override public CloseableTribbleIterator<VariantContext> query(Locatable locus) throws IOException {
                    return new NoMnpIterator(reader.query(locus));
                }
                @Override public CloseableTribbleIterator<VariantContext> query(String chr, int start, int end) throws IOException {
                    return new NoMnpIterator(reader.query(chr, start, end));
                }

                @Override
                public CloseableTribbleIterator<VariantContext> iterator() throws IOException {
                    return new NoMnpIterator(reader.iterator());
                }
            };
        } catch (final TribbleException e){
            throw new UserException("Failed to create reader from " + variantURI + " because of the following error:\n\t"
                    + e.getMessage(), e);
        }
    }

    /**
     * Input argument "overwriteExistingWorkspace" defaults to false.
     * The tool creates a new workspace if it doesn't exist, or checks that
     * a workspace exists in the incremental import case. Deletes
     * an existing workspace if argument is true
     *
     * @return  The workspace directory
     */
    private String overwriteCreateOrCheckWorkspace() {
        String workspaceDir = GATKGenomicsDBUtils.genomicsDBGetAbsolutePath(workspace);
        // From JavaDoc for GATKGenomicsDBUtils.createTileDBWorkspacevid
        //   returnCode = 0 : OK. If overwriteExistingWorkspace is true and the workspace exists, it is deleted first.
        //   returnCode = -1 : path was not a directory
        //   returnCode = -2 : failed to create workspace
        //   returnCode = 1 : if overwriteExistingWorkspace is false, return 1 if directory already exists
        int returnCode = GenomicsDBUtils.createTileDBWorkspace(workspaceDir, overwriteExistingWorkspace);
        if (returnCode == -1) {
            throw new UnableToCreateGenomicsDBWorkspace("Error creating GenomicsDB workspace: " + workspace + " already exists and is not a directory");
        } else if (returnCode < 0) {
            throw new UnableToCreateGenomicsDBWorkspace("Error creating GenomicsDB workspace: " + workspace);
        } else if (!overwriteExistingWorkspace && returnCode == 1) {
            if (doIncrementalImport) {
                return workspaceDir;
            } else {
                throw new UnableToCreateGenomicsDBWorkspace("Error creating GenomicsDB workspace: " + workspace + " already exists");
            }
        } else {
            if (doIncrementalImport) {
                throw new UserException("We require an existing valid workspace when incremental import is set");
            }
            return workspaceDir;
        }
    }

    static class UnableToCreateGenomicsDBWorkspace extends UserException {
        @Serial private static final long serialVersionUID = 1L;

        UnableToCreateGenomicsDBWorkspace(final String message){
            super(message);
        }
    }

    /**
     * Loads our intervals using the best available sequence
     * dictionary (as returned by {@link #getBestAvailableSequenceDictionary})
     * to parse/verify them. Does nothing if no intervals were specified.
     */
    void initializeIntervals() {
        if (intervalArgumentCollection.intervalsSpecified()) {
            if (getIntervalsFromExistingWorkspace || doIncrementalImport) {
                logger.warn(INCREMENTAL_WORKSPACE_ARG_LONG_NAME+" was set, so ignoring specified intervals." +
                    "The tool will use the intervals specified by the initial import");
                intervals = generateIntervalListFromWorkspace();
                return;
            }
            final SAMSequenceDictionary intervalDictionary = getBestAvailableSequenceDictionary();

            if (intervalDictionary == null) {
                throw new UserException("We require at least one input source that " +
                    "has a sequence dictionary (reference or reads) when intervals are specified");
            }

            intervals = new ArrayList<>();
            final List<SimpleInterval> simpleIntervalList = intervalArgumentCollection.getIntervals(intervalDictionary);
            if (!mergeInputIntervals && simpleIntervalList.size() > INTERVAL_LIST_SIZE_WARNING_THRESHOLD) {
                logger.warn(String.format(
                        "A large number of intervals were specified. " +
                        "Using more than %d intervals in a single import is not recommended and can cause performance to suffer. " +
                        "If GVCF data only exists within those intervals, performance can be improved by aggregating intervals with the " +
                        MERGE_INPUT_INTERVALS_LONG_NAME + " argument.",
                        INTERVAL_LIST_SIZE_WARNING_THRESHOLD)
                );
            }
            intervals = mergeInputIntervals ? IntervalUtils.getSpanningIntervals(simpleIntervalList, getBestAvailableSequenceDictionary()) : simpleIntervalList;
        } else if (getIntervalsFromExistingWorkspace || doIncrementalImport) {
            // in incremental import case, we don't care if intervals were specified
            intervals = generateIntervalListFromWorkspace();
        } else {
            throw new UserException("No intervals specified");
        }
    }

    @Override
    public void onShutdown(){
        if(inputPreloadExecutorService != null) {
            inputPreloadExecutorService.shutdownNow();
        }
    }

    /**
     * Overriding getBestAvailableSequenceDictionary() to prefer the mergedVCFHeader's
     * sequence directory, if present, over any other dictionaries
     *
     * @return  Sequence directory from merged header, or super.getBestAvailableSequenceDictionary()
     *          if none
     */
    @Override
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        final SAMSequenceDictionary sequenceDictionary = mergedHeaderSequenceDictionary;
        if (sequenceDictionary == null) {
            return super.getBestAvailableSequenceDictionary();
        } else {
            return sequenceDictionary;
        }
    }

    /**
     * This class is a hack to force parallel loading of the headers and indexes of remote gvcf files.
     * It initializes a feature reader and starts a query.  This causes the header and index to be read, and also causes any
     * pre-fetching to begin if enabled.  It is very narrowly crafted and should not be used for other purposes.
     */
    private static final class InitializedQueryWrapper implements FeatureReader<VariantContext> {
        private final FeatureReader<VariantContext> reader;
        private final SimpleInterval interval;
        private CloseableTribbleIterator<VariantContext> query;

        private InitializedQueryWrapper(final FeatureReader<VariantContext> reader, final Locatable interval) throws IOException {
            this.reader = reader;
            this.interval = new SimpleInterval(interval);
            this.query = reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
        }

        @Override
        public CloseableTribbleIterator<VariantContext> query(final String chr, final int start, final int end) {
            final SimpleInterval queryInterval = new SimpleInterval(chr, start, end);
            if( !interval.equals(queryInterval)){
                throw new GATKException("Cannot call query with different interval, expected:" + this.interval + " queried with: " + queryInterval);
            }
            if( query != null ){
                final CloseableTribbleIterator<VariantContext> tmp = query;
                query = null;
                return tmp;
            } else {
                throw new GATKException("Cannot call query twice on this wrapper.");
            }
        }

        @Override
        public CloseableTribbleIterator<VariantContext> iterator() {
            throw new UnsupportedOperationException("iterator() not supported, this should not have been called and indicates an issue with GenomicsDB integration");
        }

        @Override
        public void close() throws IOException {
            reader.close();
        }

        @Override
        public List<String> getSequenceNames() {
            throw new UnsupportedOperationException("getSequenceNames() not supported, this should not have been called and indicates an issue with GenomicsDB integration");
        }

        @Override
        public Object getHeader() {
            return reader.getHeader();
        }
    }
}

