package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment.SegmentExonUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeter;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.sqlite.util.StringUtils;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGTFFieldConstants.FeatureTag.*;

/**
 * A factory to create {@link GencodeFuncotation}s.
 * This is a high-level object that interfaces with the internals of {@link Funcotator}.
 * Created by jonn on 8/30/17.
 */
public class GencodeFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Public Static Members:
    /**
     * Default name for this data source.
     */
    public static final String DEFAULT_NAME = "Gencode";

    public static final String OTHER_TRANSCRIPTS_INFO_SEP = "_";

    //==================================================================================================================
    // Private Static Members:
    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(GencodeFuncotationFactory.class);

    private static final String LOCAL_GENCODE_TRANSCRIPT_TMP_DIR_PREFIX = "localGencodeTranscriptFastaFolder";
    private static final String LOCAL_GENCODE_TRANSCRIPT_FILE_BASE_NAME = "gencodeTranscriptFastaFile";

    /**
     * Number of bases to the left and right of a variant in which to calculate the GC content.
     */
    private static final int gcContentWindowSizeBases = 200;

    /**
     * The number of leading bases to include when building the variant sequence for UTR variants.
     * Used to determine if there is a de novo start.
     */
    private static final int numLeadingBasesForUtrAnnotationSequenceConstruction = 2;

    /**
     * The number of leading bases to include when building the variant sequence for UTR variants.
     * * Used to determine if there is a de novo start.
     */
    private static final int defaultNumTrailingBasesForUtrAnnotationSequenceConstruction = AminoAcid.CODON_LENGTH;

    /**
     * The window for an indel to be within the end of a transcript to trigger padding the end of the
     * transcript with additional bases from the reference.
     */
    @VisibleForTesting
    static final int TRANSCRIPT_END_WINDOW_PADDING_THRESHOLD = 10;

    /**
     * The window around a variant to include in the reference context annotation.
     * Also used for context from which to get surrounding codon changes and protein changes.
     */
    // TODO: Make this a parameter:
    final static private int referenceWindow = 10;

    /**
     * List of valid Appris Ranks used for sorting funcotations to get the "best" one.z
     */
    private static final LinkedHashMap<String, GencodeGTFFieldConstants.FeatureTag> apprisRanks = new LinkedHashMap<>() {{
            put(APPRIS_PRINCIPAL.toString(), APPRIS_PRINCIPAL);
            put(APPRIS_PRINCIPAL_1.toString(), APPRIS_PRINCIPAL_1);
            put(APPRIS_PRINCIPAL_2.toString(), APPRIS_PRINCIPAL_2);
            put(APPRIS_PRINCIPAL_3.toString(), APPRIS_PRINCIPAL_3);
            put(APPRIS_PRINCIPAL_4.toString(), APPRIS_PRINCIPAL_4);
            put(APPRIS_PRINCIPAL_5.toString(), APPRIS_PRINCIPAL_5);
            put(APPRIS_ALTERNATIVE_1.toString(), APPRIS_ALTERNATIVE_1);
            put(APPRIS_ALTERNATIVE_2.toString(), APPRIS_ALTERNATIVE_2);
            put(APPRIS_CANDIDATE_HIGHEST_SCORE.toString(), APPRIS_CANDIDATE_HIGHEST_SCORE);
            put(APPRIS_CANDIDATE_LONGEST_CCDS.toString(), APPRIS_CANDIDATE_LONGEST_CCDS);
            put(APPRIS_CANDIDATE_CCDS.toString(), APPRIS_CANDIDATE_CCDS);
            put(APPRIS_CANDIDATE_LONGEST_SEQ.toString(), APPRIS_CANDIDATE_LONGEST_SEQ);
            put(APPRIS_CANDIDATE_LONGEST.toString(), APPRIS_CANDIDATE_LONGEST);
            put(APPRIS_CANDIDATE.toString(), APPRIS_CANDIDATE);
        }};

    /**
     * The set of {@link GencodeFuncotation.VariantClassification} types that are valid for coding regions.
     */
    private static final Set<GencodeFuncotation.VariantClassification> codingVariantClassifications =
            Sets.newHashSet(Arrays.asList(GencodeFuncotation.VariantClassification.MISSENSE,
                                          GencodeFuncotation.VariantClassification.NONSENSE,
                                          GencodeFuncotation.VariantClassification.NONSTOP,
                                          GencodeFuncotation.VariantClassification.SILENT,
                                          GencodeFuncotation.VariantClassification.IN_FRAME_DEL,
                                          GencodeFuncotation.VariantClassification.IN_FRAME_INS,
                                          GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS,
                                          GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL,
                                          GencodeFuncotation.VariantClassification.START_CODON_SNP,
                                          GencodeFuncotation.VariantClassification.START_CODON_INS,
                                          GencodeFuncotation.VariantClassification.START_CODON_DEL));
    /** Used when generating the genes field name for segment funcotations. */
    public static final String GENES_SUFFIX = "_genes";
    /** Since gencode funcotation field names are generated dynamically, we have to search the incoming funcotation fields for the proper gencode funcotations.  This suffix (and the others below)
     * help define regular expression patterns for finding the appropriate field names.
     *
     * Used when generating the start gene (gene that overlaps the starting breakpoint of a segment) field name for segment funcotations.
     */
    public static final String START_GENE_SUFFIX = "_start_gene";
    /** Used when generating the end gene (gene that overlaps the ending breakpoint of a segment) field name for segment funcotations.
     * See {@link GencodeFuncotationFactory#START_GENE_SUFFIX} for more information. */
    public static final String END_GENE_SUFFIX = "_end_gene";
    /** Used when generating start exon field name for segment funcotations.  This is the exon that overlaps the starting breakpoint.
     * A "+" means all higher number exons in the coding direction.  And "-" means lower-number.  E.g. "5+"
     * See {@link GencodeFuncotationFactory#START_GENE_SUFFIX} for more information. */
    public static final String START_EXON_SUFFIX = "_start_exon";
    /** Used when generating end exon field name for segment funcotations. This is the exon that overlaps the end breakpoint.
     * A "+" means all higher number exons in the coding direction.  And "-" means lower-number.  E.g. "7+"
     * See {@link GencodeFuncotationFactory#START_GENE_SUFFIX} for more information. */
    public static final String END_EXON_SUFFIX = "_end_exon";

    /** The separator for the genes field when funcotating segments. */
    private static final String GENES_FIELD_SEPARATOR = ",";

    /** The string for the genes field when no genes overlap a segment. */
    private static final String NO_GENE_INFO_STR = "";

    /** For legacy reasons, the alt allele column is always blank in funcotations. */
    private static final String LEGACY_SEGMENT_ALT_ALLELE_VALUE = "";

    /** For legacy reasons, the ref allele column is always blank in funcotations. */
    private static final String LEGACY_SEGMENT_REF_ALLELE_VALUE = "";


    //==================================================================================================================
    // Private Members:

    /**
     * The name of this Gencode data source.
     */
    private final String name;

    /**
     * ReferenceSequenceFile for the transcript reference file.
     */
    private final ReferenceDataSource transcriptFastaReferenceDataSource;

    /**
     * Map between transcript IDs and the IDs from the FASTA file to look up the transcript.
     * This is necessary because of the way the FASTA file contigs are named.
     */
    private final Map<String, MappedTranscriptIdInfo> transcriptIdMap;

    /**
     * The mode to select the "best" transcript (i.e. the transcript with detailed information) from the list of
     * possible transcripts.
     *
     * For more information on the specifics of the differences go here:
     * https://gatkforums.broadinstitute.org/gatk/discussion/4220/what-is-the-difference-between-tx-mode-best-effect-vs-canonical
     */
    private final TranscriptSelectionMode transcriptSelectionMode;

    /**
     * {@link List} of Transcript IDs that the user has requested that we annotate.
     * If this list is empty, will default to keeping ALL transcripts.
     * Otherwise, only transcripts on this list will be annotated.
     */
    private final Set<String> userRequestedTranscripts;

    /**
     * The {@link Path} from which we will read the sequences for the coding regions in given transcripts.
     */
    private final Path gencodeTranscriptFastaFile;

    /**
     * The ncbiBuildVersion for this {@link GencodeFuncotationFactory}.
     * Note: This value is passed in at construction time.
     */
    private String ncbiBuildVersion = null;

    /**
     * The window on either side of splice sites to mark variants as {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#SPLICE_SITE}.
     */
    private int spliceSiteVariantWindowBases = FuncotatorUtils.DEFAULT_SPLICE_SITE_WINDOW_SIZE;

    /**
     * Comparator to be used when sorting {@link Funcotation}s created by this {@link GencodeFuncotationFactory}.
     * Will be either {@link TranscriptSelectionMode.BestEffectGencodeFuncotationComparator} or {@link TranscriptSelectionMode.CanonicalGencodeFuncotationComparator}.
     */
    private final Comparator<GencodeFuncotation> gencodeFuncotationComparator;

    /**
     * Settings object containing our 5'/3' flank sizes
     */
    private final FlankSettings flankSettings;


    private final FuncotationMetadata segmentMetadata;

    /**
     * If this is false, the gencode funcotation factory will not attempt to funcotate segments.  When false, segments
     *  will be annotated with {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#COULD_NOT_DETERMINE}
     */
    private boolean isSegmentFuncotationEnabled;

    /**
     * If this is true, only MANE transcripts will be used for funcotation creation when at least one is present.
     */
    private boolean preferMANETranscripts;

    //==================================================================================================================
    // Constructors:

    /**
     * Creates a {@link GencodeFuncotationFactory} with the 5'/3' flank sizes both set to 0 and will funcotate segments.
     *
     * @param gencodeTranscriptFastaFilePath {@link Path} to the FASTA file containing the sequences of all transcripts in the Gencode data source.
     * @param version The version {@link String} of Gencode from which {@link Funcotation}s will be made.
     * @param name A {@link String} containing the name of this {@link GencodeFuncotationFactory}.
     * @param transcriptSelectionMode The {@link TranscriptSelectionMode} by which representative/verbose transcripts will be chosen for overlapping variants.
     * @param userRequestedTranscripts A {@link Set<String>} containing Gencode TranscriptIDs that the user requests to be annotated with priority over all other transcripts for overlapping variants.
     * @param annotationOverrides A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainFeatureInput The backing {@link FeatureInput} for this {@link GencodeFuncotationFactory}, from which all {@link Funcotation}s will be created.
     * @param ncbiBuildVersion The NCBI build version for this {@link GencodeFuncotationFactory} (can be found in the datasource config file)
     */
    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFilePath,
                                     final String version,
                                     final String name,
                                     final TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts,
                                     final LinkedHashMap<String, String> annotationOverrides,
                                     final FeatureInput<? extends Feature> mainFeatureInput,
                                     final String ncbiBuildVersion) {
        this(gencodeTranscriptFastaFilePath, version, name, transcriptSelectionMode, userRequestedTranscripts, annotationOverrides, mainFeatureInput, new FlankSettings(0, 0), ncbiBuildVersion);
    }

    /**
     * Create a {@link GencodeFuncotationFactory} that will funcotate segments.
     *
     * @param gencodeTranscriptFastaFilePath {@link Path} to the FASTA file containing the sequences of all transcripts in the Gencode data source.
     * @param version The version {@link String} of Gencode from which {@link Funcotation}s will be made.
     * @param name A {@link String} containing the name of this {@link GencodeFuncotationFactory}.
     * @param transcriptSelectionMode The {@link TranscriptSelectionMode} by which representative/verbose transcripts will be chosen for overlapping variants.
     * @param userRequestedTranscripts A {@link Set<String>} containing Gencode TranscriptIDs that the user requests to be annotated with priority over all other transcripts for overlapping variants.
     * @param annotationOverrides A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainFeatureInput The backing {@link FeatureInput} for this {@link GencodeFuncotationFactory}, from which all {@link Funcotation}s will be created.
     * @param flankSettings Settings object containing our 5'/3' flank sizes
     * @param ncbiBuildVersion The NCBI build version for this {@link GencodeFuncotationFactory} (can be found in the datasource config file)
     */
    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFilePath,
                                     final String version,
                                     final String name,
                                     final TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts,
                                     final LinkedHashMap<String, String> annotationOverrides,
                                     final FeatureInput<? extends Feature> mainFeatureInput,
                                     final FlankSettings flankSettings,
                                     final String ncbiBuildVersion) {
        this(gencodeTranscriptFastaFilePath, version, name, transcriptSelectionMode, userRequestedTranscripts, annotationOverrides, mainFeatureInput, flankSettings, false, ncbiBuildVersion, true);
    }

    /**
     * Create a {@link GencodeFuncotationFactory}.
     *
     * @param gencodeTranscriptFastaFilePath {@link Path} to the FASTA file containing the sequences of all transcripts in the Gencode data source.
     * @param version The version {@link String} of Gencode from which {@link Funcotation}s will be made.
     * @param name A {@link String} containing the name of this {@link GencodeFuncotationFactory}.
     * @param transcriptSelectionMode The {@link TranscriptSelectionMode} by which representative/verbose transcripts will be chosen for overlapping variants.
     * @param userRequestedTranscripts A {@link Set<String>} containing Gencode TranscriptIDs that the user requests to be annotated with priority over all other transcripts for overlapping variants.
     * @param annotationOverrides A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainFeatureInput The backing {@link FeatureInput} for this {@link GencodeFuncotationFactory}, from which all {@link Funcotation}s will be created.
     * @param flankSettings Settings object containing our 5'/3' flank sizes
     * @param isDataSourceB37 If {@code true}, indicates that the data source behind this {@link GencodeFuncotationFactory} contains B37 data.
     * @param ncbiBuildVersion The NCBI build version for this {@link GencodeFuncotationFactory} (can be found in the datasource config file)
     */
    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFilePath,
                                     final String version,
                                     final String name,
                                     final TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts,
                                     final LinkedHashMap<String, String> annotationOverrides,
                                     final FeatureInput<? extends Feature> mainFeatureInput,
                                     final FlankSettings flankSettings,
                                     final boolean isDataSourceB37,
                                     final String ncbiBuildVersion,
                                     final boolean isSegmentFuncotationEnabled) {
        this(gencodeTranscriptFastaFilePath, version, name,
                transcriptSelectionMode, userRequestedTranscripts, annotationOverrides, mainFeatureInput,
                flankSettings, isDataSourceB37, ncbiBuildVersion, isSegmentFuncotationEnabled,
                FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT);
    }

        /**
         * Create a {@link GencodeFuncotationFactory}.
         *
         * @param gencodeTranscriptFastaFilePath {@link Path} to the FASTA file containing the sequences of all transcripts in the Gencode data source.
         * @param version The version {@link String} of Gencode from which {@link Funcotation}s will be made.
         * @param name A {@link String} containing the name of this {@link GencodeFuncotationFactory}.
         * @param transcriptSelectionMode The {@link TranscriptSelectionMode} by which representative/verbose transcripts will be chosen for overlapping variants.
         * @param userRequestedTranscripts A {@link Set<String>} containing Gencode TranscriptIDs that the user requests to be annotated with priority over all other transcripts for overlapping variants.
         * @param annotationOverrides A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
         * @param mainFeatureInput The backing {@link FeatureInput} for this {@link GencodeFuncotationFactory}, from which all {@link Funcotation}s will be created.
         * @param flankSettings Settings object containing our 5'/3' flank sizes
         * @param isDataSourceB37 If {@code true}, indicates that the data source behind this {@link GencodeFuncotationFactory} contains B37 data.
         * @param ncbiBuildVersion The NCBI build version for this {@link GencodeFuncotationFactory} (can be found in the datasource config file)
         * @param minBasesForValidSegment The minimum number of bases for a segment to be considered valid.
         */
    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFilePath,
                                     final String version,
                                     final String name,
                                     final TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts,
                                     final LinkedHashMap<String, String> annotationOverrides,
                                     final FeatureInput<? extends Feature> mainFeatureInput,
                                     final FlankSettings flankSettings,
                                     final boolean isDataSourceB37,
                                     final String ncbiBuildVersion,
                                     final boolean isSegmentFuncotationEnabled,
                                     final int minBasesForValidSegment) {

        this(gencodeTranscriptFastaFilePath, version, name, transcriptSelectionMode, userRequestedTranscripts,
                annotationOverrides, mainFeatureInput, flankSettings, isDataSourceB37, ncbiBuildVersion,
                isSegmentFuncotationEnabled, minBasesForValidSegment, FuncotatorUtils.DEFAULT_SPLICE_SITE_WINDOW_SIZE, false);
    }

    /**
     * Create a {@link GencodeFuncotationFactory}.
     *
     * @param gencodeTranscriptFastaFilePath {@link Path} to the FASTA file containing the sequences of all transcripts in the Gencode data source.
     * @param version The version {@link String} of Gencode from which {@link Funcotation}s will be made.
     * @param name A {@link String} containing the name of this {@link GencodeFuncotationFactory}.
     * @param transcriptSelectionMode The {@link TranscriptSelectionMode} by which representative/verbose transcripts will be chosen for overlapping variants.
     * @param userRequestedTranscripts A {@link Set<String>} containing Gencode TranscriptIDs that the user requests to be annotated with priority over all other transcripts for overlapping variants.
     * @param annotationOverrides A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainFeatureInput The backing {@link FeatureInput} for this {@link GencodeFuncotationFactory}, from which all {@link Funcotation}s will be created.
     * @param flankSettings Settings object containing our 5'/3' flank sizes
     * @param isDataSourceB37 If {@code true}, indicates that the data source behind this {@link GencodeFuncotationFactory} contains B37 data.
     * @param ncbiBuildVersion The NCBI build version for this {@link GencodeFuncotationFactory} (can be found in the datasource config file)
     * @param minBasesForValidSegment The minimum number of bases for a segment to be considered valid.
     * @param spliceSiteWindowSize The number of bases on either side of a splice site for a variant to be a {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#SPLICE_SITE} variant.
     */
    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFilePath,
                                     final String version,
                                     final String name,
                                     final TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts,
                                     final LinkedHashMap<String, String> annotationOverrides,
                                     final FeatureInput<? extends Feature> mainFeatureInput,
                                     final FlankSettings flankSettings,
                                     final boolean isDataSourceB37,
                                     final String ncbiBuildVersion,
                                     final boolean isSegmentFuncotationEnabled,
                                     final int minBasesForValidSegment,
                                     final int spliceSiteWindowSize,
                                     final boolean preferMANETranscriptsWhereApplicable) {

        super(mainFeatureInput, minBasesForValidSegment);

        // Set up our local transcript fasta file.
        // We must localize it (if not on disk) to make read times fast enough to be manageable:
        gencodeTranscriptFastaFile = localizeGencodeTranscriptFastaFile( gencodeTranscriptFastaFilePath );
        this.flankSettings = flankSettings;

        // Initialize our transcript data source and ID map:
        transcriptFastaReferenceDataSource = ReferenceDataSource.of(gencodeTranscriptFastaFile);
        transcriptIdMap = createTranscriptIdMap(transcriptFastaReferenceDataSource);

        this.transcriptSelectionMode = transcriptSelectionMode;

        this.version = version;

        this.name = name;

        this.dataSourceIsB37 = isDataSourceB37;

        this.ncbiBuildVersion = ncbiBuildVersion;

        this.isSegmentFuncotationEnabled = isSegmentFuncotationEnabled;

        this.spliceSiteVariantWindowBases = spliceSiteWindowSize;

        // Go through each requested transcript and remove the version numbers from them if they exist:
        this.userRequestedTranscripts = new HashSet<>();
        for ( final String transcript : userRequestedTranscripts ) {
            this.userRequestedTranscripts.add( FuncotatorUtils.getTranscriptIdWithoutVersionNumber(transcript) );
        }

        // Set our comparator for outputting our funcotations in the right order with the correct "best" transcript:
        gencodeFuncotationComparator = transcriptSelectionMode.getComparator(userRequestedTranscripts);

        // Initialize segment metadata.  This will only be used if funcotating segments.  If we are not doing that create
        //  an empty set of metadata for segments.
        // TODO: Refactor this away as part of https://github.com/broadinstitute/gatk/issues/5932
        segmentMetadata = isSegmentFuncotationEnabled? createSegmentFuncotationMetadata() :
                FuncotationMetadataUtils.createWithUnknownAttributes(Collections.emptyList());

        // Initialize overrides / defaults:
        initializeAnnotationOverrides( annotationOverrides );

        this.preferMANETranscripts = preferMANETranscriptsWhereApplicable;
    }

    private Path localizeGencodeTranscriptFastaFile( final Path gencodeTranscriptFastaFilePath ) {

        // Is the path local or in the cloud:
        if ( gencodeTranscriptFastaFilePath.getFileSystem().equals(FileSystems.getDefault()) ) {
            // local path, just return it:
            return gencodeTranscriptFastaFilePath;
        }

        // Not a local path!  We must localize it!

        // Get the remote paths for the index and dictionary files:
        final Path remoteGencodeTranscriptFastaIndexFilePath = IOUtils.getPath( ReferenceUtils.getFastaIndexFileName(gencodeTranscriptFastaFilePath.toUri().toString()) );
        final Path remoteGencodeTranscriptFastaSequenceDictionaryFilePath = IOUtils.getPath( ReferenceUtils.getFastaDictionaryFileName(gencodeTranscriptFastaFilePath.toUri().toString()) );

        // Create a place for the files:
        final File tmpDir = IOUtils.createTempDir(LOCAL_GENCODE_TRANSCRIPT_TMP_DIR_PREFIX);
        tmpDir.deleteOnExit();
        final Path tmpDirPath = tmpDir.toPath();

        // Create paths to the fasta, fasta index, and the sequence dictionary:
        final Path localGencodeTranscriptFastaFilePath = tmpDirPath.resolve(LOCAL_GENCODE_TRANSCRIPT_FILE_BASE_NAME + ".fa");
        final Path localGencodeTranscriptFastaIndexFilePath = IOUtils.getPath( ReferenceUtils.getFastaIndexFileName(localGencodeTranscriptFastaFilePath.toUri().toString()) );
        final Path localGencodeTranscriptFastaSequenceDictionaryFilePath = IOUtils.getPath( ReferenceUtils.getFastaDictionaryFileName(localGencodeTranscriptFastaFilePath.toUri().toString()) );

        // Copy the files to our local machine:
        logger.info("Localizing Gencode transcript FASTA file for faster lookup times...");

        // Copy FASTA:
        NioFileCopierWithProgressMeter.create(gencodeTranscriptFastaFilePath, localGencodeTranscriptFastaFilePath, true).initiateCopy();

        // Copy Index:
        NioFileCopierWithProgressMeter.create(remoteGencodeTranscriptFastaIndexFilePath, localGencodeTranscriptFastaIndexFilePath, true).initiateCopy();

        // Copy Sequence Dictionary:
        NioFileCopierWithProgressMeter.create(remoteGencodeTranscriptFastaSequenceDictionaryFilePath, localGencodeTranscriptFastaSequenceDictionaryFilePath, true).initiateCopy();

        // Bye Bye!
        return localGencodeTranscriptFastaFilePath;
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public Class<? extends Feature> getAnnotationFeatureClass() {
        return GencodeGtfFeature.class;
    }

    @Override
    public String getInfoString() {
        return getName() + " " + getVersion() + " " + transcriptSelectionMode.toString();
    }

    @Override
    public void close() {
        transcriptFastaReferenceDataSource.close();
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {

            return new LinkedHashSet<>(Arrays.asList(
                    getName() + "_" + getVersion() + "_hugoSymbol",
                    getName() + "_" + getVersion() + "_ncbiBuild",
                    getName() + "_" + getVersion() + "_chromosome",
                    getName() + "_" + getVersion() + "_start",
                    getName() + "_" + getVersion() + "_end",
                    getName() + "_" + getVersion() + "_variantClassification",
                    getName() + "_" + getVersion() + "_secondaryVariantClassification",
                    getName() + "_" + getVersion() + "_variantType",
                    getName() + "_" + getVersion() + "_refAllele",
                    getName() + "_" + getVersion() + "_tumorSeqAllele1",
                    getName() + "_" + getVersion() + "_tumorSeqAllele2",
                    getName() + "_" + getVersion() + "_genomeChange",
                    getName() + "_" + getVersion() + "_annotationTranscript",
                    getName() + "_" + getVersion() + "_transcriptStrand",
                    getName() + "_" + getVersion() + "_transcriptExon",
                    getName() + "_" + getVersion() + "_transcriptPos",
                    getName() + "_" + getVersion() + "_cDnaChange",
                    getName() + "_" + getVersion() + "_codonChange",
                    getName() + "_" + getVersion() + "_proteinChange",
                    getName() + "_" + getVersion() + "_gcContent",
                    getName() + "_" + getVersion() + "_referenceContext",
                    getName() + "_" + getVersion() + "_otherTranscripts"
            ));
    }

    @Override
    protected List<Funcotation> createDefaultFuncotationsOnVariant( final VariantContext variant, final ReferenceContext referenceContext ) {

        if (FuncotatorUtils.isSegmentVariantContext(variant, minBasesForValidSegment)) {
            return createSegmentFuncotations(variant, Collections.emptyList(), null, null, null, null);
        } else {
            // Simply create IGR
            final List<Funcotation> funcotationList = new ArrayList<>();
            funcotationList.addAll(createIgrFuncotations(variant, referenceContext));
            return funcotationList;
        }
    }

    private List<GencodeFuncotation> createGencodeFuncotationsByAllTranscripts( final VariantContext variant, final ReferenceContext referenceContext, final Allele altAllele, final List<GencodeGtfGeneFeature> gencodeGtfGeneFeatures) {

        // If the variant overlaps more than one gene, we need to create a flat list of the transcripts in all genes.
        final List<GencodeFuncotation> gencodeFuncotationList = gencodeGtfGeneFeatures.stream()
                .map(f -> createFuncotationsHelper(variant, altAllele, f, referenceContext))
                .flatMap(List::stream).collect(Collectors.toList());

        // Sort the funcotations, and populate other transcripts. Note that we're guaranteed not to have any
        // IGR funcotations at this point due to the contract of createFuncotationsHelper().
        if (gencodeFuncotationList.size() > 0) {
            // Get our "Best Transcript" from our list.
            sortFuncotationsByTranscriptForOutput(gencodeFuncotationList);
            populateOtherTranscriptsMapForFuncotation(gencodeFuncotationList);
        }

        // Grab the best choice in the case of transcript selection modes other than ALL.  The selection will be the first
        // transcript in the list.
        if ((this.transcriptSelectionMode != TranscriptSelectionMode.ALL) && (gencodeFuncotationList.size() > 0)) {
            return Collections.singletonList(gencodeFuncotationList.get(0));
        }

        return gencodeFuncotationList;
    }

    private void populateOtherTranscriptsMapForFuncotation(final List<GencodeFuncotation> gencodeFuncotations) {
        // This method assumes unique transcript names in each GencodeFuncotation.  Do not use Function.identity() as
        //  the key or error will result.
        // First create a map that goes from each given funcotation to its condensed version.
        final Map<String, String> funcotationToCondensedString = gencodeFuncotations.stream()
                .collect(Collectors.toMap(f -> f.getAnnotationTranscript(), GencodeFuncotationFactory::condenseGencodeFuncotation, (x,y) -> y, LinkedHashMap::new));

        for (final GencodeFuncotation g: gencodeFuncotations) {
            final List<String> otherTranscriptStrings = funcotationToCondensedString.keySet().stream()
                .filter(f -> !f.equals(g.getAnnotationTranscript()))
                .map(f -> funcotationToCondensedString.get(f))
                .collect(Collectors.toList());
            g.setOtherTranscripts(otherTranscriptStrings);
        }
    }

    /**
     * Attempts to treat the given features as {@link GencodeGtfFeature} objects in order to
     * create funcotations for the given variant and reference.
     *
     * This is the entry point into the Factory
     */
    @Override
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> geneFeatureList) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        final List<Feature> nonNullFeatures = geneFeatureList.stream().filter(g -> g != null).collect(Collectors.toList());

        // If we have features we need to annotate, go through them and create annotations:
        if ( nonNullFeatures.size() > 0 ) {
            for ( final Allele altAllele : variant.getAlternateAlleles() ) {

                // At this point we know the feature list is composed of GTF Gene Features
                final List<GencodeGtfGeneFeature> gencodeGtfGeneFeatures = convertFeaturesToGencodeGtfGeneFeatures(nonNullFeatures);

                // By this point we know the feature type is correct, so we cast it:
                final List<GencodeFuncotation> gencodeFuncotationList = createGencodeFuncotationsByAllTranscripts(variant, referenceContext, altAllele, gencodeGtfGeneFeatures);

                // Add the returned funcotations here:
                if ( ! gencodeFuncotationList.isEmpty() ) {
                    outputFuncotations.addAll(gencodeFuncotationList);
                } else {
                    // If the funcotation List is empty, it means that this variant was IGR for all transcripts,
                    // or that there were no basic transcripts:
                    outputFuncotations.add(createIgrFuncotation(variant, altAllele, referenceContext));
                }
            }
        } else {
            // This is an IGR. We have to have an annotation on all variants, so we must annotate it.
            outputFuncotations.addAll(createIgrFuncotations(variant, referenceContext));
        }

        return outputFuncotations;
    }

    private static List<GencodeGtfGeneFeature> convertFeaturesToGencodeGtfGeneFeatures(final List<Feature> nonNullFeatures) {
        return nonNullFeatures.stream()
                        .filter(g -> g != null)
                        .map(g -> (GencodeGtfGeneFeature) g)
                        .collect(Collectors.toList());
    }

    /**
     * If MANE_Plus_Clinical transcripts are avalible, only return them, followed by MANE_Select transcripts, followed by only the basic transcripts if none were MANE_Plus_Clinical or MANE_Select.
     * @param transcripts of gencode transcripts to possibly filter
     * @return
     */
    @VisibleForTesting
    static List<GencodeGtfTranscriptFeature> retreiveMANESelectModeTranscriptsCriteria(final List<GencodeGtfTranscriptFeature> transcripts) {
        final List<GencodeGtfTranscriptFeature> plusClincal = transcripts.stream()
                .filter(g -> hasTag(g, MANE_PLUS_CLINICAL)).toList();
        if (plusClincal.size() > 0) {
            return plusClincal;
        }

        final List<GencodeGtfTranscriptFeature> maneSelectTranscripts = transcripts.stream()
                .filter(g -> hasTag(g, MANE_SELECT)).toList();

        if (maneSelectTranscripts.size() > 0) {
            return maneSelectTranscripts;
        }

        return transcripts.stream().filter(GencodeFuncotationFactory::isBasic).collect(Collectors.toList());
    }

    /**
     * {@inheritDoc}
     * This method should never be called on a {@link GencodeFuncotationFactory}, as that would imply a time-loop.
     */
    @Override
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        throw new GATKException("This method should never be called on a "+ this.getClass().getName());
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.GENCODE;
    }

    /**
     * {@inheritDoc}
     *
     * We override this method to request extra padding around queries on our FeatureContext to take into account
     * the 5' and 3' flank sizes in our {@link #flankSettings}. We use the max of these two values as our
     * query padding.
     */
    @Override
    protected SimpleInterval transformFeatureQueryInterval(final SimpleInterval queryInterval) {
        final int queryPadding = Math.max(flankSettings.fivePrimeFlankSize, flankSettings.threePrimeFlankSize);

        final int newStart = queryInterval.getStart() - queryPadding;
        final int newEnd   = queryInterval.getEnd() + queryPadding;

        return new SimpleInterval( queryInterval.getContig(), newStart < 1 ? 1 : newStart, newEnd);
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Filter the given list of {@link GencodeFuncotation} to only contain those funcotations that have transcriptIDs that
     * appear in the given {@code selectedTranscripts}.
     * Ignores transcript version numbers.
     * @param funcotations The {@link List} of {@link GencodeFuncotation} to filter.
     * @param selectedTranscripts The {@link Set} of transcript IDs to keep in the given {@code funcotations}.
     */
    static void filterAnnotationsByUserTranscripts( final List<GencodeFuncotation> funcotations,
                                                    final Set<String> selectedTranscripts ) {
        if ( (selectedTranscripts != null) && (!selectedTranscripts.isEmpty()) ) {
            funcotations.removeIf( f -> !FuncotatorUtils.isFuncotationInTranscriptList(f, selectedTranscripts) );
        }
    }

    /**
     * Creates a map of Transcript IDs for use in looking up transcripts from the FASTA dictionary for the GENCODE Transcripts.
     * We include the start and stop codons in the transcripts so we can handle start/stop codon variants.
     * @param fastaReference The {@link ReferenceDataSource} corresponding to the Transcript FASTA file for this GENCODE dataset.
     * @return A {@link Map} of {@link String} -> {@link MappedTranscriptIdInfo} which maps real transcript IDs to the information about that transcript in the transcript FASTA file.
     */
    @VisibleForTesting
    static Map<String, MappedTranscriptIdInfo> createTranscriptIdMap(final ReferenceDataSource fastaReference) {

        final Map<String, MappedTranscriptIdInfo> idMap = new HashMap<>();

        for ( final SAMSequenceRecord sequence : fastaReference.getSequenceDictionary().getSequences() ) {

            final MappedTranscriptIdInfo transcriptInfo = createMappedTranscriptIdInfo( sequence );

            // The names in the file are actually in a list with | between each sequence name.
            // We need to split the names and add them to the dictionary so we can resolve them to the full
            // sequence name as it appears in the file:
            for ( final String transcriptId : Utils.split(sequence.getSequenceName(), "|") ) {
                idMap.put(transcriptId, transcriptInfo);
            }
        }

        return idMap;
    }

    /**
     * Creates a {@link MappedTranscriptIdInfo} object based on the given {@link SAMSequenceRecord}.
     * This method is a helper method to get information out of a GENCODE transcript FASTA file easily.
     * This method assumes that {@code sequence} is a sequence from a GENCODE transcript FASTA file.
     * @param sequence The {@link SAMSequenceRecord} from which to create the {@link MappedTranscriptIdInfo}.
     * @return A populated {@link MappedTranscriptIdInfo} object based on the given {@link SAMSequenceRecord}.
     */
    private static MappedTranscriptIdInfo createMappedTranscriptIdInfo( final SAMSequenceRecord sequence ) {

        final MappedTranscriptIdInfo transcriptIdInfo = new MappedTranscriptIdInfo();

        final Pattern utrPattern = Pattern.compile("UTR[35]:(\\d+)-(\\d+)");
        final Pattern cdsPattern = Pattern.compile("CDS:(\\d+)-(\\d+)");

        boolean has3pUtr = false;
        boolean has5pUtr = false;

        // Now let's go through the sequence name and pull out the salient features for each field:
        for (final String field : Utils.split(sequence.getSequenceName(), "|")) {
            if ((field.length() > 4) && (field.substring(0, 5).equals("UTR5:"))) {
                final Matcher m = utrPattern.matcher(field);
                m.find();
                transcriptIdInfo.fivePrimeUtrStart = Integer.valueOf(m.group(1));
                transcriptIdInfo.fivePrimeUtrEnd = Integer.valueOf(m.group(2));
                has5pUtr = true;
            } else if ((field.length() > 4) && (field.substring(0, 5).equals("UTR3:"))) {
                final Matcher m = utrPattern.matcher(field);
                m.find();
                transcriptIdInfo.threePrimeUtrStart = Integer.valueOf(m.group(1));
                transcriptIdInfo.threePrimeUtrEnd = Integer.valueOf(m.group(2));
                has3pUtr = true;
            } else if ((field.length() > 3) && (field.substring(0, 4).equals("CDS:"))) {
                final Matcher m = cdsPattern.matcher(field);
                m.find();
                transcriptIdInfo.codingSequenceStart = Integer.valueOf(m.group(1));
                transcriptIdInfo.codingSequenceEnd = Integer.valueOf(m.group(2));
            }
        }

        // TODO: This may not be correct.  It's not clear that the whole sequence should be taken.  Perhaps it should be handled with logic downstream...
        // If no CDS was specified, we use the whole sequence:
        if ( transcriptIdInfo.codingSequenceStart == 0 ) {
            transcriptIdInfo.codingSequenceStart = 1;
        }
        if ( transcriptIdInfo.codingSequenceEnd == 0 ) {
            transcriptIdInfo.codingSequenceEnd = sequence.getSequenceLength();
        }

        transcriptIdInfo.mapKey = sequence.getSequenceName();
        transcriptIdInfo.has3pUtr = has3pUtr;
        transcriptIdInfo.has5pUtr = has5pUtr;

        return transcriptIdInfo;
    }

    /**
     * Get the coding sequence from the GENCODE Transcript FASTA file for a given {@code transcriptId}.
     * This will get ONLY the coding sequence for the given {@code transcriptId} and will not include any UTRs.
     * @param transcriptId The ID of the transcript to get from the FASTA file.
     * @param transcriptIdMap A map from transcriptId to MappedTranscriptIdInfo, which tells us how to pull information for the given {@code transcriptId} out of the given {@code transcriptFastaReferenceDataSource}.
     * @param transcriptFastaReferenceDataSource A {@link ReferenceDataSource} for the GENCODE transcript FASTA file.
     * @param transcriptTailPaddingBaseString Bases to add to the end of the transcript base string to enable processing variants that overrrun the end of the transcript.
     * @return The coding sequence for the given {@code transcriptId} as represented in the GENCODE transcript FASTA file.
     */
    private static String getCodingSequenceFromTranscriptFasta( final String transcriptId,
                                                                final Map<String, MappedTranscriptIdInfo> transcriptIdMap,
                                                                final ReferenceDataSource transcriptFastaReferenceDataSource,
                                                                final String transcriptTailPaddingBaseString) {

        final MappedTranscriptIdInfo transcriptMapIdAndMetadata = transcriptIdMap.get(transcriptId);

        if ( transcriptMapIdAndMetadata == null ) {
            throw new UserException.BadInput( "Unable to find the given Transcript ID in our transcript list for our coding sequence (not in given transcript FASTA file): " + transcriptId );
        }

        final SimpleInterval transcriptInterval = new SimpleInterval(
                transcriptMapIdAndMetadata.mapKey,
                transcriptMapIdAndMetadata.codingSequenceStart,
                transcriptMapIdAndMetadata.codingSequenceEnd
        );

        return transcriptFastaReferenceDataSource.queryAndPrefetch( transcriptInterval ).getBaseString() + transcriptTailPaddingBaseString;
    }

    /**
     * Get the 5' UTR sequence from the GENCODE Transcript FASTA file for a given {@code transcriptId}.
     * This will get ONLY the 5' UTR sequence for the given {@code transcriptId} and will NOT include the coding sequence or the 3' UTR.
     * If the given transcript has no 5' UTR, this will return an empty {@link String}.
     * @param transcriptId The ID of the transcript to get from the FASTA file.
     * @param transcriptIdMap A map from transcriptId to MappedTranscriptIdInfo, which tells us how to pull information for the given {@code transcriptId} out of the given {@code transcriptFastaReferenceDataSource}.
     * @param transcriptFastaReferenceDataSource A {@link ReferenceDataSource} for the GENCODE transcript FASTA file.
     * @param extraBases The number of extra bases from the coding region to include in the results after the 5' UTR region.
     * @return The coding sequence for the given {@code transcriptId} as represented in the GENCODE transcript FASTA file, or an empty {@link String}.
     */
    private static String getFivePrimeUtrSequenceFromTranscriptFasta( final String transcriptId,
                                                                      final Map<String, MappedTranscriptIdInfo> transcriptIdMap,
                                                                      final ReferenceDataSource transcriptFastaReferenceDataSource,
                                                                      final int extraBases, 
                                                                      final int transcriptLength ) {

        final MappedTranscriptIdInfo transcriptMapIdAndMetadata = transcriptIdMap.get(transcriptId);

        if ( transcriptMapIdAndMetadata == null ) {
            throw new UserException.BadInput( "Unable to find the given Transcript ID in our transcript list for our 5'UTR (not in given transcript FASTA file): " + transcriptId );
        }

	boolean intervalWithinTranscript = transcriptMapIdAndMetadata.fivePrimeUtrEnd + extraBases < transcriptLength;

        if ( transcriptMapIdAndMetadata.has5pUtr && intervalWithinTranscript ) {

            final SimpleInterval transcriptInterval = new SimpleInterval(
                    transcriptMapIdAndMetadata.mapKey,
                    transcriptMapIdAndMetadata.fivePrimeUtrStart,
                    transcriptMapIdAndMetadata.fivePrimeUtrEnd + extraBases
            );

            return transcriptFastaReferenceDataSource.queryAndPrefetch(transcriptInterval).getBaseString();
        }
        else {
            return "";
        }
    }

    /**
     * Returns whether a variant is in a coding region based on its primary and secondary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification}.
     * @param varClass Primary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} of a variant.
     * @param secondaryVarClass Secondary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} of a variant.
     * @return {@code true} if the corresponding variant is in a coding region; {@code false} otherwise.
     */
    @VisibleForTesting
    static boolean isVariantInCodingRegion(final GencodeFuncotation.VariantClassification varClass,
                                           final GencodeFuncotation.VariantClassification secondaryVarClass ) {

        Utils.nonNull( varClass );

        if ( varClass == GencodeFuncotation.VariantClassification.SPLICE_SITE ) {
            Utils.nonNull( secondaryVarClass );
            return secondaryVarClass != GencodeFuncotation.VariantClassification.INTRON;
        }
        else {
            return codingVariantClassifications.contains(varClass);
        }
    }

    //==================================================================================================================
    // Normal Methods:

    /**
     * Creates a {@link List} of {@link GencodeFuncotation}s based on the given {@link VariantContext}, {@link Allele}, and {@link GencodeGtfGeneFeature}.
     * Will return an empty List if the variant was in the IGR for all transcripts.
     *
     * @param variant The variant to annotate.
     * @param altAllele The allele of the given variant to annotate.
     * @param gtfFeature The GTF feature on which to base annotations.
     * @return A {@link List} of {@link GencodeFuncotation}s for the given variant, allele and gtf feature.
     *         Will be an empty List if the variant was in the IGR for all transcripts.
     */
    private List<GencodeFuncotation> createFuncotationsHelper(final VariantContext variant, final Allele altAllele, final GencodeGtfGeneFeature gtfFeature, final ReferenceContext reference) {

        List<GencodeGtfTranscriptFeature> transcriptList;

        // Only get basic transcripts if we're using data from Gencode:
        if ( gtfFeature.getGtfSourceFileType().equals(GencodeGtfCodec.GTF_FILE_TYPE_STRING) ) {
            if (preferMANETranscripts) {
                // Filter out the non-MANE_Select/Mane_Plus_Clinical transcripts if we're only using MANE transcripts:
                transcriptList = retreiveMANESelectModeTranscriptsCriteria(gtfFeature.getTranscripts());
            } else {
                transcriptList = retrieveBasicTranscripts(gtfFeature);
            }
        } else {
            transcriptList = gtfFeature.getTranscripts();
        }


        return createFuncotationsHelper(variant, altAllele, reference, transcriptList);
    }

    private static List<GencodeGtfTranscriptFeature> retrieveBasicTranscripts(final GencodeGtfGeneFeature gtfFeature) {
        return gtfFeature.getTranscripts().stream()
                .filter(GencodeFuncotationFactory::isBasic).collect(Collectors.toList());
    }

    private List<GencodeFuncotation> createFuncotationsHelper(final VariantContext variant, final Allele altAllele, final ReferenceContext reference, final List<GencodeGtfTranscriptFeature> basicTranscripts) {
        final List<GencodeFuncotation> outputFuncotations = new ArrayList<>();

        // For each applicable transcript, create an annotation.
        for ( final GencodeGtfTranscriptFeature transcript : basicTranscripts ) {

            // Try to create the annotation:
            try {
                final GencodeFuncotation gencodeFuncotation = createGencodeFuncotationOnSingleTranscript(variant, altAllele, reference, transcript);

                // Add the functotation for this transcript into our output funcotations. It will be null if this was an IGR.
                if ( gencodeFuncotation != null ) {
                    outputFuncotations.add(gencodeFuncotation);
                }
            }
            catch ( final FuncotatorUtils.TranscriptCodingSequenceException ex ) {
                 logger.error("Problem creating a GencodeFuncotation on transcript " + transcript.getTranscriptId() + " for variant: " +
                        variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + "(" + variant.getReference() + " -> " + altAllele + "): " +
                         ex.getMessage()
                 );
                logger.warn("Creating default GencodeFuncotation on transcript " + transcript.getTranscriptId() + " for problem variant: " +
                                variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + "(" + variant.getReference() + " -> " + altAllele + ")");
                outputFuncotations.add( createDefaultFuncotationsOnProblemVariant( variant, altAllele, reference, transcript, version, getName() ) );
            }
        }
        return outputFuncotations;
    }

    /**
     * Create a placeholder funcotation on a given {@code variant} and {@code allele} pair for a case that funcotator
     * cannot yet handle, or would currently get wrong.
     * Primarily this occurs when a variant is long and spans multiple types of {@link GencodeGtfFeature}s
     * (i.e. it starts in an intron and ends in an exon or visa-versa). or is
     * long and begins in a transcript and extends beyond a given transcript's end point.
     * There are two such cases right now as manifested in the following issues:
     *     https://github.com/broadinstitute/gatk/issues/3749
     *     https://github.com/broadinstitute/gatk/issues/4307
     * As noted in the above issues, other functional annotation tools also get these kinds of cases wrong.
     * @param variant The {@link VariantContext} to annotate.
     * @param altAllele The alternate {@link Allele} to annotate.
     * @param reference The {@link ReferenceContext} for the given {@code variant}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} which is being used to annotate the given {@code variant}.
     * @param version A {@link String} representing the version of the {@link GencodeFuncotationFactory} being used to annotate the given {@code variant}.
     * @param dataSourceName A {@link String} containing the name of the data source instance.
     * @return A placeholder {@link GencodeFuncotation} for the given {@code variant}.
     */
    private GencodeFuncotation createDefaultFuncotationsOnProblemVariant( final VariantContext variant,
                                                                               final Allele altAllele,
                                                                               final ReferenceContext reference,
                                                                               final GencodeGtfTranscriptFeature transcript,
                                                                               final String version,
                                                                               final String dataSourceName) {
        return createDefaultFuncotationsOnProblemVariant(variant, altAllele, reference, transcript,
                version, dataSourceName, this.ncbiBuildVersion);
    }

    /**
     * Create a placeholder funcotation on a given {@code variant} and {@code allele} pair for a case that funcotator
     * cannot yet handle, or would currently get wrong.
     * Primarily this occurs when a variant is long and spans multiple types of {@link GencodeGtfFeature}s
     * (i.e. it starts in an intron and ends in an exon or visa-versa). or is
     * long and begins in a transcript and extends beyond a given transcript's end point.
     * There are two such cases right now as manifested in the following issues:
     *     https://github.com/broadinstitute/gatk/issues/3749
     *     https://github.com/broadinstitute/gatk/issues/4307
     * As noted in the above issues, other functional annotation tools also get these kinds of cases wrong.
     * @param variant The {@link VariantContext} to annotate.
     * @param altAllele The alternate {@link Allele} to annotate.
     * @param reference The {@link ReferenceContext} for the given {@code variant}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} which is being used to annotate the given {@code variant}.
     * @param version A {@link String} representing the version of the {@link GencodeFuncotationFactory} being used to annotate the given {@code variant}.
     * @param dataSourceName A {@link String} containing the name of the data source instance.
     * @param ncbiBuildVersion NCBI build version
     * @return A placeholder {@link GencodeFuncotation} for the given {@code variant}.
     */
    @VisibleForTesting
    static final GencodeFuncotation createDefaultFuncotationsOnProblemVariant( final VariantContext variant,
                                                                               final Allele altAllele,
                                                                               final ReferenceContext reference,
                                                                               final GencodeGtfTranscriptFeature transcript,
                                                                               final String version,
                                                                               final String dataSourceName,
                                                                               final String ncbiBuildVersion) {
        // Create basic annotation information:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, transcript, ncbiBuildVersion);

        // Set our version:
        gencodeFuncotationBuilder.setVersion(version);

        // Get the reference bases for our current variant:
        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, reference, transcript.getGenomicStrand(), referenceWindow);

        // Set the reference context with the bases from the sequence comparison
        // NOTE: The reference context is ALWAYS from the + strand, so we need to reverse our bases back in the - case:
        gencodeFuncotationBuilder.setReferenceContext(referenceBases.getBaseString(Strand.POSITIVE));

        // Get GC Content:
        gencodeFuncotationBuilder.setGcContent( calculateGcContent( variant.getReference(), altAllele, reference, gcContentWindowSizeBases ) );

        gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE);

        gencodeFuncotationBuilder.setDataSourceName(dataSourceName);

        return gencodeFuncotationBuilder.build();
    }

    private static boolean isBasic(final GencodeGtfTranscriptFeature transcript) {
        // Check if this transcript has the `basic` tag:
        return hasTag(transcript, GencodeGTFFieldConstants.FeatureTag.BASIC);
    }

    private static boolean hasTag(final GencodeGtfTranscriptFeature transcript, final GencodeGTFFieldConstants.FeatureTag tag) {
        // Check if this transcript has the given tag:
        return transcript.getOptionalFields().stream()
                .filter( f -> f.getName().equals("tag") )
                .filter( f -> f.getValue().equals(tag.toString()) )
                .count() > 0;
    }

    /**
     * Create a {@link GencodeFuncotation} for a given variant and transcript. Returns null for IGR variants.
     *
     * @param variant The {@link VariantContext} to annotate.
     * @param altAllele The alternate {@link Allele} to annotate.
     * @param reference The {@link ReferenceContext} for the given {@link VariantContext}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation}, or null if the variant was IGR.
     */
    private GencodeFuncotation createGencodeFuncotationOnSingleTranscript(final VariantContext variant,
                                                                          final Allele altAllele,
                                                                          final ReferenceContext reference,
                                                                          final GencodeGtfTranscriptFeature transcript) {
        final GencodeFuncotation gencodeFuncotation;

        // TODO: check for complex INDEL and warn and skip (https://github.com/broadinstitute/gatk/issues/5411).

        // Find the sub-feature of transcript that contains our variant:
        final GencodeGtfFeature containingSubfeature = getContainingGtfSubfeature(variant, transcript);

        // First handle flanking and IGR variants
        if ( containingSubfeature == null ) {

            // TODO: If the transcript overlaps the variant, but getContainingGtfSubfeature() returned null, it means
            // TODO: that the variant spans across a transcript boundary, such as the boundary between the UTR and the
            // TODO: flank region. We cannot currently handle this case, so we throw an exception to cause us to annotate
            // TODO: with VariantClassification.COULD_NOT_DETERMINE. See https://github.com/broadinstitute/gatk/issues/4307
            if ( transcript.overlaps(variant) ) {
                throw new FuncotatorUtils.TranscriptCodingSequenceException(String.format(
                        "Variant overlaps transcript but is not completely contained within it. " +
                        "Funcotator cannot currently handle this case. Transcript: %s  Variant: %s",
                        transcript.getTranscriptId(), variant));
            }
            else if ( isFivePrimeFlank(variant, transcript, flankSettings.fivePrimeFlankSize) ) {
                return createFlankFuncotation(variant, altAllele, transcript, reference, GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK);
            }
            else if ( isThreePrimeFlank(variant, transcript, flankSettings.threePrimeFlankSize) ) {
                return createFlankFuncotation(variant, altAllele, transcript, reference, GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK);
            }
            else {
                // This is an IGR, so we return nothing here. If we end up with only IGRs at the end of annotation,
                // we'll create the actual IGR funcotations at a higher level.
                return null;
            }
        }

        // Make sure that we have accounted for "special" allele cases:
        final GencodeFuncotation specialAlleleFuncotation = checkForAndCreateSpecialAlleleFuncotation(variant, altAllele, reference, transcript);
        if ( specialAlleleFuncotation != null ) {
            return specialAlleleFuncotation;
        }

        // Make sure the sub-regions in the transcript actually contain the variant:
        // TODO: this is slow, and repeats work that is done later in the process (we call getSortedCdsAndStartStopPositions when creating the sequence comparison)
        final int startPosInTranscript =  FuncotatorUtils.getStartPositionInTranscript(variant, getSortedCdsAndStartStopPositions(transcript), transcript.getGenomicStrand() );

        if ( GencodeGtfExonFeature.class.isAssignableFrom(containingSubfeature.getClass()) ) {

            if ( startPosInTranscript == -1 ) {
                // we overlap an exon but we don't start in one.  Right now this case cannot be handled.
                // Bubble up an exception and let the caller handle this case.
                // TODO: fix this case, issue #4307 (https://github.com/broadinstitute/gatk/issues/4307)
                throw new FuncotatorUtils.TranscriptCodingSequenceException("Cannot yet handle indels starting outside an exon and ending within an exon.");
            }
            else {
                // We have a coding region variant
                gencodeFuncotation = createExonFuncotation(variant, altAllele, reference, transcript, (GencodeGtfExonFeature) containingSubfeature);
            }
        }
        else if ( GencodeGtfUTRFeature.class.isAssignableFrom(containingSubfeature.getClass()) ) {
            // We have a UTR variant
            gencodeFuncotation = createUtrFuncotation(variant, altAllele, reference, transcript, (GencodeGtfUTRFeature) containingSubfeature);
        }
        else if ( GencodeGtfTranscriptFeature.class.isAssignableFrom(containingSubfeature.getClass()) ) {
            // We have an intron variant
            gencodeFuncotation = createIntronFuncotation(variant, altAllele, reference, transcript);
        }
        else {
            // Uh-oh!  Problemz.
            throw new GATKException.ShouldNeverReachHereException("Unable to determine type of variant-containing subfeature: " + containingSubfeature.getClass().getName());
        }

        return gencodeFuncotation;
    }

    /**
     * Checks the given variant alt allele for whether it's a special case.
     * If so, will create the appropriate funcotation for that case.
     * @param variant {@link VariantContext} for the given {@code altAllele}.
     * @param altAllele The alternate {@link Allele} to be checked for special cases.
     * @param reference {@link ReferenceContext} supporting the given {@code variant}.
     * @param transcript {@link GencodeGtfTranscriptFeature} containing the given {@code variant}.
     * @return A {@link GencodeFuncotation} for the given special case alt allele, or {@code null} if the allele is not a special case.
     */
    private GencodeFuncotation checkForAndCreateSpecialAlleleFuncotation(final VariantContext variant,
                                                                         final Allele altAllele,
                                                                         final ReferenceContext reference,
                                                                         final GencodeGtfTranscriptFeature transcript) {
        // If the alt allele is a spanning deletion or a symbolic allele, create an unknown funcotation:
        if (altAllele.isSymbolic() || altAllele.equals(Allele.SPAN_DEL) ) {
            return createFuncotationForSymbolicAltAllele(variant, altAllele, transcript, reference);
        }

        // If the reference allele or alternate allele contains any masked bases:
        if ( variant.getReference().getBaseString().contains(FuncotatorConstants.MASKED_ANY_BASE_STRING) || altAllele.getBaseString().contains(FuncotatorConstants.MASKED_ANY_BASE_STRING) ) {
            return createFuncotationForMaskedBases(variant, altAllele, transcript, reference);
        }

        return null;
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in a given {@code exon}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the current data set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param exon The {@link GencodeGtfExonFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code exon}.
     */
    private GencodeFuncotation createExonFuncotation(final VariantContext variant,
                                                     final Allele altAllele,
                                                     final ReferenceContext reference,
                                                     final GencodeGtfTranscriptFeature transcript,
                                                     final GencodeGtfExonFeature exon) {

        // Before we get started, check to see if this is a non-protein-coding feature.
        // If it is, we must handle it differently:
        if ( GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString().equals(transcript.getGeneType()) ) {
            return createCodingRegionFuncotationForProteinCodingFeature(variant, altAllele, reference, transcript, exon);
        }
        else {
            return createCodingRegionFuncotationForNonProteinCodingFeature(variant, altAllele, reference, transcript, exon);
        }
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in a coding region in a given {@code exon}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the current data set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param exon The {@link GencodeGtfExonFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code exon}.
     */
    private GencodeFuncotation createCodingRegionFuncotationForNonProteinCodingFeature(final VariantContext variant,
                                                                                       final Allele altAllele,
                                                                                       final ReferenceContext reference,
                                                                                       final GencodeGtfTranscriptFeature transcript,
                                                                                       final GencodeGtfExonFeature exon) {

        // Get the list of exons by their locations so we can use them to determine our location in the transcript and get
        // the transcript code itself:
        final List<? extends Locatable> exonPositionList = getSortedCdsAndStartStopPositions(transcript);

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, transcript);

        // Set the exon number:
        gencodeFuncotationBuilder.setTranscriptExonNumber(exon.getExonNumber());

        // Set our version:
        gencodeFuncotationBuilder.setVersion(version);

        // Set our transcript positions:
        setTranscriptPosition(variant, altAllele, transcript, gencodeFuncotationBuilder);

        // Get the reference bases for our current variant:
        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, reference, exon.getGenomicStrand(), referenceWindow);

        // Set the reference context with the bases from the sequence comparison
        // NOTE: The reference context is ALWAYS from the + strand
        gencodeFuncotationBuilder.setReferenceContext(referenceBases.getBaseString(Strand.POSITIVE));

        // Get our exon positions:
        final SimpleInterval exonOverlapInterval = FuncotatorUtils.getOverlappingExonPositions(variant.getReference(), altAllele, variant.getContig(), variant.getStart(), variant.getEnd(), transcript.getGenomicStrand(), exonPositionList);

        // Set the GC content:
        // Set the cDNA change:
        gencodeFuncotationBuilder
                .setGcContent(calculateGcContent(variant.getReference(), altAllele, reference, gcContentWindowSizeBases))
                .setcDnaChange(
                    FuncotatorUtils.getCodingSequenceChangeString(
                        FuncotatorUtils.getStartPositionInTranscript(variant, exonPositionList, transcript.getGenomicStrand()),
                        (exon.getGenomicStrand() == Strand.FORWARD ? variant.getReference().getBaseString() : ReadUtils.getBasesReverseComplement(variant.getReference().getBases())),
                        (exon.getGenomicStrand() == Strand.FORWARD ? altAllele.getBaseString() : ReadUtils.getBasesReverseComplement(altAllele.getBases())),
                        exon.getGenomicStrand(),
                        (exonOverlapInterval != null ? exonOverlapInterval.getStart() : null),
                        (exonOverlapInterval != null ? exonOverlapInterval.getStart() : null),
                        (exonOverlapInterval != null ? variant.getStart() : null)
                    )
                );

        // Set the VariantClassification through a simple equivalency on the gene type (since we have no transcript info):
        gencodeFuncotationBuilder.setVariantClassification( convertGeneTranscriptTypeToVariantClassification(exon.getGeneType()) );

        // Set our data source name:
        gencodeFuncotationBuilder.setDataSourceName(getName());

        //==============================================================================================================

        return gencodeFuncotationBuilder.build();
    }

    /**
     * Set the transcript start and end position in the given {@link GencodeFuncotationBuilder}.
     * @param variant The {@link VariantContext} for which to set the position.
     * @param altAllele The alternate {@link Allele} for which to set the position.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the {@code variant} occurs.
     * @param gencodeFuncotationBuilder The {@link GencodeFuncotationBuilder}  in which to set the transcript start and end positions.
     */
    private void setTranscriptPosition(final VariantContext variant, final Allele altAllele, final GencodeGtfTranscriptFeature transcript, final GencodeFuncotationBuilder gencodeFuncotationBuilder) {
        final int txStart = FuncotatorUtils.getTranscriptAlleleStartPosition(variant, transcript.getExons(), transcript.getGenomicStrand());
        setTranscriptPosition(variant, altAllele, txStart, gencodeFuncotationBuilder);
    }

    /**
     * Set the transcript start and end position in the given {@link GencodeFuncotationBuilder}.
     * @param variant The {@link VariantContext} for which to set the position.
     * @param altAllele The alternate {@link Allele} for which to set the position.
     * @param txStart The start position of the {@code variant} in the transcript.
     * @param gencodeFuncotationBuilder The {@link GencodeFuncotationBuilder}  in which to set the transcript start and end positions.
     */
    private void setTranscriptPosition(final VariantContext variant, final Allele altAllele, final int txStart, final GencodeFuncotationBuilder gencodeFuncotationBuilder) {

        // Set the transcript start position:
        gencodeFuncotationBuilder.setTranscriptStartPos(txStart);

        // For SNPs, start == end:
        if ( variant.getReference().length() == 1 ) {
            gencodeFuncotationBuilder.setTranscriptEndPos( txStart );
        }
        else {
            // For indels, we have to adjust by the leading placeholder base.
            // For MNPs we can use the raw length.
            int indelOffset = 0;
            if ( GATKVariantContextUtils.isIndel(variant.getReference(), altAllele) ) {
                // Account for the leading base in the variant:
                indelOffset = 1;
            }
            // Subtract 1 here for the inclusive nature of the positions.
            // NOTE: this is in addition to the indel offset.
            gencodeFuncotationBuilder.setTranscriptEndPos(txStart + variant.getReference().length() - indelOffset -1);
        }
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in a coding region in a given {@code exon}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the current data set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param exon The {@link GencodeGtfExonFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code exon}.
     */
    private GencodeFuncotation createCodingRegionFuncotationForProteinCodingFeature(final VariantContext variant,
                                                                                    final Allele altAllele,
                                                                                    final ReferenceContext reference,
                                                                                    final GencodeGtfTranscriptFeature transcript,
                                                                                    final GencodeGtfExonFeature exon) {

        // Get the list of exons by their locations so we can use them to determine our location in the transcript and get
        // the transcript code itself:
        final List<? extends Locatable> exonPositionList = getSortedCdsAndStartStopPositions(transcript);

        // NOTE: Regardless of strandedness, we always report the alleles as if they appeared in the forward direction.
        final GencodeFuncotation.VariantType variantType =
                getVariantType(variant.getReference(),
                        altAllele);

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, transcript);

        // Set the exon number:
        gencodeFuncotationBuilder.setTranscriptExonNumber(exon.getExonNumber());

        // Set our version:
        gencodeFuncotationBuilder.setVersion(version);

        // Set up our SequenceComparison object so we can calculate some useful fields more easily
        // These fields can all be set without knowing the alternate allele:
        final SequenceComparison sequenceComparison = createSequenceComparison(variant, altAllele, reference, transcript, exonPositionList, transcriptIdMap, transcriptFastaReferenceDataSource, true);

        // Set our transcript positions:
        setTranscriptPosition(variant, altAllele, sequenceComparison.getTranscriptAlleleStart(), gencodeFuncotationBuilder);

        // Set the reference context with the bases from the sequence comparison
        // NOTE: The reference context is ALWAYS from the + strand, so we need to reverse our bases back in the - case:
        if ( sequenceComparison.getStrand() == Strand.POSITIVE ) {
            gencodeFuncotationBuilder.setReferenceContext(sequenceComparison.getReferenceBases());
        }
        else {
            gencodeFuncotationBuilder.setReferenceContext(ReadUtils.getBasesReverseComplement(sequenceComparison.getReferenceBases().getBytes()));
        }

        // Set the GC content
        // Set the cDNA change:
        gencodeFuncotationBuilder.setGcContent(sequenceComparison.getGcContent())
                                 .setcDnaChange(FuncotatorUtils.getCodingSequenceChangeString(
                                         sequenceComparison.getCodingSequenceAlleleStart(),
                                         sequenceComparison.getReferenceAllele(),
                                         sequenceComparison.getAlternateAllele(),
                                         sequenceComparison.getStrand(),
                                         sequenceComparison.getExonStartPosition(),
                                         sequenceComparison.getExonEndPosition(),
                                         sequenceComparison.getAlleleStart()
                                 ));

        //==============================================================================================================
        // Set the Codon and Protein changes and the Variant Classification
        // but only if we have the sequence information to do so.
        // NOTE: This should always be true in this method, but we need to have this if statement just in case it does.
        //       A warning will have been generated in createSequenceComparison if the sequenceComparison does not have
        //       coding sequence information.
        if ( sequenceComparison.hasSequenceInfo() ) {
            final String codonChange = FuncotatorUtils.getCodonChangeString(sequenceComparison, exon.getStartCodon());
            final String proteinChange = FuncotatorUtils.renderProteinChangeString(sequenceComparison, exon.getStartCodon());

            gencodeFuncotationBuilder.setCodonChange(codonChange)
                                     .setProteinChange(proteinChange);

            // Set the Variant Classification:
            final GencodeFuncotation.VariantClassification varClass = createVariantClassification(
                    variant, altAllele, variantType, exon, transcript.getExons().size(), sequenceComparison, spliceSiteVariantWindowBases
            );
            final GencodeFuncotation.VariantClassification secondaryVarClass;
            gencodeFuncotationBuilder.setVariantClassification(varClass);
            if ( varClass == GencodeFuncotation.VariantClassification.SPLICE_SITE ) {
                secondaryVarClass = getVariantClassificationForCodingRegion(variant, altAllele, variantType, sequenceComparison);
                gencodeFuncotationBuilder.setSecondaryVariantClassification(secondaryVarClass);
            }
        }
        else {
            // Set the variant classification here.
            // We should have sequence information but we don't... this is not good, but we have to put something here:
            gencodeFuncotationBuilder.setVariantClassification( convertGeneTranscriptTypeToVariantClassification(exon.getGeneType()) );
        }

        // Set our data source name:
        gencodeFuncotationBuilder.setDataSourceName(getName());

        return gencodeFuncotationBuilder.build();
    }

    /**
     * Gets a list of locatables representing the start codon, cds, and stop codon containing coding regions within the given {@code transcript}.
     * These locatables are sorted by exon-number order.
     * @param transcript A {@link GencodeGtfTranscriptFeature} from which to pull the exons.
     * @return A list of {@link Locatable} objects representing the regions in the exons in the given {@code transcript} in the order in which the appear in the expressed protein.
     */
    @VisibleForTesting
    static List<? extends Locatable> getSortedCdsAndStartStopPositions(final GencodeGtfTranscriptFeature transcript) {

        // Sort by exon number first:
        transcript.getExons().sort((lhs, rhs) -> lhs.getExonNumber() < rhs.getExonNumber() ? -1 : (lhs.getExonNumber() > rhs.getExonNumber() ) ? 1 : 0 );

        final List<GencodeGtfFeature> regionList = new ArrayList<>(transcript.getExons().size());
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {

            // Add in a CDS region:
            if ( exon.getCds() != null ) {

                // If we have a start codon that is not in the CDS for some reason,
                // we need to add it to our list:
                if (exon.getStartCodon() != null) {
                    if ( !exon.getCds().contains(exon.getStartCodon()) ) {
                        regionList.add( exon.getStartCodon() );
                    }
                }

                regionList.add( exon.getCds() );

                // If we have a stop codon that is not in the CDS for some reason,
                // we need to add it to our list:
                if (exon.getStopCodon() != null) {
                    if ( !exon.getCds().contains(exon.getStopCodon()) ) {
                        regionList.add( exon.getStopCodon() );
                    }
                }

            }
            else if (exon.getStartCodon() != null) {
                regionList.add( exon.getStartCodon() );
            }
            else if ( exon.getStopCodon() != null ) {
                regionList.add( exon.getStopCodon() );
            }
        }

        // Now we have to sort by the transcript direction because of assumptions in the code later.
        // If the transcript is in the + direction, we don't do anything.
        // For - transcripts, we sort in order by reverse starting position.
        // This is a fix for issue https://github.com/broadinstitute/gatk/issues/7051
        if (transcript.getGenomicStrand() == Strand.NEGATIVE) {
            regionList.sort(Comparator.comparingInt(GencodeGtfFeature::getStart).reversed());
        }

        return regionList;
    }

    /**
     * Gets the {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} of the given {@code altAllele} for the given {@code variant}.
     * @param variant The {@link VariantContext} to classify.
     * @param altAllele The {@link Allele} of the given {@code variant} to classify.
     * @param variantType The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType} of the given {@code variant}.
     * @param exon The {@link GencodeGtfExonFeature} in which the given {@code variant} occurs.
     * @param numberOfExonsInTranscript The number of exons in the transcript in which the given {@code variant} occurs. (Must be > 0).
     * @param sequenceComparison The {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} for the given {@code variant}.
     * @param spliceSiteWindowBases The window on either side of splice sites to mark variants as {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#SPLICE_SITE}.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} based on the given {@code allele}, {@code variant}, {@code exon}, and {@code sequenceComparison}.
     */
    @VisibleForTesting
    static GencodeFuncotation.VariantClassification createVariantClassification(final VariantContext variant,
                                                                                final Allele altAllele,
                                                                                final GencodeFuncotation.VariantType variantType,
                                                                                final GencodeGtfExonFeature exon,
                                                                                final int numberOfExonsInTranscript,
                                                                                final SequenceComparison sequenceComparison,
                                                                                final int spliceSiteWindowBases){
        Utils.nonNull(variant);
        Utils.nonNull(altAllele);
        Utils.nonNull(variantType);
        Utils.nonNull(exon);
        ParamUtils.isPositive(numberOfExonsInTranscript, "Number of exons in transcript must be positive (given: " +numberOfExonsInTranscript + ")");
        Utils.nonNull(sequenceComparison);

        GencodeFuncotation.VariantClassification varClass = null;

        boolean hasBeenClassified = false;

        // Adjust the variant interval for the overlap check, specifically to properly test for the indel cases.
        // That is, since we always have a preceding base with indels, we adjust the variant interval to remove this
        // base so we can use the `overlaps` methods and get the correct results:
        final SimpleInterval realVariationInterval = getBasesChangedIntervalIgnoringLeadingVcfContextBase(variant, altAllele);

        // Check for non-stop first:
        if ( (exon.getStopCodon() != null) && (exon.getStopCodon().overlaps(realVariationInterval)) ) {

            boolean foundStop = false;

            // The -1 here is to account for the exclusive second argument to `String::substring`.
            // This will allow for all potential codons to be checked here.
            for (int i = 0; (i+AminoAcid.CODON_LENGTH-1) < sequenceComparison.getAlignedCodingSequenceAlternateAllele().length(); i+=AminoAcid.CODON_LENGTH ){
                final String codon = sequenceComparison.getAlignedCodingSequenceAlternateAllele().substring(i, i+AminoAcid.CODON_LENGTH);
                if (FuncotatorUtils.getEukaryoticAminoAcidByCodon(codon) == AminoAcid.STOP_CODON) {
                    foundStop = true;
                    break;
                }
            }

            if ( !foundStop ) {
                varClass = GencodeFuncotation.VariantClassification.NONSTOP;
                hasBeenClassified = true;
            }
        }

        // Now check all other cases:
        if ( !hasBeenClassified ) {

            // First check for splice site:

            final boolean isInternalExon = (exon.getExonNumber() != 1) && (exon.getExonNumber() != numberOfExonsInTranscript);
            final boolean doLeftOverlapCheck   = isInternalExon ||
                    ((sequenceComparison.getStrand() == Strand.NEGATIVE) && (exon.getExonNumber() == 1)) ||
                    ((sequenceComparison.getStrand() == Strand.POSITIVE) && (exon.getExonNumber() == numberOfExonsInTranscript));
            final boolean doRightOverlapCheck  = isInternalExon ||
                    ((sequenceComparison.getStrand() == Strand.POSITIVE) && (exon.getExonNumber() == 1)) ||
                    ((sequenceComparison.getStrand() == Strand.NEGATIVE) && (exon.getExonNumber() == numberOfExonsInTranscript));

            // Flags to check if the variant overlaps the ends of this exon:
            boolean overlapsLeft  = false;
            boolean overlapsRight = false;

            // Adjust the exon interval if we have an insertion because everything needs to be adjusted to account
            // for the newly inserted bases:
            final int adjustedExonStart = adjustLocusForInsertion(exon.getStart(), variant, altAllele, realVariationInterval);
            final int adjustedExonEnd = adjustLocusForInsertion(exon.getEnd(), variant, altAllele, realVariationInterval);

            // If we have 0 padding, we want to make sure the variant overlaps the exon start itself,
            // not within a base of the exon start, so we have to create this adjustment here so that we can
            // not subtract one and create an invalid interval.
            final int intervalEndCoordAdjuster = spliceSiteWindowBases > 0 ? 1 : 0;

            if ( doLeftOverlapCheck ) {
                final SimpleInterval leftSideInterval = new SimpleInterval(exon.getContig(), adjustedExonStart - spliceSiteWindowBases, adjustedExonStart + (spliceSiteWindowBases-intervalEndCoordAdjuster));
                overlapsLeft = leftSideInterval.overlaps(realVariationInterval);
            }
            if ( doRightOverlapCheck ) {
                final SimpleInterval rightSideInterval = new SimpleInterval(exon.getContig(), adjustedExonEnd - spliceSiteWindowBases + 1, adjustedExonEnd + (spliceSiteWindowBases-intervalEndCoordAdjuster) + 1);
                overlapsRight = rightSideInterval.overlaps(realVariationInterval);
            }

            // Check for splice site variants.
            if ( overlapsLeft || overlapsRight ) {
                varClass = GencodeFuncotation.VariantClassification.SPLICE_SITE;
            }
            else if ((exon.getStartCodon() != null) && (exon.getStartCodon().overlaps(realVariationInterval))) {
                switch (variantType) {
                    case INS:
                        varClass = GencodeFuncotation.VariantClassification.START_CODON_INS;
                        break;
                    case DEL:
                        varClass = GencodeFuncotation.VariantClassification.START_CODON_DEL;
                        break;
                    default:
                        varClass = GencodeFuncotation.VariantClassification.START_CODON_SNP;
                        break;
                }
            }
            else {
                varClass = getVariantClassificationForCodingRegion(variant, altAllele, variantType, sequenceComparison);
            }
        }

        return varClass;
    }

    /**
     * Shifts a given genomic locus by the number of inserted bases in the given {@code variant}/{@code altAllele} pair
     * if the position occurs after the start of the insertion (and if the variant is an insertion).
     * Assumes that the position and the variant are on the same contig.
     * @param genomicLocus The locus to be potentially adjusted.
     * @param variant The {@link VariantContext} to use for the reference allele.
     * @param altAllele The alternate {@link Allele} to use to compare against the reference allele.
     * @param changedBasesInterval The interval representing the bases actually changed in the given {@code variant} and {@code altAllele} pair.
     * @return A genomic locus adjusted for the bases inserted before it, if any.
     */
    private static int adjustLocusForInsertion(final int genomicLocus, final VariantContext variant,
                                               final Allele altAllele, final SimpleInterval changedBasesInterval) {

        int adjustedPosition = genomicLocus;

        // Is our variant / alt pair an insertion AND Is our position after the changed bases have started:
        if ( (altAllele.length() > variant.getReference().length()) && (genomicLocus > changedBasesInterval.getStart()) ) {
            // Adjust our position by the number of inserted bases:
            adjustedPosition += (changedBasesInterval.getEnd() - changedBasesInterval.getStart() + 1);
        }

        return adjustedPosition;

    }

    /**
     * Get the {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} for a given {@code variant}/{@code allele} in a coding region.
     * @param variant The {@link VariantContext} to classify.
     * @param altAllele The {@link Allele} of the given {@code variant} to classify.
     * @param variantType The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType} of the given {@code variant}.
     * @param sequenceComparison The {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} for the given {@code variant}.  Must have a non-null proteinChangeInfo Object.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} based on the given {@code allele}, {@code variant}, {@code exon}, and {@code sequenceComparison}.
     */
    private static GencodeFuncotation.VariantClassification getVariantClassificationForCodingRegion(final VariantContext variant,
                                                                                             final Allele altAllele,
                                                                                             final GencodeFuncotation.VariantType variantType,
                                                                                             final SequenceComparison sequenceComparison) {
        final GencodeFuncotation.VariantClassification varClass;

        if (variantType == GencodeFuncotation.VariantType.INS) {
            if ( GATKVariantContextUtils.isFrameshift(variant.getReference(), altAllele)) {
                varClass = GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS;
            }
            else {
                varClass = GencodeFuncotation.VariantClassification.IN_FRAME_INS;
            }
        }
        else if (variantType == GencodeFuncotation.VariantType.DEL) {
            if (GATKVariantContextUtils.isFrameshift(variant.getReference(), altAllele)) {
                varClass = GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL;
            }
            else {
                varClass = GencodeFuncotation.VariantClassification.IN_FRAME_DEL;
            }
        }
        else {
            // This is a SNP/MNP
            // We just check to see what the protein change is to check for MISSENSE/NONSENSE/SILENT:
            varClass = getVarClassFromEqualLengthCodingRegions( sequenceComparison );
        }

        return varClass;
    }

    /**
     * Gets the {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} for a {@code variant} where the reference and alternate
     * alleles are the same length and the variant appears in a coding region.
     * This essentially compares the amino acid sequences of both alleles and returns a value based on the differences between them.
     * @return The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} corresponding to the given variant / reference allele / alternate allele.
     */
    private static GencodeFuncotation.VariantClassification getVarClassFromEqualLengthCodingRegions(final SequenceComparison sequenceComparison) {

        GencodeFuncotation.VariantClassification varClass = GencodeFuncotation.VariantClassification.SILENT;

        boolean foundErroneousStop = false;

        final String refAaSeq = sequenceComparison.getProteinChangeInfo().getRefAaSeq();
        final String altAaSeq = sequenceComparison.getProteinChangeInfo().getAltAaSeq();

        for ( int i = 0; i < altAaSeq.length(); ++i ) {
            final char altAminoAcid = altAaSeq.charAt(i);

            if ( altAminoAcid != refAaSeq.charAt(i) ) {
                if ( FuncotatorUtils.getAminoAcidByLetter(altAminoAcid) == AminoAcid.STOP_CODON ) {
                    foundErroneousStop = true;
                    break;
                }
                varClass = GencodeFuncotation.VariantClassification.MISSENSE;
            }
        }

        if ( foundErroneousStop ) {
            varClass = GencodeFuncotation.VariantClassification.NONSENSE;
        }

        return varClass;
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in an untranslated region in a given {@code transcript}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the current data set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param utr The {@link GencodeGtfUTRFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code utr}.
     */
    private GencodeFuncotation createUtrFuncotation(final VariantContext variant,
                                                    final Allele altAllele,
                                                    final ReferenceContext reference,
                                                    final GencodeGtfTranscriptFeature transcript,
                                                    final GencodeGtfUTRFeature utr) {

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, transcript);

        // Find which exon this UTR is in:
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
            if ( exon.contains( utr ) ) {
                gencodeFuncotationBuilder.setTranscriptExonNumber( exon.getExonNumber() );
                break;
            }
        }

        // Set GC Content:
        gencodeFuncotationBuilder.setGcContent( calculateGcContent( variant.getReference(), altAllele, reference, gcContentWindowSizeBases ) );

        // Get the strand:
        final Strand strand = transcript.getGenomicStrand();

        // Get the strand-corrected alleles from the inputs.
        // Also get the reference sequence for the variant region.
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele                        strandCorrectedAltAllele = FuncotatorUtils.getStrandCorrectedAllele(altAllele, strand);
        final StrandCorrectedReferenceBases referenceBases           = FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, reference, strand, referenceWindow);

        // Set our reference sequence in the Gencode Funcotation Builder:
        // NOTE: The reference context is ALWAYS from the + strand
        gencodeFuncotationBuilder.setReferenceContext(referenceBases.getBaseString(Strand.POSITIVE));

        // Set whether it's the 5' or 3' UTR:
        if ( is5PrimeUtr(utr, transcript) ) {

            // We're 5' to the coding region.
            // Default our variant classification to 5'UTR and try to refine it later:
            gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR);

            // Now we can check for de novo starts:

            // Only try to get the sequence if our transcript occurs in the FASTA file:
            if ( transcriptIdMap.containsKey(transcript.getTranscriptId()) ) {

                // Get the 5' UTR sequence here.
                // Note: We grab 3 extra bases at the end (from the coding sequence) so that we can check for denovo starts
                //       even if the variant occurs in the last base of the UTR.
                final int numExtraTrailingBases = variant.getReference().length() < defaultNumTrailingBasesForUtrAnnotationSequenceConstruction ? defaultNumTrailingBasesForUtrAnnotationSequenceConstruction : variant.getReference().length() + 1;
		int transcriptLength = getTranscriptLength( transcriptFastaReferenceDataSource, transcriptIdMap, transcript );
		
                final String fivePrimeUtrCodingSequence =
                        getFivePrimeUtrSequenceFromTranscriptFasta( transcript.getTranscriptId(), transcriptIdMap, transcriptFastaReferenceDataSource, numExtraTrailingBases, transcriptLength );

                // Get our start position in our coding sequence:
                final int codingStartPos = FuncotatorUtils.getStartPositionInTranscript(variant, transcript.getExons(), strand);

                // We can use the `referenceBases` to get the alt subsequence, but we need to adjust it for indels.
                // We make this adjustment with knowledge of how the referenceBases sequence was created in the first
                // place.
                final int indelOffset = variant.getReference().length() != altAllele.length() ? 1 : 0;

                final int frontOffset = strand == Strand.POSITIVE ? indelOffset : 0;
                final int backOffset  = strand == Strand.NEGATIVE ? indelOffset : 0;

                final String rawAltUtrSubSequence = (referenceBases.getBaseString().substring(referenceWindow-numLeadingBasesForUtrAnnotationSequenceConstruction + frontOffset, referenceWindow) +
                        strandCorrectedAltAllele +
                        referenceBases.getBaseString().substring(referenceWindow + variant.getReference().length(), referenceWindow + numExtraTrailingBases + backOffset));

                // Check for de novo starts in the raw sequence:
                boolean startFound = false;
                int codingRegionOffset = frontOffset-numLeadingBasesForUtrAnnotationSequenceConstruction;
                for ( int i = 0; (i+AminoAcid.CODON_LENGTH < rawAltUtrSubSequence.length()) ; ++i ) {
                    startFound = FuncotatorUtils.getEukaryoticAminoAcidByCodon( rawAltUtrSubSequence.substring(i, i+AminoAcid.CODON_LENGTH) ) == AminoAcid.METHIONINE;
                    if (startFound) {
                        codingRegionOffset += i;
                        break;
                    }
                }
                // If we found a start codon, we should set the variant classification as such:
                if ( startFound ) {
                    if ( FuncotatorUtils.isInFrameWithEndOfRegion(codingStartPos + codingRegionOffset, fivePrimeUtrCodingSequence.length()) ) {
                        gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME);
                    }
                    else {
                        gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME);
                    }
                }
            }
            else {
                logger.warn("Attempted to process transcript information for transcript WITHOUT sequence data.  Ignoring sequence information for Gencode Transcript ID: " + transcript.getTranscriptId());
            }
        }
        else {
            gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR);
        }

        // Set our version:
        gencodeFuncotationBuilder.setVersion(version);

        // Set our data source name:
        gencodeFuncotationBuilder.setDataSourceName(getName());

        return gencodeFuncotationBuilder.build();
    }

    /**
     * Get the length of the given transcript by querying the reference for the transcript ID
     * @return the length of the requested {@code transcript} in the {@code transcriptFastaReferenceDataSource}
     */
    private static int getTranscriptLength(
	    ReferenceDataSource transcriptFastaReferenceDataSource,
	    Map<String, MappedTranscriptIdInfo> transcriptIdMap,
	    GencodeGtfTranscriptFeature transcript
    ){
		SAMSequenceDictionary dict = transcriptFastaReferenceDataSource.getSequenceDictionary();
		final MappedTranscriptIdInfo transcriptMapIdAndMetadata = transcriptIdMap.get( transcript.getTranscriptId() );

		if ( transcriptMapIdAndMetadata == null ) {
			String transcriptId = transcript.getTranscriptId();
			throw new UserException.BadInput( "Unable to find the given Transcript ID in our transcript list for ou     r 5'UTR (not in given transcript FASTA file): " + transcriptId );
		}

		String transcriptReferenceKey = transcriptMapIdAndMetadata.mapKey;
		SAMSequenceRecord transcriptRecord = dict.getSequence( transcriptReferenceKey );
		int transcriptLength = transcriptRecord.getSequenceLength();

		return transcriptLength;
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in an intron in the given {@code transcript}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the given {@code variant}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code transcript}.
     */
    private GencodeFuncotation createIntronFuncotation(final VariantContext variant,
                                                       final Allele altAllele,
                                                       final ReferenceContext reference,
                                                       final GencodeGtfTranscriptFeature transcript) {

        // Get the strand-corrected alleles from the inputs.
        // Also get the reference sequence for the variant region.
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele strandCorrectedRefAllele = FuncotatorUtils.getStrandCorrectedAllele(variant.getReference(), transcript.getGenomicStrand());
        final Allele strandCorrectedAltAllele = FuncotatorUtils.getStrandCorrectedAllele(altAllele, transcript.getGenomicStrand());

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, transcript);

        // Set our reference sequence in the Gencode Funcotation Builder:

        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, reference, transcript.getGenomicStrand(), referenceWindow);

        // NOTE: The reference context is ALWAYS from the + strand:
        gencodeFuncotationBuilder.setReferenceContext(referenceBases.getBaseString(Strand.POSITIVE));

        // Set the VariantClassification:
        if ( GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString().equals(transcript.getGeneType()) ) {
            gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.INTRON);
        }
        else {
            gencodeFuncotationBuilder.setVariantClassification(convertGeneTranscriptTypeToVariantClassification(transcript.getGeneType()));
        }

        // Set GC Content:
        gencodeFuncotationBuilder.setGcContent( calculateGcContent( variant.getReference(), altAllele, reference, gcContentWindowSizeBases ) );

        // Need to check if we're within the window for splice site variants:
        final GencodeGtfExonFeature spliceSiteExon = getExonWithinSpliceSiteWindow(variant, altAllele, transcript, spliceSiteVariantWindowBases);
        if ( spliceSiteExon != null ) {
            // Set the variant classification:
            gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.SPLICE_SITE)
                                     .setSecondaryVariantClassification(GencodeFuncotation.VariantClassification.INTRON);

            // In deletions we have added a base to the front because of VCF requirements, thus we add an
            // offset of 1 to account for that:
            // (TODO: come to think of it this is really bad, because we're tying our parsing / computations to a data format).
            int offsetIndelAdjustment = 0;
            if ( GATKVariantContextUtils.isDeletion(strandCorrectedRefAllele, strandCorrectedAltAllele) ) {
                offsetIndelAdjustment = 1;
            }

            gencodeFuncotationBuilder.setCodonChange(
                    FuncotatorUtils.createSpliceSiteCodonChange(variant.getStart(), spliceSiteExon.getExonNumber(), spliceSiteExon.getStart(), spliceSiteExon.getEnd(), transcript.getGenomicStrand(), offsetIndelAdjustment)
            );
        }

        // Set our cDNA string:
        gencodeFuncotationBuilder.setcDnaChange(
                FuncotatorUtils.createIntronicCDnaString(
                        variant.getStart(),
                        transcript.getExons(),
                        strandCorrectedRefAllele.getBaseString(),
                        strandCorrectedAltAllele.getBaseString()
                )
        );

        // Set our version:
        gencodeFuncotationBuilder.setVersion(version);

        // Set our data source name:
        gencodeFuncotationBuilder.setDataSourceName(getName());

        return gencodeFuncotationBuilder.build();
    }

    private static SimpleInterval getBasesChangedIntervalIgnoringLeadingVcfContextBase(final VariantContext variant,
                                                                                       final Allele altAllele) {

        // Adjust the variant interval for the overlap check, specifically to properly test for the indel cases:
        final SimpleInterval changedBasesInterval;
        if ( GATKVariantContextUtils.typeOfVariant(variant.getReference(), altAllele).equals(VariantContext.Type.INDEL) ) {

            final int adjustedStart;
            final int end;
            // Insertion:
            if ( variant.getReference().length() < altAllele.length() ) {
                // We use the variant start here because by inserting bases we're not making the actually changed
                // bases any closer to the boundaries of the exon.
                // That is, the inserted bases shouldn't count towards the extents of the variant in genomic space
                // with respect to the exon boundaries.
                adjustedStart = FuncotatorUtils.getIndelAdjustedAlleleChangeStartPosition(variant, altAllele);
                end = variant.getEnd() + (altAllele.length() - (adjustedStart - variant.getStart()));
            }
            // Deletion:
            else {
                // Because we're deleting bases from the reference allele, we need to adjust the start position of the
                // variant to reflect that the leading base(s) do(es) not matter.
                adjustedStart = FuncotatorUtils.getIndelAdjustedAlleleChangeStartPosition(variant, altAllele);

                // The original end position should be correct:
                end = variant.getEnd();
            }

            changedBasesInterval = new SimpleInterval( variant.getContig(), adjustedStart, end );
        }
        else {
            changedBasesInterval = new SimpleInterval(variant);
        }

        return changedBasesInterval;
    }

    /**
     * Gets the {@link GencodeGtfExonFeature} that is within {@code spliceSiteVariantWindowBases} bases of the given {@code variant}.
     * @param variant The {@link VariantContext} to check for position within {@code spliceSiteVariantWindowBases} bases of each {@link GencodeGtfExonFeature} in {@code transcript}.
     * @param altAllele The alternate {@link Allele} to check against the reference allele in the given {@code variant}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} containing the given {@code variant}.
     * @return The {@link GencodeGtfExonFeature} that is within {@code spliceSiteVariantWindowBases} bases of the given {@code variant}; {@code null} if no such {@link GencodeGtfExonFeature} exists in the given {@code transcript}.
     */
    private static GencodeGtfExonFeature getExonWithinSpliceSiteWindow( final VariantContext variant,
                                                                        final Allele altAllele,
                                                                        final GencodeGtfTranscriptFeature transcript,
                                                                        final int spliceSiteVariantWindowBases ) {
        GencodeGtfExonFeature spliceSiteExon = null;

        // Adjust the variant interval for the overlap check, specifically to properly test for the indel cases:
        final SimpleInterval changedBasesInterval = getBasesChangedIntervalIgnoringLeadingVcfContextBase(variant, altAllele);

        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {

            // Adjust the exon interval if we have an insertion because everything needs to be adjusted to account
            // for the newly inserted bases:
            final int adjustedExonStart = adjustLocusForInsertion(exon.getStart(), variant, altAllele, changedBasesInterval);
            final int adjustedExonEnd = adjustLocusForInsertion(exon.getEnd(), variant, altAllele, changedBasesInterval);

            final SimpleInterval exonStartInterval = new SimpleInterval(exon.getContig(), adjustedExonStart, adjustedExonStart);
            final SimpleInterval exonEndInterval   = new SimpleInterval(exon.getContig(), adjustedExonEnd, adjustedExonEnd);

            if ( changedBasesInterval.overlapsWithMargin(exonStartInterval, spliceSiteVariantWindowBases) ||
                 changedBasesInterval.overlapsWithMargin(exonEndInterval, spliceSiteVariantWindowBases) ) {
                spliceSiteExon = exon;
                break;
            }
        }

        return spliceSiteExon;
    }

    /**
     * Get the subfeature contained in {@code transcript} that contains the given {@code Locatable}.
     * The returned subfeature will be of type {@link GencodeGtfFeature} with concrete type based on the type of region
     * in which the variant is found:
     *      Found in coding region -> {@link GencodeGtfExonFeature}
     *      Found in UTR ->{@link GencodeGtfUTRFeature}
     *      Found in intron ->{@link GencodeGtfTranscriptFeature}
     *      Not Found in transcript ->{@code null}
     * @param variant A {@link Locatable} of which to determine the containing subfeature.
     * @param transcript A {@link GencodeGtfTranscriptFeature} in which to find the subfeature containing the given {@code variant}.
     * @return The {@link GencodeGtfFeature} corresponding to the subfeature of {@code transcript} in which the given {@code variant} was found.
     */
    private static GencodeGtfFeature getContainingGtfSubfeature(final Locatable variant, final GencodeGtfTranscriptFeature transcript) {

        boolean determinedRegionAlready = false;
        GencodeGtfFeature subFeature = null;

        if ( transcript.contains(variant) ) {

            if ( transcript.getUtrs().size() > 0 ) {
                for ( final GencodeGtfUTRFeature utr : transcript.getUtrs() ) {
                    if ( utr.overlaps(variant) ) {
                        subFeature = utr;
                        determinedRegionAlready = true;
                    }
                }
            }

            // Even though we may have an overlapping UTR already, we may be able to find a spot in the transcript
            // where this overlaps something more meaningful.
            // For example, see HG19 - chr19:8959608
            for (final GencodeGtfExonFeature exon : transcript.getExons()) {
                // TODO: This `contains` is here for issue #4307 - https://github.com/broadinstitute/gatk/issues/4307
                if ((exon.getCds() != null) && (exon.getCds().contains(variant))) {
                    subFeature = exon;
                    determinedRegionAlready = true;
                }
                else if ((exon.getStartCodon() != null) && (exon.getStartCodon().overlaps(variant))) {
                    subFeature = exon;
                    determinedRegionAlready = true;
                }
                else if ((exon.getStopCodon() != null) && (exon.getStopCodon().overlaps(variant))) {
                    subFeature = exon;
                    determinedRegionAlready = true;
                }
            }

            if ( !determinedRegionAlready ) {
                subFeature = transcript;
            }
        }

        return subFeature;
    }

    /**
     * Creates a {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} object with the fields populated.
     * @param variant The {@link VariantContext} for the current variant.
     * @param alternateAllele The current alternate {@link Allele} for the variant.
     * @param reference The {@link ReferenceContext} for the current sample set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} for the current gene feature / alt allele.
     * @param exonPositionList A {@link List} of {@link htsjdk.samtools.util.Locatable} objects representing exon positions in the transcript.
     * @param transcriptIdMap The {@link Map} of TranscriptID to {@link MappedTranscriptIdInfo} for all transcripts in the current Gencode data source.
     * @param transcriptFastaReferenceDataSource The {@link ReferenceDataSource} of the transcript FASTA file containing the sequence information for all Transcripts in the current Gencode data source.
     * @param processSequenceInformation If {@code true} will attempt to process and create sequence information for the given {@code variant}.
     * @return A populated {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} object.
     */
    @VisibleForTesting
    static SequenceComparison createSequenceComparison(final VariantContext variant,
                                                       final Allele alternateAllele,
                                                       final ReferenceContext reference,
                                                       final GencodeGtfTranscriptFeature transcript,
                                                       final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                                       final Map<String, MappedTranscriptIdInfo> transcriptIdMap,
                                                       final ReferenceDataSource transcriptFastaReferenceDataSource,
                                                       final boolean processSequenceInformation) {

        // TODO: Somewhere down the line we should adjust the positions at creation-time to account for the leading bases in VCF input files.  (issue 5349 - https://github.com/broadinstitute/gatk/issues/5349)
        // This will have ramifications down the line for all fields that get rendered.

        final SequenceComparison sequenceComparison = new SequenceComparison();

        // Get the contig:
        sequenceComparison.setContig(variant.getContig());

        // Get the strand:
        sequenceComparison.setStrand( transcript.getGenomicStrand() );

        // Get the alleles from the inputs
        // Also get the reference sequence for the variant region
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        // Get the strand-corrected alleles from the inputs.
        // Also get the reference sequence for the variant region.
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele                        refAllele      = FuncotatorUtils.getStrandCorrectedAllele(variant.getReference(), transcript.getGenomicStrand());
        final Allele                        altAllele      = FuncotatorUtils.getStrandCorrectedAllele(alternateAllele, transcript.getGenomicStrand());
        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.createReferenceSnippet(variant.getReference(), alternateAllele, reference, transcript.getGenomicStrand(), referenceWindow);

        // Set our reference sequence in the SequenceComparison:
        sequenceComparison.setReferenceWindow(referenceWindow);
        sequenceComparison.setReferenceBases(referenceBases.getBaseString());

        // Set our GC content:
        sequenceComparison.setGcContent(calculateGcContent(variant.getReference(), altAllele, reference, gcContentWindowSizeBases));

        // Get the ref allele:
        sequenceComparison.setReferenceAllele(refAllele.getBaseString());

        // Get the allele genomic start position:
        sequenceComparison.setAlleleStart(variant.getStart());

        // Get the allele transcript start position:
        sequenceComparison.setTranscriptAlleleStart(
                FuncotatorUtils.getTranscriptAlleleStartPosition(variant, transcript.getExons(), transcript.getGenomicStrand())
        );

        // Get the coding region start position (in the above computed transcript coding region):
        sequenceComparison.setCodingSequenceAlleleStart(
                FuncotatorUtils.getStartPositionInTranscript(variant, exonPositionList, transcript.getGenomicStrand())
        );

        // Get the overlapping exon start / stop as an interval from the given variant:
        //TODO: See the overlap detector for this:
        sequenceComparison.setExonPosition(
                FuncotatorUtils.getOverlappingExonPositions(refAllele, altAllele, variant.getContig(), variant.getStart(), variant.getEnd(), transcript.getGenomicStrand(), exonPositionList)
        );

        // Get the in-frame start position of the codon containing the given variant:
        sequenceComparison.setAlignedCodingSequenceAlleleStart(FuncotatorUtils.getAlignedPosition(sequenceComparison.getCodingSequenceAlleleStart()));

        // Get the in-frame stop position of the codon containing the given variant:
        sequenceComparison.setAlignedReferenceAlleleStop(
                FuncotatorUtils.getAlignedEndPosition(
                        // Subtract 1 because of the 1-based/inclusive nature of genetic coordinates:
                        sequenceComparison.getCodingSequenceAlleleStart() + refAllele.length() - 1
                )
        );

        // Get the in-frame/codon-aligned region containing the reference allele:
        sequenceComparison.setAlignedReferenceAllele(
                FuncotatorUtils.getAlignedRefAllele(
                        referenceBases,
                        referenceWindow,
                        refAllele,
                        altAllele,
                        sequenceComparison.getCodingSequenceAlleleStart(),
                        sequenceComparison.getAlignedCodingSequenceAlleleStart(),
                        sequenceComparison.getStrand(),
                        variant)
        );

        // Set our alternate allele:
        sequenceComparison.setAlternateAllele(altAllele.getBaseString());

        // Set our stop position:
        sequenceComparison.setAlignedAlternateAlleleStop(
                FuncotatorUtils.getAlignedEndPosition(
                        sequenceComparison.getCodingSequenceAlleleStart() + altAllele.length() - 1
                )
        );

        // Get the aligned alternate allele:
        final int alignedRefAlleleStartPos = sequenceComparison.getCodingSequenceAlleleStart() - sequenceComparison.getAlignedCodingSequenceAlleleStart() + 1;
        sequenceComparison.setAlignedAlternateAllele(
                FuncotatorUtils.getAlternateSequence(
                        new StrandCorrectedReferenceBases(sequenceComparison.getAlignedReferenceAllele(), transcript.getGenomicStrand()),
                        alignedRefAlleleStartPos,
                        refAllele,
                        altAllele,
                        sequenceComparison.getStrand())
        );

        //==============================================================================================================
        // Get the coding sequence for the transcript if we have a transcript sequence for this variant:

        if ( processSequenceInformation ) {
            if ( transcriptIdMap.containsKey(transcript.getTranscriptId()) ) {

                // Get padding bases just in case this variant is an indel and trails off the end of our transcript:
                final String transcriptTailPaddingBaseString = getTranscriptEndPaddingBases(variant, altAllele, exonPositionList, reference);

                // NOTE: This can't be null because of the Funcotator input args.
                final String rawCodingSequence = getCodingSequenceFromTranscriptFasta(
                        transcript.getTranscriptId(),
                        transcriptIdMap,
                        transcriptFastaReferenceDataSource,
                        transcriptTailPaddingBaseString
                );

                // Now that we have our transcript sequence, we must make sure that our reference allele is in it
                // correctly.
                // This is because if the user specifies a ref allele that is NOT the same as what is in the reference,
                // their specified allele takes precedence and overrides the allele from the reference genome.
                final String correctedCodingSequence;

                // We can't yet handle sequences that overrun the end of the coding sequence (Issue 4307 - https://github.com/broadinstitute/gatk/issues/4307):
                if ( (sequenceComparison.getCodingSequenceAlleleStart() - 1 + refAllele.getBaseString().length()) > rawCodingSequence.length() ) {
                    throw new FuncotatorUtils.TranscriptCodingSequenceException("Reference allele runs off end of coding sequence.  Cannot yet handle this case.");
                }
                else {
                    correctedCodingSequence = rawCodingSequence.substring(0, sequenceComparison.getCodingSequenceAlleleStart() - 1) +
                            refAllele.getBaseString() +
                            rawCodingSequence.substring(sequenceComparison.getCodingSequenceAlleleStart() + refAllele.length() - 1);
                }
                // Get the transcript sequence as described by the given exonPositionList:
                sequenceComparison.setTranscriptCodingSequence(new ReferenceSequence(transcript.getTranscriptId(), transcript.getStart(), correctedCodingSequence.getBytes()));

                // Get the in-frame/codon-aligned CODING region containing the reference allele:
                // NOTE: We are calling this with Strand.POSITIVE because we have already reverse complemented the reference sequence.
                sequenceComparison.setAlignedCodingSequenceReferenceAllele(
                        FuncotatorUtils.getAlignedCodingSequenceAllele(
                                sequenceComparison.getTranscriptCodingSequence().getBaseString(),
                                sequenceComparison.getAlignedCodingSequenceAlleleStart(),
                                sequenceComparison.getAlignedReferenceAlleleStop(),
                                refAllele,
                                sequenceComparison.getCodingSequenceAlleleStart(),
                                Strand.POSITIVE)
                );

                // Get the aligned coding sequence alternate allele:
                sequenceComparison.setAlignedCodingSequenceAlternateAllele(
                        FuncotatorUtils.getAlternateSequence(
                                new StrandCorrectedReferenceBases(sequenceComparison.getAlignedCodingSequenceReferenceAllele(), transcript.getGenomicStrand()),
                                alignedRefAlleleStartPos,
                                refAllele,
                                altAllele,
                                sequenceComparison.getStrand())
                );

                final ProteinChangeInfo proteinChangeInfo = ProteinChangeInfo.create(
                        refAllele,
                        altAllele,
                        sequenceComparison.getCodingSequenceAlleleStart(),
                        sequenceComparison.getAlignedCodingSequenceAlleleStart(),
                        correctedCodingSequence,
                        sequenceComparison.getStrand(),
                        // Figure out if we are in a mitochondrial contig:
                        // TODO: Make this more robust by detecting the mito contig based on the reference used.  (issue https://github.com/broadinstitute/gatk/issues/5364).
                        FuncotatorConstants.MITOCHONDRIAL_CONTIG_NAMES.contains(variant.getContig())
                );

                // Set our protein change:
                sequenceComparison.setProteinChangeInfo( proteinChangeInfo );
            }
            else {
                logger.warn("Attempted to process transcript information for transcript WITHOUT sequence data.  Ignoring sequence information for Gencode Transcript ID: " + transcript.getTranscriptId());
            }
        }

        //=============================================================================================================

        return sequenceComparison;
    }

    @VisibleForTesting
    static String getTranscriptEndPaddingBases(final VariantContext variant,
                                               final Allele altAllele,
                                               final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                               final ReferenceContext reference) {

        // We pad the end of the transcript to allow for the case where an indel runs off the end of a transcript
        // and needs to be annotated.
        // One variant where this happens is:
        //   <B37 Ref>:  1:178514560 A->AT
        //
        // We need only do this when the variant is close to the end of the transcript.
        //
        // Unfortunately we need to pad the end by the number of bases beyond the boundary, rounded up to the
        // next codon end position.
        // This corresponds to (with an extra codon for safety):
        //        (Math.ceil(<number of inserted bases>/AminoAcid.CODON_LENGTH)+1)*AminoAcid.CODON_LENGTH
        // This is a problem because transcriptFastaReferenceDataSource has only the bases in a given transcript.
        // Because of this we need to go to the real reference sequence and grab additional bases to pad onto the
        // end of the transcript coding sequence.

        final int transcriptEndGenomicPosition = exonPositionList.get(exonPositionList.size()-1).getEnd();
        // Add one because of inclusive positions:
        final int basesToTranscriptEnd = transcriptEndGenomicPosition - variant.getStart() + 1;

        final byte[] transcriptTailPaddingBases;
        if ( (variant.getType() == VariantContext.Type.INDEL) &&
                (basesToTranscriptEnd < TRANSCRIPT_END_WINDOW_PADDING_THRESHOLD) ) {
            final int numIndelBases = Math.abs(variant.getReference().length() - altAllele.length());
            final int numPaddingBases = (int)((Math.ceil(numIndelBases/((double)AminoAcid.CODON_LENGTH))+1)*AminoAcid.CODON_LENGTH);

            // Get extra bases from the reference:
            transcriptTailPaddingBases = reference.getBases(new SimpleInterval(reference.getWindow().getContig(), transcriptEndGenomicPosition+1, transcriptEndGenomicPosition + numPaddingBases));
        }
        else {
            // No bases needed:
            transcriptTailPaddingBases = new byte[]{};
        }

        return new String(transcriptTailPaddingBases);
    }

    /**
     * Calculates the fraction of Guanine and Cytosine bases in a window of a given size around a variant.
     * Note: Since Guanine and Cytosine are complementary bases, strandedness makes no difference.
     * @param refAllele Reference {@link Allele} for the locus in question.
     * @param altAllele Alternate {@link Allele} for the locus in question.
     * @param referenceContext The {@link ReferenceContext} for a variant.  Assumed to already be centered on the variant of interest.  Must not be {@code null}.
     * @param windowSize The number of bases to the left and right of the given {@code variant} to calculate the GC Content.  Must be >=1.
     * @return The fraction of Guanine and Cytosine bases / total bases in a window of size {@code windowSize} around a variant.
     */
    public static double calculateGcContent( final Allele refAllele,
                                             final Allele altAllele,
                                             final ReferenceContext referenceContext,
                                             final int windowSize ) {

        // TODO: this seems to do something similar to FuncotatorUtils::getBasesInWindowAroundReferenceAllele - should this method call into that?

        Utils.nonNull( referenceContext );
        ParamUtils.isPositive( windowSize, "Window size must be >= 1." );

        final int leadingWindowSize;
        final int trailingWindowSize = windowSize;

        if ( GATKVariantContextUtils.isInsertion(refAllele, altAllele) ||
                GATKVariantContextUtils.isDeletion(refAllele, altAllele)) {
            // If we have an insertion, we take 1 less base from the front
            // because the insertion happens between two codons.
            // The preceding padding base is there as a convenience in VCF files.
            // Thus the prior <windowSize> bases will contain this leading padding base.

            // If we have a deletion, the convention in VCF files is to include a
            // padding base at the front prior to the deleted bases so the alternate
            // allele can be non-empty.
            // Because of this we subtract 1 from the leading window size.

            leadingWindowSize = windowSize - 1;
        }
        else {
            leadingWindowSize = windowSize;
        }

        // Get the bases:
        byte[] bases = referenceContext.getBases(leadingWindowSize, trailingWindowSize);

        // Get the gcCount:
        long gcCount = 0;
        for ( final byte base : bases ) {
            if ( BaseUtils.basesAreEqual(base, BaseUtils.Base.G.base) || BaseUtils.basesAreEqual(base, BaseUtils.Base.C.base) ) {
                ++gcCount;
            }
        }

        // Calculate the ratio itself:
        return ((double)gcCount) / ((double) bases.length);
    }

    /**
     * Creates a {@link GencodeFuncotationBuilder} with some of the fields populated.
     *
     * @param variant The {@link VariantContext} for the current variant.
     * @param altAllele The alternate {@link Allele} we are currently annotating.
     * @param transcript The current {@link GencodeGtfTranscriptFeature} containing our {@code alternateAllele}.
     * @return A trivially populated {@link GencodeFuncotationBuilder} object.
     */
     private GencodeFuncotationBuilder createGencodeFuncotationBuilderWithTrivialFieldsPopulated(final VariantContext variant,
                                                                                                final Allele altAllele,
                                                                                                final GencodeGtfTranscriptFeature transcript) {
        return createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, transcript, this.ncbiBuildVersion);
     }

    /**
     * Creates a {@link GencodeFuncotationBuilder} with some of the fields populated.
     *
     * @param variant The {@link VariantContext} for the current variant.
     * @param altAllele The alternate {@link Allele} we are currently annotating.
     * @param transcript The current {@link GencodeGtfTranscriptFeature} containing our {@code alternateAllele}.
     * @param ncbiBuildVersion NCBI build version
     * @return A trivially populated {@link GencodeFuncotationBuilder} object.
     */
     @VisibleForTesting
     static GencodeFuncotationBuilder createGencodeFuncotationBuilderWithTrivialFieldsPopulated(final VariantContext variant,
                                                                                                final Allele altAllele,
                                                                                                final GencodeGtfTranscriptFeature transcript,
                                                                                                final String ncbiBuildVersion) {

         //TODO: Add gtfFeature.getGeneType() as an annotation field in the GencodeFuncotation - Issue #4408

         final GencodeFuncotationBuilder gencodeFuncotationBuilder = new GencodeFuncotationBuilder();

         gencodeFuncotationBuilder
                 .setRefAllele(variant.getReference())
                 .setStrand(transcript.getGenomicStrand())
                 .setHugoSymbol(transcript.getGeneName())
                 .setNcbiBuild(ncbiBuildVersion)
                 .setChromosome(transcript.getChromosomeName())
                 .setStart(variant.getStart())
                 .setGeneTranscriptType(transcript.getTranscriptType());

         // The end position is inclusive, so we need to make sure we don't double-count the start position (so we subtract 1):
         gencodeFuncotationBuilder
                 .setEnd(variant.getEnd())
                 .setVariantType(getVariantType(variant.getReference(), altAllele))
                 .setTumorSeqAllele2(altAllele.getBaseString())
                 .setGenomeChange(getGenomeChangeString(variant, altAllele))
                 .setAnnotationTranscript(transcript.getTranscriptId());

         // Check for the optional non-serialized values for sorting:
         // NOTE: For ENSEMBL gtf files, this field is optional.
         // NOTE 2: This is kind of a kludge:
         if ( transcript.getLocusLevel() == null ) {
             gencodeFuncotationBuilder.setLocusLevel(0);
         }
         else {
             gencodeFuncotationBuilder.setLocusLevel(Integer.valueOf(transcript.getLocusLevel().toString()));
         }

        // Check for existence of Appris Rank and set it:
         gencodeFuncotationBuilder.setApprisRank( getApprisRank( transcript ) );

         // Get the length of the transcript:
         // NOTE: We add 1 because of genomic cordinates:
        gencodeFuncotationBuilder.setTranscriptLength(
                transcript.getExons().stream()
                .mapToInt(Locatable::getLengthOnReference)
                .sum()
        );

         return gencodeFuncotationBuilder;
    }

    /**
     * Determines if the given UTR is 5' of the given transcript.
     * Assumes the UTR is part of the given transcript.
     * Essentially this will check to see if the UTR comes before or after the start codon.
     * @param utr The {@link GencodeGtfUTRFeature} to check for relative location in the given {@link GencodeGtfTranscriptFeature}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which to check for the given {@code utr}.
     * @return {@code true} if the given {@code utr} is 5' for the given {@code transcript}; {@code false} otherwise.
     */
    @VisibleForTesting
    static boolean is5PrimeUtr(final GencodeGtfUTRFeature utr, final GencodeGtfTranscriptFeature transcript) {
        final boolean isBefore;

        GencodeGtfStartCodonFeature startCodon = null;
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
            if ( exon.getStartCodon() != null ) {
                startCodon = exon.getStartCodon();
                break;
            }
        }

        if ( transcript.getGenomicStrand() == Strand.POSITIVE ) {
            isBefore = (startCodon != null) &&
                    ((utr.getStart() < startCodon.getStart()) &&
                    (utr.getEnd() < startCodon.getStart()));
        }
        else {
            isBefore = (startCodon != null) &&
                    ((utr.getStart() > startCodon.getEnd()) &&
                    (utr.getEnd() > startCodon.getEnd()));
        }


        return isBefore;
    }

    /**
     * Determines whether the provided variant is in the 5' flanking region of the provided transcript. A variant
     * is in the 5' flanking region if it overlaps any of the fivePrimeFlankSize bases before the 5' end of the
     * transcript, but does not overlap any part of the transcript itself.
     *
     * @param variant variant to check
     * @param transcript transcript whose flanking regions to consider
     * @param fivePrimeFlankSize size in bases of the 5' flanking region
     * @return true if the variant overlaps the 5' flanking region and does not overlap the transcript itself,
     *         otherwise false
     */
    @VisibleForTesting
    static boolean isFivePrimeFlank(final VariantContext variant, final GencodeGtfTranscriptFeature transcript, final int fivePrimeFlankSize) {
        return (transcript.getGenomicStrand() == Strand.POSITIVE && isInTranscriptLeftFlank(variant, transcript, fivePrimeFlankSize)) ||
               (transcript.getGenomicStrand() == Strand.NEGATIVE && isInTranscriptRightFlank(variant, transcript, fivePrimeFlankSize));
    }

    /**
     * Determines whether the provided variant is in the 3' flanking region of the provided transcript. A variant
     * is in the 3' flanking region if it overlaps any of the threePrimeFlankSize bases after the 3' end of the
     * transcript, but does not overlap any part of the transcript itself.
     *
     * @param variant variant to check
     * @param transcript transcript whose flanking regions to consider
     * @param threePrimeFlankSize size in bases of the 3' flanking region
     * @return true if the variant overlaps the 3' flanking region and does not overlap the transcript itself,
     *         otherwise false
     */
    @VisibleForTesting
    static boolean isThreePrimeFlank(final VariantContext variant, final GencodeGtfTranscriptFeature transcript, final int threePrimeFlankSize) {
        return (transcript.getGenomicStrand() == Strand.POSITIVE && isInTranscriptRightFlank(variant, transcript, threePrimeFlankSize)) ||
               (transcript.getGenomicStrand() == Strand.NEGATIVE && isInTranscriptLeftFlank(variant, transcript, threePrimeFlankSize));
    }

    private static boolean isInTranscriptLeftFlank(final VariantContext variant, final GencodeGtfTranscriptFeature transcript, final int flankSize) {
        return variant.getContig().equals(transcript.getContig()) &&
                variant.getEnd() < transcript.getStart() &&
                transcript.getStart() - variant.getEnd() <= flankSize;
    }

    private static boolean isInTranscriptRightFlank(final VariantContext variant, final GencodeGtfTranscriptFeature transcript, final int flankSize) {
        return variant.getContig().equals(transcript.getContig()) &&
                variant.getStart() > transcript.getEnd() &&
                variant.getStart() - transcript.getEnd() <= flankSize;
    }

    /**
     * Creates a string representing the genome change given the variant, allele, and gene feature for this variant.
     * NOTE: The genome change will only reflect positions and bases within an exon.
     *       Positions beyond the ends of exons will be changed to the last position in the exon.
     *       Bases beyond the ends of exons will be truncated from the resulting string.
     * @param variant {@link VariantContext} of which to create the change.
     * @param altAllele {@link Allele} representing the alternate allele for this variant.
     * @return A short {@link String} representation of the genomic change for the given variant, allele, and feature.
     */
    @VisibleForTesting
    static String getGenomeChangeString(final VariantContext variant,
                                        final Allele altAllele) {

        // Check for insertion:
        if ( variant.getReference().length() < altAllele.length() ) {
            final String cleanAltAlleleString = FuncotatorUtils.getNonOverlappingAltAlleleBaseString( variant.getReference(), altAllele, false);

            return "g." + variant.getContig() +
                    ":" + variant.getStart() + "_" + (variant.getStart() + 1) + "ins" +
                    cleanAltAlleleString;
        }
        // Check for deletion:
        else if ( variant.getReference().length() > altAllele.length() ) {

            final String cleanAltAlleleString = FuncotatorUtils.getNonOverlappingAltAlleleBaseString(variant.getReference(), altAllele, true);

            final int startPos = variant.getStart() + 1;
            final int endPos = variant.getStart() + variant.getReference().length() - 1;

            if ( startPos == endPos ) {
                return "g." + variant.getContig() +
                        ":" + startPos + "del" + cleanAltAlleleString;
            }
            else {
                return "g." + variant.getContig() +
                        ":" + startPos + "_" + endPos +
                        "del" + cleanAltAlleleString;
            }
        }
        // Check for SNP:
        else if ( variant.getReference().length() == 1 ) {
            return "g." + variant.getContig() +
                    ":" + variant.getStart() +
                    variant.getReference().getBaseString() + ">" + altAllele.getBaseString();
        }
        else {
            return "g." + variant.getContig() +
                    ":" + variant.getStart() + "_" + ( variant.getStart() + variant.getReference().length() - 1) +
                    variant.getReference().getBaseString() + ">" + altAllele.getBaseString();
        }
    }

    /**
     * Sort the given list of funcotations such that the list becomes sorted in "best"->"worst" order by each funcotation's
     * transcript.
     * @param funcotationList The {@link List} of {@link GencodeFuncotation} to sort.
     */
    private void sortFuncotationsByTranscriptForOutput( final List<GencodeFuncotation> funcotationList ) {
        //TODO: Make this sort go from "worst" -> "best" so we can just pop the last element off and save some time.
        //TODO: Make the Comparator object a private instance variable.
        funcotationList.sort(gencodeFuncotationComparator);
    }

    /**
     * Creates a {@link List} of {@link GencodeFuncotation}s based on the given {@link VariantContext} with type
     * {@link GencodeFuncotation.VariantClassification#IGR}.
     * @param variant The {@link VariantContext} for which to create {@link Funcotation}s.
     * @param reference The {@link ReferenceContext} against which to compare the given {@link VariantContext}.
     * @return A list of IGR annotations for the given variant.
     */
    private List<GencodeFuncotation> createIgrFuncotations(final VariantContext variant, final ReferenceContext reference) {
        // for each allele, create an annotation.

        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>();

        for ( final Allele altAllele : variant.getAlternateAlleles() ) {
            gencodeFuncotations.add( createIgrFuncotation(variant, altAllele, reference) );
        }

        return gencodeFuncotations;
    }

    /**
     * Condenses a given {@link GencodeFuncotation} into a string for the `other transcripts` annotation.
     * @param funcotation The {@link GencodeFuncotation} to condense.
     * @return A {@link String} representing the given {@link GencodeFuncotation}.
     */
    private static String condenseGencodeFuncotation( final GencodeFuncotation funcotation ) {
        Utils.nonNull( funcotation );

        final StringBuilder condensedFuncotationStringBuilder = new StringBuilder();

        if ( !funcotation.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR) ) {
            condensedFuncotationStringBuilder.append(funcotation.getHugoSymbol());
            condensedFuncotationStringBuilder.append(OTHER_TRANSCRIPTS_INFO_SEP);
            condensedFuncotationStringBuilder.append(funcotation.getAnnotationTranscript());
            condensedFuncotationStringBuilder.append(OTHER_TRANSCRIPTS_INFO_SEP);
            condensedFuncotationStringBuilder.append(funcotation.getVariantClassification());

            if ( (funcotation.getProteinChange() != null) &&
                 !(funcotation.getVariantClassification().equals(GencodeFuncotation.VariantClassification.INTRON) ||
                 ((funcotation.getSecondaryVariantClassification() != null) && funcotation.getSecondaryVariantClassification().equals(GencodeFuncotation.VariantClassification.INTRON))) ) {
                condensedFuncotationStringBuilder.append(OTHER_TRANSCRIPTS_INFO_SEP);
                condensedFuncotationStringBuilder.append(funcotation.getProteinChange());
            }
        }
        else {
            //TODO: This is known issue #3849:
            condensedFuncotationStringBuilder.append("IGR_ANNOTATON");
        }

        return condensedFuncotationStringBuilder.toString();
    }

    /**
     * Creates a {@link GencodeFuncotation}s based on the given {@link Allele} with type
     * {@link GencodeFuncotation.VariantClassification#IGR}.
     * Reports reference bases as if they are on the {@link Strand#POSITIVE} strand.
     * @param variant The {@link VariantContext} associated with this annotation.
     * @param altAllele The alternate {@link Allele} to use for this {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} in which the given {@link Allele}s appear.
     * @return An IGR funcotation for the given allele.
     */
    private GencodeFuncotation createIgrFuncotation(final VariantContext variant,
                                                    final Allele altAllele,
                                                    final ReferenceContext reference){

        final GencodeFuncotationBuilder funcotationBuilder = new GencodeFuncotationBuilder();

        // Get GC Content:
        funcotationBuilder.setGcContent( calculateGcContent( variant.getReference(), altAllele, reference, gcContentWindowSizeBases ) );

        final String alleleString = altAllele.getBaseString().isEmpty() ? altAllele.toString() : altAllele.getBaseString();

        funcotationBuilder.setVariantClassification( GencodeFuncotation.VariantClassification.IGR )
                          .setRefAllele( variant.getReference() )
                          .setTumorSeqAllele2( alleleString )
                          .setStart(variant.getStart())
                          .setEnd(variant.getEnd())
                          .setVariantType(getVariantType(variant.getReference(), altAllele))
                          .setChromosome(variant.getContig())
                          .setAnnotationTranscript(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY);

        if ( (!altAllele.isSymbolic()) && (!altAllele.equals(Allele.SPAN_DEL)) ) {
            funcotationBuilder.setGenomeChange(getGenomeChangeString(variant, altAllele));
        }

        funcotationBuilder.setNcbiBuild(ncbiBuildVersion);

        // Set our reference context in the the FuncotatonBuilder:
        funcotationBuilder.setReferenceContext(FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, reference, Strand.POSITIVE, referenceWindow).getBaseString());

        // Set our version:
        funcotationBuilder.setVersion(version);

        // Set our data source name:
        funcotationBuilder.setDataSourceName(getName());

        return funcotationBuilder.build();
    }

    /**
     * Creates a {@link GencodeFuncotation}s based on the given {@link Allele} with type
     * {@link GencodeFuncotation.VariantClassification#FIVE_PRIME_FLANK} or
     * {@link GencodeFuncotation.VariantClassification#THREE_PRIME_FLANK}
     *
     * @param variant The {@link VariantContext} associated with this annotation.
     * @param altAllele The alternate {@link Allele} to use for this {@link GencodeFuncotation}.
     * @param transcript {@link GencodeGtfTranscriptFeature} associated with this flank funcotation.
     * @param reference The {@link ReferenceContext} in which the given {@link Allele}s appear.
     * @param flankType {@link GencodeFuncotation.VariantClassification#FIVE_PRIME_FLANK} or
     *                  {@link GencodeFuncotation.VariantClassification#THREE_PRIME_FLANK}
     * @return A flank funcotation for the given allele
     */
    private GencodeFuncotation createFlankFuncotation(final VariantContext variant, final Allele altAllele, final GencodeGtfTranscriptFeature transcript, final ReferenceContext reference, final GencodeFuncotation.VariantClassification flankType) {

        // Create a symbolic allele funcotation for the flank type, if appropriate:
        if (altAllele.isSymbolic() || altAllele.equals(Allele.SPAN_DEL) ) {
            return createFuncotationForSymbolicAltAllele(variant, altAllele, transcript, reference, flankType);
        }

        final GencodeFuncotationBuilder funcotationBuilder = new GencodeFuncotationBuilder();

        // NOTE: The reference context is ALWAYS from the + strand
        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, reference, Strand.POSITIVE, referenceWindow);
        final String referenceBasesString = referenceBases.getBaseString(Strand.POSITIVE);

        final double gcContent = calculateGcContent(variant.getReference(), altAllele, reference, gcContentWindowSizeBases);

        funcotationBuilder.setHugoSymbol(transcript.getGeneName())
                .setChromosome(variant.getContig())
                .setStart(variant.getStart())
                .setEnd(variant.getEnd())
                .setStrand(transcript.getGenomicStrand())
                .setVariantClassification(flankType)
                .setVariantType(getVariantType(variant.getReference(), altAllele))
                .setRefAllele(variant.getReference())
                .setTumorSeqAllele2(altAllele.getBaseString())
                .setGenomeChange(getGenomeChangeString(variant, altAllele))
                .setAnnotationTranscript(transcript.getTranscriptId())
                .setReferenceContext(referenceBasesString)
                .setGcContent(gcContent)
                .setNcbiBuild(ncbiBuildVersion)
                .setGeneTranscriptType(transcript.getTranscriptType());

        // Set our version:
        funcotationBuilder.setVersion(version);

        // Set our data source name:
        funcotationBuilder.setDataSourceName(getName());

        return funcotationBuilder.build();
    }

    /**
     * Creates a {@link GencodeFuncotation} with a variant classification of {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#COULD_NOT_DETERMINE}
     * based on a given symbolic alternate {@link Allele}.
     *
     * Reports reference bases as if they are on the {@link Strand#POSITIVE} strand.
     * @param variant The {@link VariantContext} associated with this annotation.
     * @param altSymbolicAllele The alternate symbolic {@link Allele} associated with this annotation.
     * @param annotationTranscript The transcript in which this variant occurs.
     * @param reference The {@link ReferenceContext} in which the given {@link Allele}s appear.
     * @return An appropriate funcotation for the given symbolic allele.
     */
    private GencodeFuncotation createFuncotationForSymbolicAltAllele(final VariantContext variant,
                                                                     final Allele altSymbolicAllele,
                                                                     final GencodeGtfTranscriptFeature annotationTranscript,
                                                                     final ReferenceContext reference) {
        return createFuncotationForSymbolicAltAllele( variant, altSymbolicAllele, annotationTranscript, reference, GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE );
    }

    /**
     * Creates a {@link GencodeFuncotation} based on a given symbolic alternate {@link Allele}.
     *
     * Reports reference bases as if they are on the {@link Strand#POSITIVE} strand.
     * @param variant The {@link VariantContext} associated with this annotation.
     * @param altSymbolicAllele The alternate symbolic {@link Allele} associated with this annotation.
     * @param annotationTranscript The transcript in which this variant occurs.
     * @param reference The {@link ReferenceContext} in which the given {@link Allele}s appear.
     * @param variantClassification A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} to give to the resulting {@link GencodeFuncotation}.
     * @return An appropriate funcotation for the given symbolic allele.
     */
    private GencodeFuncotation createFuncotationForSymbolicAltAllele(final VariantContext variant,
                                                                     final Allele altSymbolicAllele,
                                                                     final GencodeGtfTranscriptFeature annotationTranscript,
                                                                     final ReferenceContext reference,
                                                                     final GencodeFuncotation.VariantClassification variantClassification) {

        logger.warn("Cannot create complete funcotation for variant at " + variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + " due to alternate allele: " + altSymbolicAllele.toString());

        final GencodeFuncotationBuilder funcotationBuilder = new GencodeFuncotationBuilder();

        // Get GC Content:
        funcotationBuilder.setGcContent( calculateGcContent( variant.getReference(), altSymbolicAllele, reference, gcContentWindowSizeBases ) );

        final String alleleString = altSymbolicAllele.getBaseString().isEmpty() ? altSymbolicAllele.toString() : altSymbolicAllele.getBaseString();

        funcotationBuilder.setVariantClassification( variantClassification )
                .setRefAllele( variant.getReference() )
                .setStrand(Strand.POSITIVE)
                .setTumorSeqAllele2( alleleString )
                .setStart(variant.getStart())
                .setEnd(variant.getEnd())
                .setVariantType(getVariantType(variant.getReference(), altSymbolicAllele))
                .setChromosome(variant.getContig())
                .setOtherTranscripts(Collections.emptyList())
                .setAnnotationTranscript(annotationTranscript.getTranscriptId());

        funcotationBuilder.setHugoSymbol( annotationTranscript.getGeneName() );

        funcotationBuilder.setNcbiBuild(ncbiBuildVersion);

        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(variant.getReference(), reference, Strand.POSITIVE, referenceWindow);

        // Set our reference context in the the FuncotatonBuilder:
        funcotationBuilder.setReferenceContext( referenceBases.getBaseString() );

        // Set our version:
        funcotationBuilder.setVersion(version);

        // Set our data source name:
        funcotationBuilder.setDataSourceName(getName());

        return funcotationBuilder.build();
    }

    /**
     * Creates a {@link GencodeFuncotation} based on a given variant and alt allele.
     * Either the reference or alternate allele contains at least one masked base.
     *
     * Reports reference bases as if they are on the {@link Strand#POSITIVE} strand.
     * @param variant The {@link VariantContext} associated with this annotation.
     * @param altAllele The alternate {@link Allele} associated with this annotation.
     * @param annotationTranscript The transcript in which this variant occurs.
     * @param reference The {@link ReferenceContext} in which the given {@link Allele}s appear.
     * @return An IGR funcotation for the given allele.
     */
    private GencodeFuncotation createFuncotationForMaskedBases(final VariantContext variant,
                                                             final Allele altAllele,
                                                             final GencodeGtfTranscriptFeature annotationTranscript,
                                                             final ReferenceContext reference) {

        final StringBuilder alleleLogStringBuilder = new StringBuilder();
        if ( variant.getReference().getBaseString().contains(FuncotatorConstants.MASKED_ANY_BASE_STRING) ) {
            alleleLogStringBuilder.append( " ref allele: " + variant.getReference().toString() );
        }
        if ( altAllele.getBaseString().contains(FuncotatorConstants.MASKED_ANY_BASE_STRING) ) {
            alleleLogStringBuilder.append( " alt allele: " + altAllele.toString() );
        }
        logger.warn("Cannot create complete funcotation for variant at " + variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + " due to" + alleleLogStringBuilder.toString() );

        final GencodeFuncotationBuilder funcotationBuilder = new GencodeFuncotationBuilder();

        // Get GC Content:
        funcotationBuilder.setGcContent( calculateGcContent( variant.getReference(), altAllele, reference, gcContentWindowSizeBases ) );

        final String alleleString = altAllele.getBaseString().isEmpty() ? altAllele.toString() : altAllele.getBaseString();

        funcotationBuilder.setVariantClassification( GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE )
                .setRefAllele( variant.getReference() )
                .setStrand(Strand.POSITIVE)
                .setTumorSeqAllele2( alleleString )
                .setStart(variant.getStart())
                .setEnd(variant.getEnd())
                .setVariantType(getVariantType(variant.getReference(), altAllele))
                .setChromosome(variant.getContig())
                .setOtherTranscripts(Collections.emptyList())
                .setAnnotationTranscript(annotationTranscript.getTranscriptId());

        funcotationBuilder.setHugoSymbol( annotationTranscript.getGeneName() );

        funcotationBuilder.setNcbiBuild(ncbiBuildVersion);

        final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(variant.getReference(), reference, Strand.POSITIVE, referenceWindow);

        // Set our reference context in the the FuncotatonBuilder:
        funcotationBuilder.setReferenceContext(referenceBases.getBaseString(Strand.POSITIVE));

        // Set our version:
        funcotationBuilder.setVersion(version);

        // Set our data source name:
        funcotationBuilder.setDataSourceName(getName());

        return funcotationBuilder.build();
    }

    /**
     * Determines the variant type based on the given reference allele and alternate allele.
     * @param refAllele The reference {@link Allele} for this variant.
     * @param altAllele The alternate {@link Allele} for this variant.
     * @return A {@link GencodeFuncotation.VariantType} representing the variation type between the given reference and alternate {@link Allele}.
     * Spanning deletions and no calls will get a type of {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType#NA}
     */
    @VisibleForTesting
    static GencodeFuncotation.VariantType getVariantType( final Allele refAllele, final Allele altAllele ) {

        if ( altAllele.length() > refAllele.length() ) {
            return GencodeFuncotation.VariantType.INS;
        }
        else if (altAllele.equals(Allele.SPAN_DEL) || altAllele.equals(Allele.NO_CALL) || altAllele.isSymbolic() ) {
            return GencodeFuncotation.VariantType.NA;
        }
        else if (altAllele.length() < refAllele.length()) {
            return GencodeFuncotation.VariantType.DEL;
        }
        else {
            // We know they are the same length, now we just need to check one of them:
            switch (refAllele.length()) {
                case 1:  return GencodeFuncotation.VariantType.SNP;
                case 2:  return GencodeFuncotation.VariantType.DNP;
                case 3:  return GencodeFuncotation.VariantType.TNP;
                default: return GencodeFuncotation.VariantType.ONP;
            }
        }
    }

    /**
     * Get the Appris Rank from the given {@link GencodeGtfGeneFeature}.
     * Appris ranks are specified as annotations using {@link GencodeGTFFieldConstants.FeatureTag}s.
     * @param gtfFeature The {@link GencodeGtfTranscriptFeature} from which to get the Appris Rank.
     * @return The highest Appris Rank found in the given {@code gtfFeature}; if no Appris Rank exists, {@code null}.
     */
    @VisibleForTesting
    static GencodeGTFFieldConstants.FeatureTag getApprisRank(final GencodeGtfTranscriptFeature gtfFeature ) {

        // Get the Appris Rank tags and convert them to Sortable Enums:
        final List<GencodeGTFFieldConstants.FeatureTag> gtfApprisTags = gtfFeature.getOptionalFields().stream()
                .filter( f -> f.getName().equals("tag") )
                .filter( f -> apprisRanks.containsKey( f.getValue() ) )
                .map( f -> apprisRanks.get(f.getValue()) ).collect(Collectors.toList());

        if ( gtfApprisTags.isEmpty() ) {
            return null;
        }
        else if ( gtfApprisTags.size() == 1 ) {
            return gtfApprisTags.get(0);
        }
        else {
            // This case should never happen, but just in case we take the highest Appris Rank:
            gtfApprisTags.sort( Comparator.naturalOrder() );
            return gtfApprisTags.get(0);
        }
    }

    /**
     * Converts a given GeneTranscriptType {@link String} to a {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification}.
     * Assumes the given {@code type} is not {@link GencodeGTFFieldConstants.KnownGeneBiotype#PROTEIN_CODING}.
     * If no type can be assessed, returns {@code null}.
     * @param type A {@link String} representing a GeneTranscriptType to convert to a {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification}.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} representing the given GeneTranscriptType {@link String}, or {@code null}.
     */
    private static GencodeFuncotation.VariantClassification convertGeneTranscriptTypeToVariantClassification (final String type ) {

        //TODO: This all needs to be fixed so there is a 1:1 mapping of GencodeGtfFeature.KnownGeneBiotype->VariantClassification - Issue #4405
        if (GencodeGTFFieldConstants.KnownGeneBiotype.LINCRNA.toString().equals(type) ||
                GencodeGTFFieldConstants.KnownGeneBiotype.MACRO_LNCRNA.toString().equals(type)) {
			return GencodeFuncotation.VariantClassification.LINCRNA;
		}
        return GencodeFuncotation.VariantClassification.RNA;
    }

    //==================================================================================================================
    // Helper Data Types:

    /**
     * A simple data object class to hold information about the transcripts in the
     * GENCODE transcript FASTA file.
     */
    @VisibleForTesting
    static class MappedTranscriptIdInfo {
        /**
         * The key in the GENCODE transcript FASTA file to use to get the coding sequence associated with this Transcript.
         */
        String mapKey;

        /**
         * The start position (1-based, inclusive) of the coding sequence in this transcript.
         */
        int codingSequenceStart;

        /**
         * The start position (1-based, inclusive) of the coding sequence in this transcript.
         */
        int codingSequenceEnd;

        /**
         * Whether or not the transcript has a 3' UTR.
         */
        boolean has3pUtr;

        /**
         * The start position (1-based, inclusive) of the 3' UTR in this transcript.
         */
        int threePrimeUtrStart;

        /**
         * The end position (1-based, inclusive) of the 3' UTR in this transcript.
         */
        int threePrimeUtrEnd;

        /**
         * Whether or not the transcript has a 3' UTR.
         */
        boolean has5pUtr;

        /**
         * The start position (1-based, inclusive) of the 5' UTR in this transcript.
         */
        int fivePrimeUtrStart;

        /**
         * The end position (1-based, inclusive) of the 5' UTR in this transcript.
         */
        int fivePrimeUtrEnd;
    }

    //==================================================================================================================
    // Segment annotation specific Methods

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFieldsForSegments() {
        return segmentMetadata.retrieveAllHeaderInfo().stream().map(VCFInfoHeaderLine::getID)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    private List<String> getSupportedFuncotationFieldsForSegmentsAsList() {
        return segmentMetadata.retrieveAllHeaderInfo().stream().map(VCFInfoHeaderLine::getID)
                .collect(Collectors.toList());
    }


    /**
     * Return the metadata for segment annotations.
     *
     * @return Never {@code null}
     */
    private FuncotationMetadata createSegmentFuncotationMetadata() {
        return VcfFuncotationMetadata.create(
                Arrays.asList(
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + GENES_SUFFIX,1, VCFHeaderLineType.String, "The genes overlapping the segment.  Blank if none."),
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + START_GENE_SUFFIX,1, VCFHeaderLineType.String, "The genes overlapping the start of the segment.  Blank if none."),
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + END_GENE_SUFFIX,1, VCFHeaderLineType.String, "The genes overlapping the end of the segment.  Blank if none."),
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + START_EXON_SUFFIX,1, VCFHeaderLineType.String, "The genes overlapping the start of the segment.  Blank if none."),
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + END_EXON_SUFFIX,1, VCFHeaderLineType.String, "The genes overlapping the end of the segment.  Blank if none."),
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + "_alt_allele",1, VCFHeaderLineType.String, "Always blank.  Included for legacy reasons."),
                        new VCFInfoHeaderLine(getName() + "_" + getVersion() + "_ref_allele",1, VCFHeaderLineType.String, "Always blank.  Included for legacy reasons.")
                )
        );
    }


    @Override
    public boolean isSupportingSegmentFuncotation() {
        return isSegmentFuncotationEnabled;
    }

    /**
     * Create Funcotations for a given segment variant context.  Will not support the {@link TranscriptSelectionMode#ALL} mode,
     * since that is undefined for segment annotations.
     *
     * @param segmentVariantContext Should be created with {@link AnnotatedIntervalToSegmentVariantContextConverter#convert(AnnotatedInterval, ReferenceContext)},
     *                              since this will adhere to some assumed conventions.  Never {@code null}
     * @param referenceContext Never {@code null}
     * @param gencodeGtfGeneFeaturesAsFeatures These must be GencodeGeneGtfFeatures.  Assumed to overlap the entire segment, not just the endpoints.  Never {@code null}
     * @return Newly created funcotations for Gencode information on the segment.  Never {@code null}
     */
    @Override
    public List<Funcotation> createFuncotationsOnSegment(final VariantContext segmentVariantContext, final ReferenceContext referenceContext, final List<Feature> gencodeGtfGeneFeaturesAsFeatures) {

        Utils.validateArg(this.transcriptSelectionMode != TranscriptSelectionMode.ALL, "Cannot create funcotations on segments if the selection mode is " + TranscriptSelectionMode.ALL);

        final List<GencodeGtfGeneFeature> geneFeatures = gencodeGtfGeneFeaturesAsFeatures.stream().map(g -> (GencodeGtfGeneFeature) g).collect(Collectors.toList());
        return createSegmentFuncotations(segmentVariantContext, referenceContext, geneFeatures, gencodeFuncotationComparator);
    }

    // Create funcotations for the start and endpoints of a segment and then create the segment funcotations from those.
    //  Assumes the reference context here is the same as the segment even if the contigs are different (ala hg19 vs b37)
    //   and uses the one in the reference context for segment breakpoints.
    private List<Funcotation> createSegmentFuncotations(final VariantContext segmentVariantContext, final ReferenceContext referenceContext, final List<GencodeGtfGeneFeature> geneFeatures, final Comparator<GencodeFuncotation> comparator) {
        // Create a funcotation for the start position of the segment.  Using the ref allele as the ref and the alt.
        final List<GencodeGtfTranscriptFeature> allBasicOverlappingTranscripts = geneFeatures.stream().flatMap(gf -> gf.getTranscripts().stream())
                .filter(GencodeFuncotationFactory::isBasic)
                .filter(t -> t.overlaps(segmentVariantContext)).collect(Collectors.toList());

        // Get the genes funcotation field
        final List<String> genes = retrieveGeneNamesFromTranscripts(allBasicOverlappingTranscripts).stream().sorted().collect(Collectors.toList());

        // Get the segment endpoints as variant contexts (assume the alternate allele is a dummy)
        final VariantContext segStartAsVariant = createSubSegmentAsVariantContext(segmentVariantContext, segmentVariantContext.getStart(), segmentVariantContext.getStart());
        final SimpleInterval segStartAsVariantInterval = new SimpleInterval(referenceContext.getInterval().getContig(), segStartAsVariant.getStart(), segStartAsVariant.getEnd());
        final List<GencodeGtfTranscriptFeature> transcriptsOverlappingStart = subsetToOverlappingTranscripts(segStartAsVariant, allBasicOverlappingTranscripts);

        final VariantContext segEndAsVariant = createSubSegmentAsVariantContext(segmentVariantContext, segmentVariantContext.getEnd(), segmentVariantContext.getEnd());
        final SimpleInterval segEndAsVariantInterval = new SimpleInterval(referenceContext.getInterval().getContig(), segEndAsVariant.getStart(), segEndAsVariant.getEnd());
        final List<GencodeGtfTranscriptFeature> transcriptsOverlappingEnd = subsetToOverlappingTranscripts(segEndAsVariant, allBasicOverlappingTranscripts);

        // Create funcotations for start of segment
        //  Make sure to alter the reference context as well to speed up the funcotation rendering.
        final List<GencodeFuncotation> startGencodeFuncotations = createFuncotationsHelper(segStartAsVariant,
                segStartAsVariant.getReference(),
                new ReferenceContext(referenceContext, segStartAsVariantInterval), transcriptsOverlappingStart);
        startGencodeFuncotations.sort(comparator);
        final GencodeFuncotation startFuncotation = startGencodeFuncotations.size() == 0 ? null: startGencodeFuncotations.get(0);

        // Create funcotations for end of segment
        //  Make sure to alter the reference context as well to speed up the funcotation rendering.
        final List<GencodeFuncotation> endGencodeFuncotations = createFuncotationsHelper(segEndAsVariant,
                segEndAsVariant.getReference(),
                new ReferenceContext(referenceContext, segEndAsVariantInterval), transcriptsOverlappingEnd);
        endGencodeFuncotations.sort(comparator);
        final GencodeFuncotation endFuncotation = endGencodeFuncotations.size() == 0 ? null: endGencodeFuncotations.get(0);

        // Remember that the start funcotation could be null (or otherwise not have a transcript)
        final GencodeGtfTranscriptFeature chosenTranscriptStart = startFuncotation != null ?
                findFirstTranscriptMatch(transcriptsOverlappingStart, startGencodeFuncotations.get(0).getAnnotationTranscript()) :
                null;
        final GencodeGtfTranscriptFeature chosenTranscriptEnd = endFuncotation != null ?
                findFirstTranscriptMatch(transcriptsOverlappingEnd, endGencodeFuncotations.get(0).getAnnotationTranscript()) :
                null;

        return createSegmentFuncotations(segmentVariantContext, genes, startFuncotation, endFuncotation, chosenTranscriptStart, chosenTranscriptEnd);
    }

    private static List<GencodeGtfTranscriptFeature> subsetToOverlappingTranscripts(final VariantContext variant, final List<GencodeGtfTranscriptFeature> allBasicOverlappingTranscripts) {
        return allBasicOverlappingTranscripts.stream()
                .filter(tx -> tx.overlaps(variant)).collect(Collectors.toList());
    }

    private static VariantContext createSubSegmentAsVariantContext(final VariantContext segmentVariantContext, final int start, final int end) {

        return new VariantContextBuilder()
                .chr(segmentVariantContext.getContig())
                .start(start)
                .stop(end)
                .alleles(segmentVariantContext.getAlleles())
                .make();
    }

    private List<Funcotation> createSegmentFuncotations(final VariantContext segmentVariantContext, final List<String> genes, final GencodeFuncotation startFuncotation, final GencodeFuncotation endFuncotation,
                                                        final GencodeGtfTranscriptFeature chosenTranscriptStart,  final GencodeGtfTranscriptFeature chosenTranscriptEnd) {

        final String genesValue = StringUtils.join(genes, GENES_FIELD_SEPARATOR);
        final String startGeneValue = startFuncotation == null ? NO_GENE_INFO_STR : startFuncotation.getHugoSymbol();
        final String endGeneValue = endFuncotation == null ? NO_GENE_INFO_STR : endFuncotation.getHugoSymbol();
        final String startExonValue = chosenTranscriptStart == null ? NO_GENE_INFO_STR : SegmentExonUtils.determineSegmentExonPosition(chosenTranscriptStart, segmentVariantContext).getSegmentStartExonOverlap();
        final String endExonValue = chosenTranscriptEnd == null ? NO_GENE_INFO_STR : SegmentExonUtils.determineSegmentExonPosition(chosenTranscriptEnd, segmentVariantContext).getSegmentEndExonOverlap();

        // Ordering of the fields as list is in createSegmentFuncotationMetadata().  Note that this will not dictate
        //  the rendering in the final output of funcotator, since that is the purview of the OutputRenderer.
        return segmentVariantContext.getAlternateAlleles().stream().map(a -> TableFuncotation.create(
                getSupportedFuncotationFieldsForSegmentsAsList(),
                Arrays.asList(
                    genesValue, startGeneValue, endGeneValue, startExonValue, endExonValue,
                        LEGACY_SEGMENT_REF_ALLELE_VALUE, LEGACY_SEGMENT_ALT_ALLELE_VALUE
                ), a, getName(), segmentMetadata)).collect(Collectors.toList());
    }

    private static Set<String> retrieveGeneNamesFromTranscripts(final List<GencodeGtfTranscriptFeature> txs) {
        return txs.stream().map(GencodeGtfTranscriptFeature::getGeneName).collect(Collectors.toSet());
    }

    private static GencodeGtfTranscriptFeature findFirstTranscriptMatch(final List<GencodeGtfTranscriptFeature> transcripts, final String txId) {
        return transcripts.stream().filter(tx -> tx.getTranscriptId().equals(txId)).findFirst().orElse(null);
    }
}
