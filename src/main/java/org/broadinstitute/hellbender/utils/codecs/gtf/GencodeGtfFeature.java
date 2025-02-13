package org.broadinstitute.hellbender.utils.codecs.gtf;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * A {@link GencodeGtfFeature} represents data in a GENCODE GTF file.
 *
 * Features are grouped logically by related data.
 * While the abstract class {@link GencodeGtfFeature} represents a single line
 * of a GENCODE GTF File, the concrete instantiations represent at least one line,
 * and often more than one.
 *
 * For example, a {@link GencodeGtfGeneFeature} represents all lines in the given
 * data file with information on a particular gene.  This includes all transcripts,
 * exons, coding regions, etc. in that gene.
 *
 * Similarly, a {@link GencodeGtfTranscriptFeature} represents all lines in the given
 * data file with information on a particular transcript.
 *
 * However, a {@link GencodeGtfSelenocysteineFeature} represents a particular line
 * in the given data file that contains information on a specific selenocysteine.
 *
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Currently only supports GENCODE versions 19-26.
 *
 * Created by jonn on 7/21/17.
 */
public abstract class GencodeGtfFeature implements Feature, Locatable, Comparable<GencodeGtfFeature> {

    private static final Logger logger = LogManager.getLogger(GencodeGtfFeature.class);

    // ===========================================================================


    public static final String ANNOTATION_SOURCE_ENSEMBL = "ENSEMBL";
    public static final String ANNOTATION_SOURCE_HAVANA = "HAVANA";
    public static final String ANNOTATION_SOURCE_ENA = "ena";

    // ===========================================================================

    // Metadata fields:

    private static final String FIELD_DELIMITER                 = "\t";

    public static final int NO_FEATURE_ORDER                    = -1;
    public static final int NO_EXON_NUMBER                      = -1;

    private static final int NUM_FIELDS                         = 9;

    private static final int CHROMOSOME_NAME_INDEX              = 0;
    private static final int ANNOTATION_SOURCE_INDEX            = 1;
    private static final int FEATURE_TYPE_INDEX                 = 2;
    private static final int START_LOCATION_INDEX               = 3;
    private static final int END_LOCATION_INDEX                 = 4;
    private static final int GENOMIC_STRAND_INDEX               = 6;
    private static final int GENOMIC_PHASE_INDEX                = 7;
    private static final int EXTRA_FIELDS_INDEX                 = 8;

    private static final String EXTRA_FIELD_DELIMITER           = ";";

    private static final int EXTRA_FIELD_KEY_INDEX              = 0;
    private static final int EXTRA_FIELD_VALUE_INDEX            = 1;
    public static final String EXTRA_FIELD_KEY_VALUE_SPLITTER   = " ";

    private static final Pattern NUMBER_PATTERN                 = Pattern.compile("\\d\\d*");

    private String ucscGenomeVersion =  null;
    @VisibleForTesting
    final GencodeGtfFeatureBaseData baseData;

    // ================================================================================================

    /**
     * Populate this GencodeGtfFeature with the given data.
     * @param gtfFields {@link String[]} containing an ordered list of fields to use to populate this {@link GencodeGtfFeature}.
     * @param gtfFileType A {@link String} containing the file type of the GTF data that created this {@link GencodeGtfFeature}.
     */
    protected GencodeGtfFeature(final String[] gtfFields, final String gtfFileType) {

        Utils.validateArg(gtfFields.length == NUM_FIELDS, "Unexpected number of fields: " + gtfFields.length + " != " + NUM_FIELDS);

        baseData = new GencodeGtfFeatureBaseData();

        try {
            baseData.genomicPosition = new SimpleInterval(
                    gtfFields[CHROMOSOME_NAME_INDEX],
                    Integer.valueOf(gtfFields[START_LOCATION_INDEX]),
                    Integer.valueOf(gtfFields[END_LOCATION_INDEX])
            );
        }
        catch (final NumberFormatException ex) {
            throw new UserException.MalformedFile("Cannot read integer value for start/end position!");
        }

        baseData.gtfSourceFileType       = gtfFileType;

        baseData.annotationSource        = gtfFields[ANNOTATION_SOURCE_INDEX];
        baseData.featureType             = GencodeGtfFeature.FeatureType.getEnum( gtfFields[FEATURE_TYPE_INDEX].toLowerCase() );
        baseData.genomicStrand           = convertStringToStrand( gtfFields[GENOMIC_STRAND_INDEX] );
        baseData.genomicPhase            = GenomicPhase.getEnum( gtfFields[GENOMIC_PHASE_INDEX] );

        // Get the extra fields from the last column:
        final String[] extraFields    = gtfFields[EXTRA_FIELDS_INDEX].split(EXTRA_FIELD_DELIMITER, -1);

        final StringBuilder anonymousOptionalFieldBuilder = new StringBuilder();

        // Now there are "optional" fields to go through (some actually required, some actually optional),
        // But we need to match up the field names to the fields themselves:
        for ( final String extraField : extraFields ) {

            final String trimmedExtraField = extraField.trim();
            if (trimmedExtraField.isEmpty()) {
                continue;
            }

            final int splitPoint = trimmedExtraField.indexOf(EXTRA_FIELD_KEY_VALUE_SPLITTER);
            if( splitPoint == -1 ) {
                throw new UserException.MalformedFile("Extraneous optional field data - not in a key/value pair: " + extraField);
            }

            final String fieldName = trimmedExtraField.substring(0, splitPoint).trim();

            // The value of the field may be between two quotes.
            // We remove them here.
            final String rawFieldValue = trimmedExtraField.substring(splitPoint + 1, trimmedExtraField.length());
            final String fieldValue = StringUtils.remove(rawFieldValue.trim(), '"');

            if( fieldValue.contains(EXTRA_FIELD_KEY_VALUE_SPLITTER) ){
                throw new UserException("Expected a key/value pair but found several values " + fieldName + "/" + fieldValue);
            }

            OptionalField<?> optionalField = null;

            switch (fieldName) {
                // Find the right field to set:
                case GencodeGTFFieldConstants.GENE_ID:
                    baseData.geneId = fieldValue;
                    break;
                case GencodeGTFFieldConstants.TRANSCRIPT_ID:
                    baseData.transcriptId = fieldValue;
                    break;
                case GencodeGTFFieldConstants.GENE_TYPE:
                    baseData.geneType = fieldValue;
                    break;
                // For ENSEMBL GTF files:
                case GencodeGTFFieldConstants.GENE_BIOTYPE:
                    baseData.geneType = fieldValue;
                    break;
                case GencodeGTFFieldConstants.GENE_STATUS:
                    baseData.geneStatus = fieldValue;
                    break;
                case GencodeGTFFieldConstants.GENE_NAME:
                    baseData.geneName = fieldValue;
                    break;
                case GencodeGTFFieldConstants.TRANSCRIPT_TYPE:
                    baseData.transcriptType = fieldValue;
                    break;
                case GencodeGTFFieldConstants.TRANSCRIPT_BIOTYPE:
                    baseData.transcriptType = fieldValue;
                    break;
                case GencodeGTFFieldConstants.TRANSCRIPT_STATUS:
                    baseData.transcriptStatus = fieldValue;
                    break;
                case GencodeGTFFieldConstants.TRANSCRIPT_NAME:
                    baseData.transcriptName = fieldValue;
                    break;
                case GencodeGTFFieldConstants.EXON_NUMBER:
                    try {
                        baseData.exonNumber = Integer.valueOf(fieldValue);
                    }
                    catch (final NumberFormatException ex) {
                        throw new UserException.MalformedFile("Could not convert field value into integer: " + fieldValue);
                    }
                    break;
                case GencodeGTFFieldConstants.EXON_ID:
                    baseData.exonId = fieldValue;
                    break;
                case GencodeGTFFieldConstants.LEVEL:
                    baseData.locusLevel = fieldValue;
                    break;
                default:
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
            }

            // If the optional field was good, we add it:
            if ( optionalField != null ) {
                baseData.optionalFields.add(optionalField);
            }
        }
    }

    /**
     * Converts the given {@link String} into a {@link Strand}.
     * @param s {@link String} to convert into a {@link Strand}.
     * @return The {@link Strand} corresponding to {@code s}.
     */
    private static Strand convertStringToStrand( final String s ) {
        if ( s.equals("+") ) {
            return Strand.POSITIVE;
        }
        else if ( s.equals("-") ) {
            return Strand.NEGATIVE;
        }
        else {
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }

    /**
     * Populate this GencodeGtfFeature with the given data.
     */
    protected GencodeGtfFeature(final GencodeGtfFeatureBaseData baseData) {
        this.baseData = baseData;
    }

    // ================================================================================================

    /**
     * Create the appropriate {@link GencodeGtfFeature} object based on the given {@code baseData}
     * @param baseData A {@link GencodeGtfFeatureBaseData} object containing all data for a single line in a GENCODE GTF File.
     * @return A {@link GencodeGtfFeature} containing the data in {@code baseData}
     */
    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        Utils.nonNull(baseData);

        // Create our feature:
        return baseData.featureType.create(baseData);
    }

    /**
     * Create a {@link GencodeGtfFeature} based on a line from a Gencode GTF File.
     * @param gtfLine A line from a Gencode GTF File to convert into a {@link GencodeGtfFeature} object.
     * @param gtfFileType A {@link String} containing the file type of the GTF data that created this {@link GencodeGtfFeature}.
     * @return The {@link GencodeGtfFeature} representing the information in {@code gtfLine}
     */
    public static GencodeGtfFeature create(final String gtfLine, final String gtfFileType) {
        Utils.nonNull(gtfLine);
        return create(gtfLine.split(FIELD_DELIMITER), gtfFileType);
    }

    /**
     * Create a {@link GencodeGtfFeature} based on a line from a Gencode GTF File.
     * @param gtfFields A line from a Gencode GTF File split on the {@link #FIELD_DELIMITER} character.
     * @param gtfFileType A {@link String} containing the file type of the GTF data that created this {@link GencodeGtfFeature}.
     * @return The {@link GencodeGtfFeature} representing the information in {@code gtfLine}
     */
    public static GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
        Utils.nonNull(gtfFields);

        // Ensure that the input data are superficially well-formed:
        if ( gtfFields.length != GencodeGtfCodec.NUM_COLUMNS ) {
            throw new UserException.MalformedFile("Invalid number of fields in the given GENCODE line " +
                    " - Given: " + gtfFields.length + " Expected: " + GencodeGtfCodec.NUM_COLUMNS);
        }

        final FeatureType featureType = FeatureType.getEnum( gtfFields[FEATURE_TYPE_INDEX] );

        // Return our feature:
        return featureType.create(gtfFields, gtfFileType);
    }

    // ================================================================================================

    @Override
    public String getContig() {
        return baseData.genomicPosition.getContig();
    }

    @Override
    public int getStart() {
        return baseData.genomicPosition.getStart();
    }

    @Override
    public int getEnd() {
        return baseData.genomicPosition.getEnd();
    }

    // ================================================================================================

    /**
     * Get all the features from this {@link GencodeGtfFeature} itself.
     * This is useful to get any subfeatures included in this {@link GencodeGtfFeature}.
     * @return A {@link List} of the features represented in this {@link GencodeGtfFeature}.
     */
    @VisibleForTesting
    public List<GencodeGtfFeature> getAllFeatures() {
        final List<GencodeGtfFeature> list = new ArrayList<>();
        list.add(this);
        return list;
    }


    /**
     * Get all the {@link GencodeGtfFeatureBaseData} objects from this {@link GencodeGtfFeature} itself.
     * This is useful to get any subfeatures included in this {@link GencodeGtfFeature}.
     * @return A {@link List} of the {@link GencodeGtfFeatureBaseData} objects represented in this {@link GencodeGtfFeature}.
     */
    @VisibleForTesting
    List<GencodeGtfFeatureBaseData> GencodeGtfFeatureBaseData() {
        final List<GencodeGtfFeatureBaseData> baseDataParts = new ArrayList<>();
        baseDataParts.add(this.baseData);
        getAllFeatures().forEach(feature -> baseDataParts.add(feature.baseData));
        return baseDataParts;
    }

    /**
     * Serializes the base data in {@link GencodeGtfFeature} to a string.
     * @return a {@link String} representing this {@link GencodeGtfFeature}
     */
    private String serializeToStringHelper() {

        final StringBuilder stringBuilder = new StringBuilder();

        stringBuilder.append( baseData.genomicPosition.getContig() );
        stringBuilder.append( '\t' );
        stringBuilder.append( baseData.annotationSource );
        stringBuilder.append( '\t' );
        stringBuilder.append( baseData.featureType );
        stringBuilder.append( '\t' );
        stringBuilder.append( baseData.genomicPosition.getStart() );
        stringBuilder.append( '\t' );
        stringBuilder.append( baseData.genomicPosition.getEnd() );
        stringBuilder.append( "\t.\t" );
        stringBuilder.append( baseData.genomicStrand );
        stringBuilder.append( '\t' );
        stringBuilder.append( baseData.genomicPhase );
        stringBuilder.append( '\t' );

        if ( baseData.geneId != null ) {
            stringBuilder.append("gene_id \"");
            stringBuilder.append(baseData.geneId);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.transcriptId != null) {
            stringBuilder.append("transcript_id \"");
            stringBuilder.append(baseData.transcriptId);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.geneType != null ) {
            stringBuilder.append("gene_type \"");
            stringBuilder.append(baseData.geneType);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.geneStatus != null ) {
            stringBuilder.append("gene_status \"");
            stringBuilder.append(baseData.geneStatus);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.geneName != null ) {
            stringBuilder.append("gene_name \"");
            stringBuilder.append(baseData.geneName);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.transcriptType != null ) {
            stringBuilder.append("transcript_type \"");
            stringBuilder.append(baseData.transcriptType);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.transcriptStatus != null ) {
            stringBuilder.append("transcript_status \"");
            stringBuilder.append(baseData.transcriptStatus);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.transcriptName != null ) {
            stringBuilder.append("transcript_name \"");
            stringBuilder.append(baseData.transcriptName);
            stringBuilder.append( "\"; " );
        }
        if ( baseData.exonNumber != NO_EXON_NUMBER ) {
            stringBuilder.append("exon_number ");
            stringBuilder.append(baseData.exonNumber);
            stringBuilder.append( "; " );
        }
        if ( baseData.exonId != null) {
            stringBuilder.append("exon_id \"");
            stringBuilder.append(baseData.exonId);
            stringBuilder.append( "\"; ");
        }
        if (baseData.locusLevel != null) {
            stringBuilder.append("level ");
            stringBuilder.append(baseData.locusLevel);
            stringBuilder.append("; ");
        }

        // = = = = = = = = = = = = = = = = = = = = = = =

        // Output our optional fields:
        stringBuilder.append(
                baseData.optionalFields.stream().map(Object::toString).collect(Collectors.joining(" "))
        );

        return stringBuilder.toString().trim();
    }

    /**
     * Serializes all data in {@link GencodeGtfFeature} to a string.
     * This includes all subfields of child classes.
     * @return a {@link String} representing this {@link GencodeGtfFeature}
     */
    public String serializeToString() {
        final StringBuilder stringBuilder = new StringBuilder();

        final List<GencodeGtfFeature> features = getAllFeatures();
        Collections.sort( features );

        for ( final GencodeGtfFeature feature : features ) {
            stringBuilder.append( feature.serializeToStringHelper() );
            stringBuilder.append("\n");
        }

        return stringBuilder.toString().trim();
    }

    @Override
    public String toString() {
        return serializeToString();
    }

    // ================================================================================================

    public String getGtfSourceFileType() { return baseData.gtfSourceFileType; }

    public String getUcscGenomeVersion() {
        return ucscGenomeVersion;
    }

    public void setUcscGenomeVersion(final String ucscGenomeVersion) {
        this.ucscGenomeVersion = ucscGenomeVersion;
    }

    public SimpleInterval getGenomicPosition() { return baseData.genomicPosition; }

    public int getFeatureOrderNumber() { return baseData.featureOrderNumber; }

    public String getChromosomeName() {
        return baseData.genomicPosition.getContig();
    }

    public String getAnnotationSource() {
        return baseData.annotationSource;
    }

    public FeatureType getFeatureType() {
        return baseData.featureType;
    }

    public int getGenomicStartLocation() {
        return baseData.genomicPosition.getStart();
    }

    public int getGenomicEndLocation() {
        return baseData.genomicPosition.getEnd();
    }

    public Strand getGenomicStrand() {
        return baseData.genomicStrand;
    }

    public GenomicPhase getGenomicPhase() {
        return baseData.genomicPhase;
    }

    public String getGeneId() {
        return baseData.geneId;
    }

    public String getTranscriptId() {
        return baseData.transcriptId;
    }

    public String getGeneType() {
        return baseData.geneType;
    }

    public String getGeneName() {
        return baseData.geneName;
    }

    public String getTranscriptType() {
        return baseData.transcriptType;
    }

    public String getTranscriptName() {
        return baseData.transcriptName;
    }

    public String getGeneStatus() {
        return baseData.geneStatus;
    }

    public String getTranscriptStatus() {
        return baseData.transcriptStatus;
    }

    public int getExonNumber() {
        return baseData.exonNumber;
    }

    public String getExonId() {
        return baseData.exonId;
    }

    public String getLocusLevel() {
        return baseData.locusLevel;
    }

    public List<OptionalField<?>> getOptionalFields() {
        return baseData.optionalFields;
    }

    // Certian optional fields support multiple values (for example "tag"), so we need to return a list of them.
    public List<OptionalField<?>> getOptionalField(final String key) {
        final List<OptionalField<?>> optionalFields = new ArrayList<>();
        for (final OptionalField<?> optionalField : baseData.optionalFields) {
            if ( optionalField.getName().equals(key) ) {
                optionalFields.add(optionalField);
            }
        }
        return optionalFields;
    }

    /**
     * Comparable interface implementation for {@link GencodeGtfFeature}.
     *
     * Order is determined by {@link GencodeGtfFeatureBaseData#featureOrderNumber}
     *
     * @param other {@link GencodeGtfFeature} to which to compare
     * @return -1 if this < other; 0 if this == other; 1 if this > other
     */
    @Override
    public int compareTo(final GencodeGtfFeature other) {
        Utils.nonNull(other);
        return (baseData.featureOrderNumber - other.baseData.featureOrderNumber);
    }

    @Override
    public boolean equals(final Object that) {
        if (that == null) {
            return false;
        }
        else if ( this == that ) {
            return true;
        }

        boolean isEqual = that instanceof GencodeGtfFeature;
        if (isEqual) {
            final GencodeGtfFeature thatFeature = (GencodeGtfFeature) that;
            isEqual = Objects.equals(baseData, thatFeature.baseData);

            if ( isEqual ) {
                isEqual = ucscGenomeVersion.equals( thatFeature.getUcscGenomeVersion() );
            }
        }

        return isEqual;
    }

    @Override
    public int hashCode() {
        return baseData != null ? baseData.hashCode() : 0;
    }

    /**
     * Checks if {@code other} is contained within this {@link GencodeGtfFeature}.
     * Comparison is made using {@link SimpleInterval#contains(Locatable)} ala {@link GencodeGtfFeatureBaseData#genomicPosition}
     * @param other {@link Locatable} of which to check the bounds.
     * @return true if {@code other} is contained within the bounds of this {@link GencodeGtfFeature}, false otherwise.
     */
    public boolean contains(final Locatable other) {
        return baseData.genomicPosition.contains(other);
    }

    /**
     * Checks if {@code other} overlaps with this {@link GencodeGtfFeature}.
     * Comparison is made using {@link SimpleInterval#overlaps(Locatable)} ala {@link GencodeGtfFeatureBaseData#genomicPosition}
     * @param other {@link Locatable}-derived class of which to check the bounds.
     * @return true if {@code other} overlaps the bounds of this {@link GencodeGtfFeature}, false otherwise.
     */
    public boolean overlaps(final Locatable other) {
        return baseData.genomicPosition.overlaps(other);
    }

    public void setFeatureOrderNumber(final int featureOrderNumber) {
        this.baseData.featureOrderNumber = featureOrderNumber;
    }

    // ================================================================================================

    static public class OptionalField<T> {

        private String name;
        private T value;

        public OptionalField(final String name, final T value) {
            this.name = name;
            this.value = value;
        }

        public String getName() {
            return name;
        }

        public void setName(final String name) {
            this.name = name;
        }

        public T getValue() {
            return value;
        }

        public void setValue(final T value) {
            this.value = value;
        }

        @Override
        public String toString() {

            final StringBuilder sb = new StringBuilder();

            sb.append(name);
            sb.append(" ");

            // We need to do some formatting for the numbers / non-numbers in the field:
            final String valueString = value.toString();
            if ( NUMBER_PATTERN.matcher(valueString).matches() ) {
                sb.append(valueString);
                sb.append(";");
            }
            else {
                sb.append("\"");
                sb.append(valueString);
                sb.append("\";");
            }

            return sb.toString();
        }

        @Override
        public int hashCode() {
            int result = name != null ? name.hashCode() : 0;
            result = 31 * result + (value != null ? value.hashCode() : 0);
            return result;
        }

        @Override
        public boolean equals(final Object other) {

            if (other == null) {
                return false;
            }
            else if ( this == other ) {
                return true;
            }

            if ( !(other instanceof OptionalField) ) {
                return false;
            }

            final OptionalField<?> otherOptionalField = (OptionalField<?>) other;

            return (name.equals(otherOptionalField.name)) &&
                    (value.equals(otherOptionalField.value));
        }
    }

    // ================================================================================================



    // ================================================================================================

    /**
     * Keyword identifying the source of the feature, like a program
     * (e.g. Augustus or RepeatMasker) or an organization (like TAIR).
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum AnnotationSource {
        ENSEMBL,
        HAVANA,
        ena // From ENSEMBLE GTFs
    }

    /**
     * Type of the feature represented in a single line of a GENCODE GTF File.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum FeatureType {
        GENE("gene"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfGeneFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfGeneFeature.create(gtfFields, gtfFileType);
            }
        },
        TRANSCRIPT("transcript"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfTranscriptFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfTranscriptFeature.create(gtfFields, gtfFileType);
            }
        },
        SELENOCYSTEINE("Selenocysteine"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfSelenocysteineFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfSelenocysteineFeature.create(gtfFields, gtfFileType);
            }
        },
        EXON("exon"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfExonFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfExonFeature.create(gtfFields, gtfFileType);
            }
        },
        CDS("CDS"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfCDSFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfCDSFeature.create(gtfFields, gtfFileType);
            }
        },
        START_CODON("start_codon"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfStartCodonFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfStartCodonFeature.create(gtfFields, gtfFileType);
            }
        },
        STOP_CODON("stop_codon"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfStopCodonFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfStopCodonFeature.create(gtfFields, gtfFileType);
            }
        },
        UTR("UTR"){
            public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
                return GencodeGtfUTRFeature.create(baseData);
            }
            public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
                return GencodeGtfUTRFeature.create(gtfFields, gtfFileType);
            }
        };

        @SuppressWarnings("unchecked")
        private static final Map<String, FeatureType> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        FeatureType(final String serializedValue) { serialized = serializedValue; }

        @Override
        public String toString() { return serialized; }

        public static FeatureType getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }

        /**
         * Create a {@link GencodeGtfFeature} of this type given {@code baseData}
         * @param baseData The data to use to create a {@link GencodeGtfFeature}
         * @return The {@link GencodeGtfFeature} represented by the given {@code baseData}
         */
        abstract public GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData);

        /**
         * Create a {@link GencodeGtfFeature} of this type given {@code gtfFields}
         * @param gtfFields The data to use to create a {@link GencodeGtfFeature}
         * @param gtfFileType A {@link String} containing the file type of the GTF data that created this {@link GencodeGtfFeature}.
         * @return The {@link GencodeGtfFeature} represented by the given {@code gtfFields}
         */
        abstract public GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType);
    }

    /**
     * Whether the first base of the CDS segment is the first (frame 0), second (frame 1) or third (frame 2) \
     * in the codon of the ORF.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum GenomicPhase {
        ZERO("0"),
        ONE ("1"),
        TWO ("2"),
        DOT (".");

        @SuppressWarnings("unchecked")
        private static final Map<String, GenomicPhase> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        GenomicPhase(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static GenomicPhase getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }


}
