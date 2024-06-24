package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorConstants;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.SamplePairExtractor;
import org.broadinstitute.hellbender.tools.funcotator.metadata.TumorNormalPair;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * A Funcotator output renderer for writing to MAF files.
 * Complies with version 2.4 of the MAF standard found here:
 * https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4
 * Created by jonn on 12/5/17.
 */
public class MafOutputRenderer extends OutputRenderer {

    //==================================================================================================================
    // Public Static Members:

    /**
     * Version of the MAF standard that this {@link MafOutputRenderer} writes.
     */
    public static String VERSION = "2.4";

    //==================================================================================================================
    // Private Static Members:

    private static final Logger logger = LogManager.getLogger(MafOutputRenderer.class);

    private static final Set<String> HG_19_CHR_SET = new HashSet<>(Arrays.asList("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"));

    private static final List<String> ORDERED_GENCODE_VARIANT_CLASSIFICATIONS = new ArrayList<> (Arrays.asList(
            GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
            GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
            GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
            GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
            GencodeFuncotation.VariantClassification.MISSENSE.toString(),
            GencodeFuncotation.VariantClassification.NONSENSE.toString(),
            GencodeFuncotation.VariantClassification.SILENT.toString(),
            GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
            GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString(),
            GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString(),
            GencodeFuncotation.VariantClassification.START_CODON_SNP.toString(),
            GencodeFuncotation.VariantClassification.START_CODON_INS.toString(),
            GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
            GencodeFuncotation.VariantClassification.NONSTOP.toString(),
            GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
            GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
            GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
            GencodeFuncotation.VariantClassification.INTRON.toString(),
            GencodeFuncotation.VariantClassification.LINCRNA.toString()
    ));
    private static final List<String> ORDERED_MAF_VARIANT_CLASSIFICATIONS = new ArrayList<> (Arrays.asList(
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_INS.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString()),
        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString())
    ));

    //==================================================================================================================
    // Private Members:

    /**
     * Default set of columns to include in this {@link MafOutputRenderer}.
     * Order of the columns is preserved by the {@link LinkedHashMap}, while still being able to access each field via
     * the associated key.
     */
    private final LinkedHashMap<String, String> defaultMap = new LinkedHashMap<>();

    /**
     * Map for: outputFieldName -> dataSourceFieldName1, dataSourceFieldName2 ...
     *
     * This map informs how to fill out outputFields using the input funcotations.
     *
     * Strings in the value list are in order of priority for use in the outputField.
     * That is, if both dataSourceFieldName1 and dataSourceFieldName2 are present as funcotation fields, then only
     * dataSourceFieldName1 will be used as the outputField.
     */
    private final Map<String, List<String>> outputFieldNameMap = new LinkedHashMap<>();

    /** Flag to see if the header has been written to the output file yet. */
    private boolean hasWrittenHeader = false;

    /**
     * {@link Path} to the resulting MAF output file for this {@link MafOutputRenderer}.
     */
    private final Path outputFilePath;

    /**
     * {@link java.io.PrintWriter} to which to write the output MAF file.
     */
    private Writer writer;

    /**
     * Tool header information to go into the header.
     */
    private final LinkedHashSet<String> toolHeaderLines;

    /**
     * Existing header from the input file to preserve in the output.
     */
    private final VCFHeader inputFileHeader;

    /** Override annotation list from the user. */
    private final LinkedHashMap<String, String> overrideAnnotations;

    /** The tumor normal pairs discovered in the input */
    private final List<TumorNormalPair> tnPairs;

    /** The version of the reference used to create annotations that will be output by this {@link MafOutputRenderer}.*/
    private final String referenceVersion;

    /** Fields that should be removed in the final MAF file. */
    private final Set<String> excludedOutputFields;

    //==================================================================================================================
    // Constructors:

    /**
     * Create a {@link MafOutputRenderer}.  Usage for germline use cases is unsupported.
     *
     * @param outputFilePath {@link Path} to output file (must not be null).
     * @param dataSources {@link List} of {@link DataSourceFuncotationFactory} to back our annotations (must not be null).
     * @param inputFileHeader {@link VCFHeader} of input VCF file to preserve (must not be null).
     * @param unaccountedForDefaultAnnotations {@link LinkedHashMap} of default annotations that must be added (must not be null).
     * @param unaccountedForOverrideAnnotations {@link LinkedHashMap} of override annotations that must be added (must not be null).
     * @param toolHeaderLines Lines to add to the header with information about the tool (must not be null).
     * @param referenceVersion Version of the reference we're using (must not be null).
     * @param excludedOutputFields Fields that should not be rendered in the final output. Only exact name matches will be excluded (must not be null).
     * @param toolVersion The version number of the tool used to produce the MAF file (must not be null).
     */
    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final VCFHeader inputFileHeader,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                             final Set<String> toolHeaderLines,
                             final String referenceVersion,
                             final Set<String> excludedOutputFields,
                             final String toolVersion) {
        super(toolVersion);

        Utils.nonNull(outputFilePath);
        Utils.nonNull(dataSources);
        Utils.nonNull(inputFileHeader);
        Utils.nonNull(unaccountedForDefaultAnnotations);
        Utils.nonNull(unaccountedForOverrideAnnotations);
        Utils.nonNull(toolHeaderLines);
        Utils.nonNull(referenceVersion);
        Utils.nonNull(excludedOutputFields);

        // Set our internal variables from the input:
        this.outputFilePath = outputFilePath;
        this.toolHeaderLines = new LinkedHashSet<>(toolHeaderLines);
        this.inputFileHeader = inputFileHeader;
        this.dataSourceFactories = dataSources;
        this.referenceVersion = referenceVersion;

        this.tnPairs = SamplePairExtractor.extractPossibleTumorNormalPairs(this.inputFileHeader);
        if (tnPairs.size() == 0) {
            logger.warn("No tumor/normal pairs were seen, cannot populate the some of the MAF fields (e.g. t_alt_count).  Please add '##tumor_sample=<tumor_sample_name>' and (if applicable) '##normal_sample=<normal_sample_name>' to the input VCF header");
        }

        // TODO: Make this check unnecessary and use the contained information to populate the MAF entries correctly.  (https://github.com/broadinstitute/gatk/issues/4912)
        if (this.tnPairs.size() > 1) {
            throw new UserException.BadInput("Input files with more than one tumor normal pair are currently not supported.  Found: " + tnPairs.stream()
                    .map(tn -> tn.toString())
                    .collect(Collectors.joining("; ")) );
        }

        // Merge the default annotations into our manualAnnotations:
        manualAnnotations = new LinkedHashMap<>();
        manualAnnotations.putAll(unaccountedForDefaultAnnotations);

        // Handle our override annotations a little differently:
        this.overrideAnnotations = unaccountedForOverrideAnnotations;

        // Fill in our default output map:
        initializeDefaultMapWithKeys();

        // Fill in our default map for outputField -> funcotation name:
        initializeOutputFieldNameMap();

        // TODO: Make this FASTER!
        // Update our defaultMap with the dataSourceFactories:
        for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {
            for ( final String field : funcotationFactory.getSupportedFuncotationFields() ) {

                // Check if it's in our output map already:
                if ( !defaultMap.containsKey(field) ) {
                    boolean found = false;

                    // Now check our field assignment/translation map:
                    for ( final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet() ) {
                        if ( entry.getValue().contains(field) ) {
                            found = true;
                            break;
                        }
                    }
                    if ( !found ) {
                        defaultMap.put(field, FuncotatorConstants.UNKNOWN_VALUE_STRING);
                    }
                }
            }
        }

        // Update the Default map with manual annotations that we have
        // and remove those annotations from the manual map:
        // NOTE: We must check the outputFieldNameMap in case we override an aliased output field:
        final Iterator<Map.Entry<String, String>> it = manualAnnotations.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<String, String> manualAnnotation = it.next();
            if ( defaultMap.containsKey(manualAnnotation.getKey()) ) {
                defaultMap.put(manualAnnotation.getKey(), manualAnnotation.getValue());
                it.remove();
            }
            else {
                // Check the aliased output fields.
                // NOTE: A field may be aliased many times, so we should check ALL fields before
                //       removing the iterator and moving to the next value:
                boolean isAliasedName = false;
                for ( final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet() ) {
                    if ( entry.getValue().contains(manualAnnotation.getKey())) {
                        isAliasedName = true;
                        defaultMap.put( entry.getKey(), manualAnnotation.getValue() );
                    }
                }
                if (isAliasedName) {
                    it.remove();
                }
            }
        }

        // Set values for unused fields:
        defaultMap.put(MafOutputRendererConstants.FieldName_Score, MafOutputRendererConstants.UNUSED_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_BAM_File, MafOutputRendererConstants.UNUSED_STRING);

        // Open the output object:
        try {
            writer = new BufferedWriter(new OutputStreamWriter(Files.newOutputStream(outputFilePath)));
        }
        catch (final IOException ex) {
            throw new UserException("Error opening output file path: " + outputFilePath.toUri().toString(), ex);
        }

        this.excludedOutputFields = excludedOutputFields;
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void close() {
        if (!hasWrittenHeader) {
            // The alt allele can be anything here.  We just need to write the header, not any actual funcotations.
            final String dummyAltAllele = "AT";
            final LinkedHashMap<String, String> dummyMafCompliantOutputMap = createMafCompliantOutputMap(Allele.create(dummyAltAllele), Collections.emptyList());
            writeHeader(new ArrayList<>(dummyMafCompliantOutputMap.keySet()));
        }
        if ( writer != null ) {
            try {
                writer.flush();
                writer.close();
            } catch (IOException e){
                throw new UserException.CouldNotCreateOutputFile("Failed while closing the maf writer, the output may be corrupted.", e);
            }
        }
    }

    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

        if (txToFuncotationMap.getTranscriptList().size() > 1) {
            logger.warn("MAF typically does not support multiple transcripts per variant, though this should be able to render (grouped by transcript).  No user action needed.");
        }

        // Add the generated count funcotations necessary for a MAF output rendering.
        final List<Funcotation> customMafCountFuncotations = CustomMafFuncotationCreator.createCustomMafCountFields(variant, tnPairs);
        txToFuncotationMap.getTranscriptList().forEach(txId -> txToFuncotationMap.add(txId, customMafCountFuncotations));

        // Add the custom dbSNP funcotations (e.g. dbSNP validation field) to populate the MAF correctly
        for (final String txId : txToFuncotationMap.getTranscriptList()) {
            final List<Funcotation> customDbSnpFuncotations = CustomMafFuncotationCreator.createCustomMafDbSnpFields(txToFuncotationMap.get(txId));
            if (customDbSnpFuncotations.size() == 0) {
                logger.warn("No dbSNP annotations exist for this variant.  Cannot render the dbSNP fields in the MAF.  These fields will not be correct.  " + variant);
            } else {
                txToFuncotationMap.add(txId, customDbSnpFuncotations);
            }
        }

        // Loop through each alt allele in our variant:
        for ( final Allele altAllele : variant.getAlternateAlleles() ) {

            // Ignore spanning deletions.  Those have no meaning in a MAF.
            if (altAllele.equals(Allele.SPAN_DEL)) {
                continue;
            }

            for (final String txId : txToFuncotationMap.getTranscriptList()) {

                final List<Funcotation> funcotations = txToFuncotationMap.get(txId);
                final LinkedHashMap<String, String> mafCompliantOutputMap = createMafCompliantOutputMap(altAllele, funcotations);

                // Write our header if we have to:
                if (!hasWrittenHeader) {
                    // Please note that we are implicitly using the ordering of a LinkedHashMap under the hood.
                    writeHeader(new ArrayList<>(mafCompliantOutputMap.keySet()));
                }

                try {
                    // Write the output (with manual annotations at the end):
                    final Iterator<Map.Entry<String, String>> entryIterator =
                            mafCompliantOutputMap.entrySet().iterator();

                    // We should always have at least 1 field here:
                    Map.Entry<String, String> entry = entryIterator.next();
                    writeString(entry.getValue());

                    while (entryIterator.hasNext()) {
                        entry = entryIterator.next();
                        writeString(MafOutputRendererConstants.FIELD_DELIMITER);
                        writeString(entry.getValue());
                    }
                    writeLine("");
                } catch (IOException e){
                    throw new UserException.CouldNotCreateOutputFile("Error while writing maf file, cause by: " + e.getMessage(), e);
                }
            }
        }
    }

    @VisibleForTesting
    LinkedHashMap<String, String> createMafCompliantOutputMap(final Allele altAllele, final List<Funcotation> funcotations) {
        // Create our output maps:
        final LinkedHashMap<String, Object> outputMap = new LinkedHashMap<>(defaultMap);
        final LinkedHashMap<String, Object> extraFieldOutputMap = new LinkedHashMap<>();

        // Get our funcotations for this allele and add them to the output maps:
        for (final Funcotation funcotation : funcotations) {
            if (funcotation.getAltAllele().equals(altAllele)) {
                // Add all the fields from the other funcotations into the extra field output:
                for (final String field : funcotation.getFieldNames()) {
                    setField(extraFieldOutputMap, field, funcotation.getField(field));
                }
            }
        }

        // Now add in our annotation overrides so they can be aliased correctly with the outputFieldNameMap:
        extraFieldOutputMap.putAll(overrideAnnotations);

        // Go through all output fields and see if any of the names in the value list are in our extraFieldOutputMap.
        // For any that match, we remove them from our extraFieldOutputMap and add them to the outputMap with the
        // correct key.
        for (final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet()) {
            for (final String fieldName : entry.getValue()) {
                if (extraFieldOutputMap.containsKey(fieldName)) {
                    outputMap.put(entry.getKey(), extraFieldOutputMap.remove(fieldName));
                    break;
                }
            }
        }

        // Merge our output maps together:
        outputMap.putAll(extraFieldOutputMap);

        // Now translate fields/values to the field names/values that MAF likes:
        final LinkedHashMap<String, String> mafCompliantMap = replaceFuncotationValuesWithMafCompliantValues(outputMap);

        // Remove any fields that are excluded and sanitize any field values.
        return mafCompliantMap.keySet().stream()
            .filter(k -> !excludedOutputFields.contains(k))
            .collect(Collectors.toMap(
                    Function.identity(), k -> FuncotatorUtils.sanitizeFuncotationFieldForMaf(mafCompliantMap.get(k)),
                    (u, v) -> {
                        throw new GATKException.ShouldNeverReachHereException("Found duplicate keys for MAF output");
                        },
                    LinkedHashMap::new));
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * @return The {@link #defaultMap} with currently populated fields.
     */
    @VisibleForTesting
    LinkedHashMap<String, String> getDefaultMap() {
        return defaultMap;
    }

    /**
     * Ensures each value in the given {@code outputMap} is compliant with MAF-expected values.
     * @param outputMap The {@link Map} of output field -> values to be checked for MAF compliance.
     * @return A {@link LinkedHashMap} of output field strings -> values to be written to the MAF file.
     */
    @VisibleForTesting
    LinkedHashMap<String, String> replaceFuncotationValuesWithMafCompliantValues( final Map<String, Object> outputMap ) {

        final LinkedHashMap<String, String> finalOutMap = new LinkedHashMap<>(outputMap.size());

        // Massage individual Key/Value pairs:
        for ( final String key : outputMap.keySet() ) {
            finalOutMap.put(key, mafTransform(key, outputMap.get(key).toString(), referenceVersion) );
        }

        // Massage the OtherTranscripts field:
        if ( finalOutMap.containsKey(MafOutputRendererConstants.FieldName_Other_Transcripts) ) {
            finalOutMap.put(
                    MafOutputRendererConstants.FieldName_Other_Transcripts, finalOutMap.get(MafOutputRendererConstants.FieldName_Other_Transcripts)
                            .replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)
            );
        }

        // Fix the alleles in the case of INDELS:
        adjustIndelAlleleInformationForMafOutput(finalOutMap);

        return finalOutMap;
    }

    /**
     * Checks the given {@code outputMap} has MAF-correct values for information relating to INDEL alleles.
     * NOTE: The output map is modified in place.
     * @param outputMap The {@link Map} of output field -> values to be checked for MAF compliance.
     */
    @VisibleForTesting
    void adjustIndelAlleleInformationForMafOutput(final LinkedHashMap<String, String> outputMap) {
        // Massage the start/end/alleles in the case of INDELs
        // (Because MAF has different conventions from VCF for start/end positions of INDELs)
        if ( outputMap.containsKey(MafOutputRendererConstants.FieldName_Variant_Type) &&
            (outputMap.get(MafOutputRendererConstants.FieldName_Variant_Type).equals(MafOutputRendererConstants.FieldValue_Variant_Type_Insertion) ||
             outputMap.get(MafOutputRendererConstants.FieldName_Variant_Type).equals(MafOutputRendererConstants.FieldValue_Variant_Type_Deletion)) ) {

            final int refAlleleLength = outputMap.get(MafOutputRendererConstants.FieldName_Reference_Allele).length();
            final int altAlleleLength = outputMap.get(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2).length();

            // TODO: port these changes to GencodeFuncotationFactory (issue: https://github.com/broadinstitute/gatk/issues/4378)
            // Check to see if it's an insertion:
            if ( refAlleleLength < altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the ALT_allele where N = length(ref_allele)
                //    Replace the ref_allele with "-"
                //    Replace the Tumor_Seq_Allele1 with "-"
                //    Set the End_Position to be Start_Position + 1 (All Insertions should have length 1 to represent the bases between which the insertion occurs).
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, outputMap.get(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2).substring(refAlleleLength));
                outputMap.put(MafOutputRendererConstants.FieldName_Reference_Allele,  MafOutputRendererConstants.EmptyAllele);
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, MafOutputRendererConstants.EmptyAllele);
                outputMap.put(MafOutputRendererConstants.FieldName_End_Position, String.valueOf(Integer.valueOf(outputMap.get(MafOutputRendererConstants.FieldName_Start_Position)) + 1));
            }
            // Check to see if it's a deletion:
            else if ( refAlleleLength > altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the REF_allele where N = length(alt_allele)
                //    Remove the first N bases from the Tumor_Seq_Allele1
                //    Replace the alt_allele with "-"
                //    Increment the Start_Position by 1 (start position should be inclusive of the first base deleted)
                //    Set the End_Position to Start_Position + M - 2, where M = length(ref_allele) (end position should be inclusive of the last base deleted, explanation below)
                outputMap.put(MafOutputRendererConstants.FieldName_Reference_Allele,  outputMap.get(MafOutputRendererConstants.FieldName_Reference_Allele).substring(altAlleleLength));
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, outputMap.get(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1).substring(altAlleleLength));
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, MafOutputRendererConstants.EmptyAllele);
                outputMap.put(MafOutputRendererConstants.FieldName_Start_Position, String.valueOf(Integer.valueOf(outputMap.get(MafOutputRendererConstants.FieldName_Start_Position)) + 1));

                // Use the new start position we just calculated for this end position.
                // Then subtract 2 from ref allele length:
                //    -1 for the removed first base
                //    -1 for the inclusive nature of positions
                outputMap.put(MafOutputRendererConstants.FieldName_End_Position, String.valueOf(Integer.valueOf(outputMap.get(MafOutputRendererConstants.FieldName_Start_Position)) + refAlleleLength - 2));
            }
        }
    }

    /**
     * Transforms a given {@code value} to the equivalent MAF-valid value based on the given {@code key}.
     * @param key The {@code key} off of which to base the transformation.  This key is the final (transformed) key for output (i.e. the column name in the MAF file).
     * @param value The {@code value} to transform into a MAF-valid value.
     * @param referenceVersion The version of the reference used to create these annotations.
     * @return The MAF-valid equivalent of the given {@code value}.
     */
    public static String mafTransform(final String key, final String value, final String referenceVersion) {

        switch (key) {
            case MafOutputRendererConstants.FieldName_Variant_Classification:
                if ( MafOutputRendererConstants.VariantClassificationMap.containsKey(value)) {
                    return MafOutputRendererConstants.VariantClassificationMap.get(value);
                }
                break;
                
            case MafOutputRendererConstants.FieldName_Chromosome:
                if ( value.equals(MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito) ) {
                    return MafOutputRendererConstants.FieldValue_Chromosome_Mito;
                }
                else if ( value.toLowerCase().startsWith("chr") && (referenceVersion.toLowerCase().equals("hg19") || referenceVersion.toLowerCase().equals("b37"))) {
                    final String trimVal = value.substring(3);
                    if ( HG_19_CHR_SET.contains(trimVal)) {
                        return trimVal;
                    }
                }
                break;
            case MafOutputRendererConstants.FieldName_Other_Transcripts:
                // Use apache commons string utils because it's much, much faster to do this replacement:
                return StringUtils.replaceEachRepeatedly(value, ORDERED_GENCODE_VARIANT_CLASSIFICATIONS.toArray(new String[]{}), ORDERED_MAF_VARIANT_CLASSIFICATIONS.toArray(new String[]{}));
        }

        return value;
     }

    /**
     * Transforms a given {@code value} from a MAF-valid value to a general-purpose value based on the given {@code key}.
     * @param key The {@code key} off of which to base the transformation.  This key is the transformed key for output (i.e. the column name in the MAF file).
     * @param value The {@code value} to transform from a MAF-valid value into a general-purpose value.
     * @param referenceVersion The version of the reference used to create these annotations.
     * @return The general-purpose equivalent of the given MAF-valid {@code value}.
     */
    public static String mafTransformInvert(final String key, final String value, final String referenceVersion ) {
        switch (key) {
            case MafOutputRendererConstants.FieldName_Variant_Classification:
                if ( MafOutputRendererConstants.VariantClassificationMapInverse.containsKey(value)) {
                    return MafOutputRendererConstants.VariantClassificationMapInverse.get(value);
                }
                break;

            case MafOutputRendererConstants.FieldName_Chromosome:
                if ( value.equals(MafOutputRendererConstants.FieldValue_Chromosome_Mito) ) {
                    return MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito;
                }
                else if (referenceVersion.toLowerCase().equals("hg19") || referenceVersion.toLowerCase().equals("b37")) {
                    if ( value.length() <= 2 ) {
                        if ( HG_19_CHR_SET.contains(value) ) {
                            return "chr" + value;
                        }
                    }
                }
                break;
            case MafOutputRendererConstants.FieldName_Other_Transcripts:

                // Use apache commons string utils because it's much, much faster to do this replacement:
                // But we have to do the LINCRNA -> RNA conversion separately, so exclude those indices:
                String replacement = StringUtils.replaceEachRepeatedly(
                        value,
                        ORDERED_MAF_VARIANT_CLASSIFICATIONS.subList(0, ORDERED_MAF_VARIANT_CLASSIFICATIONS.size()-1).toArray(new String[]{}),
                        ORDERED_GENCODE_VARIANT_CLASSIFICATIONS.subList(0, ORDERED_MAF_VARIANT_CLASSIFICATIONS.size()-1).toArray(new String[]{})
                );

                // Handle the RNA/LINCRNA case specially:
                final String needle = MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString());
                int indx = replacement.indexOf(needle);
                while ( indx != -1 ) {
                    // make sure the previous letters are not LINC:
                    if ( (indx <= 3 ) || ((indx > 3) && (!replacement.substring(indx - 4, indx).equals("LINC")))) {
                        replacement = replacement.substring(0, indx) + GencodeFuncotation.VariantClassification.LINCRNA.toString() + replacement.substring(indx + needle.length());
                        indx += needle.length();
                    }
                    else {
                        ++indx;
                    }

                    indx = replacement.indexOf(needle, indx);
                }

                return replacement;
        }

        return value;
    }

    /**
     * Set the field in the given {@code outputMap} specified by the given {@code key} to the given {@code value}.
     * Checks the given {@code key} to see if the given {@code value} must be transformed to be a valid MAF field.
     * Will only set this field if the given {@code value} is not {@code null}.
     * @param outputMap {@link Map} of {@link String} to {@link Object} to hold output annotations.  Must not be {@code null}.
     * @param key {@link String} key to add to the output map.  Must not be {@code null}.
     * @param value {@link Object} value for the output map.
     */
    private void setField(final Map<String, Object> outputMap, final String key, final Object value) {
        Utils.nonNull(outputMap);
        Utils.nonNull(key);

        if ( value != null ) {
            outputMap.put(key, value);
        }
        else {
            outputMap.put(key, FuncotatorConstants.UNKNOWN_VALUE_STRING);
        }
    }

    /**
     * Write the given line to the {@link #writer} and append a newline.
     * @param line The {@link String} to write as a line to the {@link #writer}.
     */
    private void writeLine(final String line) throws IOException {
        writeString(line);
        writeString(System.lineSeparator());
    }

    /**
     * Write the given {@link String} to the {@link #writer}.
     * @param s The {@link String} to write to the {@link #writer}.
     */
    private void writeString(final String s) throws IOException {
        writer.write(s);
    }

    /**
     * Write the header to the output file.
     * @param outputFields Ordered list of the header columns.  These will be written as presented.
     */
    protected void writeHeader(final List<String> outputFields) {
        try {
            // Write out version:
            writeLine(MafOutputRendererConstants.COMMENT_STRING + "version " + VERSION);
            writeLine(MafOutputRendererConstants.COMMENT_STRING + MafOutputRendererConstants.COMMENT_STRING);

            // Write previous header info:
            for (final VCFHeaderLine line : inputFileHeader.getMetaDataInInputOrder()) {
                writer.write(MafOutputRendererConstants.COMMENT_STRING);
                writer.write(MafOutputRendererConstants.COMMENT_STRING);
                writer.write(" ");
                writeLine(line.toString());
            }

            // Write any default tool header lines:
            for (final String line : toolHeaderLines) {
                writer.write(MafOutputRendererConstants.COMMENT_STRING);
                writer.write(MafOutputRendererConstants.COMMENT_STRING);
                writer.write(" ");
                writeLine(line);
            }

            // Write tool name and the data sources with versions:
            writer.write(MafOutputRendererConstants.COMMENT_STRING);
            writer.write(MafOutputRendererConstants.COMMENT_STRING);
            writer.write(" ");
            writer.write(" Funcotator ");
            writer.write(toolVersion);
            writer.write(" | Date ");
            writer.write(new SimpleDateFormat("yyyymmdd'T'hhmmss").format(new Date()));
            writer.write(" | ");
            writer.write(getDataSourceInfoString());
            writeLine("");

            // Write the column headers for our output set and our manual annotations:
            writer.write(outputFields.stream().collect(Collectors.joining(MafOutputRendererConstants.FIELD_DELIMITER)));
            if (manualAnnotations.size() > 0) {
                writeLine(MafOutputRendererConstants.FIELD_DELIMITER + manualAnnotations.keySet()
                        .stream()
                        .collect(Collectors.joining(MafOutputRendererConstants.FIELD_DELIMITER)));
            }
            else {
                writer.write("\n");
            }
            // Make sure we keep track of the fact that we've now written the header:
            hasWrittenHeader = true;
        } catch (IOException e){
            throw new UserException.CouldNotCreateOutputFile("Exception while writing maf header, caused by " + e.getMessage(), e);
        }
    }

    /** No field ordering is preserved in the output of this method.
     * @return Immutable copy of the reverse output field map used by this MafOutputRenderer.
     * In other words, return a mapping for the input field name to the MAF column this will produce. Never {@code null}
     */
    public ImmutableMap<String, Set<String>> getReverseOutputFieldNameMap() {
        return  ImmutableMap.copyOf(Utils.getReverseValueToListMap(outputFieldNameMap));
    }

    /**
     * Initializes {@link MafOutputRenderer#defaultMap} with the default keys for the columns in a MAF file.
     */
    protected void initializeDefaultMapWithKeys() {
        
        // Baseline required fields:
        defaultMap.put(MafOutputRendererConstants.FieldName_Hugo_Symbol                   ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Entrez_Gene_Id                ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Center                        ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_NCBI_Build                    ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Chromosome                    ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Start_Position                ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_End_Position                  ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Strand                        ,               MafOutputRendererConstants.FieldValue_Strand );
        defaultMap.put(MafOutputRendererConstants.FieldName_Variant_Classification        ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Variant_Type                  ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Reference_Allele              ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1             ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2             ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_dbSNP_RS                      ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_dbSNP_Val_Status              ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Sample_Barcode          ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Matched_Norm_Sample_Barcode   ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele1        ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele2        ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Validation_Allele1      ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Validation_Allele2      ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele1 ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele2 ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Verification_Status           ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Validation_Status             ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Mutation_Status               ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Sequencing_Phase              ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Sequence_Source               ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Validation_Method             ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Score                         ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_BAM_File                      ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Sequencer                     ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Sample_UUID             ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Matched_Norm_Sample_UUID      ,               FuncotatorConstants.UNKNOWN_VALUE_STRING );

        // Required "optional" fields:
        defaultMap.put(MafOutputRendererConstants.FieldName_Genome_Change                          ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Annotation_Transcript                  ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Transcript_Strand                      ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Transcript_Exon                        ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Transcript_Position                    ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_cDNA_Change                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Codon_Change                           ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Protein_Change                         ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Other_Transcripts                      ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Refseq_mRNA_Id                         ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Refseq_prot_Id                         ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_SwissProt_acc_Id                       ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_SwissProt_entry_Id                     ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Description                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_AApos                          ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Region                         ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Site                           ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Natural_Variations             ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Experimental_Info              ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_GO_Biological_Process                  ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_GO_Cellular_Component                  ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_GO_Molecular_Function                  ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_overlapping_mutations           ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_fusion_genes                    ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_tissue_types_affected           ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_total_alterations_in_gene       ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumorscape_Amplification_Peaks         ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumorscape_Deletion_Peaks              ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_TCGAscape_Amplification_Peaks          ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_TCGAscape_Deletion_Peaks               ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_DrugBank                               ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_ref_context                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_gc_content                             ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_overlapping_mutations     ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_total_mutations_in_gene   ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Mutation_Type                      ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Translocation_Partner              ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Somatic                ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Germline               ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Other_Diseases                     ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_DNARepairGenes_Activity_linked_to_OMIM ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_FamilialCancerDatabase_Syndromes       ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_MUTSIG_Published_Results               ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_OREGANNO_ID                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_OREGANNO_Values                        ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_tumor_f                                ,      FuncotatorConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_t_alt_count                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_t_ref_count                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_n_alt_count                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_n_ref_count                            ,      FuncotatorConstants.UNKNOWN_VALUE_STRING );

    }

    /**
     * Initializes the output map
     */
    private void initializeOutputFieldNameMap() {

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Hugo_Symbol                            , MafOutputRendererConstants.OutputFieldNameMap_Hugo_Symbol );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Entrez_Gene_Id                         , MafOutputRendererConstants.OutputFieldNameMap_Entrez_Gene_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Center                                 , MafOutputRendererConstants.OutputFieldNameMap_Center );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_NCBI_Build                             , MafOutputRendererConstants.OutputFieldNameMap_NCBI_Build );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Chromosome                             , MafOutputRendererConstants.OutputFieldNameMap_Chromosome );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Start_Position                         , MafOutputRendererConstants.OutputFieldNameMap_Start_Position );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_End_Position                           , MafOutputRendererConstants.OutputFieldNameMap_End_Position );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Strand                                 , MafOutputRendererConstants.OutputFieldNameMap_Strand );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Variant_Classification                 , MafOutputRendererConstants.OutputFieldNameMap_Variant_Classification );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Variant_Type                           , MafOutputRendererConstants.OutputFieldNameMap_Variant_Type );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Reference_Allele                       , MafOutputRendererConstants.OutputFieldNameMap_Reference_Allele );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1                      , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Seq_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2                      , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Seq_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_dbSNP_RS                               , MafOutputRendererConstants.OutputFieldNameMap_dbSNP_RS );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_dbSNP_Val_Status                       , MafOutputRendererConstants.OutputFieldNameMap_dbSNP_Val_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Sample_Barcode                   , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Sample_Barcode );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Matched_Norm_Sample_Barcode            , MafOutputRendererConstants.OutputFieldNameMap_Matched_Norm_Sample_Barcode );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele1                 , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Seq_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele2                 , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Seq_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Validation_Allele1               , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Validation_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Validation_Allele2               , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Validation_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele1          , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Validation_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele2          , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Validation_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Verification_Status                    , MafOutputRendererConstants.OutputFieldNameMap_Verification_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Validation_Status                      , MafOutputRendererConstants.OutputFieldNameMap_Validation_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Mutation_Status                        , MafOutputRendererConstants.OutputFieldNameMap_Mutation_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Sequencing_Phase                       , MafOutputRendererConstants.OutputFieldNameMap_Sequencing_Phase );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Sequence_Source                        , MafOutputRendererConstants.OutputFieldNameMap_Sequence_Source );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Validation_Method                      , MafOutputRendererConstants.OutputFieldNameMap_Validation_Method );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Score                                  , MafOutputRendererConstants.OutputFieldNameMap_Score );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_BAM_File                               , MafOutputRendererConstants.OutputFieldNameMap_BAM_File );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Sequencer                              , MafOutputRendererConstants.OutputFieldNameMap_Sequencer );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Sample_UUID                      , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Sample_UUID );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Matched_Norm_Sample_UUID               , MafOutputRendererConstants.OutputFieldNameMap_Matched_Norm_Sample_UUID );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Genome_Change                          , MafOutputRendererConstants.OutputFieldNameMap_Genome_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Annotation_Transcript                  , MafOutputRendererConstants.OutputFieldNameMap_Annotation_Transcript );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Transcript_Strand                      , MafOutputRendererConstants.OutputFieldNameMap_Transcript_Strand );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Transcript_Exon                        , MafOutputRendererConstants.OutputFieldNameMap_Transcript_Exon );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Transcript_Position                    , MafOutputRendererConstants.OutputFieldNameMap_Transcript_Position );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_cDNA_Change                            , MafOutputRendererConstants.OutputFieldNameMap_cDNA_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Codon_Change                           , MafOutputRendererConstants.OutputFieldNameMap_Codon_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Protein_Change                         , MafOutputRendererConstants.OutputFieldNameMap_Protein_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Other_Transcripts                      , MafOutputRendererConstants.OutputFieldNameMap_Other_Transcripts );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Refseq_mRNA_Id                         , MafOutputRendererConstants.OutputFieldNameMap_Refseq_mRNA_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Refseq_prot_Id                         , MafOutputRendererConstants.OutputFieldNameMap_Refseq_prot_Id );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_SwissProt_acc_Id                       , MafOutputRendererConstants.OutputFieldNameMap_SwissProt_acc_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_SwissProt_entry_Id                     , MafOutputRendererConstants.OutputFieldNameMap_SwissProt_entry_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Description                            , MafOutputRendererConstants.OutputFieldNameMap_Description );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_AApos                          , MafOutputRendererConstants.OutputFieldNameMap_UniProt_AApos );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Region                         , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Region );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Site                           , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Site );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Natural_Variations             , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Natural_Variations );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Experimental_Info              , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Experimental_Info );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_GO_Biological_Process                  , MafOutputRendererConstants.OutputFieldNameMap_GO_Biological_Process );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_GO_Cellular_Component                  , MafOutputRendererConstants.OutputFieldNameMap_GO_Cellular_Component );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_GO_Molecular_Function                  , MafOutputRendererConstants.OutputFieldNameMap_GO_Molecular_Function );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_overlapping_mutations           , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_overlapping_mutations );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_fusion_genes                    , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_fusion_genes );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_tissue_types_affected           , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_tissue_types_affected );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_total_alterations_in_gene       , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_total_alterations_in_gene );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumorscape_Amplification_Peaks         , MafOutputRendererConstants.OutputFieldNameMap_Tumorscape_Amplification_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumorscape_Deletion_Peaks              , MafOutputRendererConstants.OutputFieldNameMap_Tumorscape_Deletion_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_TCGAscape_Amplification_Peaks          , MafOutputRendererConstants.OutputFieldNameMap_TCGAscape_Amplification_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_TCGAscape_Deletion_Peaks               , MafOutputRendererConstants.OutputFieldNameMap_TCGAscape_Deletion_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_DrugBank                               , MafOutputRendererConstants.OutputFieldNameMap_DrugBank );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_ref_context                            , MafOutputRendererConstants.OutputFieldNameMap_ref_context );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_gc_content                             , MafOutputRendererConstants.OutputFieldNameMap_gc_content );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_overlapping_mutations     , MafOutputRendererConstants.OutputFieldNameMap_CCLE_ONCOMAP_overlapping_mutations );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_total_mutations_in_gene   , MafOutputRendererConstants.OutputFieldNameMap_CCLE_ONCOMAP_total_mutations_in_gene );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Mutation_Type                      , MafOutputRendererConstants.OutputFieldNameMap_CGC_Mutation_Type );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Translocation_Partner              , MafOutputRendererConstants.OutputFieldNameMap_CGC_Translocation_Partner );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Somatic                , MafOutputRendererConstants.OutputFieldNameMap_CGC_Tumor_Types_Somatic );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Germline               , MafOutputRendererConstants.OutputFieldNameMap_CGC_Tumor_Types_Germline );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Other_Diseases                     , MafOutputRendererConstants.OutputFieldNameMap_CGC_Other_Diseases );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_DNARepairGenes_Activity_linked_to_OMIM , MafOutputRendererConstants.OutputFieldNameMap_DNARepairGenes_Activity_linked_to_OMIM );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_FamilialCancerDatabase_Syndromes       , MafOutputRendererConstants.OutputFieldNameMap_FamilialCancerDatabase_Syndromes );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_MUTSIG_Published_Results               , MafOutputRendererConstants.OutputFieldNameMap_MUTSIG_Published_Results );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_OREGANNO_ID                            , MafOutputRendererConstants.OutputFieldNameMap_OREGANNO_ID );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_OREGANNO_Values                        , MafOutputRendererConstants.OutputFieldNameMap_OREGANNO_Values );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_tumor_f                                , MafOutputRendererConstants.OutputFieldNameMap_tumor_f );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_t_alt_count                            , MafOutputRendererConstants.OutputFieldNameMap_t_alt_count );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_t_ref_count                            , MafOutputRendererConstants.OutputFieldNameMap_t_ref_count );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_n_alt_count                            , MafOutputRendererConstants.OutputFieldNameMap_n_alt_count );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_n_ref_count                            , MafOutputRendererConstants.OutputFieldNameMap_n_ref_count );
    }

    //==================================================================================================================
    // Helper Data Types:

    //------------------------------------------------------------------------------------------------------------------
    // Columns:
    
// Required:
// 1	Hugo_Symbol
// 2	Entrez_Gene_Id
// 3	Center
// 4	NCBI_Build
// 5	Chromosome
// 6	Start_Position
// 7	End_Position
// 8	Strand
// 9	Variant_Classification
//10	Variant_Type
//11	Reference_Allele
//12	Tumor_Seq_Allele1
//13	Tumor_Seq_Allele2
//14	dbSNP_RS
//15	dbSNP_Val_Status
//16	Tumor_Sample_Barcode
//17	Matched_Norm_Sample_Barcode
//18	Match_Norm_Seq_Allele1
//19	Match_Norm_Seq_Allele2
//20	Tumor_Validation_Allele1
//21	Tumor_Validation_Allele2
//22	Match_Norm_Validation_Allele1
//23	Match_Norm_Validation_Allele2
//24	Verification_Status
//25	Validation_Status
//26	Mutation_Status
//27	Sequencing_Phase
//28	Sequence_Source
//29	Validation_Method
//30	Score
//31	BAM_File
//32	Sequencer
//33    Tumor_Sample_UUID
//34    Matched_Norm_Sample_UUID

// "Optional" Required:
//35    Genome_Change
//36    Annotation_Transcript
//37    Transcript_Strand
//38    Transcript_Exon
//39    Transcript_Position
//40    cDNA_Change
//41    Codon_Change
//42    Protein_Change
//43    Other_Transcripts
//44    Refseq_mRNA_Id
//45    Refseq_prot_Id
//46    SwissProt_acc_Id
//47    SwissProt_entry_Id
//48    Description
//49    UniProt_AApos
//50    UniProt_Region
//51    UniProt_Site
//52    UniProt_Natural_Variations
//53    UniProt_Experimental_Info
//54    GO_Biological_Process
//55    GO_Cellular_Component
//56    GO_Molecular_Function
//57    COSMIC_overlap
//58    ping_mutations
//59    COSMIC_fusion_genes
//60    COSMIC_tissue_types_affected
//61    COSMIC_total_alterations_in_gene
//62    Tumorscape_Amplification_Peaks
//63    Tumorscape_Deletion_Peaks
//64    TCGAscape_Amplification_Peaks
//65    TCGAscape_Deletion_Peaks
//66    DrugBank
//67    ref_context
//68    gc_content
//69    CCLE_ONCOMAP_overlapping_mutations
//70    CCLE_ONCOMAP_total_mutations_in_gene
//71    CGC_Mutation_Type
//72    CGC_Translocation_Partner
//73    CGC_Tumor_Types_Somatic
//74    CGC_Tumor_Types_Germline
//75    CGC_Other_Diseases
//76    DNARepairGenes_Role
//77    FamilialCancerDatabase_Syndromes
//78    MUTSIG_Published_Results
//79    OREGANNO_ID
//80    OREGANNO_Values
//81    tumor_f

}
