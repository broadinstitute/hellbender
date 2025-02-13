package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.*;
import java.util.function.Function;

import static org.broadinstitute.hellbender.utils.Utils.split;

/**
 * Extract fields from a VCF file to a tab-delimited table
 *
 * <p>
 *     This tool extracts specified fields for each variant in a VCF file to a tab-delimited table, which may be easier
 *     to work with than a VCF. By default, the tool only extracts PASS or . (unfiltered) variants in the VCF file. Filtered variants may be
 *     included in the output by adding the --show-filtered flag. The tool can extract both INFO (i.e. site-level) fields and
 *     FORMAT (i.e. sample-level) fields. If the tool is run without specifying any fields, it defaults to include all fields
 *     declared in the VCF header.
 * </p>
 *
 * <h4>INFO/site-level fields</h4>
 * <p>
 *     Use the `-F` argument to extract INFO fields; each field will occupy a single column in the output
 * file.  The field can be any standard VCF column (e.g. CHROM, ID, QUAL) or any annotation name in the INFO field (e.g. AC, AF).
 * The tool also supports the following additional fields:
 * <ul>
 *     <li> EVENTLENGTH (length of the event) </li>
 *     <li> TRANSITION (1 for a bi-allelic transition (SNP), 0 for bi-allelic transversion (SNP), -1 for INDELs and multi-allelics) </li>
 *     <li> HET (count of het genotypes) </li>
 *     <li> HOM-REF (count of homozygous reference genotypes) </li>
 *     <li> HOM-VAR (count of homozygous variant genotypes) </li>
 *     <li> NO-CALL (count of no-call genotypes) </li>
 *     <li> TYPE (type of variant, possible values are NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, and MIXED </li>
 *     <li> VAR (count of non-reference genotypes) </li>
 *     <li> NSAMPLES (number of samples) </li>
 *     <li> NCALLED (number of called samples) </li>
 *     <li> MULTI-ALLELIC (is this variant multi-allelic? true/false) </li>
 * </ul>
 * <p>
 *     Use the `-ASF` argument to extract allele-specific/per allele INFO fields and split them appropriately when
 *     splitting multi-allelic variants.
 * </p>
 *
 * <h4>FORMAT/sample-level fields</h4>
 * <p>
 *     Use the `-GF` argument to extract FORMAT/sample-level fields. The tool will create a new column per sample
 *     with the name "SAMPLE_NAME.FORMAT_FIELD_NAME" e.g. NA12877.GQ, NA12878.GQ.
 * </p>
 * <p>
 *     Use the `-ASGF` argument to extract allele-specific/per allele FORMAT fields and split them appropriately
 *     when splitting multi-allelic variants.  If AD is specified as an allele-specific genotype field the ref and alt
 *     counts will be given for each alt.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li> A VCF file to convert to a table </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 *     A tab-delimited file containing the values of the requested fields in the VCF file.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 *     gatk VariantsToTable \
 *     -V input.vcf \
 *     -F CHROM -F POS -F TYPE -GF AD \
 *     -O output.table
 * </pre>
 * <p>would produce a file that looks like:</p>
 * <pre>
 *     CHROM  POS        TYPE   HSCX1010N.AD  HSCX1010T.AD
 *     1      31782997   SNP    77,0          53,4
 *     1      40125052   SNP    97,0          92,7
 *     1      65068538   SNP    49,0          35,4
 *     1      111146235  SNP    69,1          77,4
 * </pre>
 * <pre>
 *     gatk VariantsToTable \
 *     -V input.vcf \
 *     -O output.table
 * </pre>
 * <p>would produce a file that includes all fields declared in the VCF header.</p>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li> It is common for certain annotations to be absent for some variants. By default, this tool will emit an NA for a missing annotation. If you prefer that the tool fail upon encountering a missing annotation, use the --error-if-missing-data flag. </li>
 *     <li> If multiple samples are present in the VCF, the genotype fields will be ordered alphabetically by sample name. </li>
 *     <li> Filtered sites are ignored by default. To include them in the output, use the --show-filtered flag. </li>
 *     <li> Allele-specific filtering is not yet supported.  For PASS sites, all alleles will be given, regardless of their AS_FilterStatus.</li>
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier to work with than a VCF",
        oneLineSummary = "Extract fields from a VCF file to a tab-delimited table",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public final class VariantsToTable extends VariantWalker {
    public final static String SPLIT_MULTI_ALLELIC_LONG_NAME = "split-multi-allelic";
    public final static String SPLIT_MULTI_ALLELIC_SHORT_NAME = "SMA";
    public static final String NUMERIC_GT_FULLNAME = "use-numeric-gt";

    static final Logger logger = LogManager.getLogger(VariantsToTable.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which the tab-delimited table is written")
    private String out = null;

    /**
     * Any standard VCF column (CHROM, ID, QUAL) or any annotation name in the INFO field (e.g., -F AC) to include in
     * the output table. To capture FORMAT field values, see the -GF argument. This argument accepts any number
     * of inputs e.g. -F CHROM -F POS
     */
    @Argument(fullName="fields",
            shortName="F",
            doc="The name of a standard VCF field or an INFO field to include in the output table", optional=true)
    protected List<String> fieldsToTake = new ArrayList<>();

    /**
     * Any annotation name in the FORMAT field (e.g., GQ, PL) to include in the output table.
     * This argument accepts any number of inputs e.g. -GF GQ -GF PL
     */
    @Argument(fullName="genotype-fields",
            shortName="GF",
            doc="The name of a genotype field to include in the output table", optional=true)
    private List<String> genotypeFieldsToTake = new ArrayList<>();

    @Argument(shortName="ASF", doc="The name of an allele-specific INFO field to be split if present", optional=true)
    private List<String> asFieldsToTake = new ArrayList<>();

    @Argument(shortName="ASGF", doc="The name of an allele-specific FORMAT field to be split if present", optional=true)
    private List<String> asGenotypeFieldsToTake = new ArrayList<>();

    /**
     * By default this tool only emits values for records where the FILTER field is either PASS or . (unfiltered).
     * Turn on this flag to emit values regardless of the value of the FILTER field.
     */
    @Advanced
    @Argument(fullName="show-filtered",
            shortName="raw",
            doc="Include filtered records in the output", optional=true)
    private boolean showFiltered = false;

    /**
     * By default, a variant record with multiple ALT alleles will be summarized in one line, with per alt-allele fields
     * (e.g. allele depth) separated by commas. This may cause difficulty when the table is loaded by an R script, for example.
     * Use this flag to write multi-allelic records on separate lines of output. Fields that are not allele-specific will be duplicated.
     */
    @Argument(fullName=SPLIT_MULTI_ALLELIC_LONG_NAME,
            shortName=SPLIT_MULTI_ALLELIC_SHORT_NAME,
            doc="Split multi-allelic records into multiple lines", optional=true)
    protected boolean splitMultiAllelic = false;

    /**
     * Use this flag to emit each field within a variant on a separate line. The resulting table will have
     * four columns: RecordID, Sample, Variable, and Value. Variable refers to the field name, Value to the value of the
     * field. The tool prints "site" under Sample column for an INFO or standard field.
     *
     * Example: -F CHROM -GF AD will print the following table
     * RecordID  Sample   Variable  Value
     * 1         site     CHROM     20
     * 1         NA12878  AD        36,28
     * 2         site     CHROM     20
     * 2         NA12878  AD        26,27
     * 3         site     CHROM     20
     */
    @Advanced
    @Argument(fullName="moltenize",
            shortName="moltenize",
            doc="Produce molten output", optional=true)
    private boolean moltenizeOutput = false;

    /**
     * By default, this tool will write NA for missing data.
     * Turn on this flag, and the tool will throw an error and exit if it encounters missing data.
     */
    @Advanced
    @Argument(fullName="error-if-missing-data",
            shortName="EMD",
            doc="Fail on missing data", optional=true)
    public boolean errorIfMissingData = false;

    @Argument(fullName = NUMERIC_GT_FULLNAME,
              doc = "write the GT field the way it appears in a VCF ( ex. 0/1 instead of A/T )",
              optional = true)
    public boolean useNumericGT = false;

    private static final String MISSING_DATA = "NA";

    private SortedSet<String> samples;
    private long nRecords = 0L;
    private PrintStream outputStream = null;
    private VCFHeader inputHeader;

    @Override
    public void onTraversalStart() {
        inputHeader = getHeaderForVariants();
        outputStream = createPrintStream();

        // if no fields specified, default to include all fields listed in header into table
        if(fieldsToTake.isEmpty() && genotypeFieldsToTake.isEmpty() && asFieldsToTake.isEmpty() && asGenotypeFieldsToTake.isEmpty()){
            logger.warn("No fields were specified. All fields declared in the VCF header will be included in the output table.");

            // add all mandatory VCF fields (except INFO)
            for(VCFHeader.HEADER_FIELDS headerField : VCFHeader.HEADER_FIELDS.values()){
                if(!headerField.name().equals(VCFHeader.HEADER_FIELDS.INFO.name())) {
                    fieldsToTake.add(headerField.name());
                }
            }

            // add all INFO fields present in VCF header
            for (final VCFInfoHeaderLine infoLine : inputHeader.getInfoHeaderLines()) {
                fieldsToTake.add(infoLine.getID());
            }

            // add all FORMAT fields present in VCF header
            for (final VCFFormatHeaderLine formatLine : inputHeader.getFormatHeaderLines()) {
                // ensure GT field listed as first FORMAT field
                if(formatLine.getID().equals(VCFConstants.GENOTYPE_KEY)) {
                    genotypeFieldsToTake.add(0, formatLine.getID());
                }
                else {
                    genotypeFieldsToTake.add(formatLine.getID());
                }
            }
        }

        // if fields specified, but none are genotype fields, set samples to empty
        if (genotypeFieldsToTake.isEmpty() && asGenotypeFieldsToTake.isEmpty()) {
            samples = Collections.emptySortedSet();
        }
        else {
            final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
            samples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

            // if there are no samples, we don't have to worry about any genotype fields
            if (samples.isEmpty()) {
                genotypeFieldsToTake.clear();
                asGenotypeFieldsToTake.clear();
                logger.warn("There are no samples - the genotype fields will be ignored");
                if (fieldsToTake.isEmpty() && asFieldsToTake.isEmpty()){
                    throw new UserException("There are no samples and no fields - no output will be produced");
                }
            }
        }

        if (asGenotypeFieldsToTake.isEmpty() && asFieldsToTake.isEmpty() && !splitMultiAllelic) {
            logger.warn("Allele-specific fields will only be split if splitting multi-allelic variants is specified (`--" + SPLIT_MULTI_ALLELIC_LONG_NAME + "` or `-" + SPLIT_MULTI_ALLELIC_SHORT_NAME + "`");
        }

        // print out the header
        if ( moltenizeOutput ) {
            outputStream.println("RecordID\tSample\tVariable\tValue");
        } else {
            final List<String> fields = new ArrayList<>();

            fields.addAll(fieldsToTake);
            fields.addAll(asFieldsToTake);
            fields.addAll(createGenotypeFields());
            final String header = Utils.join("\t", fields);
            outputStream.println(header);
        }
    }

    private PrintStream createPrintStream() {
        try {
            return out != null ? new PrintStream(out) : System.out;
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(out, e);
        }
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if ( showFiltered || vc.isNotFiltered() ) {
            nRecords++;
            final List<List<String>> records = extractFields(vc);
            if (moltenizeOutput){
                records.forEach(record -> emitMoltenizedOutput(record));
            } else {
                records.forEach(record -> outputStream.println(Utils.join("\t", record)));
            }
        }
    }

    private static boolean isWildCard(final String s) {
        return s.endsWith("*");
    }

    private List<String> createGenotypeFields() {
        final List<String> allGenotypeFieldsToTake = new ArrayList<>(genotypeFieldsToTake);
        allGenotypeFieldsToTake.addAll(asGenotypeFieldsToTake);

        final List<String> genotypeFields = new ArrayList<>();
        for ( final String sample : samples ) {
            for ( final String gf : allGenotypeFieldsToTake ) {
                // spaces in sample names are legal but wreak havoc in R data frames
                final StringBuilder sb = new StringBuilder(sample.replace(" ","_"))
                        .append(".")
                        .append(gf);
                genotypeFields.add(sb.toString());
            }
        }
        return genotypeFields;
    }

    private void emitMoltenizedOutput(final List<String> record) {
        int index = 0;
        for ( final String field : fieldsToTake ) {
            outputStream.println(String.format("%d\tsite\t%s\t%s", nRecords, field, record.get(index++)));
        }
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                outputStream.println(String.format("%d\t%s\t%s\t%s", nRecords, sample.replace(" ","_"), gf, record.get(index++)));
            }
        }
    }

    /**
     * Utility function that returns the list of values for each field in fields from vc.
     *
     * @param vc                the VariantContext whose field values we can to capture
     * @return List of lists of field values
     */
    private List<List<String>> extractFields(final VariantContext vc) {

        final int numRecordsToProduce = splitMultiAllelic ? vc.getAlternateAlleles().size() : 1;
        final List<List<String>> records = new ArrayList<>(numRecordsToProduce);

        for ( int i = 0; i < numRecordsToProduce; i++ ) {
            records.add(new ArrayList<>());
        }

        for ( final String field : fieldsToTake ) {
            if ( splitMultiAllelic && field.equals("ALT") ) { // we need to special case the ALT field when splitting out multi-allelic records
                addFieldValue(splitAltAlleles(vc), records);
            } else if ( getters.containsKey(field) ) {
                addFieldValue(getters.get(field).apply(vc), records);
            } else if ( vc.hasAttribute(field) ) {
                addFieldValue(vc.getAttribute(field, null), records);
            } else if ( isWildCard(field) ) {
                final SortedSet<String> wildVals = new TreeSet<>();
                for ( final Map.Entry<String,Object> elt : vc.getAttributes().entrySet()) {
                    if ( elt.getKey().startsWith(field.substring(0, field.length() - 1)) ) {
                        wildVals.add(elt.getValue().toString());
                    }
                }

                final String val = wildVals.isEmpty() ? MISSING_DATA : Utils.join(",", wildVals);

                addFieldValue(val, records);
            } else {
                handleMissingData(errorIfMissingData, field, records, vc);
            }
        }

        for ( final String field : asFieldsToTake) {
            if (vc.hasAttribute(field)) {
                if (splitMultiAllelic) {
                    addAlleleSpecificFieldValue(Arrays.asList(vc.getAttributeAsString(field, ".").replace("[", "").replace("]", "").split(",")), records, inputHeader.getInfoHeaderLine(field).getCountType());

                } else {
                    addFieldValue(vc.getAttributeAsString(field, ".").replace("[","").replace("]",""), records);
                }
            } else {
                handleMissingData(errorIfMissingData, field, records, vc);
            }
        }

        if ( !genotypeFieldsToTake.isEmpty() || !asGenotypeFieldsToTake.isEmpty() ) {
            addGenotypeFieldsToRecords(vc, records, errorIfMissingData);
        }

        return records;
    }

    private void addGenotypeFieldsToRecords(final VariantContext vc, final List<List<String>> records, final boolean errorIfMissingData) {
        for ( final String sample : samples ) {
            final Genotype genotype = vc.getGenotype(sample);
            for ( final String gf : genotypeFieldsToTake ) {
                if ( genotype != null && genotype.hasAnyAttribute(gf) ) {
                    if (VCFConstants.GENOTYPE_KEY.equals(gf)) {
                        final String gtRepresentation = useNumericGT
                                ? VCFEncoder.encodeGtField(vc, genotype) : genotype.getGenotypeString(true);
                        addFieldValue(gtRepresentation, records);
                    } else {
                        /*
                         * TODO - If gf == "FT" and the GT record is not filtered, Genotype.getAnyAttribute == null. Genotype.hasAnyAttribute should be changed so it
                         * returns false for this condition. Presently, it always returns true. Once this is fixed, then only the "addFieldValue" statement will
                         * remain in the following logic block.
                         */
                        if (genotype.getAnyAttribute(gf) != null) {
                            addFieldValue(genotype.getAnyAttribute(gf), records);
                        } else {
                            handleMissingData(errorIfMissingData, gf, records, vc);
                        }
                    }
                } else {
                    handleMissingData(errorIfMissingData, gf, records, vc);
                }
            }

            for ( final String field : asGenotypeFieldsToTake) {
                if ( genotype != null && genotype.hasAnyAttribute(field) ) {
                    if (splitMultiAllelic) {
                        if (VCFConstants.GENOTYPE_ALLELE_DEPTHS.equals(field)) {
                            final List<String> altDepths = new ArrayList<>();
                            final int[] allDepths = genotype.getAD();
                            for (int i = 1; i < allDepths.length; i++) {
                                altDepths.add(allDepths[0] + "," + allDepths[i]);
                            }
                            addFieldValue(altDepths, records);
                        } else {
                            addAlleleSpecificFieldValue(split(genotype.getExtendedAttribute(field).toString(), ','),
                                    records, inputHeader.getFormatHeaderLine(field).getCountType());
                        }
                    } else {
                        final String value = genotype.getAnyAttribute(field).toString();
                        if (field.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
                            addFieldValue(value.replace("[","").replace("]","").replaceAll("\\s",""),records);
                        } else {
                            addFieldValue(value, records);
                        }
                    }
                } else {
                    handleMissingData(errorIfMissingData, field, records, vc);
                }
            }
        }
    }

    private static void handleMissingData(final boolean errorIfMissingData, final String field, final List<List<String>> records, final VariantContext vc) {
        if (errorIfMissingData) {
            throw new UserException(String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc));
        } else {
            addFieldValue(MISSING_DATA, records);
        }
    }

    private static void addFieldValue(final Object val, final List<List<String>> result) {
        final int numResultRecords = result.size();

        // if we're trying to create a single output record, add it
        if ( numResultRecords == 1 ) {
            result.get(0).add(prettyPrintObject(val));
        }
        // if this field is a list of the proper size, add the appropriate entry to each record
        else if ( (val instanceof List<?> list) && list.size() == numResultRecords ) {
            for ( int i = 0; i < numResultRecords; i++ ) {
                result.get(i).add(list.get(i).toString());
            }
        }
        // otherwise, add the original value to all of the records
        else {
            final String valStr = prettyPrintObject(val);
            for ( final List<String> record : result ) {
                record.add(valStr);
            }
        }
    }

    /**
     * Handle per-allele/allele-specific annotations as described in the header
     * @param val the annotation value
     * @param result the cummulative output
     * @param alleleCount scalar, R-type or A-type values
     */
    private static void addAlleleSpecificFieldValue(final Object val, final List<List<String>> result, final VCFHeaderLineCount alleleCount) {
        if (val instanceof final List<?> list && alleleCount.equals(VCFHeaderLineCount.R)) {
            addFieldValue(new ArrayList<>(list.subList(1, list.size())), result);
        }
        else {
            addFieldValue(val, result);
        }
    }

    private static String prettyPrintObject(final Object val) {
        if ( val == null ) {
            return "";
        }

        if ( val instanceof List<?> list) {
            return prettyPrintObject(list.toArray());
        }

        if ( !val.getClass().isArray() ) {
            return val.toString();
        }

        final int length = Array.getLength(val);
        if ( length == 0 ) {
            return "";
        }

        final StringBuilder sb = new StringBuilder(prettyPrintObject(Array.get(val, 0)));
        for ( int i = 1; i < length; i++ ) {
            sb.append(",");
            sb.append(prettyPrintObject(Array.get(val, i)));
        }
        return sb.toString();
    }

    // ----------------------------------------------------------------------------------------------------
    //
    //  Getting values from VC by name.
    //
    // ----------------------------------------------------------------------------------------------------

    private final Map<String, Function<VariantContext, String>> getters = new LinkedHashMap<>();
    {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", vc -> vc.getContig());
        getters.put("POS", vc -> Integer.toString(vc.getStart()));
        getters.put("REF", vc -> vc.getReference().getDisplayString());
        getters.put("ALT", vc -> {
            final StringBuilder x = new StringBuilder();
            final int n = vc.getAlternateAlleles().size();
            if ( n == 0 ) {
                return ".";
            }

            for ( int i = 0; i < n; i++ ) {
                if ( i != 0 ) {
                    x.append(",");
                }
                x.append(vc.getAlternateAllele(i));
            }
            return x.toString();
        });
        getters.put("EVENTLENGTH", vc -> {
            int maxLength = 0;
            for ( final Allele a : vc.getAlternateAlleles() ) {
                final int length = a.length() - vc.getReference().length();
                if( Math.abs(length) > Math.abs(maxLength) ) { maxLength = length; }
            }
            return Integer.toString(maxLength);
        });
        getters.put("QUAL", vc -> Double.toString(vc.getPhredScaledQual()));
        getters.put("TRANSITION", vc -> {
            if ( vc.isSNP() && vc.isBiallelic() ) {
                return GATKVariantContextUtils.isTransition(vc) ? "1" : "0";
            } else {
                return "-1";
            }
        });
        getters.put("FILTER", vc -> vc.isNotFiltered() ? "PASS" : Utils.join(",", vc.getFilters()));
        getters.put("ID", vc -> vc.getID());
        getters.put("HET", vc -> Integer.toString(vc.getHetCount()));
        getters.put("HOM-REF", vc -> Integer.toString(vc.getHomRefCount()));
        getters.put("HOM-VAR", vc -> Integer.toString(vc.getHomVarCount()));
        getters.put("NO-CALL", vc -> Integer.toString(vc.getNoCallCount()));
        getters.put("TYPE", vc -> vc.getType().toString());
        getters.put("VAR", vc -> Integer.toString(vc.getHetCount() + vc.getHomVarCount()));
        getters.put("NSAMPLES", vc -> Integer.toString(vc.getNSamples()));
        getters.put("NCALLED", vc -> Integer.toString(vc.getNSamples() - vc.getNoCallCount()));
        getters.put("MULTI-ALLELIC", vc -> Boolean.toString(vc.getAlternateAlleles().size() > 1));
        getters.put("SAMPLE_NAME", vc -> vc.getGenotype(0).getSampleName());
    }

    private static Object splitAltAlleles(final VariantContext vc) {
        final int numAltAlleles = vc.getAlternateAlleles().size();
        if ( numAltAlleles == 1 ) {
            return vc.getAlternateAllele(0);
        }

        return vc.getAlternateAlleles();
    }
}
