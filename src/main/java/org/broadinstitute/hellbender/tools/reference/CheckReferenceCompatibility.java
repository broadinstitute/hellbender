package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.util.*;

/**
 * Check a BAM/VCF for compatibility against specified references.
 *
 * <p>This tool generates a table analyzing the compatibility of a sequence file against provided references. The tool works to compare
 * BAM/CRAMs (specified using the -I argument) as well as VCFs (specified using the -V argument) against provided
 * reference(s), specified using the -references-to-compare argument. When MD5s are present, the tool decides compatibility based on all sequence
 * information (MD5, name, length); when MD5s are missing, the tool makes compatibility calls based only on sequence name
 * and length. The table can be directed to a file or standard output using provided command-line arguments.
 *</p>
 *
 * <h3>Input</h3>
 * <p>
 * A BAM/CRAM/VCF and reference(s) for comparison.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A TSV file or output stream displaying the compatibility table.
 * <pre>
 * #Current Reference: reads_data_source_test1_withmd5s_missingchr1.bam
 * Reference	Compatibility	Summary
 * hg19mini.fasta	COMPATIBLE_SUBSET	The sequence dictionary in reads_data_source_test1_withmd5s_missingchr1.bam is a subset of the hg19mini.fasta reference sequence dictionary. Missing sequence(s): [1]
 * hg19mini_1renamed.fasta	COMPATIBLE_SUBSET	The sequence dictionary in reads_data_source_test1_withmd5s_missingchr1.bam is a subset of the hg19mini_1renamed.fasta reference sequence dictionary. Missing sequence(s): [chr1]
 * hg19mini_missingchr1.fasta	COMPATIBLE	The sequence dictionaries exactly match
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CheckReferenceCompatibility \
 *   -refcomp reference.fasta \
 *   -I reads.bam \
 *   -O output.table
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = "Check a BAM/VCF for compatibility against specified references and output a tab-delimited table detailing the compatibility status and relevant information about the status.",
        oneLineSummary = "Check a BAM/VCF for compatibility against specified references.",
        programGroup = ReferenceProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public class CheckReferenceCompatibility extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", optional = true)
    public GATKPath vcfPath;

    /**
     * Output file will be written here.
     *
     * Note: If no output file provided, table will print to standard output.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "If specified, output reference sequence table in TSV format to this file. Otherwise print it to stdout.", optional = true)
    private GATKPath output;

    private SAMSequenceDictionary queryDictionary;
    private GATKPath queryDictionaryPath;
    private String queryDictionaryName;
    private Map<GATKPath, SAMSequenceDictionary> dictionaries = new LinkedHashMap<>();
    private List<CompatibilityRecord> referenceCompatibilities = new ArrayList<>();
    private boolean md5sPresent;

    // sequence dictionary validation disabled since dictionaries may not be valid before running tool
    @Override
    protected SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection() {
        return new SequenceDictionaryValidationArgumentCollection.NoValidationCollection();
    }

    @Override
    public void onTraversalStart() {
        initializeSequenceDictionaryForInput();
    }

    private void initializeSequenceDictionaryForInput(){
        if(hasReads() && vcfPath != null){
            throw new UserException.BadInput("Both BAM and VCF specified. Tool analyzes one input at a time.");
        }
        else{
            // BAMs/CRAMs
            if (hasReads()) {
                if (readArguments.getReadPathSpecifiers().size() > 1) {
                    throw new UserException.BadInput("Tool analyzes one reads input at a time.");
                }
                queryDictionary = getHeaderForReads().getSequenceDictionary();
                queryDictionaryPath = readArguments.getReadPathSpecifiers().get(0);
                queryDictionaryName = queryDictionaryPath.toPath().getFileName().toString();
                md5sPresent = dictionaryHasMD5s(queryDictionary);
            }
            // VCFs
            else if(vcfPath!=null){
                try (final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(vcfPath.toString())) {
                    VCFHeader header = (VCFHeader) vcfReader.getHeader();
                    queryDictionary = header.getSequenceDictionary();
                    queryDictionaryPath = vcfPath;
                    queryDictionaryName = queryDictionaryPath.toPath().getFileName().toString();
                    md5sPresent = dictionaryHasMD5s(queryDictionary);
                }
            }
            dictionaries.put(queryDictionaryPath, queryDictionary);
            for (GATKPath path : references) {
                dictionaries.put(path, ReferenceDataSource.of(path.toPath()).getSequenceDictionary());
            }
            if (!hasReads() && vcfPath == null) {
                throw new UserException.BadInput("No input provided.");
            }
        }
    }

    @Override
    public void traverse() {
        if(md5sPresent){
            ReferenceSequenceTable table = new ReferenceSequenceTable(dictionaries);
            table.build();

            List<ReferencePair> refPairs = table.compareAgainstKeyReference(queryDictionaryPath);
            for(ReferencePair pair : refPairs){
                referenceCompatibilities.add(evaluateCompatibilityWithMD5Table(pair, table));
            }
        }
        else{
            logger.warn("************************************************************************************************************************");
            logger.warn("* Comparison lacking MD5. All comparisons based on sequence name and length, which could hide mismatching references.  *");
            logger.warn("************************************************************************************************************************");

            for(Map.Entry<GATKPath, SAMSequenceDictionary> entry : dictionaries.entrySet()){
                if(!entry.getValue().equals(queryDictionary)){
                    referenceCompatibilities.add(evaluateCompatibilityWithoutMD5(entry.getValue(), entry.getKey()));
                }
            }
        }
        writeOutput();
    }

    /**
     * Write compatibility information to output in table format.
     *
     * Note: if no output file specified, displays to standard output.
     *
     */
    private void writeOutput(){
        TableColumnCollection columns = new TableColumnCollection(Arrays.asList("Reference", "Compatibility", "Summary"));
        try(CheckReferenceCompatibilityTableWriter writer = output == null
                ? new CheckReferenceCompatibilityTableWriter(new OutputStreamWriter(System.out), columns)
                : new CheckReferenceCompatibilityTableWriter(output.toPath(), columns)
        ){
            writer.writeComment(String.format("Current Reference: %s", queryDictionaryName));
            writer.writeHeaderIfApplies();
            for(CompatibilityRecord record : referenceCompatibilities){
                writer.writeRecord(record);
            }
        }
        catch(IOException exception){
            throw new UserException.CouldNotCreateOutputFile(output, "Failed to write output table.", exception);
        }
    }

    private boolean dictionaryHasMD5s(SAMSequenceDictionary dictionary){
        for (SAMSequenceRecord record : dictionary.getSequences()) {
            String md5FromDict = record.getMd5();
            if (md5FromDict == null) {
                return false;
            }
        }
        return true;
    }

    /** Evaluate the compatibility of a given pair of references when MD5s are present.
     *
     * @param refPair ReferencePair containing the dictionary file and a reference it is being compared against
     * @param table
     */
    private CompatibilityRecord evaluateCompatibilityWithMD5Table(ReferencePair refPair, ReferenceSequenceTable table){
        EnumSet<ReferencePair.Status> status = refPair.getAnalysis();
        if(status.contains(ReferencePair.Status.EXACT_MATCH)){
            return new CompatibilityRecord(refPair.getRef2(), CompatibilityRecord.Compatibility.COMPATIBLE, "The sequence dictionaries exactly match");
        }
        else if(status.contains(ReferencePair.Status.SUBSET) && status.size() == 1){
            return new CompatibilityRecord(refPair.getRef2(), CompatibilityRecord.Compatibility.COMPATIBLE_SUBSET, String.format("The sequence dictionary in %s is a subset of the %s reference sequence dictionary. Missing sequence(s): %s", refPair.getRef1AsString(), refPair.getRef2AsString(), getMissingSequencesIfSubset(dictionaries.get(refPair.getRef2()))));
        }
        else{
            return new CompatibilityRecord(refPair.getRef2(), CompatibilityRecord.Compatibility.NOT_COMPATIBLE, String.format("Status: %s. Run CompareReferences tool for more information on reference differences.", status));
        }
    }

    /** Evaluate the compatibility of a given pair of references when MD5s are not present.
     *
     * Note: compatibility calls are made based on sequence name and length, alone, and are therefore imprecise.
     *
     * @param referenceDict SAMSequenceDictionary to be compared against the key sequence file
     * @param referenceDictPath path to provided dictionary
     */
    private CompatibilityRecord evaluateCompatibilityWithoutMD5(SAMSequenceDictionary referenceDict, GATKPath referenceDictPath){
        String referenceDictName = referenceDictPath.toPath().getFileName().toString();
        SequenceDictionaryUtils.SequenceDictionaryCompatibility compatibilityStatus = SequenceDictionaryUtils.compareDictionaries(referenceDict, queryDictionary, false);
        if(compatibilityStatus.equals(SequenceDictionaryUtils.SequenceDictionaryCompatibility.IDENTICAL)){
            return new CompatibilityRecord(referenceDictPath, CompatibilityRecord.Compatibility.MAYBE_COMPATIBLE, "All sequence names and lengths match in the sequence dictionaries. Since the MD5s are lacking, we can't confirm there aren't mismatching bases in the references.");
        }
        else if(compatibilityStatus.equals(SequenceDictionaryUtils.SequenceDictionaryCompatibility.SUPERSET)){
            return new CompatibilityRecord(referenceDictPath, CompatibilityRecord.Compatibility.COMPATIBLE_SUBSET, String.format("All sequence names and lengths present in the sequence dictionaries match, but %s is a subset of %s. Missing sequence(s): %s. Since the MD5s are lacking, we can't confirm there aren't mismatching bases in the references.", queryDictionaryName, referenceDictName, getMissingSequencesIfSubset(referenceDict)));
        }
        else{
            return new CompatibilityRecord(referenceDictPath, CompatibilityRecord.Compatibility.NOT_COMPATIBLE, String.format("Status: %s. Run CompareReferences tool for more information on reference differences.", compatibilityStatus));
        }
    }

    private List<String> getMissingSequencesIfSubset(SAMSequenceDictionary referenceDict){
        Set<String> commonSequences = SequenceDictionaryUtils.getCommonContigsByName(queryDictionary, referenceDict);
        List<String> missingSequences = SequenceDictionaryUtils.getContigNamesList(referenceDict);
        missingSequences.removeAll(commonSequences);

        return missingSequences;
    }

    /**
     * Minimal class representing a record in the compatibility output table.
     */
    private static class CompatibilityRecord {
        private final GATKPath ref;
        private final Compatibility compatibilityStatus;

        private enum Compatibility{
            COMPATIBLE,
            COMPATIBLE_SUBSET,
            NOT_COMPATIBLE,
            MAYBE_COMPATIBLE
        }
        private final String summary;

        CompatibilityRecord(GATKPath ref, Compatibility compatibilityStatus, String summary){
            this.ref = ref;
            this.compatibilityStatus = compatibilityStatus;
            this.summary = summary;
        }

        public String getRefAsString(){
            return ref.toPath().getFileName().toString();
        }

        public Compatibility getCompatibilityStatus(){
            return compatibilityStatus;
        }

        public String getSummary(){
            return summary;
        }
    }

    /**
     * TableWriter to format and write the table output.
     */
    public static class CheckReferenceCompatibilityTableWriter extends TableWriter<CompatibilityRecord> {

        public CheckReferenceCompatibilityTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
        }

        public CheckReferenceCompatibilityTableWriter(final Writer writer, TableColumnCollection columns) throws IOException {
            super(writer, columns);
        }

        @Override
        protected void composeLine(final CompatibilityRecord record, final DataLine dataLine){
            dataLine.set("Reference", record.getRefAsString())
                    .set("Compatibility", record.getCompatibilityStatus().toString())
                    .set("Summary", record.getSummary());
        }
    }
}
