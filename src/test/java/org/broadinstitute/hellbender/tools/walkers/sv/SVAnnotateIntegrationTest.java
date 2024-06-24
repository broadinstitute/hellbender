package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;


public class SVAnnotateIntegrationTest extends CommandLineProgramTest {
    final String INPUT_VCF_PATH = getToolTestDataDir() + "integration.vcf.gz";
    final File inputVCF = new File(INPUT_VCF_PATH);
    private final String LARGE_FILE_DIR = GATKBaseTest.largeFileTestDir + "SVAnnotate/";
    final File GTF_FILE = new File(LARGE_FILE_DIR + "MANE.selected.GRCh38.v0.95.select_ensembl_genomic.gtf");
    final File NONCODING_ELEMENTS_FILE = new File(LARGE_FILE_DIR + "noncoding.selected.hg38.bed.gz");

    final List<String> allAnnotationInfoKeys = Arrays.asList(GATKSVVCFConstants.LOF, GATKSVVCFConstants.PROMOTER,
            GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.NONCODING_BREAKPOINT, GATKSVVCFConstants.NONCODING_SPAN,
            GATKSVVCFConstants.INV_SPAN, GATKSVVCFConstants.PROMOTER, GATKSVVCFConstants.COPY_GAIN,
            GATKSVVCFConstants.INTERGENIC, GATKSVVCFConstants.NEAREST_TSS, GATKSVVCFConstants.INT_EXON_DUP,
            GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.UTR,
            GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.BREAKEND_EXON,
            GATKSVVCFConstants.PARTIAL_DISPERSED_DUP);

    private void assertVariantAnnotatedAsExpected(final List<VariantContext> vcf, final String variantID,
                                                  Map<String, Object> expectedAnnotations) {
        for (VariantContext variant : vcf) {
            // Iterate until specified variant ID
            if (variant.getID().equals(variantID)) {
                // check for presence of expected annotations
                for (String key : expectedAnnotations.keySet()) {
                    if (!variant.hasAttribute(key)) {
                        throw new AssertionError("Variant " + variantID + " is missing annotation " + key);
                    } else if (!variant.getAttributeAsStringList(key, null).equals(expectedAnnotations.get(key))) {
                        throw new AssertionError("Variant " + variantID + " has incorrect value for " +
                                key + "annotation. Expected: " + expectedAnnotations.get(key) + ". Actual: " +
                                variant.getAttribute(key));
                    }
                }
                // check for absence of unexpected annotations
                for (String key : allAnnotationInfoKeys) {
                    if (variant.hasAttribute(key) && !expectedAnnotations.containsKey(key)) {
                        throw new AssertionError("Variant " + variantID + " has unexpected annotation " + key);
                    }
                }
                break;
            }
        }
    }

    private void assertEqualVariantsWithIgnoredAttributes(final List<VariantContext> v1, final List<VariantContext> v2,
                                                          final List<String> attributesToIgnore) {
        // Check existence of both VCFs
        Utils.nonNull(v1, "v1");
        Utils.nonNull(v2, "v2");
        // Check input and output VCFs have same number of variants
        if (v1.size() != v2.size()){
            throw new AssertionError("different sizes " + v1.size()+ " vs " + v2.size());
        }

        boolean equalityCheckPassed = true;
        int numFailed = 0;

        // Check input and output variants are identical except for any annotation INFO keys
        for (int i = 0; i < v1.size(); i++) {
            try {
                VariantContextTestUtils.assertVariantContextsAreEqual(v1.get(i), v2.get(i), attributesToIgnore,
                        Collections.emptyList());
            } catch (AssertionError e) {
                logger.error("Variant contexts differ: " + i + ":\n" + v1.get(i) + "\n" + v2.get(i));
                equalityCheckPassed = false;
                ++numFailed;
            }
        }
        if (!equalityCheckPassed) {
            throw new AssertionError("Variant comparison failed!  Num non-matching variant pairs: " + numFailed);
        }
    }

    @DataProvider(name = "argumentSets")
    public Object[][] getArgumentSets() {
        return new Object[][]{
                // Run with GTF and noncoding BED file. Check a variant that has both
                {new ArgumentsBuilder()
                        .addVCF(inputVCF)
                        .add(SVAnnotate.PROTEIN_CODING_GTF_NAME, GTF_FILE)
                        .add(SVAnnotate.NON_CODING_BED_NAME, NONCODING_ELEMENTS_FILE),
                        "ref_panel_1kg_v1_INV_chr21_1",
                        SVAnnotateEngineUnitTest.createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.NONCODING_BREAKPOINT),
                                Arrays.asList("APP", "Tommerup_TADanno"))},
                // Noncoding only. Check same variant for just the noncoding annotation
                {new ArgumentsBuilder()
                        .addVCF(inputVCF)
                        .add(SVAnnotate.NON_CODING_BED_NAME, NONCODING_ELEMENTS_FILE),
                        "ref_panel_1kg_v1_INV_chr21_1",
                        SVAnnotateEngineUnitTest.createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_BREAKPOINT),
                                Arrays.asList("Tommerup_TADanno"))},
                // Protein-coding GTF only
                {new ArgumentsBuilder()
                        .addVCF(inputVCF)
                        .add(SVAnnotate.PROTEIN_CODING_GTF_NAME, GTF_FILE),
                        "ref_panel_1kg_v1_BND_chr22_7",
                        SVAnnotateEngineUnitTest.createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.INTRONIC), Arrays.asList("BCL2L13"))},
                // Toggle promoter window and check for addition of promoter annotation to previous variant
                {new ArgumentsBuilder()
                        .addVCF(inputVCF)
                        .add(SVAnnotate.PROTEIN_CODING_GTF_NAME, GTF_FILE)
                        .add(SVAnnotate.PROMOTER_WINDOW_NAME, 8000),
                        "ref_panel_1kg_v1_BND_chr22_7",
                        SVAnnotateEngineUnitTest.createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.PROMOTER),
                                Arrays.asList("BCL2L13", "ATP6V1E1"))},
                // Toggle BND annotation and check for change from default
                {new ArgumentsBuilder()
                        .addVCF(inputVCF)
                        .add(SVAnnotate.PROTEIN_CODING_GTF_NAME, GTF_FILE)
                        .add(SVAnnotate.MAX_BND_LEN_NAME, 12000),
                        "ref_panel_1kg_v1_BND_chr22_7",
                        SVAnnotateEngineUnitTest.createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.LOF), Arrays.asList("BCL2L13"))},
                // No reference arguments, no annotations
                {new ArgumentsBuilder()
                        .addVCF(inputVCF),
                        "ref_panel_1kg_v1_BND_chr22_7",
                        SVAnnotateEngineUnitTest.createAttributesMap( Collections.emptyList(), Collections.emptyList())}
        };
    }

    @Test(dataProvider = "argumentSets")
    public void testSVAnnotate(
            ArgumentsBuilder commandLineArgs,
            String variantID,
            Map<String, Object> expectedAnnotations
    ) {
        // Run each set of command line arguments to create a tmp annotated VCF
        final File output = createTempFile("annotated",".vcf");
        final ArgumentsBuilder args = commandLineArgs.addOutput(output);
        runCommandLine(args, SVAnnotate.class.getSimpleName());

        // Load input and output VCFs into memory for comparisons
        Pair<VCFHeader, List<VariantContext>> inputVCFHeaderAndVariants =
                VariantContextTestUtils.readEntireVCFIntoMemory(INPUT_VCF_PATH);
        Pair<VCFHeader, List<VariantContext>> outputVCFHeaderAndVariants =
                VariantContextTestUtils.readEntireVCFIntoMemory(output.getPath());
        List<VariantContext> inputVariants = inputVCFHeaderAndVariants.getRight();
        List<VariantContext> outputVariants = outputVCFHeaderAndVariants.getRight();

        // Check input and output VCFs have same variants (ignore differences in annotation INFO keys)
        assertEqualVariantsWithIgnoredAttributes(inputVariants, outputVariants, allAnnotationInfoKeys);

        // Check one variant for expected annotations
        assertVariantAnnotatedAsExpected(outputVariants, variantID, expectedAnnotations);
    }
}
