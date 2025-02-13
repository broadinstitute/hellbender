package org.broadinstitute.hellbender.tools.walkers.variantutils;

import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.beta.io.bundle.BundleJSON;
import htsjdk.beta.plugin.variants.VariantsBundle;
import htsjdk.io.IOPath;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.writers.IntervalFilteringVcfWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import com.google.common.collect.Comparators;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class SelectVariantsIntegrationTest extends CommandLineProgramTest {

    private static String baseTestString(String args, String testFile) {
        return " --variant " + testFile
                    + " -O %s "
                    + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false "
                    + args;
    }

    @Test
    public void testSampleSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + hg19MiniReference
                        + " --variant " + testFile
                        + " -sn NA11918 "
                        + " --suppress-reference-path " // suppress reference file path in output for test differencing
                        + " -O %s "
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SimpleSelection.vcf")
        );

        spec.executeTest("testSampleSelection--" + testFile, this);
    }

    @Test
    public void testExpressionSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteringDepthInFormat.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + hg19MiniReference
                        + " --variant " + testFile
                        + " -select 'DP < 7' "
                        + " --suppress-reference-path " // suppress reference file path in output for test differencing
                        + " -O %s  --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SimpleExpressionSelection.vcf")
        );

        spec.executeTest("testSimpleExpressionSelection--" + testFile, this);
    }

    @Test(expectedExceptions = UserException.ValidationFailure.class)
    public void testResortingFileWarning() throws IOException {
        final String testFile = getToolTestDataDir() + "unsortedGenotypeFieldsTestFile.vcf";
        final File output = File.createTempFile("test_unsortedGenotypeField", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(testFile)
                .addOutput(output)
                .addFlag("fail-on-unsorted-genotype");

        runCommandLine(args);

    }

    @Test
    public void testRepeatedLineSelectionAndExludeFiltered() throws IOException {
        final String testFile = getToolTestDataDir() + "test.dup.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn A -sn B -sn C -exclude-filtered ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RepeatedLineSelection.vcf")
        );

        spec.executeTest("testRepeatedLineSelection--" + testFile, this);
    }

    @Test
        public void testComplexSelection()  throws IOException {
            final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
            final String samplesFile = getToolTestDataDir() + "samples.args";

            final IntegrationTestSpec spec = new IntegrationTestSpec(
                    baseTestString(" -sn NA11894 -se 'NA069*' -sn " + samplesFile + " -select 'RMSMAPQ < 170.0'", testFile),
                    Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ComplexSelection.vcf")
            );

            spec.executeTest("testComplexSelection--" + testFile, this);
        }

    /**
     * When input variants are untrimmed, they can be trimmed by select variants, which may change their order.
     * This test confirms that this case is handled correctly, and the resulting variants are ouput correctly sorted.
     */
    @Test
    public void testUntrimmedVariants() {
        final File testFile = new File(getToolTestDataDir() + "untrimmed.vcf");
        final File output = IOUtils.createTempFile("test_untrimmed", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(testFile)
                .addOutput(output)
                .add(StandardArgumentDefinitions.SAMPLE_NAME_SHORT_NAME, "SAMPLE_01");

        runCommandLine(args);

        final List<VariantContext> vcs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getPath()).getRight();

        Assert.assertTrue(Comparators.isInOrder(vcs, Comparator.comparingInt(VariantContext::getStart)));
    }

    @Test
    public void testUntrimmedVariantsWithSetFilteredGtToNocall() {
        final File testFile = new File(getToolTestDataDir() + "untrimmed.vcf");
        final File output = IOUtils.createTempFile("test_untrimmed", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(testFile)
                .addOutput(output)
                .add(StandardArgumentDefinitions.SAMPLE_NAME_SHORT_NAME, "SAMPLE_01")
                .addFlag("set-filtered-gt-to-nocall");

        runCommandLine(args);

        final List<VariantContext> vcs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getPath()).getRight();

        Assert.assertTrue(Comparators.isInOrder(vcs, Comparator.comparingInt(VariantContext::getStart)));
    }

    @Test
    public void testComplexSelectionWithNonExistingSamples()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --allow-nonoverlapping-command-line-samples  -select 'RMSMAPQ < 170.0' -sn Z -sn " // non existent samples on command line
                        + samplesFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ComplexSelectionWithNonExistingSamples.vcf")
        );
        spec.executeTest("testComplexSelectionWithNonExistingSamples--" + testFile, this);
    }

    @Test
    public void testNonExistentSampleFile() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final File nonExistentFile = GATKBaseTest.getSafeNonExistentFile("nonexistentSamples.args");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn A -sn Z -sn Q -sn " + nonExistentFile, testFile),
                1,
                CommandLineException.class
        );
        spec.executeTest("testNonExistentSampleFile--" + testFile, this);
    }

    @Test
    public void testNonExistingFieldSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --exclude-non-variants -select 'foo!=0 || RMSMAPQ < 170.0' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_NonExistingSelection.vcf")
        );

        spec.executeTest("testNonExistingSelection--" + testFile, this);
    }

    /**
     * Test excluding samples from file and sample name
     */
    @Test
    public void testSampleExclusionFromFileAndSeparateSample()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl-sn NA11894 -xl-sn " + samplesFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionFromFileAndSeparateSample.vcf")
        );

        spec.executeTest("testSampleExclusionFromFileAndSeparateSample--" + testFile, this);
    }

    /**
     * Test excluding samples from file
     */
    @Test
    public void testSampleExclusionJustFromFile()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl-sn " + samplesFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionJustFromFile.vcf")
        );

        spec.executeTest("testSampleExclusionJustFromFile--" + testFile, this);
    }

    /**
     * Test excluding samples from expression
     */
    @Test
    public void testSampleExclusionJustFromExpression()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl-se 'NA069*' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionJustFromExpression.vcf")
        );

        spec.executeTest("testSampleExclusionJustFromExpression--" + testFile, this);
    }

    /**
     * Test excluding samples from negation expression
     */
    @Test
    public void testSampleExclusionJustFromNegationExpression()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -se 'NA[0-9]{4}[^1-9]' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionJustFromRegexExpression.vcf")
        );

        spec.executeTest("testSampleExclusionJustFromRegexExpression--" + testFile, this);
    }

    /**
     * Test including samples that are not in the VCF
     */

    @Test
    public void testSampleInclusionWithNonexistingSamples()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn A -sn Z -sn Q -sn " + samplesFile, testFile),
                1,
                UserException.BadInput.class
        );

        spec.executeTest("testSampleInclusionWithNonexistingSamples--" + testFile, this);
    }

    @Test
    public void testDiscordance() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String discordanceFile = getToolTestDataDir() + "vcfexample2DiscordanceConcordance.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11992 " // not present in discordance track
                                + " -disc " + discordanceFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Discordance.vcf")
        );

        spec.executeTest("testDiscordance--" + testFile, this);
    }

    @Test
    public void testConcordance()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String concordanceFile = getToolTestDataDir() + "vcfexample2DiscordanceConcordance.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -conc " + concordanceFile  + " --lenient ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Concordance.vcf")
        );

        spec.executeTest("testConcordance--" + testFile, this);
    }

    /**
     * Test including variant types.
     */
    @Test
    public void testVariantTypeSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --restrict-alleles-to MULTIALLELIC --select-type-to-include MIXED ",testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_VariantTypeSelection.vcf")
        );

        spec.executeTest("testVariantTypeSelection--" + testFile, this);
    }

    /**
     * Test including variant types in GVCF files.
     */
    @Test
    public void testVariantTypeSelectionForGVCF() throws IOException {
        final String testFile = getToolTestDataDir() + "gvcfExample.g.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --select-type-to-include SNP --ignore-non-ref-in-types ",testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_VariantTypeSelectionForGVCF.vcf")
        );

        spec.executeTest("testVariantTypeSelectionForGVCF--" + testFile, this);
    }

    /**
     * Test excluding indels that are larger than the specified size
     */
    @Test
    public void testMaxIndelLengthSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --select-type-to-include INDEL --max-indel-size 2 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MaxIndelLengthSelection.vcf")
        );

        spec.executeTest("testMaxIndelLengthSelection--" + testFile, this);
    }

    /**
     * Test excluding indels that are smaller than the specified size
     */
    @Test
    public void testMinIndelLengthSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --select-type-to-include INDEL --min-indel-size 2 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MinIndelLengthSelection.vcf")
        );

        spec.executeTest("testMinIndelLengthSelection--" + testFile, this);
    }

    @Test
    public void testRemoveMLE() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.withMLE.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12892 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RemoveMLE.vcf")
        );

        spec.executeTest("testRemoveMLE--" + testFile, this);
    }

    @Test
    public void testKeepOriginalAC() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.loseAlleleInSelection.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --keep-original-ac -sn NA12892 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepOriginalAC.vcf")
        );

        spec.executeTest("testKeepOriginalAC--" + testFile, this);
    }

    @Test
    public void testKeepOriginalACAndENV() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.loseAlleleInSelection.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --keep-original-ac -sn NA12892 --exclude-non-variants --remove-unused-alternates ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepOriginalACAndENV.vcf")
        );

        spec.executeTest("testKeepOriginalACAndENV--" + testFile, this);
    }

    @Test
    public void testKeepOriginalDP() throws IOException {
        final String testFile = getToolTestDataDir() + "CEUtrioTest.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --keep-original-dp -sn NA12892 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepOriginalDP.vcf")
        );

        spec.executeTest("testKeepOriginalDP--" + testFile, this);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() throws IOException {
        final String testFile = getToolTestDataDir() + "selectVariants.onePosition.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -select 'KG_FREQ < 0.5' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MultipleRecordsAtOnePosition.vcf")
        );

        spec.executeTest("testMultipleRecordsAtOnePosition--" + testFile, this);
    }

    @Test
    public void testNoGTs() throws IOException {
        final String testFile = getToolTestDataDir() + "vcf4.1.example.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec (
                " --variant " + testFile + " -O %s  --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_NoGTs.vcf")
        );

        spec.executeTest("testNoGTs--" + testFile, this);
    }

    @Test
    public void testRemoveSingleSpanDelAlleleNoSpanDel() throws IOException {
        final String testFile = getToolTestDataDir() + "spanning_deletion.vcf";
        final String sampleName = "NA1";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn " + sampleName + " --remove-unused-alternates --exclude-non-variants", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RemoveSingleSpanDelAlleleNoSpanDel.vcf")
        );
        spec.executeTest("test encounter no instance of '*' as only ALT allele and ensure line is removed when only monomorphic allele exists" + testFile, this);
    }

    @Test
    public void testRemoveSingleSpanDelAlleleExNonVar() throws IOException {
        final String testFile = getToolTestDataDir() + "spanning_deletion.vcf";
        final String sampleName = "NA2";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn " + sampleName + " --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RemoveSingleSpanDelAlleleExNoVar.vcf")
        );
        spec.executeTest("test will not remove variant line where '*' is only ALT allele because --exclude-non-variants not called --" + testFile, this);
    }

    @Test
    public void testRemoveSingleSpanDelAllele() throws IOException {
        final String testFile = getToolTestDataDir() + "spanning_deletion.vcf";
        final String sampleName = "NA2";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn " + sampleName + " --exclude-non-variants --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RemoveSingleSpanDelAllele.vcf")
        );
        spec.executeTest("test removes variant line where '*' is only ALT allele --" + testFile, this);
    }

    @Test
    public void testSelectFromMultiAllelic() throws IOException {
        final String testFile = getToolTestDataDir() + "multi-allelic.bi-allelicInGIH.vcf";
        final String sampleName = getToolTestDataDir() + "GIH.samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn " + sampleName + " --exclude-non-variants --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MultiAllelicExcludeNonVar.vcf")
        );
        spec.executeTest("test select from multi allelic with exclude-non-variants --" + testFile, this);
    }

    @Test
    public void testMultiAllelicAnnotationOrdering() throws IOException {
        final String testFile = getToolTestDataDir() + "multi-allelic-ordering.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn SAMPLE-CC -sn SAMPLE-CT -sn SAMPLE-CA --exclude-non-variants", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MultiAllelicAnnotationOrdering.vcf")
        );
        spec.executeTest("test multi allelic annotation ordering --" + testFile, this);
    }

    @Test
    public void testFileWithoutInfoLineInHeader() throws IOException {
        testFileWithoutInfoLineInHeader("testSelectVariants_FileWithoutInfoLineInHeader", IllegalStateException.class);
    }

    @Test
    public void testFileWithoutInfoLineInHeaderWithOverride() throws IOException {
        testFileWithoutInfoLineInHeader("testSelectVariants_FileWithoutInfoLineInHeaderWithOverride", null);
    }

    private void testFileWithoutInfoLineInHeader(final String name, final Class<? extends Exception> expectedException) throws IOException {
        final String testFile = getToolTestDataDir() + "missingHeaderLine.vcf";
        final String outFile = getToolTestDataDir() + "expected/" + name + ".vcf";

        final String cmd = baseTestString(" -sn NA12892 " + (expectedException == null ? " --lenient" : ""), testFile);

        IntegrationTestSpec spec =
                expectedException != null
                        ? new IntegrationTestSpec(cmd, 1, expectedException)
                        : new IntegrationTestSpec(cmd, Collections.singletonList(outFile));

        spec.executeTest(name, this);
    }

    private void beforeJexlTest(){
        // NOTE: JexlEngine singleton construction in VariantContextUtils sets silent to false.
        // However VariantFiltration.initialize() sets setSilent(true) on the shared instance.
        // Just in case this test runs after a VariantFiltration in the same VM, always set silent back to false.
        htsjdk.variant.variantcontext.VariantContextUtils.engine.get().setSilent(false);
    }

    final GATKPath svHGDBMultiSampleVcf = new GATKPath(getToolTestDataDir() + "hgdb_sv_multi_sample_mini.vcf");
    final Pair<VCFHeader, List<VariantContext>> svHGDBMultiSampleVariants = VariantContextTestUtils.readEntireVCFIntoMemory(svHGDBMultiSampleVcf.toString());
    final int gqThreshold = 60; // For JEXL filtering by sample GQ

    @Test
    public void testInvalidJexl() throws IOException {
        beforeJexlTest();

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -select 'vc.getGenotype(\"FAKE_SAMPLE\").isHomRef()' ", svHGDBMultiSampleVcf.toString()),
                1,
                UserException.class);
        spec.executeTest("InvalidJexl", this);
    }

    // If two separate -select command lines are used e.g. --select "'AF > 0.05'" --select "'MQRankSum > 3'"
    // then the booleans should be combined with the logical-or
    @Test
    public void testMultipleINFOJexls() {
        beforeJexlTest();
        final File testOutput = createTempFile("multiple_jexl_test", "vcf");
        final Predicate<VariantContext> filter = vc -> vc.getAttributeAsInt("AC", -1) > 0 ||
                vc.getAttributeAsString("SVTYPE", "").equals("DEL");

        final int expectedCount = (int) svHGDBMultiSampleVariants.getRight().stream().filter(filter).count();
        Assert.assertTrue(expectedCount > 0); // Otherwise the test is not informative

        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "-select", "SVTYPE == 'DEL'",
                        "-select", "AC > 0"),
                SelectVariants.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualCount = result.getRight().size();
        Assert.assertEquals(expectedCount, actualCount);
    }

    // Test JEXL expressions like "AC > 0 && AF > 0.001"
    // The recommended usage for logical-or is to place individual expressions as separate --select command
    // e.g. -select AC > 0 -select AF > 0.001 for "AC > 0 || AF > 0.001"
    @Test
    public void testLogicalAndOrInJexl() {
        beforeJexlTest();
        final File testOutput = createTempFile("multiple_jexl_test", "vcf");

        final Predicate<VariantContext> andFilter = vc -> vc.getAttributeAsInt("AC", -1) > 0 &&
                vc.getAttributeAsString("SVTYPE", "").equals("DEL");
        final int expectedCount = (int) svHGDBMultiSampleVariants.getRight().stream().filter(andFilter).count();
        Assert.assertTrue(expectedCount > 0);

        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "--select", "SVTYPE == 'DEL' && AC > 0"),
                SelectVariants.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualCount = result.getRight().size();
        Assert.assertEquals(expectedCount, actualCount);
    }


    /** One of the five variants does not have the GQ field. It is silently removed from the output. **/
    @Test
    public void testJexlGenotypeFilter() {
        final File testOutput = createTempFile("jexl_genotype_test", "vcf");
        beforeJexlTest();
        final int threshold = 1;
        final Predicate<VariantContext> filter = vc -> vc.getGenotypes().stream().anyMatch(g -> g.getGQ() > threshold);
        final int expectedNumVCs = (int) svHGDBMultiSampleVariants.getRight().stream().filter(filter).count();
        Assert.assertTrue(expectedNumVCs > 0);

        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "-" + SelectVariants.GENOTYPE_SELECT_SHORT_NAME,"GQ > " + threshold),
                SelectVariants.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualNumVCs = result.getRight().size();

        Assert.assertEquals(expectedNumVCs, actualNumVCs);
    }

    @Test
    public void testJexlDirectlyAccessGenotype() {
        final File testOutput = createTempFile("jexl_genotype_test", "vcf");
        beforeJexlTest();
        final String sampleName = "__HGDP00029__";
        final int threshold = 1;
        final Predicate<VariantContext> filter = vc -> vc.getGenotype(sampleName).getGQ() > threshold;
        final int expectedNumVCs = (int) svHGDBMultiSampleVariants.getRight().stream().filter(filter).count();
        Assert.assertTrue(expectedNumVCs > 0);

        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "-" + SelectVariants.SELECT_NAME, "vc.getGenotype('__HGDP00029__').getGQ() > " + threshold),
                SelectVariants.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualNumVCs = result.getRight().size();

        Assert.assertEquals(expectedNumVCs, actualNumVCs);
    }

    @Test
    public void testUseCorrectDPInJexl() {
        beforeJexlTest();
        final GATKPath haploidMultiSampleVcf = new GATKPath(getToolTestDataDir() + "haploid-multisample.vcf");
        final Pair<VCFHeader, List<VariantContext>> haploidMultiSampleVcfVariants = VariantContextTestUtils.readEntireVCFIntoMemory(haploidMultiSampleVcf.toString());
        final File testOutputInfo = createTempFile("jexl_info", "vcf");
        final File testOutputGenotype = createTempFile("jexl_genotype", "vcf");

        final int threshold = 20;
        final Predicate<VariantContext> filterInfoDep = vc -> vc.getAttributeAsInt("DP", -1) > threshold;
        final Predicate<VariantContext> filterGenotypeDep = vc -> vc.getGenotypes().stream().anyMatch(g -> g.getDP() > threshold);

        final int expectedNumVCsInfo = (int) haploidMultiSampleVcfVariants.getRight().stream().filter(filterInfoDep).count();
        final int expectedNumVCsGenotype = (int) haploidMultiSampleVcfVariants.getRight().stream().filter(filterGenotypeDep).count();


        runCommandLine(Arrays.asList("-V", haploidMultiSampleVcf.toString(),
                        "-O", testOutputInfo.getAbsolutePath(),
                        "-select","DP > 20"),
                SelectVariants.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> resultInfo = VariantContextTestUtils.readEntireVCFIntoMemory(testOutputInfo.getAbsolutePath());
        final int actualNumVCsInfo = resultInfo.getRight().size();
        Assert.assertEquals(expectedNumVCsInfo, actualNumVCsInfo);

        runCommandLine(Arrays.asList("-V", haploidMultiSampleVcf.toString(),
                        "-O", testOutputGenotype.getAbsolutePath(),
                        "-select-genotype","DP > 20"),
                SelectVariants.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> resultGenotype = VariantContextTestUtils.readEntireVCFIntoMemory(testOutputGenotype.getAbsolutePath());
        final int actualNumVCsGenotype = resultGenotype.getRight().size();
        Assert.assertEquals(expectedNumVCsGenotype, actualNumVCsGenotype);
    }

    // Combining -select and -select-genotype is done by logical-or by default.
    @Test
    public void testCombineInfoAndGenotypeJexl() {
        final File testOutput = createTempFile("jexl_genotype_info_test", "vcf");
        beforeJexlTest();
        final Predicate<VariantContext> filter = vc -> vc.getAttributeAsString("SVTYPE", "").equals("DEL") ||
                vc.getGenotypes().stream().anyMatch(g -> g.getGQ() > gqThreshold);
        final int expectedNumVCs = (int) svHGDBMultiSampleVariants.getRight().stream().filter(filter).count();

        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "-" + SelectVariants.SELECT_NAME, "SVTYPE == 'DEL'",
                        "-" + SelectVariants.GENOTYPE_SELECT_SHORT_NAME, "GQ > " + gqThreshold),
                SelectVariants.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualNumVCs = result.getRight().size();
        Assert.assertEquals(expectedNumVCs, actualNumVCs);
    }

    // Erroneously give an INFO field to -select-genotype
    // This (unfortunately) works, since we give the whole VariantContext object to the JEXL parser (so it can find info fields)
    // when we filter by the genotype field.
    // Adding the capability to detect such a mix-up should be added in the future.
    @Test
    public void testINFOFieldGivenToSelectGenotype() {
        final Predicate<VariantContext> filter = vc -> vc.getAttributeAsString("SVTYPE", "").equals("DEL");
        final int expectedNumVCs = (int) svHGDBMultiSampleVariants.getRight().stream().filter(filter).count();
        Assert.assertTrue(expectedNumVCs > 0);

        final File testOutput = createTempFile("jexl_genotype_info_mixup_test", "vcf");
        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "-" + SelectVariants.GENOTYPE_SELECT_SHORT_NAME, "SVTYPE == 'DEL'"), // INFO filter expression given to genotype filter
                SelectVariants.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualNumVCs2 = result.getRight().size();
        Assert.assertEquals(expectedNumVCs, actualNumVCs2);
    }

    // Erroneously give a genotype field expression --select
    // The tool silently filters everything (i.e. output an empty vcf), because we don't give the genotype object to the jexl parser,
    // and so the parser cannot find the requested genotype fields.
    // Adding the capability to detect such a mix-up should be added in the future.
    @Test
    public void testFORMATFieldGivenToSelect() {
        final File testOutput = createTempFile("jexl_genotype_info_mixup_test", "vcf");
        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutput.getAbsolutePath(),
                        "-" + SelectVariants.SELECT_NAME,"GQ > " + gqThreshold), // genotype filter expression given to INFO filter
                SelectVariants.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath());
        final int actualNumVCs3 = result.getRight().size();
        Assert.assertEquals(0, actualNumVCs3);
    }

    @Test
    public void testSelectRandomFraction() {
        beforeJexlTest();
        final File testOutput = createTempFile("jexl_genotype_info_test", "vcf");
        // There are 80 variants in this vcf
        final File inputVcf = getTestFile("tetraploid-multisample-sac.g.vcf");
        final int numInputVariants = VariantContextTestUtils.readEntireVCFIntoMemory(inputVcf.getAbsolutePath()).getRight().size();
        final double fractionToKeep = 0.5;

        runCommandLine(Arrays.asList("-V", inputVcf.getAbsolutePath(),
                        "-O", testOutput.getAbsolutePath(),
                        "-" + SelectVariants.FRACTION_TO_KEEP_SHORT_NAME, Double.toString(fractionToKeep)),
                SelectVariants.class.getSimpleName());

        final int numOutputVariants = VariantContextTestUtils.readEntireVCFIntoMemory(testOutput.getAbsolutePath()).getRight().size();
        final double fractionKept = (double) numOutputVariants/numInputVariants;
        final double epsilon = 0.12; // Not chosen rigorously.
        Assert.assertEquals(MathUtils.compareDoubles(fractionToKeep, fractionKept, epsilon), 0);
    }

    // Check that the applying JEXL filter before the sample subsetting does not affect the final output
    @Test
    public void testJEXLFilterFirst() {
        beforeJexlTest();
        final File testOutputJexlBefore = createTempFile("jexl_before_test", "vcf");
        final File testOutputJexlAfter = createTempFile("jexl_after_test", "vcf");

        final Predicate<VariantContext> filter = vc -> vc.getAttributeAsInt("AC", -1) > 0 ||
                vc.getAttributeAsString("SVTYPE", "").equals("DEL");

        final List<VariantContext> expectedVariants = svHGDBMultiSampleVariants.getRight().stream().filter(filter).collect(Collectors.toList());
        Assert.assertTrue(expectedVariants.size() > 0);

        // Run the JEXL filters before subsetting samples
        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutputJexlBefore.getAbsolutePath(),
                        "-select", "SVTYPE == 'DEL'",
                        "-select", "AC > 0",
                        "-" + SelectVariants.APPLY_JEXL_FIRST_SHORT_NAME, "true"),
                SelectVariants.class.getSimpleName());

        // Run then after subsetting samples
        runCommandLine(Arrays.asList("-V", svHGDBMultiSampleVcf.toString(),
                        "-O", testOutputJexlAfter.getAbsolutePath(),
                        "-select", "SVTYPE == 'DEL'",
                        "-select", "AC > 0"),
                SelectVariants.class.getSimpleName());
        final List<VariantContext> resultJexlBefore = VariantContextTestUtils.readEntireVCFIntoMemory(testOutputJexlBefore.getAbsolutePath()).getRight();
        final List<VariantContext> resultJexlAfter = VariantContextTestUtils.readEntireVCFIntoMemory(testOutputJexlAfter.getAbsolutePath()).getRight();

        Assert.assertEquals(expectedVariants.size(), resultJexlBefore.size());
        Assert.assertEquals(resultJexlBefore.size(), resultJexlAfter.size());
        for (int i = 0; i < expectedVariants.size(); i++){
            final VariantContext expectedVariant = expectedVariants.get(i);
            final VariantContext variantBefore = resultJexlBefore.get(i);
            final VariantContext variantAfter = resultJexlAfter.get(i);
            Assert.assertEquals(expectedVariant.getStart(), variantBefore.getStart());
            Assert.assertEquals(variantBefore.getStart(), variantAfter.getStart());
        }
    }

    @Test
    public void testAlleleTrimming() throws IOException {
        final String testFile = getToolTestDataDir() + "forHardLeftAlignVariantsTest.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12878 --exclude-non-variants --remove-unused-alternates ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_AlleleTrimming.vcf"));
        spec.executeTest("testAlleleTrimming", this);
    }

    @DataProvider(name="unusedAlleleTrimmingProvider")
    public Object[][] unusedAlleleTrimmingProvider() {
        final String expectedPath = getToolTestDataDir() + "expected/";
        return new Object[][] {
                {
                        getToolTestDataDir() + "forHardLeftAlignVariantsTest.vcf",
                        "--remove-unused-alternates",
                        expectedPath + "testSelectVariants_UnusedAlleleHardLeftTrim.vcf"
                },
                {
                        getToolTestDataDir() + "forHardLeftAlignVariantsTest.vcf",
                        null,
                        expectedPath + "testSelectVariants_UnusedAlleleHardLeft.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCT.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT --exclude-non-variants",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCTEnv.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT --remove-unused-alternates",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCTTrim.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT --exclude-non-variants --remove-unused-alternates",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCTTrimAltEnv.vcf"
                }
        };
    }

    @Test(dataProvider="unusedAlleleTrimmingProvider")
    public void testUnusedAlleleTrimming(final String vcf, final String extraArgs, final String expectedOutput) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(extraArgs == null ? "" : extraArgs, vcf),
                Collections.singletonList(expectedOutput)
        );

        spec.executeTest(
                String.format("testUnusedAlleleTrimming: (%s,%s)", new File(vcf).getName(), extraArgs == null ? "(none)" : extraArgs),
                this);
    }

    /**
     *  Test with an empty VCF file
     */
    @Test
    public void testEmptyVcfException() throws IOException {
        final String testFile = getToolTestDataDir() + "reallyEmpty.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("", testFile),
                1,
                UserException.NoSuitableCodecs.class
        );

        spec.executeTest("testEmptyVcfException--" + testFile, this);
    }

    /**
     * Test with a VCF file that is not a file
     */
    @Test
    public void testNotFileVcfException() throws IOException {
        final String testFile = getToolTestDataDir();

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("", testFile),
                1,
                UserException.CouldNotReadInputFile.class
        );

        spec.executeTest("testNotFileVcfException--" + testFile, this);
    }

    /**
     * Test with a VCF file that does not exist
     */
    @Test
    public void testMissingVcfException() throws IOException {
        final String testFile = "test.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("", testFile),
                1,
                UserException.CouldNotReadInputFile.class
        );

        spec.executeTest("testMissingVcfException--" + testFile, this);
    }

    /**
     * Test inverting the variant selection criteria by the -invertSelect argument
     */
    @Test
    public void testInvertSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -sn " + samplesFile +
                                    " -select 'RMSMAPQ < 170.0' --invert-select ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_InvertSelection.vcf")
        );

        spec.executeTest("testInvertSelection--" + testFile, this);
    }

    /**
     * Test inverting the variant selection criteria by inverting the JEXL expression logic following -select
     */
    @Test
    public void testInvertJexlSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -sn " + samplesFile +
                        " -select 'RMSMAPQ > 170.0' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_InvertJexlSelection.vcf")
        );

        spec.executeTest("testInvertJexlSelection--" + testFile, this);
    }

    /**
     * Test selecting variants with rsIDs from a .list file
     */
    @Test
    public void testKeepSelectionIDFromFile() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";
        final String idFile = getToolTestDataDir() + "complexExample1.vcf.id.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -ids " + idFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepSelectionID.vcf")
        );

        spec.executeTest("testKeepSelectionIDFile--" + testFile, this);
    }

    /**
     * Test selecting variants with literal rsIDs
     */
    @Test
    public void testKeepSelectionIDLiteral() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -ids testid1", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepSelectionID.vcf")
        );

        spec.executeTest("testKeepSelectionIDLiteral--" + testFile, this);
    }

    /**
     * Test excluding variants with rsIDs from a file
     */
    @Test
    public void testExcludeSelectionIDFromFile() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";
        final String idFile = getToolTestDataDir() + "complexExample1.vcf.id.args";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl-ids " + idFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ExcludeSelectionID.vcf")
        );

        spec.executeTest("testExcludeSelectionIDFile--" + testFile, this);
    }

    /**
     * Test excluding variants with literal rsIDs
     */
    @Test
    public void testExcludeSelectionIDLiteral() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl-ids testid1", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ExcludeSelectionID.vcf")
        );

        spec.executeTest("testExcludeSelectionIDLiteral--" + testFile, this);
    }

    /**
     * Test excluding variant types
     */
    @Test
    public void testExcludeSelectionType() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --select-type-to-exclude SNP ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ExcludeSelectionType.vcf")
        );

        spec.executeTest("testExcludeSelectionType--" + testFile, this);
    }

    @Test
    public void testMendelianViolationSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "CEUtrioTest.vcf";
        final String pedFile = getToolTestDataDir() + "CEUtrio.ped";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -ped " + pedFile + " --mendelian-violation --mendelian-violation-qual-threshold 0 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MendelianViolationSelection.vcf")
        );

        spec.executeTest("testMendelianViolationSelection--" + testFile, this);
    }

    @Test
    public void testInvertMendelianViolationSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "CEUtrioTest.vcf";
        final String pedFile = getToolTestDataDir() + "CEUtrio.ped";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --mendelian-violation --mendelian-violation-qual-threshold 0 --invert-mendelian-violation -ped " + pedFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_InvertMendelianViolationSelection.vcf")
        );

        spec.executeTest("testInvertMendelianViolationSelection--" + testFile, this);
    }

    @Test
    public void testMaxFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --max-filtered-genotypes 1 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MaxFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMaxFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testMinFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --min-filtered-genotypes 2 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MinFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMinFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testMaxFractionFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --max-fraction-filtered-genotypes 0.4 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MaxFractionFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMaxFractionFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testMinFractionFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final  IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --min-fraction-filtered-genotypes 0.6 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MinFractionFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMinFractionFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testSetFilteredGtoNocall() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --set-filtered-gt-to-nocall ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SetFilteredGtoNocall.vcf")
        );

        spec.executeTest("testSetFilteredGtoNocall--" + testFile, this);
    }

    @Test
    public void testMaxNoCall1() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.forNoCallFiltering.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --max-nocall-number 1", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_maxNOCALLnumber1.vcf")
        );

        spec.executeTest("testMaxNoCall1--" + testFile, this);
    }

    @Test
    public void testMaxNoCall0_25() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.forNoCallFiltering.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --max-nocall-fraction 0.25", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_maxNOCALLnumber1.vcf")
        );

        spec.executeTest("testMaxNoCall0_25--" + testFile, this);
    }

    @Test
    public void testMaxNoCall2() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.forNoCallFiltering.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --max-nocall-number 2", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_maxNOCALLnumber2.vcf")
        );

        spec.executeTest("testMaxNoCall2--" + testFile, this);
    }

    @Test
    public void testMaxNoCall0_5() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.forNoCallFiltering.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --max-nocall-fraction 0.5", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_maxNOCALLnumber2.vcf")
        );

        spec.executeTest("testMaxNoCall0_5--" + testFile, this);
    }

    @Test
    public void testHaploid() throws IOException {
        final String testFile = getToolTestDataDir() + "haploid-multisample.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn HG00610 -select 'DP > 7' --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Haploid.vcf")
        );

        spec.executeTest("testHaploid--" + testFile, this);
    }

    @Test
    public void testTetraploid() throws IOException {
        final String testFile = getToolTestDataDir() + "tetraploid-multisample.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA18486 -select 'DP > 57' --remove-unused-alternates ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Tetraploid.vcf")
        );

        spec.executeTest("testTetraploid--" + testFile, this);
    }

    @Test
    public void testTetraDiploid() throws IOException {
        final String testFile = getToolTestDataDir() + "tetra-diploid.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12878 -select 'DP > 48' --remove-unused-alternates ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_TetraDiploid.vcf")
        );

        spec.executeTest("testTetraDiploid--" + testFile, this);
    }
    
    @Test
    public void testSACSimpleDiploid() throws IOException {
        final String testFile = getToolTestDataDir() + "261_S01_raw_variants_gvcf.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SimpleDiploid.vcf")
        );

        spec.executeTest("testSACSimpleDiploid" + testFile, this);
    }

    @Test
    public void testSACDiploid() throws IOException {
        final String testFile = getToolTestDataDir() + "diploid-multisample-sac.g.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12891 --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SACDiploid.vcf")
        );

        spec.executeTest("testSACDiploid" + testFile, this);
    }

    @Test
    public void testSACNonDiploid() throws IOException {
        final String testFile = getToolTestDataDir() + "tetraploid-multisample-sac.g.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12891 --remove-unused-alternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SACNonDiploid.vcf")
        );

        spec.executeTest("testSACNonDiploid" + testFile, this);
    }

    @Test
    public void testSetFilteredGtoNocallUpdateInfo() throws IOException {
        final String testFile = getToolTestDataDir() + "selectVariantsInfoField.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --set-filtered-gt-to-nocall --remove-unused-alternates --exclude-non-variants", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SetFilteredGtoNocallUpdateInfo.vcf")
        );

        spec.executeTest("testSetFilteredGtoNocallUpdateInfo--" + testFile, this);
    }

    @DataProvider(name = "dropAnnotationsDataProvider")
    Object[][] dropAnnotationsDataProvider() {
        return new Object[][]{
                {"-DA FisherStrand -DA OnOffGenotype -DGA RD -sn NA11894", "testSelectVariants_DropAnnotations.vcf", "standard"},
                {"-DA FisherStrand -DA OnOffGenotype -DGA RD -sn NA11894 -DA NotAnAnnotation -DGA AlsoNotAnAnnotation", "testSelectVariants_DropAnnotations.vcf", "unused_annotations"},
                {"-DA FisherStrand -DA OnOffGenotype -DGA RD -sn NA11894 -select 'FisherStrand > 10.0'", "testSelectVariants_DropAnnotationsSelectFisherStrand.vcf", "select_on_dropped_annotation"},
                {"-DA FisherStrand -DA OnOffGenotype -DGA RD -sn NA11894 -select 'RMSMAPQ > 175.0'", "testSelectVariants_DropAnnotationsSelectRMSMAPQ.vcf", "select_on_kept_annotation"},
                {"-DA FisherStrand -DA OnOffGenotype -DGA RD -sn NA11894 -select 'vc.getGenotype(\"NA11894\").getExtendedAttribute(\"RD\")>6'", "testSelectVariants_DropAnnotationsSelectRD.vcf", "select_on_dropped_genotype_annotation"},
                {"-DA FisherStrand -DA OnOffGenotype -DGA RD -sn NA11894 -select 'vc.getGenotype(\"NA11894\").getGQ()==1'", "testSelectVariants_DropAnnotationsSelectGQ.vcf", "select_on_kept_genotype_annotation"}
        };
    }

    @Test(dataProvider = "dropAnnotationsDataProvider")
    public void testDropAnnotations(String args, String expectedFile, String testName) throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(args, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + expectedFile)
        );
        spec.executeTest("testDropAnnotations--" + testName, this);
    }

    @Test(groups = "bucket")
    public void testSampleSelectionOnNio() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final String out = BucketUtils.getTempFilePath(
            getGCPTestStaging() +"testSelectVariants_SimpleSelection", ".vcf");

        final String[] args = new String[]{
            "SelectVariants",
            "-R", hg19MiniReference
            , "--variant", testFile
            , "-sn", "NA11918"
            , "--suppress-reference-path" // suppress reference file path in output for test differencing
            , "-O", out
            , "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"};

        final String expectedFile = getToolTestDataDir() + "expected/" + "testSelectVariants_SimpleSelection.vcf";

        new Main().instanceMain(args);

        IntegrationTestSpec.assertEqualTextFiles(IOUtils.getPath(out), IOUtils.getPath(expectedFile), null);
    }

    @DataProvider(name="cloudBucketBundles")
    public Object[][] cloudBucketBundles() {
        return new Object[][]{
                {
                    // cloud vcf with a .idx
                        new GATKPath("gs://hellbender/test/resources/large/dbsnp_138.b37.1.1-65M.vcf"),
                        new GATKPath("gs://hellbender/test/resources/large/dbsnp_138.b37.1.1-65M.vcf.idx"),
                        "1:20547",
                        "expected/testSelectVariants_BundleWith_idx.vcf"
                },
                {
                    // cloud vcf with a .tbi
                        new GATKPath("gs://hellbender/test/resources/large/dbsnp_138.b37.20.21.vcf.blockgz.gz"),
                        new GATKPath("gs://hellbender/test/resources/large/dbsnp_138.b37.20.21.vcf.blockgz.gz.tbi"),
                        "21:10000128",
                        "expected/testSelectVariants_BundleWith_tbi.vcf"
                }
        };
    }

    @Test(dataProvider="cloudBucketBundles", groups = "bucket")
    public void testBundleRemoteQueryWithVCFAndIndexInSeparateBuckets(
            final IOPath vcfPath,
            final IOPath indexPath,
            final String queryInterval,
            final String expectedOutputFile
    ) throws IOException {
        // test that a query using a bundle containing remote nio paths to a vcf works when the index is
        // in a different bucket than the vcf
        final String targetBucketName = BucketUtils.randomRemotePath(getGCPTestStaging(), "testSelectVariantsBundleOnNio", "") + "/";
        IOUtils.deleteOnExit(IOUtils.getPath(targetBucketName));
        final GATKPath targetVCFPath = new GATKPath(targetBucketName + vcfPath.toPath().getFileName());

        // copy our test vcf (which is already in a bucket with it's accompanying index) to a different bucket where the
        // index does not reside, and then create a bundle referencing the copy of the vcf and the original index, so
        // the vcf and the index are in separate buckets
        Files.copy(vcfPath.toPath(), targetVCFPath.toPath(), StandardCopyOption.REPLACE_EXISTING);

        // create a bundle with the copied vcf and the original index
        final VariantsBundle vcfBundle = new VariantsBundle(new GATKPath(targetVCFPath), indexPath);
        final GATKPath bundleFilePath = new GATKPath(createTempFile("testSelectVariantsBundleOnNio", ".json").toString());
        IOPathUtils.writeStringToPath(bundleFilePath, BundleJSON.toJSON(vcfBundle));

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + bundleFilePath
                        + " -L " + queryInterval
                        + " --suppress-reference-path " // suppress reference file path in output for test differencing
                        + " -O %s "
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(getToolTestDataDir() + expectedOutputFile)
        );

        spec.executeTest("testBundle", this);
    }

    @Test(groups = "bucket")
    public void testBundleRemoteWithQueryParams() throws IOException {
        final VariantsBundle vcfBundle = new VariantsBundle(
                new GATKPath(BucketUtils.createSignedUrlToGcsObject("gs://hellbender/test/resources/large/1000G.phase3.broad.withGenotypes.chr20.10100000.vcf", 1)),
                new GATKPath(BucketUtils.createSignedUrlToGcsObject("gs://hellbender/test/resources/large/1000G.phase3.broad.withGenotypes.chr20.10100000.vcf.idx", 1)));
        final IOPath bundleFile = new GATKPath(createTempFile("testSelectVariantsThroughBundle", ".json").toString());
        IOPathUtils.writeStringToPath(bundleFile, BundleJSON.toJSON(vcfBundle));

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + bundleFile.toString()
                        + " -L 20:10001365 "
                        + " --suppress-reference-path " // suppress reference file path in output for test differencing
                        + " -O %s "
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Bundle.vcf")
        );

        spec.executeTest("testBundle", this);
    }

    // the input test file is a somatic VCF with several many-allelic sites and no PLs.  This tests that the tool does not attempt
    // to create a PL-to-alleles cache, which would cause the tool to freeze.  See https://github.com/broadinstitute/gatk/issues/6291
    @Test
    public void testManyAllelicWithoutPLsDoesntFreeze() {
        final File input = new File(getToolTestDataDir(), "many-allelic-somatic.vcf");
        final File output = createTempFile("output", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(input)
                .addReference(b37Reference)
                .addOutput(output);
        runCommandLine(args);
    }

    @Test
    public void testOutputOnlyStartingInInterval() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("-L 1:10045604-10046982 " +
                        "--"+ StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME + " " + IntervalFilteringVcfWriter.Mode.STARTS_IN, testFile) ,
                Collections.singletonList(getToolTestDataDir() + "expected/testSelectVariants_onlyOutputInInterval.vcf")
        );
        spec.executeTest("test that variant interval output filter mode is applied", this);
    }
}
