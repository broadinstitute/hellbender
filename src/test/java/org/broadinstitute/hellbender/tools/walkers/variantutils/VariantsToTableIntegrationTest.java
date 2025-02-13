package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class VariantsToTableIntegrationTest extends CommandLineProgramTest {
    private String variantsToTableCmd(final String moreArgs) {
        return  " --variant " + getToolTestDataDir() + "soap_gatk_annotated.noChr_lines.vcf" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F TRANSITION -F DP -F SB -F set -F RankSumP -F refseq.functionalClass*" +
                " -O %s " + moreArgs;
    }

    private String variantsToTableMultiAllelicCmd(final String moreArgs) {
        return  " --variant " + getToolTestDataDir() + "multiallelic.vcf" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC -F AC -F AF" +
                " -O %s" + moreArgs;
    }

    private String variantsToTableCmdNoSamples(final String moreArgs) {
        return  " --variant " + getToolTestDataDir() + "vcfexample.noSamples.vcf" +
                " -O %s" + moreArgs;
    }

    @Test
    public void testInputFileFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant does_not_exist.vcf -O %s",
                1, UserException.CouldNotReadInputFile.class);
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testComplexVariantsToTable-FAIL", this);
    }

    @Test
    public void testOutputFileFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "soap_gatk_annotated.noChr_lines.vcf" +
                " -O /does_not_exists/txt.table",
                1, UserException.CouldNotCreateOutputFile.class);
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testComplexVariantsToTable-FAIL", this);
    }

    @Test
    public void testComplexVariantsToTableFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmd("--error-if-missing-data"),
                1, UserException.class);
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testComplexVariantsToTable-FAIL", this);
    }

    @Test
    public void testUnfilteredGenotypeFieldsFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -GF RD -GF FT --error-if-missing-data" +
                        " -O %s",
                1,
                UserException.class);
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testUnfilteredGenotypeFields-FAIL", this);
    }

    @Test
    public void testNoSamples() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmdNoSamples(" -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F TRANSITION -F EVENTLENGTH"),
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample.noSamples.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testNoSamples", this);
    }

    @Test
    public void testNoSamplesSoNoGenotypes() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmdNoSamples(" -GF DP"),
                1,
                UserException.class);
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testNoSamples", this);
    }

    @Test
    public void testComplexVariantsToTable() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmd(""),
                Arrays.asList(getToolTestDataDir() + "expected.soap_gatk_annotated.noChr_lines.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testComplexVariantsToTable", this);
    }

    @Test
    public void testMultiAllelicOneRecord() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableMultiAllelicCmd(""),
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testMultiAllelicOneRecord", this);
    }

    @Test
    public void testMultiAllelicSplitRecords() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableMultiAllelicCmd(" -SMA"),
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic.SMA.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testMultiAllelicSplitRecords", this);
    }

    @Test
    public void testGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -GF RD" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testGenotypeFields", this);
    }

    @Test
    public void testUnfilteredGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -GF RD -GF FT" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.FT.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testUnfilteredGenotypeFields", this);
    }

    @Test
    public void testMultiallelicGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "multiallelic_gt.vcf" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC" +
                        " -GF PL -GF AD" +
                        " -SMA" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic_gt.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testMultiallelicGenotypeFields", this);
    }

    @Test
    public void testGenotypeFieldsWithInline() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -GF RD -GF GT -GF GQ" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.GF_GT.GF_GT.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testGenotypeFieldsWithInline", this);
    }

    @Test
    public void testSplitMultiallelicFields() throws IOException {
        //missing AS INFO and FORMAT fields are handled
        //R-type FORMAT annotations work (MMQ)
        //A-type FORMAT annotations wotk (TLOD)
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "../../GenotypeGVCFs/threeSamples.2alts.vcf" +
                        " -SMA -F CHROM -F POS -F REF -F ALT -F FOO -ASF TLOD -ASGF TLOD -ASGF AD -ASGF MMQ -ASGF BAR -raw" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.threeSamples.2alts.MT.txt"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testGenotypeFieldsWithInline", this);

        //asking for allele-specific fields without splitting produces reasonable output
        final IntegrationTestSpec spec2 = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "../../GenotypeGVCFs/threeSamples.2alts.vcf" +
                        " -F CHROM -F POS -F REF -F ALT -ASGF TLOD -ASGF AD -ASGF MMQ -raw" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.threeSamples.2alts.MT.noSplit.txt"));
        spec2.setTrimWhiteSpace(false);
        spec2.executeTest("testGenotypeFieldsWithInline", this);

        //A-type INFO annotations work
        final IntegrationTestSpec spec4 = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "../../../VQSR/expected/applyIndelAlleleSpecificResult.vcf" +
                        " -SMA -F CHROM -F POS -F REF -F ALT -ASF AS_BaseQRankSum -ASGF AD -raw -ASF AS_FilterStatus" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.ASindelVQSR.txt"));
        spec4.setTrimWhiteSpace(false);
        spec4.executeTest("testGenotypeFieldsWithInline", this);
    }

    @Test
    public void testListFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample.withMLE.vcf" +
                        " -GF PL" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample.withMLE.GF_PL.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testGenotypeFields", this);
    }

    @Test
    public void testMoltenOutput() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER" +
                        " --moltenize" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.moltenize.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testMoltenOutput", this);
    }

    @Test
    public void testMoltenOutputWithGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -GF RD" +
                        " --moltenize" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.moltenize.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testMoltenOutputWithGenotypeFields", this);
    }

    @Test
    public void testMoltenOutputWithMultipleAlleles() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "multiallelic.vcf" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC -F AC -F AF" +
                        " --moltenize -SMA" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic.moltenize.SMA.table"));
        spec.setTrimWhiteSpace(false);
        spec.executeTest("testMoltenOutputWithMultipleAlleles", this);
    }

    @Test
    public void testNoFieldsSpecifiedNoSamples() throws IOException {
        final File inputFile = new File(getToolTestDataDir(), "VCFWithoutGenotypes_dbsnp_138.snippet.vcf");
        final File outputFile = createTempFile("noFieldsSpecifiedOutput", ".table");
        final File expectedFile = new File(getToolTestDataDir(), "expected.noFieldsSpecifiedNoSamples.table");

        final String[] args = new String[] {"--variant", inputFile.getAbsolutePath(),
                "-O", outputFile.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
    }

    @Test
    public void testNoFieldsSpecifiedWithSamples() throws IOException {
        final File inputFile = new File(getToolTestDataDir(), "VCFWithGenotypes_1000G.phase3.snippet.vcf");
        final File outputFile = createTempFile("noFieldsSpecifiedWithSamplesOutput", ".table");
        final File expectedFile = new File(getToolTestDataDir(), "expected.noFieldsSpecifiedWithSamples.table");

        final String[] args = new String[] {"--variant", inputFile.getAbsolutePath(),
                "-O", outputFile.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
    }

    @Test
    public void testNoFieldsSpecifiedFormatFieldInHeaderNoSamples() throws IOException {
        final File inputFile = new File(getToolTestDataDir(), "VCFWithoutGenotypesWithFormatField_dbsnp_138.snippet.vcf");
        final File outputFile = createTempFile("noFieldsSpecifiedNoSamplesOutput", ".table");
        final File expectedFile = new File(getToolTestDataDir(), "expected.noFieldsSpecifiedNoSamples.table");

        final String[] args = new String[] {"--variant", inputFile.getAbsolutePath(),
                "-O", outputFile.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
    }

    @Test
    public void testNumericGTFlag() throws IOException {
        final File inputFile = new File(getToolTestDataDir(), "VCFWithGenotypes_1000G.phase3.snippet.vcf");
        final File outputFile = createTempFile("numericGT", ".table");
        final File expectedFile = new File(getToolTestDataDir(), "expected.numericGT.table");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(inputFile)
                        .addOutput(outputFile)
                        .addFlag(VariantsToTable.NUMERIC_GT_FULLNAME);
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
    }
    
}
