package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CompareReferencesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testCompareReferencesRenamedSequence() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File output = createTempFile("testCompareReferencesRenamedSequence", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesRenamedSequence.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMultipleReferences() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesMultipleReferences", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesMultipleReferences.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMissingValue() throws IOException{
        // SNP on chr2 causes the sequence to have a different MD5 than the hg19mini sequence of the same name
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesMissingValue", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @DataProvider(name = "md5ArgumentData")
    public Object[][] md5ArgumentData() {
        return new Object[][]{
                // ref1, ref2, expected output, md5 calculation mode
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table"),
                        "USE_DICT"
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesRecalculateIfMissing.table"),
                        "RECALCULATE_IF_MISSING"
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesAlwaysRecalculateMD5.table"),
                        "ALWAYS_RECALCULATE"
                },
        };
    }

    @Test(dataProvider = "md5ArgumentData")
    public void testCompareReferencesMD5CalculationModes(File ref1, File ref2, File expected, String md5Mode) throws IOException {
        final File output = createTempFile("testCompareReferencesMD5CalculationModes", ".table");
        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", md5Mode};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCompareReferencesUseDictMD5MissingValue() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesUseDictMD5", ".table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", "USE_DICT"};
        runCommandLine(args);
    }

    @DataProvider(name = "cloudInputData")
    public Object[][] cloudInputData() {
        String testBucket = getGCPTestInputPath() + "org/broadinstitute/hellbender/tools/reference/";
        return new Object[][]{
                // list of refs, expected output
                new Object[]{ Arrays.asList(testBucket + "hg19mini.fasta",
                        testBucket + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table")
                },
                new Object[]{ Arrays.asList(testBucket + "hg19mini.fasta",
                        testBucket + "hg19mini_1renamed.fasta",
                        testBucket + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesMultipleReferences.table")
                },
        };
    }

    @Test(dataProvider = "cloudInputData", groups = "bucket")
    public void testCompareReferencesCloudInputs(List<String> fastas, File expected) throws IOException{
        final File output = createTempFile("testCompareReferencesCloudInputs", ".table");
        List<String> arguments = new ArrayList<>();
        for(int i = 0; i < fastas.size(); i++){
            if(i == 0){
                arguments.add("-R");

            }
            else{
                arguments.add("-refcomp");
            }
            arguments.add(fastas.get(i));
        }
        arguments.add("-O");
        arguments.add(output.getAbsolutePath());

        final String[] args = arguments.toArray(new String[0]);
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    // no assertions made, testing tool runs successfully without -O argument
    @Test
    public void testCompareReferencesToStdOutput() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath()};
        runCommandLine(args);
    }

    // The following test runs the tool on different combinations of reference files
    // and produce output to stdout for the sake of manually inspecting outputs.
    // Disabled, as no actual assertions made.
    @Test(enabled = false)
    public void testCompareReferencesMultipleReferencesStdOut() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File ref4 = new File(getToolTestDataDir() + "hg19mini_missingchr1.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref3.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath(),
                "-refcomp", ref4.getAbsolutePath(), "-display-sequences-by-name", "-display-only-differing-sequences"};
        runCommandLine(args);
    }

    // FIND_SNPS_ONLY tests:
    @DataProvider(name = "findSNPsData")
    public Object[][] findSNPsData() {
        return new Object[][]{
                // ref1, ref2, output name, expected output
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2multiplesnps.fasta"),
                        "hg19mini.fasta_hg19mini_chr2multiplesnps.fasta_snps.tsv",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr2multiplesnps.fasta_snps.tsv"),
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2iupacsnps.fasta"),
                        "hg19mini.fasta_hg19mini_chr2iupacsnps.fasta_snps.tsv",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr2iupacsnps.fasta_snps.tsv"),
                },
                // produces empty output file bc FIND_SNPS_ONLY doesn't work on indels
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr1indel.fasta"),
                        "hg19mini.fasta_hg19mini_chr1indel.fasta_snps.tsv",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr1indel.fasta_snps.tsv"),
                },
        };
    }
    @Test(dataProvider = "findSNPsData")
    public void testFindSNPs(File ref1, File ref2, String actualOutputFileName, File expectedOutput) throws IOException{
        final File output = IOUtils.createTempDir("tempFindSNPs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FIND_SNPS_ONLY", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File actualOutput = new File(output, actualOutputFileName);
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @DataProvider(name = "baseComparisonIllegalArgumentsData")
    public Object[][] baseComparisonIllegalArgumentsData() {
        return new Object[][]{
                // base comparison mode
                new Object[]{"FULL_ALIGNMENT"},
                new Object[]{"FIND_SNPS_ONLY"},
        };
    }

    // tool throws a UserException.CouldNotCreateOutputFile but gets wrapped in a NullPointerException
    @Test(dataProvider = "baseComparisonIllegalArgumentsData", expectedExceptions = NullPointerException.class)
    public void testFindSNPsNoOutputDir(String baseComparisonMode) throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2iupacsnps.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", baseComparisonMode};
        runCommandLine(args);
    }

    @Test(dataProvider = "baseComparisonIllegalArgumentsData", expectedExceptions = UserException.CouldNotCreateOutputFile.class)
    public void testFindSNPsNonexistentOutputDir(String baseComparisonMode) throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2iupacsnps.fasta");
        final File nonexistentOutputDir = new File("/tempreferences");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", baseComparisonMode, "-base-comparison-output", nonexistentOutputDir.getAbsolutePath()};
        runCommandLine(args);
    }

    @Test(dataProvider = "baseComparisonIllegalArgumentsData", expectedExceptions = UserException.BadInput.class)
    public void testFindSNPsMoreThanTwoReferences(String baseComparisonMode) throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2iupacsnps.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_chr2multiplesnps.fasta");
        final File output = IOUtils.createTempDir("tempFindSNPs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath(), "-base-comparison", baseComparisonMode, "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);
    }

    // FULL_ALIGNMENT tests:
    @DataProvider(name = "fullAlignmentData")
    public Object[][] fullAlignmentData() {
        return new Object[][]{
                // ref1, ref2, output name, expected output
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2multiplesnps.fasta"),
                        "hg19mini.fasta_hg19mini_chr2multiplesnps.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr2multiplesnps.fasta.vcf"),
                },
                // order that these 2 references are specified determines whether it is an insertion or deletion
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr1indel.fasta"),
                        "hg19mini.fasta_hg19mini_chr1indel.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.testDeletion.hg19mini.fasta_hg19mini_chr1indel.fasta.vcf"),
                },
                // order that these 2 references are specified determines whether it is an insertion or deletion
                new Object[]{ new File(getToolTestDataDir() + "hg19mini_chr1indel.fasta"),
                        new File(getToolTestDataDir() + "hg19mini.fasta"),
                        "hg19mini_chr1indel.fasta_hg19mini.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.testInsertion.hg19mini_chr1indel.fasta_hg19mini.fasta.vcf"),
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_snpsmultiplecontigs.fasta"),
                        "hg19mini.fasta_hg19mini_snpsmultiplecontigs.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_snpsmultiplecontigs.fasta.vcf"),
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_snpandindel.fasta"),
                        "hg19mini.fasta_hg19mini_snpandindel.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.SNPandINDEL.hg19mini.fasta_hg19mini_snpandindel.fasta.vcf"),
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_manySNPsINDELs.fasta"),
                        "hg19mini.fasta_hg19mini_manySNPsINDELs.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_manySNPsINDELs.fasta.vcf"),
                },
        };
    }

    @Test(dataProvider = "fullAlignmentData")
    public void testFullAlignmentMode(File ref1, File ref2, String actualOutputFileName, File expectedOutput) throws IOException{
        final File output = IOUtils.createTempDir("tempDirFullAlignment");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File actualOutput = new File(output, actualOutputFileName);
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @DataProvider(name = "findSNPsEquivalentSequencesData")
    public Object[][] findSNPsEquivalentSequencesData() {
        return new Object[][]{
                // ref1, ref2, equivalency file, output name, expected output
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_1renamed_snp1and2.fasta"),
                        new File(getToolTestDataDir() + "equivalent_sequences.tsv"),
                        "hg19mini.fasta_hg19mini_1renamed_snp1and2.fasta_snps.tsv",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_1renamed_snp1and2.fasta_snps.tsv"),
                }
        };
    }
    @Test(dataProvider = "findSNPsEquivalentSequencesData")
    public void testFindSNPsEquivalentSequences(File ref1, File ref2, File equivalentSequences, String actualOutputFileName, File expectedOutput) throws IOException{
        final File output = IOUtils.createTempDir("tempFindSNPs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FIND_SNPS_ONLY", "-base-comparison-output", output.getAbsolutePath(), "-sequences-to-align", equivalentSequences.getAbsolutePath()};
        runCommandLine(args);

        final File actualOutput = new File(output, actualOutputFileName);
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @DataProvider(name = "fullAlignmentEquivalentSequencesData")
    public Object[][] findFullAlignmentEquivalentSequencesData() {
        return new Object[][]{
                // ref1, ref2, equivalency file, output name, expected output
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_1renamed_snp1and2.fasta"),
                        new File(getToolTestDataDir() + "equivalent_sequences.tsv"),
                        "hg19mini.fasta_hg19mini_1renamed_snp1and2.fasta.vcf",
                        new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_1renamed_snp1and2.fasta.vcf"),
                }
        };
    }
    @Test(dataProvider = "fullAlignmentEquivalentSequencesData")
    public void testFullAlignmentEquivalentSequences(File ref1, File ref2, File equivalentSequences, String actualOutputFileName, File expectedOutput) throws IOException{
        final File output = IOUtils.createTempDir("tempFullAlignment");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath(), "-sequences-to-align", equivalentSequences.getAbsolutePath()};
        runCommandLine(args);

        final File actualOutput = new File(output, actualOutputFileName);
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testCasePreservation() throws IOException{
        File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        File ref2 = new File(getToolTestDataDir() + "hg19mini_snp1_casechange1.fasta");
        final File output = IOUtils.createTempDir("tempIgnoreCaseDifferences");
        String actualOutputFileName =  "hg19mini.fasta_hg19mini_snp1_casechange1.fasta_snps.tsv";
        File expectedOutput = new File(getToolTestDataDir() + "expected.case_sensitive_hg19mini.fasta_hg19mini_snp1_casechange1.fasta_snps.tsv");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FIND_SNPS_ONLY", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File actualOutput = new File(output, actualOutputFileName);
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testIgnoreCaseDifferences() throws IOException{
        File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        File ref2 = new File(getToolTestDataDir() + "hg19mini_snp1_casechange1.fasta");
        String actualOutputFileName =  "hg19mini.fasta_hg19mini_snp1_casechange1.fasta_snps.tsv";
        File expectedOutput = new File(getToolTestDataDir() + "expected.ignore_case_differences_hg19mini.fasta_hg19mini_snp1_casechange1.fasta_snps.tsv");
        final File output = IOUtils.createTempDir("tempIgnoreCaseDifferences");


        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FIND_SNPS_ONLY", "-base-comparison-output", output.getAbsolutePath(), "-ignore-case-level-differences"};
        runCommandLine(args);

        final File actualOutput = new File(output, actualOutputFileName);
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

}
