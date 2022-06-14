package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

/**
 * These tests are scaled down versions of GATK3 tests (the reduction coming from smaller query intervals),
 * and the expected results files were generated by running GATK on the same files with the same arguments.
 * In most cases, the GATK3 runs were done using modifications of the integration tests, as opposed to running
 * from the command line, in order to ensure that the initial state of the random number generator is the same
 * between versions. In all cases, the variants in the expected files are identical to those produced by GATK3,
 * though the VCF headers were hand modified to account for small differences in the metadata lines.
 *
 * UPDATE: The expected results from GATK3 were updated in https://github.com/broadinstitute/gatk/pull/7709.
 * However, we left the original comments referencing GATK3 below untouched.
 */
public class VariantRecalibratorIntegrationTest extends CommandLineProgramTest {

    private final String[] VQSRSNPParamsWithResources =
        new String[] {
            "--variant",
            getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf",
            "-L","20:1,000,000-10,000,000",
            "--resource:known,known=true,prior=10.0",
            getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf",
            "--resource:truth_training1,truth=true,training=true,prior=15.0",
            getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf",
            "--resource:truth_training2,training=true,truth=true,prior=12.0",
            getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf",
            "-an", "QD", "-an", "HaplotypeScore", "-an", "HRun",
            "--trust-all-polymorphic", // for speed
            "-mode", "SNP",
            "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

    private final String[] VQSRSNPsWithAnnotationDupe =
            new String[] {
                    "--variant",
                    getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf",
                    "-L","20:1,000,000-10,000,000",
                    "--resource:known,known=true,prior=10.0",
                    getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf",
                    "--resource:truth_training1,truth=true,training=true,prior=15.0",
                    getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf",
                    "--resource:truth_training2,training=true,truth=true,prior=12.0",
                    getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf",
                    "-an", "QD", "-an", "HaplotypeScore", "-an", "HRun", "-an", "HRun",
                    "--trust-all-polymorphic", // for speed
                    "-mode", "SNP",
                    "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
            };

    private final String[] VQSRBothParamsWithResources =
            new String[] {
                    "--variant",
                    getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf",
                    "--variant",
                    getLargeVQSRTestDataDir() + "g94982_20_1m_10m_python_2dcnn.indels.vcf.gz",
                    "-L","20:1,000,000-10,000,000",
                    "--resource:known,known=true,prior=10.0",
                    getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf",
                    "--resource:truth_training1,truth=true,training=true,prior=15.0",
                    getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf",
                    "--resource:truth_training2,training=true,truth=true,prior=12.0",
                    getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf",
                    "-an", "QD", "-an", "HaplotypeScore", "-an", "HRun",
                    "--trust-all-polymorphic", // for speed
                    "-mode", "BOTH",
                    "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                    "--max-gaussians", "6"
            };

    private final String[] VQSRBothAggregateParamsWithResources =
            new String[] {
                    "--variant",
                    getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf",
                    "--aggregate",
                    getLargeVQSRTestDataDir() + "g94982_20_1m_10m_python_2dcnn.indels.vcf.gz",
                    "-L","20:1,000,000-10,000,000",
                    "--resource:known,known=true,prior=10.0",
                    getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf",
                    "--resource:truth_training1,truth=true,training=true,prior=15.0",
                    getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf",
                    "--resource:truth_training2,training=true,truth=true,prior=12.0",
                    getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf",
                    "-an", "QD", "-an", "HaplotypeScore", "-an", "HRun",
                    "--trust-all-polymorphic", // for speed
                    "-mode", "BOTH",
                    "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                    " --output %s" + " -tranches-file %s",
                    "--max-gaussians", "6"
            };

    private final String[] alleleSpecificVQSRParams =
        new String[] {
            "--variant",
            getLargeVQSRTestDataDir() + "chr1snippet.doctoredMQ.sites_only.vcf.gz",
            "-L","chr1:1-10,000,000",
            "-resource:same,known=false,training=true,truth=true,prior=15",
            getLargeVQSRTestDataDir() + "chr1snippet.doctoredMQ.sites_only.vcf.gz",
            "-an", "AS_QD", "-an", "AS_ReadPosRankSum", "-an", "AS_MQ", "-an", "AS_SOR", //AS_MQRankSum has zero variance and AS_FS is nearly constant; also different annotation orders may not converge
            "--trust-all-polymorphic", // for speed
            "--use-allele-specific-annotations",
            "-mode", "SNP",
            "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

    private final String[] alleleSpecificVQSRParamsTooManyAnnotations =
            new String[] {
                    "--variant",
                    getLargeVQSRTestDataDir() + "chr1snippet.doctoredMQ.sites_only.vcf.gz",
                    "-L","chr1:1-10,000,000",
                    "-resource:same,known=false,training=true,truth=true,prior=15",
                    getLargeVQSRTestDataDir() + "chr1snippet.doctoredMQ.sites_only.vcf.gz",
                    "-an", "AS_QD", "-an", "AS_ReadPosRankSum", "-an", "AS_MQ", "-an", "AS_SOR", "-an", "AS_FS", "-an", "AS_MQRankSum", //AS_MQRankSum has zero variance and AS_FS is nearly constant; also different annotation orders may not converge
                    "--trust-all-polymorphic", // for speed
                    "--use-allele-specific-annotations",
                    "-mode", "SNP",
                    "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
            };

    private final String[] alleleSpecificVQSRParamsNoNegativeData =
            new String[] {
                    "--variant",
                    getLargeVQSRTestDataDir() + "chr1snippet.doctoredMQ.sites_only.vcf.gz",
                    "-L","chr1:1-10,000,000", "--bad-lod-score-cutoff", "-5.5",
                    "-resource:same,known=false,training=true,truth=true,prior=15",
                    getLargeVQSRTestDataDir() + "chr1snippet.doctoredMQ.sites_only.vcf.gz",
                    "-an", "AS_QD", "-an", "AS_ReadPosRankSum", "-an", "AS_MQ", "-an", "AS_SOR", //AS_MQRankSum has zero variance and AS_FS is nearly constant; also different annotation orders may not converge
                    "--trust-all-polymorphic", // for speed
                    "--use-allele-specific-annotations",
                    "-mode", "SNP",
                    "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
            };

    @Override
    public String getToolTestDataDir(){
        return toolsTestDir + "walkers/VQSR/";
    }

    private String getLargeVQSRTestDataDir(){
        return largeFileTestDir + "VQSR/";
    }

    @BeforeMethod
    public void initializeVariantRecalTests() {
        //Reset the RNG in order to align it with the initial state of the GATK3 RNG at the time
        //the tests start running, which we want to do in order to get the same results produced by
        //GATK3. Note that this means the results of running these tests will be different if they
        //are run manually outside of the test framework.
        logger.info("Initializing VQSR tests/resetting random number generator");
        Utils.resetRandomGenerator();
    }

    @DataProvider(name="VarRecalSNP")
    public Object[][] getVarRecalSNPData() {
        return new Object[][] {
                {VQSRSNPsWithAnnotationDupe,
                        getLargeVQSRTestDataDir() + "expected/SNPDefaultTranches.txt",
                        getLargeVQSRTestDataDir() + "snpRecal.vcf"
                },
                {
                    VQSRSNPParamsWithResources,
                getLargeVQSRTestDataDir() + "expected/SNPDefaultTranches.txt",
                getLargeVQSRTestDataDir() + "snpRecal.vcf"
            },
            {
                alleleSpecificVQSRParams,
                getToolTestDataDir() + "expected.AS.tranches",
                getLargeVQSRTestDataDir() + "expected/expected.AS.recal.vcf"
            }

        };
    }

    @DataProvider(name="VarRecalSNPAlternateTranches")
    public Object[][] getVarRecalSNPAlternateTranchesData() {
        return new Object[][] {
                {
                        VQSRSNPParamsWithResources,
                    getLargeVQSRTestDataDir() + "expected/SNPAlternateTranches.txt",
                    getLargeVQSRTestDataDir() + "snpRecal.vcf"
                },
                {
                    alleleSpecificVQSRParams,
                    getToolTestDataDir() + "expected.AS.alternate.tranches",
                    getLargeVQSRTestDataDir() + "expected/expected.AS.recal.vcf"
                }

        };
    }

    @DataProvider(name="SNPRecalCommand")
    public Object[][] getSNPRecalCommand() {
        return new Object[][] {
                {
                        VQSRBothParamsWithResources
                }
        };
    }

    private void doSNPTest(final String[] params, final String expectedTranchesFile, final String expectedRecalFile) throws IOException {
        //NOTE: The number of iterations required to ensure we have enough negative training data to proceed,
        //as well as the test results themselves, are both very sensitive to the state of the random number
        //generator at the time the tool starts to execute. Sampling a single integer from the RNG at the
        //start aligns the initial state of the random number generator with the initial state of the GATK3
        //random number generator at the time the tests start executing (both RNGs use the same seed, but
        //the pathway to the test in GATK3 results in an additional integer being sampled), thus allowing the
        //results to be easily compared against those produced by GATK3. This also happens to allow this
        //test to succeed on the first iteration (max_attempts=1) which we want for test performance reasons.
        //Failing to call nextInt here would result in the model failing on the first 3 attempts (it would
        //succeed on the 4th), and the results would not match those produced by GATK3 on these same inputs.
        //
        //Also note that due to this RNG conditioning, this test will produce different results when manually
        //outside of the test framework.
        @SuppressWarnings("unused")
        final int hack = Utils.getRandomGenerator().nextInt();

        // use an ArrayList - ArgumentBuilder tokenizes using the "=" in the resource args
        List<String> args = new ArrayList<>(params.length);
        Stream.of(params).forEach(arg -> args.add(arg));

        File recalOut = createTempFile("testVarRecalSnp", ".vcf");
        File tranchesOut = createTempFile("testVarRecalSnp", ".txt");
        args.addAll(addTempFileArgs(recalOut, tranchesOut));

        final VariantRecalibrator varRecalTool = new VariantRecalibrator();
        Assert.assertEquals(varRecalTool.instanceMain(args.toArray(new String[args.size()])), true);

        // the expected vcf is not in the expected dir because its used
        // as input for the ApplyVQSR test
        IntegrationTestSpec.assertEqualTextFiles(recalOut, new File(expectedRecalFile));
        IntegrationTestSpec.assertEqualTextFiles(tranchesOut, new File(expectedTranchesFile));
    }

    @Test(dataProvider = "VarRecalSNP")
    public void testVariantRecalibratorSNP(final String[] params, final String tranchesPath, final String recalPath) throws IOException {
        doSNPTest(params, tranchesPath, recalPath);
    }

    @Test(dataProvider = "VarRecalSNPAlternateTranches")
    public void testVariantRecalibratorSNPAlternateTranches(final String[] params, final String tranchesPath, final String recalPath) throws IOException {
        // same as testVariantRecalibratorSNP but with specific tranches
        List<String> args = new ArrayList<>(params.length);
        Stream.of(params).forEach(arg -> args.add(arg));
        args.addAll(
                Arrays.asList(
                        "-tranche", "100.0",
                        "-tranche", "99.95",
                        "-tranche", "99.9",
                        "-tranche", "99.8",
                        "-tranche", "99.6",
                        "-tranche", "99.5",
                        "-tranche", "99.4",
                        "-tranche", "99.3",
                        "-tranche", "99.0",
                        "-tranche", "98.0",
                        "-tranche", "97.0",
                        "-tranche", "90.0"
                )
        );
        doSNPTest(args.toArray(new String[args.size()]), tranchesPath, recalPath);
    }

    @Test(dataProvider = "SNPRecalCommand")
    public void testVariantRecalibratorSNPMaxAttempts(final String[] params) throws IOException {
        // For this test, we deliberately *DON'T* sample a single random int as above; this causes
        // the tool to require 4 attempts to acquire enough negative training data to succeed

        // use an ArrayList - ArgumentBuilder tokenizes using the "=" in the resource args
        List<String> args = new ArrayList<>(params.length);
        Stream.of(params).forEach(arg -> args.add(arg));
        File recalOut = createTempFile("testVarRecalMaxAttempts", ".vcf");
        File tranchesOut = createTempFile("testVarRecalMaxAttempts", ".txt");
        args.addAll(addTempFileArgs(recalOut, tranchesOut));

        args.add("--max-attempts");
        args.add("4"); // it takes for for this test to wind up with enough training data

        runCommandLine(args);
        final VariantRecalibrator varRecalTool = new VariantRecalibrator();
        Assert.assertEquals(varRecalTool.instanceMain(args.toArray(new String[args.size()])), true);
        Assert.assertEquals(varRecalTool.max_attempts, 4);
    }

    private List<String> addTempFileArgs(final File recalOutFile, final File tranchesOutFile) {
        List<java.lang.String> args = new ArrayList<>(2);
        args.add("--output");
        args.add(recalOutFile.getAbsolutePath());
        args.add("--tranches-file");
        args.add(tranchesOutFile.getAbsolutePath());
        return args;
    }

    @Test
    public void testVariantRecalibratorIndel() throws IOException {
        @SuppressWarnings("unused")
        final int hack = Utils.getRandomGenerator().nextInt();
        final String inputFile = getLargeVQSRTestDataDir() + "combined.phase1.chr20.raw.indels.filtered.sites.1M-10M.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --resource:known,known=true,prior=10.0 " + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf" +
                " --resource:truth_training,training=true,truth=true,prior=15.0 " + getLargeVQSRTestDataDir() + "ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf" +
                " --variant " + inputFile +
                " -L 20:1,000,000-10,000,000" +
                " -an QD -an ReadPosRankSum -an HaplotypeScore" +
                " -mode INDEL -max-gaussians 3" +
                " --trust-all-polymorphic" + // for speed
                " --output %s" +
                " -tranches-file %s" +
                " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Arrays.asList(
                        // the "expected" vcf is not in the expected dir because its used
                        // as input for the ApplyVQSR test
                        getLargeVQSRTestDataDir() + "indelRecal.vcf",
                        getLargeVQSRTestDataDir() + "expected/indelTranches.txt"));
        spec.executeTest("testVariantRecalibratorIndel"+  inputFile, this);
    }

    @Test
    public void testBothRecalMode() throws IOException {
        final String args = StringUtils.join(VQSRBothParamsWithResources, " ");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                args + " --output %s -tranches-file %s",
                Arrays.asList(
                // the "expected" vcf is not in the expected dir because it [should be] used
                // as input for a ApplyVQSR test
                getLargeVQSRTestDataDir() + "bothRecal.vcf",
                getLargeVQSRTestDataDir() + "expected/bothTranches.txt"));
        spec.executeTest("testBothRecalMode", this);
    }

    @Test
    public void testBothAggregateRecalMode() throws IOException {
        final String args = StringUtils.join(VQSRBothAggregateParamsWithResources, " ");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                args, Arrays.asList(
                // the "expected" vcf is not in the expected dir because it [should be] used
                // as input for a ApplyVQSR test
                getLargeVQSRTestDataDir() + "bothRecalWithAggregate.vcf",
                getLargeVQSRTestDataDir() + "expected/bothTranchesWithAggregate.txt"));
        spec.executeTest("testBothRecalMode", this);
    }

    private final String tmpDir = createTempDir(this.getTestedClassName()).getAbsolutePath();
    private final String modelReportFilename = tmpDir + "/snpSampledModel.report";
    private final String modelReportRecal = getLargeVQSRTestDataDir() + "expected/snpSampledRecal.vcf";
    private final String modelReportTranches = getLargeVQSRTestDataDir() + "expected/snpSampledTranches.txt";

    @Test
    public void testVariantRecalibratorSampling() throws IOException {
        final String inputFile = getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + inputFile +
                " -L 20:1,000,000-10,000,000" +
                " --resource:known,known=true,prior=10.0 " + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf" +
                " --resource:truth_training1,truth=true,training=true,prior=15.0 " + getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf" +
                " --resource:truth_training2,training=true,truth=true,prior=12.0 " + getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf" +
                " -an QD -an HaplotypeScore -an HRun" +
                " --trust-all-polymorphic" + // for speed
                " --output %s" +
                " -tranches-file %s" +
                " --output-model " + modelReportFilename +
                " -mode SNP --max-gaussians 3" +  //reduce max gaussians so we have negative training data with the sampled input
                " -sample-every 2" +
                " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Arrays.asList(
                        modelReportRecal,
                        modelReportTranches));
        spec.executeTest("testVariantRecalibratorSampling"+  inputFile, this);
    }

    // Expected exception is a UserException but gets wrapped in a RuntimeException
    @Test(expectedExceptions = RuntimeException.class)
    public void testVariantRecalibratorFailedRscriptOutput() throws IOException {
        final String inputFile = getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + inputFile +
                        " -L 20:1,000,000-10,000,000" +
                        " --resource:known,known=true,prior=10.0 " + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf" +
                        " --resource:truth_training1,truth=true,training=true,prior=15.0 " + getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf" +
                        " --resource:truth_training2,training=true,truth=true,prior=12.0 " + getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf" +
                        " -an QD -an HaplotypeScore -an HRun" +
                        " --trust-all-polymorphic" + // for speed
                        " --output %s" +
                        " -tranches-file %s" +
                        " --output-model " + modelReportFilename +
                        " -mode SNP --max-gaussians 3" +  //reduce max gaussians so we have negative training data with the sampled input
                        " -sample-every 2" +
                        " --rscript-file " + createTempFile("rscriptOutput", ".R") +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Arrays.asList(
                        modelReportRecal,
                        modelReportTranches));
        spec.executeTest("testVariantRecalibratorFailedRscriptOutput"+  inputFile, this);
    }


    @Test
    public void testVariantRecalibratorRScriptOutput() throws IOException{
        final String inputFile = getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf";

        File unrunRscript = createTempFile("rscriptOutput", ".R");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + inputFile +
                        " -L 20:1,000,000-10,000,000" +
                        " --resource:known,known=true,prior=10.0 " + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf" +
                        " --resource:truth_training1,truth=true,training=true,prior=15.0 " + getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf" +
                        " --resource:truth_training2,training=true,truth=true,prior=12.0 " + getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf" +
                        " -an QD -an HaplotypeScore -an HRun" +
                        " --trust-all-polymorphic" + // for speed
                        " --output %s" +
                        " -tranches-file %s" +
                        " --output-model " + modelReportFilename +
                        " -mode SNP --max-gaussians 3" +  //reduce max gaussians so we have negative training data with the sampled input
                        " -sample-every 2" +
                        " --disable-rscriptexecutor " +
                        " --rscript-file " + unrunRscript +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Arrays.asList(
                        modelReportRecal,
                        modelReportTranches));
        spec.executeTest("testVariantRecalibratorRscriptOutput"+  inputFile, this);
        Assert.assertTrue(Files.exists(IOUtils.fileToPath(unrunRscript)), "Rscript file was not generated.");
        Assert.assertTrue(Files.size(IOUtils.fileToPath(unrunRscript))>0, "Rscript file was empty.");
    }



    @Test(dependsOnMethods = {"testVariantRecalibratorSampling"})
    public void testVariantRecalibratorModelInput() throws IOException {
        final String inputFile = getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + inputFile +
                        " -L 20:1,000,000-10,000,000" +
                        " --resource:known,known=true,prior=10.0 " + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf" +
                        " --resource:truth_training1,truth=true,training=true,prior=15.0 " + getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf" +
                        " --resource:truth_training2,training=true,truth=true,prior=12.0 " + getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf" +
                        " -an QD -an HaplotypeScore -an HRun" +
                        " --trust-all-polymorphic" + // for speed
                        " --output %s" +
                        " -tranches-file %s" +
                        " --input-model " + modelReportFilename +
                        " -mode SNP -max-gaussians 3" +  //reduce max gaussians so we have negative training data with the sampled input
                        " -sample-every 2" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Arrays.asList(
                        modelReportRecal,
                        modelReportTranches));
        spec.executeTest("testVariantRecalibratorModelInput"+  inputFile, this);
    }

    private final String annoOrderRecal = getLargeVQSRTestDataDir() + "expected/anno_order.recal";
    private final String annoOrderTranches = getLargeVQSRTestDataDir() + "expected/anno_order.tranches";
    private final String exacModelReportFilename = publicTestDir + "/subsetExAC.snps_model.report";

    @Test
    public void testVQSRAnnotationOrder() throws IOException {
        final String inputFile = publicTestDir + "/oneSNP.vcf";

        // We don't actually need resources because we are using a serialized model,
        // so we just pass input as resource to prevent a crash
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + inputFile +
                        " -L 1:110201699" +
                        " --resource:hapmap,known=false,training=true,truth=true,prior=15 " + inputFile +
                        " -an FS -an ReadPosRankSum -an MQ -an MQRankSum -an QD -an SOR" +
                        " --output %s" +
                        " -tranches-file %s" +
                        " --input-model " + exacModelReportFilename +
                        " --add-output-vcf-command-line false" +
                        " -ignore-all-filters -mode SNP",
                Arrays.asList(
                        annoOrderRecal,
                        annoOrderTranches));
        spec.executeTest("testVariantRecalibratorModelInput"+  inputFile, this);

        Utils.resetRandomGenerator();
        // Change annotation order and assert consistent outputs
        final IntegrationTestSpec spec2 = new IntegrationTestSpec(
                " --variant " + inputFile +
                        " -L 1:110201699" +
                        " --resource:hapmap,known=false,training=true,truth=true,prior=15 " + inputFile +
                        " -an ReadPosRankSum -an MQ -an MQRankSum -an QD -an SOR -an FS" +
                        " --output %s" +
                        " -tranches-file %s" +
                        " --input-model " + exacModelReportFilename +
                        " --add-output-vcf-command-line false" +
                        " -ignore-all-filters -mode SNP",
                Arrays.asList(
                        annoOrderRecal,
                        annoOrderTranches));
        spec2.executeTest("testVariantRecalibratorModelInput"+  inputFile, this);

    }


    @DataProvider(name="VarRecalSNPScattered")
    public Object[][] getVarRecalSNPScatteredData() {
        return new Object[][] {
                {
                        new String[] {
                                "--variant",
                                getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf",
                                "-L","20:1,000,000-10,000,000",
                                "--resource:known,known=true,prior=10.0",
                                getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf",
                                "--resource:truth_training1,truth=true,training=true,prior=15.0",
                                getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf",
                                "--resource:truth_training2,training=true,truth=true,prior=12.0",
                                getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf",
                                "-an", "QD", "-an", "HaplotypeScore", "-an", "HRun",
                                "-trust-all-polymorphic", // for speed
                                "-mode", "SNP",
                                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                                "--output-tranches-for-scatter",
                                "--vqslod-tranche", "10.0",
                                "--vqslod-tranche", "8.0",
                                "--vqslod-tranche", "6.0",
                                "--vqslod-tranche", "4.0",
                                "--vqslod-tranche", "2.0",
                                "--vqslod-tranche", "0.0",
                                "--vqslod-tranche", "-2.0",
                                "--vqslod-tranche", "-4.0",
                                "--vqslod-tranche", "-6.0",
                                "--vqslod-tranche", "-8.0",
                                "--vqslod-tranche", "-10.0",
                                "--vqslod-tranche", "-12.0"
                        }
                },
        };
    }

    @Test(dataProvider = "VarRecalSNPScattered")
    //the only way the recal file will match here is if we use the doSNPTest infrastructure -- as an IntegrationTestSpec it doesn't match for some reason
    public void testVariantRecalibratorSNPscattered(final String[] params) throws IOException {
        doSNPTest(params, getLargeVQSRTestDataDir() + "/snpTranches.scattered.txt", getLargeVQSRTestDataDir() + "snpRecal.vcf"); //tranches file isn't in the expected/ directory because it's input to GatherTranchesIntegrationTest
    }

    // One of the Gaussians has a covariance matrix with a determinant of zero,
    // (can be confirmed that the entries of sigma for the row and column with the index of the constant annotation are zero)
    // which leads to all +Inf LODs if we don't throw.
    //
    // UPDATE: Originally, this test checked for expectedExceptions = {UserException.VQSRPositiveModelFailure.class}.
    // However, the convergence failure was fixed in https://github.com/broadinstitute/gatk/pull/7709,
    // so we now expect the test to complete and can instead treat it as a regression test.
    @Test
    public void testAnnotationsWithNoVarianceSpecified() throws IOException {
        // use an ArrayList - ArgumentBuilder tokenizes using the "=" in the resource args
        List<String> args = new ArrayList<>(alleleSpecificVQSRParamsTooManyAnnotations.length);
        Stream.of(alleleSpecificVQSRParamsTooManyAnnotations).forEach(arg -> args.add(arg));

        File recalOut = createTempFile("testVarRecalSnp", ".vcf");
        File tranchesOut = createTempFile("testVarRecalSnp", ".txt");
        args.addAll(addTempFileArgs(recalOut, tranchesOut));

        final VariantRecalibrator varRecalTool = new VariantRecalibrator();
        Assert.assertEquals(varRecalTool.instanceMain(args.toArray(new String[args.size()])), true);
    }

    @Test(expectedExceptions = {UserException.VQSRNegativeModelFailure.class})
    public void testNoNegativeTrainingData() throws IOException {
        // use an ArrayList - ArgumentBuilder tokenizes using the "=" in the resource args
        List<String> args = new ArrayList<>(alleleSpecificVQSRParamsNoNegativeData.length);
        Stream.of(alleleSpecificVQSRParamsNoNegativeData).forEach(arg -> args.add(arg));

        File recalOut = createTempFile("testVarRecalSnp", ".vcf");
        File tranchesOut = createTempFile("testVarRecalSnp", ".txt");
        args.addAll(addTempFileArgs(recalOut, tranchesOut));

        final VariantRecalibrator varRecalTool = new VariantRecalibrator();
        Assert.assertEquals(varRecalTool.instanceMain(args.toArray(new String[args.size()])), true);
    }
}

