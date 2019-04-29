package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class ReblockGVCFIntegrationTest extends CommandLineProgramTest {

    private static final String hg38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";
    private static final String b37_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.fasta";

    @Test  //covers inputs with "MQ" annotation
    public void testJustOneSample() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-L chr20:69485-69791 -O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf -rgq-threshold 20" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testJustOneSample.expected.g.vcf"));
        spec.executeTest("testJustOneSample", this);
    }

    @Test
    public void testGVCFReblockingIsContiguous() throws Exception {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final File expected = new File(largeFileTestDir + "testProductionGVCF.expected.g.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addArgument("V", largeFileTestDir + "NA12878.prod.chr20snippet.g.vcf.gz")
                .addArgument("rgq-threshold", "20")
                .addArgument("L", "20:60001-1000000")
                .addArgument("A", "Coverage")
                .addArgument("A", "RMSMappingQuality")
                .addArgument("A", "ReadPosRankSumTest")
                .addArgument("A", "MappingQualityRankSumTest")
                .addBooleanArgument("disable-tool-default-annotations", true)
                .addOutput(output);
        runCommandLine(args);

        final CommandLineProgramTester validator = ValidateVariants.class::getSimpleName;
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.addArgument("R", b37_reference_20_21);
        args2.addArgument("V", output.getAbsolutePath());
        args2.addArgument("L", "20:60001-1000000");
        args2.add("-gvcf");
        validator.runCommandLine(args2);  //will throw a UserException if GVCF isn't contiguous

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList()));
        }
    }

    @Test  //absolute minimal output
    public void testOneSampleAsForGnomAD() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-drop-low-quals -do-qual-approx -L chr20:69485-69791 -O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " -A Coverage -A RMSMappingQuality -A ReadPosRankSumTest -A MappingQualityRankSumTest --disable-tool-default-annotations true",
                Arrays.asList(getToolTestDataDir() + "testOneSampleAsForGnomAD.expected.g.vcf"));
        spec.executeTest("testOneSampleDropLows", this);
    }

    //TODO: this isn't actually correcting non-ref GTs because I changed some args around -- separate out dropping low qual alleles and low qual sites?
    @Test  //covers non-ref AD and non-ref GT corrections
    public void testNonRefADCorrection() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "nonRefAD.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testNonRefADCorrection.expected.g.vcf"));
        spec.executeTest("testNonRefADCorrection", this);
    }

    @Test //covers inputs with "RAW_MQ" annotation
    public void testRawMQInput() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "prod.chr20snippet.withRawMQ.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "prod.chr20snippet.withRawMQ.expected.g.vcf"));
        spec.executeTest("testRawMQInput", this);
    }

    @Test
    public void testASAnnotationsAndSubsetting() throws Exception {
        //some subsetting, but never dropping the first alt
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.NA12878.AS.chr20snippet.reblocked.g.vcf"));
        spec.executeTest("testASAnnotationsAndSubsetting", this);

        //one case where first alt is dropped
        final IntegrationTestSpec spec2 = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12892.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.NA12892.AS.chr20snippet.reblocked.g.vcf"));
        spec2.executeTest("testASAnnotationsAndSubsetting2", this);
    }

    @Test
    public void testNewCompressionScheme() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60",
                Arrays.asList(getToolTestDataDir() + "expected.NA12878.AS.chr20snippet.reblocked.hiRes.g.vcf"));
        spec.executeTest("testNewCompressionScheme", this);
    }
}