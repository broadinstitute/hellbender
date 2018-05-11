package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithVariants;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public final class FeatureSupportIntegrationTest extends CommandLineProgramTest {
    private static final String FEATURE_INTEGRATION_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Override
    public String getTestedClassName() {
        return ExampleReadWalkerWithVariants.class.getSimpleName();
    }

    @Test
    public void testFeatureSupportUsingVCF() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + "TestFeatures:" + FEATURE_INTEGRATION_TEST_DIRECTORY + "feature_data_source_test.vcf" +
                " --groupVariantsBySource" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeatureSupportUsingVCF_output.txt")
        );
        testSpec.executeTest("testFeatureSupportUsingVCF", this);
    }

    @Test
    public void testFeaturesAsIntervals() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/intervals_from_features_test.vcf" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeaturesAsIntervals_output.txt")
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testFeaturesAsIntervalsWithExclusion() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/intervals_from_features_test.vcf" +
                " -XL " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/intervals_from_features_test_exclude.vcf" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeaturesAsIntervalsWithExclusion_output.txt")
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testFeaturesAsIntervalsNonExistentFile() throws IOException {
        // Non-existent files should be interpreted as interval strings and fail interval parsing
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "non_existent_file.vcf" +
                " -O %s",
                1,
                UserException.MalformedGenomeLoc.class
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testFeaturesAsIntervalsUnrecognizedFormatFile() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/unrecognized_format_file.xyz" +
                " -O %s",
                1,
                UserException.CouldNotReadInputFile.class
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testSitesOnlyMode() {
        File out = createTempFile("GTStrippedOutput", "vcf");
        String[] args = new String[] {
                "-V",  FEATURE_INTEGRATION_TEST_DIRECTORY + "vcf_with_genotypes.vcf",
                "--" + StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
                "-O",
                out.getAbsolutePath()};
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList(args), SelectVariants.class.getSimpleName()));

        // Assert that the genotype field has been stripped from the file
        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(out.getAbsolutePath());

        Assert.assertFalse(results.getLeft().hasGenotypingData());
        for (VariantContext v: results.getRight()) {
            Assert.assertFalse(v.hasGenotypes());
        }
    }

    @Test
    // this test asserts that a helpful exception is thrown for blockZipped files lacking an index as they may not be fully supported
    //TODO this is a temporary fix until https://github.com/broadinstitute/gatk/issues/4224 has been resolved
    public void testUnindexedBZippedFile() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                        " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                        " -V " + toolsTestDir + "IndexFeatureFile/4featuresHG38Header.unindexed.vcf.gz" +
                        " -O %s",
                1,
                UserException.MissingIndex.class
        );
        testSpec.executeTest("testMissingIndexFeatureFile", this);
    }
}
