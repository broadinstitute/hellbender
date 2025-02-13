package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithVariants;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

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
                " -V:TestFeatures " + FEATURE_INTEGRATION_TEST_DIRECTORY + "feature_data_source_test.vcf" +
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
}
