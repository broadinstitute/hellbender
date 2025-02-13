package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ExampleVariantLocusWalkerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleVariantLocusWalker() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -L 1" +
                " -R " + hg19MiniReference +
                " -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                " -O %s",
            Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleVariantLocusWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleVariantLocusWalker", this);
    }

}
