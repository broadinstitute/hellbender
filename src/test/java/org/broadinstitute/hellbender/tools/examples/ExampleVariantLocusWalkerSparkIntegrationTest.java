package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class ExampleVariantLocusWalkerSparkIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleVariantWalker() throws IOException {
        final File out = File.createTempFile("out", ".txt");
        Assert.assertTrue(out.delete(), "failed to perform necessary deletion during test setup");
        out.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-L 1");
        args.add("-V");
        args.add(TEST_DATA_DIRECTORY + "feature_data_source_test.vcf");
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(hg19MiniReference);
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleVariantLocusWalkerIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

}
