package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ShardingArgumentCollection;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ExampleReadSliderWalkerIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @DataProvider
    public Object[][] testWindowArguments() {
        return new Object[][] {
                // default values (overlapping windows without padding)
                {ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE, ExampleReadSliderWalker.DEFAULT_WINDOW_STEP, ExampleReadSliderWalker.DEFAULT_WINDOW_PAD},
                // non-overlapping windows (without padding)
                {ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE, ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE, 0},
                // non-overlapping windows (but padding adds overlaps)
                {ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE, ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE, ExampleReadSliderWalker.DEFAULT_WINDOW_STEP},
                // jumping windows (without padding)
                {ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE, ExampleReadSliderWalker.DEFAULT_WINDOW_SIZE * 2, 0},
        };
    }

    @Test(dataProvider = "testWindowArguments")
    public void testExampleReadSliderWalker(final int windowSize, final int stepSize, final int windowPad) throws IOException {
        final String expectedOutput = String.format("%sexpected_%s.%sw_%ss_%sp.txt",
                TEST_OUTPUT_DIRECTORY, getTestedToolName(),
                windowSize, stepSize, windowPad);
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                        String.format(" --%s %s ", ShardingArgumentCollection.WINDOW_SIZE_NAME, windowSize) +
                        String.format(" --%s %s ", ShardingArgumentCollection.WINDOW_STEP_NAME, stepSize) +
                        String.format(" --%s %s ", ShardingArgumentCollection.WINDOW_PAD_NAME, windowPad) +
                        " -L 1:200-1125" + // region with variants/reads
                        " -L 4:15951-16000" + // region with reference bases
                        " -R " + hg19MiniReference +
                        " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                        " -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                        " -O %s",
                Arrays.asList(expectedOutput)
        );
        testSpec.executeTest("testExampleReadSliderWalker", this);
    }
}