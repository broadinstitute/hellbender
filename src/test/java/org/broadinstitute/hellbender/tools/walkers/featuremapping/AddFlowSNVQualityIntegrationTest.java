package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;

public class AddFlowSNVQualityIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir + FlowTestConstants.ADD_FLOW_SNVQ_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testAddFlowSNVQTest");
        final String filename = "add_flow_snvq_output.sam";
        final File expectedFile = new File(testDir + "/" + filename);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + filename);

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", publicTestDir + FlowTestConstants.ADD_FLOW_SNVQ_DATA_DIR + "/add_flow_snvq_input.bam",
                "-L", "chr1:1-15000",
                "--max-phred-score", "50",
                "--debug-read-name", "30020185_2-UGAv3-182-1989782468",
                "--verbosity", "INFO"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@");
        }
    }

    @Test
    public void testBasicBQ() throws IOException {

        final File outputDir = createTempDir("testAddFlowSNVQTest");
        final String filename = "add_flow_snvq_output_bq.sam";
        final File expectedFile = new File(testDir + "/" + filename);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + filename);

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", publicTestDir + FlowTestConstants.ADD_FLOW_SNVQ_DATA_DIR + "/add_flow_snvq_input.bam",
                "-L", "chr1:1-15000",
                "--max-phred-score", "50",
                "--debug-read-name", "30020185_2-UGAv3-182-1989782468",
                "--verbosity", "INFO",
                "--output-quality-attribute", "BQ"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@");
        }
    }
}
