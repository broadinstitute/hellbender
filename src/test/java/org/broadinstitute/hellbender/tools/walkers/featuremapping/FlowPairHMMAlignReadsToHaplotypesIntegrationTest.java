package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import org.broadinstitute.hellbender.CommandLineProgramTest;

public class FlowPairHMMAlignReadsToHaplotypesIntegrationTest extends CommandLineProgramTest{
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static final String testDir = publicTestDir + FlowTestConstants.FEATURE_MAPPING_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testExpanded() throws IOException {

        final File outputDir = createTempDir("testFlowAlignReads");
        final File expectedFile = new File(testDir + "/read_to_hap_expanded.txt");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/read_to_hap_expanded.txt");

        final String[] args = new String[] {
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/alignReadsToHaplotypesTest.bam",
                "-H", testDir + "/alignReadsToHaplotypesTest.fa",
                "--flow-use-t0-tag",
                "-E", "FlowBased",
                "--flow-fill-empty-bins-value", "0.00001",
                "--flow-likelihood-optimized-comp"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "#");
        }
    }

    @Test
    public void testConcise() throws IOException {

        final File outputDir = createTempDir("testFlowAlignReads");
        final File expectedFile = new File(testDir + "/read_to_hap_concise.txt");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/read_to_hap_concise.txt");

        final String[] args = new String[] {
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/alignReadsToHaplotypesTest.bam",
                "-H", testDir + "/alignReadsToHaplotypesTest.fa",
                "--ref-haplotype", "Hap_2",
                "--concise-output-format",
                "--flow-use-t0-tag",
                "-E", "FlowBased",
                "--flow-fill-empty-bins-value", "0.00001",
                "--flow-likelihood-optimized-comp"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "#");
        }
    }
}