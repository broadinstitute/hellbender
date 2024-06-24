package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;

import java.io.File;
import java.io.IOException;

public class FlowFeatureMapperIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir + FlowTestConstants.FEATURE_MAPPING_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "RG",
                "--copy-attr", "AS,Integer,AS attribute, as copied",
                "--copy-attr", "rq,Float",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755"
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
    public void testSurroundingMedianQuality() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_smq_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_smq_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "tr",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755",
                "--surrounding-median-quality-size", "1000" // use a high value to stress the code
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
    public void testSurroundingAllQuality() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_smq_all_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_smq_all_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "tr",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755",
                "--surrounding-median-quality-size", "1000",
                "--surrounding-mean-quality-size", "1000"
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
    public void testQCFailed() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_qc_failed_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_qc_failed_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "tr",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755",
                "--include-qc-failed-reads", "true"
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
    public void testReportAllAlts() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_report_all_alts_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_report_all_alts_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "RG",
                "--copy-attr", "AS,Integer,AS attribute, as copied",
                "--copy-attr", "rq,Float",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755",
                "--report-all-alts",
                "-L", "chr20:1099787-1101000"
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
    public void testTagBasesWithAdjacentRefDiff() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_tag_adjs_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_tag_adjs_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "RG",
                "--copy-attr", "AS,Integer,AS attribute, as copied",
                "--copy-attr", "rq,Float",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755",
                "--tag-bases-with-adjacent-ref-diff",
                "-L", "chr20:1099787-1101000"
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
