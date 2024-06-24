package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class GroundTruthScorerIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String OUTPUT_FILENAME = "ground_truth_scorer_output.csv";
    public static final String OUTPUT_FILENAME_MEANCALL = "ground_truth_scorer_output_meancall.csv";
    public static final String OUTPUT_FILENAME_UNCLIPPED = "ground_truth_scorer_output_unclipped.csv";
    public static final String OUTPUT_FILENAME_GENOME_PRIOR = "ground_truth_scorer_output_genome_prior.csv";
    public static final String OUTPUT_FILENAME_FILTER = "ground_truth_scorer_output_filter.csv";
    public static final String OUTPUT_FILENAME_REPORT = "ground_truth_scorer_output_report.txt";
    public static final String GT_SCORER_INPUT_BAM = "gt_scorer_input.bam";
    public static final String REF_38_FASTA = "Homo_sapiens_assembly38.fasta.gz";
    public static final String CHR9_SMALL_INTERVAL = "chr9:71000-74000";
    public static final String CHR9_ALL_INTERVAL = "chr9";

    private static String testDir = publicTestDir + FlowTestConstants.GROUND_TRUTH_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME);

        final String[] args = buildCommonArgs(outputFile, GT_SCORER_INPUT_BAM, true);

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
        }
    }

    @Test
    public void testMeanCall() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME_MEANCALL);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME_MEANCALL);

        final String[] args = ArrayUtils.addAll(buildCommonArgs(outputFile, GT_SCORER_INPUT_BAM, true),
                "--add-mean-call"
        );
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
        }
    }

    @Test
    public void testUseSoftclippedBases() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME_UNCLIPPED);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME_UNCLIPPED);

        final String[] args = ArrayUtils.addAll(buildCommonArgs(outputFile, GT_SCORER_INPUT_BAM, true),
                            "--use-softclipped-bases", "true",
                            "--features-file", testDir + "/" + "dbsnp.chr9.vcf.gz"
        );

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
        }
    }

    @Test
    public void testGenomePrior() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME_GENOME_PRIOR);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME_GENOME_PRIOR);

        final String[] args = ArrayUtils.addAll(buildCommonArgs(outputFile, GT_SCORER_INPUT_BAM, true),
                new String[] {"--genome-prior", testDir + "/" + "genome_prior.csv"});

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
        }
    }

    @Test
    public void testFilter() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME_FILTER);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME_FILTER);

        final String[] args = ArrayUtils.addAll(buildCommonArgs(outputFile, GT_SCORER_INPUT_BAM, true),
                new String[] {
                        "--features-file", testDir + "/" + "dbsnp.chr9.vcf.gz"
                });

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
        }
    }

    @Test
    public void testReports() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME);
        final File reportExpectedFile = new File(testDir + "/" + OUTPUT_FILENAME_REPORT);

        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME);
        final File reportFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? reportExpectedFile : new File(outputDir + "/" + OUTPUT_FILENAME_REPORT);

        final String[] args = ArrayUtils.addAll(buildCommonArgs(outputFile, GT_SCORER_INPUT_BAM, true),
                new String[] {
                        "--report-file", reportFile.getAbsolutePath(),
                        "--omit-zeros-from-report",
                        "--quality-percentiles", "0,20,50,90,99,100"
                });

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
            IntegrationTestSpec.assertEqualTextFiles(reportFile, reportExpectedFile);
        }
    }

    private String[] buildCommonArgs(final File outputFile, final String inputFile, boolean smallInterval) {

        return new String[] {
                "-R", largeFileTestDir + "/" + REF_38_FASTA,
                "-I", testDir + "/" + inputFile,
                "--output-csv", outputFile.getAbsolutePath(),
                "--intervals", smallInterval ? CHR9_SMALL_INTERVAL : CHR9_ALL_INTERVAL,
                "--likelihood-calculation-engine", "FlowBased"
        };

    }
}
