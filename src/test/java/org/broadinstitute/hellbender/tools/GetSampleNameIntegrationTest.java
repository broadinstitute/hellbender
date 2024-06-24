package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class GetSampleNameIntegrationTest extends CommandLineProgramTest {
    private static final File SINGLE_SAMPLE_BAM_FILE = new File(toolsTestDir, "valid.bam");
    private static final File BAD_MULTI_SAMPLE_BAM_FILE = new File(toolsTestDir, "multi_sample_bam_header.bam");

    @Test
    public void testBasicUsage() throws IOException {

        final File outputFile = createTempFile("get-sample-name", ".txt");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, SINGLE_SAMPLE_BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().count() == 1);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().filter(n -> n.equals("Hi,Mom!")).count() == 1);
    }

    @Test
    public void testUrlEncoding() throws IOException {

        final File outputFile = createTempFile("get-sample-name", ".txt");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, SINGLE_SAMPLE_BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-encode",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().count() == 1);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().filter(n -> n.equals("Hi%2CMom%21")).count() == 1);
    }

    @Test
    public void testMultiSampleBam() throws IOException {
        final File outputFile = createTempFile("get-sample-name-ms", ".txt");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAD_MULTI_SAMPLE_BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().count() == 3);
    }
}
