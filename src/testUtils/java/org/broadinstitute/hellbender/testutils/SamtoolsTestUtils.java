package org.broadinstitute.hellbender.testutils;

import htsjdk.beta.plugin.IOUtils;
import htsjdk.io.IOPath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.ProcessExecutor;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Test utilities for running samtools from GATK tests.
 */
public class SamtoolsTestUtils {
    private static final String SAMTOOLS_BINARY_ENV_VARIABLE = "GATK_SAMTOOLS_BIN";
    public final static String expectedSamtoolsVersion = "1.21";

    /**
     * @return true if samtools is available, otherwise false
     */
    public static boolean isSamtoolsAvailable() {
        final String binPath = getSamtoolsBin();
        final Path binFile = Paths.get(binPath);
        return Files.exists(binFile);
    }

    /**
     * @return true if a local samtools executable is available, otherwise throws a runtimeException
     */
    public static void assertSamtoolsAvailable() {
        if (!isSamtoolsAvailable()) {
            throw new RuntimeException(
                    String.format(
                            "No samtools executable can be found." +
                                    " The %s environment variable must be set to the name of the local samtools executable.",
                            SAMTOOLS_BINARY_ENV_VARIABLE));
        }
    }

    /**
     * @return the name and location of the local samtools executable as specified by the environment
     * variable SAMTOOLS_BINARY_ENV_VARIABLE, or the default value of "/usr/bin/samtools" if the environment
     * variable is not set
     */
    public static String getSamtoolsBin() {
        final String samtoolsPath = System.getenv(SAMTOOLS_BINARY_ENV_VARIABLE);
        return samtoolsPath == null ? "/usr/bin/samtools" : samtoolsPath;
    }

    /**
     * Execute a samtools command line if a local samtools executable is available see {@link #isSamtoolsAvailable()}.
     *
     * @param commandLine samtools command line string, excluding the "samtools" prefix. For example:
     *                    {@code "view -h -b my.sam -o my.bam"}
     * @return the {@link ProcessExecutor.ExitStatusAndOutput} resulting from the command execution, if
     * the command succeeds
     * @throws RuntimeException if the command fails, or if a local samtools executable is not available.
     */
    public static ProcessExecutor.ExitStatusAndOutput executeSamToolsCommand(final String commandLine) {
        assertSamtoolsAvailable();
        final String commandString = String.format("%s %s", getSamtoolsBin(), commandLine);
        final ProcessExecutor.ExitStatusAndOutput processStatus =
                ProcessExecutor.executeAndReturnInterleavedOutput(commandString);
        if (processStatus.exitStatus != 0) {
            // samtools seems to write some errors to stdout
            throw new RuntimeException(
                    String.format("Failure code %d returned from samtools command %s\n (stderr: %.500s)\n (stdout: %.500s)\n",
                            processStatus.exitStatus,
                            commandString,
                            processStatus.stderr == null ? "" : processStatus.stderr,
                            processStatus.stdout == null ? "" : processStatus.stdout));
        }
        return processStatus;
    }

    /**
     * Convert an input sam/bam/cram file to a temporary CRAM file using the samtools "view" command. The temp
     * file will be deleted when the process exits. Use {@link #isSamtoolsAvailable()} to determine if its safe
     * to use this method.
     *
     * @param inputSAMBAMCRAMFile input file to convert
     * @param referenceFile a valid reference file
     * @param commandLineOptions additional command line options (--input-fmt-option or --output-fmt-option)
     * @return a temporary file containing the samtools-generated results.
     */
    public static final IOPath convertToCRAM(
            final IOPath inputSAMBAMCRAMFile,
            final IOPath referenceFile,
            final String commandLineOptions) {
        assertSamtoolsAvailable();
        final IOPath tempCRAMPath = IOUtils.createTempPath("samtoolsTemporaryCRAM", FileExtensions.CRAM);
        tempCRAMPath.toPath().toFile().deleteOnExit();
        final String commandString = String.format("view -h -C -T %s %s %s -o %s",
                referenceFile.toPath().toAbsolutePath(),
                commandLineOptions == null ? "" : commandLineOptions,
                inputSAMBAMCRAMFile.toPath().toAbsolutePath(),
                tempCRAMPath.toPath().toAbsolutePath());
        executeSamToolsCommand(commandString);
        return tempCRAMPath;
}

    public static final void convertToCRAM(
            final IOPath inputSAMBAMCRAMPath,
            final IOPath outputPath,
            final IOPath referencePath,
            final String commandLineOptions) {
        assertSamtoolsAvailable();
        final String commandString = String.format("view -h -C -T %s %s %s -o %s",
                referencePath.toPath().toAbsolutePath(),
                commandLineOptions == null ? "" : commandLineOptions,
                inputSAMBAMCRAMPath.toPath().toAbsolutePath(),
                outputPath.toPath().toAbsolutePath());
        executeSamToolsCommand(commandString);
    }
}