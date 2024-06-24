package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

public final class SystemCommandUtilsTest extends GATKBaseTest {

    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/extract");
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    /**
     * Run the h5diff utility to compare the specified files. If h5diff exits with a non-zero
     * status code (indicating that there were differences), the differences will get printed
     * in the exception message and show up in the TestNG logs on github.
     *
     * @param expected HDF5 file with the expected results
     * @param actual HDF5 file with the actual results
     */
    static void runH5Diff(final String expected, final String actual) {
        final ProcessController controller = ProcessController.getThreadLocal();

        // -r: Report mode. Print the differences.
        // --use-system-epsilon: Return a difference if and only if the difference between two data values exceeds
        //                       the system value for epsilon.
        final String[] command = new String[] { "h5diff", "-r", "--use-system-epsilon", expected, actual };

        runProcessAndCaptureOutputInExceptionMessage(controller, command);
    }

    /**
     * Run the standard diff utility to compare the specified files. If diff exits with a non-zero
     * status code (indicating that there were differences), the differences will get printed
     * in the exception message and show up in the TestNG logs on github.
     *
     * @param expected file with the expected results
     * @param actual file with the actual results
     */
    static void runDiff(final String expected, final String actual) {
        final ProcessController controller = ProcessController.getThreadLocal();
        final String[] command = new String[] { "diff", expected, actual };
        runProcessAndCaptureOutputInExceptionMessage(controller, command);
    }

    @Test(groups = {"python"}) // python environment is required to use h5diff
    public void testRunH5DiffIdenticalFile() {
        runH5Diff(String.format("%s/extract.AS.indel.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR),
                  String.format("%s/extract.AS.indel.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR));
    }

    @Test
    public void testRunDiffIdenticalFile() {
        runDiff(String.format("%s/extract.AS.indel.pos.vcf", EXPECTED_TEST_FILES_DIR),
                String.format("%s/extract.AS.indel.pos.vcf", EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class, groups = {"python"}) // python environment is required to use h5diff
    public void testRunH5DiffException() {
        runH5Diff(String.format("%s/extract.AS.indel.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR),
                  String.format("%s/extract.AS.snp.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testRunDiffException() {
        runDiff(String.format("%s/extract.AS.indel.pos.vcf", EXPECTED_TEST_FILES_DIR),
                String.format("%s/extract.AS.snp.pos.vcf", EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testRunDiffNoSuchFile() {
        runDiff(String.format("%s/blahblah", EXPECTED_TEST_FILES_DIR),
                String.format("%s/extract.AS.snp.pos.vcf", EXPECTED_TEST_FILES_DIR));
    }
}
