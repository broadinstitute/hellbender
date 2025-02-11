package org.broadinstitute.hellbender.testutils;

import htsjdk.io.HtsPath;
import htsjdk.io.IOPath;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.ProcessExecutor;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SamtoolsTestUtilsTest extends GATKBaseTest {
    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/");

    @Test
    public void testSamtoolsIsAvailable() {
        Assert.assertTrue(SamtoolsTestUtils.isSamtoolsAvailable());
    }

    @Test
    public void testSamtoolsVersion() {
        if (!SamtoolsTestUtils.isSamtoolsAvailable()) {
            throw new SkipException("Samtools not available on local device");
        }
        // If this test runs, but fails because version validation fails, then the local samtools version is
        // not the one expected by the htsjdk tests
        final ProcessExecutor.ExitStatusAndOutput processStatus = SamtoolsTestUtils.executeSamToolsCommand("--version");
        Assert.assertTrue(processStatus.stdout.contains(SamtoolsTestUtils.expectedSamtoolsVersion));
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testSamtoolsPresentButCommandFails() {
        if (!SamtoolsTestUtils.isSamtoolsAvailable()) {
            throw new SkipException("Samtools not available on local device");
        }
        SamtoolsTestUtils.executeSamToolsCommand("--notASamtoolsCommand");
    }

    @Test
    public void testCRAMConversion()throws IOException {
        if (!SamtoolsTestUtils.isSamtoolsAvailable()) {
            throw new SkipException("Samtools not available on local device");
        }

        // Validates CRAM 3.1 conversion.
        final File sourceFile = new File(TEST_DATA_DIR, "print_reads.cram");
        final File cramReference = new File(TEST_DATA_DIR, "print_reads.fasta");
        // This also validates that any extra command line arguments are passed through to samtools by requesting
        // that NM/MD values are synthesized in the output file (which is required for the output records to match).
        final IOPath tempSamtoolsPath = SamtoolsTestUtils.convertToCRAM(
                new HtsPath(sourceFile.getAbsolutePath()),
                new HtsPath(cramReference.getAbsolutePath()),
                "--output-fmt cram,version=3.0,fast");
        final SamReaderFactory factory = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.LENIENT)
                .referenceSequence(cramReference);
        try (final SamReader originalReader = factory.open(sourceFile);
             final SamReader samtoolsCopyReader = factory.open(tempSamtoolsPath.toPath());
             final CloseableIterator<SAMRecord> originalIt = originalReader.iterator();
             final CloseableIterator<SAMRecord> samtoolsIt = samtoolsCopyReader.iterator()) {
            while (originalIt.hasNext() && samtoolsIt.hasNext()) {
                Assert.assertEquals(originalIt.next(), samtoolsIt.next());
            }
            Assert.assertEquals(samtoolsIt.hasNext(), originalIt.hasNext());
        }
    }
}
