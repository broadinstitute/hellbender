package org.broadinstitute.hellbender.testutils;

import com.google.common.io.Files;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import org.aeonbits.owner.util.Collections;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;

import java.io.*;
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public final class IntegrationTestSpec {

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(IntegrationTestSpec.class);

    public static final String DEFAULT_TEMP_EXTENSION = ".tmp";
    public static final String DEFAULT_TEMP_PREFIX = "walktest.tmp_param";
    public static final Set<String> INDEX_EXTENSIONS = Collections.set(
            FileExtensions.BAI_INDEX,
            FileExtensions.COMPRESSED_VCF_INDEX,
            FileExtensions.CRAM_INDEX,
            FileExtensions.FASTA_INDEX,
            FileExtensions.TABIX_INDEX,
            FileExtensions.TRIBBLE_INDEX,
            FileExtensions.VCF_INDEX
    );

    private final String args;
    private final Class<?> expectedException;
    private final int nOutputFiles;
    private final List<String> expectedFileNames;
    private String tempExtension;
    private boolean trimWhiteSpace;

    //If this field is set to true, bam files will be compared after they get sorted.
    //This is needed as a workaround because Spark tools don't respect a pre-ordered BAMs
    // and so may create BAMs that are sorted differently than the input (though both orders will be valid).
    private boolean compareBamFilesSorted;

    //Stringency for validation of bams.
    private ValidationStringency validationStringency;

    //callback interface for determining if a file is binary or text (by filename). used by EqualZipFiles
    public interface EqualZipFilesAssist {
        boolean isFilenameText(String filename);
    }

    static public EqualZipFilesAssist EqualZipFilesAssist_AllText = new EqualZipFilesAssist() {
        @Override
        public boolean isFilenameText(String filename) {
            return true;
        }
    };

    public IntegrationTestSpec(String args, List<String> expectedFileNames) {
        this.args = args;
        this.nOutputFiles = expectedFileNames.size();
        this.expectedException = null;
        this.expectedFileNames = expectedFileNames;
        this.compareBamFilesSorted = false;
        this.validationStringency = ValidationStringency.DEFAULT_STRINGENCY;
        this.tempExtension = DEFAULT_TEMP_EXTENSION;
        this.trimWhiteSpace = true;
    }

    public IntegrationTestSpec(String args, int nOutputFiles, Class<?> expectedException) {
        if (expectedException == null) {
            throw new IllegalArgumentException("expected exception is null");
        }
        this.args = args;
        this.nOutputFiles = nOutputFiles;
        this.expectedException = expectedException;
        this.expectedFileNames = null;
        this.compareBamFilesSorted = false;
        this.validationStringency = ValidationStringency.DEFAULT_STRINGENCY;
        this.tempExtension = DEFAULT_TEMP_EXTENSION;
        this.trimWhiteSpace = true;
    }

    public boolean expectsException() {
        return expectedException != null;
    }

    public final Class<?> getExpectedException() {
        if (!expectsException())
            throw new GATKException("Tried to get exception for walker test that doesn't expect one");
        return expectedException;
    }

    public void setOutputFileExtension(final String ext) {
        tempExtension = ext;
    }

    public void setCompareBamFilesSorted(final boolean compareBamFilesSorted) {
        this.compareBamFilesSorted = compareBamFilesSorted;
    }

    public void setTrimWhiteSpace(final boolean trimWhiteSpace) {
        this.trimWhiteSpace = trimWhiteSpace;
    }

    public void setValidationStringency(final ValidationStringency validationStringency) {
        this.validationStringency = validationStringency;
    }

    public String getArgs() {
        return args;
    }

    public Collection<String> expectedFileNames() {
        return expectedFileNames;
    }

    public void executeTest(final String name, final CommandLineProgramTester test) throws IOException {
        executeTest(name, test, null);
    }

    public void executeTest(final String name, CommandLineProgramTester test, final String expectedIndexExtension) throws IOException {
        List<File> tmpFiles = new ArrayList<>();
        for (int i = 0; i < nOutputFiles; i++) {
            File fl = BaseTest.createTempFile(String.format(DEFAULT_TEMP_PREFIX + ".%d", i), tempExtension);
            tmpFiles.add(fl);
        }

        final String preFormattedArgs = getArgs();
        final String formattedArgs = String.format(preFormattedArgs, tmpFiles.toArray());
        System.out.println(StringUtils.repeat('-', 80));

        if (expectsException()) {
            // this branch handles the case were we are testing that a walker will fail as expected
            executeTest(name, test, null, null, tmpFiles, formattedArgs, getExpectedException(), expectedIndexExtension);
        } else {
            final List<String> expectedFileNames = new ArrayList<>(expectedFileNames());
            if (!expectedFileNames.isEmpty() && preFormattedArgs.equals(formattedArgs)) {
                throw new GATKException("Incorrect test specification - you're expecting " + expectedFileNames.size() + " file(s) the specified arguments do not contain the same number of \"%s\" placeholders");
            }

            executeTest(name, test, null, expectedFileNames, tmpFiles, formattedArgs, null, expectedIndexExtension);
        }
    }

    /**
     * execute the test, given the following:
     *
     * @param testName               the name of the test
     * @param testClass              the object that contains the test
     * @param expectedFileNames      the list of expectedFileNames
     * @param tmpFiles               the temp file corresponding to the expectedFileNames list
     * @param args                   the argument list
     * @param expectedException      the expected exception or null
     * @param expectedIndexExtension the extension of output indexes or null
     */
    private void executeTest(String testName, CommandLineProgramTester testClass, File outputFileLocation, List<String> expectedFileNames, List<File> tmpFiles, String args, Class<?> expectedException, String expectedIndexExtension) throws IOException {
        if (outputFileLocation != null) {
            args += " -O " + outputFileLocation.getAbsolutePath();
        }
        executeTest(testName, testClass, args, expectedException);

        if (expectedException == null && !expectedFileNames.isEmpty()) {
            assertMatchingFiles(tmpFiles, expectedFileNames, compareBamFilesSorted, validationStringency, trimWhiteSpace);
            if (expectedIndexExtension != null) {
                for (final File f : tmpFiles) {
                    final String indexPath = f.getAbsolutePath() + expectedIndexExtension;
                    Assert.assertTrue(new File(indexPath).exists(), "Index expected at " + indexPath);
                }
            }
        }
    }

    /**
     * execute the test, given the following:
     *
     * @param testName          the name of the test
     * @param testClass         the object that contains the test
     * @param args              the argument list
     * @param expectedException the expected exception or null
     */
    private void executeTest(String testName, CommandLineProgramTester testClass, String args, Class<?> expectedException) {
        String[] command = Utils.escapeExpressions(args);
        // run the executable
        boolean gotAnException = false;
        try {
            final String now = LocalDateTime.now().format(DateTimeFormatter.ofPattern("HH:mm:ss"));
            System.out.println(String.format("[%s] Executing test %s:%s", now, testClass.getClass().getSimpleName(), testName));
            testClass.runCommandLine(command);
        } catch (Exception e) {
            gotAnException = true;
            if (expectedException == null) {
                // we didn't expect an exception but we got one :-(
                throw new RuntimeException(e);
            }
            // we expect an exception
            if (!expectedException.isInstance(e)) {
                final String message = String.format("Test %s:%s expected exception %s but instead got %s with error message %s",
                        testClass, testName, expectedException, e.getClass(), e.getMessage());
                if (e.getCause() != null) {
                    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
                    final PrintStream ps = new PrintStream(baos);
                    e.getCause().printStackTrace(ps);
                    BaseTest.log(message);
                    BaseTest.log(baos.toString());
                }
                Assert.fail(message);
            }
        }

        if (expectedException != null && !gotAnException) {
            // we expected an exception but didn't see it
            Assert.fail(String.format("Test %s:%s expected exception %s but none was thrown", testClass.getClass().getSimpleName(), testName, expectedException.toString()));
        }
    }

    public static void assertMatchingFiles(final List<File> resultFiles, final List<String> expectedFiles, final boolean compareBamFilesSorted, final ValidationStringency stringency) throws IOException {
        assertMatchingFiles(resultFiles, expectedFiles, compareBamFilesSorted, stringency, false);
    }

    public static void assertMatchingFiles(final List<File> resultFiles, final List<String> expectedFiles, final boolean compareBamFilesSorted, final ValidationStringency stringency, final boolean trimWhiteSpace) throws IOException {
        Assert.assertEquals(resultFiles.size(), expectedFiles.size());
        for (int i = 0; i < resultFiles.size(); i++) {
            final File resultFile = resultFiles.get(i);
            final String expectedFileName = expectedFiles.get(i);
            final File expectedFile = new File(expectedFileName);
            final boolean isIndex = INDEX_EXTENSIONS.stream().anyMatch(ext -> expectedFileName.endsWith(ext));
            if (isIndex) {
                Assert.assertTrue(Files.equal(expectedFile, resultFile), "Resulting index file different from expected");
            } else if (expectedFileName.endsWith(".bam")) {
                SamAssertionUtils.assertEqualBamFiles(resultFile, expectedFile, compareBamFilesSorted, stringency);
            } else {
                assertEqualTextFiles(resultFile, expectedFile, null, trimWhiteSpace);
            }
        }
    }

    public static void assertMatchingFiles(final List<Path> resultFiles, final ValidationStringency stringency, final boolean trimWhiteSpace, final List<Path> expectedFiles) throws IOException {
        Assert.assertEquals(resultFiles.size(), expectedFiles.size());
        for (int i = 0; i < resultFiles.size(); i++) {
            final Path resultFile = resultFiles.get(i);
            final String expectedFileName = expectedFiles.get(i).toString();
            final Path expectedFile = expectedFiles.get(i);
            final boolean isIndex = INDEX_EXTENSIONS.stream().anyMatch(ext -> expectedFileName.endsWith(ext));
            if (isIndex) {
                Assert.assertEquals(java.nio.file.Files.readAllBytes(resultFile), java.nio.file.Files.readAllBytes(expectedFile), String.format("Resulting index file %s different from expected",expectedFileName));
            } else if (expectedFileName.endsWith(".bam")) {
                SamAssertionUtils.assertEqualBamFiles(resultFile, expectedFile, null, stringency);
            } else {
                assertEqualTextFiles(resultFile, expectedFile, null, trimWhiteSpace);
            }
        }
    }

    public static void assertEqualTextFiles(final File resultFile, final File expectedFile) throws IOException {
        assertEqualTextFiles(resultFile, expectedFile, null, true);
    }

    public static void assertEqualTextFiles(final File resultFile, final File expectedFile, final String commentPrefix) throws IOException {
        assertEqualTextFiles(resultFile.toPath(), expectedFile.toPath(), commentPrefix, true);
    }

    public static void assertEqualTextFiles(final File resultFile, final File expectedFile, final String commentPrefix, final boolean doTrimWhitespace) throws IOException {
        assertEqualTextFiles(resultFile.toPath(), expectedFile.toPath(), commentPrefix, doTrimWhitespace);
    }

    public static void assertEqualTextFiles(final Path resultFile, final Path expectedFile, final String commentPrefix) throws IOException {
        assertEqualTextFiles(resultFile, expectedFile, commentPrefix, true);
    }

    public static void assertEqualTextFiles(final Path resultFile, final Path expectedFile, final String commentPrefix, final boolean doTrimWhitespace) throws IOException {
        assertEqualTextFiles(IOUtils.makeReaderMaybeGzipped(resultFile), IOUtils.makeReaderMaybeGzipped(expectedFile)
                ,resultFile.toString(), expectedFile.toString(), commentPrefix, doTrimWhitespace);
    }

    /**
     * Compares two text files and ignores all lines that start with the comment prefix.
     */
    public static void assertEqualTextFiles(final Reader resultReader, final Reader expectedReader, final String resultName, final String expectedName, final String commentPrefix, final boolean doTrimWhitespace) throws IOException {

        XReadLines actual = new XReadLines(resultReader, doTrimWhitespace, commentPrefix);
        XReadLines expected = new XReadLines(expectedReader, doTrimWhitespace, commentPrefix);

        // For ease of debugging, we look at the lines first and only then check their counts.
        // For performance, we stream the lines through instead of loading everything first.
        int numUnequalLines = 0;
        int i = 0;
        while (actual.hasNext() && expected.hasNext()) {
          final String expectedLine = expected.next();
          final String actualLine = actual.next();
          if ( !actualLine.equals(expectedLine) ) {
            logger.error( "Line number " + i + " (not counting comments) expected " +
                expectedLine + " actual " + actualLine + '\n' +
                "Expected :" + expectedLine  + '\n' +
                "Actual   :" + actualLine  + '\n'
            );
            ++numUnequalLines;
          }
          i++;
        }

        final boolean sizeMatches = (actual.hasNext() == expected.hasNext());

        // Check our error cases:
        if ( (numUnequalLines != 0) && (!sizeMatches) ) {
            throw new AssertionError("File sizes are unequal - actual = " + actual.readLines().size() + ", expected = " + expected.readLines().size() + " AND detected unequal lines: " + numUnequalLines);
        }
        else if ( numUnequalLines != 0 ) {
            throw new AssertionError("Detected unequal lines: " + numUnequalLines + " between files actual file = "+resultName+", expected file = "+expectedName);
        }
        else if (!sizeMatches) {
            throw new AssertionError("File sizes are unequal - actual = " + (i + actual.readLines().size()) + ", expected = " + (i + expected.readLines().size()) + " between files actual file = "+resultName+", expected file = "+expectedName);
        }
    }

    public static void assertEqualZipFiles(final File resultFile, final File expectedFile, final EqualZipFilesAssist assist) throws IOException {

        // open both as zip files
        final ZipFile             resultZip = new ZipFile(resultFile);
        final ZipFile             expectedZip = new ZipFile(expectedFile);
//
//        File expectedTMP = BaseTest.createTempDir("expected."+expectedFile.getName());
//        File resultTMP = BaseTest.createTempDir("expected."+expectedFile.getName());

        try (FileSystem expectedFS = FileSystems.newFileSystem(expectedFile.toPath(), null);
             FileSystem resultFS = FileSystems.newFileSystem(resultFile.toPath(), null)
        ) {
            // walk the streams, comparing name and content
            final List<Path> resultEntries = StreamSupport.stream(resultFS.getRootDirectories().spliterator(), false).flatMap(root -> {
                try {
                    return java.nio.file.Files.walk(root);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return Stream.empty();
            }).filter(java.nio.file.Files::isRegularFile).collect(Collectors.toList());

            final List<Path> expectedEntries = StreamSupport.stream(expectedFS.getRootDirectories().spliterator(), false).flatMap(root -> {
                try {
                    return java.nio.file.Files.walk(root);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return Stream.empty();
            }).filter(java.nio.file.Files::isRegularFile).collect(Collectors.toList());

            assertMatchingFiles(resultEntries, ValidationStringency.DEFAULT_STRINGENCY, false, expectedEntries);
        }
//            while (resultEntries.hasNext() || expectedEntries.hasNext()) {
//
//                // access entries
//                final Path resultEntry = resultEntries.hasNext() ? resultEntries.next() : null;
//                final Path expectedEntry = expectedEntries.hasNext() ? expectedEntries.next() : null;
//
//                // both should have valid entries
//                if (resultEntry == null) {
//                    throw new AssertionError("result zip missing entries, near " + expectedEntry);
//                } else if (expectedEntry == null) {
//                    throw new AssertionError("result zip has extra extries near " + resultEntry);
//                }
//
//                // name should match
//                if (!resultEntry.getName(-1).equals(expectedEntry.getName(-1))) {
//                    throw new AssertionError("result entry name differs from expected: " + resultEntry.getName() + " != " + expectedEntry.getName());
//                }
//
//                // determine text/binary (names already verified to be the same)
//                final boolean isText = (assist != null) ? assist.isFilenameText(resultEntry.getName()) : false;
//                if (isText) {
//                    IntegrationTestSpec.assertMatchingFiles();
//
//                    final Reader resultReader = new InputStreamReader(new BlockCompressedInputStream(resultZip.getInputStream(resultEntry)));
//                    final Reader expectedReader = new InputStreamReader(new BlockCompressedInputStream(expectedZip.getInputStream(expectedEntry)));
//
//                    assertEqualTextFiles(resultReader, expectedReader, resultEntry.getName(), expectedEntry.getName(), null, true);
//                } else {
//
//                    // binary file must match in size
//                    if (resultEntry.getSize() != expectedEntry.getSize()) {
//                        throw new AssertionError("result entry " + resultEntry.getName() + " differs from expected in size: "
//                                + resultEntry.getSize() + " != " + expectedEntry.getSize());
//                    }
//                }
//            }
//        } finally {
//            resultZip.close();
//            expectedZip.close();
//        }
    }



}
