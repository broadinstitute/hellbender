package org.broadinstitute.hellbender.utils.io;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.util.FileExtensions;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.logging.log4j.core.util.FileUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class IOUtilsUnitTest extends GATKBaseTest {

    private FileSystem jimfs;

    @BeforeClass
    public void setUp() {
        jimfs = Jimfs.newFileSystem(Configuration.unix());
    }

    @AfterClass
    public void tearDown() {
        try {
            jimfs.close();
        } catch (final IOException e) {
            log("Unable to close JimFS");
        }
    }

    private List<String> getSortedRecursiveDirectoryContentsFileNames(final Path directoryPath) throws IOException {
        return Files.walk(directoryPath)
                .map(f -> f.toFile().getName())
                .sorted(Comparator.naturalOrder())
                .collect(Collectors.toList());
    }

    private void assertContentsTheSame(final Path baseActualPath, final Path baseExpectedPath) {

        try {
            final List<String> expectedFiles = getSortedRecursiveDirectoryContentsFileNames(baseExpectedPath);
            final List<String> actualFiles = getSortedRecursiveDirectoryContentsFileNames(baseActualPath);

            // Make sure the list of output files is what we expect:
            if (!actualFiles.equals(expectedFiles)) {
                final String actualFilesString = String.join(",", actualFiles);
                final String expectedFilesString = String.join(",", expectedFiles);

                // Duplicating a little work here for a good error message:
                Assert.assertEquals(actualFiles, expectedFiles, "Actual file list differs from expected file list: ACTUAL=[" + actualFilesString + "] , EXPECTED=[" + expectedFilesString + "]");
            }

            // Check that the files and directories are the same:
            Files.find(baseActualPath, Integer.MAX_VALUE,
                    (actualPath, fileAttributes) -> {

                        // First check that the corresponding file exists in our expected unzipped archive:
                        final Path expectedPath = baseExpectedPath.resolve(baseActualPath.relativize(actualPath));

                        final Path actualFileName = actualPath.getFileName();
                        final Path expectedFileName = expectedPath.getFileName();

                        Assert.assertEquals(actualFileName, expectedFileName);

                        // For regular files, make sure the contents of the corresponding files are the same:
                        if ( fileAttributes.isRegularFile() && !fileAttributes.isSymbolicLink()) {
                            try {
                                final byte[] actualContents   = Files.readAllBytes(actualPath);
                                final byte[] expectedContents = Files.readAllBytes(expectedPath);

                                Assert.assertEquals(actualContents, expectedContents);
                            }
                            catch (final IOException ex) {
                                throw new GATKException("Could not verify contents of files: " + actualPath.toUri() + " AND " + expectedPath.toUri(), ex);
                            }
                        }

                        return true;
                    }
            );
        }
        catch (final IOException ex) {
            throw new GATKException("Could not verify identical contents of test directories!");
        }
    }

    @DataProvider
    private Object[][] provideForTestExtractTarGz() {
        return new Object[][] {
                { null },
                { createTempDir("provideForTestExtractTarGz_Case2").toPath() },
        };
    }

    @Test(dataProvider = "provideForTestExtractTarGz")
    public void testExtractTarGz(final Path destDir) {

        final Path outputDataSourcesPath;

        // Prepare our output directory:
        if ( destDir == null ) {
            final File tmpDir = createTempDir("IOUtilsUnitTest_testExtractTarGz");
            final Path tmpDirPath            = tmpDir.toPath();
            outputDataSourcesPath = tmpDirPath.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());

            // Copy our data sources to the destination folder:
            try {
                Files.copy(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ), outputDataSourcesPath);
            }
            catch (final IOException ex) {
                throw new GATKException("Could not copy files for testing!", ex);
            }

            // Extract the files:
            IOUtils.extractTarGz( outputDataSourcesPath );
        }
        else {
            // Extract the files:
            IOUtils.extractTarGz( IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ), destDir );

            outputDataSourcesPath = destDir.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());
        }

        // Get the base expected data sources output path:
        final Path baseExpectedOutputDataSourcesPath = IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_FOLDER);

        // Get the base file path for the output:
        final Path baseActualOutputDataSourcesPath = outputDataSourcesPath.getParent().resolve(
                baseExpectedOutputDataSourcesPath.getFileName()
        );

        // Now compare the files with the files in the known output folder:
        assertContentsTheSame(baseActualOutputDataSourcesPath, baseExpectedOutputDataSourcesPath);
    }

    @Test(expectedExceptions = UserException.class)
    public void testExtractTarGzThrowsExceptionOnExistingOutput() {

        // Prepare our output directory:
        final File tmpDir = createTempDir("IOUtilsUnitTest_testExtractTarGz");
        final Path tmpDirPath            = tmpDir.toPath();
        final Path outputDataSourcesPath = tmpDirPath.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());

        // Copy our data sources to the destination folder:
        try {
            Files.copy(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ), outputDataSourcesPath);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not copy files for testing!", ex);
        }

        // Extract the files:
        IOUtils.extractTarGz( outputDataSourcesPath );

        // Extract them again and explode:
        IOUtils.extractTarGz( outputDataSourcesPath );
    }

    @Test
    public void testExtractTarGzOverwrite() {

        // Prepare our output directory:
        final File tmpDir = createTempDir("IOUtilsUnitTest_testExtractTarGz");
        final Path tmpDirPath            = tmpDir.toPath();
        final Path outputDataSourcesPath = tmpDirPath.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());

        // Copy our data sources to the destination folder:
        try {
            Files.copy(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ), outputDataSourcesPath);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not copy files for testing!", ex);
        }

        // Extract the files:
        IOUtils.extractTarGz( outputDataSourcesPath );

        // Get the base expected data sources output path:
        final Path baseExpectedOutputDataSourcesPath = IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_FOLDER);

        // Get the base file path for the output:
        final Path baseActualOutputDataSourcesPath = outputDataSourcesPath.getParent().resolve(
                baseExpectedOutputDataSourcesPath.getFileName()
        );

        // Now compare the files with the files in the known output folder:
        assertContentsTheSame(baseActualOutputDataSourcesPath, baseExpectedOutputDataSourcesPath);

        // Extract them again and overwrite them:
        IOUtils.extractTarGz( outputDataSourcesPath, outputDataSourcesPath.getParent(), true );

        // Now compare the files with the files in the known output folder AGAIN:
        assertContentsTheSame(baseActualOutputDataSourcesPath, baseExpectedOutputDataSourcesPath);
    }

    @Test
    public void testWriteTarGz() throws Exception {
        final File tarFile = createTempFile("tar", ".tar");

        // Prepare our output directory:
        final File tmpDir = createTempDir("IOUtilsUnitTest_testExtractTarGz");
        final Path tmpDirPath            = tmpDir.toPath();
        final Path outputDataSourcesPath = tmpDirPath.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());

        // Copy our data sources to the destination folder:
        try {
            Files.copy(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ), outputDataSourcesPath);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not copy files for testing!", ex);
        }

        // Extract the files:
        IOUtils.extractTarGz( outputDataSourcesPath );

        IOUtils.writeTarGz(tarFile.getAbsolutePath(), tmpDir.listFiles());
    }

    @Test
    public void testCreateFifoFile() {

        // Create an output location for the test files to go:
        final File tmpDir = createTempDir("IOUtilsUnitTest_testCreateFifoFile");
        final Path tmpDirPath = tmpDir.toPath();

        // Create a FIFO file:
        final Path fifoFilePath = tmpDirPath.resolve(IOUtils.getPath("FIFOFILE"));
        IOUtils.createFifoFile( fifoFilePath);

        // Verify the FIFO file exists and is a FIFO file:
        Assert.assertTrue( Files.exists(fifoFilePath) );
        try {
            Assert.assertTrue( Files.readAttributes(fifoFilePath, BasicFileAttributes.class).isOther() );
        }
        catch (final IOException ex) {
            throw new GATKException("Could not verify file is a FIFO file!", ex);
        }

        // Create the same FIFO file again and overwrite it:
        IOUtils.createFifoFile( fifoFilePath, true );

        // Verify the FIFO file exists and is a FIFO file:
        Assert.assertTrue( Files.exists(fifoFilePath) );
        try {
            Assert.assertTrue( Files.readAttributes(fifoFilePath, BasicFileAttributes.class).isOther() );
        }
        catch (final IOException ex) {
            throw new GATKException("Could not verify file is a FIFO file!", ex);
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testCreateFifoFileThrows() {

        // Create an output location for the test files to go:
        final Path tmpDirPath = createTempDir("IOUtilsUnitTest_testCreateFifoFileThrows").toPath();

        // Create a FIFO file:
        final Path fifoFilePath = tmpDirPath.resolve(IOUtils.getPath("FIFOFILE"));
        IOUtils.createFifoFile( fifoFilePath );

        // Verify the FIFO file exists and is a FIFO file:
        Assert.assertTrue( Files.exists(fifoFilePath) );
        try {
            Assert.assertTrue( Files.readAttributes(fifoFilePath, BasicFileAttributes.class).isOther() );
        }
        catch (final IOException ex) {
            throw new GATKException("Could not verify file is a FIFO file!", ex);
        }

        // Create the same FIFO file again and overwrite it:
        // This will throw:
        IOUtils.createFifoFile( fifoFilePath );
    }

    @Test
    public void testTempDir() {
        File tempDir = IOUtils.createTempDir("Q-Unit-Test");
        Assert.assertTrue(tempDir.exists());
        Assert.assertFalse(tempDir.isFile());
        Assert.assertTrue(tempDir.isDirectory());
        boolean deleted = IOUtils.tryDelete(tempDir);
        Assert.assertTrue(deleted);
        Assert.assertFalse(tempDir.exists());
    }

    @Test
    public void testIsSpecialFile() {
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/null")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/full")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/stdout")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/stderr")));
        Assert.assertFalse(IOUtils.isSpecialFile(null));
        Assert.assertFalse(IOUtils.isSpecialFile(new File("/home/user/my.file")));
        Assert.assertFalse(IOUtils.isSpecialFile(new File("/devfake/null")));
    }

    @DataProvider( name = "ByteArrayIOTestData")
    public Object[][] byteArrayIOTestDataProvider() {
        return new Object[][] {
                // file size, read buffer size
                { 0,     4096 },
                { 1,     4096 },
                { 2000,  4096 },
                { 4095,  4096 },
                { 4096,  4096 },
                { 4097,  4096 },
                { 6000,  4096 },
                { 8191,  4096 },
                { 8192,  4096 },
                { 8193,  4096 },
                { 10000, 4096 }
        };
    }

    @Test( dataProvider = "ByteArrayIOTestData" )
    public void testWriteThenReadFileIntoByteArray ( int fileSize, int readBufferSize ) {
        File tempFile = createTempFile(String.format("testWriteThenReadFileIntoByteArray_%d_%d", fileSize, readBufferSize), "tmp");

        byte[] dataWritten = getDeterministicRandomData(fileSize);
        IOUtils.writeByteArrayToFile(dataWritten, tempFile);
        byte[] dataRead = IOUtils.readFileIntoByteArray(tempFile, readBufferSize);

        Assert.assertEquals(dataRead.length, dataWritten.length);
        Assert.assertEquals(dataRead, dataWritten);
    }

    @Test( dataProvider = "ByteArrayIOTestData" )
    public void testWriteThenReadStreamIntoByteArray ( int fileSize, int readBufferSize ) throws Exception {
        File tempFile = createTempFile(String.format("testWriteThenReadStreamIntoByteArray_%d_%d", fileSize, readBufferSize), "tmp");

        byte[] dataWritten = getDeterministicRandomData(fileSize);
        IOUtils.writeByteArrayToStream(dataWritten, new FileOutputStream(tempFile));
        byte[] dataRead = IOUtils.readStreamIntoByteArray(new FileInputStream(tempFile), readBufferSize);

        Assert.assertEquals(dataRead.length, dataWritten.length);
        Assert.assertEquals(dataRead, dataWritten);
    }

    @Test( expectedExceptions = UserException.CouldNotReadInputFile.class )
    public void testReadNonExistentFileIntoByteArray() {
        File nonExistentFile = GATKBaseTest.getSafeNonExistentFile("djfhsdkjghdfk");
        Assert.assertFalse(nonExistentFile.exists());

        IOUtils.readFileIntoByteArray(nonExistentFile);
    }

    @Test( expectedExceptions = IllegalArgumentException.class )
    public void testReadStreamIntoByteArrayInvalidBufferSize() throws Exception {
        IOUtils.readStreamIntoByteArray(new FileInputStream(createTempFile("testReadStreamIntoByteArrayInvalidBufferSize", "tmp")),
                -1);
    }

    private byte[] getDeterministicRandomData ( int size ) {
        Utils.resetRandomGenerator();
        Random rand = Utils.getRandomGenerator();

        byte[] randomData = new byte[size];
        rand.nextBytes(randomData);

        return randomData;
    }

    @Test
    public void testDeleteDirOnExit() throws IOException {
        //This just tests that the code runs without crashing.
        //It runs at jvm shutdown so there isn't a good way to test it properly.
        //If you see a directory in the hellbender main folder called
        final File dir = new File(GATKBaseTest.publicTestDir + "I_SHOULD_HAVE_BEEN_DELETED");
        IOUtils.deleteOnExit(dir.toPath());

        FileUtils.mkdir(dir, true);
        final File subdir = new File(dir, "subdir");
        FileUtils.mkdir(subdir, true);
        File someFile = new File(dir, "someFile");
        Assert.assertTrue(someFile.createNewFile());
        File anotherFile = new File(subdir, "anotherFile");
        Assert.assertTrue(anotherFile.createNewFile());
    }

    @DataProvider(name = "extensionsToReplace")
    public Object[][] getExtensionsToReplace(){
        return new Object[][] {
                {"file.old", "file.new"},
                {"file.something.old", "file.something.new"},
                {"src/test.something/file", "src/test.something/file.new"},
                {"/.src.folder/some/thing/.secret/file.old", "/.src.folder/some/thing/.secret/file.new" }
        };
    }

    @Test(dataProvider = "extensionsToReplace")
    public void testReplaceExtension(String input, String expected){
        Assert.assertEquals(IOUtils.replaceExtension(input, "new"), expected);
        Assert.assertEquals(IOUtils.replaceExtension(input, "new"), IOUtils.replaceExtension(input,"..new"));
        Assert.assertEquals(IOUtils.replaceExtension(new File(input), "new"), new File(expected));
    }

    @Test(groups={"bucket"})
    public void testGetPath() throws IOException {
        innerTestGetPath(getGCPTestInputPath() + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        innerTestGetPath("file://" + NA12878_20_21_WGS_bam);
        innerTestGetPath(NA12878_20_21_WGS_bam);
    }

    @Test
    public void testGetPathHandlesIntervals() {
        // Make sure we don't crash if passing intervals to getPath.
        // Also, it shouldn't crash when we check whether the file exists.
        Assert.assertFalse(Files.exists(IOUtils.getPath("chr1:10-11")));
        Assert.assertFalse(Files.exists(IOUtils.getPath("chr1:10")));
        Assert.assertFalse(Files.exists(IOUtils.getPath("1:10-11")));
        Assert.assertFalse(Files.exists(IOUtils.getPath("1:10")));
    }

    private void innerTestGetPath(String s) throws IOException {
        Path p = IOUtils.getPath(s);
        long size = Files.size(p);
        Assert.assertTrue(size>0);
    }

    @DataProvider
    public Object[][] absoluteNames() {
        return new Object[][] {
                {"relative/example.txt", new File("relative/example.txt").getAbsolutePath()},
                {"/local/example.txt", "/local/example.txt"},
                // note that path normalization removes the extra / in file://
                {"/local/file://example.txt", "/local/file:/example.txt"},
                {"file:///local/example.txt", "/local/example.txt"},
                {"gs://dir/example.txt", "gs://dir/example.txt"}
        };
    }

    @Test(dataProvider = "absoluteNames")
    public void testGetAbsolutePathWithoutFileProtocol(final String uriString, final String expected) {
        Assert.assertEquals(IOUtils.getAbsolutePathWithoutFileProtocol(IOUtils.getPath(uriString)), expected);
    }

    @Test
    public void testAppendPathToDir() throws Exception {
        Assert.assertEquals(IOUtils.appendPathToDir("dir", "file"), "dir/file");
        Assert.assertEquals(IOUtils.appendPathToDir("dir/", "file"), "dir/file");
        Assert.assertEquals(IOUtils.appendPathToDir("dir", "/file"), "/file");
        Assert.assertEquals(IOUtils.appendPathToDir("dir/", "/file"), "/file");
        Assert.assertEquals(IOUtils.appendPathToDir("/path/to/dir", "anotherdir/file"), "/path/to/dir/anotherdir/file");

        // see https://github.com/eclipse/jetty.project/issues/8549
        if (isGATKDockerContainer()) {
            // for the docker tests, the test dependencies are in a separate jar
            throw new SkipException("skipping due to jetty jar parsing issues (https://github.com/eclipse/jetty.project/issues/8549)");
        }
        // hdfs: URI
        MiniDFSCluster cluster = null;
        try {
            cluster = MiniClusterUtils.getMiniCluster();
            Path tempPath = IOUtils.getPath(MiniClusterUtils.getWorkingDir(cluster).toUri().toString());
            Assert.assertEquals(IOUtils.appendPathToDir(tempPath.toString(), "temp"), tempPath.toString() + "/temp");
        }
        finally {
            MiniClusterUtils.stopCluster(cluster);
        }

        // gs: URI
        Assert.assertEquals(IOUtils.appendPathToDir("gs://abucket/dir", "file"), "gs://abucket/dir/file");
        Assert.assertEquals(IOUtils.appendPathToDir("gs://abucket/dir/", "file"), "gs://abucket/dir/file");

        // file: URI
        Assert.assertEquals(IOUtils.appendPathToDir("file:///dir", "file"), "file:///dir/file");
        Assert.assertEquals(IOUtils.appendPathToDir("file:///dir/", "file"), "file:///dir/file");
    }

    @Test
    public void testSuccessfulCanReadFileCheck() {
        final File expectedFile = createTempFile("Utils-can-read-test",".txt");
        IOUtils.canReadFile(expectedFile);
    }

    @Test
    public void testSuccessfulCanReadFilesCheck() {
        final File file1 = createTempFile("Utils-can-read-test1",".txt");
        final File file2 = createTempFile("Utils-can-read-test2",".txt");
        IOUtils.canReadFile(file1, file2);
    }

    @DataProvider(name = "unsuccessfulCanReadFileCheckData")
    public Object[][] unsuccessfulCanReadFileCheckData() {
        final File directory = createTempDir("Utils-can-read-file-Dir");
        final File nonExistingFile = createTempFile("Utils-cant-read-NoFile", ".file");
        Assert.assertTrue(nonExistingFile.delete());
        final File nonReadable = createTempFile("Utils-cant-read-NotReadable", ".file");

        //Note: if this test suite is run as root (eg in a Docker image), setting readable to false may fail.
        // So we check it here and skip this test if we can't make a file non readable.
        final boolean successSetNotReadable = nonReadable.setReadable(false) && !nonReadable.canRead();
        if (successSetNotReadable) {
            return new Object[][]{{directory}, {nonExistingFile}, {nonReadable}};
        } else {
            //pass in a special value to be able to throw a SkipException and make the same number of tests
            logger.debug("cannot make a file unreadable (maybe you're running as root)");
            return new Object[][]{{directory}, {nonExistingFile}, {null}};
        }
    }

    @Test(dataProvider = "unsuccessfulCanReadFileCheckData",
            expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testUnsuccessfulCanReadFileCheck(final File file) {
        if (file == null){
            throw new SkipException("cannot make a file unreadable (maybe you're running as root)");
        }
        IOUtils.canReadFile(file);
    }

    @Test
    public void testCreateDirectory() throws IOException {

        // see https://github.com/eclipse/jetty.project/issues/8549
        if (isGATKDockerContainer()) {
            // for the docker tests, the test dependencies are in a separate jar
            throw new SkipException("skipping due to jetty jar parsing issues (https://github.com/eclipse/jetty.project/issues/8549)");
        }
        // hdfs
        Path tempPath = IOUtils.getPath(MiniClusterUtils.getWorkingDir(MiniClusterUtils.getMiniCluster()).toUri().toString())
                .resolve("temp");
        IOUtils.createDirectory(tempPath.toUri().toString());
        Assert.assertTrue(java.nio.file.Files.exists(tempPath));

        // local
        tempPath = IOUtils.getPath(BaseTest.createTempDir("xxx").getAbsolutePath().concat("/sub"));
        IOUtils.createDirectory(tempPath.toUri().toString());
        Assert.assertTrue(java.nio.file.Files.exists(tempPath));
    }

    @DataProvider(name="urlEncodeDecode")
    public Object[][] urlEncodeDecode() {
        return new Object[][]{
                // string, url encoding
                { "string", "string"},
                { "string.", "string."},
                { "string1", "string1"},
                { "string with space", "string+with+space"},
                { "string://", "string%3A%2F%2F"},
        };
    }

    @Test(dataProvider = "urlEncodeDecode")
    public void testUrlEncodeDecode(final String string, final String encoded) {
        Assert.assertEquals(IOUtils.urlEncode(string), encoded);
        Assert.assertEquals(string, IOUtils.urlDecode(encoded));
    }

    @DataProvider(name="resourcePaths")
    public Object[][] getResourcePaths() {
        final String testResourceContents = "this is a test resource file";
        return new Object[][] {
                { Resource.LARGE_RUNTIME_RESOURCES_PATH + "/testResourceFile.txt", null , testResourceContents},
                { "testResourceFile.txt", IOUtils.class , testResourceContents},
        };
    }

    @Test(dataProvider = "resourcePaths")
    public void testWriteTempResource(
            final String resourcePath, final Class<?> relativeClass, final String expectedFirstLine) throws IOException
    {
        final Resource largeResource = new Resource(resourcePath, relativeClass);
        final File resourceFile = IOUtils.writeTempResource(largeResource);
        final String resourceContentsFirstLine = getFirstLineAndDeleteTempFile(resourceFile);
        Assert.assertEquals(resourceContentsFirstLine, expectedFirstLine);
    }

    @Test(dataProvider = "resourcePaths")
    public void testWriteTempResourceFromPath(
            final String resourcePath, final Class<?> relativeClass, final String expectedFirstLine) throws IOException
    {
        final File resourceFile = IOUtils.writeTempResourceFromPath(resourcePath, relativeClass);
        final String resourceContentsFirstLine = getFirstLineAndDeleteTempFile(resourceFile);
        Assert.assertEquals(resourceContentsFirstLine, expectedFirstLine);
    }

    @Test
    public void testCreateTempFileInDirectory() {
        final File tempDir = createTempDir("testCreateTempFileInDirectory");
        final File tempFile = IOUtils.createTempFileInDirectory("testCreateTempFileInDirectory", ".txt", tempDir);
        Assert.assertTrue(tempFile.exists(), "file was not written to temp file: " + tempFile);
        Assert.assertEquals(tempFile.getParentFile().getAbsolutePath(), (tempDir.getAbsolutePath()),
                "file was not written to temp file: " + tempFile + " in dir: " + tempDir);
    }

    @Test
    public void testCreateTempFileWithContent() throws IOException {
        final String prefix = "prefix";
        final String suffix = ".txt";
        final List<String> inputStrings = Arrays.asList("first line", "second line");

        final File tempFile = IOUtils.writeTempFile(inputStrings,prefix, suffix);

        Assert.assertTrue(tempFile.exists(), "file was not written to temp file: " + tempFile);
        Assert.assertTrue(tempFile.getName().startsWith("prefix"));
        Assert.assertTrue(tempFile.getName().endsWith(".txt"));

        final List<String> actualStrings = new ArrayList<>();
        try (Stream<String> lines = Files.lines( tempFile.toPath())) {
            lines.forEach(actualStrings::add);
        }

        Assert.assertEquals(actualStrings, inputStrings);
    }

    @DataProvider
    public Object[][] tmpPathDirs() throws Exception {
        return new Object[][] {
                {createTempDir("local").toPath()},
                {Files.createDirectory(jimfs.getPath("tmp"))}
        };
    }

    @Test(dataProvider = "tmpPathDirs", singleThreaded = true)
    public void testTempPath(final Path tempDir) throws Exception {
        // store the previous tmp.dir to check that it is working
        final String previousTmpDir = System.getProperty("java.io.tmpdir");
        try {
            System.setProperty("java.io.tmpdir", IOUtils.getAbsolutePathWithoutFileProtocol(tempDir));
            final Path tempFile = IOUtils.createTempPath(UUID.randomUUID().toString(), ".txt");
            Assert.assertTrue(Files.exists(tempFile),
                    "file was not written to temp file: " + tempFile);
            Assert.assertEquals(tempFile.getParent().toUri().toString(), tempDir.toUri().toString(),
                    "file was not written to temp file: " + tempFile + " in dir: " + tempDir);
        } finally {
            System.setProperty("java.io.tmpdir", previousTmpDir);
        }
    }

    private String getFirstLineAndDeleteTempFile(final File tempResourceFile) throws IOException {
        try (final FileReader fr = new FileReader(tempResourceFile);
             final BufferedReader br = new BufferedReader(fr)) {
            return br.readLine();
        } finally {
            tempResourceFile.delete();
        }
    }

    @DataProvider
    private Object[][] getIsHDF5TestFiles() {
        return new Object[][] {
                { getToolTestDataDir() + "/isValidHDF5.hdf5", true },
                { getToolTestDataDir() + "isValidHDF5.ext", true },
                { getToolTestDataDir() + "/isTSV.tsv", false },
                { getToolTestDataDir() + "/isTSV.hdf", false },
                { getToolTestDataDir() + "/isTSV.hdf5", false },
                { getToolTestDataDir() + "/isEmpty.txt", false },
        };
    }

    @Test(dataProvider = "getIsHDF5TestFiles")
    public void testIsHDF5File(final String filePath, final boolean expected) {
        final Path testPath = Paths.get(filePath);
        Assert.assertEquals(IOUtils.isHDF5File(testPath), expected);
    }

    @DataProvider(name = "GenomicsDBTestPathData")
    public Object[][] genomicsDBTestPathData() {
        return new Object[][]{
                //path, getGenomicsDBPath, getAbsolutePathWithGenomicsDBURIScheme, isGenomicsDBPath
                {null, null, null, false},
                {"", null, null, false},
                {"dfdfdf://fdfdf", null, null, false},
                {"fdfdf", null, null, false},
                {"gendbdfdfdf://fdfdf", null, null, false},
                {"gendb-dfdfdf://fdfdf", null, null, false},
                {"gendb", null, null, false},
                {"gendbdfdf", null, null, false},
                {"agendb://dfdfd", null, null, false},
                {"gendb.://fdfdf", null, null, false},
                {"gendb.", null, null, false},

                {"gendb.dfdfdf://fdfdf", "dfdfdf://fdfdf", "gendb.dfdfdf://fdfdf", true},
                {"gendb://fdfdf", "fdfdf", "gendb://" + new File("fdfdf").getAbsolutePath(), true},
                {"gendb:///fdfd", "/fdfd", "gendb:///fdfd", true},
                {"gendb:///", "/", "gendb:///", true},
                {"gendb.hdfs://this-node:9000/dir", "hdfs://this-node:9000/dir", "gendb.hdfs://this-node:9000/dir", true},
                {"gendb.gs://my-bucket/dir", "gs://my-bucket/dir", "gendb.gs://my-bucket/dir", true},
                {"gendb.s3://my-bucket/dir", "s3://my-bucket/dir", "gendb.s3://my-bucket/dir", true},

                {"gendb-hdfs://this-node:9000/dir", null, null, false},
                {"gendb-gs://this-node:9000/dir", null, null, false}
        };
    }

    @Test(dataProvider = "GenomicsDBTestPathData")
    public void testGenomicsDBPathParsing(String path, String expectedPath, String gendbExpectedAbsolutePath, boolean expectedComparison) {
        Assert.assertEquals(IOUtils.getGenomicsDBPath(path), expectedPath, "getGenomicsDBPath() returned the wrong value");
        if (path != null) {
            // GATKPath can't take a null path
            Assert.assertEquals(IOUtils.getAbsolutePathWithGenomicsDBURIScheme(new GATKPath(path)), gendbExpectedAbsolutePath, "getAbsolutePathWithGenDBScheme() returned the wrong value");
        }
        Assert.assertEquals(IOUtils.isGenomicsDBPath(path), expectedComparison, "isGenomicsDBPath() returned the wrong value");
    }

    @Test
    public void testDeleteRecursively() throws IOException {
        final Path dir = Files.createTempDirectory("test-dir").normalize();
        final Path file = dir.resolve("new-file");
        Files.createFile(file);
        Assert.assertTrue(Files.exists(file));
        IOUtils.deleteRecursively(dir);
        Assert.assertFalse(Files.exists(dir));
        Assert.assertFalse(Files.exists(file));
    }

    @Test(groups={"bucket"})
    public void testDeleteRecursivelyGCS() throws IOException {
        final String gcsFolder = BucketUtils.randomRemotePath(getGCPTestStaging(), "test-dir", "");
        final Path gcsFolderPath = IOUtils.getPath(gcsFolder+"/");
        Files.createDirectory(gcsFolderPath);
        Assert.assertTrue(Files.exists(gcsFolderPath));
        Assert.assertTrue(Files.isDirectory(gcsFolderPath));
        final Path gcsFilePath = Files.createFile(IOUtils.getPath(gcsFolder+"/"+"new-file"));
        Assert.assertTrue(Files.exists(gcsFilePath));
        IOUtils.deleteRecursively(gcsFolderPath);
        Assert.assertFalse(Files.exists(gcsFilePath));
        Assert.assertFalse(Files.exists(IOUtils.getPath(gcsFolder)));
    }

    @DataProvider(name = "ExtensionRemoverProvider")
    public Object[][] getExtensionData() {
        // Files, expected to support serial iteration (false if any input is a .sam)
        return new Object[][] {
                { "" },
                { "/path/to" },
                { "gs://bucket/subdir" }
        };
    }

    @Test(dataProvider = "ExtensionRemoverProvider")
    public void testRemoveExtension(final String basePathString) {
        final Path baseFilePath = Paths.get(basePathString, "test_file");
        final List<String> extensionsList = new ArrayList<>(FileExtensions.VCF_LIST.size() + 1);
        extensionsList.add("");
        extensionsList.addAll(FileExtensions.VCF_LIST);
        for (final String extension : extensionsList) {
            final Path inputPath = Paths.get(baseFilePath.toString() + extension);
            final Path outputPath = IOUtils.removeExtension(inputPath, FileExtensions.VCF_LIST);
            Assert.assertEquals(outputPath, baseFilePath);
        }
    }
}
