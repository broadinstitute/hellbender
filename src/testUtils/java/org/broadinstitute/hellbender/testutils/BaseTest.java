package org.broadinstitute.hellbender.testutils;

import com.google.common.base.Strings;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.BeforeSuite;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;


/**
 * This is the base test class for all of our test cases.  All test cases should extend from this
 * class; it sets up the logger, and resolves the location of directories that we rely on.
 */
public abstract class BaseTest {

    static {
        // set properties for the local Spark runner
        System.setProperty("dataflow.spark.test.reuseSparkContext", "true");
        SparkContextFactory.enableTestSparkContext();
        System.setProperty("picard.useLegacyParser", "false");
    }

    /**
     * run an external process and assert that it finishes with exit code 0
     * @param processController a ProcessController to use
     * @param command command to run, the 0th element must be the executable name
     */
    public static void runProcess(final ProcessController processController, final String[] command) {
        runProcess(processController, command, "Process exited with non-zero value. Command: "+ Arrays.toString(command) + "\n");
    }

    /**
     * run an external process and assert that it finishes with exit code 0
     * if it exits with non-zero, include the process stdout/stderr in the exception message
     *
     * @param processController a ProcessController to use
     * @param command command to run, the 0th element must be the executable name
     */
    public static void runProcessAndCaptureOutputInExceptionMessage(final ProcessController processController, final String[] command) {
        runProcess(processController, command, null,"Process exited with non-zero value. Command: "+ Arrays.toString(command) + "\n", true);
    }

    /**
     * run an external process and assert that it finishes with exit code 0
     * @param processController a ProcessController to use
     * @param command command to run, the 0th element must be the executable name
     * @param message error message to display on failure
     */
    public static void runProcess(final ProcessController processController, final String[] command, final String message) {
        runProcess(processController, command, null, message);
    }

    /**
     * run an external process and assert that it finishes with exit code 0
     * @param processController a ProcessController to use
     * @param command command to run, the 0th element must be the executable name
     * @param environment what to use as the process environment variables
     * @param message error message to display on failure
     */
    public static void runProcess(final ProcessController processController, final String[] command, final Map<String, String> environment, final String message) {
        runProcess(processController, command, environment, message, false);
    }

    /**
     * run an external process and assert that it finishes with exit code 0
     * @param processController a ProcessController to use
     * @param command command to run, the 0th element must be the executable name
     * @param environment what to use as the process environment variables
     * @param message error message to display on failure
     * @param includeProcessOutputInExceptionMessage if true, include the process's stdout/stderr (if any) in the exception
     *                                               message thrown upon exit with a non-zero value
     */
    public static void runProcess(final ProcessController processController, final String[] command, final Map<String, String> environment, final String message, final boolean includeProcessOutputInExceptionMessage) {
        final ProcessSettings prs = new ProcessSettings(command);
        prs.getStderrSettings().printStandard(true);
        prs.getStdoutSettings().printStandard(true);

        // Capture stdout/stderr to a buffer if we need to include it in the exception message
        if ( includeProcessOutputInExceptionMessage ) {
            prs.getStderrSettings().setBufferSize(65536);
            prs.getStdoutSettings().setBufferSize(65536);
        }

        prs.setEnvironment(environment);
        final ProcessOutput output = processController.exec(prs);

        final String fullMessage = includeProcessOutputInExceptionMessage ?
                message + "\n***Process stdout:***\n" + output.getStdout() + "\n***Process stderr:***\n" + output.getStderr() :
                message;

        Assert.assertEquals(output.getExitValue(), 0, fullMessage);
    }

    /**
     * Spawn a new jvm with the same classpath as this one and run a gatk CommandLineProgram
     * This is useful for running tests that require changing static state that is not allowed to change during
     * a tool run but which needs to be changed to test some condition.
     *
     * @param toolName CommandLineProgram to run
     * @param arguments arguments to provide to the tool
     */
    public static void runToolInNewJVM(String toolName, ArgumentsBuilder arguments){
        runToolInNewJVM(toolName, arguments, System.getenv());
    }

    /**
     * Spawn a new jvm with the same classpath as this one and run a gatk CommandLineProgram
     * This is useful for running tests that require changing static state that is not allowed to change during
     * a tool run but which needs to be changed to test some condition.
     *
     * @param toolName CommandLineProgram to run
     * @param arguments arguments to provide to the tool
     * @param environment a map of key-value pairs which will be used as to set the System environment
     */
    public static void runToolInNewJVM(String toolName, ArgumentsBuilder arguments, Map<String, String> environment){
        final String javaHome = System.getProperty("java.home");
        final String javaBin = javaHome + File.separator + "bin" + File.separator + "java";
        final String classpath = System.getProperty("java.class.path");;
        final List<String> baseCommand = new ArrayList<>(Arrays.asList(
                javaBin,
                "-cp", classpath,
                Main.class.getName(),
                toolName));
        baseCommand.addAll(arguments.getArgsList());
        runProcess(ProcessController.getThreadLocal(), baseCommand.toArray(new String[0]), environment,
                "Java exited with non-zero value. Command: "+ String.join(" ",baseCommand) + "\n");
    }

    @BeforeSuite
    public void setTestVerbosity(){
        LoggingUtils.setLoggingLevel(Log.LogLevel.WARNING);
    }

    public static final Logger logger = LogManager.getLogger("org.broadinstitute.gatk");

    /**
     *  This is a public requester pays bucket owned by the broad-gatk-test project.
     *  It must be owned by a different project than the service account doing the testing or the test may fail because it can access the
     *  file directly through alternative permissions.
     */
    public static final String REQUESTER_PAYS_BUCKET = "gs://hellbender-requester-pays-test/";

    /**
     * A publicly readable GCS bucket set as requester pays, this should not be owned by the same project that is set
     * as {@link BaseTest#getGCPTestProject()} or the tests for requester pays access may be invalid.
     * @return HELLBENDER_REQUESTER_PAYS_BUCKET env. var if defined, {@value BaseTest#REQUESTER_PAYS_BUCKET}.
     */
    public static final String getGCPRequesterPaysBucket(){
        final String valueFromEnvironment = System.getenv("HELLBENDER_REQUESTER_PAYS_BUCKET");
        return Strings.isNullOrEmpty(valueFromEnvironment) ? REQUESTER_PAYS_BUCKET : valueFromEnvironment;
    }

    /**
     * name of the google cloud project that stores the data and will run the code
     * @return HELLBENDER_TEST_PROJECT env. var if defined, throws otherwise.
     */
    public static String getGCPTestProject() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_PROJECT");
    }

    /**
     * A writable GCS path where java files can be cached and temporary test files can be written,
     * of the form gs://bucket/, or gs://bucket/path/.
     * @return HELLBENDER_TEST_STAGING env. var if defined, throws otherwise.
     */
    public static String getGCPTestStaging() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_STAGING");
    }

    /**
     *  A GCS path where the test inputs are stored.
     *
     *  The value of HELLBENDER_TEST_INPUTS should end in a "/" (for example, "gs://hellbender/test/resources/")
     *  
     *  @return HELLBENDER_TEST_INPUTS env. var if defined, throws otherwise.
     */
    public static String getGCPTestInputPath() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_INPUTS");
    }

    /**
     *  A path where the test inputs for the Funcotator LargeDataValidationTest are stored.
     *
     *  The value of FUNCOTATOR_LARGE_TEST_INPUTS should end in a "/" (for example, "gs://hellbender/funcotator/test/resources/")
     *
     *  @return FUNCOTATOR_LARGE_TEST_INPUTS env. var if defined, throws otherwise.
     */
    public static String getFuncotatorLargeDataValidationTestInputPath() {
        return getNonNullEnvironmentVariable("FUNCOTATOR_LARGE_TEST_INPUTS");
    }

    /**
     * A local path where the service account credentials are stored
     * @return GOOGLE_APPLICATION_CREDENTIALS env. var if defined, throws otherwise.
     */
    public static String getGoogleServiceAccountKeyPath() {
      return getNonNullEnvironmentVariable("GOOGLE_APPLICATION_CREDENTIALS");
    }

    protected static String getNonNullEnvironmentVariable(String envVarName) {
        String value = System.getenv(envVarName);
        if (null == value) {
            throw new UserException("For this test, please define environment variable \""+envVarName+"\"");
        }
        return value;
    }

    /**
     * Returns the location of the resource directory for the command line program.
     */
    public String getToolTestDataDir(){
        return "src/test/resources/" +this.getClass().getPackage().getName().replace(".","/") +"/" + getTestedClassName() + "/";
    }

    /**
     * Returns the name of the class being tested.
     * The default implementation takes the simple name of the test class and removes the trailing "Test".
     * Override if needed.
     */
    public String getTestedClassName(){
        if (getClass().getSimpleName().contains("IntegrationTest"))
            return getClass().getSimpleName().replaceAll("IntegrationTest$", "");
        else if (getClass().getSimpleName().contains("UnitTest"))
            return getClass().getSimpleName().replaceAll("UnitTest$", "");
        else
            return getClass().getSimpleName().replaceAll("Test$", "");
    }

    /**
     * @param fileName the name of a file
     * @return a File resolved using getToolTestDataDir as the parent and fileName
     */
    public File getTestFile(String fileName) {
        return new File(getToolTestDataDir(), fileName);
    }

    /**
     * @param fileName the name of a file
     * @return a File resolved using getToolTestDataDir as the parent and fileName
     */
    public GATKPath getTestFileGATKPath(String fileName) {
        return new GATKPath(String.format("%s/%s", getToolTestDataDir(), fileName));
    }

    /**
     * Simple generic utility class to creating TestNG data providers:
     *
     * 1: inherit this class, as in
     *
     *      private class SummarizeDifferenceTest extends TestDataProvider {
     *         public SummarizeDifferenceTest() {
     *           super(SummarizeDifferenceTest.class);
     *         }
     *         ...
     *      }
     *
     * Provide a reference to your class to the TestDataProvider constructor.
     *
     * 2: Create instances of your subclass.  Return from it the call to getTests, providing
     * the class type of your test
     *
     * <code>
     * {@literal @}DataProvider(name = "summaries")
     * public Object[][] createSummaries() {
     *   new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
     *   new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
     *   return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
     * }
     * </code>
     *
     * This class magically tracks created objects of this
     */
    public static class TestDataProvider {
        private static final Map<Class<?>, List<Object>> tests = new LinkedHashMap<>();
        protected String name;

        /**
         * Create a new TestDataProvider instance bound to the class variable C
         */
        public TestDataProvider(Class<?> c, String name) {
            if ( ! tests.containsKey(c) )
                tests.put(c, new ArrayList<>());
            tests.get(c).add(this);
            this.name = name;
        }

        public TestDataProvider(Class<?> c) {
            this(c, "");
        }

        public void setName(final String name) {
            this.name = name;
        }

        /**
         * Return all of the data providers in the form expected by TestNG of type class C
         */
        public static Object[][] getTests(Class<?> c) {
            List<Object[]> params2 = new ArrayList<>();
            for ( Object x : tests.get(c) ) params2.add(new Object[]{x});
            return params2.toArray(new Object[][]{});
        }

        @Override
        public String toString() {
            return "TestDataProvider("+name+")";
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(final String name, final String extension) {
        return IOUtils.createTempFile(name, extension);
    }

    /**
     * Create a temp file with an arbitrary name and extension
     */
    public static File createTempFile(){
        return createTempFile("default", ".tmp");
    }

    /**
     * Creates a temp path that will be deleted on exit after tests are complete.
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the path.
     * @param extension Extension to concat to the end of the path.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static Path createTempPath(final String name, final String extension) {
        return IOUtils.createTempPath(name, extension);
    }

    /**
     * Return a File object representing a file with the given name and extension that is guaranteed not to exist.
     * @param fileNameWithExtension
     * @return File object representing a file that is guaranteed not to exist
     */
    public static File getSafeNonExistentFile(final String fileNameWithExtension) {
        final File tempDir = createTempDir("nonExistentFileHolder");
        final File nonExistingFile = new File(tempDir, fileNameWithExtension);
        return nonExistingFile;
    }

    /**
     * Return a Path object representing a file with the given name and extension that is guaranteed not to exist.
     * @param fileNameWithExtension
     * @return Path object representing a file that is guaranteed not to exist
     */
    public static Path getSafeNonExistentPath(final String fileNameWithExtension) {
        return getSafeNonExistentGATKPath(fileNameWithExtension).toPath();
    }

    /**
     * Return a GATKPath object representing a file with the given name and extension that is guaranteed not to exist.
     * @param fileNameWithExtension
     * @return Path object representing a file that is guaranteed not to exist
     */
    public static GATKPath getSafeNonExistentGATKPath(final String fileNameWithExtension) {
        final File tempDir = createTempDir("nonExistentFileHolder");
        return new GATKPath(String.format("%s/%s", tempDir.getAbsolutePath(), fileNameWithExtension));
    }

    /**
     * Creates an empty temp directory which will be deleted on exit after tests are complete
     *
     * @param prefix prefix for the directory name
     * @return an empty directory starting with prefix which will be deleted after the program exits
     */
    public static File createTempDir(final String prefix){
        return IOUtils.createTempDir(prefix);
    }

    /**
     * Creates an empty temp directory which will be deleted on exit after tests are complete
     *
     * @return an empty directory will be deleted after the program exits
     */
    public static File createTempDir(){
        return createTempDir("tmp");
    }

    /**
     * Log this message so that it shows up inline during output as well as in html reports
     */
    public static void log(final String message) {
        Reporter.log(message, true);
    }

    private static final double DEFAULT_FLOAT_TOLERANCE = 1e-1;


    /**
     * Checks whether two double array contain the same values or not.
     * @param actual actual produced array.
     * @param expected expected array.
     * @param tolerance maximum difference between double value to be consider equivalent.
     */
    protected static void assertEqualsDoubleArray(final double[] actual, final double[] expected, final double tolerance) {
        if (expected == null)
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual.length,expected.length,"array length");
        }
        for (int i = 0; i < actual.length; i++)
            Assert.assertEquals(actual[i],expected[i],tolerance,"array position " + i);
    }

    public static void assertEqualsDoubleSmart(final Object actual, final Double expected, final double tolerance) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double) (Double) actual, (double) expected, tolerance);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected) {
        assertEqualsDoubleSmart(actual, expected, DEFAULT_FLOAT_TOLERANCE);
    }

    public static <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new LinkedHashSet<>(actual);
        final Set<T> expectedSet = new LinkedHashSet<>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance) {
        assertEqualsDoubleSmart(actual, expected, tolerance, null);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance, final String message) {
        if ( Double.isNaN(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isNaN(actual), "expected is nan, actual is not");
        else if ( Double.isInfinite(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isInfinite(actual), "expected is infinite, actual is not");
        else {
            final double delta = Math.abs(actual - expected);
            final double ratio = Math.abs(actual / expected - 1.0);
            Assert.assertTrue(delta < tolerance || ratio < tolerance, "expected = " + expected + " actual = " + actual
                    + " not within tolerance " + tolerance
                    + (message == null ? "" : "message: " + message));
        }
    }

    public static void assertEqualsIntSmart(final int actual, final int expected, final int tolerance, final String message) {
        final double delta = Math.abs(actual - expected);
        final double ratio = Math.abs(actual / expected - 1.0);
        Assert.assertTrue(delta < tolerance || ratio < tolerance, "expected = " + expected + " actual = " + actual
                + " not within tolerance " + tolerance
                + (message == null ? "" : "message: " + message));
    }


    /**
     * captures {@link System#out} while runnable is executing
     * @param runnable a code block to execute
     * @return everything written to {@link System#out} by runnable
     */
    public static String captureStdout(Runnable runnable){
        return captureSystemStream(runnable, System.out, System::setOut);
    }


    /**
     * captures {@link System#err} while runnable is executing
     * @param runnable a code block to execute
     * @return everything written to {@link System#err} by runnable
     */
    public static String captureStderr(Runnable runnable){
        return captureSystemStream(runnable, System.err, System::setErr);
    }

    private static String captureSystemStream(Runnable runnable,  PrintStream stream, Consumer<? super PrintStream> setter){
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        setter.accept(new PrintStream(out));
        try {
            runnable.run();
        } finally{
            setter.accept(stream);
        }
        return out.toString();
    }

    public static void assertContains(String actual, String expectedSubstring){
        Assert.assertTrue(actual.contains(expectedSubstring),  expectedSubstring +" was not found in " + actual + ".");
    }

    public static <T> void assertCondition(Iterable<T> actual, Iterable<T> expected, BiConsumer<T,T> assertion){
        final Iterator<T> iterActual = actual.iterator();
        final Iterator<T> iterExpected = expected.iterator();
        while(iterActual.hasNext() && iterExpected.hasNext()){
            assertion.accept(iterActual.next(), iterExpected.next());
        }
        if (iterActual.hasNext()){
            Assert.fail("actual is longer than expected with at least one additional element: " + iterActual.next());
        }
        if (iterExpected.hasNext()){
            Assert.fail("actual is shorter than expected, missing at least one element: " + iterExpected.next());
        }
    }

    /**
     * assert that the iterable is sorted according to the comparator
     */
    public static <T> void assertSorted(Iterable<T> iterable, Comparator<T> comparator){
        final Iterator<T> iter = iterable.iterator();
        assertSorted(iter, comparator);
    }

    /**
     * assert that the iterator is sorted according to the comparator
     */
    public static <T> void assertSorted(Iterator<T> iterator, Comparator<T> comparator){
        assertSorted(iterator, comparator, null);
    }


    /**
     * assert that the iterator is sorted according to the comparator
     */
    public static <T> void assertSorted(Iterable<T> iterable, Comparator<T> comparator, String message){
       assertSorted(iterable.iterator(), comparator, message);
    }


    /**
     * assert that the iterator is sorted according to the comparator
     */
    public static <T> void assertSorted(Iterator<T> iterator, Comparator<T> comparator, String message){
        T previous = null;
        while(iterator.hasNext()){
            T current = iterator.next();
            if( previous != null) {
                Assert.assertTrue(comparator.compare(previous, current) <= 0, "Expected " + previous + " to be <= " + current + (message == null ? "" : "\n"+message));
            }
            previous = current;
        }
    }

    /**
     * Get a FileSystem that uses the explicit credentials instead of the default
     * credentials.
     *
     * @param bucket the name of the bucket to access.
     *
     * @return A FileSystem for that bucket on GCS, using explicit credentials.
     *
     * @throws IOException
     */
    protected FileSystem getAuthenticatedGcs(final String bucket) throws IOException {
        final byte[] creds = Files.readAllBytes(Paths.get(getGoogleServiceAccountKeyPath()));
        return BucketUtils.getAuthenticatedGcs(getGCPTestProject(), bucket, creds);
    }

    /**
     * Print information to the console that helps understand which credentials are in use
     * (are they the ones you meant to use?).
     */
    protected void helpDebugAuthError() {
        final String key = "GOOGLE_APPLICATION_CREDENTIALS";
        String credsFile = System.getenv(key);
        if (null == credsFile) {
            System.err.println("$" + key + " is not defined.");
            return;
        }
        System.err.println("$" + key + " = " + credsFile);
        Path credsPath = Paths.get(credsFile);
        boolean exists = Files.exists(credsPath);
        System.err.println("File exists: " + exists);
        if (exists) {
            try {
                System.err.println("Key lines from file:");
                printKeyLines(credsPath, "\"type\"", "\"project_id\"", "\"client_email\"");
            } catch (IOException x2) {
                System.err.println("Unable to read: " + x2.getMessage());
            }
        }
    }

    private void printKeyLines(Path path, String... keywords) throws IOException {
        for (String line : Files.readAllLines(path)) {
            for (String keyword : keywords) {
                if (line.contains(keyword)) {
                    System.err.println(line);
                }
            }
        }
    }
}

