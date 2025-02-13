package org.broadinstitute.hellbender.tools;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.beta.io.bundle.Bundle;
import htsjdk.beta.io.bundle.BundleJSON;
import htsjdk.beta.plugin.registry.HaploidReferenceResolver;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadNameReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;

public abstract class AbstractPrintReadsIntegrationTest extends CommandLineProgramTest {

    protected static final File TEST_DATA_DIR = getTestDataDir();

    public void doFileToFile(String fileIn, String extOut, String reference, boolean testMD5) throws Exception {
        String samFile = fileIn;
        final File outFile = GATKBaseTest.createTempFile(samFile + ".", extOut);
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final GATKPath refFile;

        final ArrayList<String> args = new ArrayList<>();
        args.add("--input"); args.add(ORIG_BAM.getAbsolutePath());
        args.add("--output"); args.add(outFile.getAbsolutePath());
        if (reference != null) {
            if (reference.endsWith(BundleJSON.BUNDLE_EXTENSION)) {
                // the test json files are temporary files, not files in TEST_DATA_DIR
                refFile = new GATKPath(reference);
                args.add("-R"); args.add(reference);
            } else {
                refFile = new GATKPath(new File(TEST_DATA_DIR, reference).getAbsolutePath());
                args.add("-R"); args.add(refFile.toString());
            }
        }
        else {
            refFile = null;
        }
        if (testMD5) {
            args.add("--" + StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME);
            args.add("true");
        }
        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM, refFile == null ? null : refFile.toPath().toFile());

        if (testMD5) {
            checkMD5asExpected(outFile);
        }
    }

    public void doFileToFileUsingReferenceBundle(String fileIn, String extOut, String reference, boolean testMD5) throws Exception {
        final String referenceToUse;
        if (reference != null) {
            // create the bundle, using inference to find the sibling files, then write the bundle out to a temp file
            final Bundle referenceBundle = HaploidReferenceResolver.referenceBundleFromFastaPath(
                    new GATKPath(new File(TEST_DATA_DIR, reference).toPath().toString()),
                    GATKPath::new);
            final GATKPath tempBundlePath = new GATKPath(
                    IOUtils.createTempFile("printReadsRefBundle", ".json").getAbsolutePath()
            );
            IOPathUtils.writeStringToPath(tempBundlePath, BundleJSON.toJSON(referenceBundle));
            referenceToUse = tempBundlePath.toString();
        } else {
            referenceToUse = reference;
        }

        // no run the regular test, but using the reference bundle
        doFileToFile(fileIn, extOut, referenceToUse, testMD5);
    }

    private void checkMD5asExpected(final File outFile) throws IOException {
        final File md5File = new File(outFile.getAbsolutePath() + ".md5");
        if (md5File.exists()) {
            md5File.deleteOnExit();
        }
        Assert.assertTrue(md5File.exists(), md5File + " does not exist");
        final String expectedMD5 = Utils.calculateFileMD5(outFile);
        final String actualMD5 = FileUtils.readFileToString(md5File, StandardCharsets.UTF_8);
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test(dataProvider="testingData")
    public void testFileToFileWithReferenceBundle(String fileIn, String extOut, String reference) throws Exception {
        doFileToFileUsingReferenceBundle(fileIn, extOut, reference, false);
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sam", ".sam", null},
                {"print_reads.sam", ".bam", null},
                {"print_reads.sam", ".cram", "print_reads.fasta"},
                {"print_reads.bam", ".sam", null},
                {"print_reads.bam", ".bam", null},
                {"print_reads.bam", ".cram", "print_reads.fasta"},
                {"print_reads.cram", ".sam", "print_reads.fasta"},
                {"print_reads.cram", ".bam", "print_reads.fasta"},
                {"print_reads.cram", ".cram", "print_reads.fasta"},

                {"print_reads.sorted.sam", ".sam", null},
                {"print_reads.sorted.sam", ".bam", null},
                {"print_reads.sorted.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.bam", ".sam", null},
                {"print_reads.sorted.bam", ".bam", null},
                {"print_reads.sorted.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".cram", "print_reads.fasta"},

                {"print_reads.sorted.queryname.sam", ".sam", null},
                {"print_reads.sorted.queryname.sam", ".bam", null},
                {"print_reads.sorted.queryname.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.bam", ".sam", null},
                {"print_reads.sorted.queryname.bam", ".bam", null},
                {"print_reads.sorted.queryname.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".cram", "print_reads.fasta"},

                //test queryname-sorted crams with multiref containers in GATK:
                //print_reads.sorted.queryname_htsjdk_2.1.0.cram was generated from print_reads.sam
                //using gatk4 PrintReads/htsjdk.2.1.0, which includes changes to support
                //multireference containers
                {"print_reads.sorted.queryname.htsjdk-2.1.0.cram", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.htsjdk-2.1.0.cram", ".sam", "print_reads.fasta"}
        };
    }

    @Test(dataProvider="testingData")
    public void testFileToFileUsingReferenceBundle(String fileIn, String extOut, String reference) throws Exception {
        if (reference != null) {
            doFileToFileUsingReferenceBundle(fileIn, extOut, reference, false);
        }
    }

    @Test(dataProvider="testingData")
    public void testFileToFile(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, false);
    }

    @Test
    public void testReadThatConsumesNoReferenceBases() throws IOException {
        final File zeroRefBasesReadBam = new File(TEST_DATA_DIR, "read_consumes_zero_ref_bases.bam");
        final File outFile = GATKBaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input" , zeroRefBasesReadBam.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        // Make sure no exception is thrown given an input containing a read that consumes no reference bases
        runCommandLine(args);

        //Make sure we print the read, ie not lose it.
        SamAssertionUtils.assertSamsEqual(outFile, zeroRefBasesReadBam);
    }

    @DataProvider(name = "UnmappedReadInclusionTestData")
    public Object[][] unmappedReadInclusionTestData() {
        // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
        final File unmappedBam = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_with_unmapped.bam");

        // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
        // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
        // reads with no mapped mate.
        final File ceuSnippet = new File(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");
        final File ceuSnippetCram = new File(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.cram");

        return new Object[][] {
                { unmappedBam, null, Arrays.asList("unmapped"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // The same interval as above in an intervals file
                { unmappedBam, null, Arrays.asList(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped.intervals"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                { unmappedBam, null, Arrays.asList("1:200-300", "unmapped"), Arrays.asList("a", "b", "c", "u1", "u2", "u3", "u4", "u5") },
                { unmappedBam, null, Arrays.asList("1:200-300", "4:700-701", "unmapped"), Arrays.asList("a", "b", "c", "k", "u1", "u2", "u3", "u4", "u5") },
                // The same intervals as above in an intervals file
                { unmappedBam, null, Arrays.asList(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped2.intervals"), Arrays.asList("a", "b", "c", "k", "u1", "u2", "u3", "u4", "u5") },
                { ceuSnippet, null, Arrays.asList("unmapped"), Arrays.asList("g", "h", "h", "i", "i") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000011"), Arrays.asList("a", "b", "c", "d", "e") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000011", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000013", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i") },
                { ceuSnippetCram, b37_reference_20_21, Arrays.asList("unmapped"), Arrays.asList("g", "h", "h", "i", "i") },
                { ceuSnippetCram, b37_reference_20_21, Arrays.asList("20:10000009-10000011", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i") },
                { ceuSnippetCram, b37_reference_20_21, Arrays.asList("20:10000009-10000013", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i") }
        };
    }

    @Test(dataProvider = "UnmappedReadInclusionTestData")
    public void testUnmappedReadInclusion( final File input, final String reference, final List<String> intervalStrings, final List<String> expectedReadNames ) {
        final File outFile = createTempFile("testUnmappedReadInclusion", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("-I"); args.addRaw(input.getAbsolutePath());
        args.addRaw("-O"); args.addRaw(outFile.getAbsolutePath());
        for ( final String intervalString : intervalStrings ) {
            args.addRaw("-L"); args.addRaw(intervalString);
        }
        if ( reference != null ) {
            args.addRaw("-R"); args.addRaw(reference);
        }

        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsPathDataSource(outFile.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                actualReads.add(read);
            }

            Assert.assertEquals(actualReads.size(), expectedReadNames.size(), "Wrong number of reads output");

            for ( int readNumber = 0; readNumber < actualReads.size(); ++readNumber ) {
                Assert.assertEquals(actualReads.get(readNumber).getName(), expectedReadNames.get(readNumber), "Unexpected read name");
            }
        }
    }

    @DataProvider(name="readFilterTestData")
    public static Object[][] testReadFilterData() {
        return new Object[][]{
                {"print_reads_one_malformed_read.sam", null, ".sam", Collections.emptyList(), 7},
                {"print_reads_one_malformed_read.sam", null, ".sam", Arrays.asList("--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS), 8},
                {"print_reads_one_malformed_read.sam", null, ".sam",
                        Arrays.asList("--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, "WellformedReadFilter"), 8},
                {"print_reads.sorted.sam", null, ".sam", Arrays.asList("--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS), 8},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadNameReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "both_reads_align_clip_adapter"),
                        2},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "100",
                                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "200"),
                        8},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "1",
                                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "10"),
                        0},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadNameReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "both_reads_align_clip_adapter",
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "100",
                                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "101"),
                        2},
                {"print_reads.sorted.bam", null, ".sam", Arrays.asList("--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS), 8},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadNameReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "both_reads_align_clip_adapter"),
                        2},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "100",
                                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "101"),
                        8},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadNameReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "both_reads_align_clip_adapter",
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "100",
                                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "101"),
                        2},
                {"print_reads.sorted.cram", "print_reads.fasta", ".sam",
                        Arrays.asList(
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadNameReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "both_reads_align_clip_adapter",
                                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "100",
                                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "101"),
                        2},
        };
    }

    @Test(dataProvider = "readFilterTestData")
    public void testReadFilters(
            final String input,
            final String reference,
            final String extOut,
            final List<String> inputArgs,
            final int expectedCount) throws IOException
    {
        final File outFile = createTempFile("testReadFilter", extOut);

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("-I"); args.addRaw(new File(TEST_DATA_DIR, input).getAbsolutePath());
        args.addRaw("-O"); args.addRaw(outFile.getAbsolutePath());
        if ( reference != null ) {
            args.addRaw("-R"); args.addRaw(new File(TEST_DATA_DIR, reference).getAbsolutePath());
        }
        for (final String filter : inputArgs) {
            args.addRaw(filter);
        }

        runCommandLine(args);


        SamReaderFactory factory = SamReaderFactory.makeDefault();
        if (reference != null) {
            factory = factory.referenceSequence(new File(TEST_DATA_DIR, reference));
        }
        int count = 0;
        try (final SamReader reader = factory.open(outFile)) {
            Iterator<SAMRecord> it = reader.iterator();
            while (it.hasNext()) {
                SAMRecord rec = it.next();
                count++;
            }
        }
        Assert.assertEquals(count, expectedCount);
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testNonExistentReference() throws Exception {
        final File inCram = new File(TEST_DATA_DIR, "print_reads.sorted.cram");
        final File outCram = GATKBaseTest.createTempFile("print_reads_bad_reference", ".cram");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.addRaw(inCram.getCanonicalPath());
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.addRaw(outCram.getCanonicalPath());
        args.addRaw("-R");
        args.addRaw(GATKBaseTest.getSafeNonExistentFile("Nonexistent.fasta").getCanonicalPath());

        runCommandLine(args.getArgsArray());
    }

    @DataProvider
    public Object[][] gcsTestingData() {
        return new Object[][]{
                {"org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam",
                        ".bam", false, null},
                {"org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam",
                        ".bam", true, null},
        };
    }

    /**
     * Test GCS access.
     *
     * For this to work, the settings in src/main/resources/core-site.xml must be correct,
     * and the project name and credential file it points to must be present.
     */
    @Test(dataProvider = "gcsTestingData", groups = "bucket")
    public void testGCSInputsAndOutputs(final String gcsInput, final String outputExtension,
                                        final boolean outputToGCS, final File expectedOutput) {
        final String gcsInputPath = getGCPTestInputPath() + gcsInput;
        final String outputPrefix = outputToGCS ? getGCPTestStaging() : "testGCSInputsAndOutputs";
        final String outputPath = BucketUtils.getTempFilePath(outputPrefix, outputExtension);

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.add("input", gcsInputPath)
                .add("output", outputPath);
        runCommandLine(argBuilder);
    }

    @Test()
    public void testCoordinateSortedInRegion() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.sorted.bam");
        final File expectedBam = new File(getTestDataDir(), "print_reads.sorted.chr1_1.bam");
        final File outBam = GATKBaseTest.createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.addRaw(inBam.getCanonicalPath());
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.addRaw(outBam.getCanonicalPath());
        args.addRaw("-L chr7:1-100 -XL chr7:2-100");

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, expectedBam);
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testSequenceDictionaryValidation() throws Exception {
        final File inCram = new File(getTestDataDir(), "print_reads.sorted.cram");
        final File inRef = new File(getTestDataDir(), "print_reads.chr1only.fasta");
        final File outBam = GATKBaseTest.createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.addRaw(inCram.getCanonicalPath());
        args.addRaw("-R");
        args.addRaw(inRef.getCanonicalPath());
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.addRaw(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());
    }

    @Test()
    public void testUnSorted() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.unsorted.bam");
        try (ReadsDataSource ds = new ReadsPathDataSource(inBam.toPath())){
            Assert.assertEquals(ds.getHeader().getSortOrder(), SAMFileHeader.SortOrder.unsorted);
        }
        final File outBam = GATKBaseTest.createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.addRaw(inBam.getCanonicalPath());
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.addRaw(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

    @Test()
    public void testReadFiltering() throws IOException {
        final File samWithOneMalformedRead = new File(getTestDataDir(), "print_reads_one_malformed_read.sam");
        final File outBam = GATKBaseTest.createTempFile("print_reads_testReadFiltering", ".bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.addRaw(samWithOneMalformedRead.getCanonicalPath());
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.addRaw(outBam.getCanonicalPath());

        runCommandLine(args.getArgsArray());
        SamAssertionUtils.assertSamsEqual(outBam, new File(getTestDataDir(), "expected.print_reads_one_malformed_read.bam"));
    }

    @Test()
    public void testReadFiltering_asIntegrationTestSpec() throws IOException {
        final File samWithOneMalformedRead = new File(getTestDataDir(), "print_reads_one_malformed_read.sam");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --" + StandardArgumentDefinitions.INPUT_LONG_NAME + " " + samWithOneMalformedRead.getCanonicalPath() +
                        " --" + StandardArgumentDefinitions.OUTPUT_LONG_NAME + " " + "%s",
                Arrays.asList(new File(getTestDataDir(), "expected.print_reads_one_malformed_read.bam").getCanonicalPath())
        );
        spec.executeTest("testReadFiltering_asIntegrationTestSpec", this);
    }
}