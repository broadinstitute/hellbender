package org.broadinstitute.hellbender.tools.walkers.bqsr;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.SamReaderFactory;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public final class ApplyBQSRIntegrationTest extends CommandLineProgramTest {
    private static class ABQSRTest {
        final String bam;
        final String reference;
        final String outputExtension;
        final String args[];
        final String expectedFile;

        private ABQSRTest(String bam, String reference, String outputExtension, String args[], String expectedFile) {
            this.bam= bam;
            this.reference = reference;
            this.outputExtension = outputExtension;
            this.args = args;
            this.expectedFile = expectedFile;
        }

        @Override
        public String toString() {
            return String.format("ApplyBQSR(args='%s')", args == null ? "" : StringUtils.join(args));
        }
    }

    @Override
    public String getTestedClassName() {
        return ApplyBQSR.class.getSimpleName();
    }

    final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
    final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.bam";
    final String hiSeqCram = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.cram";
    final String hiSeqBamAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    final String hiSeqCramAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.cram";

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--quantize-quals", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--quantize-quals", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30", "--round-down-quantized"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030RDQ.bam")});

        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--quantize-quals", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--quantize-quals", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30", "--round-down-quantized"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030RDQ.bam")});

        //CRAM - input and output crams generated by direct conversion of the corresponding BAM test files with samtools 1.3
        tests.add(new Object[]{new ABQSRTest(hiSeqCram, hg18Reference, ".cram", new String[] {"--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.cram")});
        tests.add(new Object[]{new ABQSRTest(hiSeqCramAligned, hg18Reference, ".cram", new String[] {"--quantize-quals", "6", "--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.cram")});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "MiniApplyBQSRTest")
    public Object[][] createMiniABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest")
    public void testApplyBQSRFile(ABQSRTest params) throws IOException {
        File outFile = GATKBaseTest.createTempFile("applyBQSRTest", params.outputExtension);
        final ArrayList<String> args = new ArrayList<>();
        File refFile = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O");
        args.add(outFile.getAbsolutePath());
        if (params.reference != null) {
            refFile = new File(params.reference);
            args.add("-R");
            args.add(refFile.getAbsolutePath());
            if (params.args != null) {
                Stream.of(params.args).forEach(arg -> args.add(arg));
            }

            runCommandLine(args);

            SamAssertionUtils.assertSamsEqual(outFile, new File(params.expectedFile), refFile);
        }
    }

    @Test(dataProvider = "MiniApplyBQSRTest")
    public void testApplyBQSRPath(ABQSRTest params) throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path outPath = jimfs.getPath("applyBQSRTest"+params.outputExtension);

            final ArrayList<String> args = new ArrayList<>();
            Path refPath = null;

            args.add("-I");
            args.add(new File(params.bam).getAbsolutePath());
            args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
            args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
            args.add("-O"); args.add(outPath.toUri().toString());
            if (params.reference != null) {
                File refFile = new File(params.reference);
                args.add("-R"); args.add(refFile.getAbsolutePath());
                refPath = refFile.toPath();
            }
            if (params.args != null) {
                Stream.of(params.args).forEach(arg -> args.add(arg));
            }

            runCommandLine(args);

            SamAssertionUtils.assertSamsEqual(outPath, new File(params.expectedFile).toPath(), refPath);
        }
    }

    @Test(dataProvider = "ApplyBQSRTest", groups={"bucket"})
    public void testApplyBQSRCloud(ABQSRTest params) throws IOException {
        // getTempFilePath also deletes the file on exit.
        final String outString = BucketUtils.getTempFilePath(getGCPTestStaging() + "tmp/testApplyBQSRCloud",  params.outputExtension);
        final Path outPath = BucketUtils.getPathOnGcs(outString);
        final ArrayList<String> args = new ArrayList<>();
        Path refPath = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O");
        args.add(outString);
        if (params.reference != null) {
            File refFile = new File(params.reference);
            args.add("-R");
            args.add(refFile.getAbsolutePath());
            refPath = refFile.toPath();
        }
        if (params.args != null) {
            Stream.of(params.args).forEach(arg -> args.add(arg));
        }

        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outPath, new File(params.expectedFile).toPath(), refPath);
    }

    @Test
    public void testMissingReadGroup() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.20mb.1RG.table.missingRG.gz" +
                        " -O /dev/null", 0,
                IllegalStateException.class);
        spec.executeTest("testMissingReadGroup", this);
    }

    @Test
    public void testemptyBqsrRecalFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + createTempFile("emptyBqsrRecal", "").toString() +
                        " -O /dev/null", 0,
                UserException.class);
        spec.executeTest("testemptyBqsrRecalFile", this);
    }

    @Test
    public void testPRNoFailWithHighMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O /dev/null",
                Arrays.<String>asList());
        spec.executeTest("testPRNoFailWithHighMaxCycle", this);      //this just checks that the tool does not blow up
    }


    @Test
    public void testHelp() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --help --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O /dev/null",
                Arrays.<String>asList());
        spec.executeTest("testHelp", this);      //this just checks that the tool does not blow up
    }

    @Test
    public void testPRFailWithLowMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + hiSeqBamAligned +
                        " --"  + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME+ " " + resourceDir + "HiSeq.1mb.1RG.lowMaxCycle.table.gz" +
                        " -O /dev/null",
                0,
                UserException.class);
        spec.executeTest("testPRFailWithLowMaxCycle", this);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndSQQ() throws IOException {
        // --quantize-quals and --static-quantized-quals shouldn't be able to be run in the same command
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBam +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.20mb.1RG.table.gz" +
                        " --static-quantized-quals 9 --quantize-quals 4 " +
                        " -O /dev/null",
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndRDQ() throws IOException {
        // --quantize-quals and --static-quantized-quals shouldn't be able to be run in the same command
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBam +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " "  + resourceDir + "HiSeq.20mb.1RG.table.gz" +
                        " --round-down-quantized --quantize-quals 4 " +
                        " -O /dev/null",
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
    }

    @Test
    public void testOverfiltering() throws IOException {
        final File zeroRefBasesReadBam = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = GATKBaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input", zeroRefBasesReadBam.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--use-original-qualities",
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
        //The expected output is actually the same as inputs for this read
        SamAssertionUtils.assertSamsEqual(outFile, zeroRefBasesReadBam);
    }

    @Test
    public void testAddingPG() throws IOException {
        final File inFile = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = GATKBaseTest.createTempFile("testAddingPG", ".bam");
        final String[] args = new String[] {
                "--input", inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--use-original-qualities",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //The expected output is actually the same as inputs for this read (this ignores the PGs header)
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has no GATK ApplyBQSR in headers
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));

        //output has a GATK ApplyBQSR in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));
    }
}
