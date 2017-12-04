package org.broadinstitute.hellbender.tools.spark;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public final class ApplyBQSRSparkIntegrationTest extends CommandLineProgramTest {
    private final static String THIS_TEST_FOLDER = "org/broadinstitute/hellbender/tools/BQSR/";

    private static class ABQSRTest {
        final String bam;
        final String args;
        final String reference;
        final String outputExtension;
        final String expectedFile;

        private ABQSRTest(String bam, String reference, String outputExtension, String args, String expectedFile) {
            this.bam= bam;
            this.reference = reference;
            this.outputExtension = outputExtension;
            this.args = args;
            this.expectedFile = expectedFile;
        }

        @Override
        public String toString() {
            return String.format("ApplyBQSR(args='%s')", args);
        }
    }

    final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
    final String hiSeqCram = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.cram";

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", "-qq -1", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", "-qq 6", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", "-OQ", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", "-SQQ 10 -SQQ 20 -SQQ 30", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", "-SQQ 10 -SQQ 20 -SQQ 30 -RDQ", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030RDQ.bam")});

        //CRAM - input and output crams generated by direct conversion of the corresponding BAM test files with samtools 1.3
        tests.add(new Object[]{new ABQSRTest(hiSeqCram, hg18Reference, ".cram", "--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME + " true", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.cram")});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "ApplyBQSRTestGCS")
    public Object[][] createABQSRTestDataGCS() {
        final String resourceDirGCS = getGCPTestInputPath() + THIS_TEST_FOLDER;
        final String hiSeqBamGCS = resourceDirGCS + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";

        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ABQSRTest(hiSeqBamGCS, null, ".bam", "", resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});

        // TODO: add test inputs with some unaligned reads

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest", groups = "spark")
    public void testApplyBQSR(ABQSRTest params) throws IOException {
        File outFile = GATKBaseTest.createTempFile("applyBQSRTest", params.outputExtension);
        final ArrayList<String> args = new ArrayList<>();
        File refFile = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O"); args.add(outFile.getAbsolutePath());
        if (params.reference != null) {
            refFile = new File(params.reference);
            args.add("-R"); args.add(refFile.getAbsolutePath());
        }
        if (params.args != null) {
            Stream.of(params.args.split(" ")).forEach(arg -> args.add(arg));
        }

        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outFile, new File(params.expectedFile), refFile);
    }

    //TODO: This is disabled because we can't read a google bucket as a hadoop file system outside of the dataproc environment yet
    //Renable when we've figured out how to setup the google hadoop fs connector
    @Test(dataProvider = "ApplyBQSRTestGCS", groups = {"spark", "bucket"}, enabled = false)
    public void testPR_GCS(ABQSRTest params) throws IOException {
        String args =
                " -I " + params.bam +
                        " --apiKey " + getGCPTestApiKey() +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.20mb.1RG.table.gz " +
                        params.args +
                        " -O %s";
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFile));
        spec.executeTest("testPrintReads-" + params.args, this);
    }

    @Test(dataProvider = "ApplyBQSRTestGCS", groups = {"spark", "bucket"}, enabled = false)
    public void testPR_Cloud(ABQSRTest params) throws IOException {
        String args =
                " -I " + params.bam +
                        " --apiKey " + getGCPTestApiKey() +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + getGCPTestInputPath() + THIS_TEST_FOLDER + "HiSeq.20mb.1RG.table.gz " +
                        params.args +
                        " -O %s";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args,
                Arrays.asList(params.expectedFile));
        spec.executeTest("testPrintReads-" + params.args, this);
    }

    @Test(groups = "spark")
    public void testPRNoFailWithHighMaxCycle() throws IOException {
        String args = " -I " + hiSeqBam +
                " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                "" +
                " -O " + createTempFile("ignore",".me");
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString() ,
                Arrays.<String>asList());
        spec.executeTest("testPRNoFailWithHighMaxCycle", this);      //this just checks that the tool does not blow up
    }

    @Test(groups = "spark")
    public void testPRFailWithLowMaxCycle() throws IOException {
        String args =  " -I " + hiSeqBam +
                " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.lowMaxCycle.table.gz" +
                " -O /dev/null";
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                0,
                UserException.class);
        spec.executeTest("testPRFailWithLowMaxCycle", this);
    }
}
