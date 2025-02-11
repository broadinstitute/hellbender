package org.broadinstitute.hellbender;

import htsjdk.beta.plugin.IOUtils;
import htsjdk.io.IOPath;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
//import org.broadinstitute.hellbender.testutils.SamtoolsTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class PrintFileDiagnosticsIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "fileDiagnosticsTestCases")
    public Object[][] getFileDiagnosticsTestCases() {
        //the pathnames used by the diagnostics tool wind up embedded in the diagnostics output file, so use a
        // relative pathname instead of the named constants (i.e., NA12878_20_21_WGS_cram) in order to avoid test
        // failures in CI caused by the full pathname varying in different environments
        return new Object[][]{
                {
                        "src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.0.samtools.cram",
                        List.of(Pair.of("count-limit", "10")),
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.txt"
                },
                {
                        "src/test/resources/org/broadinstitute/hellbender/metrics/analysis/MeanQualityByCycle/example_pfFail_reads.bam.bai",
                        null,
                        "src/test/resources/filediagnostics/example_pfFail_reads.bam.bai.txt"
                },
                {
                        "src/test/resources/org/broadinstitute/hellbender/engine/cram_with_crai_index.cram.crai",
                        null,
                        "src/test/resources/filediagnostics/cram_with_crai_index.cram.crai.txt"
                },
                {
                    // cram file that uses all the new 3.1 codecs (fqzcomp, name tok, ransNx15, and adaptive arithmetic)
                        "src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.1.samtools.archive.cram",
                        null,
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.1.samtools.archive.cram.txt"
                },
        };
    }

    @Test(dataProvider = "fileDiagnosticsTestCases")
    public void testFileDiagnostics(
            final String inputPath,
            final List<Pair<String, String>> extraArgs,
            final String expectedOutputPath) throws IOException {
        final IOPath outFile = IOUtils.createTempPath("testFileDiagnostics", ".txt");
        runFileDiagnosticsTool(new GATKPath(inputPath), extraArgs, outFile);
        IntegrationTestSpec.assertEqualTextFiles(outFile.toPath().toFile(), new File(expectedOutputPath));
    }

    private void runFileDiagnosticsTool(
            final IOPath inputPath,
            final List<Pair<String, String>> extraArgs,
            final IOPath outputPath) {
        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addInput(inputPath.getRawInputString());
        argBuilder.addOutput(outputPath.getRawInputString());
        if (extraArgs != null) {
            extraArgs.forEach(argPair -> argBuilder.add(argPair.getKey(), argPair.getValue()));
        }
        runCommandLine(argBuilder.getArgsList());
    }

    // TODO: we can't rely on samtools until we reconcile the version used in CI with the version used on the gatk docker,
    // which is currently installed as part of the gatk base image
//    @DataProvider(name = "cram31DiagnosticsTestCases")
//    public Object[][] getCRAM31DiagnosticsTestCases() {
//        return new Object[][]{
//                {
//                        new GATKPath(NA12878_20_21_WGS_cram_31),
//                        new GATKPath(b37_reference_20_21),
//                        List.of(Pair.of("count-limit", "20")),
//                        "--output-fmt cram,version=3.1,archive",
//                        new GATKPath("expectedOutputPath")
//                },
//        };
//    }
//
//    @Test(dataProvider = "cram31DiagnosticsTestCases")
//    public void testCRAM31Diagnostics(
//            final IOPath testInput,
//            final IOPath testReference,
//            final List<Pair<String,String>> diagnosticsArgs,
//            final String samtoolsCommandLineArgs,
//            final IOPath expectedOutputPath) throws IOException {
//        // use samtools to convert the input to CRAM 3.1
//        final IOPath cram31Path = IOUtils.createTempPath("cram31Test", ".cram");
//        SamtoolsTestUtils.convertToCRAM(
//                testInput,
//                cram31Path,
//                testReference,
//                samtoolsCommandLineArgs);
//
//        final IOPath cram31DiagnosticsPath = IOUtils.createTempPath("cram31DiagnosticsTest", ".txt");
//        runFileDiagnosticsTool(cram31Path, diagnosticsArgs, cram31DiagnosticsPath);
//        IntegrationTestSpec.assertEqualTextFiles(cram31DiagnosticsPath.toPath().toFile(), expectedOutputPath.toPath().toFile());
//    }

}
