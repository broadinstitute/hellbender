package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class PasteDepthEvidenceIntegrationTest extends CommandLineProgramTest {
    public static final String testDir = toolsTestDir + "sv/PasteDepthEvidence/";

    @Override public String getTestedToolName() { return "PasteDepthEvidence"; }

    @Test
    public void testAllOK() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00096.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00129.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00140.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(testDir + "merged.rd.txt"));
        testSpec.setOutputFileExtension("rd.txt");
        testSpec.executeTest("test1: a-ok", this);
    }

    @Test
    public void testSmallBin() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00096.binsize.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00129.binsize.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00140.binsize.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(testDir + "merged.binsize.rd.txt"));
        testSpec.setOutputFileExtension("rd.txt");
        testSpec.executeTest("test2: small bin", this);
    }

    @Test
    public void testDifferentIntervals() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00096.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00129.binsize.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, testDir + "HG00140.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(), 1, UserException.class);
        testSpec.setOutputFileExtension("rd.txt");
        testSpec.executeTest("test3: different intervals", this);
    }
}
