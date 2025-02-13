package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;


public class CoveredIntervalFinderIntegrationTest extends GatkToolIntegrationTest {

    @Test
    // Baseline test asserting that histogram window boundaries work and that coverage counts are correct when all base quality scores are included.
    public void testBaseOutput10Coverage() throws IOException {
        final File output = createTempFile("testBaseOutput10Coverage", ".bed");
        final File expected = getTestFile("coverageInterval10.bed");

        final String[] args = {
                "-I", largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam",
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "--minimum-depth", "10"
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }
}