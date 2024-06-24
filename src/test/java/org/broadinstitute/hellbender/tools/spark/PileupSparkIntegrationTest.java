package org.broadinstitute.hellbender.tools.spark;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration tests for {@link PileupSpark}.
 */
public final class PileupSparkIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "walkers/qc/pileup");

    @DataProvider(name = "shuffle")
    public Object[][] shuffleParameters() {
        return new Object[][] { { false }, { true } };
    }

    private File createAndDeleteTempFile() {
        final File out = IOUtils.createTempFile("out", ".txt");
        out.delete();
        return out;
    }
    
    @Test(dataProvider = "shuffle")
    public void testSimplePileup(boolean useShuffle) throws Exception {
        final File out = createAndDeleteTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(NA12878_20_21_WGS_bam);
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(b37_reference_20_21);
        args.addRaw("-L 20:9999900-10000000");
        if (useShuffle) {
            args.addRaw("--shuffle");
            args.addRaw("--num-reducers").addRaw("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedSimplePileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testVerbosePileup(boolean useShuffle) throws Exception {
        final File out = createAndDeleteTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(NA12878_20_21_WGS_bam);
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(b37_reference_20_21);
        args.addRaw("-L 20:9999990-10000000");
        args.addRaw("-verbose");
        if (useShuffle) {
            args.addRaw("--shuffle");
            args.addRaw("--num-reducers").addRaw("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedVerbosePileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testFeaturesPileup(boolean useShuffle) throws Exception {
        final File out = createAndDeleteTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(NA12878_20_21_WGS_bam);
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(b37_reference_20_21);
        args.addRaw("-L 20:10000092-10000112");
        args.addRaw("-metadata " + dbsnp_138_b37_20_21_vcf);
        if (useShuffle) {
            args.addRaw("--shuffle");
            args.addRaw("--num-reducers").addRaw("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedFeaturesPileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testInsertLengthPileup(boolean useShuffle) throws Exception {
        final File out = createAndDeleteTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(NA12878_20_21_WGS_bam);
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(b37_reference_20_21);
        args.addRaw("-L 20:10000092-10000112");
        args.addRaw("--output-insert-length");
        if (useShuffle) {
            args.addRaw("--shuffle");
            args.addRaw("--num-reducers").addRaw("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedInsertLengthPileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testFeaturesPileupHdfs(boolean useShuffle) throws Exception {
        if (isGATKDockerContainer()) {
            // see https://github.com/eclipse/jetty.project/issues/8549
            // for the docker tests, the test dependencies are in a separate jar
            throw new SkipException("skipping due to jetty jar parsing issues (https://github.com/eclipse/jetty.project/issues/8549)");
        }
        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path vcfPath = new Path(workingDirectory, "dbsnp_138.b37.20.21.vcf");
            final Path idxPath = new Path(workingDirectory, "dbsnp_138.b37.20.21.vcf.idx");
            cluster.getFileSystem().copyFromLocalFile(new Path(dbsnp_138_b37_20_21_vcf), vcfPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(dbsnp_138_b37_20_21_vcf + ".idx"), idxPath);

            final File out = createAndDeleteTempFile();
            final ArgumentsBuilder args = new ArgumentsBuilder();
            args.addRaw("--input");
            args.addRaw(NA12878_20_21_WGS_bam);
            args.addRaw("--output");
            args.addRaw(out.getAbsolutePath());
            args.addRaw("--reference");
            args.addRaw(b37_reference_20_21);
            args.addRaw("-L 20:10000092-10000112");
            args.addRaw("-metadata " + vcfPath.toString());
            if (useShuffle) {
                args.addRaw("--shuffle");
                args.addRaw("--num-reducers").addRaw("1");
            }
            this.runCommandLine(args.getArgsArray());
            File expected = new File(TEST_DATA_DIR, "expectedFeaturesPileup.txt");
            IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);

        });


    }

}