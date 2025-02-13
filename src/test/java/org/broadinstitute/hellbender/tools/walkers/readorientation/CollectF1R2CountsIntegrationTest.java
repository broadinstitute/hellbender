package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2TestingUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

/**
 * Created by tsato on 3/14/18.
 */
public class CollectF1R2CountsIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testOnSyntheticBam() throws IOException {
        final File outputTarGz = createTempFile("f1r2", ".tar.gz");
        final String sample = "SAMPLE";

        final int numAltReads = 30;
        final int numRefReads = 70;
        final int depth = numAltReads + numRefReads;
        final File samFile = createSyntheticSam(numRefReads, numAltReads, sample);

        final String[] args = {
                "-R", hg19_chr1_1M_Reference,
                "-I", samFile.getAbsolutePath(),
                "-O", outputTarGz.getAbsolutePath()
        };

        runCommandLine(args);

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(outputTarGz.toPath(), extractedDir.toPath());

        final File refHistogramFile = F1R2CountsCollector.getRefHistogramsFromExtractedTar(extractedDir).stream().findFirst().get();

        final MetricsFile<?, Integer> referenceSiteMetrics = new MetricsFile<>();
        final Reader in = IOUtil.openFileForBufferedReading(refHistogramFile);
        referenceSiteMetrics.read(in);
        CloserUtil.close(in);

        List<Histogram<Integer>> histograms = referenceSiteMetrics.getAllHistograms();

        /** Expected result
         *
         * (position chr1: {@code alignmentStart})
         *
         * coord.:    99,990             100,000             100,010             100,020
         *              |                   |                   |                   |
         * reference:   T C A T C A C A C T C A C T A A G C A C A C A G A G A A T A A T
         * alt read:    T C A T C A C A C T C T C T G A C C A A A G A G A A T A A T A A
         *                                    *     *   *     *
         *
         * At each alt site, we have 75 ref and 25 alt reads
         *
         * context, alt counts, alt f1r2 counts, depth
         * -------------------------------------------
         * CAC, {70, 0, 0, 30}, {35, 0, 0, 15}, 100
         * TAA, {70, 0, 30, 0}, {35, 0, 15, 0}, 100
         * AGC, {0, 30, 70, 0}, {0, 15, 35, 0}, 100
         * ACA, {30, 70, 0, 0}, {15, 35, 0, 0}, 100
         *
         *
         * Ref: 6 sites at depth 100
         * TCA, CAT, ATC, TCA, CAC, ACA, CAC, ACT, CTC, TCA, ACT, CTA, AAG, GCA, CAC, CAC, ACA, CAG, AGA, GAG, AGA, GAA, AAT, ATA, TAA, AAT
         * TCA*3, CAT, ATC, CAC*4, ACA*2, ACT*2, CTC, CTA, AAG, GCA, CAG, AGA*2, GAG, GAA, AAT*2, ATA, TAA
         *
         **/

        // alt site tests
        final File altTableFileCAC = F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir, "CAC").stream().findFirst().get();
        List<AltSiteRecord> altDesignMatrixCAC = AltSiteRecord.readAltSiteRecords(altTableFileCAC.toPath()).getRight();
        AltSiteRecord recordCAC = altDesignMatrixCAC.get(0);
        Assert.assertEquals(recordCAC.getRefCount(), 70);
        Assert.assertEquals(recordCAC.getAltCount(), 30);
        Assert.assertEquals(recordCAC.getRefF1R2(), 35);
        Assert.assertEquals(recordCAC.getAltF1R2(), 15);

        final File altTableFileTAA = F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir, "TAA").stream().findFirst().get();
        List<AltSiteRecord> altDesignMatrixTAA = AltSiteRecord.readAltSiteRecords(altTableFileCAC.toPath()).getRight();
        AltSiteRecord recordTAA = altDesignMatrixCAC.get(0);
        Assert.assertEquals(recordTAA.getRefCount(), 70);
        Assert.assertEquals(recordTAA.getAltCount(), 30);
        Assert.assertEquals(recordTAA.getRefF1R2(), 35);
        Assert.assertEquals(recordTAA.getAltF1R2(), 15);

        final File altTableFileAGC = F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir, "AGC").stream().findFirst().get();
        List<AltSiteRecord> altDesignMatrixAGC = AltSiteRecord.readAltSiteRecords(altTableFileCAC.toPath()).getRight();
        AltSiteRecord recordAGC = altDesignMatrixCAC.get(0);
        Assert.assertEquals(recordAGC.getRefCount(), 70);
        Assert.assertEquals(recordAGC.getAltCount(), 30);
        Assert.assertEquals(recordAGC.getRefF1R2(), 35);
        Assert.assertEquals(recordAGC.getAltF1R2(), 15);

        final File altTableFileACA = F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir, "ACA").stream().findFirst().get();
        List<AltSiteRecord> altDesignMatrixACA = AltSiteRecord.readAltSiteRecords(altTableFileCAC.toPath()).getRight();
        AltSiteRecord recordACA = altDesignMatrixCAC.get(0);
        Assert.assertEquals(recordACA.getRefCount(), 70);
        Assert.assertEquals(recordACA.getAltCount(), 30);
        Assert.assertEquals(recordACA.getRefF1R2(), 35);
        Assert.assertEquals(recordACA.getAltF1R2(), 15);

        // check ref histograms
        // ImmutableMap.of can take up to 5 elements, and ImmutableMap.Builder is ugly
        final ImmutableMap<String, Integer> expectedReferenceCounts1 = ImmutableMap.of(
                "TCA", 3,
                "CAT", 1,
                "ATC", 1,
                "CAC", 4,
                "ACA", 2);
        final ImmutableMap<String, Integer>  expectedReferenceCounts2 = ImmutableMap.of(
                "ACT", 2,
                "CTC", 1,
                "CTA", 1,
                "AAG", 1,
                "GCA", 1);
        final ImmutableMap<String, Integer>  expectedReferenceCounts3 = ImmutableMap.of(
                "CAG", 1,
                "AGA", 2,
                "GAG", 1,
                "GAA", 1,
                "AAT", 2);
        final ImmutableMap<String, Integer>  expectedReferenceCounts4 = ImmutableMap.of("ATA", 1, "TAA", 1);

        for (ImmutableMap<String, Integer> expectedContextCountPairs : Arrays.asList(
                expectedReferenceCounts1, expectedReferenceCounts2, expectedReferenceCounts3, expectedReferenceCounts4)){
            expectedContextCountPairs.entrySet().stream().forEach(entrySet -> {
                final String context = entrySet.getKey();
                final double expectedCount = entrySet.getValue();
                Histogram<Integer> histogram = histograms.stream()
                        .filter(hist -> hist.getValueLabel().equals(context))
                        .findFirst().get();
                Assert.assertEquals(histogram.get(depth).getValue(), expectedCount);
                Assert.assertEquals(histogram.getSumOfValues(), expectedCount );
            });
        }
    }

    @Test
    public void testHistograms() throws IOException {
        final File outputTarGz = createTempFile("f1r2", ".tar.gz");
        final String sample = "SAMPLE";
        final File sam = createSyntheticSam(30, 1, sample);

        final String[] args = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam.getAbsolutePath(),
                "-O", outputTarGz.getAbsolutePath()
        };

        runCommandLine(args);

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(outputTarGz.toPath(), extractedDir.toPath());

        final File altHistogramFile = F1R2CountsCollector.getAltHistogramsFromExtractedTar(extractedDir).stream().findFirst().get();

        final MetricsFile<?, Integer> referenceSiteMetrics = new MetricsFile<>();
        final Reader in = IOUtil.openFileForBufferedReading(altHistogramFile);
        referenceSiteMetrics.read(in);
        CloserUtil.close(in);

        final MetricsFile<?, Integer> altSiteMetrics = new MetricsFile<>();
        final Reader altMetricsReader = IOUtil.openFileForBufferedReading(altHistogramFile);
        altSiteMetrics.read(altMetricsReader);
        CloserUtil.close(altMetricsReader);

        final List<Histogram<Integer>> histograms = altSiteMetrics.getAllHistograms();

        // TODO: should there be 64*3*2 = 384 histograms or just the non-zero ones?
        final String[] expectedTransitions = new String[]{"CAC_T_F2R1", "TAA_G_F2R1", "AGC_C_F2R1", "ACA_A_F2R1"};
        // Assert.assertEquals(histograms.size(), expectedTransitions.length, "alt histogram must contain the expected number of contexts");

        for (String transition : expectedTransitions ){
            Optional<Histogram<Integer>> histogram = histograms.stream()
                    .filter(h -> h.getValueLabel().equals(transition))
                    .findFirst();
            Assert.assertTrue(histogram.isPresent(), "histogram must exist");
            Assert.assertEquals((int) histogram.get().getSumOfValues(), 1, "histogram must only contain one read");
        }
    }

    @Test
    public void testTwoSamples() throws IOException {
        final File outputTarGZ1 = createTempFile("f1r21", ".tar.gz");
        final File outputTarGZ2 = createTempFile("f1r22",".tar.gz");
        final File outputTarGZ12 = createTempFile("f1r212",".tar.gz");
        final String sample1 = "SAMPLE1";
        final String sample2 = "SAMPLE2";
        final File sam1 = createSyntheticSam(10, 1, sample1);
        final File sam2 = createSyntheticSam(20, 2, sample2);

        final String[] args1 = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam1.getAbsolutePath(),
                "-O", outputTarGZ1.getAbsolutePath()
        };

        final String[] args2 = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam2.getAbsolutePath(),
                "-O", outputTarGZ2.getAbsolutePath()
        };

        final String[] args12 = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam1.getAbsolutePath(),
                "-O", outputTarGZ12.getAbsolutePath()
        };

        runCommandLine(args1);
        runCommandLine(args2);
        runCommandLine(args12);
        final File extractedDir1 = createTempDir("extracted1");
        IOUtils.extractTarGz(outputTarGZ1.toPath(), extractedDir1.toPath());

        final File extractedDir12 = createTempDir("extracted12");
        IOUtils.extractTarGz(outputTarGZ12.toPath(), extractedDir12.toPath());

        for (final String context : List.of("ACC", "TGC", "CAA")) {
            IntegrationTestSpec.assertEqualTextFiles(new File(extractedDir1, sample1 + "." + context + F1R2CountsCollector.ALT_TABLE_EXTENSION),
                    new File(extractedDir12, sample1 + "." + context + F1R2CountsCollector.ALT_TABLE_EXTENSION));
        }

        for (final String extension : List.of(F1R2CountsCollector.REF_HIST_EXTENSION, F1R2CountsCollector.ALT_HIST_EXTENSION)) {
            IntegrationTestSpec.assertEqualTextFiles(new File(extractedDir1, sample1 + extension),
                    new File(extractedDir12, sample1 + extension));
        }
    }

    static File createSyntheticSam(final int refDepth, final int altDepth, final String sampleName) throws IOException {
        final File samFile = IOUtils.createTempFile("synthetic", ".bam");
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader(sampleName);
        final SAMFileGATKReadWriter writer = M2TestingUtils.getBareBonesSamWriter(samFile, samHeader);
        //            Ref Sequence: "CATCACACTCACTAAGCACACAGAGAATAAT".getBytes();
        //          SNPs positions:            *  * *  *
        final byte[] altReadBases = "CATCACACTCTCTGACCAAACAGAGAATAAT".getBytes();
        final List<GATKRead> refReads = M2TestingUtils.createReads(refDepth, M2TestingUtils.DEFAULT_REF_BASES, samHeader, (byte)30, "ref");
        final List<GATKRead> alt1Reads = M2TestingUtils.createReads(altDepth, altReadBases, samHeader, (byte)30, "alt");
        refReads.forEach(writer::addRead);
        alt1Reads.forEach(writer::addRead);
        writer.close(); // closing the writer writes to the file

        return samFile;
    }

}