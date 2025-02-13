package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.PerBaseCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.PerBaseCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Integration test for {@link CollectPerBaseCounts}.  Uses a BAM with sites generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Robert Klein &lt;rklein@broadinstitute.org&gt;
 */
public final class CollectPerBaseCountsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "collect-per-base-counts-normal.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "collect-per-base-counts-tumor.bam");
    private static final File SITES_FILE = new File(TEST_SUB_DIR, "collect-per-base-counts-sites.interval_list");

    private static final String NORMAL_SAMPLE_NAME_EXPECTED = "HCC1143 BL";
    private static final String TUMOR_SAMPLE_NAME_EXPECTED = "HCC1143";

    private SAMSequenceRecord getRecordName() {
        final SAMSequenceRecord recordName = new SAMSequenceRecord("20", 1000000);
        recordName.setAttribute("UR", "http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta");
        recordName.setAssembly("GRCh37");
        recordName.setMd5("0dec9660ec1efaaf33281c0d5ea2560f");
        recordName.setSpecies("Homo Sapiens");
        return(recordName);
    }

    private SampleLocatableMetadata getNormalMetadataExpected() {
        final SAMSequenceRecord recordName = this.getRecordName();
        return new SimpleSampleLocatableMetadata(
                NORMAL_SAMPLE_NAME_EXPECTED,
                new SAMSequenceDictionary(Arrays.asList(
                        recordName)));
    }

    private SampleLocatableMetadata getTumorMetadataExpected() {
        final SAMSequenceRecord recordName = this.getRecordName();
        return new SimpleSampleLocatableMetadata(
                TUMOR_SAMPLE_NAME_EXPECTED, new SAMSequenceDictionary(Arrays.asList(
                recordName)));
    }


    @DataProvider(name = "testData")
    public Object[][] testData() {
        //counts from IGV with minMQ = 30 and minBQ = 20
        final PerBaseCountCollection normalCountsExpected = new PerBaseCountCollection(
                this.getNormalMetadataExpected(),
                Arrays.asList(
                        new PerBaseCount(new SimpleInterval("20", 167218, 167218), PerBaseCount.getPerBaseCountFromCounts(0, 0, 10, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167310, 167310), PerBaseCount.getPerBaseCountFromCounts(7, 0, 4, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167325, 167325), PerBaseCount.getPerBaseCountFromCounts(0, 0, 10, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167431, 167431), PerBaseCount.getPerBaseCountFromCounts(5, 0, 3, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167441, 167441), PerBaseCount.getPerBaseCountFromCounts(0, 0, 8, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167455, 167455), PerBaseCount.getPerBaseCountFromCounts(8, 0, 0, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167647, 167647), PerBaseCount.getPerBaseCountFromCounts(0, 0, 0, 12, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167774, 167774), PerBaseCount.getPerBaseCountFromCounts(0, 0, 13, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167833, 167833), PerBaseCount.getPerBaseCountFromCounts(0, 0, 0, 10, 0))));

        final PerBaseCountCollection tumorCountsExpected = new PerBaseCountCollection(
                this.getTumorMetadataExpected(),
                Arrays.asList(
                        new PerBaseCount(new SimpleInterval("20", 167218, 167218), PerBaseCount.getPerBaseCountFromCounts(0, 0, 15, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167310, 167310), PerBaseCount.getPerBaseCountFromCounts(5, 0, 12, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167325, 167325), PerBaseCount.getPerBaseCountFromCounts(0, 0, 20, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167431, 167431), PerBaseCount.getPerBaseCountFromCounts(18, 0, 5, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167441, 167441), PerBaseCount.getPerBaseCountFromCounts(0, 0, 21, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167455, 167455), PerBaseCount.getPerBaseCountFromCounts(19, 0, 0, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167647, 167647), PerBaseCount.getPerBaseCountFromCounts(0, 0, 0, 17, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167774, 167774), PerBaseCount.getPerBaseCountFromCounts(0, 0, 16, 0, 0)),
                        new PerBaseCount(new SimpleInterval("20", 167833, 167833), PerBaseCount.getPerBaseCountFromCounts(0, 0, 0, 10, 0))));

        return new Object[][]{
                {NORMAL_BAM_FILE, normalCountsExpected},
                {TUMOR_BAM_FILE, tumorCountsExpected}
        };
    }

    @Test(dataProvider = "testData")
    public void test(final File inputBAMFile,
                     final PerBaseCountCollection countsExpected) {
        final File outputFile = createTempFile("collect-per-base-counts-test-output", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBAMFile.getAbsolutePath(),
                "-L", SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        final PerBaseCountCollection countsResult = new PerBaseCountCollection(outputFile);
        Assert.assertEquals(countsExpected.getRecords(), countsResult.getRecords());
        Assert.assertEquals(countsExpected.getMetadata(), countsResult.getMetadata());
    }
}
