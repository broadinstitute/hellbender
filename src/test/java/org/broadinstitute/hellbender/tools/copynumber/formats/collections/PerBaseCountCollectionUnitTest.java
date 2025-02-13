package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.PerBaseCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Unit tests for {@link PerBaseCountCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Robert Klein &lt;rklein@broadinstitute.org&gt;
 */
public final class PerBaseCountCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir + "copynumber/formats/collections");
    private static final File PER_BASE_COUNTS_FILE = new File(TEST_SUB_DIR, "per-base-count-collection-normal.tsv");
    private static final File PER_BASE_COUNTS_NO_N = new File(TEST_SUB_DIR, "per-base-count-collection-normal-no-n.tsv");

    private static final String NORMAL_SAMPLE_NAME_EXPECTED = "HCC1143 BL";

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

    private final PerBaseCountCollection PER_BASE_COUNTS_EXPECTED = new PerBaseCountCollection(
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

    @Test
    public void testRead() {
        final PerBaseCountCollection perBaseCounts = new PerBaseCountCollection(PER_BASE_COUNTS_FILE);
        Assert.assertEquals(perBaseCounts, PER_BASE_COUNTS_EXPECTED);
    }

    // files must have all columns in order to be read in as a PerBaseCountCollection
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadMissingNucleotides() {
        new PerBaseCountCollection(PER_BASE_COUNTS_NO_N);
    }

    @Test
    public void testWrite() throws IOException {
        final File outputFile = createTempFile("per-base-count-collection-test-output", ".tsv");
        PER_BASE_COUNTS_EXPECTED.write(outputFile);
        Assert.assertTrue(FileUtils.contentEquals(outputFile, PER_BASE_COUNTS_FILE));
    }
}
