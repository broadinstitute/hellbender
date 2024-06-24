package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.SparkConf;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ReferenceMultiSparkSourceUnitTest extends GATKBaseTest {
    private String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @Test
    public void testSerializeRoundTrip2Bit() {
        ReferenceMultiSparkSource referenceMultiSource = new ReferenceMultiSparkSource(new GATKPath(twoBitRefURL), ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final ReferenceMultiSparkSource roundTrippedReference = SparkTestUtils.roundTripInKryo(referenceMultiSource, ReferenceMultiSparkSource.class, new SparkConf());

        assertSequenceDictionariesAreEqual(roundTrippedReference.getReferenceSequenceDictionary(null), referenceMultiSource.getReferenceSequenceDictionary(null));
        Assert.assertNotNull(roundTrippedReference.getReferenceWindowFunction());
    }

    // The HTSJDK equality checks rely on interning of Strings, which is difficult to make work with Spark.
    // So we check for equality using String.equals() instead
    private void assertSequenceDictionariesAreEqual(final SAMSequenceDictionary firstDict, final SAMSequenceDictionary secondDict) {
        Assert.assertEquals(firstDict.size(), secondDict.size(), "Mismatch in sizes of sequence dictionaries after serialization");

        for ( int i = 0; i < firstDict.size(); i++ ) {
            final SAMSequenceRecord firstRecord = firstDict.getSequence(i);
            final SAMSequenceRecord secondRecord = secondDict.getSequence(i);

            Assert.assertEquals(firstRecord.getSequenceName(), secondRecord.getSequenceName(),
                    "Mismatching contig names in sequence dictionaries after serialization");
            Assert.assertEquals(firstRecord.getSequenceLength(), secondRecord.getSequenceLength(),
                    "Mismatching contig lengths in sequence dictionaries after serialization");
            Assert.assertEquals(firstRecord.getAttributes(), secondRecord.getAttributes(),
                    "Mismatching attributes in sequence dictionaries after serialization");
        }
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testBadReferenceFile() {
        new ReferenceMultiSparkSource(
                new GATKPath(GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath()),
                ReferenceWindowFunctions.IDENTITY_FUNCTION);
    }

}
