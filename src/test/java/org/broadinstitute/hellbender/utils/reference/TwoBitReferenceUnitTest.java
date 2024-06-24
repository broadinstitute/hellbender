package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.charset.StandardCharsets;

public class TwoBitReferenceUnitTest extends GATKBaseTest {

    // This FASTA contains some ambiguity codes such as M, which 2bit does not support.
    // Therefore, in the tests that compare the fasta against the corresponding 2bit file,
    // preservation of ambiguity codes is disabled.
    private static final GATKPath largeFasta = new GATKPath(b38_reference_20_21);
    private static final GATKPath largeTwoBit = new GATKPath(b38_2bit_reference_20_21);

    // A small reference with many masked bases, interspersed with N's. The complete sequence is:
    // >chrMaskTest LN:180
    // TTCCAttgTTGTGATTTTGTGctaTTAAAATGATCAAAACANNNCCCTTAAAAATCTTAT
    // TCTAACCTCTCAANNNCTTTTAAAaatgaNNNATTTCAGTACAGTCGGATGCATCTGTAA
    // AAGATAAAAAtaTaACATTGATTAGTTTgCAAAAATAATTGTTTGACCCCAGTTAAGNga
    //
    // Created by running faToTwoBit (https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToTwoBit)
    // on a fasta file with the above contents.
    private static final GATKPath maskedTwoBitTest = new GATKPath(packageRootTestDir + "utils/reference/masked.2bit");

    @Test
    public void testGetSequenceDictionary() {
        try ( final ReferenceDataSource fastaReader = ReferenceDataSource.of(largeFasta.toPath(), false);
              final TwoBitReference twoBitReader = new TwoBitReference(largeTwoBit) ) {

            final SAMSequenceDictionary twoBitDictionary = twoBitReader.getSequenceDictionary();
            final SAMSequenceDictionary fastaDictionary = fastaReader.getSequenceDictionary();

            // Can't use SAMSequenceDictionary.equals(), as it's too strict (checks optional attributes)
            Assert.assertEquals(twoBitDictionary.size(), fastaDictionary.size(), "Size of 2bit dictionary does not match corresponding fasta dictionary");
            for ( int i = 0; i < fastaDictionary.size(); i++ ) {
                final SAMSequenceRecord twoBitRecord = twoBitDictionary.getSequence(i);
                final SAMSequenceRecord fastaRecord = fastaDictionary.getSequence(i);

                Assert.assertEquals(twoBitRecord.getSequenceName(), fastaRecord.getSequenceName(),
                        "2bit sequence dictionary record at index " + i + " has the wrong sequence name");
                Assert.assertEquals(twoBitRecord.getSequenceLength(), fastaRecord.getSequenceLength(),
                        "2bit sequence dictionary record at index " + i + " has the wrong sequence length");
            }
        }
    }

    @Test
    public void testReadEntireReference() {
        try ( final ReferenceDataSource fastaReader = ReferenceDataSource.of(largeFasta.toPath(), false);
              final TwoBitReference twoBitReader = new TwoBitReference(largeTwoBit) ) {

            final SAMSequenceDictionary fastaDictionary = fastaReader.getSequenceDictionary();

            for ( int i = 0; i < fastaDictionary.size(); i++ ) {
                final SAMSequenceRecord fastaRecord = fastaDictionary.getSequence(i);
                final SimpleInterval wholeContigInterval = new SimpleInterval(fastaRecord.getSequenceName(), 1, fastaRecord.getSequenceLength());

                compareTwoBitAndFastaQueryResults(fastaReader, twoBitReader, wholeContigInterval);
            }
        }
    }

    @DataProvider(name = "queryIntervals")
    public Object[][] queryIntervals() {
        return new Object[][] {
                // Just ACGT bases
                { new SimpleInterval("chr20", 60100, 60200)},

                // A stretch of consecutive N's, and then ACGT bases
                { new SimpleInterval("chr20", 59950, 60100) },

                // Multiple nearby blocks of N's separated by ACGT bases
                { new SimpleInterval("chr20", 28646100, 28648200) },

                // A different contig
                { new SimpleInterval("chr21", 20000000, 20001000) },

                // An alt contig
                { new SimpleInterval("chr21_KI270873v1_alt", 100000, 101000) }
        };
    }

    @Test(dataProvider = "queryIntervals")
    public void testQueryValidSubinterval(final SimpleInterval queryInterval) {
        try ( final ReferenceDataSource fastaReader = ReferenceDataSource.of(largeFasta.toPath(), false);
              final TwoBitReference twoBitReader = new TwoBitReference(largeTwoBit) ) {

            compareTwoBitAndFastaQueryResults(fastaReader, twoBitReader, queryInterval);
        }
    }

    private void compareTwoBitAndFastaQueryResults( final ReferenceDataSource fastaReader, final TwoBitReference twoBitReader, final SimpleInterval queryInterval ) {
        final ReferenceSequence fastaSequence = fastaReader.queryAndPrefetch(queryInterval);
        final ReferenceSequence twoBitSequence = twoBitReader.getReferenceBases(queryInterval);

        final byte[] fastaBases = fastaSequence.getBases();
        final byte[] twoBitBases = twoBitSequence.getBases();

        Assert.assertEquals(twoBitBases.length, fastaBases.length,
                "Wrong number of bases returned for " + queryInterval + " by 2bit reader");

        for ( int baseIndex = 0; baseIndex < fastaBases.length; baseIndex++ ) {
            Assert.assertEquals(twoBitBases[baseIndex], fastaBases[baseIndex],
                    "Wrong base returned by 2bit reader at position " + (baseIndex + 1) + " for interval " + queryInterval);
        }
    }

    @DataProvider(name = "invalidIntervals")
    public Object[][] invalidIntervals() {
        return new Object[][] {
                // Contig not in reference
                { new SimpleInterval("chr1", 1, 100) },

                // Interval goes 1 base past the end of a valid contig
                { new SimpleInterval("chr20", 1, 64444168) },

                // Interval goes multiple bases past the end of a valid contig
                { new SimpleInterval("chr20", 1, 64445000) },

                // Interval starts after the end of a valid contig
                { new SimpleInterval("chr20", 64444168, 64445000) }
        };
    }

    @Test(dataProvider = "invalidIntervals", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidIntervals(final SimpleInterval invalidInterval) {
        try ( final TwoBitReference twoBitReader = new TwoBitReference(largeTwoBit) ) {
            twoBitReader.getReferenceBases(invalidInterval);
        }
    }

    @DataProvider(name = "maskedTestIntervals")
    public Object[][] maskedTestIntervals() {
        return new Object[][] {
                // Entire contig
                { new SimpleInterval("chrMaskTest", 1, 180), "TTCCAttgTTGTGATTTTGTGctaTTAAAATGATCAAAACANNNCCCTTAAAAATCTTATTCTAACCTCTCAANNNCTTTTAAAaatgaNNNATTTCAGTACAGTCGGATGCATCTGTAAAAGATAAAAAtaTaACATTGATTAGTTTgCAAAAATAATTGTTTGACCCCAGTTAAGNga" },

                // Just masked bases
                { new SimpleInterval("chrMaskTest", 6, 8), "ttg" },

                // Masked bases, preceded and followed by ACGT bases
                { new SimpleInterval("chrMaskTest", 5, 9), "AttgT" },

                // Start and end within the middle of a block of masked bases,
                // separated by a stretch of ACGT
                { new SimpleInterval("chrMaskTest", 7, 23), "tgTTGTGATTTTGTGct" },

                // Just masked bases and N's
                { new SimpleInterval("chrMaskTest", 85, 92), "aatgaNNN" }
        };
    }

    @Test(dataProvider = "maskedTestIntervals")
    public void testMaskedBases(final SimpleInterval testInterval, final String expectedBases) {
        // Create our TwoBitReference with "preserveCase" set to true for this masking test:
        try ( final TwoBitReference twoBitReader = new TwoBitReference(maskedTwoBitTest, true) ) {
            final byte[] twoBitBases = twoBitReader.getReferenceBases(testInterval).getBases();
            final byte[] expectedBaseArray = expectedBases.getBytes(StandardCharsets.US_ASCII);

            Assert.assertEquals(twoBitBases, expectedBaseArray,
                    "Wrong bases returned from query on 2bit reference. Actual: " + new String(twoBitBases) + " Expected: " + expectedBases);
        }
    }

    @Test
    public void testMaskedReferenceWithUppercasing() {
        final SimpleInterval testInterval = new SimpleInterval("chrMaskTest", 1, 180);
        final String expectedBases = "TTCCATTGTTGTGATTTTGTGCTATTAAAATGATCAAAACANNNCCCTTAAAAATCTTATTCTAACCTCTCAANNNCTTTTAAAAATGANNNATTTCAGTACAGTCGGATGCATCTGTAAAAGATAAAAATATAACATTGATTAGTTTGCAAAAATAATTGTTTGACCCCAGTTAAGNGA";

        // This time create our TwoBitReference with "preserveCase" set to false, to test auto-uppercasing:
        try ( final TwoBitReference twoBitReader = new TwoBitReference(maskedTwoBitTest, false) ) {
            final byte[] twoBitBases = twoBitReader.getReferenceBases(testInterval).getBases();
            final byte[] expectedBaseArray = expectedBases.getBytes(StandardCharsets.US_ASCII);

            Assert.assertEquals(twoBitBases, expectedBaseArray,
                    "Wrong bases returned from query on 2bit reference. Actual: " + new String(twoBitBases) + " Expected: " + expectedBases);
        }
    }
}
