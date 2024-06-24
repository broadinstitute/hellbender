package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class SimpleIntervalUnitTest extends GATKBaseTest {

    @DataProvider(name = "badIntervals")
    public Object[][] badIntervals(){
        return new Object[][]{
                {null,1,12, "null contig"},
                {"1", 0, 10, "start==0"},
                {"1", -10, 10, "negative start"},
                {"1", 10, 9, "end < start"}
        };
    }

    @Test(dataProvider = "badIntervals", expectedExceptions = IllegalArgumentException.class)
    public void badIntervals(String contig, int start, int end, String name){
        SimpleInterval interval = new SimpleInterval(contig, start, end);

    }

    @Test(dataProvider = "badIntervals", expectedExceptions = IllegalArgumentException.class)
    public void badIntervalsFromLocatable(String contig, int start, int end, String name){
        Locatable l = getLocatable(contig, start, end);
        new SimpleInterval(l);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void illegalArgumentExceptionFromNullLocatable(){
        new SimpleInterval((Locatable)null);
    }

    @DataProvider(name = "goodIntervals")
    public Object[][] goodIntervals(){
        return new Object[][]{
                {"1", "1", 1, Integer.MAX_VALUE},
                {"1:2", "1", 2, 2},
                {"1:3+", "1", 3, Integer.MAX_VALUE},
                {"1:4-5", "1", 4, 5},
                {"1:2,000", "1", 2000, 2000},
                {"1:3,000+", "1", 3000, Integer.MAX_VALUE},
                {"1:4,000-5,000", "1", 4000, 5000},
                {"1:4,0,0,0-5,0,0,0", "1", 4000, 5000}, //this is OK too, we just remove commas wherever they are
        };
    }

    @Test(dataProvider = "goodIntervals")
    public void testGoodIntervals(String str, String contig, int start, int end){
        SimpleInterval interval = new SimpleInterval(str);
        Assert.assertEquals(interval.getContig(), contig, "contig");
        Assert.assertEquals(interval.getStart(), start, "start");
        Assert.assertEquals(interval.getGA4GHStart(), start-1, "GA4GH start");
        Assert.assertEquals(interval.getEnd(), end, "end");
        Assert.assertEquals(interval.getGA4GHEnd(), end, "GA4GH start");
    }

    @Test(dataProvider = "goodIntervals")
    public void testGoodIntervalsFromLocatable(String ignoreMe, String contig, int start, int end){
        Locatable l = getLocatable(contig, start, end);

        SimpleInterval interval = new SimpleInterval(l);
        Assert.assertEquals(interval.getContig(), contig, "contig");
        Assert.assertEquals(interval.getStart(), start, "start");
        Assert.assertEquals(interval.getEnd(), end, "end");
    }

    private static Locatable getLocatable(final String contig, final int start, final int end) {
        return new Locatable() {
            @Override
            public String getContig() {
                return contig;
            }

            @Override
            public int getStart() {
                return start;
            }

            @Override
            public int getEnd() {
                return end;
            }
        };
    }

    @Test
    public void testEquality(){
        final SimpleInterval i1 = new SimpleInterval("1",1,100);
        final SimpleInterval i2 = new SimpleInterval("1",1,100);
        final SimpleInterval i3 = new SimpleInterval("chr1", 1, 100);
        final SimpleInterval i4 = new SimpleInterval("chr1", 2, 100);
        final SimpleInterval i5 = new SimpleInterval("chr1", 1, 200);
        final SimpleInterval i6 = new SimpleInterval("1:1-100");

        Assert.assertTrue(i1.equals(i1));
        Assert.assertTrue(i1.equals(i2));
        Assert.assertTrue(i2.equals(i1));
        Assert.assertTrue(i1.equals(i6));
        Assert.assertTrue(i6.equals(i1));

        //i1 i2 i6 are same
        Assert.assertEquals(i1.hashCode(), i2.hashCode());
        Assert.assertEquals(i1.hashCode(), i6.hashCode());

        //i3 i4 i5 are are different from each other and from i1 i2 i6
        Assert.assertNotEquals(i1.hashCode(), i3.hashCode());
        Assert.assertNotEquals(i1.hashCode(), i4.hashCode());
        Assert.assertNotEquals(i1.hashCode(), i5.hashCode());
        Assert.assertNotEquals(i3.hashCode(), i4.hashCode());
        Assert.assertNotEquals(i3.hashCode(), i5.hashCode());
        Assert.assertNotEquals(i4.hashCode(), i5.hashCode());

        Assert.assertFalse(i1.equals(i3));
        Assert.assertFalse(i1.equals(i4));
        Assert.assertFalse(i1.equals(i5));
    }

    @Test
    public void testToString(){
        final SimpleInterval i1 = new SimpleInterval("1",1,100);
        Assert.assertEquals(i1.toString(), "1:1-100");
    }

    @DataProvider(name = "IntervalSizeData")
    public Object[][] getIntervalSizeData() {
        // Intervals + expected sizes
        return new Object[][]{
                { new SimpleInterval("1", 1, 1), 1 },
                { new SimpleInterval("1", 1, 2), 2 },
                { new SimpleInterval("1", 1, 10), 10 },
                { new SimpleInterval("1", 2, 10), 9 }
        };
    }

    @Test(dataProvider = "IntervalSizeData")
    public void testGetSize( final SimpleInterval interval, final int expectedSize ) {
        Assert.assertEquals(interval.size(), expectedSize, "size() incorrect for interval " + interval);
    }

    @DataProvider(name = "IntervalOverlapData")
    public static Object[][] getIntervalOverlapData() {
        final SimpleInterval standardInterval = new SimpleInterval("1", 10, 20);
        final SimpleInterval oneBaseInterval = new SimpleInterval("1", 10, 10);

        return new Object[][] {
                { standardInterval, new SimpleInterval("2", 10, 20), false },
                { standardInterval, new SimpleInterval("1", 1, 5), false },
                { standardInterval, new SimpleInterval("1", 1, 9), false },
                { standardInterval, new SimpleInterval("1", 1, 10), true },
                { standardInterval, new SimpleInterval("1", 1, 15), true },
                { standardInterval, new SimpleInterval("1", 10, 10), true },
                { standardInterval, new SimpleInterval("1", 10, 15), true },
                { standardInterval, new SimpleInterval("1", 10, 20), true },
                { standardInterval, new SimpleInterval("1", 15, 20), true },
                { standardInterval, new SimpleInterval("1", 15, 25), true },
                { standardInterval, new SimpleInterval("1", 20, 20), true },
                { standardInterval, new SimpleInterval("1", 20, 25), true },
                { standardInterval, new SimpleInterval("1", 21, 25), false },
                { standardInterval, new SimpleInterval("1", 25, 30), false },
                { oneBaseInterval, new SimpleInterval("2", 10, 10), false },
                { oneBaseInterval, new SimpleInterval("1", 1, 5), false },
                { oneBaseInterval, new SimpleInterval("1", 1, 9), false },
                { oneBaseInterval, new SimpleInterval("1", 1, 10), true },
                { oneBaseInterval, new SimpleInterval("1", 10, 10), true },
                { oneBaseInterval, new SimpleInterval("1", 10, 15), true },
                { oneBaseInterval, new SimpleInterval("1", 11, 15), false },
                { oneBaseInterval, new SimpleInterval("1", 15, 20), false },
                { standardInterval, null, false },
                { standardInterval, standardInterval, true }
        };
    }

    @Test(dataProvider = "IntervalOverlapData")
    public void testOverlap( final SimpleInterval firstInterval, final SimpleInterval secondInterval, final boolean expectedOverlapResult ) {
        Assert.assertEquals(firstInterval.overlaps(secondInterval), expectedOverlapResult,
                            "overlap() returned incorrect result for intervals " + firstInterval + " and " + secondInterval);
    }

    @DataProvider(name = "overlapsWithMargin")
    public Object[][] overlapsWithMargin(){
        final SimpleInterval standardInterval = new SimpleInterval("1", 10, 20);
        final SimpleInterval middleInterval = new SimpleInterval("1", 100, 200);

        return new Object[][] {
                { standardInterval, new SimpleInterval("2", 10, 20), 100, false },
                { standardInterval, new SimpleInterval("1", 1, 15), 0, true },
                { standardInterval, new SimpleInterval("1", 30, 50), 9, false },
                { standardInterval, new SimpleInterval("1", 30, 50), 10, true },
                { middleInterval, new SimpleInterval("1", 50, 99), 0, false },
                { middleInterval, new SimpleInterval("1", 50, 90), 9, false },
                { middleInterval, new SimpleInterval("1", 50, 90), 10, true },
        };
    }
    @Test(dataProvider = "overlapsWithMargin")
    public void testOverlapWithMargin( final SimpleInterval firstInterval, final SimpleInterval secondInterval, int margin, final boolean expectedOverlapResult ) {
        Assert.assertEquals(firstInterval.overlapsWithMargin(secondInterval, margin), expectedOverlapResult,
                "overlap() returned incorrect result for intervals " + firstInterval + " and " + secondInterval);
    }

    @DataProvider()
    private Object[][] overlapsWithMarginExpectingException() {
        final SimpleInterval standardInterval = new SimpleInterval("1", 10, 20);
        final SimpleInterval middleInterval = new SimpleInterval("1", 100, 200);

        return new Object[][] {
                { standardInterval, new SimpleInterval("2", 10, 20), -100 },
                { standardInterval, new SimpleInterval("1", 30, 50), -9 },
                { standardInterval, new SimpleInterval("1", 30, 50), -10 },
                { middleInterval, new SimpleInterval("1", 50, 90), -9 },
                { middleInterval, new SimpleInterval("1", 50, 90), -10 },
        };
    }
    @Test(dataProvider = "overlapsWithMarginExpectingException", expectedExceptions = IllegalArgumentException.class)
    public void testOverlapWithMarginExpectingException( final SimpleInterval firstInterval, final SimpleInterval secondInterval,
                                                         int margin) {
        firstInterval.overlapsWithMargin(secondInterval, margin);
    }

    @DataProvider(name = "IntervalContainsData")
    public Object[][] getIntervalContainsData() {
        final SimpleInterval containingInterval = new SimpleInterval("1", 10, 20);

        return new Object[][] {
                { containingInterval, new SimpleInterval("2", 10, 20), false },
                { containingInterval, new SimpleInterval("1", 1, 5), false },
                { containingInterval, new SimpleInterval("1", 1, 10), false },
                { containingInterval, new SimpleInterval("1", 5, 15), false },
                { containingInterval, new SimpleInterval("1", 9, 10), false },
                { containingInterval, new SimpleInterval("1", 9, 20), false },
                { containingInterval, new SimpleInterval("1", 10, 10), true },
                { containingInterval, new SimpleInterval("1", 10, 15), true },
                { containingInterval, new SimpleInterval("1", 10, 20), true },
                { containingInterval, new SimpleInterval("1", 10, 21), false },
                { containingInterval, new SimpleInterval("1", 15, 25), false },
                { containingInterval, new SimpleInterval("1", 20, 20), true },
                { containingInterval, new SimpleInterval("1", 20, 21), false },
                { containingInterval, new SimpleInterval("1", 20, 25), false },
                { containingInterval, new SimpleInterval("1", 21, 25), false },
                { containingInterval, new SimpleInterval("1", 25, 30), false },
                { containingInterval, null, false },
                { containingInterval, containingInterval, true }
        };
    }


    @Test(dataProvider = "IntervalContainsData")
    public void testContains( final SimpleInterval firstInterval, final SimpleInterval secondInterval, final boolean expectedContainsResult ) {
        Assert.assertEquals(firstInterval.contains(secondInterval), expectedContainsResult,
                            "contains() returned incorrect result for intervals " + firstInterval + " and " + secondInterval);
    }

    @DataProvider(name = "subtractIntervalData")
    private Object[][] subtractIntervalData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 20, 40),
                        new SimpleInterval("chr1", 10, 20) },
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 5, 15),
                        new SimpleInterval("chr1", 15, 30) },
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 10, 20),
                        new SimpleInterval("chr1", 20, 30) },
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 20, 30),
                        new SimpleInterval("chr1", 10, 20) }
        };
    }

    @Test(dataProvider = "subtractIntervalData")
    public void testSubtractInterval( final SimpleInterval firstInterval,
                                      final SimpleInterval secondInterval,
                                      final SimpleInterval expectedResult ) {
        Assert.assertEquals(firstInterval.subtract(secondInterval), expectedResult);
    }

    @DataProvider(name = "subtractIntervalDataExpectingException")
    private Object[][] subtractIntervalDataExpectingException() {
        return new Object[][] {
                // different contigs
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr2", 20, 40) },
                // non-overlapping intervals on same contig
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 50, 150) },
                // second interval contains first
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 10, 40) }
        };
    }
    @Test(dataProvider = "subtractIntervalDataExpectingException", expectedExceptions = IllegalArgumentException.class)
    public void testSubtractIntervalExpectingException( final SimpleInterval firstInterval,
                                                        final SimpleInterval secondInterval) {
        firstInterval.subtract(secondInterval);
    }

    @DataProvider(name = "mergeWithContiguousData")
    private Object[][] mergeWithContiguousData() {
        return new Object[][] {
                // first is upstream, overlapping
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 20, 40),
                        new SimpleInterval("chr1", 10, 40) },
                // first is downstream, overlapping
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 5, 15),
                        new SimpleInterval("chr1", 5, 30) },
                // first contains second
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 15, 20),
                        new SimpleInterval("chr1", 10, 30) },
                // second contains first
                { new SimpleInterval("chr1", 20, 30),
                        new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 10, 30) },
                // first is upstream, overlapping by 1
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 30, 50),
                        new SimpleInterval("chr1", 10, 50) },
                // first is upstream, adjacent
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 31, 50),
                        new SimpleInterval("chr1", 10, 50) },
                // first is downstream, overlapping by 1
                { new SimpleInterval("chr1", 40, 60),
                        new SimpleInterval("chr1", 30, 40),
                        new SimpleInterval("chr1", 30, 60) },
                // first is downstream, adjacent
                { new SimpleInterval("chr1", 40, 60),
                        new SimpleInterval("chr1", 30, 39),
                        new SimpleInterval("chr1", 30, 60) }
        };
    }

    @Test(dataProvider = "mergeWithContiguousData")
    public void testMergeWithContiguous( final SimpleInterval firstInterval,
                                      final SimpleInterval secondInterval,
                                      final SimpleInterval expectedResult ) {
        Assert.assertEquals(firstInterval.mergeWithContiguous(secondInterval), expectedResult);
    }

    @DataProvider(name = "mergeWithContiguousDataExpectingException")
    private Object[][] mergeWithContiguousDataExpectingException() {
        return new Object[][] {
                // different contigs
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr2", 20, 40) },
                // non-contiguous intervals on same contig, first is upstream
                { new SimpleInterval("chr1", 10, 30),
                        new SimpleInterval("chr1", 50, 150) },
                // non-contiguous intervals on same contig, first is downstream
                { new SimpleInterval("chr1", 20, 30),
                        new SimpleInterval("chr1", 5, 15) }
        };
    }
    @Test(dataProvider = "mergeWithContiguousDataExpectingException", expectedExceptions = GATKException.class)
    public void testMergeWithContiguousExpectingException( final SimpleInterval firstInterval,
                                                        final SimpleInterval secondInterval) {
        firstInterval.mergeWithContiguous(secondInterval);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoNullInConstruction() throws Exception {
        new SimpleInterval((String)null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyStringInConstruction() throws Exception {
        new SimpleInterval("");
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadParsePosition() throws Exception {
        new SimpleInterval("chr1:fred");
    }

    @Test(expectedExceptions = GATKException.class)
    public void testNotContiguousLocs() {
        final SimpleInterval loc1 = new SimpleInterval("1", 10, 20);
        final SimpleInterval loc3 = new SimpleInterval("1", 31, 40);
        loc1.mergeWithContiguous(loc3);
    }

    @DataProvider(name = "ExpandWithinContigData")
    public Object[][] expandWithinContigTestData() {
        final int CONTIG_LENGTH = 10000;
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", CONTIG_LENGTH)));

        return new Object[][] {
                { new SimpleInterval("1", 5, 10), 0, CONTIG_LENGTH, dictionary, new SimpleInterval("1", 5, 10) },
                { new SimpleInterval("1", 5, 10), 1, CONTIG_LENGTH, dictionary, new SimpleInterval("1", 4, 11) },
                { new SimpleInterval("1", 1, 10), 10, CONTIG_LENGTH, dictionary, new SimpleInterval("1", 1, 20) },
                { new SimpleInterval("1", 10, 20), 10, CONTIG_LENGTH, dictionary, new SimpleInterval("1", 1, 30) },
                { new SimpleInterval("1", 10, 20), 9, CONTIG_LENGTH, dictionary, new SimpleInterval("1", 1, 29) },
                { new SimpleInterval("1", 30, 40), 5, CONTIG_LENGTH, dictionary, new SimpleInterval("1", 25, 45) },
                { new SimpleInterval("1", CONTIG_LENGTH - 10, CONTIG_LENGTH), 10, CONTIG_LENGTH, dictionary, new SimpleInterval("1", CONTIG_LENGTH - 20, CONTIG_LENGTH) },
                { new SimpleInterval("1", CONTIG_LENGTH - 20, CONTIG_LENGTH - 10), 11, CONTIG_LENGTH, dictionary, new SimpleInterval("1", CONTIG_LENGTH - 31, CONTIG_LENGTH) },
                { new SimpleInterval("1", CONTIG_LENGTH - 20, CONTIG_LENGTH - 10), 10, CONTIG_LENGTH, dictionary, new SimpleInterval("1", CONTIG_LENGTH - 30, CONTIG_LENGTH) }
        };
    }

    @Test(dataProvider = "ExpandWithinContigData")
    public void testExpandWithinContig( final SimpleInterval startingInterval, final int padding, final int contigLength, final SAMSequenceDictionary dictionary, final SimpleInterval expectedInterval ) {
        Assert.assertEquals(startingInterval.expandWithinContig(padding, contigLength), expectedInterval);
        Assert.assertEquals(startingInterval.expandWithinContig(padding, dictionary), expectedInterval);
    }

    @DataProvider(name = "ExpandWithinContigInvalidData")
    public Object[][] expandWithinContigInvalidTestData() {
        final int CONTIG_LENGTH = 10000;
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", CONTIG_LENGTH)));
        final SAMSequenceDictionary badDictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("2", CONTIG_LENGTH)));

        return new Object[][] {
                { new SimpleInterval("1", 1, 10), -1, CONTIG_LENGTH, dictionary },
                { new SimpleInterval("1", 1, 10), 1, 0, dictionary },
                { new SimpleInterval("1", 1, 10), 1, -1, dictionary },
                { new SimpleInterval("1", 1, 10), 1, CONTIG_LENGTH, null },
                { new SimpleInterval("1", 1, 10), 1, CONTIG_LENGTH, badDictionary }
        };
    }

    @Test(dataProvider = "ExpandWithinContigInvalidData", expectedExceptions = IllegalArgumentException.class)
    public void testExpandWithinContigInvalidArgs( final SimpleInterval startingInterval, final int padding, final int contigLength, final SAMSequenceDictionary dictionary ) {
        startingInterval.expandWithinContig(padding, contigLength);
        startingInterval.expandWithinContig(padding, dictionary);
    }
}
