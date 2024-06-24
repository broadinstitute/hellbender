package org.broadinstitute.hellbender.utils;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.SimpleFeature;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;


/**
 * test out the interval utility methods
 */
public final class IntervalUtilsUnitTest extends GATKBaseTest {
    public static final String INTERVAL_TEST_DATA = publicTestDir + "org/broadinstitute/hellbender/utils/interval/";
    public static final String emptyIntervals = publicTestDir + "empty_intervals.list";

    private List<GenomeLoc> hg19ReferenceLocs;
    private List<GenomeLoc> hg19exomeIntervals;


    @BeforeClass
    public void init() {
        hg19ReferenceLocs = Collections.unmodifiableList(GenomeLocSortedSet.createSetFromSequenceDictionary(hg19Header.getSequenceDictionary()).toList()) ;
        hg19exomeIntervals = Collections.unmodifiableList(IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(hg19MiniIntervalFile)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullArgSpanningInterval() throws Exception {
        IntervalUtils.getSpanningInterval(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSpanningInterval_differentConfigs() throws Exception {
        final SimpleInterval chr1_200_300 = new SimpleInterval("1:200-300");
        final SimpleInterval chr2_1_100 = new SimpleInterval("2:1-100");

        IntervalUtils.getSpanningInterval(Arrays.asList(chr1_200_300, chr2_1_100));
    }

    @Test
    public void testCompareLocatables() throws Exception {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("1", 1000));
        dict.addSequence(new SAMSequenceRecord("2", 1000));
        final SimpleInterval chr1_1_100 = new SimpleInterval("1:1-100");
        final SimpleInterval chr1_5_100 = new SimpleInterval("1:5-100");
        final SimpleInterval chr2_1_100 = new SimpleInterval("2:1-100");

        // equal intervals comparison return 0
        Assert.assertEquals(IntervalUtils.compareLocatables(chr1_1_100, chr1_1_100, dict), 0);
        Assert.assertEquals(IntervalUtils.compareLocatables(chr1_5_100, chr1_5_100, dict), 0);
        Assert.assertEquals(IntervalUtils.compareLocatables(chr2_1_100, chr2_1_100, dict), 0);
        // first < second return negative
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_1_100, chr1_5_100, dict) < 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_1_100, chr2_1_100, dict) < 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_5_100, chr2_1_100, dict) < 0);
        // first > second return positive
        Assert.assertTrue(IntervalUtils.compareLocatables(chr2_1_100, chr1_1_100, dict) > 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr2_1_100, chr1_5_100, dict) > 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_5_100, chr1_1_100, dict) > 0);
    }

    @DataProvider
    public Object[][] testLocatableOrderingData() {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("1", 1000));
        dict.addSequence(new SAMSequenceRecord("2", 1000));

        return new Object[][] {
                // first, second, dictionary, expectedIsBefore, expectedIsAfter, expectedContigComparison
                { new SimpleInterval("1", 1, 100), new SimpleInterval("1", 200, 300), dict, true, false, 0 },
                { new SimpleInterval("1", 1, 100), new SimpleInterval("1", 101, 200), dict, true, false, 0 },
                { new SimpleInterval("1", 1, 100), new SimpleInterval("1", 100, 200), dict, false, false, 0 },
                { new SimpleInterval("1", 1, 100), new SimpleInterval("1", 99, 200), dict, false, false, 0 },

                { new SimpleInterval("1", 200, 300), new SimpleInterval("1", 1, 100), dict, false, true, 0 },
                { new SimpleInterval("1", 101, 200), new SimpleInterval("1", 1, 100), dict, false, true, 0 },
                { new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1, 100), dict, false, false, 0 },
                { new SimpleInterval("1", 99, 200), new SimpleInterval("1", 1, 100), dict, false, false, 0 },

                { new SimpleInterval("1", 1, 100), new SimpleInterval("2", 1, 100), dict, true, false, -1 },
                { new SimpleInterval("1", 1, 100), new SimpleInterval("2", 101, 200), dict, true, false, -1 },
                { new SimpleInterval("1", 101, 200), new SimpleInterval("2", 1, 100), dict, true, false, -1 },
                { new SimpleInterval("2", 1, 100), new SimpleInterval("1", 1, 100), dict, false, true, 1 },
                { new SimpleInterval("2", 101, 200), new SimpleInterval("1", 1, 100), dict, false, true, 1 },
                { new SimpleInterval("2", 1, 100), new SimpleInterval("1", 101, 200), dict, false, true, 1 }
        };
    }

    @Test(dataProvider = "testLocatableOrderingData")
    public void testLocatableOrdering(final SimpleInterval first, final SimpleInterval second, final SAMSequenceDictionary dict, final boolean expectedIsBefore, final boolean expectedIsAfter, final int expectedContigComparison) {
        Assert.assertEquals(IntervalUtils.isBefore(first, second, dict), expectedIsBefore, "Wrong return value from IntervalUtils.isBefore()");
        Assert.assertEquals(IntervalUtils.isAfter(first, second, dict), expectedIsAfter, "Wrong return value from IntervalUtils.isAfter()");

        if ( expectedContigComparison < 0 ) {
            Assert.assertTrue(IntervalUtils.compareContigs(first, second, dict) < 0, "Wrong return value from IntervalUtils.compareContigs()");
        } else if ( expectedContigComparison > 0 ) {
            Assert.assertTrue(IntervalUtils.compareContigs(first, second, dict) > 0, "Wrong return value from IntervalUtils.compareContigs()");
        } else {
            Assert.assertTrue(IntervalUtils.compareContigs(first, second, dict) == 0, "Wrong return value from IntervalUtils.compareContigs()");
        }
    }

    @Test
    public void testSpanningInterval_nullIfEmptyInput() throws Exception {
        Assert.assertNull(IntervalUtils.getSpanningInterval(Collections.emptyList()));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSpanningInterval_nullElement() throws Exception {
        final SimpleInterval chr1_200_300 = new SimpleInterval("1:200-300");
        Assert.assertNull(IntervalUtils.getSpanningInterval(Arrays.asList(chr1_200_300, null)));
    }

    @DataProvider(name = "SpanningInterval")
    public Object[][] SpanningInterval() {
        final SimpleInterval chr1_1_100   = new SimpleInterval("1:1-100");
        final SimpleInterval chr1_1_200 = new SimpleInterval("1:1-200");
        final SimpleInterval chr1_1_300 = new SimpleInterval("1:1-300");
        final SimpleInterval chr1_100_200 = new SimpleInterval("1:100-200");
        final SimpleInterval chr1_200_300 = new SimpleInterval("1:200-300");
        return new Object[][]{
           {Arrays.asList(chr1_1_100),                             chr1_1_100, },
           {Arrays.asList(chr1_1_100, chr1_100_200),               chr1_1_200, },
           {Arrays.asList(chr1_1_100, chr1_100_200, chr1_200_300), chr1_1_300, },
           {Arrays.asList(chr1_1_100, chr1_200_300),               chr1_1_300, },
        };
    }

    @Test(dataProvider = "SpanningInterval")
    public void testSpanningInterval(final List<? extends Locatable> locs, final SimpleInterval expectedResult) throws Exception {
        Assert.assertEquals(IntervalUtils.getSpanningInterval(locs), expectedResult);
    }

    @Test
    public void testSpanningIntervalForIntervalList() throws Exception {
        final List<SimpleInterval> exonsTwoContigs =
                Arrays.asList(new SimpleInterval("chr2:11719-12377"),
                                new SimpleInterval("chr2:12445-14651"),
                                new SimpleInterval("chr2:14855-15188"),
                                new SimpleInterval("chr2:15646-16097"),
                                new SimpleInterval("chr2:16457-18516"),
                                new SimpleInterval("chr10:38664-41779"),
                                new SimpleInterval("chr10:42659-43102"),
                                new SimpleInterval("chr10:45288-46537"),
                                new SimpleInterval("chr10:46657-47020"),
                                new SimpleInterval("chr10:197419-202755"));
        
        //spanning intervals should be one per contig, sorted by contig according to the dictionary (not lexicographically -- I'm looking at you, Bedtools)
        final List<SimpleInterval> expected =
                Arrays.asList(new SimpleInterval("chr2:11719-18516"), new SimpleInterval("chr10:38664-202755"));

        final SAMSequenceDictionary dictionary = SAMSequenceDictionaryExtractor.extractDictionary(new File(FULL_HG38_DICT).toPath());
        final List<SimpleInterval> spanning = IntervalUtils.getSpanningIntervals(exonsTwoContigs, dictionary);
        Assert.assertTrue(spanning.size() == expected.size());
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertTrue(spanning.get(i).equals(expected.get(i)));
        }
    }

    // -------------------------------------------------------------------------------------
    //
    // tests to ensure the quality of the interval cuts of the interval cutting functions
    //
    // -------------------------------------------------------------------------------------

    private class IntervalSlicingTest extends TestDataProvider {
        public int parts;
        public double maxAllowableVariance;

        private IntervalSlicingTest(final int parts, final double maxAllowableVariance) {
            super(IntervalSlicingTest.class);
            this.parts = parts;
            this.maxAllowableVariance = maxAllowableVariance;
        }

        public String toString() {
            return String.format("IntervalSlicingTest parts=%d maxVar=%.2f", parts, maxAllowableVariance);
        }
    }

    @DataProvider(name = "intervalslicingdata")
    public Object[][] createTrees() {
        new IntervalSlicingTest(1, 0);
        new IntervalSlicingTest(2, 1);
        new IntervalSlicingTest(5, 1);
        new IntervalSlicingTest(10, 1);
        new IntervalSlicingTest(67, 1);
        new IntervalSlicingTest(100, 1);

        return IntervalSlicingTest.getTests(IntervalSlicingTest.class);
    }

    @Test(dataProvider = "intervalslicingdata")
    public void testFixedScatterIntervalsAlgorithm(IntervalSlicingTest test) {
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(hg19exomeIntervals, test.parts);

        long totalSize = IntervalUtils.intervalSize(hg19exomeIntervals);
        long idealSplitSize = totalSize / test.parts;

        long sumOfSplitSizes = 0;
        int counter = 0;
        for ( final List<GenomeLoc> split : splits ) {
            long splitSize = IntervalUtils.intervalSize(split);
            double sigma = (splitSize - idealSplitSize) / (1.0 * idealSplitSize);
            //logger.warn(String.format("Split %d size %d ideal %d sigma %.2f", counter, splitSize, idealSplitSize, sigma));
            counter++;
            sumOfSplitSizes += splitSize;
            Assert.assertTrue(Math.abs(sigma) <= test.maxAllowableVariance, String.format("Interval %d (size %d ideal %d) has a variance %.2f outside of the tolerated range %.2f", counter, splitSize, idealSplitSize, sigma, test.maxAllowableVariance));
        }

        Assert.assertEquals(totalSize, sumOfSplitSizes, "Split intervals don't contain the exact number of bases in the origianl intervals");
    }

    // -------------------------------------------------------------------------------------
    //
    // splitLocusIntervals tests
    //
    // -------------------------------------------------------------------------------------

    /** large scale tests for many intervals */
    private class SplitLocusIntervalsTest extends TestDataProvider {
        final List<GenomeLoc> originalIntervals;
        final public int parts;

        private SplitLocusIntervalsTest(final String name, List<GenomeLoc> originalIntervals, final int parts) {
            super(SplitLocusIntervalsTest.class, name);
            this.parts = parts;
            this.originalIntervals = originalIntervals;
        }

        public String toString() {
            return String.format("%s parts=%d", super.toString(), parts);
        }
    }

    @DataProvider(name = "IntervalRepartitionTest")
    public Object[][] createIntervalRepartitionTest() {
        for ( int parts : Arrays.asList(1, 2, 3, 10, 13, 100, 151, 1000, 10000) ) {
            //for ( int parts : Arrays.asList(10) ) {
            new SplitLocusIntervalsTest("hg19RefLocs", hg19ReferenceLocs, parts);
            new SplitLocusIntervalsTest("hg19ExomeLocs", hg19exomeIntervals, parts);
        }

        return SplitLocusIntervalsTest.getTests(SplitLocusIntervalsTest.class);
    }

    @Test(dataProvider = "IntervalRepartitionTest")
    public void testIntervalRepartition(SplitLocusIntervalsTest test) {
        List<List<GenomeLoc>> splitByLocus = IntervalUtils.splitLocusIntervals(test.originalIntervals, test.parts);
        Assert.assertEquals(splitByLocus.size(), test.parts, "SplitLocusIntervals failed to generate correct number of intervals");
        List<GenomeLoc> flat = IntervalUtils.flattenSplitIntervals(splitByLocus);

        // test overall size
        final long originalSize = IntervalUtils.intervalSize(test.originalIntervals);
        final long flatSize = IntervalUtils.intervalSize(flat);
        Assert.assertEquals(flatSize, originalSize, "SplitLocusIntervals locs cover an incorrect number of bases");

        // test size of each split
        final long ideal = (long)Math.floor(originalSize / (1.0 * test.parts));
        final long maxSize = ideal + (originalSize % test.parts) * test.parts; // no more than N * rounding error in size
        for ( final List<GenomeLoc> split : splitByLocus ) {
            final long splitSize = IntervalUtils.intervalSize(split);
            Assert.assertTrue(splitSize >= ideal && splitSize <= maxSize,
                    String.format("SplitLocusIntervals interval (start=%s) has size %d outside of bounds ideal=%d, max=%d",
                            split.get(0), splitSize, ideal, maxSize));
        }

        // test that every base in original is covered once by a base in split by locus intervals
        String diff = IntervalUtils.equateIntervals(test.originalIntervals, flat);
        Assert.assertNull(diff, diff);
    }

    /** small scale tests where the expected cuts are enumerated upfront for testing */
    private class SplitLocusIntervalsSmallTest extends TestDataProvider {
        final List<GenomeLoc> original;
        final public int parts;
        final public int expectedParts;
        final List<GenomeLoc> expected;

        private SplitLocusIntervalsSmallTest(final String name, List<GenomeLoc> originalIntervals, final int parts, List<GenomeLoc> expected) {
            this(name, originalIntervals, parts,  expected, parts);
        }

        private SplitLocusIntervalsSmallTest(final String name, List<GenomeLoc> originalIntervals, final int parts, List<GenomeLoc> expected, int expectedParts) {
            super(SplitLocusIntervalsSmallTest.class, name);
            this.parts = parts;
            this.expectedParts = expectedParts;
            this.original = originalIntervals;
            this.expected = expected;
        }

        public String toString() {
            return String.format("%s parts=%d", super.toString(), parts);
        }
    }

    @DataProvider(name = "SplitLocusIntervalsSmallTest")
    public Object[][] createSplitLocusIntervalsSmallTest() {
        GenomeLoc bp01_10 = hg19GenomeLocParser.createGenomeLoc("1", 1, 10);

        GenomeLoc bp1_5 = hg19GenomeLocParser.createGenomeLoc("1", 1, 5);
        GenomeLoc bp6_10 = hg19GenomeLocParser.createGenomeLoc("1", 6, 10);
        new SplitLocusIntervalsSmallTest("cut into two", Arrays.asList(bp01_10), 2, Arrays.asList(bp1_5, bp6_10));

        GenomeLoc bp20_30 = hg19GenomeLocParser.createGenomeLoc("1", 20, 30);
        new SplitLocusIntervalsSmallTest("two in two", Arrays.asList(bp01_10, bp20_30), 2, Arrays.asList(bp01_10, bp20_30));

        GenomeLoc bp1_7 = hg19GenomeLocParser.createGenomeLoc("1", 1, 7);
        GenomeLoc bp8_10 = hg19GenomeLocParser.createGenomeLoc("1", 8, 10);
        GenomeLoc bp20_23 = hg19GenomeLocParser.createGenomeLoc("1", 20, 23);
        GenomeLoc bp24_30 = hg19GenomeLocParser.createGenomeLoc("1", 24, 30);
        new SplitLocusIntervalsSmallTest("two in three", Arrays.asList(bp01_10, bp20_30), 3,
                Arrays.asList(bp1_7, bp8_10, bp20_23, bp24_30));

        GenomeLoc bp1_2 = hg19GenomeLocParser.createGenomeLoc("1", 1, 2);
        GenomeLoc bp1_1 = hg19GenomeLocParser.createGenomeLoc("1", 1, 1);
        GenomeLoc bp2_2 = hg19GenomeLocParser.createGenomeLoc("1", 2, 2);
        new SplitLocusIntervalsSmallTest("too many pieces", Arrays.asList(bp1_2), 5, Arrays.asList(bp1_1, bp2_2), 2);

        new SplitLocusIntervalsSmallTest("emptyList", Collections.<GenomeLoc>emptyList(), 5, Collections.<GenomeLoc>emptyList(), 0);

        return SplitLocusIntervalsSmallTest.getTests(SplitLocusIntervalsSmallTest.class);
    }

    @Test(dataProvider = "SplitLocusIntervalsSmallTest")
    public void splitLocusIntervalsSmallTest(SplitLocusIntervalsSmallTest test) {
        List<List<GenomeLoc>> splitByLocus = IntervalUtils.splitLocusIntervals(test.original, test.parts);
        Assert.assertEquals(splitByLocus.size(), test.expectedParts, "SplitLocusIntervals failed to generate correct number of intervals");
        List<GenomeLoc> flat = IntervalUtils.flattenSplitIntervals(splitByLocus);

        // test sizes
        final long originalSize = IntervalUtils.intervalSize(test.original);
        final long splitSize = IntervalUtils.intervalSize(flat);
        Assert.assertEquals(splitSize, originalSize, "SplitLocusIntervals locs cover an incorrect number of bases");

        Assert.assertEquals(flat, test.expected, "SplitLocusIntervals locs not expected intervals");
    }

    //
    // Misc. tests
    //

    @Test(expectedExceptions=UserException.class)
    public void testMergeListsBySetOperatorNoOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> listEveryTwoFromOne = new ArrayList<>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(hg19GenomeLocParser.createGenomeLoc("chr1",x,x));
            else
                listEveryTwoFromOne.add(hg19GenomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 100);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, null);
        Assert.assertEquals(ret.size(), 100);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 0);
    }

    @Test
    public void testMergeListsBySetOperatorAllOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
            allSites.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 150);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, null);
        Assert.assertEquals(ret.size(), 150);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 50);
    }

    @Test
    public void testMergeListsBySetOperator() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 5 == 0) {
                listEveryTwoFromTwo.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
                allSites.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
            }
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, null);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 20);
    }

    @Test
    public void testOverlappingIntervalsFromSameSourceWithIntersection() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> source1 = new ArrayList<>();
        List<GenomeLoc> source2 = new ArrayList<>();

        source1.add(hg19GenomeLocParser.createGenomeLoc("1", 10, 20));
        source1.add(hg19GenomeLocParser.createGenomeLoc("1", 15, 25));

        source2.add(hg19GenomeLocParser.createGenomeLoc("1", 16, 18));
        source2.add(hg19GenomeLocParser.createGenomeLoc("1", 22, 24));

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(source1, source2, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 2);
    }

    @Test
    public void testGetContigLengths() {
        Map<String, Integer> lengths = IntervalUtils.getContigSizes(IOUtils.getPath(GATKBaseTest.exampleReference));
        Assert.assertEquals((long)lengths.get("1"), 16000);
        Assert.assertEquals((long)lengths.get("2"), 16000);
        Assert.assertEquals((long) lengths.get("3"), 16000);
        Assert.assertEquals((long)lengths.get("4"), 16000);
    }

    @Test
    public void testParseIntervalArguments() {
        Assert.assertEquals(hg19ReferenceLocs.size(), 4);
        Assert.assertEquals(intervalStringsToGenomeLocs("1", "2", "3").size(), 3);
        Assert.assertEquals(intervalStringsToGenomeLocs("1:1-2", "1:4-5", "2:1-1", "3:2-2").size(), 4);
    }

    private SAMSequenceDictionary getAmbiguousSequenceDictionary() {
        // The following strings represent a set of valid contig names, each of which can also be interpreted
        // as a valid interval query against some *other* contig in the same list (except the last one). That
        // is, each is comprised of a prefix that appears as another contig in the list, followed by a suffix
        // that is a valid query against that prefix.
        //
        // Given a sequence dictionary containing all of these contigs, each contig name, when used as an
        // interval query, should be rejected as ambiguous.
        List<String> ambiguousContigNames = Arrays.asList(
                "HLA-A*01:02:03:04",    // prefix: HLA-A*01:02:03       suffix: :04
                "HLA-A*01:02:03:04+",   // prefix: HLA-A*01:02:03       suffix: :04+
                "HLA-A*01:02:03:04:99", // prefix: HLA-A*01:02:03:04    suffix: :99
                "HLA-A*01:02:03:04-99", // prefix: HLA-A*01:02:03       suffix: :04-99
                "HLA-A*01:02:03",       // prefix: HLA-A*01:02          suffix: :03
                "HLA-A*01:02",          // prefix: HLA-A*01             suffix: :02

                "HLA-A*01"              // no valid prefix
        );

        final SAMSequenceDictionary ambiguousContigDictionary = new SAMSequenceDictionary();
        ambiguousContigNames.forEach(s -> ambiguousContigDictionary.addSequence(new SAMSequenceRecord(s, 99)));
        return ambiguousContigDictionary;
    }

    @DataProvider(name = "ambiguousIntervalQueries")
    public Object[][] getAmbiguousIntervalQueries() {
        return new Object[][]{
                // query interval string, expected ambiguous intervals

                // query intervals with more than one interpretation
                { "HLA-A*01:02:03:04",
                        Arrays.asList(
                            new SimpleInterval("HLA-A*01:02:03:04", 1, 99),
                            new SimpleInterval("HLA-A*01:02:03", 4, 4)) },
                { "HLA-A*01:02:03:04+",
                        Arrays.asList(
                                new SimpleInterval("HLA-A*01:02:03:04+", 1, 99),
                                new SimpleInterval("HLA-A*01:02:03", 4, 99)) },
                { "HLA-A*01:02:03:04:99",
                        Arrays.asList(
                                new SimpleInterval("HLA-A*01:02:03:04:99", 1, 99),
                                new SimpleInterval("HLA-A*01:02:03:04", 99, 99)) },
                { "HLA-A*01:02:03:04-99",
                        Arrays.asList(
                                new SimpleInterval("HLA-A*01:02:03:04-99", 1, 99),
                                new SimpleInterval("HLA-A*01:02:03", 4, 99)) },
                { "HLA-A*01:02:03",
                        Arrays.asList(
                                new SimpleInterval("HLA-A*01:02:03", 1, 99),
                                new SimpleInterval("HLA-A*01:02", 3, 3)) },
                { "HLA-A*01:02",
                        Arrays.asList(
                                new SimpleInterval("HLA-A*01:02", 1, 99),
                                new SimpleInterval("HLA-A*01", 2, 2)) },

                { "HLA:02", Collections.EMPTY_LIST},  // embedded :, but no valid prefix
                { "HLA-A*01",
                        Arrays.asList(new SimpleInterval("HLA-A*01", 1, 99)) }
        };
    }

    @Test(dataProvider = "ambiguousIntervalQueries")
    public void testAmbiguousIntervalQueries(final String ambiguousQuery, final List<String> expectedAmbiguousContigs) {
        final SAMSequenceDictionary sd = getAmbiguousSequenceDictionary();
        Assert.assertEquals(new HashSet<>(IntervalUtils.getResolvedIntervals(ambiguousQuery, sd)), new HashSet<>(expectedAmbiguousContigs));
    }

    @Test(expectedExceptions = NumberFormatException.class)
    public void testMalformedQueryFail() {
        final SAMSequenceDictionary sd = getAmbiguousSequenceDictionary();
        sd.addSequence(new SAMSequenceRecord("contig:1-10", 99));
        sd.addSequence(new SAMSequenceRecord("contig", 99));
        IntervalUtils.getResolvedIntervals("contig:abc-dce", sd);
    }

    @Test
    public void testMalformedAmbiguousQuerySucceeds() {
        final SAMSequenceDictionary sd = getAmbiguousSequenceDictionary();
        final String queryString = "contig:abc-dce";

        sd.addSequence(new SAMSequenceRecord("contig", 99));
        sd.addSequence(new SAMSequenceRecord(queryString, 99));

        // the query string matches as a full contig, but also has a prefix that matches the contig "contig". When
        // viewed as the latter, however, the start-end positions fail to parse as numbers. This resolves as a query
        // against the longer contig "contig:abc-dce", and logs a warning
        HashSet<SimpleInterval> expectedIntervals = new HashSet<>();
        expectedIntervals.add(new SimpleInterval(queryString, 1, 99));
        Assert.assertEquals(
                new HashSet<>(IntervalUtils.getResolvedIntervals(queryString, sd)),
                new HashSet<>(expectedIntervals));
    }

    @Test
    public void testUseAllHG38ContigsInIntervalQuery() {
        SAMSequenceDictionary sd;

        final File testFile = new File ("src/test/resources/org/broadinstitute/hellbender/engine/Homo_sapiens_assembly38.headerOnly.vcf.gz");
        try (VCFFileReader vcfReader = new VCFFileReader(testFile, false)) {
            sd = vcfReader.getFileHeader().getSequenceDictionary();
        }

        // Pull each contig from hg38 and validate that it can be used as a interval string without
        // being rejected as ambiguous.
        sd.getSequences().stream().forEach(
                hg38Contig -> {
                    // query using just the raw contig name
                    assertValidUniqueInterval(
                            sd,
                            hg38Contig.getSequenceName(),
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength())
                    );

                    // query using the format "name:1"
                    assertValidUniqueInterval(
                            sd,
                            hg38Contig.getSequenceName() + ":1",
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, 1)
                    );

                    // query using the format "name:1+"
                    assertValidUniqueInterval(
                            sd,
                            hg38Contig.getSequenceName() + ":1+",
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength())
                    );

                    // query using the format "name:1-end"
                    assertValidUniqueInterval(
                            sd,
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength()).toString(),
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength())
                    );
                }
        );
    }

    @Test(expectedExceptions = UserException.MalformedGenomeLoc.class)
    public void testReasonableErrorWithInvalidBedFile() throws IOException {
        File outOfBoundsBedFile = new File(INTERVAL_TEST_DATA + "bad_hg19_interval_span.bed"); // this bed file has a bogus interval the extends past the end chr1 and should throw an exception
        // this should throw a user exception because of the invalid interval
        IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, outOfBoundsBedFile.getAbsolutePath());
    }

    @Test(expectedExceptions = UserException.MalformedGenomeLoc.class)
    public void testReasonableErrorWithInvalidListFile() throws IOException {
        File outOfBoundListFile = new File(INTERVAL_TEST_DATA + "bad_hg19_example_intervals.list"); // this list file has a bogus interval the extends past the end chr1 and should throw an exception
        // this should throw a user exception because of the invalid interval
        IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, outOfBoundListFile.getAbsolutePath());
    }

    //NOTE the .interval_list style argument input gets handled differently from the other -L inputs but we still get a reasonable exception in this case
    @Test(expectedExceptions = UserException.MalformedFile.class)
    public void testReasonableErrorWithInvalidPicardListFile() throws IOException {
        File outOfBoundsIntervalListFile = new File(INTERVAL_TEST_DATA + "bad_hg19.interval_list"); // this list file has a bogus interval the extends past the end chr1 and should throw an exception
        // this should throw a user exception because of the invalid interval
        IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, outOfBoundsIntervalListFile.getAbsolutePath());
    }


    @Test
    public void testResolveAmbiguousIntervalQueryUsingBEDFormat() throws IOException {
        // For the rare case where an interval query is ambiguous, make sure it can be resolved by using a bed file.
        final SAMSequenceDictionary ambiguousSequenceDictionary = getAmbiguousSequenceDictionary();

        // use a test interval that has more than one interpretation
        // "HLA-A*01:02:03:04",
        //      -> new SimpleInterval("HLA-A*01:02:03:04", 1, 99)
        //      -> new SimpleInterval("HLA-A*01:02:03", 4, 4)
        final File bedFile = createTempFile("testAmbiguousInterval", ".bed");

        final int BED_START = 1;
        final int BED_END = 99;
        final String TARGET_CONTIG = "HLA-A*01:02:03:04";
        try (FileWriter bedWriter = new FileWriter(bedFile)) {
            bedWriter.write(String.format("%s\t%d\t%d", TARGET_CONTIG, BED_START, BED_END));
        }
        final GenomeLocParser genomeLocParser = new GenomeLocParser(ambiguousSequenceDictionary);

        // first, make sure we're using a query that is actually ambiguous
        Assert.assertThrows(UserException.MalformedGenomeLoc.class, () -> IntervalUtils.parseIntervalArguments(genomeLocParser, TARGET_CONTIG));

        // now use the bed file
        final List<GenomeLoc> intervalList = IntervalUtils.parseIntervalArguments(genomeLocParser, bedFile.getAbsolutePath());

        // validate the (translated from Bed) interval coordinates
        Assert.assertEquals(intervalList.size(), 1);
        Assert.assertEquals(intervalList.get(0).getContig(), TARGET_CONTIG);
        Assert.assertEquals(intervalList.get(0).getStart(), BED_START + 1);
        Assert.assertEquals(intervalList.get(0).getEnd(), BED_END);
    }

    private void assertValidUniqueInterval(
            final SAMSequenceDictionary sequenceDictionary,
            final String queryString,
            final SimpleInterval expectedInterval) {
        final List<SimpleInterval> ambiguousIntervals = IntervalUtils.getResolvedIntervals(queryString, sequenceDictionary);
        Assert.assertNotNull(ambiguousIntervals);
        Assert.assertEquals(ambiguousIntervals.size(), 1);
        Assert.assertEquals(ambiguousIntervals.get(0), expectedInterval);
    }

    @Test
    public void testTerminatingPlusSign() {
        // special case test for a query "prefix:+", where "prefix:" is a valid contig
        final SAMSequenceDictionary sd = getAmbiguousSequenceDictionary();
        final String contigPrefix = "prefix:";
        sd.addSequence(new SAMSequenceRecord(contigPrefix, 99));
        Assert.assertEquals(
                new HashSet<>(IntervalUtils.getResolvedIntervals(contigPrefix + "+", sd)),
                new HashSet<>());

        sd.addSequence(new SAMSequenceRecord(contigPrefix + "+", 99));
        HashSet<SimpleInterval> expectedIntervals = new HashSet<>();
        expectedIntervals.add(new SimpleInterval("prefix:+", 1, 99));
        Assert.assertEquals(
                new HashSet<>(IntervalUtils.getResolvedIntervals(contigPrefix + "+", sd)),
                expectedIntervals);
    }

    @Test
    public void testParseUnmappedIntervalArgument() {
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord("1", 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final GenomeLocParser parser = new GenomeLocParser(dictionary);

        List<GenomeLoc> unmappedLoc = IntervalUtils.parseIntervalArguments(parser, "unmapped");
        Assert.assertTrue(unmappedLoc.size() == 1);
        Assert.assertTrue(unmappedLoc.get(0).isUnmapped());
    }

    @Test
    public void testParseUnmappedIntervalArgumentInIntervalFile() {
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord("1", 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final GenomeLocParser parser = new GenomeLocParser(dictionary);

        List<GenomeLoc> unmappedLoc = IntervalUtils.parseIntervalArguments(parser, publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped.intervals");
        Assert.assertTrue(unmappedLoc.size() == 1);
        Assert.assertTrue(unmappedLoc.get(0).isUnmapped());
    }

    @Test
    public void testIsGatkIntervalFile() {
        Assert.assertTrue(IntervalUtils.isGatkIntervalFile(emptyIntervals));
        Assert.assertTrue(IntervalUtils.isGatkIntervalFile(emptyIntervals, true));

        Assert.assertFalse(IntervalUtils.isGatkIntervalFile(INTERVAL_TEST_DATA + "intervals_from_features_test.vcf"));
        Assert.assertFalse(IntervalUtils.isGatkIntervalFile(INTERVAL_TEST_DATA + "intervals_from_features_test.bed"));

        List<String> extensions = Arrays.asList("intervals", "list");
        for (String extension: extensions) {
            Assert.assertTrue(IntervalUtils.isGatkIntervalFile("test_intervals." + extension, false), "Tested interval file extension: " + extension);
        }
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testMissingIntervalFile() {
        IntervalUtils.isGatkIntervalFile(GATKBaseTest.publicTestDir + "no_such_intervals.list");
    }

    @Test
    public void testParseIntervalWithPeriodInContigName() {
        // Make sure that we don't interpret contigs with periods in their name as files
        final String contigName = "GL000249.1";
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord(contigName, 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final GenomeLocParser parser = new GenomeLocParser(dictionary);

        final List<GenomeLoc> result = IntervalUtils.parseIntervalArguments(parser, contigName);
        Assert.assertEquals(result.size(), 1);
        Assert.assertEquals(result.get(0).getContig(), contigName);
        Assert.assertEquals(result.get(0).getStart(), 1);
        Assert.assertEquals(result.get(0).getEnd(), 100);
    }

    @Test
    public void testFixedScatterIntervalsBasic() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");

        List<File> files = testFiles("basic.", 3, ".interval_list");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2", "3");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterFixedIntervalsLessFiles() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");
        GenomeLoc chr4 = hg19GenomeLocParser.parseGenomeLoc("4");

        List<File> files = testFiles("less.", 3, ".interval_list");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2", "3", "4");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
        Assert.assertEquals(locs3.get(1), chr4);
    }

    @Test(expectedExceptions=CommandLineException.BadArgumentValue.class)
    public void testSplitFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2");
        IntervalUtils.splitFixedIntervals(locs, files.size());
    }

    @Test(expectedExceptions=CommandLineException.BadArgumentValue.class)
    public void testScatterFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, locs.size()); // locs.size() instead of files.size()
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);
    }
    @Test
    public void testScatterFixedIntervalsStart() {
        List<String> intervals = Arrays.asList("1:1-2", "1:4-5", "2:1-1", "3:2-2");
        GenomeLoc chr1a = hg19GenomeLocParser.parseGenomeLoc("1:1-2");
        GenomeLoc chr1b = hg19GenomeLocParser.parseGenomeLoc("1:4-5");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".interval_list");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs2.get(0), chr1b);
        Assert.assertEquals(locs3.get(0), chr2);
        Assert.assertEquals(locs3.get(1), chr3);
    }

    @Test
    public void testScatterFixedIntervalsMiddle() {
        List<String> intervals = Arrays.asList("1:1-1", "2:1-2", "2:4-5", "3:2-2");
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2a = hg19GenomeLocParser.parseGenomeLoc("2:1-2");
        GenomeLoc chr2b = hg19GenomeLocParser.parseGenomeLoc("2:4-5");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".interval_list");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs3.get(0), chr2b);
        Assert.assertEquals(locs3.get(1), chr3);
    }

    @Test
    public void testScatterFixedIntervalsEnd() {
        List<String> intervals = Arrays.asList("1:1-1", "2:2-2", "3:1-2", "3:4-5");
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:2-2");
        GenomeLoc chr3a = hg19GenomeLocParser.parseGenomeLoc("3:1-2");
        GenomeLoc chr3b = hg19GenomeLocParser.parseGenomeLoc("3:4-5");

        List<File> files = testFiles("split.", 3, ".interval_list");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3a);
        Assert.assertEquals(locs3.get(0), chr3b);
    }

    @Test
    public void testScatterFixedIntervalsFile() {
        List<File> files = testFiles("sg.", 10, ".interval_list");
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(hg19MiniIntervalFile));
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());

        int[] counts = {
                33, 29, 68, 58, 64,
                31, 27, 62, 59, 52
        };

        //String splitCounts = "";
        for (int i = 0; i < splits.size(); i++) {
            int splitCount = splits.get(i).size();
            Assert.assertEquals(splitCount, counts[i], "Num.interval_list in split " + i);
        }
        //System.out.println(splitCounts.substring(2));

        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        int locIndex = 0;
        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), counts[i], "Intervals in " + file);
            for (GenomeLoc parsedLoc: parsedLocs)
                Assert.assertEquals(parsedLoc, locs.get(locIndex), String.format("Genome loc %d from file %d", locIndex++, i));
        }
        Assert.assertEquals(locIndex, locs.size(), "Total number of GenomeLocs");
    }

    @Test
    public void testScatterFixedIntervalsMax() {
        List<File> files = testFiles("sg.", 4, ".interval_list");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(hg19ReferenceLocs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    @Test
    public void testScatterContigIntervalsOrder() {
        List<String> intervals = Arrays.asList("2:1-1", "1:1-1", "3:2-2");
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".interval_list");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr2);
        Assert.assertEquals(locs2.get(0), chr1);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsBasic() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");

        List<File> files = testFiles("contig_basic.", 3, ".interval_list");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs("1", "2", "3"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsLessFiles() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");
        GenomeLoc chr4 = hg19GenomeLocParser.parseGenomeLoc("4");

        List<File> files = testFiles("contig_less.", 3, ".interval_list");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs("1", "2", "3", "4"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3);
        Assert.assertEquals(locs3.get(0), chr4);
    }

    @Test(expectedExceptions=UserException.BadInput.class)
    public void testScatterContigIntervalsMoreFiles() {
        List<File> files = testFiles("contig_more.", 3, ".intervals");
        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs("1", "2"), files);
    }

    @Test
    public void testScatterContigIntervalsStart() {
        List<String> intervals = Arrays.asList("1:1-2", "1:4-5", "2:1-1", "3:2-2");
        GenomeLoc chr1a = hg19GenomeLocParser.parseGenomeLoc("1:1-2");
        GenomeLoc chr1b = hg19GenomeLocParser.parseGenomeLoc("1:4-5");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("contig_split_start.", 3, ".interval_list");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs1.get(1), chr1b);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsMiddle() {
        List<String> intervals = Arrays.asList("1:1-1", "2:1-2", "2:4-5", "3:2-2");
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2a = hg19GenomeLocParser.parseGenomeLoc("2:1-2");
        GenomeLoc chr2b = hg19GenomeLocParser.parseGenomeLoc("2:4-5");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("contig_split_middle.", 3, ".interval_list");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 2);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs2.get(1), chr2b);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsEnd() {
        List<String> intervals = Arrays.asList("1:1-1", "2:2-2", "3:1-2", "3:4-5");
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:2-2");
        GenomeLoc chr3a = hg19GenomeLocParser.parseGenomeLoc("3:1-2");
        GenomeLoc chr3b = hg19GenomeLocParser.parseGenomeLoc("3:4-5");

        List<File> files = testFiles("contig_split_end.", 3 ,".interval_list");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3a);
        Assert.assertEquals(locs3.get(1), chr3b);
    }

    @Test
    public void testScatterContigIntervalsMax() {
        List<File> files = testFiles("sg.", 4, ".interval_list");
        IntervalUtils.scatterContigIntervals(hg19Header, hg19ReferenceLocs, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    private List<File> testFiles(String prefix, int count, String suffix) {
        ArrayList<File> files = new ArrayList<>();
        for (int i = 1; i <= count; i++) {
            files.add(createTempFile(prefix + i, suffix));
        }
        return files;
    }

    @DataProvider(name="unmergedIntervals")
    public Object[][] getUnmergedIntervals() {
        return new Object[][] {
                new Object[] {"small_unmerged_picard_intervals.interval_list"},
                new Object[] {"small_unmerged_gatk_intervals.list"}
        };
    }

    @Test(dataProvider="unmergedIntervals")
    public void testUnmergedIntervals(String unmergedIntervals) {
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Collections.singletonList(publicTestDir + unmergedIntervals));
        Assert.assertEquals(locs.size(), 2);

        List<GenomeLoc> merged;

        merged = IntervalUtils.mergeIntervalLocations(locs, IntervalMergingRule.ALL);
        Assert.assertEquals(merged.size(), 1);

        // Test that null means the same as ALL
        merged = IntervalUtils.mergeIntervalLocations(locs, null);
        Assert.assertEquals(merged.size(), 1);
    }

    @Test
    public void testIntervalFileIntervalsAreTrimmed() throws IOException {
        final GenomeLoc expectedGenomeLoc = hg19GenomeLocParser.createGenomeLoc("1", 1, 10);
        final File gatkIntervalFile = IOUtils.writeTempFile(
                Arrays.asList(
                    String.format(
                       // add some spaces around the start/end to test that it gets trimmed down to something valid
                        " %s:%d-%d ",
                        expectedGenomeLoc.getContig(),
                        expectedGenomeLoc.getStart(),
                        expectedGenomeLoc.getEnd())),
                "testIntervalFileIntervalsAreTrimmed", ".intervals.list");

        final List<String> intervalArgs = new ArrayList<>(1);
        intervalArgs.add(gatkIntervalFile.getAbsolutePath());

        final GenomeLocSortedSet locs = IntervalUtils.loadIntervals(
                intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
        Assert.assertEquals(locs.size(), 1);
        Assert.assertEquals(locs.toList().get(0), expectedGenomeLoc);
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class, dataProvider="invalidIntervalTestData")
    public void testInvalidGATKFileIntervalHandling(GenomeLocParser genomeLocParser,
                                                    String contig, int intervalStart, int intervalEnd ) throws Exception {

        File gatkIntervalFile = IOUtils.writeTempFile(
                Arrays.asList(String.format("%s:%d-%d", contig, intervalStart, intervalEnd)),
                "testInvalidGATKFileIntervalHandling", ".intervals");

        List<String> intervalArgs = new ArrayList<>(1);
        intervalArgs.add(gatkIntervalFile.getAbsolutePath());

        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    @Test(expectedExceptions={UserException.MalformedFile.class, UserException.MalformedGenomeLoc.class}, dataProvider="invalidIntervalTestData")
    public void testLoadIntervalsInvalidPicardIntervalHandling(GenomeLocParser genomeLocParser,
                                                  String contig, int intervalStart, int intervalEnd ) throws Exception {

        SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("1"));
        // Intentionally add an extra contig not in our genomeLocParser's sequence dictionary
        // so that we can add intervals to the Picard interval file for contigs not in our dictionary
        picardFileHeader.addSequence(new SAMSequenceRecord("2", 100000));

        IntervalList picardIntervals = new IntervalList(picardFileHeader);
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname"));

        File picardIntervalFile = createTempFile("testLoadIntervalsInvalidPicardIntervalHandling", ".intervals");
        picardIntervals.write(picardIntervalFile);

        List<String> intervalArgs = new ArrayList<>(1);
        intervalArgs.add(picardIntervalFile.getAbsolutePath());

        // loadIntervals() will validate all intervals against the sequence dictionary in our genomeLocParser,
        // and should throw for all intervals in our invalidIntervalTestData set
        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    @Test(expectedExceptions={UserException.CouldNotReadInputFile.class, UserException.MalformedGenomeLoc.class, UserException.MalformedFile.class}, dataProvider="invalidIntervalTestData")
    public void testIntervalFileToListInvalidPicardIntervalHandling(GenomeLocParser genomeLocParser,
                                       String contig, int intervalStart, int intervalEnd ) throws Exception {

        final SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("1"));
        picardFileHeader.addSequence(new SAMSequenceRecord("2", 100000));

        final IntervalList picardIntervals = new IntervalList(picardFileHeader);
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname"));

        final File picardIntervalFile = createTempFile("testIntervalFileToListInvalidPicardIntervalHandling", ".interval_list");
        picardIntervals.write(picardIntervalFile);

        // loadIntervals() will validate all intervals against the sequence dictionary in our genomeLocParser,
        // and should throw for all intervals in our invalidIntervalTestData set
        IntervalUtils.parseIntervalArguments(genomeLocParser, picardIntervalFile.getAbsolutePath());
    }

    @DataProvider(name="invalidIntervalTestData")
    public Object[][] invalidIntervalDataProvider() throws Exception {
        File fastaFile = new File(publicTestDir + "exampleFASTA.fasta");
        GenomeLocParser genomeLocParser = new GenomeLocParser(new IndexedFastaSequenceFile(fastaFile));

        return new Object[][] {
                new Object[] {genomeLocParser, "1", 10000000, 20000000},
                new Object[] {genomeLocParser, "2", 1, 2},
                new Object[] {genomeLocParser, "1", -1, 50}
        };
    }

    // TODO - remove once a corrected version of the exome interval list is released.
    @DataProvider(name="negativeOneLengthIntervalTestData")
    public Object[][] negativeOneLengthIntervalDataProvider() throws Exception {
        File fastaFile = new File(publicTestDir + "exampleFASTA.fasta");
        GenomeLocParser genomeLocParser = new GenomeLocParser(new IndexedFastaSequenceFile(fastaFile));

        return new Object[][] {
                new Object[] {genomeLocParser, "1", 2, 1},
                new Object[] {genomeLocParser, "1", 500, 499},
                new Object[] {genomeLocParser, "1", 1000, 999}
        };
    }

    @DataProvider(name = "sortAndMergeIntervals")
    public Object[][] getSortAndMergeIntervals() {
        return new Object[][] {
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, intervalStringsToGenomeLocs("1:1", "1:3", "1:2"), intervalStringsToGenomeLocs("1:1", "1:2", "1:3") },
                new Object[] { IntervalMergingRule.ALL, intervalStringsToGenomeLocs("1:1", "1:3", "1:2"), intervalStringsToGenomeLocs("1:1-3") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, intervalStringsToGenomeLocs("1:1", "1:3", "2:2"), intervalStringsToGenomeLocs("1:1", "1:3", "2:2") },
                new Object[] { IntervalMergingRule.ALL, intervalStringsToGenomeLocs("1:1", "1:3", "2:2"), intervalStringsToGenomeLocs("1:1", "1:3", "2:2") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, intervalStringsToGenomeLocs("1:1", "1"), intervalStringsToGenomeLocs("1") },
                new Object[] { IntervalMergingRule.ALL, intervalStringsToGenomeLocs("1:1", "1"), intervalStringsToGenomeLocs("1") }
        };
    }

    @Test(dataProvider = "sortAndMergeIntervals")
    public void testSortAndMergeIntervals(IntervalMergingRule merge, List<GenomeLoc> unsorted, List<GenomeLoc> expected) {
        List<GenomeLoc> sorted = IntervalUtils.sortAndMergeIntervals(hg19GenomeLocParser, unsorted, merge).toList();
        Assert.assertEquals(sorted, expected);
    }

    @DataProvider(name = "testSortAndMergeIntervalsFromSameContigData")
    public Object[][] testSortAndMergeIntervalsFromSameContigData() {
        return new Object[][] {
                // All overlapping
                { Arrays.asList(new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 3, 7),
                        new SimpleInterval("1", 5, 10)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 10))
                },

                // Some overlapping
                { Arrays.asList(new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 3, 7),
                        new SimpleInterval("1", 9, 10)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 7),
                                new SimpleInterval("1", 9, 10))
                },

                // Out of order
                { Arrays.asList(new SimpleInterval("1", 100, 110),
                        new SimpleInterval("1", 105, 115),
                        new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 3, 7),
                        new SimpleInterval("1", 5, 10)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 10),
                                new SimpleInterval("1", 100, 115))
                },

                // Adjacent, IntervalMergingRule.ALL
                { Arrays.asList(new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 6, 10),
                        new SimpleInterval("1", 11, 15)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 15))
                },

                // Adjacent, IntervalMergingRule.OVERLAPPING_ONLY
                { Arrays.asList(new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 6, 10),
                        new SimpleInterval("1", 11, 15)),
                        IntervalMergingRule.OVERLAPPING_ONLY,
                        Arrays.asList(new SimpleInterval("1", 1, 5),
                                new SimpleInterval("1", 6, 10),
                                new SimpleInterval("1", 11, 15)),
                },

                // Mix of overlapping and adjacent, IntervalMergingRule.ALL
                { Arrays.asList(new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 6, 10),
                        new SimpleInterval("1", 8, 15),
                        new SimpleInterval("1", 100, 200),
                        new SimpleInterval("1", 150, 250),
                        new SimpleInterval("1", 175, 275),
                        new SimpleInterval("1", 276, 300)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 15),
                                new SimpleInterval("1", 100, 300))
                },

                // Mix of overlapping and adjacent, IntervalMergingRule.OVERLAPPING_ONLY
                { Arrays.asList(new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 6, 10),
                        new SimpleInterval("1", 8, 15),
                        new SimpleInterval("1", 100, 200),
                        new SimpleInterval("1", 150, 250),
                        new SimpleInterval("1", 175, 275),
                        new SimpleInterval("1", 276, 300)),
                        IntervalMergingRule.OVERLAPPING_ONLY,
                        Arrays.asList(new SimpleInterval("1", 1, 5),
                                new SimpleInterval("1", 6, 15),
                                new SimpleInterval("1", 100, 275),
                                new SimpleInterval("1", 276, 300))
                },

                // Mix of overlapping and adjacent, IntervalMergingRule.ALL, out of order
                { Arrays.asList(new SimpleInterval("1", 100, 200),
                        new SimpleInterval("1", 8, 15),
                        new SimpleInterval("1", 175, 275),
                        new SimpleInterval("1", 1, 5),
                        new SimpleInterval("1", 6, 10),
                        new SimpleInterval("1", 276, 300),
                        new SimpleInterval("1", 150, 250)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 15),
                                new SimpleInterval("1", 100, 300))
                },

                // Empty list
                { Collections.emptyList(), IntervalMergingRule.ALL, Collections.emptyList() },

                // List of size 1
                { Arrays.asList(new SimpleInterval("1", 1, 5)),
                        IntervalMergingRule.ALL,
                        Arrays.asList(new SimpleInterval("1", 1, 5))
                }
        };
    }

    @Test(dataProvider = "testSortAndMergeIntervalsFromSameContigData")
    public void testSortAndMergeIntervalsFromSameContig(final List<SimpleInterval> unmergedIntervals,
                                                        final IntervalMergingRule mergingRule,
                                                        final List<SimpleInterval> expectedMergedIntervals) {

        final List<SimpleInterval> actualMergedIntervals = IntervalUtils.sortAndMergeIntervalsFromSameContig(unmergedIntervals, mergingRule);
        Assert.assertEquals(actualMergedIntervals, expectedMergedIntervals, "Wrong intervals after sortAndMergeIntervals()");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSortAndMergeIntervalsFromSameContigRejectsMultipleContigs() {
        final List<SimpleInterval> multiContigIntervals = Arrays.asList(new SimpleInterval("1", 1, 10), new SimpleInterval("1", 5, 15), new SimpleInterval("2", 1, 10));
        IntervalUtils.sortAndMergeIntervalsFromSameContig(multiContigIntervals, IntervalMergingRule.ALL);
    }

    @DataProvider(name="loadIntervalsNonMerge")
    public Object[][] getloadIntervalsNonMerge(){
        final String severalIntervals = "src/test/resources/org/broadinstitute/hellbender/utils/interval/example_intervals.list";

        return new Object[][]{
                new Object[]{Arrays.asList("1:5-5","1:6","1:7","1:8","1:9","1:10"), 0, intervalStringsToGenomeLocs("1:5-5","1:6","1:7","1:8","1:9","1:10")},
                new Object[]{Arrays.asList("1:1-2", "2:1-2"), 0, intervalStringsToGenomeLocs("1:1-2", "2:1-2")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), 0, intervalStringsToGenomeLocs("1:1-10", "1:5-15")},
                new Object[]{Arrays.asList("1:5-5"), 5, intervalStringsToGenomeLocs("1:1-10")},
                new Object[]{Arrays.asList(severalIntervals), 0, intervalStringsToGenomeLocs("1:100-200", "2:20-30", "4:50")}
        };
    }

    @Test( dataProvider = "loadIntervalsNonMerge")
    public void testLoadIntervalsNonMerging(List<String> intervals, int padding, List<GenomeLoc> results) {
        List<GenomeLoc> loadedIntervals = IntervalUtils.loadIntervalsNonMerging(intervals, padding, hg19GenomeLocParser);
        Assert.assertEquals(loadedIntervals, results);
    }

    @DataProvider(name="loadintervals")
    public Object[][] getloadIntervals(){
        final String severalIntervals = "src/test/resources/org/broadinstitute/hellbender/utils/interval/example_intervals.list";

        return new Object[][]{
                new Object[]{Arrays.asList("1:1-2"), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:1-2")},
                new Object[]{Arrays.asList("1:1-2", "2:1-2"), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:1-2", "2:1-2")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:1-15")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), IntervalSetRule.INTERSECTION, 0, intervalStringsToGenomeLocs("1:5-10")},
                new Object[]{Arrays.asList("1:5-5"), IntervalSetRule.UNION, 5, intervalStringsToGenomeLocs("1:1-10")},
                new Object[]{Arrays.asList(severalIntervals), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:100-200", "2:20-30", "4:50")}
        };
    }

    @Test( dataProvider = "loadintervals")
    public void loadIntervalsTest(List<String> intervals, IntervalSetRule setRule, int padding, List<GenomeLoc> results){
        GenomeLocSortedSet loadedIntervals = IntervalUtils.loadIntervals(intervals, setRule, IntervalMergingRule.ALL, padding, hg19GenomeLocParser);
        Assert.assertEquals(loadedIntervals, results);
    }

    @DataProvider(name = "loadIntervalsFromFeatureFileData")
    public Object[][] loadIntervalsFromFeatureFileData() {
        final List<GenomeLoc> expectedIntervals = intervalStringsToGenomeLocs("1:100-100", "1:203-206", "1:1000-1003", "2:200-200", "2:548-550", "4:776-779");

        return new Object[][] {
                { new File(INTERVAL_TEST_DATA + "intervals_from_features_test.vcf"), expectedIntervals },
                { new File(INTERVAL_TEST_DATA + "intervals_from_features_test.bed"), expectedIntervals }
        };
    }

    @Test(dataProvider = "loadIntervalsFromFeatureFileData")
    public void testLoadIntervalsFromFeatureFile( final File featureFile, final List<GenomeLoc> expectedIntervals ) {
        final GenomeLocSortedSet actualIntervals = IntervalUtils.loadIntervals(Collections.singletonList(featureFile.getAbsolutePath()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
        Assert.assertEquals(actualIntervals, expectedIntervals, "Wrong intervals loaded from Feature file " + featureFile.getAbsolutePath());
    }

    @Test(dataProvider = "loadIntervalsFromFeatureFileData")
    public void testLoadIntervalsFromFeatureFileInJimfs( final File featureFile, final List<GenomeLoc> expectedIntervals ) throws IOException {
        try(final FileSystem fs = Jimfs.newFileSystem(Configuration.unix())){
            final Path jimfsRootPath = fs.getRootDirectories().iterator().next();
            final Path jimfsCopy = Files.copy(featureFile.toPath(), jimfsRootPath.resolve(featureFile.getName()));
            final String jimfsPathString = jimfsCopy.toAbsolutePath().toUri().toString();
            final GenomeLocSortedSet actualIntervals = IntervalUtils.loadIntervals(Collections.singletonList(jimfsPathString), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
            Assert.assertEquals(actualIntervals, expectedIntervals, "Wrong intervals loaded from Feature file " + jimfsPathString);
        }
    }

    // Note: because the file does not exist and all characters are allowed in contig names,
    // we will not know that this is supposed to be interpreted as a file.
    // So we'll blow up with MalformedGenomeLoc and not anything related to files
    @Test(expectedExceptions = UserException.MalformedGenomeLoc.class)
    public void testLoadIntervalsFromNonExistentFile() {
        IntervalUtils.loadIntervals(Arrays.asList(GATKBaseTest.getSafeNonExistentFile("non_existent_file.vcf").getAbsolutePath()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testLoadIntervalsFromUnrecognizedFormatFile() {
        IntervalUtils.loadIntervals(Arrays.asList(INTERVAL_TEST_DATA + "unrecognized_format_file.xyz"), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
    }

    @Test(expectedExceptions = UserException.MalformedFile.class)
    public void loadIntervalsEmptyFile(){
        IntervalUtils.loadIntervals(Arrays.asList(emptyIntervals), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData",expectedExceptions = IllegalArgumentException.class)
    public void testGenomeLocsFromLocatablesNullParser(final GenomeLocParser parser,
                                                       final List<? extends Locatable> locatables) {
        IntervalUtils.genomeLocsFromLocatables(null, locatables);
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData",expectedExceptions = IllegalArgumentException.class)
    public void testGenomeLocsFromLocatablesNullLocatables(final GenomeLocParser parser,
                                                       final List<? extends Locatable> locatables) {
        IntervalUtils.genomeLocsFromLocatables(parser,null);
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData",expectedExceptions = IllegalArgumentException.class)
    public void testGenomeLocsFromLocatablesNullContainingLocatables(final GenomeLocParser parser,
                                                           final List<? extends Locatable> locatables) {
        if (locatables.size() == 0) {
            IntervalUtils.genomeLocsFromLocatables(parser,null);
        } else {
            final List<? extends Locatable> withNull = new ArrayList<>(locatables);
            withNull.set(withNull.size() / 2,null);
            IntervalUtils.genomeLocsFromLocatables(parser,withNull);
        }
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData")
    public void testGenomeLocsFromLocatables(final GenomeLocParser parser, final List<? extends Locatable> locatables) {
        final List<GenomeLoc> result = IntervalUtils.genomeLocsFromLocatables(parser, locatables);
        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), locatables.size());
        for (int i = 0; i < result.size(); i++) {
            final GenomeLoc resultLoc = result.get(i);
            final Locatable inputLoc = locatables.get(i);
            Assert.assertEquals(resultLoc.getContig(),inputLoc.getContig());
            Assert.assertEquals(resultLoc.getStart(),inputLoc.getStart());
            Assert.assertEquals(resultLoc.getStop(),inputLoc.getEnd());
        }
        // is it real only:
        try {
            result.add(parser.createGenomeLoc("1",1,2));
            Assert.fail("the result collection should not allow to call add");
        } catch (final UnsupportedOperationException ex) {
            // ok.
        }
        if (result.size() > 0) {
            try {
                result.set(0,parser.createGenomeLoc("1",1,2));
                Assert.fail("the result collection should not allow to call set");
            } catch (final UnsupportedOperationException ex) {
                // ok.
            }
        }
    }

    @DataProvider(name="genomeLocsFromLocatablesData")
    public Object[][] genomeLocsFromLocatablesData() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
        final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());

        return new Object[][] {
                { genomeLocParser, Collections.emptyList() },
                { genomeLocParser, Collections.singletonList(new SimpleInterval("1",1,2)) },
                { genomeLocParser, Arrays.asList(new SimpleInterval("1", 1, 2), new SimpleInterval("1", 1, 1), new SimpleInterval("1", 5, 9)) },
                { genomeLocParser, Arrays.asList(new SimpleInterval("1",9,10),new SimpleInterval("1",5,9), new SimpleInterval("1",1,9)) }
        };
    }

    @Test(dataProvider = "overlapData")
    public void testOverlap(final SimpleInterval l, final SimpleInterval r, final boolean expected) {
        Assert.assertEquals(IntervalUtils.overlaps(l, r), expected);
        Assert.assertEquals(IntervalUtils.overlaps(
                new SimpleFeature(l.getContig(),l.getStart(),l.getEnd()),
                new SimpleFeature(r.getContig(),r.getStart(),r.getEnd())),expected);
        Assert.assertEquals(IntervalUtils.overlaps(
                new SimpleFeature(l.getContig(),l.getStart(),l.getEnd()),
                r),expected);
        Assert.assertEquals(IntervalUtils.overlaps(
                l,
                new SimpleFeature(r.getContig(),r.getStart(),r.getEnd())),expected);
        Assert.assertFalse(IntervalUtils.overlaps(
                new SimpleFeature(null, r.getStart(), r.getEnd()),
                new SimpleFeature(null, l.getStart(), l.getEnd())));
    }

    @DataProvider(name = "overlapData")
    public Object[][] overlapData() {
        // removing cases where the second interval is null.
        return Stream.of(SimpleIntervalUnitTest.getIntervalOverlapData()).filter(a -> a[1] != null)
                .collect(Collectors.toList()).toArray(new Object[0][]);
    }

    @Test(dataProvider = "lexicographicalOrderComparatorData")
    public void testLexicographicalOrderComparator(final List<Locatable> unsorted) {
        final List<Locatable> sorted = unsorted.stream().sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        for (int i = 0; i < sorted.size() - 1; i++) {
           final Locatable thisLoc = sorted.get(i);
           final Locatable nextLoc = sorted.get(i + 1);
           if (thisLoc.getContig() == nextLoc.getContig()) {
               Assert.assertTrue(thisLoc.getStart() <= nextLoc.getStart());
               if (thisLoc.getStart() == nextLoc.getStart()) {
                   Assert.assertTrue(thisLoc.getEnd() <= nextLoc.getEnd());
               }
           } else if (thisLoc.getContig() == null) {
               Assert.fail("Null contig must go last");
           } else if (nextLoc.getContig() != null) {
               Assert.assertTrue(thisLoc.getContig().compareTo(nextLoc.getContig()) <= 0);
               if (thisLoc.getContig().equals(nextLoc.getContig())) {
                   Assert.assertTrue(thisLoc.getStart() <= nextLoc.getStart());
                   if (thisLoc.getStart() == nextLoc.getStart()) {
                       Assert.assertTrue(thisLoc.getEnd() <= nextLoc.getEnd());
                   }
               }
           }
        }
    }

    @DataProvider(name = "lexicographicalOrderComparatorData")
    public Object[][] lexicographicalOrderComparatorData() {
        final String[] CONTIG_NAMES = new String[] {"A","AA","B","C", null};
        final int[] CONTIG_SIZES = new int[] {200,300,400,500,600};
        final int MIN_INTERVAL_SIZE = 1;
        final int MAX_INTERVAL_SIZE = 100;
        final Random rdn = new Random(1131312131);
        final int CASE_NUMBER = 100;
        final List<Object[]> result = new ArrayList<>(100);
        for (int i = 0; i < CASE_NUMBER; i++) {
            final int locatableCount = rdn.nextInt(100) + 1;
            final List<Locatable> locatables = new ArrayList<>(locatableCount);
            for (int j = 0; j < locatableCount; j++) {
                final int contigIdx = rdn.nextInt(CONTIG_NAMES.length);
                final String contig = CONTIG_NAMES[contigIdx];

                final boolean useSimpleInterval = contig == null ? false : rdn.nextBoolean();
                final int intervalSize = rdn.nextInt(MAX_INTERVAL_SIZE - MIN_INTERVAL_SIZE + 1) + MIN_INTERVAL_SIZE;
                final int start = rdn.nextInt(CONTIG_SIZES[contigIdx] - intervalSize) + 1;
                final int end = start + intervalSize - 1;
                final Locatable loc = useSimpleInterval ? new SimpleInterval(contig,start,end) : new SimpleFeature(contig,start,end);
                locatables.add(loc);
            }
            result.add(new Object[] { locatables });
        }
        return result.toArray(new Object[result.size()][]);
    }


    @DataProvider(name="intervals")
    public Object[][] intervals(){
        ArrayList<SimpleInterval> input = Lists.newArrayList(new SimpleInterval("1", 10, 100));
        ArrayList<SimpleInterval> output = input;
        ArrayList<SimpleInterval> hundred = Lists.newArrayList(new SimpleInterval("1", 1, 100));

        return new Object[][]{
                // input fits all in one shard
                new Object[]{input, 1000, output},
                // input could fit in the shard, but is cut anyways so it's aligned to integer multiples of shard boundary
                new Object[]{input, 90, Lists.newArrayList(new SimpleInterval("1", 10, 90),new SimpleInterval("1", 91, 100))},
                // shard ends just after input
                new Object[]{hundred, 101, hundred},
                // shard ends at input
                new Object[]{hundred, 100, hundred},
                // shard ends just before input does
                new Object[]{hundred, 99, Lists.newArrayList(new SimpleInterval("1", 1, 99),new SimpleInterval("1", 100, 100))},
                // input overlaps three shards
                new Object[]{hundred, 40, Lists.newArrayList(new SimpleInterval("1", 1, 40),new SimpleInterval("1", 41, 80),new SimpleInterval("1", 81, 100))},
                // input covers four shards exactly
                new Object[]{hundred, 25, Lists.newArrayList(
                        new SimpleInterval("1", 1, 25),new SimpleInterval("1", 26, 50),new SimpleInterval("1", 51, 75),new SimpleInterval("1", 76, 100))},
        };
    }

    @Test(dataProvider = "intervals")
    public void testCutToShards(ArrayList<SimpleInterval> input, int shardSize, ArrayList<SimpleInterval> expected) throws Exception {
        List<SimpleInterval> actual = IntervalUtils.cutToShards(input, shardSize);
        Assert.assertEquals(
                actual,
                expected
        );
    }

    @DataProvider(name="shardIndex")
    public Object[][] shardIndex(){
        return new Object[][]{
                // offset, shardSize, expected shardIndex
                new Object[]{1,10,0},
                new Object[]{10,10,0},
                new Object[]{11,10,1},
                new Object[]{20,10,1},
                new Object[]{21,10,2}
        };
    }

    @Test(dataProvider = "shardIndex")
    public void testShardIndex(int offset, int shardSize, int shardIndex) throws Exception {
        Assert.assertEquals(
                IntervalUtils.shardIndex(offset, shardSize),
                shardIndex
        );
    }

    @DataProvider(name="shardBegin")
    public Object[][] shardBegin(){
        return new Object[][]{
                // shard index, shard size, shard begin
                new Object[]{0,10,1},
                new Object[]{1,10,11},
                new Object[]{1,42,43},
                new Object[]{2,10,21},
                new Object[]{2,42,85}
        };
    }

    @Test(dataProvider = "shardBegin")
    public void testBeginOfShard(int index, int shardSize, int shardBegin) throws Exception {
        Assert.assertEquals(
                IntervalUtils.beginOfShard(index, shardSize),
                shardBegin
        );
    }

    @DataProvider(name="shardEnd")
    public Object[][] shardEnd(){
        return new Object[][]{
                // shard index, shard size, shard end
                new Object[]{0,10,10},
                new Object[]{1,10,20},
                new Object[]{2,10,30},
                new Object[]{0,42,42},
                new Object[]{1,42,84},
        };
    }

    @Test(dataProvider = "shardEnd")
    public void testEndOfShard(int index, int shardSize, int shardEnd) throws Exception {
        Assert.assertEquals(
                IntervalUtils.endOfShard(index, shardSize),
                shardEnd
        );
    }

    @Test
    public void testConvertSimpleIntervalToQueryInterval() {
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord("1", 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final SimpleInterval originalInterval = new SimpleInterval("1", 5, 10);

        final QueryInterval convertedInterval = IntervalUtils.convertSimpleIntervalToQueryInterval(originalInterval, dictionary);
        Assert.assertEquals(convertedInterval.referenceIndex, 0);
        Assert.assertEquals(convertedInterval.start, 5);
        Assert.assertEquals(convertedInterval.end, 10);
    }

    @DataProvider(name = "TrimIntervalToContigData")
    public Object[][] trimIntervalToContigData() {
        return new Object[][] {
                { "1", 5, 10, 100, new SimpleInterval("1", 5, 10) },
                { "1", 1, 100, 100, new SimpleInterval("1", 1, 100) },
                { "1", -1, 100, 100, new SimpleInterval("1", 1, 100) },
                { "1", -5, 100, 100, new SimpleInterval("1", 1, 100) },
                { "1", -5, 90, 100, new SimpleInterval("1", 1, 90) },
                { "1", 1, 101, 100, new SimpleInterval("1", 1, 100) },
                { "1", 1, 105, 100, new SimpleInterval("1", 1, 100) },
                { "1", 5, 105, 100, new SimpleInterval("1", 5, 100) },
                { "1", -1, 10, 100, new SimpleInterval("1", 1, 10) },
                { "1", -5, 10, 100, new SimpleInterval("1", 1, 10) },
                { "1", 90, 101, 100, new SimpleInterval("1", 90, 100) },
                { "1", 90, 105, 100, new SimpleInterval("1", 90, 100) }
        };
    }

    @Test(dataProvider = "TrimIntervalToContigData")
    public void testTrimIntervalToContig( final String contig, final int start, final int stop, final int contigLength, final SimpleInterval expectedInterval ) {
        Assert.assertEquals(IntervalUtils.trimIntervalToContig(contig, start, stop, contigLength), expectedInterval);
    }

    @DataProvider(name = "TrimIntervalToContigInvalidData")
    public Object[][] trimIntervalToContigInvalidData() {
        return new Object[][] {
                { null, 1, 10, 100 },
                { "1", 1, 10, 0 },
                { "1", 1, 10, -1 }
        };
    }

    @Test(dataProvider = "TrimIntervalToContigInvalidData", expectedExceptions = IllegalArgumentException.class)
    public void testTrimIntervalToContigInvalidArg( final String contig, final int start, final int stop, final int contigLength ) {
        IntervalUtils.trimIntervalToContig(contig, start, stop, contigLength);
    }

    @DataProvider(name = "IntervalIsOnDictionaryContigData")
    public Object[][] intervalIsOnDictionaryContigData() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 1000)));

        return new Object[][] {
                { new SimpleInterval("1", 1, 10), dictionary, true },
                { new SimpleInterval("1", 50, 100), dictionary, true },
                { new SimpleInterval("1", 1, 1000), dictionary, true },
                { new SimpleInterval("2", 1, 10), dictionary, false },
                { new SimpleInterval("1", 1, 1001), dictionary, false },
                { new SimpleInterval("1", 1, 2000), dictionary, false }
        };
    }

    @Test(dataProvider = "IntervalIsOnDictionaryContigData")
    public void testIntervalIsOnDictionaryContig( final SimpleInterval interval, final SAMSequenceDictionary dictionary, final boolean expectedResult ) {
        Assert.assertEquals(IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary), expectedResult);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalIsOnDictionaryContigNullInterval() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 1000)));
        IntervalUtils.intervalIsOnDictionaryContig(null, dictionary);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalIsOnDictionaryContigNullDictionary() {
        IntervalUtils.intervalIsOnDictionaryContig(new SimpleInterval("1", 1, 10), null);
    }

    @Test
    public void testToGatkIntervalString(){
            final SimpleInterval interval = new SimpleInterval("1",1,100);
            Assert.assertEquals(IntervalUtils.locatableToString(interval), "1:1-100");
    }

    @DataProvider(name = "combineIntervalsTesting")
    public Object [][] combineIntervalsTesting() {

        // Overlap where one fits in another and an interval that is totally disconnected
        final List<Locatable> input1_1 = new ArrayList<>();
        final List<Locatable> input1_2 = new ArrayList<>();
        final List<Locatable> output1 = new ArrayList<>();
        input1_1.add(new SimpleInterval("1", 1000, 2000));
        input1_2.add(new SimpleInterval("1", 500, 2500));
        input1_2.add(new SimpleInterval("1", 2501, 3000));
        input1_2.add(new SimpleInterval("1", 4000, 5000));
        output1.add(new SimpleInterval("1", 500, 999));
        output1.add(new SimpleInterval("1", 1000, 2000));
        output1.add(new SimpleInterval("1", 2001, 2500));
        output1.add(new SimpleInterval("1", 2501, 3000));
        output1.add(new SimpleInterval("1", 4000, 5000));

        // One is empty
        final List<Locatable> input2_1 = new ArrayList<>();
        final List<Locatable> input2_2 = new ArrayList<>();
        final List<Locatable> output2 = new ArrayList<>();
        input2_1.add(new SimpleInterval("1", 1000, 2000));
        output2.add(input2_1.get(0));

        // Other is empty
        final List<Locatable> input3_1 = new ArrayList<>();
        final List<Locatable> input3_2 = new ArrayList<>();
        final List<Locatable> output3 = new ArrayList<>();
        input3_2.add(new SimpleInterval("1", 1000, 2000));
        output3.add(input3_2.get(0));

        // Simple overlap
        final List<Locatable> input4_1 = new ArrayList<>();
        final List<Locatable> input4_2 = new ArrayList<>();
        final List<Locatable> output4 = new ArrayList<>();
        input4_1.add(new SimpleInterval("1", 500, 1000));
        input4_2.add(new SimpleInterval("1", 750, 1500));
        output4.add(new SimpleInterval("1", 500, 749));
        output4.add(new SimpleInterval("1", 750, 1000));
        output4.add(new SimpleInterval("1", 1001, 1500));

        // No overlap
        final List<Locatable> input5_1 = new ArrayList<>();
        final List<Locatable> input5_2 = new ArrayList<>();
        final List<Locatable> output5 = new ArrayList<>();
        input5_1.add(new SimpleInterval("1", 500, 1000));
        input5_2.add(new SimpleInterval("1", 2000, 3000));
        output5.add(input5_1.get(0));
        output5.add(input5_2.get(0));

        // One is null
        final List<Locatable> input6_1 = new ArrayList<>();
        final List<Locatable> input6_2 = null;
        final List<Locatable> output6 = new ArrayList<>();
        input6_1.add(new SimpleInterval("1", 1000, 2000));
        output6.add(input6_1.get(0));

        // Other is null
        final List<Locatable> input7_1 = null;
        final List<Locatable> input7_2 = new ArrayList<>();
        final List<Locatable> output7 = new ArrayList<>();
        input7_2.add(new SimpleInterval("1", 1000, 2000));
        output7.add(input7_2.get(0));

        // single bp overlap at start
        final List<Locatable> input8_1 = new ArrayList<>();
        final List<Locatable> input8_2 = new ArrayList<>();
        final List<Locatable> output8 = new ArrayList<>();
        input8_1.add(new SimpleInterval("1", 500, 500));
        input8_2.add(new SimpleInterval("1", 500, 3000));
        output8.add(new SimpleInterval("1", 500, 500));
        output8.add(new SimpleInterval("1", 501, 3000));

        // single bp overlap at end
        final List<Locatable> input9_1 = new ArrayList<>();
        final List<Locatable> input9_2 = new ArrayList<>();
        final List<Locatable> output9 = new ArrayList<>();
        input9_1.add(new SimpleInterval("1", 3000, 3000));
        input9_2.add(new SimpleInterval("1", 500, 3000));
        output9.add(new SimpleInterval("1", 500, 2999));
        output9.add(new SimpleInterval("1", 3000, 3000));

        // single bp overlap in middle
        final List<Locatable> input10_1 = new ArrayList<>();
        final List<Locatable> input10_2 = new ArrayList<>();
        final List<Locatable> output10 = new ArrayList<>();
        input10_1.add(new SimpleInterval("1", 2000, 2000));
        input10_2.add(new SimpleInterval("1", 500, 3000));
        output10.add(new SimpleInterval("1", 500, 1999));
        output10.add(new SimpleInterval("1", 2000, 2000));
        output10.add(new SimpleInterval("1", 2001, 3000));

        // One interval overlaps two entirely
        final List<Locatable> input11_1 = new ArrayList<>();
        final List<Locatable> input11_2 = new ArrayList<>();
        final List<Locatable> output11 = new ArrayList<>();
        input11_1.add(new SimpleInterval("1", 2000, 2000));
        input11_1.add(new SimpleInterval("1", 2100, 2200));
        input11_2.add(new SimpleInterval("1", 500, 3000));
        output11.add(new SimpleInterval("1", 500, 1999));
        output11.add(new SimpleInterval("1", 2000, 2000));
        output11.add(new SimpleInterval("1", 2001, 2099));
        output11.add(new SimpleInterval("1", 2100, 2200));
        output11.add(new SimpleInterval("1", 2201, 3000));

        // One interval overlaps two partly
        final List<Locatable> input12_1 = new ArrayList<>();
        final List<Locatable> input12_2 = new ArrayList<>();
        final List<Locatable> output12 = new ArrayList<>();
        input12_1.add(new SimpleInterval("1", 1000, 2000));
        input12_1.add(new SimpleInterval("1", 3000, 4000));
        input12_2.add(new SimpleInterval("1", 1500, 3500));
        output12.add(new SimpleInterval("1", 1000, 1499));
        output12.add(new SimpleInterval("1", 1500, 2000));
        output12.add(new SimpleInterval("1", 2001, 2999));
        output12.add(new SimpleInterval("1", 3000, 3500));
        output12.add(new SimpleInterval("1", 3501, 4000));

        // Two nulls going in should yield an empty list.
        final List<Locatable> input13_1 = null;
        final List<Locatable> input13_2 = null;
        final List<Locatable> output13 = new ArrayList<>();

        // One interval overlaps two partly then an empty space then a single bp interval
        final List<Locatable> input14_1 = new ArrayList<>();
        final List<Locatable> input14_2 = new ArrayList<>();
        final List<Locatable> output14 = new ArrayList<>();
        input14_1.add(new SimpleInterval("1", 1000, 2000));
        input14_1.add(new SimpleInterval("1", 3000, 4000));
        input14_1.add(new SimpleInterval("1", 5000, 5000));
        input14_2.add(new SimpleInterval("1", 1500, 3500));
        output14.add(new SimpleInterval("1", 1000, 1499));
        output14.add(new SimpleInterval("1", 1500, 2000));
        output14.add(new SimpleInterval("1", 2001, 2999));
        output14.add(new SimpleInterval("1", 3000, 3500));
        output14.add(new SimpleInterval("1", 3501, 4000));
        output14.add(new SimpleInterval("1", 5000, 5000));

        // One interval overlaps two partly, but not really, since different contig.
        final List<Locatable> input15_1 = new ArrayList<>();
        final List<Locatable> input15_2 = new ArrayList<>();
        final List<Locatable> output15 = new ArrayList<>();
        input15_1.add(new SimpleInterval("1", 1000, 2000));
        input15_1.add(new SimpleInterval("1", 3000, 4000));
        input15_1.add(new SimpleInterval("2", 5000, 5000));
        input15_2.add(new SimpleInterval("2", 1500, 3500));
        output15.addAll(input15_1.subList(0,2));
        output15.addAll(input15_2);
        output15.add(input15_1.get(2));

        // Multiple contigs
        final List<Locatable> input16_1 = new ArrayList<>();
        final List<Locatable> input16_2 = new ArrayList<>();
        final List<Locatable> output16 = new ArrayList<>();
        input16_1.add(new SimpleInterval("1", 1000, 2000));
        input16_1.add(new SimpleInterval("1", 3000, 4000));
        input16_1.add(new SimpleInterval("2", 5000, 5000));
        input16_1.add(new SimpleInterval("2", 15000, 16000));

        input16_2.add(new SimpleInterval("1", 3500, 4500));
        input16_2.add(new SimpleInterval("2", 1500, 3500));
        input16_2.add(new SimpleInterval("2", 15000, 35000));

        output16.addAll(input16_1.subList(0,1));
        output16.add(new SimpleInterval("1", 3000, 3499));
        output16.add(new SimpleInterval("1", 3500, 4000));
        output16.add(new SimpleInterval("1", 4001, 4500));

        output16.add(input16_2.get(1));
        output16.add(input16_1.get(2));
        output16.add(new SimpleInterval("2", 15000, 16000));
        output16.add(new SimpleInterval("2", 16001, 35000));


        // Trivial to combine but has a lot of contigs
        final List<Locatable> input17_1 = Lists.newArrayList(
                new SimpleInterval("1", 1000, 2000),
                new SimpleInterval("1", 5000, 6000),
                new SimpleInterval("1", 7000, 8000),
                new SimpleInterval("1", 10000, 11000),
                new SimpleInterval("1", 12000, 14000),
                new SimpleInterval("10", 1000, 2000),
                new SimpleInterval("10", 5000, 6000),
                new SimpleInterval("10", 7000, 8000),
                new SimpleInterval("10", 10000, 11000),
                new SimpleInterval("10", 12000, 14000),
                new SimpleInterval("11", 1000, 2000),
                new SimpleInterval("11", 5000, 6000),
                new SimpleInterval("11", 7000, 8000),
                new SimpleInterval("11", 10000, 11000),
                new SimpleInterval("11", 12000, 14000),
                new SimpleInterval("3", 1000, 2000),
                new SimpleInterval("3", 5000, 6000),
                new SimpleInterval("3", 7000, 8000),
                new SimpleInterval("3", 10000, 11000),
                new SimpleInterval("3", 12000, 14000),
                new SimpleInterval("2", 1000, 2000),
                new SimpleInterval("2", 5000, 6000),
                new SimpleInterval("2", 7000, 8000),
                new SimpleInterval("2", 10000, 11000),
                new SimpleInterval("2", 12000, 14000));
        final List<Locatable> input17_2 = Lists.newArrayList(
                new SimpleInterval("10", 1000, 2000),
                new SimpleInterval("10", 5000, 6000),
                new SimpleInterval("10", 7000, 8000),
                new SimpleInterval("10", 10000, 10500),
                new SimpleInterval("10", 12000, 14000),
                new SimpleInterval("2", 1000, 2000),
                new SimpleInterval("2", 5000, 6000),
                new SimpleInterval("2", 7000, 8000),
                new SimpleInterval("2", 10000, 11000),
                new SimpleInterval("2", 12000, 14000),
                new SimpleInterval("3", 1000, 2000),
                new SimpleInterval("3", 5000, 6000),
                new SimpleInterval("3", 7000, 8000),
                new SimpleInterval("3", 10000, 11000),
                new SimpleInterval("3", 12000, 14000));
        final List<Locatable> output17 = Lists.newArrayList(
                new SimpleInterval("1", 1000, 2000),
                new SimpleInterval("1", 5000, 6000),
                new SimpleInterval("1", 7000, 8000),
                new SimpleInterval("1", 10000, 11000),
                new SimpleInterval("1", 12000, 14000),
                new SimpleInterval("2", 1000, 2000),
                new SimpleInterval("2", 5000, 6000),
                new SimpleInterval("2", 7000, 8000),
                new SimpleInterval("2", 10000, 11000),
                new SimpleInterval("2", 12000, 14000),
                new SimpleInterval("3", 1000, 2000),
                new SimpleInterval("3", 5000, 6000),
                new SimpleInterval("3", 7000, 8000),
                new SimpleInterval("3", 10000, 11000),
                new SimpleInterval("3", 12000, 14000),
                new SimpleInterval("10", 1000, 2000),
                new SimpleInterval("10", 5000, 6000),
                new SimpleInterval("10", 7000, 8000),
                new SimpleInterval("10", 10000, 10500),
                new SimpleInterval("10", 10501, 11000),
                new SimpleInterval("10", 12000, 14000),
                new SimpleInterval("11", 1000, 2000),
                new SimpleInterval("11", 5000, 6000),
                new SimpleInterval("11", 7000, 8000),
                new SimpleInterval("11", 10000, 11000),
                new SimpleInterval("11", 12000, 14000)
                );

        return new Object[][]{
                {input1_1, input1_2, output1},
                {input2_1, input2_2, output2},
                {input3_1, input3_2, output3},
                {input4_1, input4_2, output4},
                {input5_1, input5_2, output5},
                {input6_1, input6_2, output6},
                {input7_1, input7_2, output7},
                {input8_1, input8_2, output8},
                {input9_1, input9_2, output9},
                {input10_1, input10_2, output10},
                {input11_1, input11_2, output11},
                {input12_1, input12_2, output12},
                {input13_1, input13_2, output13},
                {input14_1, input14_2, output14},
                {input15_1, input15_2, output15},
                {input16_1, input16_2, output16},
                {input17_1, input17_2, output17},
        };
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCombineIntervalsErrorOverlapping() {
        final List<Locatable> input1_1 = new ArrayList<>();
        final List<Locatable> input1_2 = new ArrayList<>();
        input1_1.add(new SimpleInterval("1", 1000, 2000));
        input1_1.add(new SimpleInterval("1", 2000, 3000));
        input1_2.add(new SimpleInterval("1", 500, 2500));
        input1_2.add(new SimpleInterval("1", 2501, 3000));
        input1_2.add(new SimpleInterval("1", 4000, 5000));
        IntervalUtils.combineAndSortBreakpoints(input1_1, input1_2, getSamSequenceDictionaryForCombineIntervalTests());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCombineIntervalsErrorOverlapping2() {
        final List<Locatable> input1_1 = new ArrayList<>();
        final List<Locatable> input1_2 = new ArrayList<>();
        input1_1.add(new SimpleInterval("1", 1000, 2000));
        input1_2.add(new SimpleInterval("1", 500, 2500));
        input1_2.add(new SimpleInterval("1", 2501, 3000));
        input1_2.add(new SimpleInterval("1", 4000, 5000));
        input1_2.add(new SimpleInterval("1", 4000, 5000));
        IntervalUtils.combineAndSortBreakpoints(input1_1, input1_2, getSamSequenceDictionaryForCombineIntervalTests());
    }

    /**
     * Tests the sorting as well when using combining intervals
     */
    @Test(dataProvider = "combineIntervalsTesting")
    public void testCombineIntervalsWithShuffling(List<Locatable> locatables1, List<Locatable> locatables2, List<Locatable> gtOutput) {

        if (locatables1 != null) {
            Collections.shuffle(locatables1, new Random(4040));
        }
        if (locatables2 != null) {
            Collections.shuffle(locatables2, new Random(4040));
        }

        final SAMSequenceDictionary dictionary = getSamSequenceDictionaryForCombineIntervalTests();

        final List<Locatable> outputs = IntervalUtils.combineAndSortBreakpoints(locatables1, locatables2, dictionary);
        Assert.assertEquals(outputs.size(), gtOutput.size());
        Assert.assertEquals(outputs, gtOutput);
    }

    private SAMSequenceDictionary getSamSequenceDictionaryForCombineIntervalTests() {
        return new SAMSequenceDictionary(
                Arrays.asList(new SAMSequenceRecord("1", 500000),
                        new SAMSequenceRecord("2", 500000),
                        new SAMSequenceRecord("3", 500000),
                        new SAMSequenceRecord("4", 500000),
                        new SAMSequenceRecord("10", 500000),
                        new SAMSequenceRecord("11", 500000)));
    }

    @DataProvider(name = "overlappingData")
    public Object [][] createOverlappingData() {
        final List<SimpleInterval> k = Arrays.asList(new SimpleInterval("1", 100, 500));
        final List<SimpleInterval> k3 = Arrays.asList(new SimpleInterval("1", 100, 500),
                new SimpleInterval("2", 100, 500));
        final List<SimpleInterval> k4 = Arrays.asList(new SimpleInterval("1", 100, 500),
                new SimpleInterval("1", 800, 800),
                new SimpleInterval("2", 300, 500));
        final List<SimpleInterval> k5 = Arrays.asList(new SimpleInterval("1", 100, 300),
                new SimpleInterval("1", 301, 400));

        final List<SimpleInterval> k6 = Arrays.asList(
                new SimpleInterval("1",10001,10500),
                new SimpleInterval("1", 52500, 109750),
                new SimpleInterval("1", 109751, 230500),
                new SimpleInterval("1", 230501, 258500));

        final List<SimpleInterval> vs = Arrays.asList(
                new SimpleInterval("1", 100, 200),
                new SimpleInterval("1", 300, 800)
        );
        final List<SimpleInterval> vs2 = Arrays.asList(
                new SimpleInterval("1", 100, 200),
                new SimpleInterval("1", 700, 800)
        );
        final List<SimpleInterval> vs3 = Arrays.asList(
                new SimpleInterval("1", 800, 1200),
                new SimpleInterval("2", 100, 800)
        );
        final List<SimpleInterval> vs5 = Arrays.asList(
                new SimpleInterval("1", 100, 200),
                new SimpleInterval("1", 201, 500)
        );
        final List<SimpleInterval> vs6 = Arrays.asList(
                new SimpleInterval("1",5000, 60000),
                new SimpleInterval("1",70000, 100000),
                new SimpleInterval("1",120000, 220000),
                new SimpleInterval("1",230000, 300000)
        );

        final List<SimpleInterval> vs7 = Arrays.asList(
                new SimpleInterval("1",50, 200)
        );

        return new Object[][]{
                // The new ArrayList calls are to preserve order since we may shuffle inputs.
                // Simple tests
                {k, vs, ImmutableMap.of(k.get(0), new ArrayList<>(vs))},
                {k, vs2, ImmutableMap.of(k.get(0), new ArrayList<>(vs2).subList(0,1))},
                {k, vs3, ImmutableMap.of(k.get(0), Collections.emptyList())},

                // Slightly more sophisticated tests.
                {k3, vs, ImmutableMap.of(k3.get(0), new ArrayList<>(vs), k3.get(1), Collections.emptyList())},
                {k4, vs, ImmutableMap.of(k4.get(0), new ArrayList<>(vs), k4.get(1), new ArrayList<>(vs).subList(1,2), k4.get(2), Collections.emptyList())},
                {k5, vs5, ImmutableMap.of(k5.get(0), new ArrayList<>(vs5), k5.get(1), new ArrayList<>(vs5).subList(1,2))},
                {k6, vs6, ImmutableMap.of(k6.get(0), new ArrayList<>(vs6).subList(0,1), k6.get(1), new ArrayList<>(vs6).subList(0,2),
                        k6.get(2), new ArrayList<>(vs6).subList(2,4), k6.get(3), new ArrayList<>(vs6).subList(3,4))},
                {k, vs7, ImmutableMap.of(k.get(0), new ArrayList<>(vs7))},
        };
    }

    @Test(dataProvider = "overlappingData")
    public void testCreateOverlapMap(List<Locatable> locatables1, List<Locatable> locatables2, Map<Locatable, List<Locatable>> gtOutput) {
        assertCreateOverlapMap(locatables1, locatables2, gtOutput);
    }

    @Test(dataProvider = "overlappingData")
    public void testCreateOverlapMapUnsorted(List<Locatable> locatables1, List<Locatable> locatables2, Map<Locatable, List<Locatable>> gtOutput) {

        if (locatables1 != null) {
            Collections.shuffle(locatables1, new Random(1234));
        }
        if (locatables2 != null) {
            Collections.shuffle(locatables2, new Random(4311));
        }

        assertCreateOverlapMap(locatables1, locatables2, gtOutput);
    }

    private void assertCreateOverlapMap(List<Locatable> locatables1, List<Locatable> locatables2, Map<Locatable, List<Locatable>> gtOutput) {
        final Map<Locatable, List<Locatable>> outputs = IntervalUtils.createOverlapMap(locatables1, locatables2,
                getSamSequenceDictionaryForCombineIntervalTests());
        Assert.assertEquals(outputs.size(), gtOutput.size());
        Assert.assertEquals(outputs.keySet().size(), new HashSet<>(locatables1).size());
        Assert.assertTrue(Sets.difference(outputs.keySet(), gtOutput.keySet()).size() == 0);

        for (final Map.Entry<Locatable, List<Locatable>> overlap : outputs.entrySet()) {
            Assert.assertEquals(overlap.getValue(), gtOutput.get(overlap.getKey()));
        }
    }

    @Test(dataProvider = "genomicSortingTests")
    public void testBasicGenomicSortOfLocatable(final List<Locatable> testList, final List<Locatable> gtList) throws IOException {

        final SAMSequenceDictionary dictionary = SAMSequenceDictionaryExtractor.extractDictionary(new File(FULL_HG19_DICT).toPath());
        testList.sort(IntervalUtils.getDictionaryOrderComparator(dictionary));
        Assert.assertEquals(testList, gtList);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBasicGenomicSortOfLocatableErrorContigNotInDictionary() throws IOException {
        final List<Locatable> testList = Arrays.asList(
               new SimpleInterval("2", 1500, 2000),
               new SimpleInterval("11", 1000, 2000),
               new SimpleInterval("1", 1500, 2000),
               new SimpleInterval("10", 1500, 2000),
               new SimpleInterval("10", 2500, 3000),
               new SimpleInterval("X", 2500, 3000),
               new SimpleInterval("X", 3500, 4000),
               new SimpleInterval("MT", 3500, 4000),
               new SimpleInterval("AL1234.123NOT_IN_DICT", 3500, 4000),
               new SimpleInterval("GL000207.1", 3500, 4000),
               new SimpleInterval("Y", 3500, 4000),
               new SimpleInterval("GL1234.123NOT_IN_DICT", 3700, 4000),
               new SimpleInterval("GL1234.123NOT_IN_DICT", 3500, 4000)
        );

        final SAMSequenceDictionary dictionary = SAMSequenceDictionaryExtractor.extractDictionary(new File(FULL_HG19_DICT).toPath());
        testList.sort(IntervalUtils.getDictionaryOrderComparator(dictionary));
    }

    @DataProvider(name = "genomicSortingTests")
    public Object[][] createGenomicSortTests() {
        final List<Locatable> list1 = Arrays.asList(
                new SimpleInterval("2", 1500, 2000),
                new SimpleInterval("1", 1000, 2000),
                new SimpleInterval("1", 1500, 2000),
                new SimpleInterval("10", 1500, 2000)
        );
        final List<Locatable> gtList1 = Arrays.asList(
                list1.get(1), list1.get(2), list1.get(0), list1.get(3));

        final List<Locatable> list2 = Arrays.asList(
                new SimpleInterval("2", 1500, 2000),
                new SimpleInterval("11", 1000, 2000),
                new SimpleInterval("1", 1500, 2000),
                new SimpleInterval("10", 1500, 2000),
                new SimpleInterval("10", 2500, 3000),
                new SimpleInterval("X", 2500, 3000),
                new SimpleInterval("X", 3500, 4000),
                new SimpleInterval("MT", 3500, 4000),
                new SimpleInterval("GL000207.1", 3500, 4000),
                new SimpleInterval("Y", 3500, 4000)
        );
        final List<Locatable> gtList2 = Arrays.asList(
                list2.get(2), list2.get(0), list2.get(3), list2.get(4), list2.get(1), list2.get(5),
                list2.get(6), list2.get(9), list2.get(7), list2.get(8));

        return new Object[][]{
                {list1, gtList1},
                {list2, gtList2}
        };
    }

    @DataProvider
    public Object[][] createGetIntervalsWithFlanksData() {
        final SAMSequenceDictionary dictionary = ArtificialReadUtils.createArtificialSamHeader().getSequenceDictionary();

        return new Object[][] {
                // input intervals, padding, dictionary, expected output intervals
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400), new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250)),
                  0, dictionary,
                  Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400), new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250))
                },
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400), new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250)),
                  20, dictionary,
                  Arrays.asList(new SimpleInterval("1", 1, 120), new SimpleInterval("1", 130, 220), new SimpleInterval("1", 280, 420), new SimpleInterval("2", 30, 270))
                },
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400), new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250)),
                  25, dictionary,
                  Arrays.asList(new SimpleInterval("1", 1, 225), new SimpleInterval("1", 275, 425), new SimpleInterval("2", 25, 275))
                },
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400), new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250)),
                  50, dictionary,
                  Arrays.asList(new SimpleInterval("1", 1, 450), new SimpleInterval("2", 1, 300))
                }
        };
    }

    @Test(dataProvider = "createGetIntervalsWithFlanksData")
    public void testGetIntervalsWithFlanks(final List<SimpleInterval> inputIntervals, final int padding, final SAMSequenceDictionary dictionary, final List<SimpleInterval> expectedIntervals) {
        final List<SimpleInterval> result = IntervalUtils.getIntervalsWithFlanks(inputIntervals, padding, dictionary);
        Assert.assertEquals(result, expectedIntervals);
    }

    @DataProvider
    public Object[][] createGroupIntervalsByContigTestData() {
        return new Object[][] {
                // input intervals, expected output
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400), new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250)),
                  Arrays.asList(Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200), new SimpleInterval("1", 300, 400)), Arrays.asList(new SimpleInterval("2", 50, 150), new SimpleInterval("2", 175, 250)))
                },
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("2", 150, 200), new SimpleInterval("3", 300, 400), new SimpleInterval("5", 50, 150), new SimpleInterval("6", 175, 250)),
                  Arrays.asList(Arrays.asList(new SimpleInterval("1", 1, 100)), Arrays.asList(new SimpleInterval("2", 150, 200)), Arrays.asList(new SimpleInterval("3", 300, 400)), Arrays.asList(new SimpleInterval("5", 50, 150)), Arrays.asList(new SimpleInterval("6", 175, 250)))
                },
                { Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200)),
                  Arrays.asList(Arrays.asList(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 150, 200)))
                },
                { Arrays.asList(new SimpleInterval("1", 1, 100)),
                  Arrays.asList(Arrays.asList(new SimpleInterval("1", 1, 100)))
                },
        };
    }

    @Test(dataProvider = "createGroupIntervalsByContigTestData")
    public void testGroupIntervalsByContig(final List<SimpleInterval> inputIntervals, final List<List<SimpleInterval>> expectedResult) {
        final List<List<SimpleInterval>> result = IntervalUtils.groupIntervalsByContig(inputIntervals);
        Assert.assertEquals(result, expectedResult);
    }

    @DataProvider
    public Object[][] provideReciprocalOverlapTestData(){
        return new Object[][] {
                // Simple, obvious tests
                {new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 50), 0.5, true},
                {new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 50), 0.51, false},
                {new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 50), 0.49, true},
                {new SimpleInterval("1", 1, 50), new SimpleInterval("1", 51, 100), 0., true},
                {new SimpleInterval("1", 1, 50), new SimpleInterval("1", 51, 100), 1.0, false},
                {new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100), 1.0, true},
                {new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 101), 1.0, false},

                // Contig difference -- always false, unless threshold is zero
                {new SimpleInterval("2", 1, 100), new SimpleInterval("1", 1, 50), 0.5, false},
                {new SimpleInterval("1", 1, 100), new SimpleInterval("2", 1, 50), 0.51, false},
                {new SimpleInterval("1", 1, 100), new SimpleInterval("2", 1, 50), 0.49, false},
                {new SimpleInterval("2", 1, 100), new SimpleInterval("1", 1, 50), 0., true},

                // Total overlap
                {new SimpleInterval("1", 2, 101), new SimpleInterval("1", 1, 500), 0.5, false},
                {new SimpleInterval("1", 2, 101), new SimpleInterval("1", 1, 500), 0.1, true},
                {new SimpleInterval("1", 2, 101), new SimpleInterval("1", 1, 500), 0.199999, true},
                {new SimpleInterval("1", 2, 101), new SimpleInterval("1", 1, 500), 0.2, true},
                {new SimpleInterval("1", 2, 101), new SimpleInterval("1", 1, 500), 0.21, false},

                // Partial overlap -- low side
                {new SimpleInterval("1", 51, 150), new SimpleInterval("1", 101, 600), 0.5, false},
                {new SimpleInterval("1", 51, 150), new SimpleInterval("1", 101, 600), 0.1, true},
                {new SimpleInterval("1", 51, 150), new SimpleInterval("1", 101, 600), 0.09, true},
                {new SimpleInterval("1", 51, 150), new SimpleInterval("1", 101, 600), 0.11, false},

                // Partial overlap -- high side
                {new SimpleInterval("1", 551, 650), new SimpleInterval("1", 101, 600), 0.5, false},
                {new SimpleInterval("1", 551, 650), new SimpleInterval("1", 101, 600), 0.1, true},
                {new SimpleInterval("1", 551, 650), new SimpleInterval("1", 101, 600), 0.09, true},
                {new SimpleInterval("1", 551, 650), new SimpleInterval("1", 101, 600), 0.11, false},
        };
    }

    /**
     * Also tests that the reciprocal overlap is commutative.
     */
    @Test(dataProvider = "provideReciprocalOverlapTestData")
    public void testIsReciprocalOverlap(final SimpleInterval interval1, final SimpleInterval interval2, final double reciprocalThreshold, final boolean gt) {
        Assert.assertEquals(IntervalUtils.isReciprocalOverlap(interval1, interval2, reciprocalThreshold), gt);
        Assert.assertEquals(IntervalUtils.isReciprocalOverlap(interval2, interval1, reciprocalThreshold), gt);
    }
}
