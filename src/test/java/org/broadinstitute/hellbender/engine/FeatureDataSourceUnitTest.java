package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public final class FeatureDataSourceUnitTest extends GATKBaseTest {
    private static final String FEATURE_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File QUERY_TEST_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");
    private static final File QUERY_TEST_GVCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test_gvcf.vcf");
    private static final File UNINDEXED_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "unindexed.vcf");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFile() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>((File)null);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testHandleNonExistentFile() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(
                GATKBaseTest.getSafeNonExistentFile("nonexistent.vcf"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleInvalidQueryLookahead() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, "MyName", -1);
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleQueryOverUnindexedFile() {
        try ( FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(UNINDEXED_VCF) ) {
            featureSource.query(new SimpleInterval("1", 1, 1));  // Should throw, since we have no index
        }
    }

    @Test
    public void testGetName() {
        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, "CustomName")) {
            Assert.assertEquals(featureSource.getName(), "CustomName", "Wrong name returned from getHeader()");
        }
    }

    @Test
    public void testGetHeader() {
        Object header = null;
        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, "CustomName")) {
            header = featureSource.getHeader();
        }
        Assert.assertTrue(header instanceof VCFHeader, "Header for " + QUERY_TEST_VCF.getAbsolutePath() + " not a VCFHeader");
    }

    @Test
    public void testGetSequenceDictionaryFromIndex() {
        // this vcf has no sequence dictionary, but it has an index....
        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, "CustomName")) {
            final SAMSequenceDictionary dict = featureSource.getSequenceDictionary();
            Assert.assertEquals(dict.size(), 4);
            Assert.assertEquals(dict.getSequences().stream().map(s->s.getSequenceName()).collect(Collectors.toList()), Arrays.asList("1", "2", "3", "4"));
        }
    }

    @DataProvider(name = "CompleteIterationTestData")
    public Object[][] getCompleteIterationTestData() {
        // File to iterate over + Expected Variant ID(s)
        return new Object[][] {
                { QUERY_TEST_VCF, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z") },
                { UNINDEXED_VCF, Arrays.asList("a", "b", "c") }
        };
    }

    @Test(dataProvider = "CompleteIterationTestData")
    public void testCompleteIterationOverFile( final File vcfFile, final List<String> expectedVariantIDs ) {
        try ( FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(vcfFile) ) {
            Iterator<VariantContext> iter = featureSource.iterator();

            checkTraversalResults(iter, expectedVariantIDs, vcfFile, null);
        }
    }

    @DataProvider(name = "TraversalByIntervalsTestData")
    public Object[][] getTraversalByIntervalsTestData() {
        // Intervals for traversal + expected Variant IDs
        return new Object[][] {
                // Single interval
                { Arrays.asList(new SimpleInterval("1", 100, 200)), Arrays.asList("a", "b", "c") },

                // Two non-adjacent intervals on the same contig
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Some records overlap multiple intervals, and there are gaps between intervals
                { Arrays.asList(new SimpleInterval("1", 100, 203), new SimpleInterval("1", 205, 284), new SimpleInterval("1", 286, 1000)), Arrays.asList("a", "b", "c", "d", "e", "f", "h", "i", "j", "k") },

                // Some records overlap multiple intervals, and no gaps between intervals
                { Arrays.asList(new SimpleInterval("1", 100, 203), new SimpleInterval("1", 204, 285), new SimpleInterval("1", 286, 1000)), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k") },

                // Two intervals on different contigs
                { Arrays.asList(new SimpleInterval("1", 200, 300), new SimpleInterval("2", 500, 600)), Arrays.asList("b", "c", "d", "e", "f", "g", "h", "p", "q") },

                // More than two intervals spanning different contigs, and some records overlap multiple intervals
                { Arrays.asList(new SimpleInterval("1", 200, 203), new SimpleInterval("1", 205, 285), new SimpleInterval("2", 200, 548), new SimpleInterval("2", 550, 650), new SimpleInterval("4", 700, 800)), Arrays.asList("b", "c", "d", "e", "f", "g", "o", "p", "q", "r", "y", "z") },

                // One interval with no overlapping records at the beginning of interval list
                { Arrays.asList(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Multiple intervals with no overlapping records at the beginning of interval list
                { Arrays.asList(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 60, 70), new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // One interval with no overlapping records in the middle of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 500, 600), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Multiple intervals with no overlapping records in the middle of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 500, 600), new SimpleInterval("1", 700, 800), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // One interval with no overlapping records at the end of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000), new SimpleInterval("1", 2000, 3000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Multiple intervals with no overlapping records at the end of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000), new SimpleInterval("1", 2000, 3000), new SimpleInterval("1", 4000, 5000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // No records overlap any intervals
                { Arrays.asList(new SimpleInterval("1", 1, 99), new SimpleInterval("1", 287, 290), new SimpleInterval("1", 500, 600), new SimpleInterval("2", 201, 524), new SimpleInterval("2", 1000, 2000), new SimpleInterval("4", 1, 500)), Collections.<String>emptyList() },

                // No intervals (should traverse the entire file)
                { Collections.<SimpleInterval>emptyList(), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z") }
        };
    }

    @Test(dataProvider = "TraversalByIntervalsTestData")
    public void testTraversalByIntervals( final List<SimpleInterval> intervalsForTraversal, final List<String> expectedVariantIDs ) {
        try ( FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF) ) {
            featureSource.setIntervalsForTraversal(intervalsForTraversal);
            Iterator<VariantContext> iter = featureSource.iterator();

            checkTraversalResults(iter, expectedVariantIDs, QUERY_TEST_VCF, intervalsForTraversal);
        }
    }

    private void checkTraversalResults( final Iterator<VariantContext> traversalResults, final List<String> expectedVariantIDs, final File vcfFile, final List<SimpleInterval> traversalIntervals ) {
        final String intervalString = traversalIntervals != null ? " with intervals " + traversalIntervals : "";

        int recordCount = 0;
        while ( traversalResults.hasNext() ) {
            VariantContext record = traversalResults.next();
            Assert.assertTrue(recordCount < expectedVariantIDs.size(), "Too many records returned during iteration over " + vcfFile.getAbsolutePath() + intervalString);
            Assert.assertEquals(record.getID(), expectedVariantIDs.get(recordCount),
                                "Record #" + (recordCount + 1) + " encountered in iteration over " + vcfFile.getAbsolutePath() + intervalString + " is incorrect");
            ++recordCount;
        }

        Assert.assertEquals(recordCount, expectedVariantIDs.size(), "Wrong number of records returned in iteration over " + vcfFile.getAbsolutePath() + intervalString);
    }

    @DataProvider(name = "IndependentFeatureQueryTestData")
    public Object[][] getIndependentFeatureQueryTestData() {
        // Query Interval + Expected Variant ID(s)
        return new Object[][] {
                { new SimpleInterval("1", 1, 99), Collections.<String>emptyList() },
                { new SimpleInterval("1", 100, 100), Arrays.asList("a") },
                { new SimpleInterval("1", 100, 200), Arrays.asList("a", "b", "c") },
                { new SimpleInterval("1", 200, 202), Arrays.asList("b", "c") },
                { new SimpleInterval("1", 200, 203), Arrays.asList("b", "c", "d") },
                { new SimpleInterval("1", 201, 203), Arrays.asList("d") },
                { new SimpleInterval("1", 204, 204), Arrays.asList("d") },
                { new SimpleInterval("1", 206, 206), Arrays.asList("d") },
                { new SimpleInterval("1", 207, 207), Collections.<String>emptyList() },
                { new SimpleInterval("1", 200, 300), Arrays.asList("b", "c", "d", "e", "f", "g", "h") },
                { new SimpleInterval("1", 275, 300), Arrays.asList("e", "f", "g", "h") },
                { new SimpleInterval("1", 275, 284), Arrays.asList("e", "f") },
                { new SimpleInterval("1", 284, 284), Arrays.asList("f") },
                { new SimpleInterval("1", 284, 285), Arrays.asList("f", "g") },
                { new SimpleInterval("1", 284, 286), Arrays.asList("f", "g", "h") },
                { new SimpleInterval("1", 286, 286), Arrays.asList("f", "h") },
                { new SimpleInterval("1", 287, 290), Collections.<String>emptyList() },
                { new SimpleInterval("1", 999, 1000), Arrays.asList("i", "j", "k") },
                { new SimpleInterval("1", 1000, 1001), Arrays.asList("j", "k") },
                { new SimpleInterval("1", 1002, 1005), Arrays.asList("k") },
                { new SimpleInterval("1", 1005, 1010), Collections.<String>emptyList() },
                { new SimpleInterval("1", 1075, 1175), Arrays.asList("l", "m") },
                { new SimpleInterval("1", 1075, 1176), Arrays.asList("l", "m", "n") },
                { new SimpleInterval("1", 1077, 1176), Arrays.asList("m", "n") },
                { new SimpleInterval("1", 1170, 1180), Arrays.asList("n") },
                { new SimpleInterval("1", 1003, 1175), Arrays.asList("k", "l", "m") },
                { new SimpleInterval("1", 1, 2000), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n") },
                { new SimpleInterval("1", 2000, 3000), Collections.<String>emptyList() },
                { new SimpleInterval("1", 286, 1175), Arrays.asList("f", "h", "i", "j", "k", "l", "m") },
                { new SimpleInterval("2", 200, 700), Arrays.asList("o", "p", "q", "r", "s") },
                { new SimpleInterval("2", 201, 699), Arrays.asList("p", "q", "r") },
                { new SimpleInterval("2", 550, 560), Arrays.asList("q") },
                { new SimpleInterval("2", 549, 699), Arrays.asList("q", "r") },
                { new SimpleInterval("2", 600, 700), Arrays.asList("r", "s") },
                { new SimpleInterval("2", 701, 800), Collections.<String>emptyList() },
                { new SimpleInterval("3", 1, 300), Arrays.asList("t", "u", "v") },
                { new SimpleInterval("3", 300, 400), Arrays.asList("u", "v", "w") },
                { new SimpleInterval("3", 301, 400), Arrays.asList("v", "w") },
                { new SimpleInterval("3", 301, 399), Arrays.asList("v") },
                { new SimpleInterval("3", 305, 400), Arrays.asList("w") },
                { new SimpleInterval("3", 500, 600), Collections.<String>emptyList() },
                { new SimpleInterval("4", 1, 1000), Arrays.asList("x", "y", "z") },
                { new SimpleInterval("4", 600, 775), Arrays.asList("x", "y") },
                { new SimpleInterval("4", 775, 776), Arrays.asList("y", "z") },
                { new SimpleInterval("4", 777, 780), Arrays.asList("z") },
                { new SimpleInterval("4", 777, 300000000), Arrays.asList("z") },
        };
    }

    @Test(expectedExceptions = ArithmeticException.class)
    public void testBlowUpOnOverflow() {
        final SimpleInterval queryInterval = new SimpleInterval("4", 777, Integer.MAX_VALUE);
        try (final FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF)) {
            Iterator<VariantContext> featureIterator = featureSource.query(queryInterval);
        }
    }

    /**
     * Tests correctness of the raw query operation (without any caching coming into play) by creating
     * a new FeatureDataSource for each query to eliminate caching as a factor.
     *
     * (We test the use case of re-using a single FeatureDataSource across queries in separate tests
     *  below to make sure the caching is working as expected)
     */
    @Test(dataProvider = "IndependentFeatureQueryTestData")
    public void testIndependentFeatureQuerying( final SimpleInterval queryInterval, final List<String> expectedVariantIDs ) {
        try (final FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF)) {

            // Use query() here rather than queryAndPrefetch() so that query() will have test coverage
            // (the other tests below use queryAndPrefetch())
            Iterator<VariantContext> featureIterator = featureSource.query(queryInterval);
            List<VariantContext> queryResults = new ArrayList<>();
            while (featureIterator.hasNext()) {
                queryResults.add(featureIterator.next());
            }

            checkVariantQueryResults(queryResults, expectedVariantIDs, queryInterval);
        }
    }

    @Test
    public void testQueryAllHG38Intervals() {
        SAMSequenceDictionary sd;
        final File testFile = new File (FEATURE_DATA_SOURCE_TEST_DIRECTORY, "Homo_sapiens_assembly38.headerOnly.vcf.gz");

        try (VCFFileReader vcfReader = new VCFFileReader(testFile, false)) {
            sd = vcfReader.getFileHeader().getSequenceDictionary();
        }

        // Test that we can execute a query using any hg38 contig name against a VCF with an hg38 sequence dictionary.
        // Since the query is provided as a SimpleInterval, no interval query parsing or disambiguation is executed by
        // this code path, but GenomeLocParserUnitTest IntervalUtilsUnitTest have corresponding tests that ensures that
        // no hg38 query can be ambiguous.
        try (final FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(testFile)) {
            sd.getSequences().stream().forEach(
                    hg38Contig -> featureSource.query(new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength())));
        }
    }

    private void checkVariantQueryResults( final List<VariantContext> queryResults, final List<String> expectedVariantIDs, final SimpleInterval queryInterval ) {
        Assert.assertEquals(queryResults.size(), expectedVariantIDs.size(), "Wrong number of records returned for query on interval " + queryInterval);
        for ( int recordIndex = 0; recordIndex < queryResults.size(); ++recordIndex ) {
            Assert.assertEquals(queryResults.get(recordIndex).getID(), expectedVariantIDs.get(recordIndex),
                                "Record #" + (recordIndex + 1) + " returned from query on interval " + queryInterval + " is incorrect");
        }
    }

    /**
     * This data provider serves up sets of multiple queries paired with a single FeatureDataSource.
     * The intention is to test caching behavior, including recovery from cache misses, and correctness
     * of the code that subsets the Feature cache down to Features overlapping the desired interval
     * (whereas the primary intent of the independent querying tests above was to test correctness of the
     * raw query operation without any caching involved).
     */
    @DataProvider(name = "SingleDataSourceMultipleQueriesTestData")
    @SuppressWarnings("unchecked")
    public Object[][] getSingleDataSourceMultipleQueriesTestData() {
        // Query set #1:
        // Re-use the queries + expected results from the IndependentFeatureQueryTestData DataProvider above,
        // but this time aggregated together and executed on the same FeatureDataSource
        List<Pair<SimpleInterval, List<String>>> aggregatedIndependentQueries = new ArrayList<>();
        Object[][] independentQueryTestData = getIndependentFeatureQueryTestData();
        for ( Object[] queryTest : independentQueryTestData ) {
            aggregatedIndependentQueries.add(Pair.of((SimpleInterval)queryTest[0], (List<String>)queryTest[1]));
        }

        // Query set #2:
        // Large intervals with regularly-increasing start positions. Represents typical query access patterns
        // ideal for the caching implementation in FeatureDataSource, as it minimizes cache misses.
        List<Pair<SimpleInterval, List<String>>> regularlyIncreasingQueries = Arrays.asList(
                Pair.of(new SimpleInterval("1", 1, 100), Arrays.asList("a")),
                Pair.of(new SimpleInterval("1", 50, 150), Arrays.asList("a")),
                Pair.of(new SimpleInterval("1", 100, 200), Arrays.asList("a", "b", "c")),
                Pair.of(new SimpleInterval("1", 150, 250), Arrays.asList("b", "c", "d")),
                Pair.of(new SimpleInterval("1", 200, 300), Arrays.asList("b", "c", "d", "e", "f", "g", "h")),
                Pair.of(new SimpleInterval("1", 250, 350), Arrays.asList("e", "f", "g", "h")),
                Pair.of(new SimpleInterval("1", 300, 400), Collections.<String>emptyList()),
                Pair.of(new SimpleInterval("1", 350, 450), Collections.<String>emptyList()),
                Pair.of(new SimpleInterval("1", 950, 1050), Arrays.asList("i", "j", "k")),
                Pair.of(new SimpleInterval("1", 1000, 1100), Arrays.asList("j", "k", "l")),
                // First cache miss here given the default queryLookaheadBases value of 1000
                Pair.of(new SimpleInterval("1", 1050, 1150), Arrays.asList("l", "m")),
                Pair.of(new SimpleInterval("1", 1100, 1200), Arrays.asList("m", "n")),
                Pair.of(new SimpleInterval("1", 1150, 1250), Arrays.asList("m", "n")),
                Pair.of(new SimpleInterval("1" ,1200, 1300), Collections.<String>emptyList()),
                // Second cache miss here as we change contigs
                Pair.of(new SimpleInterval("2", 1, 100), Collections.<String>emptyList()),
                Pair.of(new SimpleInterval("2", 100, 200), Arrays.asList("o")),
                Pair.of(new SimpleInterval("2", 500, 600), Arrays.asList("p", "q")),
                Pair.of(new SimpleInterval("2", 550, 650), Arrays.asList("q", "r")),
                Pair.of(new SimpleInterval("2", 600, 700), Arrays.asList("r", "s")),
                Pair.of(new SimpleInterval("2", 650, 750), Arrays.asList("s")),
                // Third cache miss (changing contigs again)
                Pair.of(new SimpleInterval("3", 1, 200), Arrays.asList("t")),
                Pair.of(new SimpleInterval("3", 300, 400), Arrays.asList("u", "v", "w")),
                Pair.of(new SimpleInterval("3", 302, 350), Arrays.asList("v"))
        );

        // Query set #3:
        // Cache miss hell: lots of cache misses in this query set due to either backing up or skipping ahead too far
        List<Pair<SimpleInterval, List<String>>> cacheMissQueries = Arrays.asList(
                Pair.of(new SimpleInterval("1", 100, 200), Arrays.asList("a", "b", "c")),
                // Cache miss due to backup
                Pair.of(new SimpleInterval("1", 99, 205), Arrays.asList("a", "b", "c", "d")),
                // Cache miss due to jumping past query lookahead
                Pair.of(new SimpleInterval("1", 2000, 3000), Collections.<String>emptyList()),
                // Cache miss due to backup
                Pair.of(new SimpleInterval("1", 205, 285), Arrays.asList("d", "e", "f", "g")),
                Pair.of(new SimpleInterval("1", 286, 400), Arrays.asList("f", "h")),
                // Cache miss due to contig change
                Pair.of(new SimpleInterval("2", 200, 600), Arrays.asList("o", "p", "q")),
                // Cache miss due to backup
                Pair.of(new SimpleInterval("2", 1, 100), Collections.<String>emptyList()),
                // Cache miss due to jumping past query lookahead
                Pair.of(new SimpleInterval("2", 3000, 4000), Collections.<String>emptyList()),
                // Cache miss due to contig change
                Pair.of(new SimpleInterval("1", 200, 300), Arrays.asList("b", "c", "d", "e", "f", "g", "h")),
                // Cache miss due to backup
                Pair.of(new SimpleInterval("1", 100, 200), Arrays.asList("a", "b", "c")),
                // Cache miss due to backup
                Pair.of(new SimpleInterval("1", 1, 1), Collections.<String>emptyList()),
                // Cache miss due to jumping past query lookahead
                Pair.of(new SimpleInterval("1", 1100, 1200), Arrays.asList("m", "n"))
        );

        return new Object[][] {
                { aggregatedIndependentQueries },
                { regularlyIncreasingQueries },
                { cacheMissQueries }
        };
    }

    /**
     * Tests correctness of Feature caching behavior by executing multiple queries on the same FeatureDataSource
     */
    @Test(dataProvider = "SingleDataSourceMultipleQueriesTestData")
    public void testSingleDataSourceMultipleQueries( final List<Pair<SimpleInterval, List<String>>> testQueries ) {
        try (final FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF)) {

            // This test re-uses the same FeatureDataSource across queries to test caching of query results.
            for ( Pair<SimpleInterval, List<String>> testQuery : testQueries ) {
                final SimpleInterval queryInterval = testQuery.getLeft();
                final List<String> expectedVariantIDs = testQuery.getRight();

                final List<VariantContext> queryResults = featureSource.queryAndPrefetch(queryInterval);
                checkVariantQueryResults(queryResults, expectedVariantIDs, queryInterval);
            }
        }
    }

    @DataProvider(name = "GVCFQueryTestData")
    public Object[][] getGVCFQueryTestData() {

        // Ensure that queries on a FeatureDataSource take GVCF blocks into account when computing overlap.

        // Query interval + expected variant ID(s)
        return new Object[][] {
                { new SimpleInterval("1", 1, 99), Collections.<String>emptyList() },
                { new SimpleInterval("1", 50, 100), Arrays.asList("aa") },
                { new SimpleInterval("1", 50, 150), Arrays.asList("aa") },
                { new SimpleInterval("1", 100, 100), Arrays.asList("aa") },
                { new SimpleInterval("1", 100, 150), Arrays.asList("aa") },
                { new SimpleInterval("1", 100, 200), Arrays.asList("aa") },
                { new SimpleInterval("1", 150, 200), Arrays.asList("aa") },
                { new SimpleInterval("1", 150, 250), Arrays.asList("aa") },
                { new SimpleInterval("1", 200, 201), Arrays.asList("aa") },
                { new SimpleInterval("1", 201, 3000), Collections.<String>emptyList() }
        };
    }

    @Test(dataProvider = "GVCFQueryTestData")
    public void testQueryGVCF( final SimpleInterval queryInterval, final List<String> expectedVariantIDs ) {
        try ( FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_GVCF) ) {
            final List<VariantContext> queryResults = featureSource.queryAndPrefetch(queryInterval);
            checkVariantQueryResults(queryResults, expectedVariantIDs, queryInterval);
        }
    }

    /**************************************************
     * Direct testing on the FeatureCache inner class
     **************************************************/

    @SuppressWarnings("overrides")  // because I don't want to implement hashCode() but do need an equals() here
    private static class ArtificialTestFeature implements Feature {
        private String chr;
        private int start;
        private int end;

        public ArtificialTestFeature( final String chr, final int start, final int end ) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }
        //suppressing deprecation; function required because it's part of the implemented class
        @Override
        @SuppressWarnings("deprecation")
        @Deprecated
        public String getChr() {
            return getContig();
        }

        @Override
        public String getContig() {
            return chr;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }

        @Override
        public boolean equals( Object other ) {
            if ( other == null || ! (other instanceof ArtificialTestFeature) ) {
                return false;
            }

            ArtificialTestFeature otherFeature = (ArtificialTestFeature)other;
            return chr.equals(otherFeature.getContig()) && start == otherFeature.getStart() && end == otherFeature.getEnd();
        }

        @Override
        public String toString() {
            return chr + ":" + start + "-" + end;   // (to improve output on test failures involving this class)
        }
    }

    private FeatureCache<ArtificialTestFeature> initializeFeatureCache( final List<ArtificialTestFeature> features, final String cacheContig, final int cacheStart, final int cacheEnd ) {
        FeatureCache<ArtificialTestFeature> cache = new FeatureCache<>();

        cache.fill(features.iterator(), new SimpleInterval(cacheContig, cacheStart, cacheEnd));
        return cache;
    }

    @DataProvider(name = "FeatureCacheFillDataProvider")
    public Object[][] getFeatureCacheFillData() {
        return new Object[][] {
                { Arrays.asList(new ArtificialTestFeature("1", 1, 100), new ArtificialTestFeature("1", 50, 150),
                                new ArtificialTestFeature("1", 200, 300), new ArtificialTestFeature("1", 350, 400)),
                  "1", 1, 400 },
                { Arrays.asList(new ArtificialTestFeature("1", 1, 100)), "1", 1, 100 },
                { Collections.<ArtificialTestFeature>emptyList(), "1", 1, 1 }
        };
    }

    @Test(dataProvider = "FeatureCacheFillDataProvider")
    public void testCacheFill( final List<ArtificialTestFeature> features, final String cacheContig, final int cacheStart, final int cacheEnd) {
        FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(features, cacheContig, cacheStart, cacheEnd);

        List<ArtificialTestFeature> cachedFeatures = cache.getCachedFeaturesUpToStopPosition(cacheEnd);
        Assert.assertEquals(cache.getContig(), cacheContig, "Wrong contig reported by cache after fill");
        Assert.assertEquals(cache.getCacheStart(), cacheStart, "Wrong start position reported by cache after fill");
        Assert.assertEquals(cache.getCacheEnd(), cacheEnd, "Wrong stop position reported by cache after fill");
        Assert.assertEquals(cachedFeatures, features, "Wrong Features in cache after fill()");
    }

    @DataProvider(name = "FeatureCacheHitDetectionDataProvider")
    public Object[][] getFeatureCacheHitDetectionData() {
        List<ArtificialTestFeature> features = Arrays.asList(new ArtificialTestFeature("1", 1, 100),
                                                             new ArtificialTestFeature("1", 50, 150),
                                                             new ArtificialTestFeature("1", 200, 300));
        FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(features, "1", 50, 250);

        return new Object[][] {
                // Exact match for cache boundaries
                { cache, new SimpleInterval("1", 50, 250), true },
                // Interval completely contained within cache boundaries
                { cache, new SimpleInterval("1", 100, 200), true },
                // Interval left-aligned with cache boundaries
                { cache, new SimpleInterval("1", 50, 100), true },
                // Interval right-aligned with cache boundaries
                { cache, new SimpleInterval("1", 200, 250), true },
                // Interval overlaps, but is off the left edge of cache boundaries
                { cache, new SimpleInterval("1", 49, 100), false },
                // Interval overlaps, but is off the right edge of cache boundaries
                { cache, new SimpleInterval("1", 200, 251), false },
                // Interval does not overlap and is to the left of cache boundaries
                { cache, new SimpleInterval("1", 1, 40), false },
                // Interval does not overlap and is to the right of cache boundaries
                { cache, new SimpleInterval("1", 300, 350), false },
                // Interval is on different contig
                { cache, new SimpleInterval("2", 50, 250), false }
        };
    }

    @Test(dataProvider = "FeatureCacheHitDetectionDataProvider")
    public void testCacheHitDetection( final FeatureCache<ArtificialTestFeature> cache,
                                       final SimpleInterval testInterval, final boolean cacheHitExpectedResult ) {
        Assert.assertEquals(cache.cacheHit(testInterval), cacheHitExpectedResult,
                            "Cache hit detection failed for interval " + testInterval);
    }

    @DataProvider(name = "FeatureCacheTrimmingDataProvider")
    public Object[][] getFeatureCacheTrimmingData() {
        // Features are required to always be sorted by start position, but stop positions need not be sorted.
        // This complicates cache trimming.
        List<ArtificialTestFeature> feats = Arrays.asList(
                new ArtificialTestFeature("1", 1, 1),     // Feature 0
                new ArtificialTestFeature("1", 1, 100),   // Feature 1
                new ArtificialTestFeature("1", 1, 1),     // Feature 2
                new ArtificialTestFeature("1", 1, 50),    // Feature 3
                new ArtificialTestFeature("1", 1, 3),     // Feature 4
                new ArtificialTestFeature("1", 1, 5),     // Feature 5
                new ArtificialTestFeature("1", 5, 5),     // Feature 6
                new ArtificialTestFeature("1", 5, 50),    // Feature 7
                new ArtificialTestFeature("1", 5, 10),    // Feature 8
                new ArtificialTestFeature("1", 50, 100),  // Feature 9
                new ArtificialTestFeature("1", 50, 50),   // Feature 10
                new ArtificialTestFeature("1", 50, 200),  // Feature 11
                new ArtificialTestFeature("1", 100, 100), // Feature 12
                new ArtificialTestFeature("1", 100, 110), // Feature 13
                new ArtificialTestFeature("1", 100, 200), // Feature 14
                new ArtificialTestFeature("1", 100, 150), // Feature 15
                new ArtificialTestFeature("1", 100, 199)  // Feature 16
        );
        FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(feats, "1", 1, 200);

        // Pairing of start position to which to trim the cache with the List of Features we expect to see
        // in the cache after trimming
        List<Pair<Integer, List<ArtificialTestFeature>>> trimOperations = Arrays.asList(
                Pair.of(1, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(2, Arrays.asList(feats.get(1), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(3, Arrays.asList(feats.get(1), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(4, Arrays.asList(feats.get(1), feats.get(3), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(5, Arrays.asList(feats.get(1), feats.get(3), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(6, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(10, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(11, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(50, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(51, Arrays.asList(feats.get(1), feats.get(9), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(100, Arrays.asList(feats.get(1), feats.get(9), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(101, Arrays.asList(feats.get(11), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(111, Arrays.asList(feats.get(11), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(151, Arrays.asList(feats.get(11), feats.get(14), feats.get(16))),
                Pair.of(151, Arrays.asList(feats.get(11), feats.get(14), feats.get(16))),
                Pair.of(200, Arrays.asList(feats.get(11), feats.get(14)))
        );

        return new Object[][] {
                { cache, trimOperations }
        };
    }

    @Test(dataProvider = "FeatureCacheTrimmingDataProvider")
    public void testCacheTrimming( final FeatureCache<ArtificialTestFeature> cache, final List<Pair<Integer, List<ArtificialTestFeature>>> trimOperations ) {
        // Repeatedly trim the cache to ever-increasing start positions, and verify after each trim operation
        // that the cache holds the correct Features in the correc order
        for ( Pair<Integer, List<ArtificialTestFeature>> trimOperation : trimOperations ) {
            final int trimPosition = trimOperation.getLeft();
            final List<ArtificialTestFeature> expectedFeatures = trimOperation.getRight();

            cache.trimToNewStartPosition(trimPosition);

            final List<ArtificialTestFeature> actualFeatures = cache.getCachedFeaturesUpToStopPosition(cache.getCacheEnd());
            Assert.assertEquals(actualFeatures, expectedFeatures, "Wrong Features in cache after trimming start position to " + trimPosition);
        }
    }

    @DataProvider(name = "FeatureCacheRetrievalDataProvider")
    public Object[][] getFeatureCacheRetrievalData() {
        List<ArtificialTestFeature> feats = Arrays.asList(
                new ArtificialTestFeature("1", 1, 1),      // Feature 0
                new ArtificialTestFeature("1", 1, 100),    // Feature 1
                new ArtificialTestFeature("1", 5, 5),      // Feature 2
                new ArtificialTestFeature("1", 10, 10),    // Feature 3
                new ArtificialTestFeature("1", 10, 100),   // Feature 4
                new ArtificialTestFeature("1", 50, 50),    // Feature 5
                new ArtificialTestFeature("1", 51, 55),    // Feature 6
                new ArtificialTestFeature("1", 52, 52),    // Feature 7
                new ArtificialTestFeature("1", 55, 60),    // Feature 8
                new ArtificialTestFeature("1", 75, 75),    // Feature 9
                new ArtificialTestFeature("1", 80, 100)    // Feature 10
        );
        FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(feats, "1", 1, 100);

        // Pairing of end position with which to bound cache retrieval with the List of Features we expect to see
        // after retrieval
        List<Pair<Integer, List<ArtificialTestFeature>>> retrievalOperations = Arrays.asList(
                Pair.of(100, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(80, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(79, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9))),
                Pair.of(80, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(75, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9))),
                Pair.of(74, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8))),
                Pair.of(54, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7))),
                Pair.of(52, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7))),
                Pair.of(51, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6))),
                Pair.of(50, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5))),
                Pair.of(49, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4))),
                Pair.of(10, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4))),
                Pair.of(9, Arrays.asList(feats.get(0), feats.get(1), feats.get(2))),
                Pair.of(5, Arrays.asList(feats.get(0), feats.get(1), feats.get(2))),
                Pair.of(4, Arrays.asList(feats.get(0), feats.get(1))),
                Pair.of(1, Arrays.asList(feats.get(0), feats.get(1)))
        );

        return new Object[][] {
                { cache, retrievalOperations }
        };
    }

    @Test(dataProvider = "FeatureCacheRetrievalDataProvider")
    public void testCacheFeatureRetrieval( final FeatureCache<ArtificialTestFeature> cache, final List<Pair<Integer, List<ArtificialTestFeature>>> retrievalOperations ) {
        for ( Pair<Integer, List<ArtificialTestFeature>> retrievalOperation: retrievalOperations ) {
            final int stopPosition = retrievalOperation.getLeft();
            final List<ArtificialTestFeature> expectedFeatures = retrievalOperation.getRight();

            final List<ArtificialTestFeature> actualFeatures = cache.getCachedFeaturesUpToStopPosition(stopPosition);
            Assert.assertEquals(actualFeatures, expectedFeatures, "Wrong Features returned in retrieval operation with stop position " + stopPosition);
        }
    }

    /**
     * Test caching a region with no Features. This should work (we should avoid going to disk
     * to look for new records when querying within such a region).
     */
    @Test
    public void testHandleCachingOfEmptyRegion() {
        FeatureCache<ArtificialTestFeature> cache = new FeatureCache<>();
        List<ArtificialTestFeature> emptyRegion = new ArrayList<>();

        cache.fill(emptyRegion.iterator(), new SimpleInterval("1", 1, 100));

        Assert.assertTrue(cache.isEmpty(), "Cache should be empty");
        Assert.assertTrue(cache.cacheHit(new SimpleInterval("1", 1, 100)), "Unexpected cache miss");
        Assert.assertTrue(cache.cacheHit(new SimpleInterval("1", 2, 99)), "Unexpected cache miss");

        Assert.assertEquals(cache.getCachedFeaturesUpToStopPosition(100), emptyRegion, "Should get back empty List for empty region");
        cache.trimToNewStartPosition(2);
        Assert.assertTrue(cache.cacheHit(new SimpleInterval("1", 2, 100)), "Unexpected cache miss");
        Assert.assertEquals(cache.getCachedFeaturesUpToStopPosition(100), emptyRegion, "Should get back empty List for empty region");
    }

    /*********************************************************
     * End of direct testing on the FeatureCache inner class
     *********************************************************/

}
