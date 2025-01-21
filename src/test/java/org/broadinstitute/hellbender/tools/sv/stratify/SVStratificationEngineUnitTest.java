package org.broadinstitute.hellbender.tools.sv.stratify;

import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class SVStratificationEngineUnitTest extends GATKBaseTest {

    private static final GATKPath CONFIG_FILE_PATH = new GATKPath(toolsTestDir + "/sv/sv_stratify_config.tsv");

    private static final String CONTEXT_1_NAME = "context1";
    private static final String CONTEXT_2_NAME = "context2";

    private static final List<Locatable> CONTEXT_1_INTERVALS = Lists.newArrayList(new SimpleInterval("chr1", 1000, 2000));
    private static final List<Locatable> CONTEXT_2_INTERVALS = Lists.newArrayList(new SimpleInterval("chr2", 1000, 2000));
    
    private static SVStratificationEngine makeDefaultEngine() {
        return new SVStratificationEngine(SVTestUtils.hg38Dict);
    }

    @Test
    public void testAddContext() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        Assert.assertNotNull(engine.getTrackIntervals(CONTEXT_1_NAME));
        Assert.assertNull(engine.getTrackIntervals(CONTEXT_2_NAME));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddDuplicateContext() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_2_INTERVALS);
    }

    @Test
    public void testNoContexts() {
        final SVStratificationEngine engine = makeDefaultEngine();
        Assert.assertTrue(engine.getStrata().isEmpty());
    }

    @Test
    public void testAddStratification() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addStratification("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, Collections.singleton(CONTEXT_1_NAME));
        final Collection<SVStratificationEngine.Stratum> stratificationCollection = engine.getStrata();
        Assert.assertNotNull(stratificationCollection);
        Assert.assertEquals(stratificationCollection.size(), 1);
        final SVStratificationEngine.Stratum stratification = stratificationCollection.iterator().next();
        Assert.assertNotNull(stratification);
        Assert.assertEquals(stratification.getSvType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertNotNull(stratification.getMinSize());
        Assert.assertEquals(stratification.getMinSize().intValue(), 50);
        Assert.assertNotNull(stratification.getMaxSize());
        Assert.assertEquals(stratification.getMaxSize().intValue(), 500);
        Assert.assertEquals(stratification.getTrackNames().size(), 1);
        Assert.assertEquals(stratification.getTrackNames().get(0), CONTEXT_1_NAME);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddStratificationBadMinSize() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addStratification("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, -1, 500, Collections.emptySet());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddStratificationBadMaxSize() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addStratification("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, -1, Collections.emptySet());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddStratificationBadMaxSizeInfinity() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addStratification("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Integer.MAX_VALUE, Collections.emptySet());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddStratificationMaxEqualToMin() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addStratification("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 50, Collections.emptySet());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddStratificationMaxLessThanMin() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addStratification("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 49, Collections.emptySet());
    }

    @Test
    public void testCreate() {
        final Map<String, List<Locatable>> map = new HashMap<>();
        map.put(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        map.put(CONTEXT_2_NAME, CONTEXT_2_INTERVALS);
        final SVStratificationEngine engine = SVStratificationEngine.create(map, CONFIG_FILE_PATH, SVTestUtils.hg38Dict);
        Assert.assertNotNull(engine);
        Assert.assertNotNull(engine.getTrackIntervals(CONTEXT_1_NAME));
        Assert.assertEquals(engine.getStrata().size(), 7);
    }

    @DataProvider(name="testGetMatchVariantsData")
    public Object[][] testGetMatchVariantsData() {
        return new Object[][] {

                // DEL

                // Outside context interval
                { "chr1", 100, "chr1", 200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                { "chr1", 2000, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                // Simple match
                { "chr1", 1100, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 900, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 900, "chr1", 1900, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 1100, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 800, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 999, "chr1", 2001, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr2", 1100, "chr2", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                // Wrong contig
                { "chr3", 1100, "chr3", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                // Barely match
                { "chr1", 1000, "chr1", 3001, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 2, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                { "chr1", 500, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                // Barely miss overlap threshold
                { "chr1", 1000, "chr1", 3002, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                // Barely large enough
                { "chr1", 1100, "chr1", 1149, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5k_both" },
                // Too small
                { "chr1", 1100, "chr1", 1148, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },

                // INV (null context)

                // Right size
                { "chr1", 1001, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.INV, null, "INV_gt1kb" },
                { "chr1", 4001, "chr1", 5000, GATKSVVCFConstants.StructuralVariantAnnotationType.INV, null, "INV_gt1kb" },
                { "chr2", 10000, "chr2", 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.INV, null, "INV_gt1kb" },
                // Too small
                { "chr1", 1001, "chr1", 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.INV, null, null },
                { "chr1", 100, "chr1", 200, GATKSVVCFConstants.StructuralVariantAnnotationType.INV, null, null },

                // INS

                // In context
                { "chr1", 1100, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 100, "INS_context1" },
                // SVLEN should not matter
                { "chr1", 1100, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 1, "INS_context1" },
                { "chr1", 1100, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 10000, "INS_context1" },
                // Out of context
                { "chr1", 100, "chr1", 100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 100, null },
                // Out of size range for context2
                { "chr2", 1100, "chr2", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 1000, null },
                { "chr2", 1100, "chr2", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 400, null },

                // BND

                // Both ends
                { "chr1", 1000, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, "BND_context1" },
                { "chr1", 2000, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, "BND_context1" },
                // One end only
                { "chr1", 500, "chr1", 900, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, null },
                { "chr1", 1500, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, null },
                // No ends
                { "chr1", 500, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, null },

                // BND (same as CTX)

                // Both ends
                { "chr1", 1000, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, "CTX_context1" },
                { "chr1", 2000, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, "CTX_context1" },
                // One end only
                { "chr1", 500, "chr1", 900, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, null },
                { "chr1", 1500, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, null },
                // No ends
                { "chr1", 500, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, null },
        };
    }

    @Test(dataProvider = "testGetMatchVariantsData")
    public void testGetMatchVariants(final String chromA, final int posA, final String chromB, final int posB,
                                     final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                                     final Integer svlen,
                                     final String expectedStratName) {
        final Map<String, List<Locatable>> map = new HashMap<>();
        map.put(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        map.put(CONTEXT_2_NAME, CONTEXT_2_INTERVALS);
        final SVStratificationEngine engine = SVStratificationEngine.create(map, CONFIG_FILE_PATH, SVTestUtils.hg38Dict);
        final SVCallRecord record;
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            record = SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates(chromA, posA, svlen);
        } else {
            record = SVTestUtils.newCallRecordWithCoordinatesAndType("record", chromA, posA, chromB, posB, svType);
        }
        final Collection<SVStratificationEngine.Stratum> result = engine.getMatches(record, 0.5, 0, 2);
        if (expectedStratName == null) {
            Assert.assertTrue(result.isEmpty());
        } else {
            Assert.assertFalse(result.isEmpty());
            Assert.assertEquals(result.iterator().next().getName(), expectedStratName);
        }
    }

    // Not supported
    @Test(expectedExceptions = GATKException.class)
    public void testGetMatchVariantsCpx() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addTrack("context3", Lists.newArrayList(new SimpleInterval("chr1", 1500, 2500)));
        engine.addStratification("strat1", GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, 50, 500, Collections.singleton("context1"));
        engine.addStratification("strat2", GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, 50, 500, Collections.singleton("context3"));
        final SVCallRecord record = SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1800, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX);
        // Should throw error
        engine.getMatches(record, 0.5, 0, 2);
    }

    @Test
    public void testGetMatchVariantsMultiple() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addTrack("context3", Lists.newArrayList(new SimpleInterval("chr1", 1500, 2500)));
        engine.addStratification("strat1", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, Collections.singleton("context1"));
        engine.addStratification("strat2", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, Collections.singleton("context3"));
        final SVCallRecord record = SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1800, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final Collection<SVStratificationEngine.Stratum> result = engine.getMatches(record, 0.5, 0, 2);
        final List<String> names = result.stream().map(SVStratificationEngine.Stratum::getName).collect(Collectors.toList());
        Assert.assertTrue(names.contains("strat1"));
        Assert.assertTrue(names.contains("strat2"));
    }

    @Test
    public void testGetMatchVariantsNullContexts() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addStratification("strat1", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, Collections.emptySet());
        final SVCallRecord record = SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr2", 1800, "chr2", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final Collection<SVStratificationEngine.Stratum> result = engine.getMatches(record, 0.5, 0, 2);
        final List<String> names = result.stream().map(SVStratificationEngine.Stratum::getName).collect(Collectors.toList());
        Assert.assertEquals(names.size(), 1);
        Assert.assertEquals(names.get(0), "strat1");
    }

    @Test
    public void testGetMatchVariantsNoEngineContexts() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addStratification("strat1", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, Collections.emptySet());
        final SVCallRecord record = SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr2", 1800, "chr2", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final Collection<SVStratificationEngine.Stratum> result = engine.getMatches(record, 0.5, 0, 2);
        final List<String> names = result.stream().map(SVStratificationEngine.Stratum::getName).collect(Collectors.toList());
        Assert.assertEquals(names.size(), 1);
        Assert.assertEquals(names.get(0), "strat1");
    }

    @Test
    public void testTestAddStratificationInnerClass() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        final SVStratificationEngine.Stratum stratification = engine.new Stratum("strat", GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, Collections.singleton(CONTEXT_1_NAME));
        engine.addStratification(stratification);
        final Collection<SVStratificationEngine.Stratum> stratificationCollection = engine.getStrata();
        Assert.assertNotNull(stratificationCollection);
        Assert.assertEquals(stratificationCollection.size(), 1);
        final SVStratificationEngine.Stratum stratificationOut = stratificationCollection.iterator().next();
        Assert.assertNotNull(stratificationOut);
        Assert.assertEquals(stratificationOut.getSvType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertNotNull(stratificationOut.getMinSize());
        Assert.assertEquals(stratificationOut.getMinSize().intValue(), 50);
        Assert.assertNotNull(stratificationOut.getMaxSize());
        Assert.assertEquals(stratificationOut.getMaxSize().intValue(), 500);
        Assert.assertEquals(stratificationOut.getTrackNames().size(), 1);
        Assert.assertEquals(stratificationOut.getTrackNames().get(0), CONTEXT_1_NAME);
    }

    @Test
    public void testMatchesType() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                100, 500,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesType(SVTestUtils.newCallRecordWithLengthAndType(null, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesType(SVTestUtils.newCallRecordWithLengthAndType(null, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP)));
    }

    @Test
    public void testMatchesSizeSimple() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                100, 500,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(499, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(50, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(500, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeNoMin() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, 500,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(499, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(500, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeNoMax() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                50, null,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(49, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(Integer.MAX_VALUE - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeNoMinOrMax() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, null,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(Integer.MAX_VALUE - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeInsertion() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                100, 500,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(499, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(50, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(500, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
    }

    @Test
    public void testMatchesSizeInsertionNullLength() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                0, Integer.MAX_VALUE - 1,
                Collections.emptySet()
        );
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(null, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
    }

    @Test
    public void testMatchesSizeInsertionNullLength2() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                null, null,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(null, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
    }

    @Test
    public void testMatchesSizeBnd() {
        final SVStratificationEngine engine = makeDefaultEngine();
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                Collections.emptySet()
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newBndCallRecordWithStrands(true, false)));
    }


    @DataProvider(name="testMatchesContextDelData")
    public Object[][] testMatchesContextDelData() {
        return new Object[][] {
                // Outside context interval
                { "chr1", 1000, 1500, 0.5, 0, true },
                { "chr1", 500, 1500, 0.5, 0, true },
                { "chr1", 499, 1499, 0.5, 0, false },
                { "chr1", 900, 1300, 0.5, 1, true },
                { "chr1", 1999, 2000000, 0, 1, true },
                { "chr1", 500, 600, 0, 2, false },
                { "chr1", 500, 1100, 0, 2, false },
                { "chr1", 1100, 1200, 0, 2, true },
                { "chr1", 1100, 1200, 1, 2, true }
        };
    }

    @Test(dataProvider = "testMatchesContextDelData")
    public void testMatchesContextDel(final String chrom, final int start, final int end,
                                      final double overlapFraction, final int numBreakpointOverlaps,
                                      final boolean expected) {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                Collections.singleton(CONTEXT_1_NAME)
        );
        Assert.assertEquals(strat.matchesTracks(SVTestUtils.newCallRecordWithCoordinatesAndType("record", chrom, start, chrom, end, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                overlapFraction, numBreakpointOverlaps, 1), expected);
    }

    @DataProvider(name="testMatchesContextInsData")
    public Object[][] testMatchesContextInsData() {
        return new Object[][] {
                // Outside context interval
                { "chr1", 1100, 100, 0.1, 0, true },
                { "chr1", 1100, 100000, 0.1, 0, true },
                { "chr1", 999, 100, 0.1, 0, false }
        };
    }

    @Test(dataProvider = "testMatchesContextInsData")
    public void testMatchesContextIns(final String chrom, final int start, final int length,
                                      final double overlapFraction, final int numBreakpointOverlaps,
                                      final boolean expected) {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                Collections.singleton(CONTEXT_1_NAME)
        );
        Assert.assertEquals(strat.matchesTracks(SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates(chrom, start, length),
                overlapFraction, numBreakpointOverlaps, 1), expected);
    }

    @DataProvider(name="testMatchesContextBndData")
    public Object[][] testMatchesContextBndData() {
        return new Object[][] {
                { "chr1", 999, "chr1", 2001, 1, false },
                { "chr1", 1000, "chr1", 1200, 1, true },
                { "chr1", 1000, "chr1", 50000, 1, true },
                { "chr1", 1000, "chr1", 1000, 1, true },
                { "chr1", 500, "chr1", 1000, 1, true },
                { "chr1", 1000, "chr1", 1999, 2, true },
                { "chr1", 1000, "chr1", 2000, 2, true },
                { "chr1", 1000, "chr2", 1000, 2, false },
                { "chr1", 1000, "chr1", 2001, 2, false },
                { "chr1", 999, "chr1", 1000, 2, false }
        };
    }

    @Test(dataProvider = "testMatchesContextBndData")
    public void testMatchesContextBnd(final String chromA, final int posA, final String chromB, final int posB,
                                      final int numBreakpointOverlapsInterchrom, final boolean expected) {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                Collections.singleton(CONTEXT_1_NAME)
        );
        Assert.assertEquals(strat.matchesTracks(SVTestUtils.newCallRecordWithCoordinatesAndType("record", chromA, posA, chromB, posB, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, numBreakpointOverlapsInterchrom), expected);
    }

    @DataProvider(name="testCountAnyContextOverlapData")
    public Object[][] testCountAnyContextOverlapData() {
        return new Object[][] {
                { "chr1", 500, 1500, 1 },
                { "chr1", 1000, 2000, 1 },
                { "chr1", 1500, 2500, 1 },
                { "chr1", 500, 2500, 1 },
                { "chr1", 1100, 1900, 1 },
                { "chr1", 999, 999, 0 },
                { "chr1", 999, 1000, 1 },
                { "chr1", 1000, 1000, 1 },
                { "chr1", 1000, 1001, 1 },
                { "chr2", 1000, 1001, 0 },
                { "chr1", 1999, 2000, 1 },
                { "chr1", 2000, 2000, 1 },
                { "chr1", 2001, 2001, 0 }
        };
    }

    @Test(dataProvider = "testCountAnyContextOverlapData")
    public void testCountAnyContextOverlap(final String chrom, final int start, final int end, final int expected) {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                Collections.singleton(CONTEXT_1_NAME)
        );
        Assert.assertEquals(strat.countAnyTrackOverlap(new SimpleInterval(chrom, start, end)), expected);
    }

    @DataProvider(name="testIsMutuallyExclusiveData")
    public Object[][] testIsMutuallyExclusiveData() {
        return new Object[][] {
                {GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null, null, null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null, null, null,
                        false},
        };
    }

    @Test
    public void testGetters() {
        final SVStratificationEngine engine = makeDefaultEngine();
        engine.addTrack(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        final SVStratificationEngine.Stratum strat = engine.new Stratum(
                "strat",
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                50, 500,
                Collections.singleton(CONTEXT_1_NAME)
        );
        Assert.assertEquals(strat.getTrackNames().size(), 1);
        Assert.assertEquals(strat.getTrackNames().get(0), CONTEXT_1_NAME);
        Assert.assertEquals(strat.getSvType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertEquals(strat.getMinSize(), Integer.valueOf(50));
        Assert.assertEquals(strat.getMaxSize(), Integer.valueOf(500));
        Assert.assertEquals(strat.getName(), "strat");
    }
}