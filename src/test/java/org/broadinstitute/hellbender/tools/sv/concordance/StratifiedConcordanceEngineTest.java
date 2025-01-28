package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.util.Combinations;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.ClusteringParameters;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.walkers.sv.SVStratify;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class StratifiedConcordanceEngineTest extends BaseTest {

    private static final SVCallRecord TEST_VARIANT_1 = SVTestUtils.newCallRecordWithCoordinatesAndType("record1", "chr1", 1000, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
    
    private static final String STRATIFICATION_1 = "strat1";
    private static final String STRATIFICATION_2 = "strat2";
    private static final String STRATIFICATION_3 = "strat3";
    private static final String STRATIFICATION_4 = "strat4";
    private static final String TRACK_A = "trackA";
    private static final String TRACK_B = "trackB";
    private static final String TRACK_C = "trackC";
    private static final String TRACK_D = "trackD";

    private static StratifiedConcordanceEngine createConcordanceEngine() {

        final SVStratificationEngine stratificationEngine = new SVStratificationEngine(SVTestUtils.hg38Dict);
        // Diagram of the interval tracks on 1000-2000 (note there are additional intervals past 2000). Track D is empty.
        //          -----|----------------|--------|--------|--------|----------------|-------
        //             1000             1200     1300     1400     1500             2000
        // TRACK A       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        // TRACK B                        XXXXXXXXXXXXXXXXXX
        // TRACK C                                  XXXXXXXXXXXXXXXXX
        // TRACK D
        stratificationEngine.addTrack(TRACK_A, Arrays.asList(new SimpleInterval("chr1", 1000, 2000), new SimpleInterval("chr1", 3000, 4000)));
        stratificationEngine.addTrack(TRACK_B, Arrays.asList(new SimpleInterval("chr1", 1200, 1400), new SimpleInterval("chr1", 5000, 6000)));
        stratificationEngine.addTrack(TRACK_C, Arrays.asList(new SimpleInterval("chr1", 1300, 1500), new SimpleInterval("chr1", 7000, 8000)));
        stratificationEngine.addTrack(TRACK_D, Collections.emptyList());
        stratificationEngine.addStratification(STRATIFICATION_1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 100, 5000, Sets.newHashSet(TRACK_A));
        stratificationEngine.addStratification(STRATIFICATION_2, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 100, 5000, Sets.newHashSet(TRACK_B));
        stratificationEngine.addStratification(STRATIFICATION_3, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 100, 5000, Sets.newHashSet(TRACK_C));
        stratificationEngine.addStratification(STRATIFICATION_4, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 100, 5000, Sets.newHashSet(TRACK_D));

        final Map<String, ClosestSVFinder> clusterEngineMap = new HashMap<>();
        final SVConcordanceLinkage linkage1 = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        linkage1.setDepthOnlyParams(ClusteringParameters.createDepthParameters(0, 0, 10, 0));
        final SVConcordanceLinkage linkage2 = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        linkage2.setDepthOnlyParams(ClusteringParameters.createDepthParameters(0, 0, 20, 0));
        final SVConcordanceLinkage linkage3 = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        linkage3.setDepthOnlyParams(ClusteringParameters.createDepthParameters(0, 0, 30, 0));
        final SVConcordanceLinkage linkage4 = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        linkage4.setDepthOnlyParams(ClusteringParameters.createDepthParameters(0, 0, 0, 0));

        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        clusterEngineMap.put(STRATIFICATION_1, new ClosestSVFinder(linkage1, collapser::annotate, SVTestUtils.hg38Dict));
        clusterEngineMap.put(STRATIFICATION_2, new ClosestSVFinder(linkage2, collapser::annotate, SVTestUtils.hg38Dict));
        clusterEngineMap.put(STRATIFICATION_3, new ClosestSVFinder(linkage3, collapser::annotate, SVTestUtils.hg38Dict));
        clusterEngineMap.put(STRATIFICATION_4, new ClosestSVFinder(linkage4, collapser::annotate, SVTestUtils.hg38Dict));

        final SVClusterEngineArgumentsCollection defaultClusteringArgs = new SVClusterEngineArgumentsCollection();
        defaultClusteringArgs.depthBreakendWindow = 0;
        defaultClusteringArgs.depthOverlapFraction = 0;
        defaultClusteringArgs.depthSizeSimilarity = 0;

        final SVStratificationEngineArgumentsCollection stratArgs = new SVStratificationEngineArgumentsCollection();
        stratArgs.numBreakpointOverlaps = 1;
        stratArgs.overlapFraction = 0;

        return new StratifiedConcordanceEngine(clusterEngineMap, stratificationEngine, stratArgs, defaultClusteringArgs, collapser, SVTestUtils.hg38Dict);
    }

    @Test
    public void testIsEmpty() {
        final StratifiedConcordanceEngine engine = createConcordanceEngine();
        Assert.assertTrue(engine.isEmpty());
        engine.addTruthVariant(TEST_VARIANT_1);
        // Empty if only contains truth variants
        Assert.assertTrue(engine.isEmpty());
        engine.addEvalVariant(TEST_VARIANT_1);
        Assert.assertFalse(engine.isEmpty());
        engine.flush(true);
        Assert.assertTrue(engine.isEmpty());
    }

    @DataProvider(name = "testStratificationData")
    public Object[][] testStratificationData() {
        return new Object[][]{
                {100, 200, Collections.emptyList()},
                {1000, 2000, Lists.newArrayList(STRATIFICATION_1)},
                {1000, 1010, Collections.emptyList()}, // too small
                {1000, 2001, Lists.newArrayList(STRATIFICATION_1)},
                {999, 2000, Lists.newArrayList(STRATIFICATION_1)},
                {999, 2001, Collections.emptyList()},
                {4500, 5500, Lists.newArrayList(STRATIFICATION_2)},
                {4500, 7500, Lists.newArrayList(STRATIFICATION_3)},
                {6500, 8500, Collections.emptyList()},
                {500, 1250, Lists.newArrayList(STRATIFICATION_1, STRATIFICATION_2)},
                {500, 1450, Lists.newArrayList(STRATIFICATION_1, STRATIFICATION_3)},
                {500, 1350, Lists.newArrayList(STRATIFICATION_1, STRATIFICATION_2, STRATIFICATION_3)},
        };
    }

    @Test(dataProvider= "testStratificationData")
    public void testSingleRecord(final int positionA, final int positionB, final List<String> expectedNonDefault) {
        final StratifiedConcordanceEngine engine = createConcordanceEngine();
        final SVCallRecord record = SVTestUtils.newDeletionRecordWithCoordinates("record1", "chr1", positionA, positionB);
        engine.addEvalVariant(record);
        final Collection<VariantContext> result = engine.flush(true);
        Assert.assertEquals(result.size(), 1);
        final VariantContext outputRecord = result.iterator().next();
        Assert.assertNotNull(outputRecord);
        final List<String> observed = outputRecord.getAttributeAsStringList(GATKSVVCFConstants.STRATUM_INFO_KEY, "");
        final List<String> expected = new ArrayList<>(expectedNonDefault.size() + 1);
        expected.addAll(expectedNonDefault);
        expected.add(SVStratify.DEFAULT_STRATUM);
        Assert.assertEquals(observed, expected.stream().sorted().collect(Collectors.toUnmodifiableList()));
    }

    @Test
    public void testSoftFlush() {
        final StratifiedConcordanceEngine engine = createConcordanceEngine();
        final SVCallRecord eval1 = SVTestUtils.newDeletionRecordWithCoordinates("eval1", "chr1", 1100, 1900);
        final SVCallRecord eval2 = SVTestUtils.newDeletionRecordWithCoordinates("eval2", "chr1", 1300, 1900);
        final SVCallRecord eval3 = SVTestUtils.newDeletionRecordWithCoordinates("eval3", "chr1", 4000, 5000);
        final SVCallRecord truth1 = SVTestUtils.newDeletionRecordWithCoordinates("truth1", "chr1", 1200, 1900);
        final SVCallRecord truth2 = SVTestUtils.newDeletionRecordWithCoordinates("truth2", "chr1", 4000, 5000);
        engine.addEvalVariant(eval1);
        engine.addTruthVariant(truth1);
        engine.addEvalVariant(eval2);
        engine.addTruthVariant(truth2);
        engine.addEvalVariant(eval3);
        final List<VariantContext> result = engine.flush(false).stream().sorted(Comparator.comparingInt(VariantContext::getStart)).collect(Collectors.toUnmodifiableList());

        Assert.assertEquals(result.size(), 2);
        final VariantContext outputRecord1 = result.get(0);
        final VariantContext outputRecord2 = result.get(1);
        Assert.assertNotNull(outputRecord1);
        Assert.assertNotNull(outputRecord2);
        Assert.assertEquals(outputRecord1.getID(), eval1.getId());
        Assert.assertEquals(outputRecord2.getID(), eval2.getId());
    }


    @DataProvider(name = "testConcordanceSimpleData")
    public Object[][] testConcordanceSimpleData() {
        return new Object[][]{
                // default strat exact
                {100, 400, 100, 400, true},
                // default strat not within tolerance
                {100, 400, 101, 400, false},
                {100, 100, 100, 399, false},
                // strat 1 exact match
                {1100, 1900, 1100, 1900, true},
                // strat 1 within 10 bp
                {1100, 1900, 1105, 1905, true},
                // strat 1 within 10 bp
                {1100, 1900, 1105, 1910, true},
                // strat 1 and default, not within tolerance
                {1100, 1900, 1105, 1911, false},
                // strat 2, not within default/strat 1 tolerance
                {1210, 1900, 1225, 1900, true},
                // strat 3, not within default/strat 1/2 tolerance
                {1310, 1900, 1335, 1900, true}
        };
    }

    @Test(dataProvider= "testConcordanceSimpleData")
    public void testConcordanceSimple(final int evalPos, final int evalEnd, final int truthPos, final int truthEnd, final boolean truePositive) {
        final StratifiedConcordanceEngine engine = createConcordanceEngine();
        final SVCallRecord eval = SVTestUtils.newDeletionRecordWithCoordinates("eval1", "chr1", evalPos, evalEnd);
        final SVCallRecord truth = SVTestUtils.newDeletionRecordWithCoordinates("truth1", "chr1", truthPos, truthEnd);
        engine.addEvalVariant(eval);
        engine.addTruthVariant(truth);
        final Collection<VariantContext> result = engine.flush(true);

        Assert.assertEquals(result.size(), 1);
        final VariantContext outputRecord = result.iterator().next();
        Assert.assertNotNull(outputRecord);
        Assert.assertEquals(outputRecord.getID(), eval.getId());
        Assert.assertEquals(outputRecord.getAttribute(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE), truePositive ? ConcordanceState.TRUE_POSITIVE.getAbbreviation() : ConcordanceState.FALSE_POSITIVE.getAbbreviation());
        Assert.assertEquals(outputRecord.getAttribute(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO), truePositive ? truth.getId() : null);
    }


    @DataProvider(name = "testConcordanceMultiData")
    public Object[][] testConcordanceMultiData() {
        return new Object[][]{
                // default strat exact
                {100, 400, 100, 400, 101, 401, true},
                {101, 401, 100, 400, 101, 401, false},
                // strat 1 match first
                {1100, 1900, 1100, 1900, 1101, 1900, true},
                // strat 1 match second
                {1101, 1900, 1100, 1900, 1101, 1900, false},
                // strat 2 match first
                {1310, 1900, 1310, 1900, 1311, 1900, true},
                // strat 2 match second
                {1310, 1900, 1310, 2000, 1311, 1900, false},
                // strat 3 match first
                {1410, 1900, 1435, 1900, 1450, 1900, true},
                // strat 3 match second
                {1410, 1900, 1410, 1850, 1410, 1925, false},
                // matches all strats but fails matching in 3, but 2 is still found
                {1410, 1900, 1410, 1920, 1410, 1950, true},
                {1410, 1900, 1410, 1920, 1410, 1921, true},
                {1410, 1900, 1410, 1950, 1410, 1920, false},
                {1410, 1900, 1410, 1921, 1410, 1920, false},
        };
    }

    @Test(dataProvider= "testConcordanceMultiData")
    public void testConcordanceMulti(final int evalPos, final int evalEnd, final int truth1Pos, final int truth1End, final int truth2Pos, final int truth2End, final boolean truth1Match) {
        final StratifiedConcordanceEngine engine = createConcordanceEngine();
        final SVCallRecord eval = SVTestUtils.newDeletionRecordWithCoordinates("eval1", "chr1", evalPos, evalEnd);
        final SVCallRecord truth1 = SVTestUtils.newDeletionRecordWithCoordinates("truth1", "chr1", truth1Pos, truth1End);
        final SVCallRecord truth2 = SVTestUtils.newDeletionRecordWithCoordinates("truth2", "chr1", truth2Pos, truth2End);
        if (evalPos < truth1Pos) {
            engine.addEvalVariant(eval);
            engine.addTruthVariant(truth1);
        } else {
            engine.addTruthVariant(truth1);
            engine.addEvalVariant(eval);
        }
        engine.addTruthVariant(truth2);
        final Collection<VariantContext> result = engine.flush(true);

        Assert.assertEquals(result.size(), 1);
        final VariantContext outputRecord = result.iterator().next();
        Assert.assertNotNull(outputRecord);
        Assert.assertEquals(outputRecord.getID(), eval.getId());
        Assert.assertEquals(outputRecord.getAttribute(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO), truth1Match ? truth1.getId() : truth2.getId());
    }

    // Tests all combinations of matching from 3 groups
    @Test
    public void testItemTracker() {
        // In order of priority (most priority to least)
        final List<String> groups = Lists.newArrayList(STRATIFICATION_1, STRATIFICATION_3, SVStratify.DEFAULT_STRATUM);
        for (int k = 0; k < groups.size(); k++) {
            final Combinations combinations = new Combinations(groups.size(), k);
            for (int[] c : combinations) {
                final StratifiedConcordanceEngine.ItemTracker itemTracker = new StratifiedConcordanceEngine.ItemTracker(groups);
                Assert.assertEquals(itemTracker.getGroups().size(), groups.size());
                Assert.assertTrue(new HashSet<>(groups).containsAll(itemTracker.getGroups()));
                Assert.assertNull(itemTracker.getOutput());
                final Set<Integer> combinationSet = Arrays.stream(c).boxed().collect(Collectors.toUnmodifiableSet());
                for (int i = 0; i < groups.size(); i++) {
                    final String group = groups.get(i);
                    final Map<String, Object> attr = new HashMap<>();
                    attr.put(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, combinationSet.contains(i) ? ConcordanceState.TRUE_POSITIVE.getAbbreviation() : ConcordanceState.FALSE_POSITIVE.getAbbreviation());
                    attr.put(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, "var_" + group);
                    final SVCallRecord record = SVTestUtils.newDeletionRecordWithAttributes(attr);
                    Assert.assertFalse(itemTracker.allEjected());
                    itemTracker.eject(group, record);
                }
                Assert.assertTrue(itemTracker.allEjected());
                final SVCallRecord out = itemTracker.getOutput();
                if (c.length == 0) {
                    Assert.assertEquals(out.getAttributes().get(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE), ConcordanceState.FALSE_POSITIVE.getAbbreviation());
                } else {
                    // Highest priority group is expected
                    final int minIndex = c[MathUtils.minElementIndex(c)];
                    Assert.assertEquals(out.getAttributes().get(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE), ConcordanceState.TRUE_POSITIVE.getAbbreviation());
                    final String minGroup = groups.get(minIndex);
                    Assert.assertEquals(out.getAttributes().get(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO), "var_" + minGroup);
                }
            }
        }
    }
}