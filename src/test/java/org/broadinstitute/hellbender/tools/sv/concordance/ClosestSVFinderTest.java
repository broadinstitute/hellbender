package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class ClosestSVFinderTest {

    private static void assertConcordanceMembers(final SVCallRecord record, final String expectedId) {
        final String closestId = (String) record.getAttributes().get(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO);
        Assert.assertEquals(closestId, expectedId);
    }

    @Test
    public void testMetrics() {
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        final ClosestSVFinder engine = new ClosestSVFinder(linkage, collapser::annotate, SVTestUtils.hg38Dict);

        final SVCallRecord eval1 = SVTestUtils.makeRecord(
                "eval1",
                "chr1",
                100000,
                true,
                "chr1",
                100001,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList(
                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.REF_N)),
                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS))
                )
        );

        final SVCallRecord truth1 = SVTestUtils.makeRecord(
                "truth1",
                "chr1",
                100010,
                true,
                "chr1",
                100011,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList(
                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS))
                )
        );

        final SVCallRecord truth2 = SVTestUtils.makeRecord(
                "truth2",
                "chr1",
                300010,
                true,
                "chr1",
                300011,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList(
                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS))
                )
        );

        engine.add(eval1, 0L, false);
        engine.add(truth1, 1L, true);
        Assert.assertEquals(engine.flush(false).size(), 0);
        engine.add(truth2, 2L, true);
        final List<ClosestSVFinder.LinkageConcordanceRecord> out1 = engine.flush(false);
        Assert.assertEquals(engine.flush(true).size(), 0);

        Assert.assertEquals(out1.size(), 1);
        final ClosestSVFinder.LinkageConcordanceRecord outEval1 = out1.get(0);
        Assert.assertEquals(outEval1.record().getId(), eval1.getId());
        Assert.assertEquals(outEval1.id(), 0L);
        assertConcordanceMembers(outEval1.record(), truth1.getId());
        // Would be 1 if missing samples were treated as hom-ref
        Assert.assertEquals(outEval1.record().getAttributes().get(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO), 0.5);
        Assert.assertTrue(outEval1.linkage().getResult());
        Assert.assertEquals(outEval1.linkage().getReciprocalOverlap(), 0.99);
        Assert.assertEquals(outEval1.linkage().getSizeSimilarity(), 1.0);
        Assert.assertEquals(outEval1.linkage().getBreakpointDistance1(), 10);
        Assert.assertEquals(outEval1.linkage().getBreakpointDistance2(), 10);

        final SVCallRecord eval2 = SVTestUtils.makeRecord(
                "eval2",
                "chr2",
                100000,
                true,
                "chr2",
                200000,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList(
                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)),
                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N))
                )
        );

        final SVCallRecord truth3 = SVTestUtils.makeRecord(
                "truth3",
                "chr2",
                200000,
                true,
                "chr2",
                300000,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList(
                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS))
                )
        );

        engine.add(eval2, 3L, false);
        engine.add(truth3, 4L, true);
        final List<ClosestSVFinder.LinkageConcordanceRecord> out2 = engine.flush(true);

        Assert.assertEquals(out2.size(), 1);

        final ClosestSVFinder.LinkageConcordanceRecord outEval2 = out2.get(0);
        Assert.assertEquals(outEval2.record().getId(), eval2.getId());
        Assert.assertEquals(outEval2.id(), 3L);
        assertConcordanceMembers(outEval2.record(), null);

        final Map<String, Object> outEval2Attr = outEval2.record().getAttributes();
        Assert.assertEquals(outEval2Attr.get(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE), ConcordanceState.FALSE_POSITIVE.getAbbreviation());
        Assert.assertEquals(outEval2Attr.get(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO), null);
        Assert.assertEquals(outEval2.linkage(), null);
    }

    @Test
    public void testClustering() {
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        final ClosestSVFinder engine = new ClosestSVFinder(linkage, collapser::annotate, SVTestUtils.hg38Dict);

        final SVCallRecord eval1 = SVTestUtils.makeRecord(
                "eval1",
                "chr1",
                100001,
                true,
                "chr1",
                100002,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        // Clusters with eval1
        final SVCallRecord truth1 = SVTestUtils.makeRecord(
                "truth1",
                "chr1",
                100001,
                true,
                "chr1",
                100002,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        // Clusters with eval1
        final SVCallRecord truth2 = SVTestUtils.makeRecord(
                "truth2",
                "chr1",
                100401,
                true,
                "chr1",
                100401,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        // Does not cluster
        final SVCallRecord truth3 = SVTestUtils.makeRecord(
                "truth3",
                "chr1",
                100601,
                true,
                "chr1",
                100601,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        final SVCallRecord eval2 = SVTestUtils.makeRecord(
                "eval2",
                "chr1",
                200001,
                true,
                "chr1",
                200002,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        // Does not cluster
        final SVCallRecord truth4 = SVTestUtils.makeRecord(
                "truth4",
                "chr1",
                201001,
                true,
                "chr1",
                201002,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        final SVCallRecord eval3 = SVTestUtils.makeRecord(
                "eval3",
                "chr1",
                301101,
                true,
                "chr1",
                301102,
                false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Lists.newArrayList()
        );

        engine.add(eval1, 0L, false);
        engine.add(truth1, 1L, true);
        engine.add(truth2, 2L, true);
        engine.add(truth3, 3L, true);
        engine.add(eval2, 4L, false);
        engine.add(truth4, 5L, true);
        engine.add(eval3, 6L, false);
        final List<SVCallRecord> softFlush1 = engine.flush(false).stream().map(ClosestSVFinder.LinkageConcordanceRecord::record).collect(Collectors.toUnmodifiableList());
        final List<SVCallRecord> softFlush2 = engine.flush(false).stream().map(ClosestSVFinder.LinkageConcordanceRecord::record).collect(Collectors.toUnmodifiableList());
        final List<SVCallRecord> hardFlush = engine.flush(true).stream().map(ClosestSVFinder.LinkageConcordanceRecord::record).collect(Collectors.toUnmodifiableList());

        // Expect eval1 to cluster with truth1 and truth2, and eval2 and eval3 do not cluster
        Assert.assertEquals(softFlush1.size(), 2);
        final SVCallRecord outEval1 = softFlush1.get(0);
        final SVCallRecord outEval2 = softFlush1.get(1);
        Assert.assertEquals(outEval1.getId(), eval1.getId());
        Assert.assertEquals(outEval2.getId(), eval2.getId());
        assertConcordanceMembers(outEval1, truth1.getId());
        assertConcordanceMembers(outEval2, null);

        Assert.assertTrue(softFlush2.isEmpty());

        Assert.assertEquals(hardFlush.size(), 1);
        final SVCallRecord outEval3 = hardFlush.get(0);
        Assert.assertEquals(outEval3.getId(), eval3.getId());
        assertConcordanceMembers(outEval3, null);
    }

    @DataProvider(name = "testCNVMatchesData")
    public Object[][] testCNVMatchesData() {
        final int maxCopyState = 5;
        final Object[][] vals = new Object[(maxCopyState + 1) * (maxCopyState + 1)][];
        int i = 0;
        for (int c1 = 0; c1 <= maxCopyState; c1++) {
            for (int c2 = 0; c2 <= maxCopyState; c2++) {
                final Double conc;
                final boolean expectMatch;
                // Otherwise, matching SV type and then check copy numbers
                expectMatch = true;
                // CNV eval case
                if (c1 == c2) {
                    conc = 1.;
                } else {
                    conc = 0.;
                }
                vals[i++] = new Object[]{c1, c2, expectMatch, conc};
            }
        }
        return vals;
    }

    @Test(dataProvider= "testCNVMatchesData")
    public void testCNVMatches(final int evalCopyNumber, final int truthCopyNumber,
                         final boolean expectMatch, final Double expectedConcordance) {
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        final ClosestSVFinder engine = new ClosestSVFinder(linkage, collapser::annotate, SVTestUtils.hg38Dict);

        final SVCallRecord eval = SVTestUtils.makeRecord(
                "eval",
                "chr1",
                100000,
                null,
                "chr1",
                100001,
                null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null,
                SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                Collections.singletonList(SVTestUtils.getDiploidCNVGenotypeBuilder("sample1", evalCopyNumber))
        );

        final SVCallRecord truth = SVTestUtils.makeRecord(
                "truth",
                "chr1",
                100000,
                null,
                "chr1",
                100001,
                null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null,
                SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                Collections.singletonList(SVTestUtils.getDiploidCNVGenotypeBuilder("sample1", truthCopyNumber))
        );

        engine.add(eval, 0L, false);
        engine.add(truth, 1L, true);

        final List<SVCallRecord> out = engine.flush(true).stream().map(ClosestSVFinder.LinkageConcordanceRecord::record).collect(Collectors.toUnmodifiableList());
        Assert.assertEquals(out.size(), 1);
        final SVCallRecord outEval = out.get(0);
        Assert.assertEquals(outEval.getId(), eval.getId());

        final String state = String.valueOf(outEval.getAttributes().get(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE));
        final Object genotypeConcordanceObj = outEval.getAttributes().get(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO);
        final Double genotypeConcordance = genotypeConcordanceObj == null ? null : Double.valueOf(String.valueOf(genotypeConcordanceObj));
        Assert.assertEquals(genotypeConcordance, expectedConcordance);
        Assert.assertEquals(state, expectMatch ? ConcordanceState.TRUE_POSITIVE.getAbbreviation() :
                ConcordanceState.FALSE_POSITIVE.getAbbreviation());
    }
}