package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class BinnedCNVDefragmenterTest {

    private static final double paddingFraction = 0.5;
    private static final double sampleOverlap = 0.9;
    private static final SVClusterEngine defaultDefragmenter = SVClusterEngineFactory.createCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, sampleOverlap);
    private static final SVClusterEngine binnedDefragmenter = SVClusterEngineFactory.createBinnedCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, 0, SVTestUtils.targetIntervals);

    @Test
    public void testCollapser() {
        final SVCallRecord call1FlattenedDefault = defaultDefragmenter.getCollapser().apply(new SVClusterEngine.OutputCluster(Collections.singletonList(SVTestUtils.call1)));
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call1, call1FlattenedDefault);

        final SVCallRecord call1FlattenedSingleSample = binnedDefragmenter.getCollapser().apply(new SVClusterEngine.OutputCluster(Collections.singletonList(SVTestUtils.call1)));
        SVTestUtils.assertEqualsExceptMembershipAndGT(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecord sameBoundsThreeSamples = binnedDefragmenter.getCollapser().apply(new SVClusterEngine.OutputCluster(Arrays.asList(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch)));
        Assert.assertEquals(sameBoundsThreeSamples.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(sameBoundsThreeSamples.getPositionB(), SVTestUtils.call1.getPositionB());

        final Genotype testGenotype1 = sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample1.make().getSampleName());
        final Genotype expectedGenotype1 = SVTestUtils.sample1.alleles(Lists.newArrayList(Allele.REF_T, Allele.SV_SIMPLE_DEL)).make();
        Assert.assertEquals(testGenotype1.getAlleles(), expectedGenotype1.getAlleles());
        Assert.assertEquals(testGenotype1.getExtendedAttributes(), expectedGenotype1.getExtendedAttributes());

        final Genotype testGenotype2 = sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample2.make().getSampleName());
        final Genotype expectedGenotype2 = SVTestUtils.sample2.alleles(Lists.newArrayList(Allele.REF_T, Allele.SV_SIMPLE_DUP)).make();
        Assert.assertEquals(testGenotype2.getAlleles(), expectedGenotype2.getAlleles());
        Assert.assertEquals(testGenotype2.getExtendedAttributes(), expectedGenotype2.getExtendedAttributes());

        final Genotype testGenotype3 = sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample3.make().getSampleName());
        final Genotype expectedGenotype3 = SVTestUtils.sample3.alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)).make();
        Assert.assertEquals(testGenotype3.getAlleles(), expectedGenotype3.getAlleles());
        Assert.assertEquals(testGenotype3.getExtendedAttributes(), expectedGenotype3.getExtendedAttributes());

        final SVCallRecord overlapping = binnedDefragmenter.getCollapser().apply(new SVClusterEngine.OutputCluster(Arrays.asList(SVTestUtils.call1, SVTestUtils.call2)));
        Assert.assertEquals(overlapping.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(overlapping.getPositionB(), SVTestUtils.call2.getPositionB());
    }

    @DataProvider
    public Object[][] clusterTogetherInputsDefault() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true, "call1 call1"},
                {SVTestUtils.call1, SVTestUtils.call2, true, "call1 call2"},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false, "call1 nonDepthOnly"},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, false, "call1 sameBoundsSampleMismatch"}
        };
    }

    @DataProvider
    public Object[][] clusterTogetherInputsSingleSample() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true, "call1 call1"},
                {SVTestUtils.call1, SVTestUtils.call2, true, "call1 call2"},  //overlapping, same samples
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false, "call1 nonDepthOnly"},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, true, "call1 sameBoundsSampleMismatch"},
                {SVTestUtils.call1_CN1, SVTestUtils.call2_CN0, false, "call1_CN1 call2_CN0"}  //overlapping, but different copy number
        };
    }

    @Test(dataProvider = "clusterTogetherInputsDefault")
    public void testClusterTogetherDefault(final SVCallRecord call1, final SVCallRecord call2,
                                           final boolean expectedResult, final String name) {
        Assert.assertEquals(defaultDefragmenter.getLinkage().areClusterable(call1, call2).getResult(), expectedResult, name);
    }

    @Test(dataProvider = "clusterTogetherInputsSingleSample")
    public void testClusterTogetherSingleSample(final SVCallRecord call1, final SVCallRecord call2,
                                                final boolean expectedResult, final String name) {
        Assert.assertEquals(binnedDefragmenter.getLinkage().areClusterable(call1, call2).getResult(), expectedResult, name);
    }

    @Test
    public void testGetMaxClusterableStartingPosition() {
        Assert.assertEquals(defaultDefragmenter.getLinkage().getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall), SVTestUtils.chr1Length);
        Assert.assertTrue(binnedDefragmenter.getLinkage().getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) == SVTestUtils.chr1Length);  //will be less than chr1length if target intervals are smaller than chr1
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVClusterEngine temp1 = SVClusterEngineFactory.createBinnedCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, 0.8, SVTestUtils.targetIntervals);
        final List<SVCallRecord> output1 = new ArrayList<>();
        output1.addAll(temp1.addAndFlush(SVTestUtils.call1));
        //force new cluster by adding a non-overlapping event
        output1.addAll(temp1.addAndFlush(SVTestUtils.call3));
        output1.addAll(temp1.flush()); //flushes all clusters
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call3, output1.get(1));

        final SVClusterEngine temp2 = SVClusterEngineFactory.createBinnedCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, 0.8, SVTestUtils.targetIntervals);
        final List<SVCallRecord> output2 = new ArrayList<>();
        output2.addAll(temp2.addAndFlush(SVTestUtils.call1));
        output2.addAll(temp2.addAndFlush(SVTestUtils.call2));  //should overlap after padding
        //force new cluster by adding a call on another contig
        output2.addAll(temp2.addAndFlush(SVTestUtils.call4_chr10));
        output2.addAll(temp2.flush());
        Assert.assertEquals(output2.size(), 2);
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.call2.getPositionB());
        SVTestUtils.assertEqualsExceptMembershipAndGT(output2.get(1), SVTestUtils.call4_chr10);

        //cohort case, checking sample set overlap
        final SVClusterEngine temp3 = SVClusterEngineFactory.createCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, CNVLinkage.DEFAULT_PADDING_FRACTION, CNVLinkage.DEFAULT_SAMPLE_OVERLAP);
        final List<SVCallRecord> output3 = new ArrayList<>();
        output3.addAll(temp3.addAndFlush(SVTestUtils.call1));
        output3.addAll(temp3.addAndFlush(SVTestUtils.sameBoundsSampleMismatch));
        output3.addAll(temp3.flush());
        Assert.assertEquals(output3.size(), 2);
    }
}