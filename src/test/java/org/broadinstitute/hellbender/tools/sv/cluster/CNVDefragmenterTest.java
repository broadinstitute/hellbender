package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

public class CNVDefragmenterTest {

    private final SAMSequenceDictionary dictionary = SVTestUtils.hg38Dict;
    private final CNVLinkage defragmenter = new CNVLinkage(dictionary, 0.5, 0.6);

    @Test
    public void testClusterTogether() {
        final SVCallRecord deletion = new SVCallRecord("test_del", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        final SVCallRecord duplication = new SVCallRecord("test_dup", "chr1", 1000, false, "chr1", 1999, true, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertFalse(defragmenter.areClusterable(deletion, duplication).getResult(), "Different sv types should not cluster");

        final SVCallRecord duplicationNonDepthOnly = new SVCallRecord("test_dup", "chr1", 1000, false, "chr1", 1999, true, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, null, Collections.emptyList(),
                1000, Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertFalse(defragmenter.areClusterable(duplication, duplicationNonDepthOnly).getResult(), "Clustered records must be depth-only");

        final SVCallRecord cnv = new SVCallRecord("test_cnv", "chr1", 1000, null, "chr1", 1999, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertFalse(defragmenter.areClusterable(deletion, cnv).getResult(), "Different sv types should not cluster");

        final SVCallRecord insertion = new SVCallRecord("test_ins", "chr1", 1000, true, "chr1", 1001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(),
                1000, Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertFalse(defragmenter.areClusterable(insertion, insertion).getResult(), "Only CNVs should be valid");

        final SVCallRecord deletion2 = new SVCallRecord("test_del2", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertTrue(defragmenter.areClusterable(deletion, deletion2).getResult(), "Valid identical records should cluster");

        final SVCallRecord deletion3 = new SVCallRecord("test_del3", "chr1", 2999, true, "chr1", 3998, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertTrue(defragmenter.areClusterable(deletion, deletion3).getResult(), "Should cluster due to overlap");

        final SVCallRecord deletion4 = new SVCallRecord("test_del3", "chr1", 3000, true, "chr1", 3999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertFalse(defragmenter.areClusterable(deletion, deletion4).getResult(), "Should barely not cluster");
    }

    @DataProvider(name = "testOverlapBuilders")
    public Object[][] testOverlapBuilders() {
        return new Object[][]{
            {
                new GenotypeBuilder()
                        .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                        .alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)),
                    new GenotypeBuilder()
                            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
                            .alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)),
                    new GenotypeBuilder()
                            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0)
                            .alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)),
                    "CN only"
            },
            {
                new GenotypeBuilder()
                        .alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N)),
                    new GenotypeBuilder()
                            .alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)),
                    new GenotypeBuilder()
                            .alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)),
                    "GT only"
            }
        };
    }

    @Test(dataProvider= "testOverlapBuilders")
    public void testSampleOverlap(final GenotypeBuilder refGenotypeMaker,
                                   final GenotypeBuilder altGenotypeMaker1,
                                   final GenotypeBuilder altGenotypeMaker2,
                                   final String name) {
        // alt - ref - ref
        final SVCallRecord deletionARR1 = SVTestUtils.makeRecord("test_del_arr1", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                        new GenotypeBuilder(altGenotypeMaker1.name("sample1").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample2").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample3").make())
                ));
        // identical to above
        final SVCallRecord deletionARR1Copy = SVTestUtils.makeRecord("test_del_arr1_copy", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                    new GenotypeBuilder(altGenotypeMaker1.name("sample1").make()),
                    new GenotypeBuilder(refGenotypeMaker.name("sample2").make()),
                    new GenotypeBuilder(refGenotypeMaker.name("sample3").make())
                ));
        // second alt genotype - single-sample case will not cluster with ARR1
        final SVCallRecord deletionARR2 = SVTestUtils.makeRecord("test_del_arr2", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                        new GenotypeBuilder(altGenotypeMaker2.name("sample1").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample2").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample3").make())
                ));

        // ref - alt - ref
        final SVCallRecord deletionRAR = SVTestUtils.makeRecord("test_del_rar", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                        new GenotypeBuilder(refGenotypeMaker.name("sample1").make()),
                        new GenotypeBuilder(altGenotypeMaker1.name("sample2").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample3").make())
                ));

        // alt - alt - ref
        final SVCallRecord deletionAAR = SVTestUtils.makeRecord("test_del_aar", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                        new GenotypeBuilder(altGenotypeMaker1.name("sample1").make()),
                        new GenotypeBuilder(altGenotypeMaker1.name("sample2").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample3").make())
                ));
        // identical to above but second alt genotype - will cluster because there are multiple samples
        final SVCallRecord deletionAAR2 = SVTestUtils.makeRecord("test_del_aar", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                        new GenotypeBuilder(altGenotypeMaker1.name("sample1").make()),
                        new GenotypeBuilder(altGenotypeMaker1.name("sample2").make()),
                        new GenotypeBuilder(refGenotypeMaker.name("sample3").make())
                ));

        // alt - alt - alt
        final SVCallRecord deletionAAA = SVTestUtils.makeRecord("test_del_aaa", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(
                        new GenotypeBuilder(altGenotypeMaker1.name("sample1").make()),
                        new GenotypeBuilder(altGenotypeMaker1.name("sample2").make()),
                        new GenotypeBuilder(altGenotypeMaker1.name("sample3").make())
                ));
        Assert.assertFalse(defragmenter.areClusterable(deletionARR1, deletionARR2).getResult());
        Assert.assertFalse(defragmenter.areClusterable(deletionARR1, deletionRAR).getResult());
        Assert.assertFalse(defragmenter.areClusterable(deletionARR1, deletionAAR).getResult());
        Assert.assertFalse(defragmenter.areClusterable(deletionARR1, deletionAAA).getResult());
        Assert.assertFalse(defragmenter.areClusterable(deletionRAR, deletionAAR).getResult());
        Assert.assertFalse(defragmenter.areClusterable(deletionRAR, deletionAAA).getResult());

        Assert.assertTrue(defragmenter.areClusterable(deletionARR1, deletionARR1Copy).getResult());
        Assert.assertTrue(defragmenter.areClusterable(deletionAAR, deletionAAA).getResult());
        Assert.assertTrue(defragmenter.areClusterable(deletionAAR, deletionAAR2).getResult());
    }

    @DataProvider(name = "maxPositionIntervals")
    public Object[][] recordPairs() {
        return new Object[][]{
                {100, 200},
                {50000, 500000},
                {1, 1},
                {1, 2},
        };
    }

    @Test(dataProvider= "maxPositionIntervals")
    public void testGetMaxClusterableStartingPosition(final int start, final int end) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start, true, "chr1", end, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                end - start + 1, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        final int maxClusterableStart = defragmenter.getMaxClusterableStartingPosition(call1);

        final int call2Start = maxClusterableStart;
        final int call2End = dictionary.getSequence(call1.getContigA()).getSequenceLength();
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", call2Start, true, "chr1", call2End, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                call2End - call2Start + 1, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertTrue(defragmenter.areClusterable(call1, call2).getResult());

        final int call3Start = maxClusterableStart + 1;
        final int call3End = dictionary.getSequence(call1.getContigA()).getSequenceLength();
        final SVCallRecord call3 = new SVCallRecord("call3", "chr1", call3Start, true, "chr1", call3End, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                call3End - call3Start + 1, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, dictionary);
        Assert.assertFalse(defragmenter.areClusterable(call1, call3).getResult());
    }
}