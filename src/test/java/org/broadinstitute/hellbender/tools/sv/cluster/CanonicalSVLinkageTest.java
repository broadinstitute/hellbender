package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CanonicalSVLinkageTest extends GATKBaseTest {

    private static final CanonicalSVLinkage<SVCallRecord> linkage = SVTestUtils.getNewDefaultLinkage();
    private static final ClusteringParameters depthOnlyParametersSizeSimilarity = ClusteringParameters.createDepthParameters(0.1, 0.5, 10000000, 0);
    private static final ClusteringParameters mixedParametersSizeSimilarity = ClusteringParameters.createMixedParameters(0.1, 0.5, 5000, 0);
    private static final ClusteringParameters evidenceParametersSizeSimilarity = ClusteringParameters.createPesrParameters(0.1, 0.5, 5000, 0);
    private static final CanonicalSVLinkage<SVCallRecord> linkageSizeSimilarity = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, false);

    // Assign length according to coordinate scheme
    private static Integer inferLength(final String contigA, final int posA, final String contigB, final int posB) {
        if (contigA.equals(contigB)) {
            if (posA == posB) {
                // no interval
                return null;
            } else {
                // intrachromosomal and intervaled
                return posB - posA;
            }
        } else {
            // interchromosomal
            return null;
        }
    }

    @BeforeTest
    public void initializeClusterEngine() {
        linkageSizeSimilarity.setDepthOnlyParams(depthOnlyParametersSizeSimilarity);
        linkageSizeSimilarity.setMixedParams(mixedParametersSizeSimilarity);
        linkageSizeSimilarity.setEvidenceParams(evidenceParametersSizeSimilarity);
    }

    @Test
    public void testClusterTogether() {
        Assert.assertTrue(linkage.areClusterable(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff).getResult());
        Assert.assertFalse(linkage.areClusterable(SVTestUtils.depthOnly, SVTestUtils.inversion).getResult());
        Assert.assertFalse(linkage.areClusterable(SVTestUtils.call1, SVTestUtils.call2).getResult());
        Assert.assertTrue(linkage.areClusterable(SVTestUtils.call1, SVTestUtils.overlapsCall1).getResult());
    }

    @DataProvider(name = "clusterTogetherVaryPositionsProvider")
    public Object[][] clusterTogetherVaryPositionsProvider() {
        return new Object[][]{
                {500, 1001, 1001, 1502, false, 1/502., 1., 501, 501},  // abutting
                {500, 1001, 500, 1001, true, 1., 1., 0, 0},  // exactly equal
                {500, 1001, 600, 1101, true, 402/502., 1., 100, 100},  // barely meet reciprocal overlap
                {500, 1001, 601, 1102, false, 401/502., 1., 101, 101}, // call2 shifted slightly up
                {500, 1000, 600, 1101, false, 401/502., 501/502., 100, 101}, // call1 slightly larger
                {500, 501, 500, 501, true, 1., 1., 0, 0}, // tiny but equal
                {500, 501, 500, 502, false, 2/3., 2/3., 0, 1}, // tiny but call2 twice as big
                {500, 500, 500, 500, true, 1., 1., 0, 0}, // 0-length and equal
                {500, 500, 501, 501, false, 0., 1., 1, 1}, // 0-length and not equal
                {1, SVTestUtils.chr1Length, 1, SVTestUtils.chr1Length, true, 1., 1., 0, 0}, // really big
                {1, 10001, 1, 10001, true, 1., 1., 0, 0}, // left contig edge
                {SVTestUtils.chr1Length - 10000, SVTestUtils.chr1Length, SVTestUtils.chr1Length - 10000, SVTestUtils.chr1Length, true, 1., 1., 0, 0}, // right contig edge
                {100000, 200000, 101001, 201001, false, 99000/100001., 1., 1001, 1001}, // window test fail
                {100000, 200000, 101000, 201000, true, 99001/100001., 1., 1000, 1000} // window test success
        };
    }

    @Test(dataProvider= "clusterTogetherVaryPositionsProvider")
    public void testClusterTogetherVaryPositions(final int start1, final int end1, final int start2, final int end2,
                                                 final boolean result, final double reciprocalOverlap, final double sizeSimilarity, final int dist1, final int dist2) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start1, true,
                "chr1", end1, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), end1 - start1 + 1, Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                SVTestUtils.threeGenotypes, Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", start2, true,
                "chr1", end2, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), end2 - start2 + 1, Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                SVTestUtils.threeGenotypes, Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final CanonicalSVLinkage.CanonicalLinkageResult linkageResult = linkage.areClusterable(call1, call2);
        Assert.assertEquals(linkageResult.getResult(), result);
        Assert.assertEquals(linkageResult.getReciprocalOverlap(), reciprocalOverlap);
        Assert.assertEquals(linkageResult.getSizeSimilarity(), sizeSimilarity);
        Assert.assertEquals(linkageResult.getBreakpointDistance1(), dist1);
        Assert.assertEquals(linkageResult.getBreakpointDistance2(), dist2);
    }

    @DataProvider(name = "testClusterTogetherWithSizeSimilarityDataProvider")
    public Object[][] testClusterTogetherWithSizeSimilarityDataProvider() {
        return new Object[][]{
                {1001, 2000, 1001, 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, true},
                {1001, 2000, 1501, 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, true},
                {1001, 2000, 1900, 2900, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, true},
                // Fails size similarity
                {1001, 2000, 1502, 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, false},
                {1001, 2000, 1901, 2399, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, false},
                // Fails reciprocal overlap
                {1001, 2000, 1902, 2900, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, false}
        };
    }

    @Test(dataProvider= "testClusterTogetherWithSizeSimilarityDataProvider")
    public void testClusterTogetherWithSizeSimilarity(final int start1, final int end1,
                                                      final int start2, final int end2,
                                                      final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                                      final boolean expectedResult) {
        // PESR
        final SVCallRecord record1 = SVTestUtils.newPESRCallRecordWithIntervalAndType(start1, end1, svtype);
        final SVCallRecord record2 = SVTestUtils.newPESRCallRecordWithIntervalAndType(start2, end2, svtype);
        Assert.assertEquals(linkageSizeSimilarity.areClusterable(record1, record2).getResult(), expectedResult);
        Assert.assertEquals(linkageSizeSimilarity.areClusterable(record2, record1).getResult(), expectedResult);

        // Depth only
        final SVCallRecord record3 = SVTestUtils.newDepthCallRecordWithIntervalAndType(start1, end1, svtype);
        final SVCallRecord record4 = SVTestUtils.newDepthCallRecordWithIntervalAndType(start2, end2, svtype);
        Assert.assertEquals(linkageSizeSimilarity.areClusterable(record3, record4).getResult(), expectedResult);
        Assert.assertEquals(linkageSizeSimilarity.areClusterable(record4, record3).getResult(), expectedResult);
    }

    @DataProvider(name = "testClusterTogetherWithSizeSimilarityInsertionsDataProvider")
    public Object[][] testClusterTogetherWithSizeSimilarityInsertionsDataProvider() {
        return new Object[][]{
                {1000, 1000, 1000, 1000, true},
                {1000, 1000, 1000, 500, true},
                {1000, 1000, 1000, 2000, true},
                {1000, 1000, 100, 1000, true},
                {1000, 1000, 1900, 1000, true},
                // Fails reciprocal overlap
                {1000, 1000, 99, 1000, false},
                {1000, 1000, 1901, 1000, false},
                // Fails size similarity
                {1000, 1000, 1000, 499, false},
                {1000, 2001, 1000, 1000, false},
        };
    }

    @Test(dataProvider= "testClusterTogetherWithSizeSimilarityInsertionsDataProvider")
    public void testClusterTogetherWithSizeSimilarityInsertions(final int start1, final int length1,
                                                                final int start2, final int length2,
                                                                final boolean expectedResult) {
        final SVCallRecord record1 = SVTestUtils.newInsertionWithPositionAndLength(start1, length1);
        final SVCallRecord record2 = SVTestUtils.newInsertionWithPositionAndLength(start2, length2);
        Assert.assertEquals(linkageSizeSimilarity.areClusterable(record1, record2).getResult(), expectedResult);
        Assert.assertEquals(linkageSizeSimilarity.areClusterable(record2, record1).getResult(), expectedResult);
    }

    @Test(expectedExceptions = { IllegalArgumentException.class })
    public void testClusterTogetherInvalidInterval() {
        // End position beyond contig end after padding
        final SVCallRecord deletion1 = new SVCallRecord("test_del", "chr1", 1000, true, "chr1", 248956423 + SVTestUtils.defaultEvidenceParameters.getWindow(), false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                null, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord deletion2 = new SVCallRecord("test_del", "chr1", 1000, true, "chr1", 248956422, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                null, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        linkage.areClusterable(deletion1, deletion2);
        Assert.fail("Expected exception not thrown");
    }

    @Test
    public void testGetClusteringIntervalEdge() {
        //edge case - end of contig
        Assert.assertTrue(linkage.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) <= SVTestUtils.chr1Length);
    }

    @DataProvider(name = "maxPositionIntervals")
    public Object[][] recordPairs() {
        return new Object[][]{
                {100, 200},
                {50000, 500000},
                {1, 1},
                {1, 2}
        };
    }

    @Test(dataProvider= "maxPositionIntervals")
    public void testGetMaxClusterableStartingPosition(final int start, final int end) {
        testGetMaxClusterableStartingPositionWithAlgorithm(start, end, GATKSVVCFConstants.DEPTH_ALGORITHM);
        testGetMaxClusterableStartingPositionWithAlgorithm(start, end, SVTestUtils.PESR_ALGORITHM);
    }

    private void testGetMaxClusterableStartingPositionWithAlgorithm(final int start, final int end, final String algorithm) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start, true, "chr1", end, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                end - start + 1, Collections.emptyList(), Collections.singletonList(algorithm),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final int maxClusterableStart = linkage.getMaxClusterableStartingPosition(call1);

        final int call2Start = maxClusterableStart;
        final SVCallRecord call2Depth = new SVCallRecord("call2", "chr1", call2Start, true, "chr1", call2Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2Pesr = new SVCallRecord("call2", "chr1", call2Start, true, "chr1", call2Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertTrue(linkage.areClusterable(call1, call2Depth).getResult() || linkage.areClusterable(call1, call2Pesr).getResult());

        final int call3Start = maxClusterableStart + 1;
        final SVCallRecord call3Depth = new SVCallRecord("call2", "chr1", call3Start, true, "chr1", call3Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call3Pesr = new SVCallRecord("call2", "chr1", call3Start, true, "chr1", call3Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertFalse(linkage.areClusterable(call1, call3Depth).getResult() || linkage.areClusterable(call1, call3Pesr).getResult());
    }

    @Test
    public void testGetClusteringIntervalLists() {
        //test overloaded function with List
        final List<SVCallRecord> pesrClusterList = new ArrayList<>();
        pesrClusterList.add(SVTestUtils.depthAndStuff);
        final int pesrCluster1 = linkage.getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + SVTestUtils.defaultEvidenceParameters.getWindow());
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + SVTestUtils.defaultMixedParameters.getWindow());
        //add an upstream variant
        pesrClusterList.add(SVTestUtils.depthAndStuff2);
        final int pesrCluster2 = linkage.getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster2 >= pesrCluster1);
        //add a downstream variant
        pesrClusterList.add(SVTestUtils.depthAndStuff3);
        final int pesrCluster3 = linkage.getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + SVTestUtils.defaultEvidenceParameters.getWindow());
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + SVTestUtils.defaultMixedParameters.getWindow());
    }

    @Test
    public void testIsDepthOnlyCall() {
        Assert.assertTrue(SVTestUtils.call1.isDepthOnly());
        Assert.assertFalse(SVTestUtils.depthAndStuff.isDepthOnly());
        Assert.assertFalse(SVTestUtils.inversion.isDepthOnly());
    }

    @DataProvider(name = "testClusterTogetherVaryPositionsBNDData")
    public Object[][] testClusterTogetherVaryPositionsBNDData() {
        return new Object[][]{
                // length 0 edge case
                {
                        "chr1", 100, "chr1", 100,
                        "chr1", 100, "chr1", 100,
                        0.5, 0.5, 100,
                        true, 1., 1., 0, 0
                },
                // simple interchromosomal
                {
                        "chr1", 100, "chr2", 100,
                        "chr1", 100, "chr2", 100,
                        0.5, 0.5, 0,
                        true, null, null, 0, 0
                },
                {
                        "chr1", 100, "chr2", 101,
                        "chr1", 100, "chr2", 100,
                        0.5, 0.5, 0,
                        false, null, null, 0, 1
                },
                {
                        "chr1", 101, "chr2", 100,
                        "chr1", 100, "chr2", 100,
                        0.5, 0.5, 0,
                        false, null, null, 1, 0
                },
                // reciprocal overlap for intrachromosomal records
                {
                        "chr1", 100, "chr1", 200,
                        "chr1", 150, "chr1", 250,
                        0.5, 0.5, 50,
                        true, 0.5, 1., 50, 50
                },
                // fails window
                {
                        "chr1", 100, "chr1", 200,
                        "chr1", 150, "chr1", 250,
                        0.5, 0.5, 49,
                        false, 0.5, 1., 50, 50
                },
                {
                        "chr1", 100, "chr1", 200,
                        "chr1", 100, "chr1", 250,
                        0.5, 0.5, 49,
                        false, 2/3., 2/3., 0, 50
                },
                {
                        "chr1", 100, "chr1", 200,
                        "chr1", 150, "chr1", 200,
                        0.5, 0.5, 49,
                        false, 0.5, 0.5, 50, 0
                },
                // fails reciprocal overlap
                {
                        "chr1", 100, "chr1", 200,
                        "chr1", 150, "chr1", 250,
                        0.51, 0.5, 50,
                        false, 0.5, 1., 50, 50
                },
                // fails size similarity
                {
                        "chr1", 100, "chr1", 200,
                        "chr1", 150, "chr1", 200,
                        0.5, 0.51, 50,
                        false, 0.5, 0.5, 50, 0
                }
        };
    }

    @Test(dataProvider= "testClusterTogetherVaryPositionsBNDData")
    public void testClusterTogetherVaryPositionsBND(final String chrom1A, final int start1, final String chrom1B, final int end1,
                                                    final String chrom2A, final int start2, final String chrom2B, final int end2,
                                                    final double reciprocalOverlapThresh, final double sizeSimilarityThresh, final int window,
                                                    final boolean result, final Double expectedReciprocalOverlap, final Double expectedSizeSimilarity,
                                                    final Integer expectedDist1, final Integer expectedDist2) {
        final SVCallRecord call1 = new SVCallRecord("call1", chrom1A, start1, true,
                chrom1B, end1, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null, Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, SVTestUtils.BND_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call2", chrom2A, start2, true,
                chrom2B, end2, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null, Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, SVTestUtils.BND_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final CanonicalSVLinkage<SVCallRecord> linkage = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, false);
        linkage.setEvidenceParams(ClusteringParameters.createPesrParameters(reciprocalOverlapThresh, sizeSimilarityThresh, window, 0));
        final CanonicalSVLinkage.CanonicalLinkageResult linkageResult = linkage.areClusterable(call1, call2);
        Assert.assertEquals(linkageResult.getResult(), result);
        Assert.assertEquals(linkageResult.getReciprocalOverlap(), expectedReciprocalOverlap);
        Assert.assertEquals(linkageResult.getSizeSimilarity(), expectedSizeSimilarity);
        Assert.assertEquals(linkageResult.getBreakpointDistance1(), expectedDist1);
        Assert.assertEquals(linkageResult.getBreakpointDistance2(), expectedDist2);
    }

    @Test
    public void testClusterTogetherVaryTypes() {
        for (final GATKSVVCFConstants.StructuralVariantAnnotationType type1 : GATKSVVCFConstants.StructuralVariantAnnotationType.values()) {
            // Pass in null strands to let them be determined automatically
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, SVTestUtils.getValidTestStrandA(type1),
                    "chr1", 2001, SVTestUtils.getValidTestStrandB(type1), type1, null, Collections.emptyList(),
                    SVTestUtils.getLength(1000, 2001, type1), Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                    Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
            for (final GATKSVVCFConstants.StructuralVariantAnnotationType type2 : GATKSVVCFConstants.StructuralVariantAnnotationType.values()) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, SVTestUtils.getValidTestStrandA(type2),
                        "chr1", 2001, SVTestUtils.getValidTestStrandB(type2), type2, null, Collections.emptyList(),
                        SVTestUtils.getLength(1000, 2001, type2), Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                // Should only cluster together if same type, except CNVs
                if ((type1 == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV && call2.isSimpleCNV()) ||
                        (type2 == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV && call1.isSimpleCNV())) {
                    Assert.assertTrue(linkage.areClusterable(call1, call2).getResult());
                } else {
                    Assert.assertEquals(linkage.areClusterable(call1, call2).getResult(), type1 == type2);
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryStrands() {
        final List<Boolean> bools = Lists.newArrayList(Boolean.TRUE, Boolean.FALSE);
        for (final Boolean strand1A : bools) {
            for (final Boolean strand1B : bools) {
                final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, strand1A,
                        "chr1", 2001, strand1B, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                        null, Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                for (final Boolean strand2A : bools) {
                    for (final Boolean strand2B : bools) {
                        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, strand2A,
                                "chr1", 2001, strand2B, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                                null, Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                        // Should only cluster if strands match
                        Assert.assertEquals(linkage.areClusterable(call1, call2).getResult(), strand1A == strand2A && strand1B == strand2B);
                    }
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryContigs() {
        final List<String> contigs = Lists.newArrayList("chr1", "chrX");
        for (int i = 0; i < contigs.size(); i++) {
            final String contig1A = contigs.get(i);
            for (int j = i; j < contigs.size(); j++) {
                final String contig1B = contigs.get(j);
                final SVCallRecord call1 = new SVCallRecord("call1", contig1A, 1000, true,
                        contig1B, 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                        null, Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                for (int k = 0; k < contigs.size(); k++) {
                    final String contig2A = contigs.get(k);
                    for (int m = k; m < contigs.size(); m++) {
                        final String contig2B = contigs.get(m);
                        final SVCallRecord call2 = new SVCallRecord("call2", contig2A, 1000, true,
                                contig2B, 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                                null, Collections.emptyList(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                        // Should only cluster if contigs match
                        Assert.assertEquals(linkage.areClusterable(call1, call2).getResult(), contig1A.equals(contig2A) && contig1B.equals(contig2B));
                    }
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryAlgorithms() {
        final List<List<String>> algorithmsList = Lists.newArrayList(
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM),
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST
        );
        for (final List<String> algorithms1 : algorithmsList) {
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                    "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                    1002, Collections.emptyList(), algorithms1, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
            for (final List<String> algorithms2 : algorithmsList) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, true,
                        "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                        1002, Collections.emptyList(), algorithms2, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                // All combinations should cluster
                Assert.assertTrue(linkage.areClusterable(call1, call2).getResult());
            }
        }
    }

    @Test
    public void testClusterTogetherCNVs() {
        final SVCallRecord cnv1 = SVTestUtils.makeRecord("cnv1", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                Collections.emptyList());

        final SVCallRecord cnv2 = SVTestUtils.makeRecord("cnv2", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                Collections.emptyList());

        final SVCallRecord del1 = SVTestUtils.makeRecord("del1", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList());

        final SVCallRecord dup1 = SVTestUtils.makeRecord("dup1", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList());

        Assert.assertTrue(linkage.areClusterable(cnv1, cnv2).getResult());
        Assert.assertTrue(linkage.areClusterable(cnv1, del1).getResult());
        Assert.assertTrue(linkage.areClusterable(cnv1, dup1).getResult());
        Assert.assertFalse(linkage.areClusterable(del1, dup1).getResult());
    }

    @DataProvider(name = "testMatchCNVNoGTData")
    public Object[][] testMatchCNVNoGTData() {
        return new Object[][]{
                // Empty
                {0, new int[]{}, new int[]{}, true},
                // Both equal
                {0, new int[]{0}, new int[]{0}, true},
                {1, new int[]{1}, new int[]{1}, true},
                {2, new int[]{2}, new int[]{2}, true},
                {2, new int[]{3}, new int[]{3}, true},
                // Unequal
                {2, new int[]{1}, new int[]{2}, false},
                {2, new int[]{2}, new int[]{1}, false},
                // Equal multiple
                {2, new int[]{2, 2}, new int[]{2, 2}, true},
                {2, new int[]{4, 2}, new int[]{4, 2}, true},
                // Unequal multiple
                {2, new int[]{2, 2}, new int[]{2, 1}, false},
                {2, new int[]{0, 2}, new int[]{1, 1}, false},
                {2, new int[]{3, 2}, new int[]{2, 2}, false},
                {2, new int[]{6, 2}, new int[]{4, 2}, false},
        };
    }

    @Test(dataProvider= "testMatchCNVNoGTData")
    public void testMatchCNVNoGT(final int ploidy, final int[] copyNumbers1, final int[] copyNumbers2, final boolean expected) {
        final List<Allele> alleles = Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_CNV);
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
        // Create genotypes with copy number attribute (and no GT)
        final SVCallRecord recordCN1 = SVTestUtils.getCNVRecordWithCN(ploidy, alleles, svtype, copyNumbers1, GATKSVVCFConstants.COPY_NUMBER_FORMAT);
        final SVCallRecord recordCN2 = SVTestUtils.getCNVRecordWithCN(ploidy, alleles, svtype, copyNumbers2, GATKSVVCFConstants.COPY_NUMBER_FORMAT);

        // With sample overlap
        final ClusteringParameters depthOnlyParams = ClusteringParameters.createDepthParameters(0.8, 0, 10000000, 1);
        final CanonicalSVLinkage<SVCallRecord> linkage = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, false);
        linkage.setDepthOnlyParams(depthOnlyParams);

        Assert.assertEquals(linkage.areClusterable(recordCN1, recordCN2).getResult(), expected);

        final SVCallRecord recordRDCN1 = SVTestUtils.getCNVRecordWithCN(ploidy, alleles, svtype, copyNumbers1, GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT);
        final SVCallRecord recordRDCN2 = SVTestUtils.getCNVRecordWithCN(ploidy, alleles, svtype, copyNumbers2, GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT);
        Assert.assertEquals(linkage.areClusterable(recordRDCN1, recordRDCN2).getResult(), expected);
    }

    @DataProvider(name = "testClusterTogetherIntervaledComplexData")
    public Object[][] testClusterTogetherIntervaledComplexData() {
        return new Object[][]{
                // exact match
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        true
                },
                // match within parameters
                {"chr1", 1100, "chr1", 1900,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1150, 1550)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1650, 1800))
                        ),
                        true
                },
                // different contigs
                {"chr2", 1000, "chr2", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // different coordinates
                {"chr1", 1600, "chr1", 2400,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // different subtypes
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINVdel,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // different number of intervals
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500))
                        ),
                        false
                },
                // different cpx intervals
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 800, 1100)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // second cpx interval type is different
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
        };
    }

    @Test(dataProvider= "testClusterTogetherIntervaledComplexData")
    public void testClusterTogetherIntervaledComplex(final String contigA, final int posA, final String contigB, final int posB,
                                                     final GATKSVVCFConstants.ComplexVariantSubtype subtype,
                                                     final List<SVCallRecord.ComplexEventInterval> cpxIntervals, final boolean expected) {
        final SVCallRecord cpx1 = new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 2000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                Arrays.asList(SVCallRecord.ComplexEventInterval.decode("DEL_chr1:1100-1500", SVTestUtils.hg38Dict), SVCallRecord.ComplexEventInterval.decode("INV_chr1:1600-1900", SVTestUtils.hg38Dict)),
                1000, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Integer length2 = inferLength(contigA, posA, contigB, posB);
        final SVCallRecord cpx2 = new SVCallRecord("cpx2", contigA, posA, null,
                contigB, posB, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                subtype,
                cpxIntervals,
                length2, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertEquals(linkage.areClusterable(cpx1, cpx2).getResult(), expected);
    }

    @DataProvider(name = "testClusterTogetherInsertedComplexData")
    public Object[][] testClusterTogetherInsertedComplexData() {
        return new Object[][]{
                // exact match
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6000, 8000))
                        ),
                        true
                },
                // close match
                {"chr1", 1010, "chr1", 1010,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6200, 7990))
                        ),
                        true
                },
                // not match by coordinates
                {"chr1", 2000, "chr1", 3000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6200, 7990))
                        ),
                        false
                },
                // not match by subtype
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6000, 8000))
                        ),
                        false
                },
                // not match by cpx interval
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 9000, 11000))
                        ),
                        false
                },
                // different cpx interval type
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 6000, 8000))
                        ),
                        false
                },
        };
    }

    @Test(dataProvider= "testClusterTogetherInsertedComplexData")
    public void testClusterTogetherInsertedComplex(final String contigA, final int posA, final String contigB, final int posB,
                                                   final GATKSVVCFConstants.ComplexVariantSubtype subtype,
                                                   final List<SVCallRecord.ComplexEventInterval> cpxIntervals, final boolean expected) {
        final SVCallRecord cpx1 = new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 1000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                Arrays.asList(new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6000, 8000))),
                2000, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Integer length2 = cpxIntervals.get(0).getInterval().size();
        final SVCallRecord cpx2 = new SVCallRecord("cpx2", contigA, posA, null,
                contigB, posB, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                subtype,
                cpxIntervals,
                length2, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertEquals(linkage.areClusterable(cpx1, cpx2).getResult(), expected);
    }

    @Test
    public void testClusterTogetherVaryParameters() {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1002, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1100, true,
                "chr1", 2101, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1002, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        // Cluster with default parameters
        Assert.assertTrue(linkage.areClusterable(call1, call2).getResult());
        final ClusteringParameters exactMatchParameters = ClusteringParameters.createDepthParameters(1.0, 0, 0, 1.0);
        final CanonicalSVLinkage<SVCallRecord> exactMatchLinkage = SVTestUtils.getNewDefaultLinkage();
        exactMatchLinkage.setDepthOnlyParams(exactMatchParameters);
        // Do not cluster requiring exact overlap
        Assert.assertFalse(exactMatchLinkage.areClusterable(call1, call2).getResult());
    }
}