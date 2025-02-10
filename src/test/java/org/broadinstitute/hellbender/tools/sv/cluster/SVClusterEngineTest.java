package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
import static org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE;

public class SVClusterEngineTest {

    private final SVClusterEngine engine = SVTestUtils.defaultSingleLinkageEngine;

    private static final String SAMPLE_1_NAME = "sample1";
    private static final String SAMPLE_2_NAME = "sample2";
    private static final String SAMPLE_3_NAME = "sample3";
    private static final String SAMPLE_4_NAME = "sample4";

    @BeforeTest
    public void initializeClusterEngine() {
        engine.addAndFlush(SVTestUtils.call1);
    }

    @Test
    public void testLinkage() {
        Assert.assertTrue(engine.getLinkage().areClusterable(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff).getResult());
        Assert.assertFalse(engine.getLinkage().areClusterable(SVTestUtils.depthOnly, SVTestUtils.inversion).getResult());
        Assert.assertFalse(engine.getLinkage().areClusterable(SVTestUtils.call1, SVTestUtils.call2).getResult());
        Assert.assertTrue(engine.getLinkage().areClusterable(SVTestUtils.call1, SVTestUtils.overlapsCall1).getResult());
    }

    @Test
    public void testCollapser() {
        //depth only and depthAndStuff have same bounds, less than call2
        final SVClusterEngine.OutputCluster testCluster = new SVClusterEngine.OutputCluster(Arrays.asList(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff, SVTestUtils.call2));
        final SVCallRecord flattened = engine.getCollapser().apply(testCluster);
        Assert.assertEquals(flattened.getPositionA(), SVTestUtils.depthAndStuff.getPositionA());
        Assert.assertEquals(flattened.getPositionB(), SVTestUtils.depthAndStuff.getPositionB());
        //should have all the algs
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthAndStuff.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthOnly.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.call2.getAlgorithms()));

        SVTestUtils.assertContainsAllIgnoreRefAlleleBase(flattened.getGenotypes(), SVTestUtils.depthAndStuff.getGenotypes(), true);
        SVTestUtils.assertContainsAllIgnoreRefAlleleBase(flattened.getGenotypes(), SVTestUtils.depthOnly.getGenotypes(), true);
        SVTestUtils.assertContainsAllIgnoreRefAlleleBase(flattened.getGenotypes(), SVTestUtils.call2.getGenotypes(), true);

        // Test subtyped insertions
        final List<SVCallRecord> subtypedRecords = Lists.newArrayList(
                SVTestUtils.newCallRecordWithAlleles(
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        2,
                        null),
                SVTestUtils.newCallRecordWithAlleles(
                        Lists.newArrayList(Allele.REF_N, Allele.create("<INS:MEI>")),
                        Lists.newArrayList(Allele.REF_N, Allele.create("<INS:MEI>")),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        2,
                        null)
        );
        final SVCallRecord subtypedFlattened = engine.getCollapser().apply(new SVClusterEngine.OutputCluster(subtypedRecords));
        Assert.assertEquals(subtypedFlattened.getAlleles().size(), 2);
        Assert.assertEquals(subtypedFlattened.getAltAlleles().size(), 1);
        final List<Allele> collapsedAlleles = subtypedFlattened.getGenotypes().stream().map(Genotype::getAlleles)
                .flatMap(Collection::stream).distinct().sorted().collect(Collectors.toList());
        Assert.assertEquals(subtypedFlattened.getAlleles(), collapsedAlleles);
    }


    @DataProvider(name = "testAddVaryPositionsProvider")
    public Object[][] testAddVaryPositionsProvider() {
        return new Object[][]{
                {10000, 20001, 12000, 22001, 14000, 24001, SINGLE_LINKAGE, 1},
                {10000, 20001, 12000, 22001, 14001, 24002, SINGLE_LINKAGE, 2},
                {10000, 20001, 12001, 22002, 14002, 24003, SINGLE_LINKAGE, 3},
                {10000, 20001, 11000, 21001, 12000, 22001, MAX_CLIQUE, 1},
                {10000, 20001, 12000, 22001, 14000, 24001, MAX_CLIQUE, 2},
                {10000, 20001, 12001, 22002, 14002, 24003, MAX_CLIQUE, 3}
        };
    }

    @Test(dataProvider= "testAddVaryPositionsProvider")
    public void testAddVaryPositions(final int positionA1, final int positionB1,
                                     final int positionA2, final int positionB2,
                                     final int positionA3, final int positionB3,
                                     final SVClusterEngine.CLUSTERING_TYPE type,
                                     final int result) {
        final SVClusterEngine engine;
        if (type == SINGLE_LINKAGE) {
            engine = SVTestUtils.getNewDefaultSingleLinkageEngine();
        } else if (type == MAX_CLIQUE) {
            engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        } else {
            throw new TestException("Unimplemented clustering type " + type.name());
        }
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", positionA1, true,
                "chr1", positionB1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                positionB1 - positionA1 + 1, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call1", "chr1", positionA2, true,
                "chr1", positionB2, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                positionB2 - positionA2 + 1, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call3 = new SVCallRecord("call1", "chr1", positionA3, true,
                "chr1", positionB3, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                positionB3 - positionA3 + 1, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final List<SVCallRecord> output = new ArrayList<>();
        output.addAll(engine.addAndFlush(call1));
        output.addAll(engine.addAndFlush(call2));
        output.addAll(engine.addAndFlush(call3));
        output.addAll(engine.flush());
        Assert.assertEquals(output.size(), result);
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVClusterEngine temp1 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        Assert.assertTrue(temp1.isEmpty());
        final List<SVCallRecord> output1 = new ArrayList<>();
        output1.addAll(temp1.addAndFlush(SVTestUtils.call1));
        Assert.assertFalse(temp1.isEmpty());
        //force new cluster by adding a non-overlapping event
        output1.addAll(temp1.addAndFlush(SVTestUtils.call3));
        output1.addAll(temp1.flush()); //flushes all clusters
        Assert.assertTrue(temp1.isEmpty());
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call3, output1.get(1));

        final SVClusterEngine temp2 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        final List<SVCallRecord> output2 = new ArrayList<>();
        output2.addAll(temp2.addAndFlush(SVTestUtils.call1));
        output2.addAll(temp2.addAndFlush(SVTestUtils.overlapsCall1));
        //force new cluster by adding a call on another contig
        output2.addAll(temp2.addAndFlush(SVTestUtils.call4_chr10));
        output2.addAll(temp2.flush());
        Assert.assertEquals(output2.size(), 2);
        //median of two items ends up being the second item here
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.call1.getPositionB());
        SVTestUtils.assertEqualsExceptMembershipAndGT(output2.get(1), SVTestUtils.call4_chr10);

        //checking insensitivity to sample set overlap
        final SVClusterEngine temp3 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        final List<SVCallRecord> output3 = new ArrayList<>();
        output3.addAll(temp3.addAndFlush(SVTestUtils.call1));
        output3.addAll(temp3.addAndFlush(SVTestUtils.sameBoundsSampleMismatch));
        output3.addAll(temp3.flush());
        Assert.assertEquals(output3.size(), 1);
        Assert.assertEquals(output3.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output3.get(0).getPositionB(), SVTestUtils.call1.getPositionB());
        Assert.assertEquals(output3.get(0).getPositionA(), SVTestUtils.sameBoundsSampleMismatch.getPositionA());
        Assert.assertEquals(output3.get(0).getPositionB(), SVTestUtils.sameBoundsSampleMismatch.getPositionB());
    }

    @Test
    public void testAddMaxCliqueLarge() {
        final int numRecords = 100;
        final SVClusterEngine engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        final int length = 5000;
        final List<SVCallRecord> result = new ArrayList<>();
        for (int i = 0; i < numRecords; i++) {
            final int start = 1000 + 10 * i;
            final int end = start + length - 1;
            result.addAll(engine.addAndFlush(SVTestUtils.newPESRCallRecordWithIntervalAndType(start, end, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        }
        result.addAll(engine.flush());
        Assert.assertEquals(result.size(), 50);
        for (final SVCallRecord resultRecord : result) {
            Assert.assertTrue(resultRecord.getAttributes().containsKey(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            final int numMembers = VariantContextGetters.attributeToList(resultRecord.getAttributes().get(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY)).size();
            Assert.assertEquals(numMembers, 51);
        }
    }

    @DataProvider(name = "testGetCarrierSamplesBiallelicData")
    public Object[][] testGetCarrierSamplesBiallelicData() {
        return new Object[][]{
                // Empty
                {0, null, null, new int[]{}, new String[]{}},
                {0, null, Allele.SV_SIMPLE_DEL, new int[]{}, new String[]{}},
                {0, Allele.REF_N, null, new int[]{}, new String[]{}},
                {0, Allele.NO_CALL, null, new int[]{}, new String[]{}},
                {0, Allele.NO_CALL, Allele.SV_SIMPLE_DEL, new int[]{}, new String[]{}},
                {0, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{}, new String[]{}},
                // Ploidy 0
                {0, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{0}, new String[]{}},
                // Ploidy 1, no carrier
                {1, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{1, 1}, new String[]{}},
                // Ploidy 1, with carrier
                {1, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{1, 0}, new String[]{"1"}},
                // Ploidy 2, no carrier
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 2, 2}, new String[]{}},
                // Ploidy 2, with carriers
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 1, 1}, new String[]{"1", "2"}},
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 1, 0}, new String[]{"1", "2"}},
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 1, 3}, new String[]{"1"}}, // third is not consistent with DEL
                // DUP
                {2, Allele.REF_N, Allele.SV_SIMPLE_DUP, new int[]{2, 3, 2}, new String[]{"1"}},
                {2, Allele.REF_N, Allele.SV_SIMPLE_DUP, new int[]{2, 1, 3}, new String[]{"2"}}, // second is not consistent with DUP
                // Ploidy 3
                {3, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{3, 3, 3}, new String[]{}},
                {3, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{0, 1, 2, 3}, new String[]{"0", "1", "2"}}
        };
    }

    @Test(dataProvider= "testGetCarrierSamplesBiallelicData")
    public void testGetCarrierSamplesBiallelic(final int ploidy, final Allele refAllele, final Allele altAllele,
                                               final int[] copyNumbers, final String[] expectedSamples) {
        final Set<String> expectedResult = Stream.of(expectedSamples).collect(Collectors.toSet());
        final List<Allele> alleles = new ArrayList<>(2);
        if (refAllele != null) {
            alleles.add(refAllele);
        }
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype;
        if (altAllele != null) {
            if (altAllele.equals(Allele.SV_SIMPLE_DEL)) {
                svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.DEL;
            } else if (altAllele.equals(Allele.SV_SIMPLE_DUP)) {
                svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.DUP;
            } else {
                throw new TestException("Invalid alt allele: " + altAllele.getDisplayString());
            }
            alleles.add(altAllele);
        } else {
            svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.DEL; // default
        }

        // Create genotypes with copy number attribute (and no GT)
        final SVCallRecord recordCN = SVTestUtils.getCNVRecordWithCN(ploidy, alleles, svtype, copyNumbers, GATKSVVCFConstants.COPY_NUMBER_FORMAT);
        final Set<String> resultWithCopyNumber =  recordCN.getCarrierSampleSet();
        Assert.assertEquals(resultWithCopyNumber, expectedResult);

        final SVCallRecord recordRDCN = SVTestUtils.getCNVRecordWithCN(ploidy, alleles, svtype, copyNumbers, GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT);
        final Set<String> resultWithRDCopyNumber =  recordRDCN.getCarrierSampleSet();
        Assert.assertEquals(resultWithRDCopyNumber, expectedResult);

        // Create genotypes with GT (and no copy number attribute)
        final List<Genotype> genotypesWithGenotype = IntStream.range(0, copyNumbers.length)
                .mapToObj(i -> new GenotypeBuilder(String.valueOf(i))
                        .alleles(SVTestUtils.buildBiallelicListWithPloidy(altAllele, refAllele, copyNumbers[i], ploidy))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy)
                        .make())
                .collect(Collectors.toList());
        final SVCallRecord recordWithGenotype = new SVCallRecord("", "chr1", 1000, SVTestUtils.getValidTestStrandA(svtype),
                "chr1", 1999, SVTestUtils.getValidTestStrandB(svtype), svtype, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                alleles, GenotypesContext.copy(genotypesWithGenotype), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Set<String> resultWithGenotype = recordWithGenotype.getCarrierSampleSet();

        Assert.assertEquals(resultWithGenotype, expectedResult);
    }

    @Test
    public void testLargeRandom() {
        final Random rand = new Random(42);
        final List<SVCallRecord> records = new ArrayList<>(1000);
        for (int i = 0; i < 1000; i++) {
            final int pos1 = rand.nextInt(9000) + 1;
            final int pos2 = rand.nextInt(10000) + 1;
            records.add(SVTestUtils.newPESRCallRecordWithIntervalAndType(Math.min(pos1, pos2), Math.max(pos1, pos2), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL));
        }
        final SVClusterEngine engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        final List<SVCallRecord> output = new ArrayList<>(
            records.stream()
                    .sorted(SVCallRecordUtils.getCallComparator(SVTestUtils.hg38Dict))
                    .map(engine::addAndFlush)
                    .flatMap(List::stream)
                    .collect(Collectors.toUnmodifiableList())
        );
        output.addAll(engine.flush());
        Assert.assertEquals(output.size(), 2926);
    }

    @DataProvider(name = "testSampleOverlapData")
    public Object[][] testSampleOverlapData() {
        return new Object[][]{
                {new String[]{}, new String[]{}, 0, true},
                {new String[]{}, new String[]{}, 1, true},
                {new String[]{SAMPLE_1_NAME}, new String[]{}, 0, true},
                {new String[]{SAMPLE_1_NAME}, new String[]{}, 0.1, false},
                {new String[]{SAMPLE_1_NAME}, new String[]{SAMPLE_1_NAME}, 1, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_1_NAME}, 0, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_3_NAME}, 0, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_1_NAME}, 0.5, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_2_NAME}, 0.5, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_3_NAME}, 0.5, false},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_1_NAME}, 1, false},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, 1, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME}, new String[]{SAMPLE_3_NAME, SAMPLE_4_NAME}, 0.1, false},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME, SAMPLE_3_NAME}, new String[]{SAMPLE_1_NAME, SAMPLE_3_NAME, SAMPLE_4_NAME}, 0.5, true},
                {new String[]{SAMPLE_1_NAME, SAMPLE_2_NAME, SAMPLE_3_NAME}, new String[]{SAMPLE_1_NAME, SAMPLE_3_NAME, SAMPLE_4_NAME}, 1, false},
        };
    }

    @Test(dataProvider= "testSampleOverlapData")
    public void testSampleOverlap(final String[] samplesA, final String[] samplesB, final double threshold, final boolean expected) {
        final Set<String> setA = Stream.of(samplesA).collect(Collectors.toUnmodifiableSet());
        final Set<String> setB = Stream.of(samplesB).collect(Collectors.toUnmodifiableSet());
        final List<String> allSamples = List.of(SAMPLE_1_NAME, SAMPLE_2_NAME, SAMPLE_3_NAME, SAMPLE_4_NAME);
        final SVCallRecord recordA = SVTestUtils.makeRecordWithCarriers(allSamples, setA);
        final SVCallRecord recordB = SVTestUtils.makeRecordWithCarriers(allSamples, setB);
        final CanonicalSVLinkage<SVCallRecord> linkage = SVTestUtils.getNewDefaultLinkage();
        linkage.setEvidenceParams(ClusteringParameters.createPesrParameters(0.5, 0.5, 100, threshold));
        Assert.assertEquals(linkage.areClusterable(recordA, recordB).getResult(), expected);
        Assert.assertEquals(linkage.areClusterable(recordB, recordA).getResult(), expected);

        // Test complex event code path
        final SVCallRecord complexA = SVTestUtils.makeComplexRecordWithCarriers(allSamples, setA);
        final SVCallRecord complexB = SVTestUtils.makeComplexRecordWithCarriers(allSamples, setB);
        Assert.assertEquals(linkage.areClusterable(complexA, complexB).getResult(), expected);
        Assert.assertEquals(linkage.areClusterable(complexB, complexA).getResult(), expected);
    }

}