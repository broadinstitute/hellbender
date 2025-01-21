package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.vcf.GenotypeConcordanceStates;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class SVConcordanceAnnotatorTest {

    @DataProvider(name = "testMatchCnvData")
    public Object[][] testMatchCnvData() {
        return new Object[][]{

                // Note we only match on CN of truth genotypes with respect to the val expected CN

                ////////////////////////////////////////
                // CNV - Ploidy 0
                ////////////////////////////////////////

                // truth equal
                {
                        0,                                              // eval expected copy number (ie ploidy)
                        0,                                              // eval copy number
                        0,                                              // truth copy number
                        true                                            // expected result
                },

                // truth non-equal
                {
                        0,
                        0,
                        1,
                        false
                },
                {
                        0,
                        0,
                        2,
                        false
                },
                {
                        0,
                        0,
                        3,
                        false
                },

                ////////////////////////////////////////
                // CNV - eval ploidy 1
                ////////////////////////////////////////

                // truth ref
                {
                        1,
                        1,
                        1,
                        true
                },

                // truth var
                {
                        1,
                        1,
                        0,
                        false
                },
                {
                        1,
                        1,
                        2,
                        false
                },
                // same var as eval
                {
                        1,
                        0,
                        0,
                        true
                },
                {
                        1,
                        2,
                        2,
                        true
                },
                // different CN as eval and therefore false
                {
                        1,
                        2,
                        3,
                        false
                },
                {
                        1,
                        2,
                        1,
                        false
                },
                {
                        2,
                        1,
                        0,
                        false
                },

                ////////////////////////////////////////
                // CNV - eval ploidy 2
                ////////////////////////////////////////

                // truth hom ref
                {
                        2,
                        2,
                        2,
                        true
                },

                // true - same types
                {
                        2,
                        1,
                        1,
                        true
                },
                {
                        2,
                        3,
                        3,
                        true
                },
                // false - different types or different CN
                {
                        2,
                        1,
                        3,
                        false
                },
                {
                        2,
                        1,
                        0,
                        false
                },
                {
                        2,
                        3,
                        1,
                        false
                },
                {
                        2,
                        3,
                        4,
                        false
                }
        };
    }

    @Test(dataProvider= "testMatchCnvData")
    public void testMatchCnvNonNull(
            final int evalExpectedCopyNumber,
            final int evalCopyNumber,
            final int truthCopyNumber,
            final boolean expected) {

        final SVCallRecord evalRecord = SVTestUtils.newCallRecordWithAlleles(
                Collections.nCopies(evalExpectedCopyNumber, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                evalExpectedCopyNumber,
                evalCopyNumber
        );

        final String sample = "sample";
        final int truthExpectedCopyNumber = 2;  // this shouldn't matter
        final SVCallRecord truthRecord = SVTestUtils.newCallRecordWithAlleles(
                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                truthExpectedCopyNumber,
                truthCopyNumber
        );

        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        final boolean actual = collapser.copyNumbersMatch(sample, evalRecord, truthRecord);
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testMatchCnvNull() {
        final String sample = "sample";
        final SVCallRecord record = SVTestUtils.newCallRecordWithAlleles(
                Arrays.asList(Allele.NO_CALL, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                2,
                2
        );

        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();

        // Null records in null call
        Assert.assertNull(collapser.copyNumbersMatch(sample, record, null));
        Assert.assertNull(collapser.copyNumbersMatch(sample, null, record));


        final SVCallRecord recordNoSample = SVTestUtils.newCallRecordWithAllelesAndSampleName(
                "sample2",
                Arrays.asList(Allele.NO_CALL, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                2,
                2
        );

        // Null call if sample not found
        Assert.assertNull(collapser.copyNumbersMatch(sample, recordNoSample, recordNoSample));
        Assert.assertNull(collapser.copyNumbersMatch(sample, record, recordNoSample));
        Assert.assertNull(collapser.copyNumbersMatch(sample, recordNoSample, record));

        final SVCallRecord recordNoCn = SVTestUtils.newCallRecordWithAllelesAndSampleName(
                sample,
                Arrays.asList(Allele.NO_CALL, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV),
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                2,
                null
        );

        // Null call if CN not found
        Assert.assertNull(collapser.copyNumbersMatch(sample, recordNoCn, recordNoCn));
        Assert.assertNull(collapser.copyNumbersMatch(sample, record, recordNoCn));
        Assert.assertNull(collapser.copyNumbersMatch(sample, recordNoCn, record));
    }

    @DataProvider(name = "testGetStateFromGenotypeData")
    public Object[][] testGetStateFromGenotypeData() {
        return new Object[][]{

                ////////////////////////////////////////
                // Minimal edge cases
                ////////////////////////////////////////

                // no genotype alleles but not a null genotype
                { new Allele[]{}, GenotypeConcordanceStates.TruthState.NO_CALL },

                ////////////////////////////////////////
                // Haploid
                ////////////////////////////////////////

                // no-call
                { new Allele[]{ Allele.NO_CALL }, GenotypeConcordanceStates.TruthState.NO_CALL },
                // hom ref
                { new Allele[]{ Allele.REF_N }, GenotypeConcordanceStates.TruthState.HOM_REF },
                // hom var
                { new Allele[]{ Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.HOM_VAR1 },

                ////////////////////////////////////////
                // Diploid
                ////////////////////////////////////////

                // hom ref
                { new Allele[]{ Allele.REF_N, Allele.REF_N }, GenotypeConcordanceStates.TruthState.HOM_REF },
                // hom var
                { new Allele[]{ Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.HOM_VAR1 },
                // het
                { new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.HET_REF_VAR1 },
                // mixed cases
                { new Allele[]{ Allele.REF_N, Allele.NO_CALL }, GenotypeConcordanceStates.TruthState.IS_MIXED },
                { new Allele[]{ Allele.NO_CALL, Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.IS_MIXED },
        };
    }

    @Test(dataProvider= "testGetStateFromGenotypeData")
    public void testGetStateFromGenotype(final Allele[] genotypeAlleles,
                                         final GenotypeConcordanceStates.TruthState expectedTruthState) {
        final Genotype genotype = alleleArrayToGenotype(genotypeAlleles, null, null);
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        final GenotypeConcordanceStates.TruthState actualTruthState = collapser.getTruthState(genotype);
        Assert.assertEquals(actualTruthState, expectedTruthState);

        // Equivalent call states
        final GenotypeConcordanceStates.CallState expectedCallState;
        if (expectedTruthState == GenotypeConcordanceStates.TruthState.HOM_REF) {
            expectedCallState = GenotypeConcordanceStates.CallState.HOM_REF;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.HET_REF_VAR1) {
            expectedCallState = GenotypeConcordanceStates.CallState.HET_REF_VAR1;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.HOM_VAR1) {
            expectedCallState = GenotypeConcordanceStates.CallState.HOM_VAR1;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.NO_CALL) {
            expectedCallState = GenotypeConcordanceStates.CallState.NO_CALL;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.IS_MIXED) {
            expectedCallState = GenotypeConcordanceStates.CallState.IS_MIXED;
        } else {
            throw new TestException("Unexpected truth state: " + expectedTruthState.name());
        }
        final GenotypeConcordanceStates.CallState actualCallState = collapser.getEvalState(genotype);
        Assert.assertEquals(actualCallState, expectedCallState);
    }

    @Test
    public void testGetTruthNullGenotypeState() {
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        final GenotypeConcordanceStates.TruthState actualTruthState = collapser.getTruthState(null);
        Assert.assertEquals(actualTruthState, GenotypeConcordanceStates.TruthState.NO_CALL);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetNullEvalGenotypeState() {
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator();
        collapser.getEvalState(null);
    }

    @Test
    public void testAnnotateEvalKeepAF() {
        final Map<String, Object> evalAttr = new HashMap<>();
        evalAttr.put(VCFConstants.ALLELE_NUMBER_KEY, 200);
        evalAttr.put(VCFConstants.ALLELE_COUNT_KEY, 3);
        evalAttr.put(VCFConstants.ALLELE_FREQUENCY_KEY, 3./200.);
        final GenotypeBuilder builder = new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL));
        final SVCallRecord evalRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("eval", Collections.singletonList(builder.make()), evalAttr);
        final SVCallRecord truthRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("truth", Collections.singletonList(builder.make()), Collections.emptyMap());
        final ClosestSVFinder.ClosestPair pair = new ClosestSVFinder.ClosestPair(0L, evalRecord, truthRecord);
        final SVCallRecord result = new SVConcordanceAnnotator().annotate(pair);
        Assert.assertTrue(result.getAttributes().containsKey(VCFConstants.ALLELE_NUMBER_KEY));
        Assert.assertTrue(result.getAttributes().containsKey(VCFConstants.ALLELE_COUNT_KEY));
        Assert.assertTrue(result.getAttributes().containsKey(VCFConstants.ALLELE_FREQUENCY_KEY));
        Assert.assertEquals(result.getAttributes().get(VCFConstants.ALLELE_NUMBER_KEY), evalAttr.get(VCFConstants.ALLELE_NUMBER_KEY));
        Assert.assertEquals(result.getAttributes().get(VCFConstants.ALLELE_COUNT_KEY), evalAttr.get(VCFConstants.ALLELE_COUNT_KEY));
        Assert.assertEquals(result.getAttributes().get(VCFConstants.ALLELE_FREQUENCY_KEY), evalAttr.get(VCFConstants.ALLELE_FREQUENCY_KEY));
    }

    @Test
    public void testAnnotateEvalComputeAF() {
        final GenotypeBuilder builder = new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL));
        final SVCallRecord evalRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("eval", Collections.singletonList(builder.make()), Collections.emptyMap());
        final SVCallRecord truthRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("truth", Collections.singletonList(builder.make()), Collections.emptyMap());
        final ClosestSVFinder.ClosestPair pair = new ClosestSVFinder.ClosestPair(0L, evalRecord, truthRecord);
        final SVCallRecord result = new SVConcordanceAnnotator().annotate(pair);
        Assert.assertTrue(result.getAttributes().containsKey(VCFConstants.ALLELE_NUMBER_KEY));
        Assert.assertTrue(result.getAttributes().containsKey(VCFConstants.ALLELE_COUNT_KEY));
        Assert.assertTrue(result.getAttributes().containsKey(VCFConstants.ALLELE_FREQUENCY_KEY));
        Assert.assertEquals(result.getAttributes().get(VCFConstants.ALLELE_NUMBER_KEY), Integer.valueOf(2));
        Assert.assertEquals((int[]) result.getAttributes().get(VCFConstants.ALLELE_COUNT_KEY), new int[]{1});
        Assert.assertEquals((double[]) result.getAttributes().get(VCFConstants.ALLELE_FREQUENCY_KEY), new double[]{1/2.0});
    }

    @Test
    public void testAnnotateCopyTruthAF() {
        final Map<String, Object> truthAttr = new HashMap<>();
        truthAttr.put(VCFConstants.ALLELE_NUMBER_KEY, 200);
        truthAttr.put(VCFConstants.ALLELE_COUNT_KEY, 3);
        truthAttr.put(VCFConstants.ALLELE_FREQUENCY_KEY, 3./200.);
        final GenotypeBuilder builder = new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL));
        final SVCallRecord evalRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("eval", Collections.singletonList(builder.make()), Collections.emptyMap());
        final SVCallRecord truthRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("truth", Collections.singletonList(builder.make()), truthAttr);
        final ClosestSVFinder.ClosestPair pair = new ClosestSVFinder.ClosestPair(0L, evalRecord, truthRecord);
        final SVCallRecord result = new SVConcordanceAnnotator().annotate(pair);
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO), truthAttr.get(VCFConstants.ALLELE_NUMBER_KEY));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO), truthAttr.get(VCFConstants.ALLELE_COUNT_KEY));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO), truthAttr.get(VCFConstants.ALLELE_FREQUENCY_KEY));
    }

    @Test
    public void testAnnotateCalculateTruthAF() {
        final GenotypeBuilder builder = new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL));
        final SVCallRecord evalRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("eval", Collections.singletonList(builder.make()), Collections.emptyMap());
        final SVCallRecord truthRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("truth", Collections.singletonList(builder.make()), Collections.emptyMap());
        final ClosestSVFinder.ClosestPair pair = new ClosestSVFinder.ClosestPair(0L, evalRecord, truthRecord);
        final SVCallRecord result = new SVConcordanceAnnotator().annotate(pair);
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO), Integer.valueOf(2));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO), new int[]{1});
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO), new double[]{1/2.0});
    }

    @Test
    public void testAnnotateOverwriteTruthAF() {
        // Should overwrite existing TRUTH_AF annotations
        final Map<String, Object> evalAttr = new HashMap<>();
        evalAttr.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, 200);
        evalAttr.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, 3);
        evalAttr.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, 1./200.);
        final GenotypeBuilder builder = new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL));
        final SVCallRecord evalRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("eval", Collections.singletonList(builder.make()), Collections.emptyMap());
        final SVCallRecord truthRecord = SVTestUtils.newNamedDeletionRecordWithAttributesAndGenotypes("truth", Collections.singletonList(builder.make()), Collections.emptyMap());
        final ClosestSVFinder.ClosestPair pair = new ClosestSVFinder.ClosestPair(0L, evalRecord, truthRecord);
        final SVCallRecord result = new SVConcordanceAnnotator().annotate(pair);
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO), Integer.valueOf(2));
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO), new int[]{1});
        Assert.assertEquals(result.getAttributes().get(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO), new double[]{1/2.0});
    }

    private Genotype alleleArrayToGenotype(final Allele[] allelesArr,
                                           final Integer expectedCopyNumber,
                                           final Integer copyNumber) {
        final GenotypeBuilder builder = new GenotypeBuilder("");
        if (expectedCopyNumber != null) {
            builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        }
        if (copyNumber != null) {
            builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
        }
        return builder.alleles(Arrays.asList(allelesArr)).make();
    }
}
