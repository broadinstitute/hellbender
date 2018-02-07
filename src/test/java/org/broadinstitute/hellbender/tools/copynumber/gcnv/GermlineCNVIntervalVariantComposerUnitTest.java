package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableIntegerCopyNumber;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link GermlineCNVIntervalVariantComposer}.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class GermlineCNVIntervalVariantComposerUnitTest extends CommandLineProgramTest {
    private static final SimpleInterval TEST_INTERVAL = new SimpleInterval("1", 1, 10000);

    /* sums up to unity in probability space */
    private static final double[] TEST_VALID_LOG_POSTERIOR_VECTOR = new double[] {
            -22.133276724683618,
            -4.3825075766712871,
            -0.012754265709645551,
            -8.6265377789688955,
            -21.918647174298602,
            -22.133276474165779};

    /* expected PLs, copy-number MAP call, and GQ for above */
    private static final List<Integer> TEST_EXPECTED_PLS = new ArrayList<>(Arrays.asList(96, 18, 0, 37, 95, 96));
    private static final int TEST_EXPECTED_MAP_CN = 2;
    private static final int TEST_EXPECTED_GQ = 18;

    /* does not sum up to unity in probability space */
    private static final double[] TEST_INVALID_LOG_POSTERIOR_VECTOR = new double[] {
            -22.133276724683618,
            -4.3825075766712871,
            -0.12754265709645551,
            -8.6265377789688955,
            -21.918647174298602,
            -22.133276474165779};

    private static final String TEST_SAMPLE_NAME = "TEST_SAMPLE_NAME";

    private static final List<IntegerCopyNumberState> TEST_COPY_NUMBER_STATE_LIST = new ArrayList<>();
    private static final IntegerCopyNumberStateCollection TEST_COPY_NUMBER_STATE_COLLECTION;
    static {
        IntStream.range(0, 6).forEach(cn -> TEST_COPY_NUMBER_STATE_LIST.add(new IntegerCopyNumberState(cn)));
        TEST_COPY_NUMBER_STATE_COLLECTION = new IntegerCopyNumberStateCollection(TEST_COPY_NUMBER_STATE_LIST.stream()
                .map(IntegerCopyNumberState::toString)
                .collect(Collectors.toList()));
    }

    private static final CopyNumberPosteriorDistribution TEST_DISTRIBUTION = new CopyNumberPosteriorDistribution(
            IntStream.range(0, TEST_VALID_LOG_POSTERIOR_VECTOR.length)
                    .boxed()
                    .collect(Collectors.toMap(IntegerCopyNumberState::new, cn -> TEST_VALID_LOG_POSTERIOR_VECTOR[cn])));

    private static final LocatableCopyNumberPosteriorDistribution TEST_LOCATABLE_DISTRIBUTION =
            new LocatableCopyNumberPosteriorDistribution(TEST_INTERVAL, TEST_DISTRIBUTION);

    /**
     * Assert that exception is thrown if normalized posterior probabilities do not add up to 1
     */
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCopyNumberPosteriorDistributionValidation() {
        final CopyNumberPosteriorDistribution posteriorRecord = new CopyNumberPosteriorDistribution(
                IntStream.range(0, TEST_VALID_LOG_POSTERIOR_VECTOR.length)
                        .boxed()
                        .collect(Collectors.toMap(IntegerCopyNumberState::new,
                                cn -> TEST_INVALID_LOG_POSTERIOR_VECTOR[cn])));
        final LocatableCopyNumberPosteriorDistribution locatablePosteriorRecord =
                new LocatableCopyNumberPosteriorDistribution(TEST_INTERVAL, posteriorRecord);
        testGenotyping(locatablePosteriorRecord, TEST_COPY_NUMBER_STATE_COLLECTION, null, 0, 0);
    }

    /**
     * Assert that exception is thrown when a copy-number state collection has less than 3 states
     */
    @Test(expectedExceptions = IllegalStateException.class)
    public void testIntegerCopyNumberStatesCollectionValidation() {
        final File outputFile = createTempFile("test", ".vcf");
        final VariantContextWriter outputWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);
        final List<IntegerCopyNumberState> shortCopyNumberStateList = new ArrayList<>(
                Arrays.asList(new IntegerCopyNumberState(0), new IntegerCopyNumberState(1)));
        final IntegerCopyNumberStateCollection shortCopyNumberStateCollection =
                new IntegerCopyNumberStateCollection(shortCopyNumberStateList.stream()
                        .map(IntegerCopyNumberState::toString)
                        .collect(Collectors.toList()));
        final IntegerCopyNumberState testRefAutosomalCopyNumberState = new IntegerCopyNumberState(2);
        final Set<String> testAllosomalContigsSet = Collections.emptySet();
        new GermlineCNVIntervalVariantComposer(outputWriter, shortCopyNumberStateCollection, TEST_SAMPLE_NAME,
                testRefAutosomalCopyNumberState, testAllosomalContigsSet);
    }

    @Test(dataProvider = "alleleDeterminationTestData")
    public void testAlleleDetermination(final int refAutosomalCopyNumber,
                                        final Set<String> allosomalContigs,
                                        final int baselineCopyNumber,
                                        final Allele expectedAllele) {
        final File outputFile = createTempFile("test", ".vcf");
        final LocatableIntegerCopyNumber baselineCopyNumberState = new LocatableIntegerCopyNumber(TEST_INTERVAL,
                new IntegerCopyNumberState(baselineCopyNumber));
        final VariantContextWriter outputWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);
        final GermlineCNVIntervalVariantComposer variantComposer = new GermlineCNVIntervalVariantComposer(
                outputWriter, TEST_COPY_NUMBER_STATE_COLLECTION, "TEST_SAMPLE_NAME",
                new IntegerCopyNumberState(refAutosomalCopyNumber), allosomalContigs);
        final VariantContext var = variantComposer.composeVariantContext(TEST_LOCATABLE_DISTRIBUTION,
                baselineCopyNumberState, "TEST_CNV");
        final List<Allele> allAlleles = var.getAlleles();
        Assert.assertEquals(allAlleles, GermlineCNVIntervalVariantComposer.ALL_ALLELES);
        final Genotype gen = var.getGenotype(TEST_SAMPLE_NAME);
        final Allele genAllele = gen.getAlleles().get(0);
        Assert.assertEquals(genAllele, expectedAllele);
    }

    @Test(dataProvider = "testCopyNumberPosteriorDistributionRecords")
    public void testGenotyping(final LocatableCopyNumberPosteriorDistribution posteriorLocatableRecord,
                               final IntegerCopyNumberStateCollection copyNumberStateCollection,
                               final List<Integer> expectedPLs,
                               final int expectedMAPCopyNumber,
                               final int expectedGQ) {
        final List<Integer> actualPLs = GermlineCNVIntervalVariantComposer.getCopyNumberPLVector(
                posteriorLocatableRecord, copyNumberStateCollection);
        final int actualMAPCopyNumber = GermlineCNVIntervalVariantComposer.calculateMAPCopyNumberState(
                posteriorLocatableRecord, copyNumberStateCollection);
        final int actualGQ = GermlineCNVIntervalVariantComposer.calculateGenotypeQuality(
                posteriorLocatableRecord, copyNumberStateCollection);
        Assert.assertEquals(actualMAPCopyNumber, expectedMAPCopyNumber);
        IntStream.range(0, expectedPLs.size())
                .forEach(i -> Assert.assertEquals(actualPLs.get(i).intValue(), expectedPLs.get(i).intValue()));
        Assert.assertEquals(actualGQ, expectedGQ);
    }

    @DataProvider(name = "testCopyNumberPosteriorDistributionRecords")
    public Object[][] getTestCopyNumberPosteriorDistributionRecords() {
        return new Object[][] {
                {TEST_LOCATABLE_DISTRIBUTION, TEST_COPY_NUMBER_STATE_COLLECTION,
                        TEST_EXPECTED_PLS, TEST_EXPECTED_MAP_CN, TEST_EXPECTED_GQ}};
    }

    @DataProvider(name = "alleleDeterminationTestData")
    public Object[][] getAlleleDeterminationTestData() {
        return new Object[][] {
                /* if the contig of the genotyped interval is not contained in the allosomal contigs set
                 * (of if the latter is empty), the baseline copy-number state is ignored and the ref autosomal
                 * copy-number is used for allele determination */
                {0, new HashSet<String>(), TEST_EXPECTED_MAP_CN - 1, GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {0, new HashSet<String>(), TEST_EXPECTED_MAP_CN + 1, GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {0, new HashSet<String>(), TEST_EXPECTED_MAP_CN, GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {1, new HashSet<String>(), TEST_EXPECTED_MAP_CN - 1, GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {1, new HashSet<String>(), TEST_EXPECTED_MAP_CN + 1, GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {1, new HashSet<String>(), TEST_EXPECTED_MAP_CN, GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {2, new HashSet<String>(), TEST_EXPECTED_MAP_CN - 1, GermlineCNVIntervalVariantComposer.REF_ALLELE},
                {2, new HashSet<String>(), TEST_EXPECTED_MAP_CN + 1, GermlineCNVIntervalVariantComposer.REF_ALLELE},
                {2, new HashSet<String>(), TEST_EXPECTED_MAP_CN, GermlineCNVIntervalVariantComposer.REF_ALLELE},
                {3, new HashSet<String>(), TEST_EXPECTED_MAP_CN - 1, GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                {3, new HashSet<String>(), TEST_EXPECTED_MAP_CN + 1, GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                {3, new HashSet<String>(), TEST_EXPECTED_MAP_CN, GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                /* if the contig of the genotyped interval is contained in the allosomal contigs set,
                 * the baseline copy-number state is used for allele determination and the ref autosomal
                 * copy-number state is ignored */
                {0, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN - 1,
                        GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {0, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN + 1,
                        GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                {0, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN,
                        GermlineCNVIntervalVariantComposer.REF_ALLELE},
                {1, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN - 1,
                        GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {1, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN + 1,
                        GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                {1, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN,
                        GermlineCNVIntervalVariantComposer.REF_ALLELE},
                {2, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN - 1,
                        GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {2, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN + 1,
                        GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                {3, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN,
                        GermlineCNVIntervalVariantComposer.REF_ALLELE},
                {3, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN - 1,
                        GermlineCNVIntervalVariantComposer.DUP_ALLELE},
                {3, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN + 1,
                        GermlineCNVIntervalVariantComposer.DEL_ALLELE},
                {3, Collections.singleton(TEST_INTERVAL.getContig()), TEST_EXPECTED_MAP_CN,
                        GermlineCNVIntervalVariantComposer.REF_ALLELE}
        };
    }
}
