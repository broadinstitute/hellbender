package org.broadinstitute.hellbender.tools.walkers.validation;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class EvaluateGenotypingPerformanceUnitTests extends GATKBaseTest {

    @Test
    void testPearsonCorrelationAggregator() {
        final Random rand = new Random(234);
        final int nEntries = 1000;
        final double[] x = new double[nEntries];
        final double[] y = new double[nEntries];

        final EvaluateGenotypingPerformance.PearsonCorrelationAggregator aggregator = new EvaluateGenotypingPerformance.PearsonCorrelationAggregator();

        for (int i = 0; i < nEntries; i++) {
            x[i] = rand.nextDouble();
            y[i] = rand.nextDouble();

            aggregator.addEntry(x[i], y[i]);
        }

        final PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();
        Assert.assertEquals(aggregator.getCorrelation(), pearsonsCorrelation.correlation(x, y), 1e-10);

    }

    @DataProvider(name = "getConcordanceStateDataProvider")
    Object[][] getConcordanceStateDataProvider() {
        return new Object[][] {
                {Arrays.asList(Allele.ALT_A, Allele.REF_G), Arrays.asList(Allele.REF_G, Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.REF_G), Arrays.asList(Allele.REF_G, Allele.ALT_A), true, ConcordanceState.FILTERED_FALSE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_G, Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_G, Allele.ALT_A), true, ConcordanceState.FILTERED_FALSE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_G, Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A), false, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A), true, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.ALT_C), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.ALT_C), true, ConcordanceState.FILTERED_TRUE_NEGATIVE},

                //differing ploidy, generally from different treatment of chromosome X in males

                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A), false, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_A), false, ConcordanceState.FALSE_NEGATIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A), false, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A), true, ConcordanceState.FILTERED_TRUE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.REF_A), false, ConcordanceState.FALSE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A), true, ConcordanceState.FILTERED_FALSE_NEGATIVE}
        };
    }

    @Test(dataProvider = "getConcordanceStateDataProvider")
    void testGetConcordanceState(final List<Allele> truthAlleles, final List<Allele> evalAlleles, final boolean evalWasFiltered, final ConcordanceState expectedConcordanceState ) {
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
        final Genotype truthGenotype = genotypeBuilder.alleles(truthAlleles).make();
        final Genotype evalGenotype = genotypeBuilder.alleles(evalAlleles).make();
        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        evaluateGenotypingPerformance.allowDifferingPloidy = true;
        final ConcordanceState observedConcordanceState = evaluateGenotypingPerformance.getConcordanceState(truthGenotype, evalGenotype, evalWasFiltered);

        Assert.assertEquals(observedConcordanceState, expectedConcordanceState);
    }

    @DataProvider(name = "getConcordanceStateDifferingPloidyExceptionDataProvider")
    Object[][] getConcordanceStateDifferingPloidyExceptionDataProvider() {
        return new Object[][] {
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A)},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_A)},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.ALT_A)},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_A)},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A)},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.ALT_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A)}
        };
    }

    @Test(dataProvider = "getConcordanceStateDifferingPloidyExceptionDataProvider", expectedExceptions = GATKException.class)
    void testGetConcordanceStateDifferingPloidyException(final List<Allele> truthAlleles, final List<Allele> evalAlleles) {
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
        final Genotype truthGenotype = genotypeBuilder.alleles(truthAlleles).make();
        final Genotype evalGenotype = genotypeBuilder.alleles(evalAlleles).make();

        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        evaluateGenotypingPerformance.getConcordanceState(truthGenotype, evalGenotype, false);
    }

    @DataProvider(name = "buildAFMapDataProvider")
    Object[][] buildAFMapDataProvider() {
        return new Object[][]{
                {
                        ImmutableMap.builder().
                                put("annotation1", Arrays.asList(0.3, 0.4)).
                                put("annotation2", Arrays.asList(0.1, 0.05)).
                                build(),
                        0,
                        ImmutableMap.builder().
                                put("annotation1", 0.3).
                                put("annotation2", 0.1).
                                build()
                },
                {
                        ImmutableMap.builder().
                                put("annotation1", Arrays.asList(0.3, 0.4)).
                                put("annotation2", Arrays.asList(0.1, 0.05)).
                                build(),
                        1,
                        ImmutableMap.builder().
                                put("annotation1", 0.4).
                                put("annotation2", 0.05).
                                build()
                }
        };
    }

    @Test(dataProvider = "buildAFMapDataProvider")
    void testBuildAFMap(final Map<String, List<Double>> annotationsMap, final int altAlleleIndex, final Map<String, Double> expectedMap) {
        final VariantContextBuilder vb = new VariantContextBuilder("snp1", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C));
        for(final Map.Entry<String, List<Double>> entry : annotationsMap.entrySet()) {
            vb.attribute(entry.getKey(), entry.getValue());
        }

        final VariantContext vc = vb.make();

        final Map<String, Double> actualMap = EvaluateGenotypingPerformance.buildAFMapForIndex(vc, annotationsMap.keySet(), altAlleleIndex);
        Assert.assertEquals(actualMap, expectedMap);
    }

    @DataProvider(name = "getBinDataProvider")
    Object[][] getBinDataProvider() {
        return new Object[][] {
                {12, 0.01, 0.075, 5},
                {12, 0.01, 0.002, 0},
                {12, 0.01, 0, 0},
                {12, 0.01, 1, 11},
                {12, 0.01,0.86, 11},
                {12, 0.01,0.5, 10},

                {6,0.3,0.7, 4},
                {21,0.001,0.024, 10}
        };
    }

    @Test(dataProvider = "getBinDataProvider")
    void testGetBinDataProvider(final int nBins, final double firstBinRightEdge, final double val, final int expectedBin) {
        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        evaluateGenotypingPerformance.nBins = nBins;
        evaluateGenotypingPerformance.firstBinRightEdge = firstBinRightEdge;
        final int actualBin = evaluateGenotypingPerformance.getBin(val);
        Assert.assertEquals(actualBin, expectedBin);
    }

}