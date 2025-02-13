package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.Histogram;
import java.nio.file.Path;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2FilterConstants.ALL_KMERS;

public class LearnReadOrientationModelEngineUnitTest extends CommandLineProgramTest {
    private static final double EPSILON = 1e-2;

    protected final Logger logger = LogManager.getLogger(this.getClass());

    @DataProvider(name = "computeResponsibilities")
    public Object[][] computeResponsibilitiesData(){
        return new Object[][] {
                {15, 5, 5,
                        Arrays.asList(pair(ArtifactState.HOM_REF, 0.999), pair(ArtifactState.F1R2_T, 1-0.999)),
                        Arrays.asList(pair(ArtifactState.HOM_REF, 0), pair(ArtifactState.F1R2_T, 1)), 1e-3},
                {15, 5, 5,
                        Collections.EMPTY_LIST,
                        Arrays.asList(pair(ArtifactState.HOM_REF, 0), pair(ArtifactState.F1R2_T, 0.7)), 2e-1},
                {50, 10, 10,
                        Collections.EMPTY_LIST,
                        Arrays.asList(pair(ArtifactState.F1R2_T, 1.0)), 1e-2},
                {50, 10, 5,
                        Collections.EMPTY_LIST,
                        Arrays.asList(pair(ArtifactState.SOMATIC_HET, 1.0)), 2e-1},
                {50, 25, 12,
                        Collections.EMPTY_LIST,
                        Arrays.asList(pair(ArtifactState.GERMLINE_HET, 0.5)), 2e-1}
        };
    }

    @Test(dataProvider="computeResponsibilities")
    public void testComputeResponsibilities(final int refDepth, final int altDepth, final int altF1R2,
                                                 final List<ImmutablePair<ArtifactState, Double>> pairs,
                                                 final List<ImmutablePair<ArtifactState, Double>> expectedPairs,
                                                 final double epsilon){
        double[] prior = new double[ArtifactState.values().length];
        for (ImmutablePair<ArtifactState, Double> pair : pairs){
            final ArtifactState state = pair.getLeft();
            final double probability = pair.getRight();
            prior[state.ordinal()] = probability;
        }

        final Nucleotide refAllele = Nucleotide.C;
        final Nucleotide altAllele = Nucleotide.T;

        if (pairs.isEmpty()) {
            prior = LearnReadOrientationModelEngine.getFlatPrior(refAllele);
        }

        final int depth = altDepth + refDepth;
        final double[] normalizedPosterior = LearnReadOrientationModelEngine.computeResponsibilities(refAllele, altAllele, altDepth, altF1R2, depth, prior, false);
        for (ImmutablePair<ArtifactState, Double> expectedPair : expectedPairs){
            final ArtifactState state = expectedPair.getLeft();
            final double probability = expectedPair.getRight();
            Assert.assertEquals(normalizedPosterior[state.ordinal()], probability, epsilon);
        }
    }

    private ImmutablePair<ArtifactState, Double> pair(final ArtifactState state, final double probability){
        return new ImmutablePair<>(state, probability);
    }

    /**
     * 100 sites for a single context AGT. G -> T transition observed in NUM_F1R2_EXAMPLES sites
     * All of the alt sites have 100% Alt F1R2 reads. The rest of the examples are hom ref.
     */
    @Test
    public void testSimpleCase() {
        final int numRefExamples = 80;
        final int numAltExamples = 20;
        final int numExamples = numRefExamples + numAltExamples;
        final int refCount = 50;
        final int altCount = 10;
        final int refF1R2 = refCount/2;
        final int depth = refCount + altCount;

        final String refContext = "ACT";
        final Nucleotide altAllele = Nucleotide.T;
        final ArtifactState artifactState = ArtifactState.F1R2_T;


        final List<AltSiteRecord> altDesignMatrix = new ArrayList<>();
        final Histogram<Integer> refSiteHistogram = F1R2FilterUtils.createRefHistogram(refContext, F1R2FilterConstants.DEFAULT_MAX_DEPTH);

        IntStream.range(0, numRefExamples).forEach(i -> refSiteHistogram.increment(depth));
        IntStream.range(0, numAltExamples).forEach(i -> altDesignMatrix.add(new AltSiteRecord(refContext, refCount, altCount, refF1R2, altCount, altAllele)));

        final LearnReadOrientationModelEngine engine = new LearnReadOrientationModelEngine(refSiteHistogram, Collections.emptyList(), altDesignMatrix,
                LearnReadOrientationModel.DEFAULT_CONVERGENCE_THRESHOLD, LearnReadOrientationModel.DEFAULT_MAX_ITERATIONS,
                F1R2FilterConstants.DEFAULT_MAX_DEPTH, logger);
        final ArtifactPrior artifactPrior = engine.learnPriorForArtifactStates();

        IntStream.range(0, F1R2FilterConstants.DEFAULT_MAX_DEPTH).mapToDouble(i -> MathUtils.sum(engine.getRefResonsibilities(i)))
                .forEach(s -> Assert.assertEquals(s, 1.0, EPSILON));
        IntStream.range(0, numAltExamples).mapToDouble(i -> MathUtils.sum(engine.getAltResonsibilities(i)))
                .forEach(s -> Assert.assertEquals(s, 1.0, EPSILON));
        Assert.assertEquals(engine.getEffectiveCounts().getL1Norm(), numExamples, EPSILON);

        Assert.assertEquals(engine.getEffectiveCounts(artifactState), (double) numAltExamples, EPSILON);
        Assert.assertEquals(engine.getEffectiveCounts(ArtifactState.HOM_REF), (double) numRefExamples, EPSILON);

        Assert.assertEquals(artifactPrior.getPi(artifactState), (double) numAltExamples/numExamples, EPSILON);
        Assert.assertEquals(artifactPrior.getPi(ArtifactState.HOM_REF), (double) numRefExamples/numExamples, EPSILON);
    }


    /**
     * Now test the case where not all of the transitions have orientation bias. And for transitions that do sometimes exhibit
     * orientation bias filter, sill assumes single context.
     */
    @Test
    public void testMoreComplicatedCase() {
        final double epsilon = 0.2;
        final String refContext = "ACT";
        final Nucleotide refAllele = Nucleotide.C;
        final Nucleotide altAllele = Nucleotide.T;
        final int numRefExamples = 10_000;
        final int numArtifactExamples = 10;

        final int randomSeed = 42;
        final RandomGenerator rdg = RandomGeneratorFactory.createRandomGenerator(new Random(randomSeed));

        final BinomialDistribution refDistribution = new BinomialDistribution(rdg, F1R2FilterConstants.DEFAULT_MAX_DEPTH, 0.5);

        final List<AltSiteRecord> altDesignMatrix = new ArrayList<>();
        final Histogram<Integer> refSiteHistogram = F1R2FilterUtils.createRefHistogram(refContext, F1R2FilterConstants.DEFAULT_MAX_DEPTH);

        // Create ref examples
        IntStream.range(0, numRefExamples).forEach(i -> refSiteHistogram.increment(refDistribution.sample()));

        // Create alt examples with signature ACT -> ATT F1R2: assume 25% allele fraction
        final double artfiactAF = 0.25;
        IntStream.range(0, numArtifactExamples).forEach(i -> {
            final int refCount = refDistribution.sample();
            final int altCount = new BinomialDistribution(rdg, refCount, artfiactAF).sample();
            altDesignMatrix.add(new AltSiteRecord(refContext, refCount, altCount, refCount/2, altCount, altAllele));
        });

        // This is a different artifact entirely: ACT -> AAT F1R2
        final Nucleotide altAllele2 = Nucleotide.A;
        IntStream.range(0, numArtifactExamples).forEach(i -> {
            final int refCount = refDistribution.sample();
            final int altCount = new BinomialDistribution(rdg, refCount, artfiactAF).sample();
            altDesignMatrix.add(new AltSiteRecord(refContext, refCount, altCount, refCount/2, altCount, altAllele2));
        });

        // Inject some somatic hets
        final Nucleotide[] possibleAlts = new Nucleotide[]{ Nucleotide.A, Nucleotide.C, Nucleotide.T };
        final int numSomaticHetExamples = 100;
        IntStream.range(0, numSomaticHetExamples).forEach(i -> {
            final int refCount = refDistribution.sample();
            final int altCount = new BinomialDistribution(rdg, refCount, artfiactAF).sample();
            final Nucleotide somaticAlt = possibleAlts[new UniformIntegerDistribution(rdg, 0, possibleAlts.length-1).sample()];
            altDesignMatrix.add(new AltSiteRecord(refContext, refCount, altCount, refCount/2, altCount/2, somaticAlt));
        });

        final LearnReadOrientationModelEngine engine = new LearnReadOrientationModelEngine(refSiteHistogram, Collections.emptyList(), altDesignMatrix,
                LearnReadOrientationModel.DEFAULT_CONVERGENCE_THRESHOLD, LearnReadOrientationModel.DEFAULT_MAX_ITERATIONS,
                F1R2FilterConstants.DEFAULT_MAX_DEPTH, logger);
        final ArtifactPrior artifactPrior = engine.learnPriorForArtifactStates();

        Assert.assertEquals(engine.getEffectiveCounts(ArtifactState.F1R2_T), (double) numArtifactExamples, epsilon);
        Assert.assertEquals(engine.getEffectiveCounts(ArtifactState.F1R2_A), (double) numArtifactExamples, epsilon);
        Assert.assertEquals(engine.getEffectiveCounts(ArtifactState.SOMATIC_HET), (double) numSomaticHetExamples, epsilon);

        final int numExamples = numRefExamples + 2*numArtifactExamples + numSomaticHetExamples;
        Assert.assertEquals(artifactPrior.getPi(ArtifactState.SOMATIC_HET), (double) numSomaticHetExamples/numExamples, 1e-3);
        Assert.assertEquals(MathUtils.sum(artifactPrior.getPi()), 1.0, 1e-3);

        // Make sure ref->ref transitions get 0 probability
        for (ArtifactState state : ArtifactState.getRefToRefArtifacts(refAllele)){
            Assert.assertEquals(artifactPrior.getPi(state), 0.0);
        }
    }

    @Test
    public void testReverseComplement() throws IOException {
        final File f1r2Dir = createTempDir("f1r2");
        final String sample = "sample";
        final File refHistogramOutput = new File(f1r2Dir, sample + F1R2CountsCollector.REF_HIST_EXTENSION);

        final MetricsFile<?, Integer> refMetricsFile = new MetricsFile<>();

        final double epsilon = 1e-3;
        final int numRefExamples = 10_000;
        final int refDepth = 150;
        final int numArtifactExamples = 10;
        final int altDepth = 20;

        final int expectedNumUniqueContexts = 4;
        final List<Triple<String, Nucleotide, ReadOrientation>> transitions = Arrays.asList(
                new ImmutableTriple<>("TCT", Nucleotide.G, ReadOrientation.F1R2), // equivalent to AGA->C F2R1
                new ImmutableTriple<>("AGA", Nucleotide.C, ReadOrientation.F1R2), // this isn't the same as AGA->C F2R1 above
                new ImmutableTriple<>("CAT", Nucleotide.T, ReadOrientation.F2R1), // equivalent to ATG->A F1R2
                new ImmutableTriple<>("ATG", Nucleotide.A, ReadOrientation.F1R2),
                new ImmutableTriple<>("TGG", Nucleotide.C, ReadOrientation.F1R2), // should be recorded as CCA->G F2R1
                new ImmutableTriple<>("GAT", Nucleotide.G, ReadOrientation.F2R1)); // should be recorded as ATC->C F1R2

        final Map<String, List<AltSiteRecord>> altDesignMatricesByContext = ALL_KMERS.stream()
                .collect(Collectors.toMap(kmer -> kmer, kmer -> new ArrayList<>()));

        final Map<String, Histogram<Integer>> refSiteHistograms = new HashMap<>(ALL_KMERS.size());
        ALL_KMERS.forEach(context -> {
            Histogram<Integer> emptyRefHistogram = F1R2FilterUtils.createRefHistogram(context, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
            refSiteHistograms.put(context, emptyRefHistogram);
        });

        for (final Triple<String, Nucleotide, ReadOrientation> transition : transitions) {
            final String refContext = transition.getLeft();
            final Nucleotide altAllele = transition.getMiddle();
            final ReadOrientation f1r2 = transition.getRight();

            refSiteHistograms.get(refContext).increment(refDepth, numRefExamples);

            final int refF1R2 = refDepth / 2;
            final int altF1R2 = f1r2 == ReadOrientation.F1R2 ? altDepth : 0;

            IntStream.range(0, numArtifactExamples).forEach(i -> altDesignMatricesByContext.get(refContext).add(new AltSiteRecord(refContext, refDepth, altDepth, refF1R2, altF1R2, altAllele)));
        }

        refSiteHistograms.values().forEach(hist -> refMetricsFile.addHistogram(hist));
        refMetricsFile.addHeader(new StringHeader(sample));
        refMetricsFile.write(refHistogramOutput);

        for (final String context : ALL_KMERS) {
            final List<AltSiteRecord> records = altDesignMatricesByContext.get(context);
            if (!records.isEmpty()) {
                final File altMatrixOutput = new File(f1r2Dir, sample + "." + context + F1R2CountsCollector.ALT_TABLE_EXTENSION);
                final AltSiteRecord.AltSiteRecordTableWriter altTableWriter =
                        new AltSiteRecord.AltSiteRecordTableWriter(altMatrixOutput.toPath(), sample);
                altTableWriter.writeAllRecords(records);
                altTableWriter.close();
            }
        }

        final File inputTarGz = createTempFile("input", ".tar.gz");
        IOUtils.writeTarGz(inputTarGz.getAbsolutePath(), f1r2Dir.listFiles());

        final File artifactPriorTarGz = createTempFile("prior", ".tar.gz");
        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-I", inputTarGz.getAbsolutePath(),
                        "-O", artifactPriorTarGz.getAbsolutePath()),
                LearnReadOrientationModel.class.getSimpleName()));

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(artifactPriorTarGz.toPath(), extractedDir.toPath());


        final ArtifactPriorCollection artifactPriorCollection = ArtifactPriorCollection.readArtifactPriors(Arrays.stream(extractedDir.listFiles()).findFirst().get());

        Assert.assertEquals(artifactPriorCollection.getNumUniqueContexts(), expectedNumUniqueContexts);

        final double expectedRefFraction = (double) numRefExamples/(numArtifactExamples + numRefExamples);
        final double expectedArtifactFraction = (double) numArtifactExamples/(numArtifactExamples + numRefExamples);

        for (final Triple<String, Nucleotide, ReadOrientation> transition : transitions) {
            final String refContext = transition.getLeft();
            final Nucleotide altAllele = transition.getMiddle();
            final ReadOrientation f1r2 = transition.getRight();
            final ArtifactState state = f1r2 == ReadOrientation.F1R2 ? ArtifactState.getF1R2StateForAlt(altAllele) :
                    ArtifactState.getF2R1StateForAlt(altAllele);
            final ArtifactPrior artifactPrior = artifactPriorCollection.get(refContext).get();
            Assert.assertEquals(artifactPrior.getPi(state), expectedArtifactFraction, epsilon);
            Assert.assertEquals(artifactPrior.getPi(ArtifactState.HOM_REF), expectedRefFraction, epsilon);
        }
    }

    @Test
    public void testMergeHistograms(){
        final int numExamples1 = 10000;
        final int numExamples2 = 5000;
        final int numExamples = numExamples1 + numExamples2;
        final int depth1 = 100;
        final int depth2 = 50;

        // Test Merging Reference Histograms
        final Histogram<Integer> refF1R2TCT = createRefHistograms("TCT", depth1, numExamples1);
        final Histogram<Integer> refF1R2AGA = createRefHistograms("AGA", depth2, numExamples2);

        final Histogram<Integer> combinedRefAGA = LearnReadOrientationModel.combineRefHistogramWithRC(
                "AGA", refF1R2AGA, refF1R2TCT, F1R2FilterConstants.DEFAULT_MAX_DEPTH);

        Assert.assertEquals((int) combinedRefAGA.get(depth1).getValue(), numExamples1);
        Assert.assertEquals((int) combinedRefAGA.get(depth2).getValue(), numExamples2);
        Assert.assertEquals((int) combinedRefAGA.getSumOfValues(), numExamples);

        // Test Merging Alt Histograms i.e. histograms for alt sites with alt depth = 1
        final List<Histogram<Integer>> altTCT = createDepthOneAltHistograms("TCT", depth1, numExamples1);
        final List<Histogram<Integer>> altAGA = createDepthOneAltHistograms("AGA", depth2, numExamples2);

        final List<Histogram<Integer>> combinedDepthOneAGA = LearnReadOrientationModel.combineAltDepthOneHistogramWithRC(altAGA, altTCT, F1R2FilterConstants.DEFAULT_MAX_DEPTH);

        // Feeding the input the other way should fail
        try {
            LearnReadOrientationModel.combineAltDepthOneHistogramWithRC(altTCT, altAGA, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
        } catch (IllegalArgumentException e){
            // Good
        }

        Assert.assertEquals(combinedDepthOneAGA.size(), F1R2FilterConstants.numAltHistogramsPerContext);
        combinedDepthOneAGA.stream().forEach(h -> {
            Assert.assertEquals((int) h.getSumOfValues(), numExamples1 + numExamples2);
            Assert.assertEquals((int) h.get(depth1).getValue(), numExamples1);
            Assert.assertEquals((int) h.get(depth2).getValue(), numExamples2);
        });
    }


    @Test
    public void testMergeDesignMatrices(){
        final int numExamples = 1000;
        final int refDepth = 150;
        final int altDepth = 20;
        final int refF1R2 = refDepth/2;

        final Triple<String, Nucleotide, ReadOrientation> transitions1a = new ImmutableTriple<>("AGA", Nucleotide.C, ReadOrientation.F1R2);
        final int altF1R2a = 10;
        final Triple<String, Nucleotide, ReadOrientation> transitions1b = new ImmutableTriple<>("TCT", Nucleotide.G, ReadOrientation.F2R1);
        final int altF1R2b = 20;
        final Triple<String, Nucleotide, ReadOrientation> transitions1c = new ImmutableTriple<>("TCT", Nucleotide.T, ReadOrientation.F1R2);
        final int altF1R2c = 0;
        final Triple<String, Nucleotide, ReadOrientation> transitions1d = new ImmutableTriple<>("TGG", Nucleotide.C, ReadOrientation.F1R2);

        final List<AltSiteRecord> altDesignMatrix1a = createDesignMatrixOfSingleContext(transitions1a, refDepth, altDepth, refF1R2, altF1R2a, numExamples);
        final List<AltSiteRecord> altDesignMatrix1b = createDesignMatrixOfSingleContext(transitions1b, refDepth, altDepth, refF1R2, altF1R2b, numExamples);
        final List<AltSiteRecord> altDesignMatrix1c = createDesignMatrixOfSingleContext(transitions1c, refDepth, altDepth, refF1R2, altF1R2c, numExamples);
        final List<AltSiteRecord> altDesignMatrix1d = createDesignMatrixOfSingleContext(transitions1d, refDepth, altDepth, refF1R2, altF1R2c, numExamples);


        // Should be all AGA->C
        LearnReadOrientationModel.mergeDesignMatrices(altDesignMatrix1a, altDesignMatrix1b);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> a.getReferenceContext().equals(transitions1a.getLeft())).count(), 2*numExamples);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> a.getAltAllele() == Nucleotide.C).count(), 2*numExamples);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> a.getAltAllele() == Nucleotide.C).count(), 2*numExamples);

        // Now add the third, distinct transition
        LearnReadOrientationModel.mergeDesignMatrices(altDesignMatrix1a, altDesignMatrix1c);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> a.getAltAllele() == Nucleotide.C).count(), 2*numExamples);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> a.getAltAllele() == Nucleotide.A).count(), numExamples);
        // Add more tests here?

        try {
            // Merging the wrong direction should throw an error
            LearnReadOrientationModel.mergeDesignMatrices(altDesignMatrix1b, altDesignMatrix1a);
            Assert.fail();
        } catch (IllegalArgumentException e){
            // Good
        }

        try {
            // Merging non-matching contexts should throw an error
            LearnReadOrientationModel.mergeDesignMatrices(altDesignMatrix1a, altDesignMatrix1d);
            Assert.fail();
        } catch (IllegalArgumentException e){
            // Good
        }

    }

    private List<AltSiteRecord> createDesignMatrixOfSingleContext(final Triple<String, Nucleotide, ReadOrientation> transition,
                                                                  final int refDepth, final int altDepth,
                                                                  final int refF1R2, final int altF1R2,
                                                                  final int numExamples) {
        final List<AltSiteRecord> altDesignMatrix = new ArrayList<>(numExamples);

        final String refContext = transition.getLeft();
        final Nucleotide altAllele = transition.getMiddle();

        IntStream.range(0, numExamples).forEach(i ->
                altDesignMatrix.add(new AltSiteRecord(refContext, refDepth, altDepth, refF1R2, altF1R2, altAllele)));

        return altDesignMatrix;
    }

    private Histogram<Integer> createRefHistograms(final String refContext, final int refDepth, final int numExamples) {
        final Histogram<Integer> refSiteHistogram = F1R2FilterUtils.createRefHistogram(refContext, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
        refSiteHistogram.increment(refDepth, numExamples);
        return refSiteHistogram;
    }


    /**
     * For a given reference context, create alt depth-1 histograms for all possible transitions
     */
    private List<Histogram<Integer>> createDepthOneAltHistograms(final String refContext, final int depth, final int numExamples) {
        final List<Histogram<Integer>> altComputationalHistograms = new ArrayList<>(F1R2FilterConstants.numAltHistogramsPerContext);
        for (final Nucleotide altAllele : Nucleotide.STANDARD_BASES){
            if (altAllele == F1R2FilterUtils.getMiddleBase(refContext)){
                continue;
            }

            for (final ReadOrientation f1r2 : ReadOrientation.values()){
                final Histogram<Integer> altComputationalHistogram = F1R2FilterUtils.createAltHistogram(refContext, altAllele, f1r2, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
                altComputationalHistogram.increment(depth, numExamples);
                altComputationalHistograms.add(altComputationalHistogram);
            }
        }

        return altComputationalHistograms;
    }
}