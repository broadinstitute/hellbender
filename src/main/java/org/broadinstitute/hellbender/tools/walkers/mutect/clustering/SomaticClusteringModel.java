package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import htsjdk.variant.variantcontext.Allele;
import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.VariantContext;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.MathArrays;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A model for the allele fraction spectrum of somatic variation.  There are several clusters in this model:
 * 1) a "background" cluster with a broad Beta distribution of allele fractions; 2) a high-AF cluster with a Beta distribution
 * of allele fractions with essentially all of its support at high AFs; 3) several clusters with discrete allele fractions.
 *
 * Clusters of type 4) correspond to subclones in the absence of CNV.  The type 2) cluster accounts for CNV, minor subclones, and any
 * other factor that makes a discrete clustering model inapplicable.  Type 3) absorbs probability density from homozygous alts and
 * hets in pure tumors.  Such calls are obvious and we don't want to spend a bunch of type 4 clusters on them, spoiling our parsimony.
 * Type 1) is appropriate to include in the model because the likelihoods of sequencing error and real variation are modeled together
 * in the log odds from the somatc likelihoods model of Mutect2.  In contrast, artifactual error such as from sample prep and
 * mapping are not part of the clustering.
 */
public class SomaticClusteringModel {
    protected final Logger logger = LogManager.getLogger(this.getClass());

    private boolean clustersHaveBeenInitialized;

    private static final int MAX_INDEL_SIZE_IN_PRIOR_MAP = 10;
    private static final int NUM_INITIALIZATION_QUANTILES = 50;
    private static final double MIN_QUANTILE_FOR_MAKING_CLUSTER = 0.1;
    private static final int MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER = (int) (MIN_QUANTILE_FOR_MAKING_CLUSTER * NUM_INITIALIZATION_QUANTILES);
    private final Map<Integer, Double> logVariantPriors = new Int2DoubleArrayMap();
    private double logVariantVsArtifactPrior;
    private final OptionalDouble callableSites;

    private static final double INITIAL_HIGH_AF_WEIGHT = 0.01;
    public static final double MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF = 0.9;

    private double REGULARIZING_PSEUDOCOUNT = 1;

    private double[] logClusterWeights;

    private static final int NUM_ITERATIONS = 5;
    private static final int MAX_BINOMIAL_CLUSTERS = 5;

    private static final BetaDistributionShape INITIAL_HIGH_AF_BETA = new BetaDistributionShape(10, 1);
    private static final BetaDistributionShape INITIAL_BACKGROUND_BETA = BetaDistributionShape.FLAT_BETA;

    final List<Datum> data = new ArrayList<>();
    final List<AlleleFractionCluster> clusters = new ArrayList<>();
    private final MutableInt obviousArtifactCount = new MutableInt(0);
    private static final double OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD = 0.9;

    public SomaticClusteringModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).forEach(n -> logVariantPriors.put(n, MTFAC.getLogIndelPrior()));
        logVariantPriors.put(0, MTFAC.getLogSnvPrior());

        logVariantVsArtifactPrior = MTFAC.initialLogPriorOfVariantVersusArtifact;

        final OptionalDouble callableSitesFromStats = mutectStats.stream().filter(stat -> stat.getStatistic().equals(Mutect2Engine.CALLABLE_SITES_NAME))
                .mapToDouble(MutectStats::getValue).findFirst();

        final boolean noCallableSites = callableSitesFromStats.isPresent() && callableSitesFromStats.getAsDouble() < 1;
        if (noCallableSites) {
            logger.warn("No callable sites found in Mutect stats.  Running without the full somatic clustering model.  Something is seriously wrong!");
        }
        callableSites = noCallableSites ? OptionalDouble.empty() : callableSitesFromStats;


        clusters.add(new BetaBinomialCluster(INITIAL_BACKGROUND_BETA));
        clusters.add(new BetaBinomialCluster(INITIAL_HIGH_AF_BETA));
        logClusterWeights = new double[] {Math.log1p(INITIAL_HIGH_AF_WEIGHT), Math.log(INITIAL_HIGH_AF_WEIGHT)};
    }

    /**
     * Adds data to the model for every alternate allele
     * @param tumorADs for all alleles, summed over samples
     * @param tumorLogOdds for alt alleles only
     * @param artifactProbabilities by alt allele, specifically technical artifact probabilities not including sequencing error, contamination, or germline variation
     * @param nonSomaticProbabilities by alt allele, probabilities that the variants are real but not somatic ie germline or contamination
     * @param vc the variant context the data apply to
     */
    public void record(int[] tumorADs, final double[] tumorLogOdds, final List<Double> artifactProbabilities, final List<Double> nonSomaticProbabilities, final VariantContext vc) {
        Utils.validateArg(tumorADs.length == vc.getNAlleles(), "tumorADs must have one entry per allele including the ref allele");
        // get all alt allele indexes for symbolic alleles
        List<Integer> symIndexes = new IndexRange(0, vc.getNAlleles()-1).filter(n -> vc.getAlternateAllele(n).isSymbolic());

        // set tumorAD to 0 for symbolic alleles so it won't contribute to overall AD, ADs include reference allele so add one to the index
        symIndexes.forEach(i -> tumorADs[i+1] = 0);
        final int totalAD = (int) MathUtils.sum(tumorADs);

        // split into separate alt alleles for clustering
        for (int i = 0; i < tumorLogOdds.length; i++) {
            if (!vc.getAlternateAllele(i).isSymbolic()) {
                if (artifactProbabilities.get(i) > OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD) {
                    obviousArtifactCount.increment();
                    continue;
                } else if (nonSomaticProbabilities.get(i) > OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD) {
                    continue;
                }
                data.add(new Datum(tumorLogOdds[i], artifactProbabilities.get(i), nonSomaticProbabilities.get(i), tumorADs[i + 1], totalAD, indelLength(vc, i)));
            }
        }
    }

    public double getLogPriorOfSomaticVariant(final VariantContext vc, final int altIndex) {
        return getLogPriorOfSomaticVariant(indelLength(vc, altIndex));
    }

    /**
     * Given that something is not an error from the sequencing machine itself, the log prior probability that it is a
     * technical artifact (such as from sample prep or alignment) and not a real variant.
     */
    public double getLogPriorOfVariantVersusArtifact() {
        return logVariantVsArtifactPrior;
    }

    public double probabilityOfSequencingError(final Datum datum) {
        // calculate the log likelihood of real variation relative to an sequencing error log likelihood of zero
        final double[] logClusterLikelihoods = new IndexRange(0, clusters.size())
                .mapToDouble(c -> logClusterWeights[c] + clusters.get(c).correctedLogLikelihood(datum));
        final double variantLogLikelihood = NaturalLogUtils.logSumExp(logClusterLikelihoods);

        return Mutect2FilteringEngine.posteriorProbabilityOfError(variantLogLikelihood, getLogPriorOfSomaticVariant(datum.getIndelLength()));
    }

    private double probabilityOfSomaticVariant(final Datum datum) {
        final double artifactProb = datum.getArtifactProb();
        final double nonSomaticProb = datum.getNonSequencingErrorProb();
        final double sequencingErrorProb = probabilityOfSequencingError(datum);
        return (1 - artifactProb) * ( 1 - nonSomaticProb) * (1 - sequencingErrorProb);
    }

    private void initializeClusters() {
        Utils.validate(!clustersHaveBeenInitialized, "Clusters have already been initialized.");

        final double[] somaticProbs = data.stream().mapToDouble(this::probabilityOfSomaticVariant).toArray();

        double previousBIC = Double.NEGATIVE_INFINITY;

        for (int cluster = 0; cluster < MAX_BINOMIAL_CLUSTERS; cluster++) {
            final double[] oldLogClusterWeights = Arrays.copyOf(logClusterWeights, logClusterWeights.length);

            final double[] backgroundProbsGivenSomatic = data.stream().mapToDouble(datum -> backgroundProbGivenSomatic(datum.getTotalCount(), datum.getAltCount())).toArray();
            final double[] backGroundProbs = MathArrays.ebeMultiply(somaticProbs, backgroundProbsGivenSomatic);
            final double[] alleleFractionQuantiles = calculateAlleleFractionQuantiles();

            // calculate how much total probability is assigned to the background cluster, then split off a peak from the background
            final double[] totalQuantileBackgroundResponsibilities = calculateQuantileBackgroundResponsibilities(alleleFractionQuantiles, backGroundProbs);

            final List<Pair<Double, Double>> peaksAndMasses = calculatePeaksAndMasses(alleleFractionQuantiles, totalQuantileBackgroundResponsibilities);

            if (peaksAndMasses.isEmpty()) {
                break;
            }

            final Pair<Double, Double> biggestPeakAndMass = peaksAndMasses.stream().sorted(Comparator.comparingDouble(Pair<Double, Double>::getRight).reversed()).findFirst().get();

            if (biggestPeakAndMass.getLeft() < alleleFractionQuantiles[Math.min(MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER, alleleFractionQuantiles.length-1)]) {
                break;
            }

            final double totalMass = peaksAndMasses.stream().mapToDouble(Pair::getRight).sum();
            final double fractionOfBackgroundToSplit = Math.min(MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF, biggestPeakAndMass.getRight() / totalMass);
            final double newClusterLogWeight = Math.log(fractionOfBackgroundToSplit) + logClusterWeights[0];
            final double newBackgroundWeight = Math.log1p(fractionOfBackgroundToSplit) + logClusterWeights[0];

            clusters.add(new BinomialCluster(biggestPeakAndMass.getLeft()));

            final List<Double> newLogWeights = new ArrayList<>(Doubles.asList(logClusterWeights));
            newLogWeights.add(newClusterLogWeight);
            newLogWeights.set(0, newBackgroundWeight);
            logClusterWeights = Doubles.toArray(newLogWeights);

            for (int n = 0; n < NUM_ITERATIONS; n++) {
                performEMIteration(false);
            }

            final double[] logLikelihoodsGivenSomatic = data.stream().mapToDouble(datum -> logLikelihoodGivenSomatic(datum.getTotalCount(), datum.getAltCount())).toArray();

            final double weightedLogLikelihood = MathUtils.sum(MathArrays.ebeMultiply(somaticProbs, logLikelihoodsGivenSomatic));
            final double effectiveSomaticCount = MathUtils.sum(somaticProbs);
            final double numParameters = 2 * clusters.size();


            // if splitting off the peak worsened the BIC score, remove the new peak and we're done
            final double currentBIC = weightedLogLikelihood - numParameters * Math.log(effectiveSomaticCount);
            if (currentBIC < previousBIC) {
                clusters.remove(clusters.size() - 1);
                logClusterWeights = oldLogClusterWeights;
                break;
            }
            previousBIC = currentBIC;
        }

        clustersHaveBeenInitialized = true;
    }

    /**
     * Given the current model, order data by allele fraction and calculate the probability that each datum is somatic
     * Then compute probability-weighted quantiles of allele fraction.  That is, each quantile range of allele fraction
     * contains 1 / NUM_INITIALIZATION_QUANTILES of the total somatic posterior probability.
     * @return
     */
    private double[] calculateAlleleFractionQuantiles() {
        // sort data with probabilities by allele fraction
        final List<Pair<Double, Double>> alleleFractionsAndSomaticProbs = data.stream()
                .map(d -> ImmutablePair.of((double) d.getAltCount() / d.getTotalCount(), probabilityOfSomaticVariant(d)))
                .sorted(Comparator.comparingDouble(p -> p.getLeft()))
                .collect(Collectors.toList());

        final double totalSomaticProb = alleleFractionsAndSomaticProbs.stream().mapToDouble(p -> p.getRight()).sum();

        double cumulativeProb = 0;
        final double quantileStep = totalSomaticProb / NUM_INITIALIZATION_QUANTILES;
        double quantileProb = quantileStep;

        final List<Double> alleleFractionQuantilesList = new ArrayList<>(NUM_INITIALIZATION_QUANTILES);

        for (int n = 0; n < data.size(); n++) {
            cumulativeProb += alleleFractionsAndSomaticProbs.get(n).getRight();

            if (cumulativeProb > quantileProb ||  n == data.size()-1) {
                alleleFractionQuantilesList.add(alleleFractionsAndSomaticProbs.get(n).getLeft());
                while (cumulativeProb > quantileProb) {
                    quantileProb += quantileStep;
                }
            }
        }

        return Doubles.toArray(alleleFractionQuantilesList.stream().distinct().collect(Collectors.toList()));
    }

    /**
     * given the current model and an array of allele fractions, calculate the sum over all data of posterior probability
     * densities that each datum is a somatic variant at a given allele fraction, weighted by the probability that data
     * come from the background cluster.
     * @param alleleFractionQuantiles
     * @return
     */
    private double[] calculateQuantileBackgroundResponsibilities(final double[] alleleFractionQuantiles, final double[] backgroundProbs) {
        final double[] totalQuantileResponsibilities = new double[alleleFractionQuantiles.length];
        for (int n = 0; n < data.size(); n++) {
            final Datum datum = data.get(n);
            final double backgroundProb = backgroundProbs[n];

            // responsibilities for each AF quantile will be the posterior probability density of binomial likelihood
            // with a flat prior on AF from 0 to 1.  That is, the unnormalized posterior density at allele fraction f
            // is the binomial emission likelihood f^alt * (1-f)^ref, which we normalize by the integral of this likelihood
            // over f from 0 to 1, which is the beta function B(alt + 1, ref + 1).  Thus the posterior density at f is
            // f^alt * (1-f)^ref * Gamma(alt + ref + 2) / [Gamma(alt + 1) * Gamma(ref + 1)]
            // since alt and ref counts are integers this can be written in terms of factorials as
            // f^alt * (1-f)^ref * (alt + ref + 1)! / [alt! * ref!], which is equal to
            // (alt + ref + 1) times the binomial probability Binom(alt | ref, f)
            final double[] quantileResponsibilities = MathUtils.applyToArray(alleleFractionQuantiles, f -> MathUtils.binomialProbability(datum.getTotalCount(), datum.getAltCount(), f));
            MathUtils.applyToArrayInPlace(quantileResponsibilities, x -> x * backgroundProb * (datum.getTotalCount() + 1));
            MathUtils.addToArrayInPlace(totalQuantileResponsibilities, quantileResponsibilities);
        }
        return totalQuantileResponsibilities;
    }

    private List<Pair<Double, Double>> calculatePeaksAndMasses(double[] alleleFractionQuantiles, double[] totalQuantileResponsibilities) {
        // peak masses will be the trapezoid rule quadratures of responsibility between consecutive minima.
        final List<Pair<Double, Double>> peaksAndMasses = new ArrayList<>();
        double currentPeakMass = 0;
        double currentPeak = 0;
        double currentPeakResponsibility = 0;

        for (int q = 0; q < alleleFractionQuantiles.length; q++) {
            //  0 and 1 have probability density 0.
            final double leftResponsibility = q == 0? 0 : totalQuantileResponsibilities[q - 1];
            final double responsibility = totalQuantileResponsibilities[q];
            final double rightResponsibility = q == alleleFractionQuantiles.length - 1 ? 0 : totalQuantileResponsibilities[q + 1];

            final double leftAlleleFraction = q == 0 ? 0 : alleleFractionQuantiles[q-1];
            final double alleleFraction = alleleFractionQuantiles[q];

            currentPeakMass += (alleleFraction - leftAlleleFraction) * (leftResponsibility + responsibility) / 2;   //trapezoid rule
            if (responsibility > currentPeakResponsibility) {
                currentPeak = alleleFraction;
                currentPeakResponsibility = responsibility;
            }

            final int leftCompare = Double.compare(responsibility, leftResponsibility);
            final int rightCompare = Double.compare(responsibility, rightResponsibility);
            final boolean localMin = leftCompare < 0 && rightCompare <= 0 || leftCompare <= 0 && rightCompare < 0;

            // record the current peak and start a new one
            if (localMin && q > 0 || q == alleleFractionQuantiles.length - 1) {
                peaksAndMasses.add(ImmutablePair.of(currentPeak, currentPeakMass));
                currentPeakMass = 0;
                currentPeak = alleleFraction;
                currentPeakResponsibility = responsibility;
            }
        }
        return peaksAndMasses;
    }

    private double getLogPriorOfSomaticVariant(final int indelLength) {
        if (!logVariantPriors.containsKey(indelLength)) {
            logVariantPriors.put(indelLength, logVariantPriors.values().stream().mapToDouble(d -> d).min().getAsDouble());
        }

        return logVariantPriors.get(indelLength) + (indelLength == 0 ? MathUtils.LOG_ONE_THIRD : 0);
    }

    public void learnAndClearAccumulatedData() {
        if (!clustersHaveBeenInitialized) {
            initializeClusters();
        }

        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {

            performEMIteration(true);

        }

        data.clear();
        obviousArtifactCount.setValue(0);
    }

    private void performEMIteration(final boolean updateSomaticPriors) {
        final Map<Integer, MutableDouble> variantCountsByIndelLength = IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).boxed()
                .collect(Collectors.toMap(l -> l, l -> new MutableDouble(0)));
        final List<double[]> responsibilities = new ArrayList<>(data.size());

        final double[] totalClusterResponsibilities = new double[clusters.size()];
        for (final Datum datum : data) {
            final double somaticProb = probabilityOfSomaticVariant(datum);

            final int indelLength = datum.getIndelLength();
            variantCountsByIndelLength.putIfAbsent(indelLength, new MutableDouble(0));
            variantCountsByIndelLength.get(indelLength).add(somaticProb);

            final double[] clusterLogLikelihoods =  new IndexRange(0, clusters.size())
                    .mapToDouble(c -> logClusterWeights[c] + clusters.get(c).logLikelihood(datum.getTotalCount(), datum.getAltCount()));
            final double[] clusterResponsibilitiesIfSomatic = NaturalLogUtils.normalizeFromLogToLinearSpace(clusterLogLikelihoods);
            final double[] clusterResponsibilities = MathArrays.scale(somaticProb, clusterResponsibilitiesIfSomatic);
            MathUtils.addToArrayInPlace(totalClusterResponsibilities, clusterResponsibilities);
            responsibilities.add(clusterResponsibilities);
        }

        MathUtils.applyToArrayInPlace(totalClusterResponsibilities, x -> x + REGULARIZING_PSEUDOCOUNT);
        logClusterWeights = MathUtils.applyToArrayInPlace(MathUtils.normalizeSumToOne(totalClusterResponsibilities), Math::log);
        final double technicalArtifactCount = obviousArtifactCount.getValue() + data.stream().mapToDouble(Datum::getArtifactProb).sum();
        final double variantCount = variantCountsByIndelLength.values().stream().mapToDouble(MutableDouble::doubleValue).sum();

        if (updateSomaticPriors) {
            logVariantVsArtifactPrior = Math.log((variantCount + REGULARIZING_PSEUDOCOUNT) / (variantCount + technicalArtifactCount + REGULARIZING_PSEUDOCOUNT * 2));

            if (callableSites.isPresent()) {
                IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).forEach(n -> {
                    final double empiricalRatio = variantCountsByIndelLength.getOrDefault(n, new MutableDouble(0)).doubleValue() / callableSites.getAsDouble();
                    logVariantPriors.put(n, Math.log(Math.max(empiricalRatio, n == 0 ? 1.0e-8 : 1.0e-9)));
                });
            }
        }

        new IndexRange(0, clusters.size()).forEach(n -> {
            final double[] responsibilitiesForThisCluster = responsibilities.stream().mapToDouble(array -> array[n]).toArray();
            clusters.get(n).learn(data, responsibilitiesForThisCluster);
        });
    }

    /**
     * emission likelihood of given allele counts given that a variant is somatic, which is the weighted sum of cluster likelihoods
     */
    public double logLikelihoodGivenSomatic(final int totalCount, final int altCount) {
        final double[] logClusterLikelihoods = new IndexRange(0, clusters.size())
                .mapToDouble(c -> logClusterWeights[c] + clusters.get(c).logLikelihood(totalCount, altCount));
        return NaturalLogUtils.logSumExp(logClusterLikelihoods);
    }

    private double backgroundProbGivenSomatic(final int totalCount, final int altCount) {
        final double[] logClusterLikelihoods = new IndexRange(0, clusters.size())
                .mapToDouble(c -> logClusterWeights[c] + clusters.get(c).logLikelihood(totalCount, altCount));
        final double[] clusterProbabilities = NaturalLogUtils.normalizeFromLogToLinearSpace(logClusterLikelihoods);
        return clusterProbabilities[0];
    }

    public List<Pair<String, String>> clusteringMetadata() {
        final List<Pair<String, String>> result = new ArrayList<>();
        IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).forEach(n -> {
            final double logPrior = logVariantPriors.get(n);
            final String type = n == 0 ? "SNV" :
                    (n < 0 ? "deletion" : "insertion") + " of length " + Math.abs(n);
            result.add(ImmutablePair.of("Ln prior of " + type, Double.toString(logPrior)));
        });

        result.add(ImmutablePair.of("Background beta-binomial cluster",
                String.format("weight = %.4f, %s", Math.exp(logClusterWeights[0]), clusters.get(0).toString())));
        result.add(ImmutablePair.of("High-AF beta-binomial cluster",
                String.format("weight = %.4f, %s", Math.exp(logClusterWeights[1]), clusters.get(1).toString())));

        IntStream.range(2, clusters.size()).boxed()
                .sorted(Comparator.comparingDouble(c -> -logClusterWeights[c]))
                .forEach(c -> result.add(ImmutablePair.of("Binomial cluster",
                        String.format("weight = %.4f, %s", Math.exp(logClusterWeights[c]), clusters.get(c).toString()))));
        return result;
    }

    public static int indelLength(final VariantContext vc, final int altIndex) {
        return vc.getAlternateAllele(altIndex).length() - vc.getReference().length();
    }
}
