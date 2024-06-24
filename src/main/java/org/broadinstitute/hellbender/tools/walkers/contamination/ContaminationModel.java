package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.qc.Pileup;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This is the probabilistic contamination model that we use to distinguish homs from hets
 * The model is similar to that of ContEst, in that it assumes that each contaminant read is independently drawn from
 * the population.  However, it also accounts for allelic CNVs and doesn't just apply to hom alt sites.
 *
 * Because this model makes assumptions about the form of contamination and doesn't model CNVs precisely (although it
 * is robust to them in that they don't lead to bad solutions), we don't actually use it for the final contamination
 * estimate.  Instead, we learn the model, then use it to predict which sites are hom alts.  It's important to do this right
 * because contamination makes hom alts look more like hets, and allelic CNVs make hets look more like hom alts, hence
 * distinguishing correctly is not trivial.
 */
public class ContaminationModel {
    public static final double INITIAL_MAF_THRESHOLD = 0.40;
    public static final double MAF_TO_SWITCH_TO_HOM_REF = 0.25;
    public static final double MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF = 0.20;
    public static final double UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION = 0.15;
    public static final double UNSCRUPULOUS_HOM_REF_FRACTION_TO_REMOVE_FOR_POSSIBLE_LOH = 0.1;
    public static final double UNSCRUPULOUS_HOM_REF_PERCENTILE = 100 * ( 1 - UNSCRUPULOUS_HOM_REF_FRACTION_TO_REMOVE_FOR_POSSIBLE_LOH);
    public static final double MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD = 0.1;

    public static final double MAF_STEP_SIZE = 0.04;
    private static final double PRECISION_FOR_STANDARD_ERROR = 0.000001;
    private final double contamination;
    private final double errorRate;
    private final List<Double> minorAlleleFractions;
    private final List<List<PileupSummary>> segments;

    public static final int HOM_REF = 0;
    public static final int HOM_ALT = 3;
    private static final int NUM_ITERATIONS = 3;
    private static final double MIN_FRACTION_OF_SITES_TO_USE = 0.25;
    private static final double MIN_RELATIVE_ERROR = 0.2;
    private static final double MIN_ABSOLUTE_ERROR = 0.001;
    private static final List<Double> CONTAMINATION_INITIAL_GUESSES = Arrays.asList(0.02, 0.05, 0.1, 0.2);

    public ContaminationModel(List<PileupSummary> sites) {
        errorRate = calculateErrorRate(sites);

        // partition genome into minor allele fraction (MAF) segments to better distinguish hom alts from LoH hets.
        segments = ContaminationSegmenter.findSegments(sites);
        final int numSegments = segments.size();

        final List<Double> minorAlleleFractionsGuess = new ArrayList<>(Collections.nCopies(segments.size(), 0.5));
        final MutableDouble contaminationGuess = new MutableDouble(0);
        for (int n = 0; n < NUM_ITERATIONS; n++) {
            IntStream.range(0, numSegments).forEach(s -> minorAlleleFractionsGuess.set(s, calculateMinorAlleleFraction(contaminationGuess.doubleValue(), errorRate, segments.get(s))));
            final Pair<List<List<PileupSummary>>, List<Double>> nonLOHSegmentsAndMafs = getNonLOHSegments(segments, minorAlleleFractionsGuess);
            contaminationGuess.setValue(calculateContamination(errorRate, nonLOHSegmentsAndMafs.getLeft(), nonLOHSegmentsAndMafs.getRight()));
        }

        minorAlleleFractions = minorAlleleFractionsGuess;
        contamination = contaminationGuess.doubleValue();
    }

    private static Pair<List<List<PileupSummary>>, List<Double>> getNonLOHSegments(final List<List<PileupSummary>> segments, final List<Double> mafs) {
        final int numSites = segments.stream().mapToInt(List::size).sum();
        for (double minMaf = INITIAL_MAF_THRESHOLD; minMaf > 0; minMaf -= MAF_STEP_SIZE) {
            final double threshold = minMaf;
            final int[] nonLOHIndices = IntStream.range(0, segments.size()).filter(s -> mafs.get(s) > threshold).toArray();
            final List<List<PileupSummary>> nonLOHSegments = Arrays.stream(nonLOHIndices).mapToObj(segments::get).collect(Collectors.toList());
            final List<Double> nonLOHMafs = Arrays.stream(nonLOHIndices).mapToObj(mafs::get).collect(Collectors.toList());

            final int numNonLOHSites = nonLOHSegments.stream().mapToInt(List::size).sum();
            if ((double) numNonLOHSites / numSites > MIN_FRACTION_OF_SITES_TO_USE) {
                return ImmutablePair.of(nonLOHSegments, nonLOHMafs);
            }
        }
        return ImmutablePair.of(segments, mafs);
    }

    /**
     * Calculate the contamination in a tumor sample using the hom alt and hom ref sites inferred from this model, which
     * could be derived from the tumor itself or a matched normal.
     * @return
     */
    public Pair<Double, Double> calculateContaminationFromHoms(final List<PileupSummary> tumorSites) {
        for (double minMaf = INITIAL_MAF_THRESHOLD; minMaf >= 0; minMaf -= MAF_STEP_SIZE) {

            final Pair<Double, Double> result;
            if (minMaf > MAF_TO_SWITCH_TO_HOM_REF) {
                result = calculateContamination(Strategy.HOM_ALT, tumorSites, minMaf);
            } else if (minMaf > MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF) {
                result = calculateContamination(Strategy.HOM_REF, tumorSites, minMaf);
            } else {
                result = calculateContamination(Strategy.UNSCRUPULOUS_HOM_REF, tumorSites, minMaf);
            }
            if (!Double.isNaN(result.getLeft()) && result.getRight() < (result.getLeft() * MIN_RELATIVE_ERROR + MIN_ABSOLUTE_ERROR)) {
                return result;
            }
        }

        final Pair<Double, Double> result = calculateContamination(Strategy.UNSCRUPULOUS_HOM_REF, tumorSites, 0);
        return Double.isNaN(result.getLeft()) ? Pair.of(0.0, 1.0) : result;
    }

    /**
     * Depending on how many sites we have found for use in our estimate, we are more or less strict as to which calculation
     * to use.  In ideal cases (exomes and genomes) we have enough hom alt sites to base the calculation on ref contamination
     * in hom alt sites.  If this fails -- that is, if the hom alt calculation has too much uncertainty -- we calculate based
     * on the alt in hom ref signal, which is less reliable because it is more subject to error (or population-dependence) in
     * the population allele frequencies.  Finally, we may resort to a hom ref calculation that uses a heuristic instead of
     * principled genotyping based on the minor allele fraction segmentation.  This brings in additionally hom ref sites that
     * fall outside the span of the het segmentation and ought to be necessary only for very small gene panels.
     */
    private enum Strategy {
        HOM_ALT, HOM_REF, UNSCRUPULOUS_HOM_REF
    }


    private Pair<Double, Double> calculateContamination(final Strategy strategy, final List<PileupSummary> tumorSites, final double minMaf) {
        final boolean useHomAlt = strategy == Strategy.HOM_ALT;
        final List<PileupSummary> genotypingHoms;
        if (strategy == Strategy.HOM_ALT) {
            genotypingHoms = homAlts(minMaf);
        } else if (strategy == Strategy.HOM_REF) {
            genotypingHoms = homRefs(minMaf);
        } else {
            final List<PileupSummary> candidateHomRefs = tumorSites.stream()
                    .filter(site -> site.getAltFraction() < UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION)
                    .collect(Collectors.toList());
            final double altFractionThreshold = Math.max(MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD,
                    new Percentile(UNSCRUPULOUS_HOM_REF_PERCENTILE).evaluate(candidateHomRefs.stream().mapToDouble(PileupSummary::getAltFraction).toArray()));
            genotypingHoms = candidateHomRefs.stream().filter(site -> site.getAltFraction() <= altFractionThreshold).collect(Collectors.toList());
        }
        final List<PileupSummary> homs = subsetSites(tumorSites, genotypingHoms);
        final double tumorErrorRate = calculateErrorRate(tumorSites);

        // depth of ref in hom alt or alt in hom ref
        final ToIntFunction<PileupSummary> oppositeCount = useHomAlt ? PileupSummary::getRefCount : PileupSummary::getAltCount;
        final ToDoubleFunction<PileupSummary> oppositeAlleleFrequency = useHomAlt ? PileupSummary::getRefFrequency : PileupSummary::getAlleleFrequency;

        final long totalDepth = homs.stream().mapToLong(PileupSummary::getTotalCount).sum();

        // total read count of ref in hom alt or alt in hom ref, as the case may be
        final long oppositeDepth = homs.stream().mapToLong(oppositeCount::applyAsInt).sum();
        final long errorDepth = Math.round(totalDepth * tumorErrorRate / 3);
        final long contaminationOppositeDepth = Math.max(oppositeDepth - errorDepth, 0);

        final double totalDepthWeightedByOppositeFrequency = homs.stream()
                .mapToDouble(ps -> ps.getTotalCount() * oppositeAlleleFrequency.applyAsDouble(ps))
                .sum();

        final double contaminationEstimate = contaminationOppositeDepth / totalDepthWeightedByOppositeFrequency;

        // coeff1 and coeff2 are the sums of (1-f)*d and f*(1-f)*d^2, respectively, over all hom alts, where f is
        // population allele frequency and d is total read depth.  When using hom refs as a backup strategy we replace
        // 1-f with f in coeff1.
        final double coeff1 = homs.stream().mapToDouble(ps -> oppositeAlleleFrequency.applyAsDouble(ps) * ps.getTotalCount()).sum();
        final double coeff2 = homs.stream().mapToDouble(ps ->
                oppositeAlleleFrequency.applyAsDouble(ps)*(1 - oppositeAlleleFrequency.applyAsDouble(ps)) * MathUtils.square(ps.getTotalCount())
        ).sum();

        final DoubleUnaryOperator errorFunc = c -> homs.isEmpty() ? 1 : Math.sqrt(coeff1*c*(1-c) + coeff2*c*c) / totalDepthWeightedByOppositeFrequency;

        // {@code errorFunc} is a formula for the standard error (see docs/mutect/mutect.pdf) of the contamination estimate,
        // but it depends on the true (unknown) contamination, not the estimate.  For low contaminations simply plugging in
        // the estimated contamination to the formula does not work.  In particular, an estimate of 0 yields a standard error
        // of 0.  A more principled approach is to find a range of contaminations that is consistent with the error estimate.
        //  For example, suppose our error estimate is 0.0.  If the standard error as a function of contamination = 0.02 is 0.01,
        // then 0.02 is not consistent with an estimate of 0.0.  If the standard error as a function of contamination = 0.005
        // is 0.005, then it is consistent with our estimate of 0.0, and furthermore because {@code errorFunc} is monotonic we
        // know it is the largest true contamination that would be consistent.
        // More precisely, if our contamination estimate is e then the highest consistent contamination is c such that
        // c = e + errorFunc(c), which we find via binary search.  The lowest consistent contmination is found similarly,
        // and our improved standard error is max(highest - e, e - lowest).  The binary search is far from optimized (in fact we could solve this explicitly with a messy closed-form solution)
        // but it converges to a precision of 1e-6 in 20 iterations of the square root of a linear function.  This is fast enough.

        final OptionalDouble upperEstimate = MathUtils.binarySearchFindZero(c -> c - errorFunc.applyAsDouble(c) - contaminationEstimate,
                contaminationEstimate, 1.0, PRECISION_FOR_STANDARD_ERROR);
        final OptionalDouble lowerEstimate = MathUtils.binarySearchFindZero(c -> c + errorFunc.applyAsDouble(c) - contaminationEstimate,
                0.0, contaminationEstimate, PRECISION_FOR_STANDARD_ERROR);

        // if for some reason the binary search failed -- and if my math is to believed it should never fail -- we use
        // the standard error formula with the estimated contamination
        final double stdError = (upperEstimate.isPresent() && lowerEstimate.isPresent()) ?
                Math.max(upperEstimate.getAsDouble() - contaminationEstimate, contaminationEstimate - lowerEstimate.getAsDouble()) :
                errorFunc.applyAsDouble(contaminationEstimate);

        return Pair.of(Math.min(contaminationEstimate, 1.0), stdError);
    }

    private List<PileupSummary> getType(final int genotype, final double minMaf) {
        final int[] nonLOHIndices = IntStream.range(0, segments.size()).filter(s -> minorAlleleFractions.get(s) > minMaf).toArray();
        final List<List<PileupSummary>> nonLOHSegments = Arrays.stream(nonLOHIndices).mapToObj(segments::get).collect(Collectors.toList());
        final List<Double> nonLOHMafs = Arrays.stream(nonLOHIndices).mapToObj(minorAlleleFractions::get).collect(Collectors.toList());

        return IntStream.range(0, nonLOHSegments.size())
                .mapToObj(n -> nonLOHSegments.get(n).stream().filter(site -> probability(site, contamination, errorRate, nonLOHMafs.get(n), genotype) > 0.5))
                .flatMap(stream -> stream)
                .collect(Collectors.toList());
    }

    private List<PileupSummary> homAlts(final double minMaf) { return getType(HOM_ALT, minMaf); }
    private List<PileupSummary> homRefs(final double minMaf) { return getType(HOM_REF, minMaf); }

    public List<MinorAlleleFractionRecord> segmentationRecords() {
        return IntStream.range(0, segments.size()).mapToObj(n -> {
                    final List<PileupSummary> segment = segments.get(n);
                    final String contig = segment.get(0).getContig();
                    final int start = segment.get(0).getStart();
                    final int end = segment.get(segment.size() - 1).getEnd();
                    final double maf = minorAlleleFractions.get(n);
                    return new MinorAlleleFractionRecord(new SimpleInterval(contig, start, end), maf);
                }).collect(Collectors.toList());
    }

    private static double calculateErrorRate(final List<PileupSummary> sites) {
        final long totalBases = sites.stream().mapToInt(PileupSummary::getTotalCount).sum();
        final long otherAltBases = sites.stream().mapToInt(PileupSummary::getOtherAltCount).sum();
        return 1.5 * ((double) otherAltBases / totalBases);
    }

    private static double calculateMinorAlleleFraction(final double contamination, final double errorRate,  final List<PileupSummary> segment) {
        final DoubleUnaryOperator objective = maf -> segmentLogLikelihood(segment, contamination, errorRate, maf);
        return OptimizationUtils.max(objective, 0.1, 0.5, 0.4, 0.01, 0.01, 20).getPoint();
    }

    private static double calculateContamination(final double errorRate, final List<List<PileupSummary>> segments, final List<Double> mafs) {
        final DoubleUnaryOperator objective = c -> modelLogLikelihood(segments, c, errorRate, mafs);

        final List<UnivariatePointValuePair> optima = CONTAMINATION_INITIAL_GUESSES.stream()
                .map(initial -> OptimizationUtils.max(objective, 0, 0.5, initial, 1.0e-4, 1.0e-4, 30))
                .collect(Collectors.toList());

        return Collections.max(optima, Comparator.comparingDouble(UnivariatePointValuePair::getValue)).getPoint();
    }

    /**
     * Array of likelihoods over the genotypes hom ref, alt minor, alt major, hom alt
     * @param site
     * @param c the contamination
     * @param maf   the local minor allele fraction
     * @param errorRate the base error rate
     * @return a double[4], never null
     */
    private static double[] genotypeLikelihoods(final PileupSummary site, final double c, final double errorRate, double maf) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        // sample genotypes in order hom ref, alt minor, alt major, hom alt
        final double[] samplePriors = new double[] {(1 - f) * (1 - f), f * (1 - f), f * (1 - f), f * f};
        final double[] sampleAFs = new double[] {errorRate /3, maf, 1 - maf, 1 - errorRate};

        return new IndexRange(0, 4).mapToDouble(sg -> samplePriors[sg] *
                MathUtils.binomialProbability(n, k, (1 - c) * sampleAFs[sg] + c * f));
    }

    /**
     * Probability of a genotype hom ref, alt minor, alt major, hom alt
     * @param genotype 0 for hom ref, 3 for hom alt etc
     */
    private static double probability(final PileupSummary site, final double contamination, final double errorRate, final double minorAlleleFraction, final int genotype) {
        final double[] likelihoods = genotypeLikelihoods(site, contamination, errorRate, minorAlleleFraction);
        return likelihoods[genotype] / MathUtils.sum(likelihoods);
    }

    private static double segmentLogLikelihood(final List<PileupSummary> segment, final double contamination, final double errorRate, final double minorAlleleFraction) {
        return segment.stream().mapToDouble(site -> Math.log(MathUtils.sum(genotypeLikelihoods(site, contamination, errorRate, minorAlleleFraction)))).sum();
    }

    private static double modelLogLikelihood(final List<List<PileupSummary>> segments, final double contamination, final double errorRate, final List<Double> mafs) {
        Utils.validate(segments.size() == mafs.size(), " Must have one MAF per segment");
        return new IndexRange(0, segments.size()).sum(n -> segmentLogLikelihood(segments.get(n), contamination, errorRate, mafs.get(n)));
    }

    private static List<PileupSummary> subsetSites(final List<PileupSummary> sites, final List<PileupSummary> subsetLoci) {
        final OverlapDetector<PileupSummary> od = OverlapDetector.create(subsetLoci);
        return sites.stream().filter(od::overlapsAny).collect(Collectors.toList());
    }
}
