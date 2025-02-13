package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.utils.ValidationUtils;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingDebugger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Evidence-likelihoods container implementation based on integer indexed arrays.
 *
 * @param <A> the type of the allele the likelihood makes reference to.
 *
 * Note: this class uses FastUtil collections for speed.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class AlleleLikelihoods<EVIDENCE extends Locatable, A extends Allele> implements SampleList, AlleleList<A> {

    public static final double LOG_10_INFORMATIVE_THRESHOLD = 0.2;
    public static final double NATURAL_LOG_INFORMATIVE_THRESHOLD = MathUtils.log10ToLog(LOG_10_INFORMATIVE_THRESHOLD);
    private static final Logger logger = LogManager.getLogger(AlleleLikelihoods.class);

    private boolean isNaturalLog = false;
    private SimpleInterval subsettedGenomicLoc;
    private int filteredHaplotypeCount = 0;

    public boolean isNaturalLog() {
        return isNaturalLog;
    }
    private double getInformativeThreshold() {
        return isNaturalLog ? NATURAL_LOG_INFORMATIVE_THRESHOLD : LOG_10_INFORMATIVE_THRESHOLD;
    }

    /**
     * Index indicating that the reference allele is missing.
     */
    private static final int MISSING_INDEX = -1;

    /**
     * Evidence by sample index. Each sub array contains reference to the evidence of the ith sample.
     */
    private final List<List<EVIDENCE>> evidenceBySampleIndex;

    /**
     * Evidence disqualified by .
     */
    protected final List<List<EVIDENCE>> filteredEvidenceBySampleIndex;

    /**
     * Indexed per sample, allele and finally evidence (within sample).
     * <p>
     *     valuesBySampleIndex[s][a][r] == lnLk(R_r | A_a) where R_r comes from Sample s.
     * </p>
     */
    private final double[][][] valuesBySampleIndex;

    /**
     * Keeps track of the maximum number of evidences and likelihood values that can be stored
     * stored across all alleles.
     */
    private final int[] likelihoodsMatrixEvidenceCapacityBySampleIndex;

    /**
     * Holds the number of evidence per sample.
     */
    private final int[] numberOfEvidences;

    /**
     * Sample list.
     */
    protected final SampleList samples;

    /**
     * Allele list.
     */
    protected AlleleList<A> alleles;

    /**
     * Maps from each unit of evidence to its index within the sample.
     *
     * <p>In order to save CPU time the indices contained in this array (not the array itself) is
     * lazily initialized by invoking {@link #evidenceIndexBySampleIndex(int)}.</p>
     *
     * <p>Whenever a sample's evidence list is modified this cache must either be updated or invalidated by {@link #invalidateEvidenceToIndexCache(int)}.</p>
     */
    private final List<Object2IntMap<EVIDENCE>> evidenceIndexBySampleIndex;

    /**
     * Index of the reference allele if any, otherwise {@link #MISSING_INDEX}.
     */
    private int referenceAlleleIndex = MISSING_INDEX;

    /**
     * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
     */
    private final LikelihoodMatrix<EVIDENCE,A>[] sampleMatrices;

    /**
     * Constructs a new evidence-likelihood collection.
     *
     * <p>
     *     The initial likelihoods for all allele-evidence combinations are
     *     0.
     * </p>
     *
     * @param samples all supported samples in the collection.
     * @param alleles all supported alleles in the collection.
     * @param evidenceBySample evidence stratified per sample.
     *
     * @throws IllegalArgumentException if any of {@code allele}, {@code samples}
     * or {@code evidenceBySample} is {@code null},
     *  or if they contain null values.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    public AlleleLikelihoods(final SampleList samples,
                             final AlleleList<A> alleles,
                             final Map<String, List<EVIDENCE>> evidenceBySample) {
        Utils.nonNull(alleles);
        Utils.nonNull(samples);
        Utils.nonNull(evidenceBySample);

        this.samples = samples;
        this.alleles = alleles;

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        evidenceBySampleIndex = new ArrayList<>(sampleCount);
        valuesBySampleIndex = new double[sampleCount][][];
        likelihoodsMatrixEvidenceCapacityBySampleIndex = new int[sampleCount]; // set to 0s; what we need.
        referenceAlleleIndex = findReferenceAllele(alleles);
        numberOfEvidences = new int[sampleCount];

        evidenceIndexBySampleIndex = new ArrayList<>(Collections.nCopies(sampleCount, null));
        filteredEvidenceBySampleIndex = new ArrayList<>();
        samples().forEach(s -> this.filteredEvidenceBySampleIndex.add(new ArrayList<>(2)));

        setupIndexes(evidenceBySample, sampleCount, alleleCount);

        sampleMatrices = (LikelihoodMatrix<EVIDENCE, A>[]) new LikelihoodMatrix[sampleCount];
    }


    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    private AlleleLikelihoods(final AlleleList alleles,
                      final SampleList samples,
                      final List<List<EVIDENCE>> evidenceBySampleIndex,
                      final List<List<EVIDENCE>> filteredEvidenceBySampleIndex,
                      final double[][][] values) {
        this.samples = samples;
        this.alleles = alleles;
        this.evidenceBySampleIndex = evidenceBySampleIndex;
        this.valuesBySampleIndex = values;

        final int sampleCount = samples.numberOfSamples();

        this.evidenceIndexBySampleIndex = new ArrayList<>(Collections.nCopies(sampleCount, null));
        if (filteredEvidenceBySampleIndex != null) {
            this.filteredEvidenceBySampleIndex = filteredEvidenceBySampleIndex;
        } else {
            this.filteredEvidenceBySampleIndex = new ArrayList<>();
            samples().forEach(s -> this.filteredEvidenceBySampleIndex.add(new ArrayList<>(2)));
        }

        referenceAlleleIndex = findReferenceAllele(alleles);
        sampleMatrices = (LikelihoodMatrix<EVIDENCE,A>[]) new LikelihoodMatrix[sampleCount];
        numberOfEvidences = IntStream.range(0, sampleCount)
          .map(i -> evidenceBySampleIndex.get(i).size())
          .toArray();

        // We take the shortest allele's values array as the maximum evidence capacity for
        // each sample.
        likelihoodsMatrixEvidenceCapacityBySampleIndex = Arrays.stream(values)
                .mapToInt(sampleValues -> Arrays.stream(sampleValues)
                        .mapToInt(v -> v.length)
                        .min().orElse(0))
                .toArray();
    }

    /*
     * create an object using the private constructor. Unless absolutely required, please use public constructor
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    public static <EVIDENCE extends Locatable, A extends Allele> AlleleLikelihoods createAlleleLikelihoods(final AlleleList alleles,
                                            final SampleList samples,
                                            final List<List<EVIDENCE>> evidenceBySampleIndex,
                                            final List<List<EVIDENCE>> filteredEvidenceBySampleIndex,
                                            final double[][][] values) {
        return new AlleleLikelihoods<>(alleles, samples, evidenceBySampleIndex, filteredEvidenceBySampleIndex, values);
    }

    /**
     * Removes subset of alleles
     */
    public AlleleLikelihoods<EVIDENCE,A> removeAllelesToSubset(Collection<A> subsetOfAlleles){
        ValidationUtils.validateArg(alleles.asListOfAlleles().containsAll(subsetOfAlleles),
                () -> String.format("subsetOfAlleles must be a subset of the present alleles. Found atleast one allele that is not present: %s",
                        subsetOfAlleles.stream()
                                .filter(a -> !alleles.containsAllele(a))
                                .findFirst()
                                .map(A::toString)
                                .orElse("")));

        Map<A, List<A>> alleleMap = new HashMap<>();
        for (A al  : subsetOfAlleles){
            alleleMap.put(al, Arrays.asList(al));
        }

        return marginalize(alleleMap);
    }
    
    // Add all the indices to alleles, sample and evidence in the look-up maps.
    private void setupIndexes(final Map<String, List<EVIDENCE>> evidenceBySample, final int sampleCount, final int alleleCount) {
        for (int s = 0; s < sampleCount; s++) {
            final String sample = samples.getSample(s);
            final List<EVIDENCE> sampleEvidences = evidenceBySample.get(sample);
            numberOfEvidences[s] = sampleEvidences == null ? 0 : sampleEvidences.size();
            evidenceBySampleIndex.add(sampleEvidences == null ? new ArrayList<>() : new ArrayList<>(sampleEvidences));
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();

            final double[][] sampleValues = new double[alleleCount][sampleEvidenceCount];
            likelihoodsMatrixEvidenceCapacityBySampleIndex[s] = sampleEvidenceCount;
            valuesBySampleIndex[s] = sampleValues;
        }
    }

    // Search for the reference allele, if not found the index is {@link MISSING_REF}.
    private static int findReferenceAllele(final AlleleList<?> alleles) {
        return IntStream.range(0, alleles.numberOfAlleles()).filter(i -> alleles.getAllele(i).isReference()).findAny().orElse(MISSING_INDEX);
    }

    /**
     * Returns the index of a sample within the likelihood collection.
     *
     * @param sample the query sample.
     *
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    @Override
    public int indexOfSample(final String sample) {
        return samples.indexOfSample(sample);
    }

    /**
     * Number of samples included in the likelihood collection.
     * @return 0 or greater.
     */
    @Override
    public int numberOfSamples() {
        return samples.numberOfSamples();
    }

    /**
     * Returns sample name given its index.
     *
     * @param sampleIndex query index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is negative.
     *
     * @return never {@code null}.
     */
    @Override
    public String getSample(final int sampleIndex) {
        return samples.getSample(sampleIndex);
    }

    /**
     * Returns the index of an allele within the likelihood collection.
     *
     * @param allele the query allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    @Override
    public int indexOfAllele(final Allele allele) {
        return alleles.indexOfAllele(allele);
    }

    /**
     * Returns number of alleles in the collection.
     * @return 0 or greater.
     */
    @Override
    public int numberOfAlleles() {
        return alleles.numberOfAlleles();
    }

    /**
     * Returns the allele given its index.
     *
     * @param alleleIndex the allele index.
     *
     * @throws IllegalArgumentException the allele index is {@code null}.
     *
     * @return never {@code null}.
     */
    @Override
    public A getAllele(final int alleleIndex) {
        return alleles.getAllele(alleleIndex);
    }

    /**
     * Returns the units of evidence that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex the requested sample.
     * @return never {@code null} but perhaps a zero-length array if there is no evidence in sample. No element in
     *   the array will be null.
     */
    public List<EVIDENCE> sampleEvidence(final int sampleIndex) {
        return Collections.unmodifiableList(evidenceBySampleIndex.get(sampleIndex));
    }

    /**
     * Returns the units of evidence that have been removed by PairHMM error score filtering (and intentially not evidence
     * filtered by any other mechanism).
     *
     * @param sampleIndex the requested sample.
     * @return never {@code null} but perhaps a zero-length array if there is no filtered evidence for a sample. No element in
     *   the array will be null.
     */
    public List<EVIDENCE> filteredSampleEvidence(final int sampleIndex) {
        return Collections.unmodifiableList(filteredEvidenceBySampleIndex.get(sampleIndex));
    }

    /**
     * Returns an evidence vs allele likelihood matrix corresponding to a sample.
     *
     * @param sampleIndex target sample.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not null.
     *
     * @return never {@code null}
     */
    public LikelihoodMatrix<EVIDENCE,A> sampleMatrix(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final LikelihoodMatrix<EVIDENCE,A> extantResult = sampleMatrices[sampleIndex];
        if (extantResult == null) {
            return sampleMatrices[sampleIndex] = new SampleMatrix(sampleIndex);
        } else {
            return extantResult;
        }
    }

    public void switchToNaturalLog() {
        Utils.validate(!isNaturalLog, "Likelihoods have already been switched to natural log");
        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        for (int s = 0; s < sampleCount; s++) {
            final int evidenceCount = sampleEvidenceCount(s);
            for (int a = 0; a < alleleCount; a++) {
                for (int e = 0; e < evidenceCount; e++) {
                    valuesBySampleIndex[s][a][e] = MathUtils.log10ToLog(valuesBySampleIndex[s][a][e]);
                }
            }
        }
        isNaturalLog = true;
    }

    /**
     * Downsamples reads based on contamination fractions making sure that all alleles are affected proportionally.
     *
     * @param perSampleDownsamplingFraction contamination sample map where the sample name are the keys and the
     *                                       fractions are the values.
     *
     * @throws IllegalArgumentException if {@code perSampleDownsamplingFraction} is {@code null}.
     */
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        Utils.nonNull(perSampleDownsamplingFraction);

        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final double fraction = perSampleDownsamplingFraction.getOrDefault(samples.getSample(s), 0.0);
            if (Double.isNaN(fraction) || fraction <= 0.0) {
                continue;
            }

            final Map<A, List<EVIDENCE>> alleleEvidenceMap = evidenceByBestAlleleMap(s);
            final Collection<EVIDENCE> evidenceToRemove = fraction >= 1 ? evidenceBySampleIndex.get(s) :
                    AlleleBiasedDownsamplingUtils.selectAlleleBiasedEvidence(alleleEvidenceMap, fraction);
            removeEvidence(s, evidenceToRemove);
        }
    }

    /**
     * Adjusts likelihoods so that for each unit of evidence, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each unit of evidence based on the maximum alternative allele likelihood.
     *
     * @param maximumLikelihoodDifferenceCap maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    public void normalizeLikelihoods(final double maximumLikelihoodDifferenceCap, final boolean symmetricallyNormalizeAllelesToReference) {
        Utils.validateArg(maximumLikelihoodDifferenceCap < 0.0 && !Double.isNaN(maximumLikelihoodDifferenceCap),
                "the minimum reference likelihood fall must be negative");

        if (maximumLikelihoodDifferenceCap == Double.NEGATIVE_INFINITY) {
            return;
        }

        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0){ // trivial case there is no alleles.
            return;
        } else if (alleleCount == 1) {
            return;
        }

        for (int s = 0; s < valuesBySampleIndex.length; s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int evidenceCount = evidenceBySampleIndex.get(s).size();
            for (int r = 0; r < evidenceCount; r++) {
                normalizeLikelihoodsPerEvidence(maximumLikelihoodDifferenceCap, sampleValues, s, r, symmetricallyNormalizeAllelesToReference);
            }
        }
    }

    // Does the normalizeLikelihoods job for each piece of evidence.
    private void normalizeLikelihoodsPerEvidence(final double maximumBestAltLikelihoodDifference,
                                                 final double[][] sampleValues, final int sampleIndex, final int evidenceIndex, final boolean symmetricallyNormalizeAllelesToReference) {

        //allow the best allele to be the reference because asymmetry leads to strange artifacts like het calls with >90% alt reads
        final BestAllele bestAllele = searchBestAllele(sampleIndex,evidenceIndex,symmetricallyNormalizeAllelesToReference);

        final double worstLikelihoodCap = bestAllele.likelihood + maximumBestAltLikelihoodDifference;

        final int alleleCount = alleles.numberOfAlleles();

        // Guarantee to be the case by enclosing code.
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] < worstLikelihoodCap) {
                sampleValues[a][evidenceIndex] = worstLikelihoodCap;
            }
        }

    }

    /**
     * Returns the samples in this evidence-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable view on the evidence-likelihoods sample list.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<String> samples() {
        return samples.asListOfSamples();
    }

    /**
     * Returns the samples in this evidence-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable. It will not be updated if the collection
     *     allele list changes.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<A> alleles() {
        return alleles.asListOfAlleles();
    }


    /**
     * Search the best allele for a unit of evidence.
     *
     * @param sampleIndex including sample index.
     * @param evidenceIndex  target evidence index.
     *
     * @param priorities An array of allele priorities (higher values have higher priority) to be used, if present, to break ties for
     *                   uninformative likelihoods, in which case the evidence is assigned to the allele with the higher score.
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null}
     * if non-could be found.
     */
    private BestAllele searchBestAllele(final int sampleIndex, final int evidenceIndex, final boolean canBeReference, final double[] priorities) {
        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0 || (alleleCount == 1 && referenceAlleleIndex == 0 && !canBeReference)) {
            return new BestAllele(sampleIndex, evidenceIndex, MISSING_INDEX,
                    Double.NEGATIVE_INFINITY, MISSING_INDEX, Double.NEGATIVE_INFINITY);
        }

        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;

        int secondBestIndex = 0;
        double bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        double secondBestLikelihood = Double.NEGATIVE_INFINITY;

        for (int a = bestAlleleIndex + 1; a < alleleCount; a++) {
            if (!canBeReference && referenceAlleleIndex == a) {
                continue;
            }
            final double candidateLikelihood = sampleValues[a][evidenceIndex];
            if (candidateLikelihood > bestLikelihood) {
                secondBestIndex = bestAlleleIndex;
                bestAlleleIndex = a;
                secondBestLikelihood = bestLikelihood;
                bestLikelihood = candidateLikelihood;
            } else if (candidateLikelihood > secondBestLikelihood) {
                secondBestIndex = a;
                secondBestLikelihood = candidateLikelihood;
            }
        }

        if (priorities != null && bestLikelihood - secondBestLikelihood < getInformativeThreshold()) {
            double bestPriority = priorities[bestAlleleIndex];
            double secondBestPriority = priorities[secondBestIndex];
            for (int a = 0; a < alleleCount; a++) {
                final double candidateLikelihood = sampleValues[a][evidenceIndex];
                if (a == bestAlleleIndex || (!canBeReference && a == referenceAlleleIndex) || bestLikelihood - candidateLikelihood > getInformativeThreshold()) {
                    continue;
                }
                final double candidatePriority = priorities[a];

                if (candidatePriority > bestPriority) {
                    secondBestIndex = bestAlleleIndex;
                    bestAlleleIndex = a;
                    secondBestPriority = bestPriority;
                    bestPriority = candidatePriority;
                } else if (candidatePriority > secondBestPriority) {
                    secondBestIndex = a;
                    secondBestPriority = candidatePriority;
                }
            }
        }

        bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        secondBestLikelihood = secondBestIndex != bestAlleleIndex ? sampleValues[secondBestIndex][evidenceIndex] : Double.NEGATIVE_INFINITY;

        return new BestAllele(sampleIndex, evidenceIndex, bestAlleleIndex, bestLikelihood, secondBestIndex, secondBestLikelihood);
    }

    private BestAllele searchBestAllele(final int sampleIndex, final int evidenceIndex, final boolean canBeReference) {
        return searchBestAllele(sampleIndex, evidenceIndex, canBeReference, null);
    }

    public void changeEvidence(final Map<EVIDENCE, EVIDENCE> evidenceReplacements) {
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(s);
            final Object2IntMap<EVIDENCE> evidenceIndex = evidenceIndexBySampleIndex.get(s);
            final int sampleEvidenceCount = sampleEvidence.size();
            for (int r = 0; r < sampleEvidenceCount; r++) {
                final EVIDENCE evidence = sampleEvidence.get(r);
                final EVIDENCE replacement = evidenceReplacements.get(evidence);
                if (replacement == null) {
                    continue;
                }
                sampleEvidence.set(r, replacement);
                if (evidenceIndex != null) {
                    evidenceIndex.remove(evidence);
                    evidenceIndex.put(replacement, r);
                }
            }
        }
    }

    /**
     * Add alleles that are missing in the evidence-likelihoods collection giving all evidence a default
     * likelihood value.
     * @param candidateAlleles the potentially missing alleles.
     * @param defaultLikelihood the default evidence likelihood value for that allele.
     *
     * @return {@code true} iff the the evidence-likelihood collection was modified by the addition of the input alleles.
     *  So if all the alleles in the input collection were already present in the evidence-likelihood collection this method
     *  will return {@code false}.
     *
     * @throws IllegalArgumentException if {@code candidateAlleles} is {@code null} or there is more than
     * one missing allele that is a reference or there is one but the collection already has
     * a reference allele.
     */
    public boolean addMissingAlleles(final Collection<A> candidateAlleles, final double defaultLikelihood) {
        Utils.nonNull(candidateAlleles, "the candidateAlleles list cannot be null");
        if (candidateAlleles.isEmpty()) {
            return false;
        }
        final List<A> allelesToAdd = candidateAlleles.stream().filter(allele -> !alleles.containsAllele(allele)).collect(Collectors.toList());

        if (allelesToAdd.isEmpty()) {
            return false;
        }

        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = alleles.numberOfAlleles() + allelesToAdd.size();

        int referenceIndex = this.referenceAlleleIndex;

        @SuppressWarnings("unchecked")
        final List<A> newAlleles = ListUtils.union(alleles.asListOfAlleles(), allelesToAdd);
        alleles = new IndexedAlleleList<>(newAlleles);

        // if we previously had no reference allele, update the reference index if a reference allele is added
        // if we previously had a reference and try to add another, throw an exception
        final OptionalInt indexOfReferenceInAllelesToAdd = IntStream.range(0, allelesToAdd.size())
                .filter(n -> allelesToAdd.get(n).isReference()).findFirst();
        if (referenceIndex != MISSING_INDEX) {
            Utils.validateArg(!indexOfReferenceInAllelesToAdd.isPresent(), "there can only be one reference allele");
        } else if (indexOfReferenceInAllelesToAdd.isPresent()){
            referenceAlleleIndex = oldAlleleCount + indexOfReferenceInAllelesToAdd.getAsInt();
        }

        //copy old allele likelihoods and set new allele likelihoods to the default value
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();
            final double[][] newValuesBySampleIndex = Arrays.copyOf(valuesBySampleIndex[s], newAlleleCount);
            for (int a = oldAlleleCount; a < newAlleleCount; a++) {
                newValuesBySampleIndex[a] = new double[likelihoodsMatrixEvidenceCapacityBySampleIndex[s]];
                if (defaultLikelihood != 0.0) {
                    Arrays.fill(newValuesBySampleIndex[a], 0, sampleEvidenceCount, defaultLikelihood);
                }
                // Fill the rest with NaNs
                Arrays.fill(newValuesBySampleIndex[a], sampleEvidenceCount, newValuesBySampleIndex[a].length, Double.NaN);
            }
            valuesBySampleIndex[s] = newValuesBySampleIndex;
        }
        return true;
    }

    /**
     * Group evidence into lists of evidence -- for example group by read name to force read pairs to support a single haplotype.
     *
     * Log Likelihoods are summed over all evidence
     * in a group, corresponding to an independent evidence assumption.  Since this container's likelihoods generally pertain to
     * sequencing only (and not sample prep etc) this is usually a good assumption.
     *
     * @param groupingFunction Attribute function for grouping evidence, for example GATKRead::getName
     * @param gather Transformation applied to collections of evidence with same value of groupingFunction.  For example, Fragment::new
     *               to construct a fragment out of a pair of reads with the same name
     *
     * @return a new AlleleLikelihoods based on the grouped, transformed evidence.
     */
    public <U, NEW_EVIDENCE_TYPE extends Locatable> AlleleLikelihoods<NEW_EVIDENCE_TYPE, A> groupEvidence(final Function<EVIDENCE, U> groupingFunction, final Function<List<EVIDENCE>, NEW_EVIDENCE_TYPE> gather) {
        final int sampleCount = samples.numberOfSamples();
        final double[][][] newLikelihoodValues = new double[sampleCount][][];
        final int alleleCount = alleles.numberOfAlleles();

        final List<List<NEW_EVIDENCE_TYPE>> newEvidenceBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {


            final List<List<EVIDENCE>> evidenceGroups = new ArrayList<>(sampleEvidence(s).stream()
                    .collect(Collectors.groupingBy(groupingFunction)).values());


            final int newEvidenceCount = evidenceGroups.size();

            final double[][] oldSampleValues = valuesBySampleIndex[s];
            newLikelihoodValues[s] = new double[alleleCount][newEvidenceCount];

            // For each old allele and read we update the new table keeping the maximum likelihood.
            for (int newEvidenceIndex = 0; newEvidenceIndex < newEvidenceCount; newEvidenceIndex++) {
                for (int a = 0; a < alleleCount; a++) {
                    for (final EVIDENCE evidence : evidenceGroups.get(newEvidenceIndex)) {
                        final int oldEvidenceIndex = evidenceIndex(s, evidence);
                        newLikelihoodValues[s][a][newEvidenceIndex] += oldSampleValues[a][oldEvidenceIndex];
                    }
                }
            }
            newEvidenceBySampleIndex.add(evidenceGroups.stream().map(gather).collect(Collectors.toList()));
        }

        // Finally we create the new read-likelihood
        final AlleleLikelihoods<NEW_EVIDENCE_TYPE, A> result = new AlleleLikelihoods<>(
                alleles,
                samples,
                newEvidenceBySampleIndex,
                null, //TODO this is only currently used for Somatic and i'm alright with removing this for now but this is NOT robust and 3 of these methods is too many
                newLikelihoodValues);
        result.isNaturalLog = this.isNaturalLog;
        return result;
    }

    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each evidence in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and evidence as the original.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this evidence-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <B extends Allele> AlleleLikelihoods<EVIDENCE, B> marginalize(final Map<B, List<A>> newToOldAlleleMap) {
        Utils.nonNull(newToOldAlleleMap);

        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final List<BitSet> newToOldAlleleIndex = this.newToOldAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        // We calculate the marginal likelihoods.
        final double[][][] newLikelihoodValues = marginalLikelihoodsDirect( newAlleleCount, newToOldAlleleIndex);

        final int sampleCount = samples.numberOfSamples();

        final List<List<EVIDENCE>> newEvidenceBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {
            newEvidenceBySampleIndex.add(new ArrayList<>(evidenceBySampleIndex.get(s)));
        }

        // Finally we create the new evidence-likelihood
        final AlleleLikelihoods<EVIDENCE, B> result = new AlleleLikelihoods<>(
                new IndexedAlleleList(newAlleles),
                samples,
                newEvidenceBySampleIndex,
                filteredEvidenceBySampleIndex,
                newLikelihoodValues);
        result.isNaturalLog = isNaturalLog;
        return result;
    }

    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    private double[][][] marginalLikelihoods(final int oldAlleleCount, final int newAlleleCount,
                                             final int[] oldToNewAlleleIndexMap) {
        final int sampleCount = samples.numberOfSamples();
        final double[][][] result = new double[sampleCount][][];

        for (int s = 0; s < sampleCount; s++) {
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();
            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final double[][] newSampleValues = result[s] = new double[newAlleleCount][sampleEvidenceCount];
            // We initiate all likelihoods to -Inf.
            for (int a = 0; a < newAlleleCount; a++) {
                Arrays.fill(newSampleValues[a], Double.NEGATIVE_INFINITY);
            }
            // For each old allele and read we update the new table keeping the maximum likelihood.
            for (int r = 0; r < sampleEvidenceCount; r++) {
                for (int a = 0; a < oldAlleleCount; a++) {
                    final int newAlleleIndex = oldToNewAlleleIndexMap[a];
                    if (newAlleleIndex == MISSING_INDEX) {
                        continue;
                    }
                    final double likelihood = oldSampleValues[a][r];
                    if (likelihood > newSampleValues[newAlleleIndex][r]) {
                        newSampleValues[newAlleleIndex][r] = likelihood;
                    }
                }
            }
        }
        return result;
    }

    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    private double[][][] marginalLikelihoodsDirect(final int newAlleleCount, final List<BitSet> newToOldAlleleIndexMap) {

        final int sampleCount = samples.numberOfSamples();
        final double[][][] result = new double[sampleCount][][];

        for (int s = 0; s < sampleCount; s++) {
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();
            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final double[][] newSampleValues = result[s] = new double[newAlleleCount][sampleEvidenceCount];

            // For each old allele and unit of evidence we update the new table keeping the maximum likelihood.
            for (int newAllele = 0; newAllele < newAlleleCount; newAllele++) {
                final BitSet oldAlleleSet = newToOldAlleleIndexMap.get(newAllele);
                for (int r = 0; r < sampleEvidenceCount; r++) {
                    final int finalR = r;
                    newSampleValues[newAllele][finalR] = oldAlleleSet.stream()
                            .mapToDouble(oldA -> oldSampleValues[oldA][finalR])
                            .max()
                            .orElse(Double.NEGATIVE_INFINITY);
                }
            }
        }
        return result;
    }


    // calculates an old to new allele index map array.
    private <B extends Allele> int[] oldToNewAlleleIndexMap(final Map<B, List<A>> newToOldAlleleMap, final int oldAlleleCount, final B[] newAlleles) {
        Arrays.stream(newAlleles).forEach(Utils::nonNull);
        Utils.containsNoNull(newToOldAlleleMap.values(), "no new allele list can be null");
        newToOldAlleleMap.values().forEach(oldList -> Utils.containsNoNull(oldList,"old alleles cannot be null"));

        final int[] oldToNewAlleleIndexMap = new int[oldAlleleCount];
        Arrays.fill(oldToNewAlleleIndexMap, -1);  // -1 indicate that there is no new allele that make reference to that old one.

        for (int newIndex = 0; newIndex < newAlleles.length; newIndex++) {
            final B newAllele = newAlleles[newIndex];
            for (final A oldAllele : newToOldAlleleMap.get(newAllele)) {
                final int oldAlleleIndex = indexOfAllele(oldAllele);
                if (oldAlleleIndex == -1) {
                    throw new IllegalArgumentException("missing old allele " + oldAllele + " in likelihood collection ");
                }
                if (oldToNewAlleleIndexMap[oldAlleleIndex] != -1) {
                    // DK: PATCH: uncollapsed haplotypes might result in identicals. ignore this for now
                    //throw new IllegalArgumentException("collision: two new alleles make reference to the same old allele");
                    logger.debug("ignored! - collision: two new alleles make reference to the same old allele");
                    continue;
                }
                oldToNewAlleleIndexMap[oldAlleleIndex] = newIndex;
            }
        }
        return oldToNewAlleleIndexMap;
    }


    // calculates an old to new allele index map array.
    private <B extends Allele> List<BitSet> newToOldAlleleIndexMap(final Map<B, List<A>> newToOldAlleleMap, final int oldAlleleCount, final B[] newAlleles) {
        Arrays.stream(newAlleles).forEach(Utils::nonNull);
        Utils.containsNoNull(newToOldAlleleMap.values(), "no new allele list can be null");
        newToOldAlleleMap.values().forEach(oldList -> Utils.containsNoNull(oldList,"old alleles cannot be null"));

        final List<BitSet>  newToOldAlleleIndexMap = new ArrayList<>(newAlleles.length);
        IntStream.range(0, newAlleles.length).forEach(i -> newToOldAlleleIndexMap.add(i, new BitSet(oldAlleleCount)));

        for (int newIndex = 0; newIndex < newAlleles.length; newIndex++) {
            for (final A oldAllele : newToOldAlleleMap.get(newAlleles[newIndex])) {
                final int oldAlleleIndex = indexOfAllele(oldAllele);
                if (oldAlleleIndex == -1) {
                    throw new IllegalArgumentException("missing old allele " + oldAllele + " in likelihood collection ");
                }
                newToOldAlleleIndexMap.get(newIndex).set(oldAlleleIndex);
            }
        }
        return newToOldAlleleIndexMap;
    }


    /**
     * Add more evidence to the collection.
     *
     * @param evidenceBySample evidence to add.
     * @param initialLikelihood the likelihood for the new entries.
     *
     * @throws IllegalArgumentException if {@code evidenceBySample} is {@code null} or {@code evidenceBySample} contains
     *  {@code null} evidence, or {@code evidenceBySample} contains evidence already present in the evidence-likelihood
     *  collection.
     */
    public void addEvidence(final Map<String,List<EVIDENCE>> evidenceBySample, final double initialLikelihood) {
        for (final Map.Entry<String,List<EVIDENCE>> entry : evidenceBySample.entrySet()) {
            final String sample = entry.getKey();
            final List<EVIDENCE> newSampleEvidence = entry.getValue();
            final int sampleIndex = samples.indexOfSample(sample);

            if (sampleIndex == MISSING_INDEX) {
                throw new IllegalArgumentException("input sample " + sample + " is not part of the evidence-likelihoods collection");
            }

            if (newSampleEvidence == null || newSampleEvidence.isEmpty()) {
                continue;
            }

            final int oldEvidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
            appendEvidence(newSampleEvidence, sampleIndex);
            final int newEvidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
            extendsLikelihoodArrays(initialLikelihood, sampleIndex, oldEvidenceCount, newEvidenceCount);
        }
    }

    // Extends the likelihood arrays-matrices.
    private void extendsLikelihoodArrays(final double initialLikelihood, final int sampleIndex, final int oldEvidenceCount, final int newEvidenceCount) {
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        final int numberOfAlleles = alleles.numberOfAlleles();
        ensureLikelihoodsMatrixEvidenceCapacity(sampleIndex, newEvidenceCount, sampleValues, numberOfAlleles);
        for (int a = 0; a < numberOfAlleles; a++) {
            Arrays.fill(sampleValues[a], oldEvidenceCount, newEvidenceCount, initialLikelihood);
        }
    }

    // Resizes the lk value holding arrays to be able to handle at least "x" amount of evidence.
    private void ensureLikelihoodsMatrixEvidenceCapacity(final int sampleIndex, final int x,
                                                         final double[][] sampleValues, final int numberOfAlleles) {
        final int currentCapacity = likelihoodsMatrixEvidenceCapacityBySampleIndex[sampleIndex];
        if (currentCapacity < x) {
            final int newCapacity = Math.max(currentCapacity, x) << 1; // we double it to avoid repetitive 1-element extensions resizing.
            for (int a  = 0; a < numberOfAlleles; a++) {
                sampleValues[a] = Arrays.copyOf(sampleValues[a], newCapacity);
                // We don't need the following fill for this to work as intended but
                // with NaN we ensure some clear failure output in case there is a bug that
                // uses "zombie" likelihoods (e.g. github issue {@link https://github.com/broadinstitute/gatk/pull/7153 #7153).
                Arrays.fill(sampleValues[a], currentCapacity, newCapacity, Double.NaN);
            }
            likelihoodsMatrixEvidenceCapacityBySampleIndex[sampleIndex] = newCapacity;
        }
    }

    // Append the new evidence reference into the structure per-sample, returning the count of evidence actually added (duplicates are not added)
    // NOTE: the evidence-to-index cache is updated in place and not invalidated via {@link #invalidateEvidenceToIndexCache(int)} because adding new evidence
    // to the cache, as opposed to removing evidence, is just a matter of appending entries
    private void appendEvidence(final List<EVIDENCE> newSampleEvidence, final int sampleIndex) {

        final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(sampleIndex);
        final Object2IntMap<EVIDENCE> sampleEvidenceIndex = evidenceIndexBySampleIndex(sampleIndex);

        for (final EVIDENCE newEvidence : newSampleEvidence) {
            final int previousValue = sampleEvidenceIndex.put(newEvidence, sampleEvidence.size());
            if (previousValue == MISSING_INDEX) {
                sampleEvidence.add(newEvidence);
            } else {
                sampleEvidenceIndex.put(newEvidence, previousValue); // revert
            }
        }

        numberOfEvidences[sampleIndex] = sampleEvidence.size();
    }

    /**
     * Adds the non-reference allele to the evidence-likelihood collection setting each evidence likelihood to the second
     * best found (or best one if only one allele has likelihood).
     *
     * <p>Nothing will happen if the evidence-likelihoods collection already includes the non-ref allele</p>
     *
     * <p>
     *     <i>Implementation note: even when strictly speaking we do not need to demand the calling code to pass
     *     the reference the non-ref allele, we still demand it in order to lead the
     *     the calling code to use the right generic type for this likelihoods
     *     collection {@link Allele}.</i>
     * </p>
     *
     * @param nonRefAllele the non-ref allele.
     *
     * @throws IllegalArgumentException if {@code nonRefAllele} is anything but the designated &lt;NON_REF&gt;
     * symbolic allele {@link Allele#NON_REF_ALLELE}.
     */
    public void addNonReferenceAllele(final A nonRefAllele) {
        Utils.nonNull(nonRefAllele, "non-ref allele cannot be null");
        if (!nonRefAllele.equals(Allele.NON_REF_ALLELE)) {
            throw new IllegalArgumentException("the non-ref allele is not valid");
        } else if (!alleles.containsAllele(nonRefAllele)) {
            if (addMissingAlleles(Collections.singleton(nonRefAllele), Double.NEGATIVE_INFINITY)) {
                updateNonRefAlleleLikelihoods();
            }
        }
    }

    /**
     * Updates the likelihoods of the non-ref allele, if present, considering all non-symbolic alleles avaialble.
     */
    public void updateNonRefAlleleLikelihoods() {
        updateNonRefAlleleLikelihoods(alleles);
    }

    /**
     * Updates the likelihood of the NonRef allele (if present) based on the likelihoods of a set of non-symbolic
     */
    public void updateNonRefAlleleLikelihoods(final AlleleList<A> allelesToConsider) {
        final int nonRefAlleleIndex = indexOfAllele(Allele.NON_REF_ALLELE);
        if ( nonRefAlleleIndex < 0) {
            return;
        }
        final int alleleCount = alleles.numberOfAlleles();
        final int nonSymbolicAlleleCount = alleleCount - 1;
        // likelihood buffer reused across evidence:
        final double[] qualifiedAlleleLikelihoods = new double[nonSymbolicAlleleCount];
        final Median medianCalculator = new Median();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int evidenceCount = evidenceBySampleIndex.get(s).size();
            for (int r = 0; r < evidenceCount; r++) {
                final BestAllele bestAllele = searchBestAllele(s, r, true);
                int numberOfQualifiedAlleleLikelihoods = 0;
                for (int i = 0; i < alleleCount; i++) {
                    final double alleleLikelihood = sampleValues[i][r];
                    if (i != nonRefAlleleIndex && alleleLikelihood < bestAllele.likelihood
                            && !Double.isNaN(alleleLikelihood) && allelesToConsider.indexOfAllele(alleles.getAllele(i)) != MISSING_INDEX) {
                        qualifiedAlleleLikelihoods[numberOfQualifiedAlleleLikelihoods++] = alleleLikelihood;
                    }
                }
                final double nonRefLikelihood = medianCalculator.evaluate(qualifiedAlleleLikelihoods, 0, numberOfQualifiedAlleleLikelihoods);
                // when the median is NaN that means that all applicable likekihoods are the same as the best
                // so the evidence is not informative at all given the existing alleles. Unless there is only one (or zero) concrete
                // alleles with give the same (the best) likelihood to the NON-REF. When there is only one (or zero) concrete
                // alleles we set the NON-REF likelihood to NaN.
                sampleValues[nonRefAlleleIndex][r] = !Double.isNaN(nonRefLikelihood) ? nonRefLikelihood
                        : nonSymbolicAlleleCount <= 1 ? Double.NaN : bestAllele.likelihood;
            }
        }
    }

    /**
     * Returns the collection of best allele estimates for the evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit fo evidence in the evidence-likelihoods collection.
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final Function<A, Double> tieBreakingPriority) {
        return IntStream.range(0, numberOfSamples()).boxed().flatMap(n -> bestAllelesBreakingTies(n, tieBreakingPriority).stream()).collect(Collectors.toList());
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    public Collection<BestAllele> bestAllelesBreakingTies() {
        return IntStream.range(0, numberOfSamples()).boxed().flatMap(n -> bestAllelesBreakingTies(n).stream()).collect(Collectors.toList());
    }

    /**
     * Returns the collection of best allele estimates for one sample's evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit of evidence in the evidence-likelihoods collection.
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final String sample, final Function<A, Double> tieBreakingPriority) {
        final int sampleIndex = indexOfSample(sample);
        return bestAllelesBreakingTies(sampleIndex, tieBreakingPriority);
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final String sample) {
        final int sampleIndex = indexOfSample(sample);
        return bestAllelesBreakingTies(sampleIndex);
    }

    /**
     * Returns the collection of best allele estimates for one sample's evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit of evidence in the evidence-likelihoods collection.
     */
    private Collection<BestAllele> bestAllelesBreakingTies(final int sampleIndex, final Function<A, Double> tieBreakingPriority) {
        Utils.validIndex(sampleIndex, numberOfSamples());

        //TODO: this currently just does ref vs alt.  Really we want CIGAR complexity.
        final double[] priorities = alleles == null ? null :
                alleles.asListOfAlleles().stream().mapToDouble(tieBreakingPriority::apply).toArray();

        final int evidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
        final List<BestAllele> result = new ArrayList<>(evidenceCount);
        for (int r = 0; r < evidenceCount; r++) {
            result.add(searchBestAllele(sampleIndex, r, true, priorities));
        }

        return result;
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    private Collection<BestAllele> bestAllelesBreakingTies(final int sampleIndex) {
        return bestAllelesBreakingTies(sampleIndex, a -> a.isReference() ? 1.0 : 0);
    }

    /**
     * Returns evidence stratified by best allele.
     * @param sampleIndex the target sample.
     * @return never {@code null}, perhaps empty.
     */
    private Map<A,List<EVIDENCE>> evidenceByBestAlleleMap(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        final int alleleCount = alleles.numberOfAlleles();
        final int sampleEvidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
        final Map<A,List<EVIDENCE>> result = new LinkedHashMap<>(alleleCount);
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(sampleEvidenceCount));
        }
        evidenceByBestAlleleMap(sampleIndex,result);
        return result;
    }

    /**
     * Returns evidence stratified by best allele.
     * @return never {@code null}, perhaps empty.
     */
    @VisibleForTesting
    Map<A,List<EVIDENCE>> evidenceByBestAlleleMap() {
        final int alleleCount = alleles.numberOfAlleles();
        final Map<A,List<EVIDENCE>> result = new LinkedHashMap<>(alleleCount);
        final int totalEvidenceCount = evidenceCount();
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(totalEvidenceCount));
        }
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            evidenceByBestAlleleMap(s, result);
        }
        return result;
    }

    private void evidenceByBestAlleleMap(final int sampleIndex, final Map<A, List<EVIDENCE>> result) {
        final int evidenceCount = evidenceBySampleIndex.get(sampleIndex).size();

        for (int r = 0; r < evidenceCount; r++) {
            final BestAllele bestAllele = searchBestAllele(sampleIndex,r,true);
            if (!bestAllele.isInformative()) {
                continue;
            }
            result.get(bestAllele.allele).add(bestAllele.evidence);
        }
    }

    /**
     * Returns the index of evidence within a sample evidence-likelihood sub collection.
     * @param sampleIndex the sample index.
     * @param evidence the query evidence.
     * @return -1 if there is no such evidence in that sample, 0 or greater otherwise.
     */
    @VisibleForTesting
    int evidenceIndex(final int sampleIndex, final EVIDENCE evidence) {
        final Object2IntMap<EVIDENCE> index = evidenceIndexBySampleIndex(sampleIndex);
        return index.getInt(evidence);
    }

    /**
     * Returns the total count of evidence in the evidence-likelihood collection.
     */
    public int evidenceCount() {
        return evidenceBySampleIndex.stream().mapToInt(List::size).sum();
    }

    /**
     * Returns the quantity of evidence that belongs to a sample in the evidence-likelihood collection.
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid sample index.
     * @return 0 or greater.
     */
    public int sampleEvidenceCount(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        return evidenceBySampleIndex.get(sampleIndex).size();
    }


    /**
     * Remove those reads that do not comply with a requirement.
     *
     * @param predicate the predicate representing the requirement.
     *
     * <p>
     *     This method modifies the current read-likelihoods collection.
     * </p>
     * <p>
     *     Any exception thrown by the predicate will be propagated to the calling code.
     * </p>
     *
     * @throws IllegalArgumentException if {@code predicate} is {@code null}.
     */
    public void retainEvidence(final Predicate<? super EVIDENCE> predicate) {
        Utils.nonNull(predicate);
        final int sampleCount = samples.numberOfSamples();

        for (int s = 0; s < sampleCount; s++) {
            // Remove evidence from the primary data
            final List<EVIDENCE> sampleEvidence = this.evidenceBySampleIndex.get(s);
            final int[] removeIndices = IntStream.range(0, sampleEvidence.size())
                    .filter(i -> !predicate.test(sampleEvidence.get(i)))
                    .toArray();
            removeEvidenceByIndex(s, removeIndices);

            // If applicable also apply the predicate to the filters
            final List<EVIDENCE> sampleFiltered = filteredEvidenceBySampleIndex.get(s).stream()
                    .filter(predicate).collect(Collectors.toList());
            filteredEvidenceBySampleIndex.set(s, sampleFiltered);
        }
    }

    protected double maximumLikelihoodOverAllAlleles(final int sampleIndex, final int evidenceIndex) {
        double result = Double.NEGATIVE_INFINITY;
        final int alleleCount = alleles.numberOfAlleles();
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] > result) {
                result = sampleValues[a][evidenceIndex];
            }
        }
        return result;
    }

    public void setVariantCallingSubsetUsed(final SimpleInterval loc) {
        this.subsettedGenomicLoc = loc;
    }

    /**
     * Returns the location used for subsetting. May be null.
     */
    public SimpleInterval getVariantCallingSubsetApplied() {
        return subsettedGenomicLoc;
    }

    /**
     * Replaces the alleles in the readLikelihood matrix. Relevant for the uncollapsing code
     *
     * @param newAlleles
     */
    public void changeAlleles(List<A> newAlleles) {
        alleles = new IndexedAlleleList<>(newAlleles);
    }

    public int getFilteredHaplotypeCount() {
        return filteredHaplotypeCount;
    }
    public void setFilteredHaplotypeCount(int filteredHaplotypeCount) {
        this.filteredHaplotypeCount = filteredHaplotypeCount;
    }

    /**
     * Contains information about the best allele for a unit of evidence.
     */
    public final class BestAllele {

        /**
         * Null if there is no possible match (no allele?).
         */
        public final A allele;
        public final A second_best_allele;
        /**
         * The containing sample.
         */
        public final String sample;

        /**
         * The query evidence.
         */
        public final EVIDENCE evidence;

        /**
         * If allele != null, the indicates the likelihood of the evidence.
         */
        public final double likelihood;
        public final double secondBestLikelihood;

        /**
         * Confidence that the evidence actually was generated under that likelihood.
         * This is equal to the difference between this and the second best allele match.
         */
        public final double confidence;

        private BestAllele(final int sampleIndex, final int evidenceIndex, final int bestAlleleIndex,
                           final double likelihood, final int secondBestAlleleIndex, final double secondBestLikelihood) {
            allele = bestAlleleIndex == MISSING_INDEX ? null : alleles.getAllele(bestAlleleIndex);
            second_best_allele = secondBestAlleleIndex == -1 ? null : alleles.getAllele(secondBestAlleleIndex);
            this.likelihood = likelihood;
            sample = samples.getSample(sampleIndex);
            evidence = evidenceBySampleIndex.get(sampleIndex).get(evidenceIndex);
            confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
            this.secondBestLikelihood = secondBestLikelihood;
        }

        public boolean isInformative() {
            return confidence > getInformativeThreshold();
        }
    }

    private void removeEvidence(final int sampleIndex, final Collection<EVIDENCE> evidences) {
        final Object2IntMap<EVIDENCE> evidenceIndexes = evidenceIndexBySampleIndex(sampleIndex);
        final int[] indexesToRemove = evidences.stream().mapToInt(e -> {
            final int index = evidenceIndexes.getInt(e);
            if (index == MISSING_INDEX) {
                throw new IllegalArgumentException("evidence provided is not in sample");
            } else {
                return index;
            }
        }).sorted().distinct().toArray();
        removeEvidenceByIndex(sampleIndex, indexesToRemove);
    }

    // remove evidence and unset the {@code evidenceIndexBySampleIndex} cache for this sample
    // assumes that evidencesToRemove is sorted and without duplicates.
    private void removeEvidenceByIndex(final int sampleIndex, final int[] evidencesToRemove) {
        final int numToRemove = evidencesToRemove.length;
        if (numToRemove == 0) {
            return;
        }
        final int oldEvidenceCount = numberOfEvidences[sampleIndex];
        final int newEvidenceCount = oldEvidenceCount - numToRemove;

        // update the list of evidence and evidence count
        final List<EVIDENCE> oldEvidence = evidenceBySampleIndex.get(sampleIndex);
        final List<EVIDENCE> newEvidence = new ArrayList<>(newEvidenceCount);
        for (int n = 0, numRemoved = 0; n < oldEvidenceCount; n++) {
            if (numRemoved < numToRemove && n == evidencesToRemove[numRemoved]) {
                numRemoved++;
            } else {
                newEvidence.add(oldEvidence.get(n));

                // update the likelihoods arrays in place
                for (final double[] alleleValues : valuesBySampleIndex[sampleIndex]) {
                    alleleValues[n - numRemoved] = alleleValues[n];
                }
            }
        }
        // set to NaN lks of the deleted positions in lk value arrays.
        for (final double[] alleleValues : valuesBySampleIndex[sampleIndex]) {
            Arrays.fill(alleleValues, newEvidenceCount, alleleValues.length, Double.NaN);
        }
        evidenceBySampleIndex.set(sampleIndex, newEvidence);
        numberOfEvidences[sampleIndex] = newEvidenceCount;
        invalidateEvidenceToIndexCache(sampleIndex);
    }

    // The evidenceToIndex map becomes invalid when the evidence list is modified, for example by deleting evidence
    // When adding evidence it is simple enough to add new entries to the map, but we must be careful to do so.
    private void invalidateEvidenceToIndexCache(final int sampleIndex) {
        evidenceIndexBySampleIndex.set(sampleIndex, null);
    }

    /**
     * Removes those read that the best possible likelihood given any allele is just too low.
     *
     * <p>
     *     This is determined by a maximum error per read-base against the best likelihood possible.
     * </p>
     *
     * @param log10MinTrueLikelihood Function that returns the minimum likelihood that the best allele for a unit of evidence must have
     * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
     *
     * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
     */
    @VisibleForTesting
    public void filterPoorlyModeledEvidence(final ToDoubleFunction<EVIDENCE> log10MinTrueLikelihood) {
        Utils.validateArg(alleles.numberOfAlleles() > 0, "unsupported for read-likelihood collections with no alleles");

        final int numberOfSamples = samples.numberOfSamples();
        for (int s = 0; s < numberOfSamples; s++) {
            final int sampleIndex = s;
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(s);

            final int numberOfEvidence = sampleEvidence.size();
            final int[] indexesToRemove = IntStream.range(0, numberOfEvidence)
                    .filter(i -> maximumLikelihoodOverAllAlleles(sampleIndex, i) < log10MinTrueLikelihood.applyAsDouble(sampleEvidence.get(i)))
                    .toArray();

            // Retain the filtered evidence for later genotyping purposes
            final List<EVIDENCE> filtered = filteredEvidenceBySampleIndex.get(sampleIndex);
            Arrays.stream(indexesToRemove).forEach(idx -> {
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("disqualified read: " + idx + " "+((GATKRead)sampleEvidence.get(idx)).getName()+ " with max likelihood " +maximumLikelihoodOverAllAlleles(sampleIndex, idx) +" and threshold "+log10MinTrueLikelihood.applyAsDouble(sampleEvidence.get(idx)));
                }
                filtered.add(sampleEvidence.get(idx));
            });

            // Remove the evidence now
            removeEvidenceByIndex(s, indexesToRemove);
        }
    }



    private Object2IntMap<EVIDENCE> evidenceIndexBySampleIndex(final int sampleIndex) {
        final Object2IntMap<EVIDENCE> cached = evidenceIndexBySampleIndex.get(sampleIndex);
        return cached == null ? fillEvidenceToIndexCache(sampleIndex) : cached;
    }

    private Object2IntMap<EVIDENCE> fillEvidenceToIndexCache(int sampleIndex) {
        final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(sampleIndex);
        final int sampleEvidenceCount = sampleEvidence.size();
        final Object2IntMap<EVIDENCE> index = new Object2IntOpenHashMap<>(sampleEvidenceCount);
        index.defaultReturnValue(MISSING_INDEX);
        for (int r = 0; r < sampleEvidenceCount; r++) {
            index.put(sampleEvidence.get(r), r);
        }
        evidenceIndexBySampleIndex.set(sampleIndex, index);
        return index;
    }

    /**
     * Implements a likelihood matrix per sample given its index.
     */
    private final class SampleMatrix implements LikelihoodMatrix<EVIDENCE, A> {

        private final int sampleIndex;

        private SampleMatrix(final int sampleIndex) {
            this.sampleIndex = sampleIndex;
        }

        public List<EVIDENCE> evidence() {
            return sampleEvidence(sampleIndex);
        }

        @Override
        public List<A> alleles() {
            return AlleleLikelihoods.this.alleles();
        }

        @Override
        public void set(final int alleleIndex, final int evidenceIndex, final double value) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(evidenceIndex,  numberOfEvidences[sampleIndex]);
            valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex] = value;
        }

        @Override
        public double get(final int alleleIndex, final int evidenceIndex) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(evidenceIndex, numberOfEvidences[sampleIndex]);
            return valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex];
        }

        @Override
        public int indexOfAllele(final Allele allele) {
            Utils.nonNull(allele);
            return AlleleLikelihoods.this.indexOfAllele(allele);
        }

        @Override
        public int indexOfEvidence(final EVIDENCE evidence) {
            Utils.nonNull(evidence);
            return AlleleLikelihoods.this.evidenceIndex(sampleIndex, evidence);
        }

        @Override
        public int numberOfAlleles() {
            return alleles.numberOfAlleles();
        }

        @Override
        public int evidenceCount() {
            return numberOfEvidences[sampleIndex];
        }

        @Override
        public A getAllele(final int alleleIndex) {
            return AlleleLikelihoods.this.getAllele(alleleIndex);
        }

        @Override
        public EVIDENCE getEvidence(final int evidenceIndex) {
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(sampleIndex);
            Utils.validIndex(evidenceIndex, sampleEvidence.size());
            return sampleEvidence.get(evidenceIndex);
        }

        @Override
        public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
            Utils.nonNull(dest);
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            System.arraycopy(valuesBySampleIndex[sampleIndex][alleleIndex], 0, dest, offset, numberOfEvidences[sampleIndex]);
        }
    }
}
