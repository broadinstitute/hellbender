package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;

public class FilteredHaplotypeFilter extends Mutect2VariantFilter {
    private static final double GERMLINE_PROBABILITY_TO_IGNORE_NORMAL_ARTIFACT = 0.25;
    private final double maxIntraHaplotypeDistance;

    // for each pgt + pid phasing string, a list of loci-error probability pairs
    private Map<String, List<Pair<Integer, Double>>> accumulatingPhasedProbabilities = new HashMap<>();

    private Map<String, List<Pair<Integer, Double>>> phasedProbabilities = new HashMap<>();

    public FilteredHaplotypeFilter(final double maxIntraHaplotypeDistance) {
        this.maxIntraHaplotypeDistance = maxIntraHaplotypeDistance;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        // use phasing of tumor genotype with greatest allele fraction
        final Genotype tumorGenotype = vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .max(Comparator.comparingDouble(g -> MathUtils.arrayMax(VariantContextGetters.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                        () -> new double[] {0.0}, 0.0)))).get();

        final Optional<String> phasingString = makePhasingString(tumorGenotype);
        if (!phasingString.isPresent()) {
            return 0.0;
        }

        // note that we use the learned probabilities from the previous pass
        final List<Pair<Integer, Double>> phasedProbs = phasedProbabilities.get(phasingString.get());

        if (phasedProbs == null) {
            return 0.0;
        }

        return phasedProbs.stream()
                .filter(pair -> Math.abs(pair.getLeft() - vc.getStart()) <= maxIntraHaplotypeDistance)
                .mapToDouble(Pair::getRight)
                .max().orElse(0.0);
    }

    @Override
    protected void accumulateDataForLearning(final VariantContext vc, final ErrorProbabilities errorProbabilities, final Mutect2FilteringEngine filteringEngine) {
        // we record the maximum non-sequencing, non-germline, artifact probability that is not from this filter itself
        final Map<Mutect2Filter, List<Double>> probabilitiesByFilter = errorProbabilities.getProbabilitiesByFilter();

        final double germlineProbability = probabilitiesByFilter.entrySet().stream()
                .filter(e -> e.getKey().filterName() == GATKVCFConstants.GERMLINE_RISK_FILTER_NAME)
                .flatMap(e -> e.getValue().stream())  // the value is a list of double, we need the max of all the lists
                .max(Double::compareTo).orElse(0.0);

        // the normal artifact filter often lights up when there's a non-artifactual germline event, which we don't want here
        final boolean ignoreNormalArtifact = germlineProbability > GERMLINE_PROBABILITY_TO_IGNORE_NORMAL_ARTIFACT;

        final double artifactProbability = probabilitiesByFilter.entrySet().stream()
                .filter(e -> e.getKey().errorType() != ErrorType.NON_SOMATIC)
                .filter(e -> !(ignoreNormalArtifact && e.getKey().filterName() == GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME))
                .filter(e -> !e.getKey().filterName().equals(filterName())) // exclude the haplotype filter itself, which would be circular
                .flatMap(e -> e.getValue().stream())  // the value is a list of double, we need the max of all the lists
                .max(Double::compareTo).orElse(0.0);

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (!filteringEngine.isTumor(tumorGenotype)) {
                continue;
            }

            final Optional<String> phasingString = makePhasingString(tumorGenotype);

            if (!phasingString.isPresent()) {
                continue;
            }

            if (!accumulatingPhasedProbabilities.containsKey(phasingString.get())) {
                accumulatingPhasedProbabilities.put(phasingString.get(), new ArrayList<>());
            }

            accumulatingPhasedProbabilities.get(phasingString.get()).add(ImmutablePair.of(vc.getStart(), artifactProbability));
        }
    }


    @Override
    protected void clearAccumulatedData() {
        accumulatingPhasedProbabilities = new HashMap<>();
    }

    @Override
    protected void learnParameters() {
        // move the accumulating probabilities to the learned probabiities
        phasedProbabilities = accumulatingPhasedProbabilities;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredInfoAnnotations() { return Collections.emptyList(); }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() { return Optional.empty(); }

    // concatenate the PGT and PID strings, if present
    private static Optional<String> makePhasingString(final Genotype genotype) {
        final String pgt = (String) genotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, null);
        final String pid = (String) genotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, null);
        return (pgt == null || pid == null) ? Optional.empty() : Optional.of(pgt + pid);
    }
}
