package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class NormalArtifactFilter extends Mutect2VariantFilter {
    private static final double MIN_NORMAL_ARTIFACT_RATIO = 0.1;    // don't call normal artifact if allele fraction in normal is much smaller than allele fraction in tumor
    private static final int IMPUTED_NORMAL_BASE_QUALITY = 30;  // only used if normal base quality annotation fails somehow

    private final double normalPileupPValueThreshold;
    public NormalArtifactFilter(final double normalPileupPValueThreshold) {
        this.normalPileupPValueThreshold = normalPileupPValueThreshold;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        final int[] tumorAlleleDepths = filteringEngine.sumADsOverSamples(vc, true, false);
        final int tumorDepth = (int) MathUtils.sum(tumorAlleleDepths);
        final int tumorAltDepth = tumorAlleleDepths[indexOfMaxTumorLod + 1];

        final int[] normalAlleleDepths = filteringEngine.sumADsOverSamples(vc, false, true);
        final int normalDepth = (int) MathUtils.sum(normalAlleleDepths);
        final int normalAltDepth = normalAlleleDepths[indexOfMaxTumorLod + 1];

        // if normal AF << tumor AF, don't filter regardless of LOD
        final double tumorAlleleFraction = (double) tumorAltDepth / tumorDepth;
        final double normalAlleleFraction = normalDepth == 0 ? 0 : (double) normalAltDepth / normalDepth;

        if (normalAlleleFraction < MIN_NORMAL_ARTIFACT_RATIO * tumorAlleleFraction)  {
            return 0.0;
        }

        final double[] normalArtifactNegativeLogOdds = MathUtils.applyToArrayInPlace(VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_ARTIFACT_LOG_10_ODDS_KEY), x -> -MathUtils.log10ToLog(x));
        final double normalArtifactProbability = filteringEngine.posteriorProbabilityOfNormalArtifact(normalArtifactNegativeLogOdds[indexOfMaxTumorLod]);

        // the normal artifact log odds misses artifacts whose support in the normal consists entirely of low base quality reads
        // Since a lot of low-BQ reads is itself evidence of an artifact, we filter these by hand via an estimated LOD
        // that uses the average base quality of *ref* reads in the normal
        final int medianRefBaseQuality = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY, IMPUTED_NORMAL_BASE_QUALITY).get(0);
        final double normalPValue = 1 - new BinomialDistribution(null, normalDepth, QualityUtils.qualToErrorProb(medianRefBaseQuality))
                .cumulativeProbability(normalAltDepth - 1);

        return normalPValue < normalPileupPValueThreshold ? 1.0 : normalArtifactProbability;
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.empty();
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME;
    }

    @Override
    protected List<String> requiredInfoAnnotations() {
        return Arrays.asList(GATKVCFConstants.NORMAL_ARTIFACT_LOG_10_ODDS_KEY, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY);
    }
}
