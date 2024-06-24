package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.readorientation.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.util.*;
import java.util.stream.IntStream;

public class ReadOrientationFilter extends Mutect2VariantFilter {
    private Map<String, ArtifactPriorCollection> artifactPriorCollections = new HashMap<>();

    public ReadOrientationFilter(final List<File> readOrientationPriorTables) {
        readOrientationPriorTables.stream()
                .forEach(file -> {
                    final ArtifactPriorCollection artifactPriorCollection = ArtifactPriorCollection.readArtifactPriors(file);
                    artifactPriorCollections.put(artifactPriorCollection.getSample(), artifactPriorCollection);
                });
    }

    public static int[] getF1R2(final Genotype g) {
        return VariantContextGetters.getAttributeAsIntArray(g, GATKVCFConstants.F1R2_KEY, () -> null, 0);
    }

    public static int[] getF2R1(final Genotype g) {
        return VariantContextGetters.getAttributeAsIntArray(g, GATKVCFConstants.F2R1_KEY, () -> null, 0);
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        if (!vc.isSNP() && !vc.isMNP()){
            return 0;
        }

        final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

        vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .forEach(g -> {
                    final double artifactPosterior = artifactProbability(referenceContext, vc, g);
                    final int[] ADs = g.getAD();
                    final int altCount = (int) MathUtils.sum(ADs) - ADs[0];

                    depthsAndPosteriors.add(ImmutablePair.of(altCount, artifactPosterior));
                });

        final double artifactPosterior = weightedMedianPosteriorProbability(depthsAndPosteriors);
        return artifactPosterior;
    }

    @Override
    public String filterName() { return GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME; }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.READ_ORIENTATION_QUAL_KEY);
    }

    @Override
    protected List<String> requiredInfoAnnotations() { return Collections.emptyList(); }


    @VisibleForTesting
    double artifactProbability(final ReferenceContext referenceContext, final VariantContext vc, final Genotype g) {
        // As of June 2018, genotype is hom ref iff we have the normal sample, but this may change in the future
        // TODO: handle MNVs
        if (g.isHomRef() || (!vc.isSNP() && !vc.isMNP()) ){
            return 0;
        } else if (!artifactPriorCollections.containsKey(g.getSampleName())) {
            return 0;
        }

        final double[] tumorLods = VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, -1);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAllele = vc.getAlternateAllele(indexOfMaxTumorLod);
        final byte[] altBases = altAllele.getBases();

        // for MNVs, treat each base as an independent substitution and take the maximum of all error probabilities
        return IntStream.range(0, altBases.length).mapToDouble(n -> {
            final Nucleotide altBase = Nucleotide.valueOf(new String(new byte[] {altBases[n]}));

            return artifactProbability(referenceContext, vc.getStart() + n, g, indexOfMaxTumorLod, altBase);
        }).max().orElse(0.0);

    }

    private double artifactProbability(final ReferenceContext referenceContext, final int refPosition, final Genotype g, final int indexOfMaxTumorLod, final Nucleotide altBase) {
        final String refContext = referenceContext.getKmerAround(refPosition, F1R2FilterConstants.REF_CONTEXT_PADDING);
        if (refContext ==  null || refContext.contains("N")){
            return 0;
        }

        Utils.validate(refContext.length() == 2 * F1R2FilterConstants.REF_CONTEXT_PADDING + 1,
                String.format("kmer must have length %d but got %d", 2 * F1R2FilterConstants.REF_CONTEXT_PADDING + 1, refContext.length()));

        final Nucleotide refAllele = F1R2FilterUtils.getMiddleBase(refContext);

        if (!(g.hasExtendedAttribute(GATKVCFConstants.F1R2_KEY) && g.hasExtendedAttribute(GATKVCFConstants.F2R1_KEY))) {
            return 0;
        }

        final int[] f1r2 = getF1R2(g);
        final int[] f2r1 = getF2R1(g);
        final int refCount =  f1r2[0] + f2r1[0];
        final int altF1R2 = f1r2[indexOfMaxTumorLod + 1];
        final int altF2R1 = f2r1[indexOfMaxTumorLod + 1];
        final int altCount = altF1R2 + altF2R1;
        final Optional<ArtifactPrior> artifactPrior = artifactPriorCollections.get(g.getSampleName()).get(refContext);

        if (! artifactPrior.isPresent()){
            return 0;
        }

        final int depth = refCount + altCount;

        final double[] posterior = LearnReadOrientationModelEngine.computeResponsibilities(refAllele, altBase, altCount, altF1R2, depth, artifactPrior.get().getPi(), true);

        final double posteriorOfF1R2 = posterior[ArtifactState.getF1R2StateForAlt(altBase).ordinal()];
        final double posteriorOfF2R1 = posterior[ArtifactState.getF2R1StateForAlt(altBase).ordinal()];

        return Math.max(posteriorOfF1R2, posteriorOfF2R1);
    }
}
