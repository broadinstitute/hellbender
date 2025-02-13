package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Likelihood-based test for the consanguinuity among samples
 *
 * <p>This annotation estimates whether there is evidence of consanguinuity in a population. The higher the score, the
 * higher the chance that some samples are related. If samples are known to be related, a pedigree file can be provided so
 * that the calculation is only performed on founders and offspring are excluded.</p>
 *
 * <h3>Details</h3>
 * <p>The output is the inbreeding coefficient 'F' (fixation) statistic, which for large sample sizes converges to the probability
 * that an individual's two alleles are identical by descent, provided that consanguinity is the only source of deviation from Hardy-Weinberg equilibrium.
 * If this assumption is not true F may be negative and the excess heterozygosity often indicates an artifactual variant.
 * It is calculated as F = 1 - (# of het genotypes)/(# of het genotypes expected under Hardy-Weinberg equilibrium).  The number of het genotypes expected under Hardy-Weinberg equilibrium
 * is 2*(# of samples)*(ref allele frequency)*(alt allele frequency), where allele frequencies are calculated from the samples' genotypes.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The Inbreeding Coefficient annotation can only be calculated for cohorts containing at least 10 founder samples.</li>
 * <li>The Inbreeding Coefficient annotation can only be calculated for diploid samples.</li>
 * <li>The implementation here uses likelihoods rather than counts, so computed InbreedingCoeff values can be sensitive to small
 * changes in the low GQ range, as might happen after "reblocking" low coverage sites, as well as GQ0 versus no-call changes</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <p><b>ExcessHet</b> also describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Likelihood-based test for the consanguinity among samples (InbreedingCoeff)")
public final class InbreedingCoeff extends PedigreeAnnotation implements InfoFieldAnnotation, StandardAnnotation {

    private static final OneShotLogger oneShotLogger = new OneShotLogger(InbreedingCoeff.class);
    private static final int MIN_SAMPLES = 10;
    private static final boolean ROUND_GENOTYPE_COUNTS = false;  //if this changes update the caveats above

    public InbreedingCoeff(){
        super((Set<String>) null);
    }

    public InbreedingCoeff(final Set<String> founderIds){
        super(founderIds);
     }

    public InbreedingCoeff(final GATKPath pedigreeFile){
        super(pedigreeFile);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        final GenotypesContext genotypes = getFounderGenotypes(vc);
        if (genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant()) {
            oneShotLogger.warn("InbreedingCoeff will not be calculated at position " + vc.getContig() + ":" + vc.getStart() +
                    " and possibly subsequent; at least " + MIN_SAMPLES + " samples must have called genotypes");
            return Collections.emptyMap();
        }

        // See: https://github.com/broadinstitute/gatk/issues/7938
        if (genotypes.stream().anyMatch(g -> !g.hasLikelihoods())) {
            oneShotLogger.warn("InbreedingCoeff will not be calculated at position " + vc.getContig() + ":" + vc.getStart() +
                    " and possibly subsequent sites. This site lacks genotype likelihoods for at least one sample");

            return Collections.emptyMap();
        }

        final Pair<Integer, Double> sampleCountCoeff = calculateIC(vc, genotypes);
        final int sampleCount = sampleCountCoeff.getLeft();
        final double F = sampleCountCoeff.getRight();
        if (sampleCount < MIN_SAMPLES) {
            oneShotLogger.warn("InbreedingCoeff will not be calculated for at least one position; at least " + MIN_SAMPLES + " samples must have called genotypes");
            return Collections.emptyMap();
        }
        return Collections.singletonMap(getKeyNames().get(0), String.format("%.4f", F));
    }

    @VisibleForTesting
    static Pair<Integer, Double> calculateIC(final VariantContext vc, final GenotypesContext genotypes) {
        Utils.validate(!genotypes.isEmpty(), "Must provide genotypes to calculate InbreedingCoeff");
        final GenotypeCounts t = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes, ROUND_GENOTYPE_COUNTS);

        final double refCount = t.getRefs();
        final double hetCount = t.getHets();
        final double homCount = t.getHoms();
        // number of samples that have likelihoods
        final int sampleCount = (int) genotypes.stream().filter(g-> GenotypeUtils.isCalledAndDiploidWithLikelihoodsOrWithGQ(g) || GenotypeUtils.isDiploidWithLikelihoods(g)).count();

        final double p = ( 2.0 * refCount + hetCount ) / ( 2.0 * (refCount + hetCount + homCount) ); // expected reference allele frequency
        final double q = 1.0 - p; // expected alternative allele frequency
        final double expectedHets = 2.0 * p * q * sampleCount; //numbers of hets that would be expected based on the allele frequency (assuming Hardy Weinberg Equilibrium)
        final double F = 1.0 - ( hetCount / expectedHets ); // inbreeding coefficient

        return Pair.of(sampleCount, F);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY); }
}