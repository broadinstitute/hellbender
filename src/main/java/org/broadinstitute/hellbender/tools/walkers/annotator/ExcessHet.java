
package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Phred-scaled p-value for exact test of excess heterozygosity.
 *
 * <p>This annotation estimates the probability of the called samples exhibiting excess heterozygosity with respect to the null hypothesis that the samples are unrelated. The higher the score, the
 * higher the chance that the variant is a technical artifact or that there is consanguinuity among the samples. In
 * contrast to Inbreeding Coefficient, there is no minimal number of samples for this annotation. If samples are known to be related, a pedigree file can be provided so
 * that the calculation is only performed on founders and offspring are excluded.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>This annotation uses the implementation from
 * <a href='https://doi.org/10.1086/429864'>Wigginton, Cutler, and Abecasis 2005.</a>
 * to calculate a right-sided p-value (the probability of getting the observed number of heterozygotes or higher).
 * <p>However, note that in releases prior to and including GATK 4.2.3.0, an additional mid-p correction
 * not discussed in that reference was applied. That is, if n_{AB} is the observed number of heterozygous genotypes,
 * N is the total number of samples, and n_A is the number of minor alleles, then the p-value was previously reported as:
 * \sum_{N_{AB} = n_{AB}}^{N} P(N_{AB} | N, n_A) - P(n_{AB} | N, n_A) / 2.
 * Similar mid-p corrections are discussed in e.g. <a href='https://doi.org/10.1515/sagmb-2012-0039'>Graffelman and Moreno 2013</a>.</p>
 * <p>This mid-p correction was removed in subsequent releases, so that the p-value is now reported as:
 * \sum_{N_{AB} = n_{AB}}^{N} P(N_{AB} | N, n_A). This p-value is thus concordant with the P_{high} p-value
 * discussed by Wigginton et al., as well as the ExcHet annotation reported by
 * <a href='https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html'>bcftools</a>.</p>
 * <p>See discussion at the corresponding GitHub
 * <a href='https://github.com/broadinstitute/gatk/issues/7392'>issue</a> and
 * <a href='https://github.com/broadinstitute/gatk/pull/7394'>pull request</a>.</p>
 *
 * <h3>Caveat</h3>
 * <p>The Excess Heterozygosity annotation can only be calculated for diploid samples.</p>
 *
 * <h3>Related annotations</h3>
 * <p><b>InbreedingCoeff</b> also describes the heterozygosity of the called samples, though without explicitly taking into account the number of samples</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Phred-scaled p-value for exact test of excess heterozygosity (ExcessHet)")
public final class ExcessHet extends PedigreeAnnotation implements InfoFieldAnnotation, StandardAnnotation {
    private static final OneShotLogger oneShotLogger = new OneShotLogger(ExcessHet.class);

    private static final double MIN_NEEDED_VALUE = 1.0E-16;
    private static final boolean ROUND_GENOTYPE_COUNTS = true;
    
    public static final double PHRED_SCALED_MIN_P_VALUE = -10.0 * Math.log10(MIN_NEEDED_VALUE);
    public static final int NUMBER_OF_GENOTYPE_COUNTS = 3;

    public ExcessHet(final Set<String> founderIds){
        super(founderIds);
    }

    public ExcessHet(final GATKPath pedigreeFile){
        super(pedigreeFile);
    }

    public ExcessHet() {
        this((Set<String>) null);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        GenotypesContext genotypes = getFounderGenotypes(vc);
        if (genotypes == null || genotypes.isEmpty() || !vc.isVariant()) {
            return Collections.emptyMap();
        }

        // See: https://github.com/broadinstitute/gatk/issues/7938
        if (genotypes.stream().anyMatch(g -> !g.hasLikelihoods())) {
            oneShotLogger.warn("ExcessHet will not be calculated at position " + vc.getContig() + ":" + vc.getStart() +
                    " and possibly subsequent sites. This site lacks genotype likelihoods for at least one sample");
            return Collections.emptyMap();
        }

        final Pair<Integer, Double> sampleCountEH = calculateEH(vc, genotypes);
        final int sampleCount = sampleCountEH.getLeft();
        final double eh =  sampleCountEH.getRight();

        if (sampleCount < 1) {
            return Collections.emptyMap();
        }
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", eh));
    }

    @VisibleForTesting
    static Pair<Integer, Double> calculateEH(final VariantContext vc, final GenotypesContext genotypes) {
        final GenotypeCounts t = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes, ROUND_GENOTYPE_COUNTS);
        // number of samples that have likelihoods
        final int sampleCount = (int) genotypes.stream().filter(GenotypeUtils::isCalledAndDiploidWithLikelihoodsOrWithGQ).count();

        return calculateEH(t, sampleCount);
    }

    @VisibleForTesting
    public static Pair<Integer, Double> calculateEH(final GenotypeCounts t, final int sampleCount) {
        final int refCount = (int) t.getRefs();
        final int hetCount = (int) t.getHets();
        final int homCount = (int) t.getHoms();

        final double pval = exactTest(hetCount, refCount, homCount);

        // If the actual phredPval would be infinity we will probably still filter out just a very large number
        // Since the method does not guarantee precision for any p-value smaller than 1e-16, we can return the phred scaled version
        if (pval < 10e-60) {
            return Pair.of(sampleCount, PHRED_SCALED_MIN_P_VALUE);
        }
        final double phredPval = -10.0 * Math.log10(pval) + 0.; // We add 0. to prevent -0.0000 from being output

        return Pair.of(sampleCount, phredPval);
    }

    /**
     * Note that this method is not accurate for very small p-values. Beyond 1.0E-16 there is no guarantee that the
     * p-value is accurate, just that it is in fact smaller than 1.0E-16 (and therefore we should filter it). It would
     * be more computationally expensive to calculate accuracy beyond a given threshold. Here we have enough accuracy
     * to filter anything below a p-value of 10E-6.
     *
     * Discussion at https://github.com/broadinstitute/gatk/issues/7392 and
     * https://github.com/broadinstitute/gatk/pull/7394 may be useful for developers.
     *
     * @param hetCount Number of observed hets (n_ab)
     * @param refCount Number of observed homRefs (n_aa)
     * @param homCount Number of observed homVars (n_bb)
     * @return Right sided p-value or the probability of getting the observed or higher number of hets given the sample
     * size (N) and the observed number of allele a (rareCopies)
     */
    @VisibleForTesting
    static double exactTest(final int hetCount, final int refCount, final int homCount) {
        Utils.validateArg(hetCount >= 0, "Het count cannot be less than 0");
        Utils.validateArg(refCount >= 0, "Ref count cannot be less than 0");
        Utils.validateArg(homCount >= 0, "Hom count cannot be less than 0");

        // Split into observed common allele and rare allele
        final int obsHomR;
        final int obsHomC;
        if (refCount < homCount) {
            obsHomR = refCount;
            obsHomC = homCount;
        } else {
            obsHomR = homCount;
            obsHomC = refCount;
        }

        final int rareCopies = 2 * obsHomR + hetCount;
        final int N = hetCount + obsHomC + obsHomR;

        final double[] probs = new double[rareCopies + 1];

        // Find (something close to the) mode for the midpoint
        int mid = (int) Math.floor(rareCopies * (2.0 * N - rareCopies) / (2.0 * N - 1.0));
        if ((mid % 2) != (rareCopies % 2)) {
            mid++;
        }

        probs[mid] = 1.0;
        double mysum = 1.0;

        // Calculate probabilities from midpoint down
        int currHets = mid;
        int currHomR = (rareCopies - mid) / 2;
        int currHomC = N - currHets - currHomR;

        while (currHets >= 2) {
            final double potentialProb = probs[currHets] * currHets * (currHets - 1.0) / (4.0 * (currHomR + 1.0) * (currHomC + 1.0));
            if (potentialProb < MIN_NEEDED_VALUE) {
                break;
            }

            probs[currHets - 2] = potentialProb;
            mysum = mysum + probs[currHets - 2];

            // 2 fewer hets means one additional homR and homC each
            currHets = currHets - 2;
            currHomR = currHomR + 1;
            currHomC = currHomC + 1;
        }

        // Calculate probabilities from midpoint up
        currHets = mid;
        currHomR = (rareCopies - mid) / 2;
        currHomC = N - currHets - currHomR;

        while (currHets <= rareCopies - 2) {
            final double potentialProb = probs[currHets] * 4.0 * currHomR * currHomC / ((currHets + 2.0) * (currHets + 1.0));
            if (potentialProb < MIN_NEEDED_VALUE) {
                break;
            }

            probs[currHets + 2] = potentialProb;
            mysum = mysum + probs[currHets + 2];

            // 2 more hets means 1 fewer homR and homC each
            currHets = currHets + 2;
            currHomR = currHomR - 1;
            currHomC = currHomC - 1;
        }

        final double rightPval = probs[hetCount] / mysum;
        // Check if we observed the highest possible number of hets
        if (hetCount == rareCopies) {
            return Math.max(0., Math.min(1., rightPval));
        }
        return Math.max(0., Math.min(1., rightPval + StatUtils.sum(Arrays.copyOfRange(probs, hetCount + 1, probs.length)) / mysum));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.EXCESS_HET_KEY);
    }

    //@Override
    public String getRawKeyName() {
        return GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY;
    }

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     *
     * @param ref         the reference context for this annotation
     * @param vc          the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    //@Override
    public Map<String, Object> annotateRawData(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        return null;
    }

    /**
     * Combine raw data, typically during the merging of raw data contained in multiple gVCFs as in CombineGVCFs and the
     * preliminary merge for GenotypeGVCFs
     *
     * @param allelesList   The merged allele list across all variants being combined/merged
     * @param listOfRawData The raw data for all the variants being combined/merged
     * @return A key, Object map containing the annotations generated by this combine operation
     */
    //@Override
    public Map<String, Object> combineRawData(List<Allele> allelesList, List<ReducibleAnnotationData<?>> listOfRawData) {
        return null;
    }

    /**
     * Calculate the final annotation value from the raw data which was generated by either annotateRawData or calculateRawData
     *
     * @param vc         -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return A key, Object map containing the finalized annotations generated by this operation to be added to the code
     */
    //@Override
    public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
        if (vc.hasAttribute(getRawKeyName())) {
            List<Integer> counts = vc.getAttributeAsIntList(getRawKeyName(), 0);
            if (counts.size() != NUMBER_OF_GENOTYPE_COUNTS) {
                throw new IllegalStateException("Genotype counts for ExcessHet (" + getRawKeyName() + ") should have three values: homozygous reference, heterozygous with one ref allele, and homozygous variant/heterozygous non-reference");
            }
            final GenotypeCounts t = new GenotypeCounts(counts.get(0), counts.get(1), counts.get(2));
            final Pair<Integer, Double> sampleCountEH = calculateEH(t, counts.get(0)+counts.get(1)+counts.get(2));
            final int sampleCount = sampleCountEH.getLeft();
            final double eh =  sampleCountEH.getRight();

            if (sampleCount < 1) {
                return Collections.emptyMap();
            }
            return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", eh));
        }
        return null;
    }
}
