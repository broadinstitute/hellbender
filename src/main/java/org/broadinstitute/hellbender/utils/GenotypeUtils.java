package org.broadinstitute.hellbender.utils;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeIndexCalculator;
import picard.util.MathUtil;

public final class GenotypeUtils {
    final static int TYPICAL_BASE_QUALITY = 30;
    //from the genotype likelihoods equations assuming the SNP ref conf model with no mismatches
    //PL[2] = GQ; scaleFactor = PL[3]/GQ ~ -10 * DP * log10(P_error) / (-10 * DP * log10(1/ploidy)) where BASE_QUALITY = -10 * log10(P_error)
    final static int PLOIDY_2_HOM_VAR_SCALE_FACTOR = (int)Math.round(TYPICAL_BASE_QUALITY /-10.0/Math.log10(.5));

    private GenotypeUtils(){}

    /**
     * Returns true if the genotype is a diploid genotype with likelihoods.
     */
    public static boolean isDiploidWithLikelihoods(final Genotype g) {
        return Utils.nonNull(g).hasLikelihoods() && g.getPloidy() == 2;
    }

    /**
     * Returns true if the genotype is a called diploid genotype with likelihoods.
     */
    public static boolean isCalledAndDiploidWithLikelihoodsOrWithGQ(final Genotype g) {
        return Utils.nonNull(g).isCalled() && g.getPloidy() == 2 && (Utils.nonNull(g).hasLikelihoods() || g.hasGQ()) ;
    }

    /**
     * Returns a triple of ref/het/hom genotype "counts".
     *
     * The exact meaning of the count is dependent on the rounding behavior.
     * if {@code roundContributionFromEachGenotype}: the counts are discrete integer counts of the most probable genotype for each {@link Genotype}
     * else: they are the sum of the normalized likelihoods of each genotype and will not be integers
     *
     * Skips non-diploid genotypes.
     *
     *
     * @param vc the VariantContext that the {@link Genotype}s originated from, non-null
     * @param genotypes a GenotypesContext containing genotypes to count, these must be a subset of {@code vc.getGenotypes()}, non-null
     * @param roundContributionFromEachGenotype if this is true, the normalized likelihood from each genotype will be rounded before
     *                                          adding to the total count
     */
    public static GenotypeCounts computeDiploidGenotypeCounts(final VariantContext vc, final GenotypesContext genotypes,
                                                              final boolean roundContributionFromEachGenotype){
        Utils.nonNull(vc, "vc");
        Utils.nonNull(genotypes, "genotypes");

        final int idxAA = 0;
        final int idxAB = 1;
        final int idxBB = 2;

        double genotypeWithTwoRefsCount = 0;  //i.e. 0/0
        double genotypesWithOneRefCount = 0;  //e.g. 0/1, 0/2, etc.
        double genotypesWithNoRefsCount = 0;  //e.g. 1/1, 1/2, 2/2, etc.

        for (final Genotype g : genotypes) {
            //if we don't have the data we need then skip this genotype (equivalent to no-call)
            if (!isDiploidWithLikelihoods(g)
                    && !isCalledAndDiploidWithLikelihoodsOrWithGQ(g)) {
                continue;
            }

            if (!g.hasLikelihoods() && g.isHomRef()) {
                if (roundContributionFromEachGenotype) {
                    genotypeWithTwoRefsCount += 1;
                    continue;
                } else if (g.getGQ() == 0) {
                    genotypeWithTwoRefsCount += 1.0/3;
                    genotypesWithOneRefCount += 1.0/3;
                    genotypesWithNoRefsCount += 1.0/3;
                    continue;
                } else {
                    genotypeWithTwoRefsCount += QualityUtils.qualToProb(g.getGQ());
                    genotypesWithOneRefCount += 1 - QualityUtils.qualToProb(g.getGQ());
                    //assume last likelihood is negligible
                    continue;
                }
            }

            // Genotype::getLikelihoods returns a new array, so modification in-place is safe
            final double[] normalizedLikelihoods;
            if (g.hasLikelihoods()) {
                normalizedLikelihoods = MathUtils.normalizeFromLog10ToLinearSpace(g.getLikelihoods().getAsVector());
            } else {
                throw new IllegalStateException("Genotype has no likelihoods: " + g.toString());
            }
            final double[] biallelicLikelihoods;

            //if there are multiple alts, use the biallelic PLs for the best alt
            if (vc.getAlternateAlleles().size() > 1 ) {
                //check for
                int maxInd = MathUtil.indexOfMax(normalizedLikelihoods);
                GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(maxInd);
                if (alleles.alleleIndex1 != 0 && alleles.alleleIndex2 != 0) {
                    //all likelihoods go to genotypesWithNoRefsCount because no ref allele is called
                    genotypesWithNoRefsCount++;
                    continue;
                }

                double maxLikelihood = normalizedLikelihoods[idxAB];
                int hetIndex = idxAB;
                int varIndex = idxBB;
                for (final Allele currAlt : vc.getAlternateAlleles()) {
                    final int[] idxVector = vc.getGLIndicesOfAlternateAllele(currAlt);
                    int tempIndex = idxVector[1];
                    if (normalizedLikelihoods[tempIndex] > maxLikelihood) {
                        maxLikelihood = normalizedLikelihoods[tempIndex];
                        hetIndex = tempIndex;
                        varIndex = idxVector[2];
                    }
                }
                biallelicLikelihoods = MathUtils.normalizeSumToOne(new double[] {normalizedLikelihoods[idxAA], normalizedLikelihoods[hetIndex], normalizedLikelihoods[varIndex]});
            }
            else {
                biallelicLikelihoods = normalizedLikelihoods;
            }

            double refLikelihood = biallelicLikelihoods[idxAA];
            double hetLikelihood = biallelicLikelihoods[idxAB];
            double varLikelihood = biallelicLikelihoods[idxBB];

            //NOTE: rounding is special cased for [0,0,X] and [X,0,0] PLs because likelihoods can come out as [0.5, 0.5, 0] and both counts round up
            if( roundContributionFromEachGenotype ) {
                genotypeWithTwoRefsCount += MathUtils.fastRound(refLikelihood);
                if (refLikelihood != hetLikelihood) {   //if GQ = 0 (specifically [0,0,X] PLs) count as homRef and don't add to the other counts
                    genotypesWithOneRefCount += MathUtils.fastRound(hetLikelihood);
                }
                if (varLikelihood != hetLikelihood) {  //if GQ = 0 (specifically [X,0,0] PLs) count as het and don't count as variant
                    genotypesWithNoRefsCount += MathUtils.fastRound(varLikelihood); //need specific genotypesWithNoRefsCount (rather than complement of the others) for PL[0,0,0] case
                }
            } else {
                genotypeWithTwoRefsCount += refLikelihood;
                genotypesWithOneRefCount += hetLikelihood;
                genotypesWithNoRefsCount += varLikelihood;
            }
        }
        return new GenotypeCounts(genotypeWithTwoRefsCount, genotypesWithOneRefCount, genotypesWithNoRefsCount);
    }

    /**
     * Do we have (or can we infer) likelihoods necessary for allele frequency calculation?
     * Some reblocked and/or DRAGEN GVCFs omit likelihoods for ref blocks, but we can estimate them
     * Ref blocks without likelihoods need to be retrieved from GenomicDB or CombineGVCFs with the --call-genotypes arg
     * If GenomicsDB max alt threshold is too low, non-reference genotypes may also be missing PLs -- we can't estimate, so reject those
     * @param g a genotype of unknown call and ploidy
     * @return  true if we have enough info for AF calculation
     */
    public static boolean genotypeIsUsableForAFCalculation(Genotype g) {
        return g.hasLikelihoods() || (g.isHomRef() && g.hasGQ() && 2 == g.getPloidy());
    }

    /**
     * Make approximate likelihoods for a diploid genotype (with arbitrary allele count) without PLs given genotype quality GQ.
     *
     * The method is as follows:
     * 1) For the biallelic diploid case with alleles A,B, the genotype likelihoods would be
     * AA: 0, AB: GQ, BB: PLOIDY_2_HOM_VAR_SCALE_FACTOR * GQ
     *
     * 2) For arbitrary allele count, set the genotype likelihoods as
     * AA: same as AA in the biallelic case
     * AB, AC, AD etc:  same as AB in the biallelic case
     * BB, BC, CC, BD etc: same as BB in the biallelic case
     *
     * WARNING: this calculation is completely bogus!  Legacy javadoc said: "For a hom-ref, as long as we have GQ we can
     * make a very accurate QUAL calculation since the hom-var likelihood should make a minuscule contribution."  Basically,
     * the bogusness of this method doesn't matter because the voodoo only involves the very small hom-var contribution.
     * That's true, but for multiallelics it incorrectly assigns the same GQ to every het genotype, essentially deflating
     * the qulity of the hom-ref call.
     *
     * In summary, the effect of this calculation is to depress the QUAL of multiallelic hom refs by a small amount for
     * no reason whatsoever.
     * @param g a diploid genotype with GQ
     * @param nAlleles number of alleles (including reference)
     * @return log10 likelihoods
     */
    public static double[] makeApproximateDiploidLog10LikelihoodsFromGQ(Genotype g, int nAlleles) {
        Utils.validate(g.getPloidy() == 2, "This method can only be used to approximate likelihoods for diploid genotypes");
        Utils.validate(g.hasGQ(), "Genotype must have GQ in order to approximate PLs");

        final int homRefLikelihood = 0;
        final int hetLikelihood = g.getGQ();
        final int homVarLikelihood = PLOIDY_2_HOM_VAR_SCALE_FACTOR * g.getGQ();

        final int[] PLs = new int[GenotypeIndexCalculator.genotypeCount(2, nAlleles)];
        //TODO: replace with GenotypesCache::iterator
        for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(2, nAlleles)) {
                PLs[gac.index()] = gac.index() == 0 ? homRefLikelihood : (gac.containsAllele(0) ? hetLikelihood : homVarLikelihood);
        }

        return GenotypeLikelihoods.fromPLs(PLs).getAsVector();  //fromPLs converts from Phred-space back to log10-space
    }

    public static boolean shouldBeCalled(final Genotype g) {
        return !g.isNonInformative() || (g.hasGQ() && g.getGQ() > 0);
    }
}
