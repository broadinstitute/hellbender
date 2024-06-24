package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections.ListUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypesCache;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.*;
import java.util.stream.Collectors;

public final class PosteriorProbabilitiesUtils {
    final private static int minSamplesToUseInputs = 10;

    private PosteriorProbabilitiesUtils(){}

    /**
     *  A class to wrangle all the various and sundry genotype posterior options,
     *  mostly from CalculateGenotypePosteriors
     */
    public final static class PosteriorProbabilitiesOptions {
        double snpPriorDirichlet;
        double indelPriorDirichlet;
        boolean useInputSamplesAlleleCounts;
        boolean useMLEAC;
        boolean ignoreInputSamplesForMissingResources;
        boolean useFlatPriorsForIndels;

        public PosteriorProbabilitiesOptions(final double snpPriorDirichlet,
                                             final double indelPriorDirichlet,
                                             final boolean useInputSamplesAlleleCounts,
                                             final boolean useMLEAC,
                                             final boolean ignoreInputSamplesForMissingResources,
                                             final boolean useFlatPriorsForIndels) {
            this.snpPriorDirichlet = snpPriorDirichlet;
            this.indelPriorDirichlet = indelPriorDirichlet;
            this.useInputSamplesAlleleCounts = useInputSamplesAlleleCounts;
            this.useMLEAC = useMLEAC;
            this.ignoreInputSamplesForMissingResources = ignoreInputSamplesForMissingResources;
            this.useFlatPriorsForIndels = useFlatPriorsForIndels;
        }
    }

         /**
         * Calculates phred-scaled posterior probabilities for genotypes given the data and allele frequency priors.
         *
         * @param vc1 input variant context
         * @param resources variants to use to calculate genotype priors (should already be checked for matching start position)
         * @param numRefSamplesFromMissingResources the number of reference samples to use for the prior (already accounted for resource missingness)
         * @return a VariantContext with PLs and PPs in the genotypes and a PG tag for the Phred-scaled prior in the INFO field if
         */
    public static VariantContext calculatePosteriorProbs(final VariantContext vc1,
                                                         final List<VariantContext> resources,
                                                         final int numRefSamplesFromMissingResources,
                                                         final PosteriorProbabilitiesOptions opts) {
        Utils.nonNull(vc1, "VariantContext vc1 is null");
        final Map<Allele,Integer> totalAlleleCounts = new HashMap<>();

        //only use discovered allele count for missing resources if there are at least 10 samples or if we have reference samples
        final boolean useDiscoveredACForMissing = !opts.ignoreInputSamplesForMissingResources && (vc1.getNSamples() >= minSamplesToUseInputs || numRefSamplesFromMissingResources != 0);

        //deletions introduce ref allele padding issues
        @SuppressWarnings("unchecked")
        List<VariantContext> allAlleles = ListUtils.union(resources, Arrays.asList(vc1));
        final Allele commonRef = GATKVariantContextUtils.determineReferenceAllele(allAlleles, new SimpleInterval(vc1));
        final List<Allele> origAllelesRemapped = ReferenceConfidenceVariantContextMerger.remapAlleles(vc1, commonRef);


        final int referenceAlleleCountForMissing = resources.isEmpty() ? HomoSapiensConstants.DEFAULT_PLOIDY*numRefSamplesFromMissingResources : 0;  //because we're dealing with diploids
        //store the allele counts for each allele in the variant priors
        for (final VariantContext r : resources) {
            if (r.getStart() == vc1.getStart()) {
                final List<Allele> remappedAlleles = ReferenceConfidenceVariantContextMerger.remapAlleles(r, commonRef);
                addAlleleCounts(totalAlleleCounts, r, remappedAlleles, !opts.useMLEAC);
            }
        }

        //add the allele counts from the input samples (if applicable)
        if ( (opts.useInputSamplesAlleleCounts && !resources.isEmpty()) || (resources.isEmpty() && useDiscoveredACForMissing)) {
            addAlleleCounts(totalAlleleCounts,vc1,origAllelesRemapped, !opts.useMLEAC);
        }

        //add zero allele counts for any reference alleles not seen in priors (if applicable)
        final int existingRefCounts = totalAlleleCounts.getOrDefault(commonRef, 0);
        totalAlleleCounts.put(commonRef, existingRefCounts + referenceAlleleCountForMissing);



        // now extract the counts of the alleles in vc1
        //if we allow deletion priors, then we need to harmonize the allele representations between vc1 and totalAlleleCounts
        final Set<Allele> allAllelesRemapped = totalAlleleCounts.keySet();
        final Set<Allele> resourceOnlyAlleles = new HashSet<>(allAllelesRemapped);
        resourceOnlyAlleles.removeAll(origAllelesRemapped);

        //this array will have the same allele ordering as the VC
        double[] alleleCounts = new double[origAllelesRemapped.size()];
        for (int i=0; i<origAllelesRemapped.size(); i++) {
            Allele a = origAllelesRemapped.get(i);
            if (a.length() == commonRef.length()) {  //use SNP prior for SNPs
                alleleCounts[i] = opts.snpPriorDirichlet + totalAlleleCounts.getOrDefault(a, 0);
            }
            else if (a.isSymbolic()) {   //use the greater of the priors for non-ref
                alleleCounts[i] = Math.max(opts.snpPriorDirichlet, opts.indelPriorDirichlet) + totalAlleleCounts.getOrDefault(a, 0);
            }
            else {   //use indel prior for non-SNPs
                alleleCounts[i] = opts.indelPriorDirichlet + totalAlleleCounts.getOrDefault(a, 0);
            }
        }

        //put resource alleles not in input VC into non-ref, if applicable
        int nonRefInd = vc1.getAlleleIndex(Allele.NON_REF_ALLELE);
        if (nonRefInd != -1) {
            //the non-ref "allele" gets a single pseudo count since in some ways it's one allele -- if we treated it as all alleles not present it would get infinite pseudocounts
            alleleCounts[nonRefInd] = Math.max(opts.snpPriorDirichlet, opts.indelPriorDirichlet) + resourceOnlyAlleles.stream().mapToDouble(a -> totalAlleleCounts.get(a)).sum();
        }

        final List<double[]> likelihoods = vc1.getGenotypes().stream().map(g -> parsePosteriorsIntoProbSpace(g)).collect(Collectors.toList());

        final boolean useFlatPriors = (!vc1.isSNP() && opts.useFlatPriorsForIndels) || (resources.isEmpty() && !useDiscoveredACForMissing && numRefSamplesFromMissingResources == 0);

        final List<double[]> posteriors = calculatePosteriorProbs(likelihoods,alleleCounts,vc1.getMaxPloidy(HomoSapiensConstants.DEFAULT_PLOIDY), useFlatPriors);

        final GenotypesContext newContext = GenotypesContext.create();
        for ( int genoIdx = 0; genoIdx < vc1.getNSamples(); genoIdx ++ ) {
            final GenotypeBuilder builder = new GenotypeBuilder(vc1.getGenotype(genoIdx));
            builder.phased(vc1.getGenotype(genoIdx).isPhased());
            if ( posteriors.get(genoIdx) != null ) {
                GATKVariantContextUtils.makeGenotypeCall(vc1.getMaxPloidy(HomoSapiensConstants.DEFAULT_PLOIDY), builder,
                        GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, posteriors.get(genoIdx), vc1.getAlleles(), vc1.getGenotype(genoIdx), null);

                builder.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY,
                            Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(posteriors.get(genoIdx)).getAsPLs()));
                //TODO: refactor me -- this is to make up for the GGVCFs changes fixing GQ0 hom-refs
                if (!builder.makeWithShallowCopy().hasPL()) {
                    final int maxPosteriorIndex = MathUtils.maxElementIndex(posteriors.get(genoIdx));
                    builder.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(maxPosteriorIndex, posteriors.get(genoIdx)));
                    builder.alleles(GenotypesCache.get(vc1.getGenotype(genoIdx).getPloidy(), maxPosteriorIndex).asAlleleList(vc1.getAlleles()));
                }
            }
            newContext.add(builder.make());
        }

        final List<Integer> priors = Utils.listFromPrimitives(
                GenotypeLikelihoods.fromLog10Likelihoods(getDirichletPrior(alleleCounts, vc1.getMaxPloidy(HomoSapiensConstants.DEFAULT_PLOIDY),useFlatPriors)).getAsPLs());

        final VariantContextBuilder builder = new VariantContextBuilder(vc1).genotypes(newContext);
        final boolean isHomRefBlock = vc1.getAlternateAlleles().size() == 1 && vc1.getAlleles().contains(Allele.NON_REF_ALLELE);
        if (!isHomRefBlock) {
            // update/add in the AC, AF, and AN attributes and genotype prior -- this is kind of a cheat because we're doing it
            // outside the annotation engine, which has the header line responsibility for these fields.  For CGP the
            // headers should already be there, but for HaplotypeCaller GVCF mode we'll add the headers in HaplotypeCallerEngine::makeVCFHeader.
            VariantContextUtils.calculateChromosomeCounts(builder.attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors), true);
        }
        return builder.make();
    }

    public static int[] parsePosteriorsIntoPhredSpace(Genotype genotype) {
        final Object PPfromVCF = genotype.getExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);

        if (PPfromVCF == null){
            return genotype.getPL();
        } else if (PPfromVCF instanceof String) {
            final String PPstring = (String) PPfromVCF;
            //samples not in trios will have PP tag like ".,.,." if family priors are applied
            return PPstring.charAt(0)=='.' ? genotype.getPL() :
                    Arrays.stream(PPstring.split(",")).mapToInt(i -> Integer.parseInt(i)).toArray();
        } else {
            return Arrays.stream(extractInts(PPfromVCF)).toArray();
        }
    }

    public static double[] parsePosteriorsIntoProbSpace(Genotype genotype) {
        final Object PPfromVCF = genotype.getExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);

        if (PPfromVCF == null){
            return getLikelihoodsVector(genotype);
        } else if (PPfromVCF instanceof String) {
            final String PPstring = (String) PPfromVCF;
            //samples not in trios will have PP tag like ".,.,." if family priors are applied
            return PPstring.charAt(0)=='.' ? getLikelihoodsVector(genotype) :
                    Arrays.stream(PPstring.split(",")).mapToDouble(s -> Double.parseDouble(s)/-10.0).toArray();
        } else {
            return Arrays.stream(extractInts(PPfromVCF)).mapToDouble(i -> i/-10.0).toArray();
        }
    }

    // return the double[] of likelihoods if available, otherwise null
    private static double[] getLikelihoodsVector(Genotype genotype) {
        return hasRealLikelihoods(genotype) ? genotype.getLikelihoods().getAsVector() : null;
    }

    /**
     * return true if the genotype contains real data for likelihoods, not just zero placeholder
     */
    private static boolean hasRealLikelihoods(Genotype genotype) {
        if (!genotype.hasLikelihoods()) {
            return false;
        }
        if (genotype.hasDP() && genotype.getDP() == 0) {
            //if there's no depth and the PLs are all zero, then there's really no data here; don't rely on no-call
            return MathUtils.arrayMax(genotype.getPL()) > 0;
        }
        return true;
    }

    /**
     * Given genotype likelihoods and known allele counts, calculate the posterior probabilities
     * over the genotype states
     * @param genotypeLikelihoods - the genotype likelihoods for the individual
     * @param knownAlleleCountsByAllele - the known allele counts in the population. For AC=2 AN=12 site, this is {10,2}
     * @param ploidy - the ploidy to assume
     * @param useFlatPriors - if true, apply flat priors to likelihoods in order to calculate posterior probabilities
     * @return - the posterior genotype likelihoods
     */
    protected static List<double[]> calculatePosteriorProbs(final List<double[]> genotypeLikelihoods,
                                                            final double[] knownAlleleCountsByAllele,
                                                            final int ploidy,
                                                            final boolean useFlatPriors) {
        if ( ploidy != 2 ) {
            throw new IllegalStateException("Genotype posteriors not yet implemented for ploidy != 2");
        }

        final double[] genotypePriorByAllele = getDirichletPrior(knownAlleleCountsByAllele,ploidy, useFlatPriors);
        final List<double[]> posteriors = new ArrayList<>(genotypeLikelihoods.size());
        for ( final double[] likelihoods : genotypeLikelihoods ) {
            double[] posteriorProbabilities = null;

            if ( likelihoods != null ) {
                if ( likelihoods.length != genotypePriorByAllele.length ) {
                    throw new IllegalStateException(String.format("Likelihoods not of correct size: expected %d, observed %d",
                            knownAlleleCountsByAllele.length*(knownAlleleCountsByAllele.length+1)/2,likelihoods.length));
                }

                posteriorProbabilities = new double[genotypePriorByAllele.length];
                for ( int genoIdx = 0; genoIdx < likelihoods.length; genoIdx ++ ) {
                    posteriorProbabilities[genoIdx] = likelihoods[genoIdx] + genotypePriorByAllele[genoIdx];
                }

                posteriorProbabilities = MathUtils.normalizeLog10(posteriorProbabilities);

            }

            posteriors.add(posteriorProbabilities);
        }

        return posteriors;
    }

    // convenience function for a single genotypelikelihoods array. Just wraps.
    @VisibleForTesting
    static double[] calculatePosteriorProbs(final double[] genotypeLikelihoods,
                                            final double[] knownAlleleCountsByAllele,
                                            final int ploidy,
                                            final boolean useFlatPriors) {
        return calculatePosteriorProbs(Arrays.asList(genotypeLikelihoods),knownAlleleCountsByAllele,ploidy, useFlatPriors).get(0);
    }


    /**
     * Given known allele counts (whether external, from the sample, or both), calculate the prior distribution
     * over genotype states. This assumes
     *   1) Random sampling of alleles (known counts are unbiased, and frequency estimate is Dirichlet)
     *   2) Genotype states are independent (Hardy-Weinberg)
     * These assumptions give rise to a Dirichlet-Multinomial distribution of genotype states as a prior
     * (the "number of trials" for the multinomial is simply the ploidy)
     * @param knownCountsByAllele - the known counts per allele. For an AC=2, AN=12 site this is {10,2}
     * @param ploidy - the number of chromosomes in the sample. For now restricted to 2.
     * @return - the Dirichlet-Multinomial distribution over genotype states
     */
    @VisibleForTesting
    static double[] getDirichletPrior(final double[] knownCountsByAllele, final int ploidy, final boolean useFlatPrior) {
        if ( ploidy != 2 ) {
            throw new IllegalStateException("Genotype priors not yet implemented for ploidy != 2");
        }

        // multi-allelic format is
        // AA AB BB AC BC CC AD BD CD DD ...
        final double sumOfKnownCounts = MathUtils.sum(knownCountsByAllele);
        final double[] priors = new double[knownCountsByAllele.length*(knownCountsByAllele.length+1)/2];
        int priorIndex = 0;
        for ( int allele2 = 0; allele2 < knownCountsByAllele.length; allele2++ ) {
            for ( int allele1 = 0; allele1 <= allele2; allele1++) {
                if (useFlatPrior) {
                    priors[priorIndex++] = 1.0;
                } else {
                    final int[] counts = new int[knownCountsByAllele.length];
                    counts[allele1] += 1;
                    counts[allele2] += 1;
                    priors[priorIndex++] = MathUtils.dirichletMultinomial(knownCountsByAllele,counts);
                }
            }
        }

        return priors;
    }

    /**
     * Parse counts for each allele
     * @param counts - Map to store and return data
     * @param context - line to be parsed from the input VCF file
     * @param useAC - use allele count annotation value from VariantContext (vs. MLEAC)
     */
    private static void addAlleleCounts(final Map<Allele,Integer> counts, final VariantContext context, final List<Allele> remappedAlleles, final boolean useAC) {
        final int[] ac;
        //use MLEAC value...
        if ( context.hasAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) && ! useAC ) {
            ac = getAlleleCounts(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, context);
        }
        //...unless specified by the user in useAC or unless MLEAC is absent
        else if ( context.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            ac = getAlleleCounts(VCFConstants.ALLELE_COUNT_KEY, context);
        }
        //if VariantContext annotation doesn't contain AC or MLEAC then get the data from direct evaluation
        else {
            ac = new int[context.getAlternateAlleles().size()];
            int idx = 0;
            for ( final Allele allele : context.getAlternateAlleles() ) {
                ac[idx++] = context.getCalledChrCount(allele);
            }
        }

        for ( int i = 0; i < context.getAlleles().size(); i++ ) {
            final Allele allele = remappedAlleles.get(i);
            final Allele origAllele = context.getAlleles().get(i);
            final int count;
            if ( allele.isReference() ) {
                //since the allele count for the reference allele is not given in the VCF format,
                //calculate it from the allele number minus the total counts for alternate alleles
                if ( context.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
                    count = Math.max(context.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY,-1) - (int) MathUtils.sum(ac),0); //occasionally an MLEAC value will sneak in that's greater than the AN
                } else {
                    count = Math.max(context.getCalledChrCount() - (int) MathUtils.sum(ac),0);
                }
            } else {
                count = ac[context.getAlternateAlleles().indexOf(origAllele)];
            }
            //if this allele isn't in the map yet, add it
            if ( ! counts.containsKey(allele) ) {
                counts.put(allele,0);
            }
            //add the count for the current allele to the existing value in the map
            counts.put(allele,count + counts.get(allele));
        }
    }

    /**
     * Retrieve allele count data from VariantContext using VCFkey, checks for correct number of values in VCF
     * @param VCFkey VariantContext annotation tag of interest (should be AC or MLEAC)
     * @param context VariantContext from which to extract the data
     * @return int[] with allele count data
     */
    private static int[] getAlleleCounts(final String VCFkey, final VariantContext context) {
        final Object alleleCountsFromVCF = context.getAttribute(VCFkey);
        if ( alleleCountsFromVCF instanceof List) {
            if ( ((List) alleleCountsFromVCF).size() != context.getAlternateAlleles().size() ) {
                throw new UserException(String.format("Variant does not contain the same number of MLE allele counts as alternate alleles for record at %s:%d", context.getContig(), context.getStart()));
            }
        }
        else if ( alleleCountsFromVCF instanceof String || alleleCountsFromVCF instanceof Integer) {//here length is 1
            if (context.getAlternateAlleles().size() != 1) {
                throw new UserException(String.format("Variant does not contain the same number of MLE allele counts as alternate alleles for record at %s:%d", context.getContig(), context.getStart()));
            }
        }
        return extractInts(alleleCountsFromVCF);
    }

    /**
     * Check the formatting on the Object returned by a call to VariantContext::getAttribute() and parse appropriately
     * @param integerListContainingVCField - Object returned by a call to VariantContext::getAttribute()
     * @return - array of ints
     *
     * //Note: if we're working with a integerListContainingVCField that's read directly out of the file it will be a String but
     * if it gets pulled from a VariantContext object built elsewhere in the code it will be an Integer or a List,
     */
    @SuppressWarnings("unchecked")
    public static int[] extractInts(final Object integerListContainingVCField) {
        List<Integer> mleList = null;
        if ( integerListContainingVCField instanceof List) {
            if ( ((List) integerListContainingVCField).get(0) instanceof String) {
                mleList = new ArrayList<>(((List) integerListContainingVCField).size());
                for ( final Object s : ((List)integerListContainingVCField)) {
                    mleList.add(Integer.parseInt((String) s));
                }
            } else {
                mleList = (List<Integer>) integerListContainingVCField;
            }
        } else if ( integerListContainingVCField instanceof Integer) {
            mleList = Arrays.asList((Integer) integerListContainingVCField);
        } else if ( integerListContainingVCField instanceof String) {
            mleList = Arrays.asList(Integer.parseInt((String)integerListContainingVCField));
        }
        Utils.nonNull( mleList, () -> String.format("VCF does not have properly formatted %s or %s.",
                    GATKVCFConstants.MLE_ALLELE_COUNT_KEY, VCFConstants.ALLELE_COUNT_KEY));

        final int[] mle = new int[mleList.size()];

        if ( ! ( mleList.get(0) instanceof Integer) ) {
            throw new IllegalStateException("BUG: The AC values should be an Integer, but was " + mleList.get(0).getClass().getCanonicalName());
        }

        for ( int idx = 0; idx < mle.length; idx++) {
            mle[idx] = mleList.get(idx);
        }

        return mle;
    }
}
