package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.genotyper.GenotypePriorCalculator;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Base class for genotyper engines.
 */
public abstract class GenotypingEngine<Config extends StandardCallerArgumentCollection> {

    private final static int TOO_LONG_PL = 100000;

    protected final AlleleFrequencyCalculator alleleFrequencyCalculator;

    protected final Config configuration;

    protected VariantAnnotatorEngine annotationEngine;

    protected Logger logger;

    protected OneShotLogger oneShotLogger;

    protected final int numberOfGenomes;

    protected final SampleList samples;

    // the top of the queue is the upstream deletion that ends first
    // note that we can get away with ordering just by the end and not the contig as long as we preserve the invariant
    // that everything in this queue belongs to the same contig
    private final PriorityQueue<Locatable> upstreamDeletionsLoc = new PriorityQueue<>(Comparator.comparingInt(Locatable::getEnd));

    private final boolean doAlleleSpecificCalcs;

    /**
     * Construct a new genotyper engine, on a specific subset of samples.
     *
     * @param configuration engine configuration object.
     * @param samples subset of sample to work on identified by their names. If {@code null}, the full toolkit
     *                    sample set will be used instead.
     * @param doAlleleSpecificCalcs Whether the AS_QUAL key should be calculated and added to newly genotyped variants.
     *
     * @throws IllegalArgumentException if any of {@code samples}, {@code configuration} is {@code null} or if {@code samples} is empty.
     */
    protected GenotypingEngine(final Config configuration,
                               final SampleList samples,
                               final boolean doAlleleSpecificCalcs) {
        this.configuration = Utils.nonNull(configuration, "the configuration cannot be null");
        Utils.validate(!samples.asListOfSamples().isEmpty(), "the sample list cannot be null or empty");
        this.samples = samples;
        this.doAlleleSpecificCalcs = doAlleleSpecificCalcs;
        logger = LogManager.getLogger(getClass());
        oneShotLogger = new OneShotLogger(logger);
        numberOfGenomes = this.samples.numberOfSamples() * configuration.genotypeArgs.samplePloidy;
        alleleFrequencyCalculator = AlleleFrequencyCalculator.makeCalculator(configuration.genotypeArgs);
    }

    /**
     * Changes the logger for this genotyper engine.
     *
     * @param logger new logger.
     *
     * @throws IllegalArgumentException if {@code logger} is {@code null}.
     */
    public void setLogger(final Logger logger) {
        this.logger = Utils.nonNull(logger, "the logger cannot be null");
    }

    public Set<VCFInfoHeaderLine> getAppropriateVCFInfoHeaders() {
        final Set<VCFInfoHeaderLine> headerInfo = new LinkedHashSet<>();
        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
        }
        return headerInfo;
    }

    /**
     * Changes the annotation engine for this genotyping-engine.
     *
     * @param annotationEngine the new annotation engine (can be {@code null}).
     */
    public void setAnnotationEngine(final VariantAnnotatorEngine annotationEngine) {
        this.annotationEngine = annotationEngine;
    }

    /**
     * Returns a reference to the engine configuration
     *
     * @return never {@code null}.
     */
    public Config getConfiguration() {
        return configuration;
    }

    /**
     * Main entry function to calculate genotypes of a given VC with corresponding GLs that is shared across
     * genotypers (namely GGVCFs and HC).
     *
     * Completes a variant context with genotype calls and associated annotations given the genotype likelihoods and
     * the model that need to be applied.  Hom-ref likelihoods can be approximated from GQs, but if no genotype has
     * likelihoods then that variant is either all-ref or contains variants with no likelihoods, and in both cases
     * we want to exit.
     *
     * @param vc                                 Input variant context to complete.
     * @return                                   VC with assigned genotypes
     */
    public VariantContext calculateGenotypes(final VariantContext vc, final GenotypePriorCalculator gpc, final List<Event> givenAlleles) {
        // if input VC can't be genotyped, exit with either null VCC or, in case where we need to emit all sites, an empty call
        if (cannotBeGenotyped(vc) || vc.getNSamples() == 0) {
            return null;
        }

        final int defaultPloidy = configuration.genotypeArgs.samplePloidy;
        final int maxAltAlleles = configuration.genotypeArgs.maxAlternateAlleles;

        VariantContext reducedVC = vc;
        if (maxAltAlleles < vc.getAlternateAlleles().size()) {
            final List<Allele> allelesToKeep = AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, defaultPloidy, maxAltAlleles, false);
            final GenotypesContext reducedGenotypes = allelesToKeep.size() == 1 ? GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy) :
                    AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), allelesToKeep, gpc,
                            GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);  //with no PLs in some reblocked GVCFs, no-calls are just going to cause problems, so keep 0/0 genotypes as such without trying to recall
            reducedVC = new VariantContextBuilder(vc).alleles(allelesToKeep).genotypes(reducedGenotypes).make();
        }

        //Calculate the expected total length of the PL arrays for this VC to warn the user in the case that they will be exceptionally large
        final long maxPLLength = GenotypeLikelihoods.numLikelihoods(reducedVC.getNAlleles(), reducedVC.getMaxPloidy(defaultPloidy));
        if(maxPLLength >= TOO_LONG_PL) {
            oneShotLogger.warn("Length of PL arrays for this VC(position:" + reducedVC.getStart() + ", alleles:" + reducedVC.getNAlleles()
                    + ", ploidy:" + reducedVC.getMaxPloidy(defaultPloidy) + ") is likely to reach " + maxPLLength + ", so processing may take a long time.");
        }

        final AFCalculationResult AFresult = alleleFrequencyCalculator.calculate(reducedVC, defaultPloidy);
        final Set<Allele> forcedAlleles = AssemblyBasedCallerUtils.allelesConsistentWithGivenAlleles(givenAlleles, vc);
        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult, vc, forcedAlleles);

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        final double log10Confidence =
                    !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                            ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0;

        // Add 0.0 removes -0.0 occurrences.
        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        // skip this if we are already looking at a vc with NON_REF as the first alt allele i.e. if we are in GenotypeGVCFs
        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic) && !emitAllActiveSites()
                && noAllelesOrFirstAlleleIsNotNonRef(outputAlternativeAlleles.alleles) && forcedAlleles.isEmpty()) {
            return null;
        }

        // return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
        if (! emitAllActiveSites() && outputAlternativeAlleles.alleles.size() == 1 && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0))) {
            return null;
        }

        // start constructing the resulting VC
        final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
        recordDeletions(vc, outputAlleles);

        final VariantContextBuilder builder = new VariantContextBuilder(callSourceString(), vc.getContig(), vc.getStart(), vc.getEnd(), outputAlleles);

        builder.log10PError(log10Confidence);
        if ( ! passesCallThreshold(phredScaledConfidence) ) {
            builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
        }

        // create the genotypes
        //TODO: omit subsetting if output alleles is not a proper subset of vc.getAlleles
        final GenotypesContext genotypes = outputAlleles.size() == 1 ? GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy) :
                AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), outputAlleles, gpc, configuration.genotypeArgs.genotypeAssignmentMethod);

        if (configuration.genotypeArgs.usePosteriorProbabilitiesToCalculateQual && hasPosteriors(genotypes)) {
            final double log10NoVariantPosterior = phredNoVariantPosteriorProbability(outputAlleles, genotypes) * -.1;
            final double qualUpdate = !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                    ? log10NoVariantPosterior + 0.0 : MathUtils.log10OneMinusPow10(log10NoVariantPosterior) + 0.0;
            if (!Double.isNaN(qualUpdate)) {
                builder.log10PError(qualUpdate);
            }
        }

        // calculating strand bias involves overwriting data structures, so we do it last
        final Map<String, Object> attributes = composeCallAttributes(vc, outputAlternativeAlleles.alternativeAlleleMLECounts(),
                AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()),genotypes);

        return builder.genotypes(genotypes).attributes(attributes).make();
    }

    protected double phredNoVariantPosteriorProbability(final List<Allele> alleles, final GenotypesContext gc) {
        return gc.stream()
                .mapToDouble(gt -> extractPNoAlt(alleles, gt))
                .filter(d -> !Double.isNaN(d))
                .reduce(Double.NaN, (a, b) -> Double.isNaN(a) ? b : (Double.isNaN(b) ? a : a + b) );
    }

    private double extractPNoAlt(final List<Allele> alleles, final Genotype gt) {
        final double[] gpArray = VariantContextGetters.getAttributeAsDoubleArray(gt, VCFConstants.GENOTYPE_POSTERIORS_KEY, () -> new double[]{Double.NaN}, Double.NaN);
        return extractPNoAlt(alleles, gt, gpArray);
    }

    private double extractPNoAlt(final List<Allele> alleles, final Genotype gt, final double[] posteriors) {
        if (!alleles.contains(Allele.SPAN_DEL)) {
            return posteriors[0] - Math.max(0, QualityUtils.phredSum(posteriors));
        } else {
            // here we need to get indices of genotypes composed of REF and * alleles
            final int ploidy = gt.getPloidy();
            final int spanDelIndex = alleles.indexOf(Allele.SPAN_DEL);
            // allele counts are in the GenotypeLikelihoodCalculator format of {ref index, ref count, span del index, span del count}
            final double[] nonVariantLog10Posteriors = IntStream.rangeClosed(0, ploidy)
                    .map(n -> GenotypeIndexCalculator.alleleCountsToIndex(0, ploidy - n, spanDelIndex, n))
                    .mapToDouble(n -> posteriors[n])
                    .toArray();

            // when the only alt allele is the spanning deletion the probability that the site is non-variant
            // may be so close to 1 that finite precision error in log10SumLog10 (called by phredSum) yields a positive value,
            // which is bogus.  Thus we cap it at 0. See AlleleFrequencyCalculator.
            return Math.max(0, QualityUtils.phredSum(nonVariantLog10Posteriors)) - Math.max(0, QualityUtils.phredSum(posteriors));
        }
    }

    //note that posteriors could be Phred-scaled or not depending on the VCF version (!) -- see also issue https://github.com/broadinstitute/gatk/issues/7390
    private boolean hasPosteriors(final GenotypesContext gc) {
        return gc.stream().anyMatch(g -> g.hasExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY));
    }

    public VariantContext calculateGenotypes(final VariantContext vc) {
        return calculateGenotypes(vc, null, Collections.emptyList());
    }

    @VisibleForTesting
    static boolean noAllelesOrFirstAlleleIsNotNonRef(List<Allele> altAlleles) {
        Utils.nonNull(altAlleles);
        return altAlleles.isEmpty() ||  altAlleles.get(0) != (Allele.NON_REF_ALLELE);
    }

    /**
     * What string to use as source of variant-context generated by this genotyper-engine.
     * @return never {@code null} nor empty.
     */
    protected abstract String callSourceString();

    /**
     * Holds information about the alternative allele subsetting based on supporting evidence, genotyping and
     * output modes.
     */
    private static class OutputAlleleSubset {
        private  final List<Allele> alleles;
        private  final boolean siteIsMonomorphic;
        private  final List<Integer> mleCounts;

        private OutputAlleleSubset(final List<Allele> alleles, final List<Integer> mleCounts, final boolean siteIsMonomorphic) {
            Utils.nonNull(alleles, "alleles");
            Utils.nonNull(mleCounts, "mleCounts");
            this.siteIsMonomorphic = siteIsMonomorphic;
            this.alleles = alleles;
            this.mleCounts = mleCounts;
        }

        private List<Allele> outputAlleles(final Allele referenceAllele) {
            return Stream.concat(Stream.of(referenceAllele), alleles.stream()).collect(Collectors.toList());
        }

        public List<Integer> alternativeAlleleMLECounts() {
            return mleCounts;
        }
    }


    /**
     * Provided the exact mode computations it returns the appropriate subset of alleles that progress to genotyping.
     * @param afCalculationResult the allele fraction calculation result.
     * @param vc the variant context
     * @param forcedAlleles alleles from the vc input that are consistent with forced alleles in the assembly region {@link AssemblyBasedCallerUtils#allelesConsistentWithGivenAlleles}
     * @return information about the alternative allele subsetting {@code null}.
     */
    private OutputAlleleSubset calculateOutputAlleleSubset(final AFCalculationResult afCalculationResult, final VariantContext vc, final Set<Allele> forcedAlleles) {
        final List<Allele> outputAlleles = new ArrayList<>();
        final List<Integer> mleCounts = new ArrayList<>();
        boolean siteIsMonomorphic = true;
        final List<Allele> alleles = afCalculationResult.getAllelesUsedInGenotyping();
        final int alternativeAlleleCount = alleles.size() - 1;
        int referenceAlleleSize = 0;

        for (final Allele allele : alleles) {
            if (allele.isReference() ) {
                referenceAlleleSize = allele.length();
            } else {
                // we want to keep the NON_REF symbolic allele but only in the absence of a non-symbolic allele, e.g.
                // if we combined a ref / NON_REF gVCF with a ref / alt gVCF
                final boolean isNonRefWhichIsLoneAltAllele = alternativeAlleleCount == 1 && allele.equals(Allele.NON_REF_ALLELE);
                final boolean isPlausible = afCalculationResult.passesThreshold(allele, configuration.genotypeArgs.standardConfidenceForCalling);

                //it's possible that the upstream deletion that spanned this site was not emitted, mooting the symbolic spanning deletion allele
                final boolean isSpuriousSpanningDeletion = GATKVCFConstants.isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc);

                final boolean toOutput = (isPlausible || forceKeepAllele(allele) || isNonRefWhichIsLoneAltAllele || forcedAlleles.contains(allele) ) && !isSpuriousSpanningDeletion;

                siteIsMonomorphic &= !(isPlausible && !isSpuriousSpanningDeletion);

                if (toOutput) {
                    outputAlleles.add(allele);
                    mleCounts.add(afCalculationResult.getAlleleCountAtMLE(allele));
                }
            }
        }

        return new OutputAlleleSubset(outputAlleles,mleCounts,siteIsMonomorphic);
    }

    void clearUpstreamDeletionsLoc() {
        upstreamDeletionsLoc.clear();
    }

    /**
     *  Record emitted deletions in order to remove downstream spanning deletion alleles that are not covered by any emitted deletion.
     *  In addition to recording new deletions, this method culls previously-recorded deletions that end before the current variant
     *  context.  This assumes that variants are traversed in order.
     *
     * @param vc                VariantContext, potentially multiallelic and potentially containing one or more deletion alleles
     * @param emittedAlleles    The subset of the variant's alt alleles that are actually emitted
     */
    @VisibleForTesting
    void recordDeletions(final VariantContext vc, final Collection<Allele> emittedAlleles) {
        while (!upstreamDeletionsLoc.isEmpty() && (!upstreamDeletionsLoc.peek().contigsMatch(vc) || upstreamDeletionsLoc.peek().getEnd() < vc.getStart())) {
            upstreamDeletionsLoc.poll();
        }

        for (final Allele allele : emittedAlleles) {
            if (allele.isNonRefAllele()) {
                continue;
            }
            final int deletionSize = vc.getReference().length() - allele.length();

            // In a deletion
            if (deletionSize > 0) {
                final SimpleInterval genomeLoc = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart() + deletionSize);
                upstreamDeletionsLoc.add(genomeLoc);
            }
        }
    }

    /**
     * Is the variant context covered by an upstream deletion?
     *
     * @param vc    variant context
     * @return  true if the location is covered by an upstream deletion, false otherwise
     */
    boolean isVcCoveredByDeletion(final VariantContext vc) {
        // note: the code below seems like it's duplicating Locatable.overlaps, but here if the upstream deletion
        // has the same start as the vc we don't want to count it
        return !upstreamDeletionsLoc.isEmpty() && upstreamDeletionsLoc.stream()
                .anyMatch(loc -> loc.getContig().equals(vc.getContig()) && loc.getStart() < vc.getStart() && vc.getStart() <= loc.getEnd());

    }

    /**
     * Checks whether even if the allele is not well supported by the data, we should keep it for genotyping.
     *
     * @param allele target allele.
     *
     * @return {@code true} iff we need to keep this alleles even if does not seem plausible.
     */
    protected abstract boolean forceKeepAllele(final Allele allele);


    /**
     * Checks whether the variant context has too many alternative alleles for progress to genotyping the site.
     * <p>
     *     AF calculation may get into trouble with too many alternative alleles.
     * </p>
     *
     * @param vc the variant context to evaluate.
     *
     * @throws NullPointerException if {@code vc} is {@code null}.
     *
     * @return {@code true} iff there is too many alternative alleles based on
     * {@link GenotypeLikelihoods#MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED}.
     */
    protected final boolean cannotBeGenotyped(final VariantContext vc) {
        if (vc.getNAlleles() <= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED
                // likelihoods may be missing when reading from GenomicsDB if there are more alts that GDB args allow
                // so ensure all genotypes (outside of 0/0 and ./.) have likelihoods
            && vc.getGenotypes().stream().filter( g -> !(g.isNoCall() || g.isHomRef()) ).allMatch(Genotype::hasLikelihoods)
                // if all sites are no calls or hom refs then keep a site where any samples have likelihoods
            && vc.getGenotypes().stream().anyMatch(Genotype::hasLikelihoods)) {
            return false;
        }
        // protect against too many alternate alleles that we can't even run AF on:
        if (vc.getNAlleles() > GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED) {
            logger.warn("Attempting to genotype more than " + GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED +
                    " alleles. Site will be skipped at location " + vc.getContig() + ":" + vc.getStart());
            return true;
        } else {
            logger.warn("Some genotypes contained insufficient data to recalculate site and allele qualities. Site will be skipped at location "
                    + vc.getContig() + ":" + vc.getStart());
            return true;
        }
    }

    /**
     * Whether or not all calls at any region over the activity threshold should be emitted regardless of confidence.
     * This does not necessarily emit calls for all sites in a region (see {@link OutputMode}).
     */
    protected boolean emitAllActiveSites() {
        return configuration.outputMode == OutputMode.EMIT_ALL_ACTIVE_SITES;
    }

    protected final boolean passesEmitThreshold(final double conf, final boolean bestGuessIsRef) {
        return (configuration.outputMode == OutputMode.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) &&
                passesCallThreshold(conf);
    }

    protected final boolean passesCallThreshold(final double conf) {
        return conf >= configuration.genotypeArgs.standardConfidenceForCalling;
    }

    protected Map<String,Object> composeCallAttributes(final VariantContext vc, final List<Integer> alleleCountsofMLE,
                                                       final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes) {
        final Map<String, Object> attributes = new LinkedHashMap<>();

        // add the MLE AC and AF annotations
        if (!alleleCountsofMLE.isEmpty()) {
            attributes.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
            final List<Double> MLEfrequencies = calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes);
            attributes.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
        }

        if (doAlleleSpecificCalcs){
            List<Integer> perAlleleQuals = new ArrayList<>();
            //Per-allele quals are not calculated for biallelic sites
            if (AFresult.getAllelesUsedInGenotyping().size() > 2) {
                for (final Allele a : allAllelesToUse) {
                    if (a.isNonReference()) {
                        //*-10 to convert from log10-scale to Phred-scale, as QUALs are typically represented
                        perAlleleQuals.add((int)Math.round(AFresult.getLog10PosteriorOfAlleleAbsent(a)*-10));
                    }
                }
            }
            else {
                //*-10 to convert from log10-scale to Phred-scale, as QUALs are typically represented
                perAlleleQuals.add((int)Math.round(AFresult.log10ProbOnlyRefAlleleExists()*-10));
            }

            attributes.put(GATKVCFConstants.AS_QUAL_KEY, perAlleleQuals.stream().mapToInt(q -> Math.round(q)).boxed().collect(Collectors.toList()));
        }

        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
            attributes.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, vc.getAlternateAlleles().size());
        }

        return attributes;
    }

    private List<Double> calculateMLEAlleleFrequencies(final List<Integer> alleleCountsofMLE, final GenotypesContext genotypes) {
        final long AN = genotypes.stream().flatMap(g -> g.getAlleles().stream()).filter(Allele::isCalled).count();
        return alleleCountsofMLE.stream().map(AC -> Math.min(1.0, (double) AC / AN)).collect(Collectors.toList());
    }

    /**
     * Calculates the active state profile value for a single sample.
     *
     * @param log10GenotypeLikelihoods the single sample genotype likelihoods.
     * @return log10 probability from 0 to -Infinity.
     */
    public double calculateSingleSampleRefVsAnyActiveStateProfileValue(final double[] log10GenotypeLikelihoods) {
        Utils.nonNull(log10GenotypeLikelihoods, "the input likelihoods cannot be null");
        return alleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior(log10GenotypeLikelihoods, true);
    }
}
