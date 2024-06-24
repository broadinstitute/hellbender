package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticGenotypingEngine implements AutoCloseable {
    protected static final Logger logger = LogManager.getLogger(SomaticGenotypingEngine.class);

    private final M2ArgumentCollection MTAC;
    private final Set<String> normalSamples;
    final boolean hasNormal;
    protected VariantAnnotatorEngine annotationEngine;
    private final Optional<Mutect3DatasetEngine> mutect3DatasetEngine;

    // If MTAC.minAF is non-zero we softly cut off allele fractions below minAF with a Beta prior of the form Beta(1+epsilon, 1); that is
    // the prior on allele fraction f is proportional to f^epsilon.  If epsilon is small this prior vanishes as f -> 0
    // and very rapidly becomes flat.  We choose epsilon such that minAF^epsilon = 0.5.
    private final double refPseudocount = 1;
    private final double altPseudocount;

    public SomaticGenotypingEngine(final M2ArgumentCollection MTAC, final Set<String> normalSamples,
                                   final VariantAnnotatorEngine annotationEngine,
                                   final SAMFileHeader header, final SAMSequenceDictionary sequenceDictionary) {
        this.MTAC = MTAC;
        altPseudocount = MTAC.minAF == 0.0 ? 1 : 1 - Math.log(2)/Math.log(MTAC.minAF);

        this.normalSamples = normalSamples;
        hasNormal = !normalSamples.isEmpty();
        this.annotationEngine = annotationEngine;

        mutect3DatasetEngine = MTAC.mutect3Dataset == null ? Optional.empty() :
                Optional.of(new Mutect3DatasetEngine(MTAC.mutect3Dataset, MTAC.mutect3TrainingDataMode, MTAC.maxRefCountForMutect3,
                        MTAC.maxAltCountForMutect3, MTAC.mutect3NonArtifactRatio, normalSamples, header, sequenceDictionary));
        Utils.validateArg(!(MTAC.mutect3Dataset == null && MTAC.mutect3TrainingDataMode), "No dataset file specified for Mutect3 training data mode.");
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param logReadLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param activeRegionWindow                     Active window
     * @param withBamOut                            whether to annotate reads in readLikelihoods for future writing to bamout
     * @param emitRefConf                           generate reference confidence (GVCF) data?
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     */
    public CalledHaplotypes callMutations(
            final AlleleLikelihoods<GATKRead, Haplotype> logReadLikelihoods,
            final AssemblyResultSet assemblyResultSet,
            final ReferenceContext referenceContext,
            final SimpleInterval activeRegionWindow,
            final FeatureContext featureContext,
            final List<Event> givenAlleles,
            final SAMFileHeader header,
            final boolean withBamOut,
            final boolean emitRefConf,
            final Set<Integer> suspiciousLocations) {
        Utils.nonNull(logReadLikelihoods);
        Utils.validateArg(logReadLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(activeRegionWindow);

        final List<Haplotype> haplotypes = logReadLikelihoods.alleles();

        EventMap.buildEventMapsForHaplotypes(haplotypes, assemblyResultSet.getFullReferenceWithPadding(),
                assemblyResultSet.getPaddedReferenceLoc(), MTAC.assemblerArgs.debugAssembly, MTAC.maxMnpDistance);

        final List<Integer> eventStarts = EventMap.getEventStartPositions(haplotypes).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        if(withBamOut){
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(logReadLikelihoods, activeRegionWindow);
        }

        if (MTAC.likelihoodArgs.phredScaledGlobalReadMismappingRate > 0) {
            logReadLikelihoods.normalizeLikelihoods(NaturalLogUtils.qualToLogErrorProb(MTAC.likelihoodArgs.phredScaledGlobalReadMismappingRate), true);
        }
        final AlleleLikelihoods<Fragment, Haplotype> logFragmentLikelihoods = logReadLikelihoods.groupEvidence(MTAC.independentMates ? read -> read : GATKRead::getName, Fragment::createAndAvoidFailure);

        final Set<Event> potentialSomaticEventsInRegion = new HashSet<>();
        for( final int loc : eventStarts ) {
            final List<VariantContext> eventsAtThisLoc = AssemblyBasedCallerUtils.getVariantsFromActiveHaplotypes(loc, haplotypes, false);
            VariantContext merged = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( merged == null ) {
                continue;
            }
            final VariantContext mergedVC = emitRefConf ? ReferenceConfidenceUtils.addNonRefSymbolicAllele(merged) : merged;

            // converting haplotype likelihoods to allele likelihoods
            final Map<Allele, List<Haplotype>> alleleMapper = AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, haplotypes, true);
            final AlleleLikelihoods<Fragment, Allele> logLikelihoods = logFragmentLikelihoods.marginalize(alleleMapper);
            final SimpleInterval variantCallingRelevantFragmentOverlap = new SimpleInterval(mergedVC).expandWithinContig(MTAC.informativeReadOverlapMargin, header.getSequenceDictionary());
            logLikelihoods.retainEvidence(variantCallingRelevantFragmentOverlap::overlaps);

            if (emitRefConf) {
                logLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }
            final List<LikelihoodMatrix<Fragment, Allele>> tumorMatrices = IntStream.range(0, logLikelihoods.numberOfSamples())
                    .filter(n -> !normalSamples.contains(logLikelihoods.getSample(n)))
                    .mapToObj(logLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final AlleleList<Allele> alleleList = tumorMatrices.get(0);
            final LikelihoodMatrix<Fragment, Allele> logTumorMatrix = combinedLikelihoodMatrix(tumorMatrices, alleleList);
            final PerAlleleCollection<Double> tumorLogOdds = somaticLogOdds(logTumorMatrix);

            final List<LikelihoodMatrix<Fragment, Allele>> normalMatrices = IntStream.range(0, logLikelihoods.numberOfSamples())
                    .filter(n -> normalSamples.contains(logLikelihoods.getSample(n)))
                    .mapToObj(logLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final LikelihoodMatrix<Fragment, Allele> logNormalMatrix = combinedLikelihoodMatrix(normalMatrices, alleleList);
            final PerAlleleCollection<Double> normalLogOdds = diploidAltLogOdds(logNormalMatrix);
            final PerAlleleCollection<Double> normalArtifactLogOdds = somaticLogOdds(logNormalMatrix);


            final Set<Allele> forcedAlleles = AssemblyBasedCallerUtils.allelesConsistentWithGivenAlleles(givenAlleles, mergedVC);
            final List<Allele> tumorAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> forcedAlleles.contains(allele) || tumorLogOdds.getAlt(allele) > MTAC.getEmissionLogOdds())
                    .collect(Collectors.toList());

            final List<Allele> allelesToGenotype = tumorAltAlleles.stream()
                    .filter(allele -> forcedAlleles.contains(allele) || !hasNormal || MTAC.genotypeGermlineSites || normalLogOdds.getAlt(allele) > MathUtils.log10ToLog(MTAC.normalLog10Odds))
                    .toList();

            // record somatic alleles for later use in the Event Count annotation
            // note that in tumor-only calling it does not attempt to detect germline events
            mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> tumorLogOdds.getAlt(allele) > MTAC.getEmissionLogOdds())
                    .filter(allele -> !hasNormal || normalLogOdds.getAlt(allele) > MathUtils.log10ToLog(MTAC.normalLog10Odds))
                    .map(allele -> new Event(mergedVC.getContig(), mergedVC.getStart(), mergedVC.getReference(), allele))
                    .forEach(potentialSomaticEventsInRegion::add);

            // if every alt allele is germline, skip this variant.  However, if some alt alleles are germline and others
            // are not we emit them all so that the filtering engine can see them
            if (allelesToGenotype.isEmpty()) {
                continue;
            }

            final List<Allele> allAllelesToEmit = ListUtils.union(Arrays.asList(mergedVC.getReference()), tumorAltAlleles);


            final Map<String, Object> negativeLogPopulationAFAnnotation =
                    getNegativeLogPopulationAFAnnotation(featureContext.getValues(MTAC.germlineResource, loc),
                            allAllelesToEmit, MTAC.getDefaultAlleleFrequency());

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC)
                    .alleles(allAllelesToEmit)
                    .attributes(negativeLogPopulationAFAnnotation)
                    .attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, tumorAltAlleles.stream().mapToDouble(a -> MathUtils.logToLog10(tumorLogOdds.getAlt(a))).toArray());

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_ARTIFACT_LOG_10_ODDS_KEY,
                        Arrays.stream(normalArtifactLogOdds.asDoubleArray(tumorAltAlleles)).map(x->MathUtils.logToLog10(x)).toArray());
                callVcb.attribute(GATKVCFConstants.NORMAL_LOG_10_ODDS_KEY,
                        Arrays.stream(normalLogOdds.asDoubleArray(tumorAltAlleles)).map(MathUtils::logToLog10).toArray());
            }

            if (!featureContext.getValues(MTAC.pon, mergedVC.getStart()).isEmpty()) {
                callVcb.attribute(GATKVCFConstants.IN_PON_KEY, true);
            }

            addGenotypes(logLikelihoods, allAllelesToEmit, callVcb);
            final VariantContext call = callVcb.make();
            final VariantContext trimmedCall = GATKVariantContextUtils.trimAlleles(call, true, true);
            final List<Allele> trimmedAlleles = trimmedCall.getAlleles();
            final List<Allele> untrimmedAlleles = call.getAlleles();
            final Map<Allele, List<Allele>> trimmedToUntrimmedAlleleMap = IntStream.range(0, trimmedCall.getNAlleles()).boxed()
                    .collect(Collectors.toMap(n -> trimmedAlleles.get(n), n -> Arrays.asList(untrimmedAlleles.get(n))));
            final AlleleLikelihoods<Fragment, Allele> trimmedLikelihoods = logLikelihoods.marginalize(trimmedToUntrimmedAlleleMap);

            // AlleleLikelihoods for annotation only
            final AlleleLikelihoods<GATKRead, Allele> logReadAlleleLikelihoods = logReadLikelihoods.marginalize(alleleMapper);
            logReadAlleleLikelihoods.retainEvidence(variantCallingRelevantFragmentOverlap::overlaps);

            if (emitRefConf) {
                logReadAlleleLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }

            final AlleleLikelihoods<GATKRead, Allele> trimmedLikelihoodsForAnnotation = logReadAlleleLikelihoods.marginalize(trimmedToUntrimmedAlleleMap);


            final VariantContext annotatedCall =  annotationEngine.annotateContext(trimmedCall, featureContext, referenceContext,
                    trimmedLikelihoodsForAnnotation, Optional.of(trimmedLikelihoods), Optional.of(logFragmentLikelihoods), Optional.empty(), a -> true);
            if(withBamOut) {
                AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(trimmedCall, trimmedLikelihoods, Fragment::getReads);
            }

            final Optional<List<VariantContext>> truthVCs = MTAC.mutect3TrainingTruth == null ? Optional.empty() :
                    Optional.of(featureContext.getValues(MTAC.mutect3TrainingTruth, mergedVC.getStart()));
            mutect3DatasetEngine.ifPresent(engine -> engine.addData(referenceContext, annotatedCall, truthVCs,
                    trimmedLikelihoodsForAnnotation, logFragmentLikelihoods, logLikelihoods, MTAC.mutect3DatasetMode));

            call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            returnCalls.add( annotatedCall );
        }

        final List<VariantContext> outputCalls = AssemblyBasedCallerUtils.phaseCalls(returnCalls, calledHaplotypes);
        final int eventCount = outputCalls.size();

        // calculate the number of somatic events in the best haplotype of each called variant
        final Map<Haplotype, MutableInt> haplotypeSupportCounts = logReadLikelihoods.alleles().stream()
                .collect(Collectors.toMap(hap -> hap, label -> new MutableInt(0)));
        logReadLikelihoods.bestAllelesBreakingTies()
                .forEach(bestHaplotype -> haplotypeSupportCounts.get(bestHaplotype.allele).increment());

        final Map<Event, List<Haplotype>> haplotypesByEvent= new HashMap<>();
        for (final Haplotype haplotype : logReadLikelihoods.alleles()) {
            for (final Event event : haplotype.getEventMap().getEvents()) {
                haplotypesByEvent.computeIfAbsent(event, e -> new ArrayList<>()).add(haplotype);
            }
        }
        final Map<VariantContext, List<Integer>> eventCountAnnotations = new HashMap<>();
        for (final VariantContext outputCall : outputCalls) {
            for (final Allele allele : outputCall.getAlternateAlleles()) {
                // note: this creates the minimal representation behind the scenes
                final Event event = new Event(outputCall.getContig(), outputCall.getStart(), outputCall.getReference(), allele);
                // haplotypesByEvent contains every *assembled* event, including events injected into the original assembly,
                // but there are some modes where we genotype events that were never in an assembly graph, in which case
                // this annotation is irrelevant
                if (haplotypesByEvent.containsKey(event)) {
                    final Haplotype bestHaplotype = haplotypesByEvent.get(event).stream()
                            .sorted(Comparator.comparingInt(h -> haplotypeSupportCounts.getOrDefault(h, new MutableInt(0)).intValue()).reversed())
                            .findFirst().get();

                    eventCountAnnotations.computeIfAbsent(outputCall, vc -> new ArrayList<>())
                            .add((int) bestHaplotype.getEventMap().getEvents().stream().filter(potentialSomaticEventsInRegion::contains).count());
                }
            }
        }
        final List<VariantContext> outputCallsWithEventCountAnnotation = outputCalls.stream()
                .map(vc -> new VariantContextBuilder(vc)
                        .attribute(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCountAnnotations.get(vc))
                        .attribute(GATKVCFConstants.EVENT_COUNT_IN_REGION_KEY, potentialSomaticEventsInRegion.size()).make())
                .collect(Collectors.toList());
        return new CalledHaplotypes(outputCallsWithEventCountAnnotation, calledHaplotypes);
    }

    public double[] makePriorPseudocounts(final int numAlleles) {
        return new IndexRange(0, numAlleles).mapToDouble(n -> n == 0 ? refPseudocount : altPseudocount);
    }

    // compute the likelihoods that each allele is contained at some allele fraction in the sample
    protected <EVIDENCE extends Locatable> PerAlleleCollection<Double> somaticLogOdds(final LikelihoodMatrix<EVIDENCE, Allele> logMatrix) {
        final int alleleListEnd = logMatrix.alleles().size()-1;
        if (logMatrix.alleles().contains(Allele.NON_REF_ALLELE) && !(logMatrix.alleles().get(alleleListEnd).equals(Allele.NON_REF_ALLELE))) {
            throw new IllegalStateException("<NON_REF> must be last in the allele list.");
        }

        final double logEvidenceWithAllAlleles = logMatrix.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(getAsRealMatrix(logMatrix), makePriorPseudocounts(logMatrix.numberOfAlleles()));

        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int refIndex = getRefIndex(logMatrix);
        IntStream.range(0, logMatrix.numberOfAlleles()).filter(a -> a != refIndex).forEach( a -> {
            final Allele allele = logMatrix.getAllele(a);
            final LikelihoodMatrix<EVIDENCE, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(logMatrix, allele);
            final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.evidenceCount() == 0 ? 0 :
                    SomaticLikelihoodsEngine.logEvidence(getAsRealMatrix(logMatrixWithoutThisAllele), makePriorPseudocounts(logMatrixWithoutThisAllele.numberOfAlleles()));
            lods.setAlt(allele, logEvidenceWithAllAlleles - logEvidenceWithoutThisAllele);
        });
        return lods;
    }

    private <EVIDENCE extends Locatable> void addGenotypes(final AlleleLikelihoods<EVIDENCE, Allele> logLikelihoods,
                              final List<Allele> allelesToEmit,
                              final VariantContextBuilder callVcb) {
        final List<Genotype> genotypes = IntStream.range(0, logLikelihoods.numberOfSamples()).mapToObj(n -> {
            final String sample = logLikelihoods.getSample(n);
            final LikelihoodMatrix<EVIDENCE, Allele> logMatrix = new SubsettedLikelihoodMatrix<>(logLikelihoods.sampleMatrix(n), allelesToEmit);
            final double[] alleleCounts = getEffectiveCounts(logMatrix);
            final double[] flatPriorPseudocounts = new IndexRange(0, logMatrix.numberOfAlleles()).mapToDouble(a -> 1);
            final double[] alleleFractionsPosterior = logMatrix.evidenceCount() == 0 ? flatPriorPseudocounts :
                    SomaticLikelihoodsEngine.alleleFractionsPosterior(getAsRealMatrix(logMatrix), flatPriorPseudocounts);
            final double[] tumorAlleleFractionsMean = MathUtils.normalizeSumToOne(alleleFractionsPosterior);

            // TODO: We shouldn't always assume that the genotype in the normal is hom ref
            final Allele ref = logMatrix.getAllele(getRefIndex(logMatrix));
            return new GenotypeBuilder(sample, normalSamples.contains(sample) ? Collections.nCopies(2, ref) : logMatrix.alleles())
                    .AD(Arrays.stream(alleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, Arrays.copyOfRange(tumorAlleleFractionsMean, 1, tumorAlleleFractionsMean.length))
                    .make();
        }).collect(Collectors.toList());

        callVcb.genotypes(genotypes);
    }

    private static <EVIDENCE> double[] getEffectiveCounts(final LikelihoodMatrix<EVIDENCE, Allele> logLikelihoodMatrix) {
        if (logLikelihoodMatrix.evidenceCount() == 0) {
            return new double[logLikelihoodMatrix.numberOfAlleles()]; // zero counts for each allele
        }
        final RealMatrix logLikelihoods = getAsRealMatrix(logLikelihoodMatrix);
        return MathUtils.sumArrayFunction(0, logLikelihoods.getColumnDimension(),
                read -> NaturalLogUtils.normalizeFromLogToLinearSpace(logLikelihoods.getColumn(read)));
    }

    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * @param matrix a matrix of log likelihoods
     */
    private <EVIDENCE extends Locatable> PerAlleleCollection<Double> diploidAltLogOdds(final LikelihoodMatrix<EVIDENCE, Allele> matrix) {
        final int refIndex = getRefIndex(matrix);
        final int numReads = matrix.evidenceCount();
        final double homRefLogLikelihood = new IndexRange(0, numReads).sum(r -> matrix.get(refIndex,r));

        final PerAlleleCollection<Double> result = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        // hom ref likelihood for the ref allele, het likelihood for alt alleles
        IntStream.range(0, matrix.numberOfAlleles()).filter(a -> a != refIndex)
                .forEach( a -> {
                    final double hetLogLikelihood = new IndexRange(0, numReads)
                            .sum(r -> NaturalLogUtils.logSumExp(matrix.get(refIndex, r), matrix.get(a,r)) + NaturalLogUtils.LOG_ONE_HALF);
                    result.setAlt(matrix.getAllele(a), homRefLogLikelihood - hetLogLikelihood);
                });
        return result;
    }

    private <EVIDENCE> int getRefIndex(LikelihoodMatrix<EVIDENCE, Allele> matrix) {
        final OptionalInt optionalRefIndex = IntStream.range(0, matrix.numberOfAlleles()).filter(a -> matrix.getAllele(a).isReference()).findFirst();
        Utils.validateArg(optionalRefIndex.isPresent(), "No ref allele found in likelihoods");
        return optionalRefIndex.getAsInt();
    }

    //convert a likelihood matrix of alleles x reads into a RealMatrix
    public static <EVIDENCE> RealMatrix getAsRealMatrix(final LikelihoodMatrix<EVIDENCE, Allele> matrix) {
        final RealMatrix result = new Array2DRowRealMatrix(matrix.numberOfAlleles(), matrix.evidenceCount());
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return matrix.get(row, column);
            }
        });
        return result;
    }

    public static <EVIDENCE extends Locatable> LikelihoodMatrix<EVIDENCE, Allele> combinedLikelihoodMatrix(final List<LikelihoodMatrix<EVIDENCE, Allele>> matrices, final AlleleList<Allele> alleleList) {
        final List<EVIDENCE> reads = matrices.stream().flatMap(m -> m.evidence().stream()).collect(Collectors.toList());
        final AlleleLikelihoods<EVIDENCE, Allele> combinedLikelihoods = new AlleleLikelihoods<>(SampleList.singletonSampleList("COMBINED"), alleleList, ImmutableMap.of("COMBINED", reads));

        int combinedReadIndex = 0;
        final LikelihoodMatrix<EVIDENCE, Allele> result = combinedLikelihoods.sampleMatrix(0);
        final int alleleCount = result.numberOfAlleles();
        for (final LikelihoodMatrix<EVIDENCE, Allele> matrix : matrices) {
            final int readCount = matrix.evidenceCount();
            for (int r = 0; r < readCount; r++) {
                for (int a = 0; a < alleleCount; a++) {
                    result.set(a, combinedReadIndex, matrix.get(a, r));
                }
                combinedReadIndex++;
            }
        }
        return result;
    }

    private <E> Optional<E> getForNormal(final Supplier<E> supplier) {
        return hasNormal ? Optional.of(supplier.get()) : Optional.empty();
    }

    /**
     *
     * @param germlineResourceVariants  Germline resource variants from which our AF INFO field annotation is drawn
     * @param allAlleles    Every emitted allele, with the reference allele first.  Only alt alleles are annotated, but we
     *                      need the ref allele in case the germline resource has a more or less parsimoniuous representation
     *                      For example, eg ref = A, alt = C; germline ref = AT, germline alt = CT
     * @param afOfAllelesNotInGermlineResource  default value of germline AF annotation
     * @return
     */
    private static Map<String, Object> getNegativeLogPopulationAFAnnotation(List<VariantContext> germlineResourceVariants,
                                                                            final List<Allele> allAlleles,
                                                                            final double afOfAllelesNotInGermlineResource) {
        final double[] populationAlleleFrequencies = getGermlineAltAlleleFrequencies(allAlleles, germlineResourceVariants, afOfAllelesNotInGermlineResource);

        return ImmutableMap.of(GATKVCFConstants.POPULATION_AF_KEY, MathUtils.applyToArray(populationAlleleFrequencies, x -> - Math.log10(x)));
    }

    /**
     *
     * @param allAlleles    Every emitted allele, with the reference allele first.  Only alt alleles are annotated, but we
     *      *               need the ref allele in case the germline resource has a more or less parsimonious representation
     *      *               For example, eg ref = A, alt = C; germline ref = AT, germline alt = CT
     * @param germlineVCs    Germline resource variant contexts from which AF INFO field is drawn
     * @param afOfAllelesNotInGermlineResource  Default value of population AF annotation
     * @return
     */
    @VisibleForTesting
    static double[] getGermlineAltAlleleFrequencies(final List<Allele> allAlleles, final List<VariantContext> germlineVCs, final double afOfAllelesNotInGermlineResource) {
        Utils.validateArg(!allAlleles.isEmpty(), "allAlleles are empty -- there is not even a reference allele.");
        final Map<Allele, Double> alleleFrequencies = new HashMap<>();
        allAlleles.forEach(a -> alleleFrequencies.put(a, afOfAllelesNotInGermlineResource));    // initialize everything to the default

        // look through every germline resource variant context at this locus and fill in the AFs
        for (final VariantContext germlineVC : germlineVCs) {
            if (! germlineVC.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                logger.warn("Germline resource variant at " + germlineVC.getContig() + ":" + germlineVC.getStart() +" missing AF attribute");
            }

            List<OptionalInt> germlineIndices = GATKVariantContextUtils.alleleIndices(allAlleles, germlineVC.getAlleles());
            final List<Double> germlineAltAFs = Mutect2Engine.getAttributeAsDoubleList(germlineVC, VCFConstants.ALLELE_FREQUENCY_KEY, afOfAllelesNotInGermlineResource);

            if (germlineAltAFs.size() == (germlineVC.getNAlleles() - 1)) {  // skip VCs with a bad AF field that got parsed as a wrong-length list
                for (int alleleIndex = 1; alleleIndex < allAlleles.size(); alleleIndex++) { // start at 1 to skip the reference, which doesn't have an AF annotation
                    final Allele allele = allAlleles.get(alleleIndex);
                    // note the -1 since germlineAltAFs do not include ref
                    germlineIndices.get(alleleIndex).ifPresent(germlineIndex -> alleleFrequencies.put(allele, germlineAltAFs.get(germlineIndex - 1)));
                }
            }
        }

        return allAlleles.stream().skip(1).mapToDouble(alleleFrequencies::get).toArray();   // skip the reference allele
    }

    /**
     * This method must be called when the client is done with genotyping.
     * It closes any open resources.
     */
    @Override
    public void close() {
        mutect3DatasetEngine.ifPresent(engine -> engine.close());
    }
}
