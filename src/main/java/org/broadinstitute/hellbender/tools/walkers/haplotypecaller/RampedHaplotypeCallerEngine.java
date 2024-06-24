package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;

/**
 * This is a specialized haplotype caller engine, designed to allow for breaking the monolithic haplotype
 * caller process into smaller discrete steps. This version allows for stopping the process at after (or before)
 * a specific step and saving a state file (an off-ramp) that can later be used to restart the process from the
 * same step (an on-ramp).
 *
 * The major steps of the haplotype caller are: assembly, filtering and genotyping
 *
 * At this off ramps are implemented before and after the assembler and before the filtering
 * On ramps are implemented after the assembler and after the filtering.
 */
public class RampedHaplotypeCallerEngine extends HaplotypeCallerEngine {

    private static final Logger logger = LogManager.getLogger(RampedHaplotypeCallerEngine.class);

    private final RampedHaplotypeCallerArgumentCollection rpArgs;

    // on-off ramp, if present
    protected PreFilterOffRamp preFilterOffRamp = null;
    protected PostFilterOnRamp postFilterOnRamp = null;
    protected AssemblerOffRamp preAssemblerOffRamp = null;
    protected AssemblerOffRamp postAssemblerOffRamp = null;
    protected PostAssemblerOnRamp postAssemblerOnRamp = null;

    // haplotype caller phases, as consumers
    final private List<Consumer<CallRegionContext>>   phases =
            Arrays.asList(
                    p -> prepare(p),
                    p -> assemble(p),
                    p -> computeReadLikelihoods(p),
                    p -> uncollapse(p),
                    p -> filter(p),
                    p -> genotype(p)
            );

    public RampedHaplotypeCallerEngine(final HaplotypeCallerArgumentCollection hcArgs, AssemblyRegionArgumentCollection assemblyRegionArgs, boolean createBamOutIndex,
                                       boolean createBamOutMD5, final SAMFileHeader readsHeader,
                                       CachingIndexedFastaSequenceFile referenceReader, VariantAnnotatorEngine annotationEngine,
                                       RampedHaplotypeCallerArgumentCollection rpArgs) {

        super(hcArgs, assemblyRegionArgs, createBamOutIndex,
                createBamOutMD5, readsHeader,
                referenceReader, annotationEngine);
        this.rpArgs = rpArgs;

        buildRamps();
    }

    /**
     * Shutdown this HC engine, closing resources as appropriate
     */
    public void shutdown() {
        super.shutdown();
        tearRamps();
    }

    public void buildRamps() {
        try {

            if ( rpArgs.offRampType != null && rpArgs.offRampType != RampedHaplotypeCallerArgumentCollection.OffRampTypeEnum.NONE) {
                if ( rpArgs.offRampFile == null )
                    throw new RuntimeException("rampFile must be specified");

                // create ramp
                switch ( rpArgs.offRampType ) {
                    case NONE:
                        break;
                    case PRE_FILTER_OFF:
                        preFilterOffRamp = new PreFilterOffRamp(rpArgs.offRampFile);
                        break;
                    case POST_ASSEMBLER_OFF:
                        postAssemblerOffRamp = new AssemblerOffRamp(rpArgs.offRampFile);
                        break;
                    case PRE_ASSEMBLER_OFF:
                        preAssemblerOffRamp = new AssemblerOffRamp(rpArgs.offRampFile);
                        break;
                }
            }

            if ( rpArgs.onRampType != null && rpArgs.onRampType != RampedHaplotypeCallerArgumentCollection.OnRampTypeEnum.NONE ) {
                if ( rpArgs.onRampFile == null )
                    throw new RuntimeException("rampFile must be specified");

                // create ramp
                switch ( rpArgs.onRampType ) {
                    case NONE:
                        break;
                    case POST_FILTER_ON:
                        if ( rpArgs.offRampType == RampedHaplotypeCallerArgumentCollection.OffRampTypeEnum.PRE_FILTER_OFF ) {
                            throw new IllegalArgumentException("an on ramp, if present, must be before the off ramp");
                        }
                        postFilterOnRamp = new PostFilterOnRamp(rpArgs.onRampFile);
                        break;
                    case POST_ASSEMBLER_ON:
                        if ( rpArgs.offRampType != null &&
                                rpArgs.offRampType != RampedHaplotypeCallerArgumentCollection.OffRampTypeEnum.NONE ) {
                            throw new IllegalArgumentException("an on ramp, if present, must be before the off ramp");
                        }
                        postAssemblerOnRamp = new PostAssemblerOnRamp(rpArgs.onRampFile);
                        break;
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void tearRamps() {
        try {
            for (RampBase ramp : Arrays.asList(preFilterOffRamp, postFilterOnRamp,
                    preAssemblerOffRamp, postAssemblerOffRamp, postAssemblerOnRamp) ) {
                if ( ramp != null ) {
                    ramp.close();
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    private class CallRegionContext {

        // params
        final AssemblyRegion region;
        final FeatureContext features;
        final ReferenceContext referenceContext;

        // prepared fields or fields we're able to regenerate
        List<VariantContext> VCpriors;
        List<Event> givenAlleles;
        LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsing;

        // assembly results
        AssemblyResultSet assemblyResult;
        Optional<AssemblyRegion> nonVariantLeftFlankRegion;
        Optional<AssemblyRegion> nonVariantRightFlankRegion;

        // computeReadLikelihoods and filter results fields
        AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods;
        Set<Integer> suspiciousLocations;

        // genotyping result
        List<VariantContext> regionVariants;

        CallRegionContext(final AssemblyRegion region, final FeatureContext features, final ReferenceContext referenceContext) {

            this.region = region;
            this.features = features;
            this.referenceContext = referenceContext;
        }

        CallRegionContext(final CallRegionContext other) {
            this.region = other.region;
            this.features = other.features;
            this.referenceContext = other.referenceContext;
            this.VCpriors = other.VCpriors;
            this.givenAlleles = other.givenAlleles;
            this.haplotypeCollapsing = other.haplotypeCollapsing;
            this.assemblyResult = other.assemblyResult;
            this.nonVariantLeftFlankRegion = other.nonVariantLeftFlankRegion;
            this.nonVariantRightFlankRegion = other.nonVariantRightFlankRegion;
            this.readLikelihoods = other.readLikelihoods;
            this.regionVariants = other.regionVariants;
        }
    }

    /**
     * Generate variant calls for an assembly region
     *
     * @param region   region to assemble and perform variant calling on
     * @param features Features overlapping the assembly region
     * @return List of variants discovered in the region (may be empty)
     */
    @Override
    public List<VariantContext> callRegion(final AssemblyRegion region, final FeatureContext features, final ReferenceContext referenceContext) {

        // dump reads upon entry
        RampUtils.logReads(rpArgs.rampsDebugReads, "callRegion entered: " + region, region.getReads());

        // create initial context
        CallRegionContext context = new CallRegionContext(region, features, referenceContext);

        // execute stages
        final Iterator<Consumer<CallRegionContext>> iter = phases.iterator();
        while ( iter.hasNext() && context.regionVariants == null ) {
            iter.next().accept(context);
        }

        // return variants
        return context.regionVariants;
    }

    private void prepare(final CallRegionContext context) {

        // no need for this step?
        if ( context.regionVariants != null ) {
            return;
        }

        if (hcArgs.justDetermineActiveRegions) {
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            context.regionVariants = NO_CALLS;
            return;
        }

        context.VCpriors = new ArrayList<>();
        if (hcArgs.standardArgs.genotypeArgs.supportVariants != null) {
            context.features.getValues(hcArgs.standardArgs.genotypeArgs.supportVariants).stream().forEach(context.VCpriors::add);
        }

        if (hcArgs.sampleNameToUse != null) {
            removeReadsFromAllSamplesExcept(hcArgs.sampleNameToUse, context.region);
        }

        if (!context.region.isActive()) {
            // Not active so nothing to do!
            context.regionVariants = referenceModelForNoVariation(context.region, true, context.VCpriors);
            return;
        }

        context.givenAlleles = context.features.getValues(hcArgs.alleles).stream()
                .filter(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())
                .flatMap(vc -> GATKVariantContextUtils.splitVariantContextToEvents(vc, false, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false).stream())
                .collect(Collectors.toList());

        if (context.givenAlleles.isEmpty() && context.region.size() == 0) {
            // No reads here so nothing to do!
            context.regionVariants = referenceModelForNoVariation(context.region, true, context.VCpriors);
        }
    }

    private void assemble(final CallRegionContext context) {

        try {
            // no need for this step?
            if (context.regionVariants != null) {
                return;
            }

            if (assemblyDebugOutStream != null) {
                assemblyDebugOutStream.write("\n\n\n\n" + context.region.getSpan() + "\nNumber of reads in region: " + context.region.getReads().size() + "     they are:");
                for (GATKRead read : context.region.getReads()) {
                    assemblyDebugOutStream.write(read.getName() + "   " + read.convertToSAMRecord(context.region.getHeader()).getFlags());
                }
            }

            // get off?
            if (preAssemblerOffRamp != null) {
                preAssemblerOffRamp.add(context.region, "assemblyResult", context.assemblyResult, context.nonVariantLeftFlankRegion, context.nonVariantRightFlankRegion, readsHeader);
                context.regionVariants = NO_CALLS;
                return;
            }

            // assembler on ramp?
            CallRegionContext debugContext = null;
            if (postAssemblerOnRamp != null) {

                // clean part
                context.nonVariantLeftFlankRegion = postAssemblerOnRamp.getOptionalAssemblyRegion(context.region, "assemblyResult.nonVariantLeftFlankRegion", readsHeader);
                context.nonVariantRightFlankRegion = postAssemblerOnRamp.getOptionalAssemblyRegion(context.region, "assemblyResult.nonVariantRightFlankRegion", readsHeader);
                context.assemblyResult = postAssemblerOnRamp.getAssembyResult(context.region, "assemblyResult", readsHeader, hcArgs, logger, aligner, referenceReader);
                if (context.assemblyResult == null) {
                    context.regionVariants = NO_CALLS;
                    return;
                }
                context.haplotypeCollapsing = context.assemblyResult.getHaplotypeCollapsingEngine();
                RampUtils.logReads(rpArgs.rampsDebugReads, "onramp: reads before trimming", context.assemblyResult.getRegionForGenotyping().getReads());

                RampUtils.logReads(rpArgs.rampsDebugReads, "onramp: BEFORE untrimmedAssemblyResult reads", context.region.getReads());
                final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(context.region, hcArgs, readsHeader, samplesList, logger, referenceReader, assemblyEngine, aligner,
                        !hcArgs.doNotCorrectOverlappingBaseQualities, hcArgs.fbargs, postFilterOnRamp != null);
                RampUtils.logReads(rpArgs.rampsDebugReads, "onramp: AFTER untrimmedAssemblyResult reads", context.region.getReads());
                context.assemblyResult.setRegionForGenotyping(untrimmedAssemblyResult.getRegionForGenotyping());


                // restore trimmed reads
                final SortedSet<Event> allVariationEvents = context.assemblyResult.getVariationEvents(hcArgs.maxMnpDistance);
                final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(context.region, allVariationEvents, context.referenceContext);
                try {
                    context.assemblyResult = context.assemblyResult.trimTo(trimmingResult.getVariantRegion());
                } catch (IllegalStateException e) {
                    logger.warn("assembler bailing out on " + context.region + " bacause " + e.getMessage());
                    context.regionVariants = NO_CALLS;
                    return;
                    // TODO: investigate further
                }
                context.haplotypeCollapsing = context.assemblyResult.getHaplotypeCollapsingEngine();
                RampUtils.logReads(rpArgs.rampsDebugReads, "onramp: reads after trimming", context.assemblyResult.getRegionForGenotyping().getReads());

                if (rpArgs.rampsDebugPostAssemblerOn) {
                    debugContext = new CallRegionContext(context);
                    logger.debug("debugContext: " + debugContext);
                } else {
                    return;
                }
            }

            // run the local assembler, getting back a collection of information on how we should proceed
            RampUtils.logReads(rpArgs.rampsDebugReads, "BEFORE untrimmedAssemblyResult reads", context.region.getReads());
            List<VariantContext> forcedPileupAlleles = Collections.emptyList(); // TODO: we currently do not support pileup alleles in RampedHaplotypeCaller, this should be added
            final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(context.region, hcArgs, readsHeader, samplesList, logger, referenceReader, assemblyEngine, aligner,
                    !hcArgs.doNotCorrectOverlappingBaseQualities, hcArgs.fbargs, postFilterOnRamp != null);
            RampUtils.logReads(rpArgs.rampsDebugReads, "AFTER untrimmedAssemblyResult reads", context.region.getReads());
            if (postFilterOnRamp != null) {
                if (!postFilterOnRamp.hasRegion(context.region.getSpan())) {
                    context.regionVariants = NO_CALLS;
                    return;
                }
                for (Haplotype h : postFilterOnRamp.getHaplotypes(context.region, "readLikelihoods.haplotypes")) {
                    if (!h.isReference())
                        untrimmedAssemblyResult.add(h);
                }
            }
            context.haplotypeCollapsing = untrimmedAssemblyResult.getHaplotypeCollapsingEngine();

            if (assemblyDebugOutStream != null) {
                assemblyDebugOutStream.write("\nThere were " + untrimmedAssemblyResult.getHaplotypeList().size() + " haplotypes found. Here they are:");
                for (final String haplotype : untrimmedAssemblyResult.getHaplotypeList().stream().map(Haplotype::toString).sorted().collect(Collectors.toList())) {
                    assemblyDebugOutStream.write(haplotype);
                }
            }

            final SortedSet<Event> allVariationEvents = untrimmedAssemblyResult.getVariationEvents(hcArgs.maxMnpDistance);

            AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(context.region, allVariationEvents, context.referenceContext);

            if (!trimmingResult.isVariationPresent() && !hcArgs.disableOptimizations && postFilterOnRamp == null) {
                context.regionVariants = referenceModelForNoVariation(context.region, false, context.VCpriors);
                return;
            }

            RampUtils.logReads(rpArgs.rampsDebugReads, "untrimmedAssemblyResult reads", untrimmedAssemblyResult.getRegionForGenotyping().getReads());
            context.assemblyResult = untrimmedAssemblyResult.trimTo(
                    postFilterOnRamp == null
                            ? trimmingResult.getVariantRegion()
                            : trimmingResult.getVariantRegion(
                            new SimpleInterval(postFilterOnRamp.getRegion(context.region.getSpan(), "readLikelihoods").getString("r4g")),
                            new SimpleInterval(postFilterOnRamp.getRegion(context.region.getSpan(), "readLikelihoods").getString("r4gp"))));
            context.nonVariantLeftFlankRegion = trimmingResult.nonVariantLeftFlankRegion();
            context.nonVariantRightFlankRegion = trimmingResult.nonVariantRightFlankRegion();
            RampUtils.logReads(rpArgs.rampsDebugReads, "assemblyResult final reads", context.assemblyResult.getRegionForGenotyping().getReads());

            if (rpArgs.rampsDebugPostAssemblerOn && debugContext != null) {
                RampUtils.compareHaplotypes(context.assemblyResult.getHaplotypeList(),
                        debugContext.assemblyResult.getHaplotypeList());

                // DEBUG:: makes it fail here?
                //context.assemblyResult.getRegionForGenotyping().clearReads();
                //debugContext.assemblyResult.getRegionForGenotyping().getReads().forEach(r -> context.assemblyResult.getRegionForGenotyping().add(r));
                RampUtils.compareReads(context.assemblyResult.getRegionForGenotyping().getReads(),
                        debugContext.assemblyResult.getRegionForGenotyping().getReads());
            }

            // get off?
            if (postAssemblerOffRamp != null) {
                postAssemblerOffRamp.add(context.region, "assemblyResult", context.assemblyResult, context.nonVariantLeftFlankRegion, context.nonVariantRightFlankRegion, readsHeader);
                context.regionVariants = NO_CALLS;
            }
        } catch (IOException e) {
            throw new GATKException("assmble failed", e);
        }
    }

    private void computeReadLikelihoods(final CallRegionContext context) {

        // no need for this step?
        if ( context.regionVariants != null ) {
            return;
        }

        final AssemblyRegion regionForGenotyping = context.assemblyResult.getRegionForGenotyping();
        final List<GATKRead> readStubs = regionForGenotyping.getReads().stream()
                .filter(r -> r.getLength() < AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING).collect(Collectors.toList());
        regionForGenotyping.removeAll(readStubs);

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (!context.assemblyResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            context.regionVariants = referenceModelForNoVariation(context.region, false, context.VCpriors);
            return;
        }

        // For sure this is not true if gVCF is on.
        if (hcArgs.dontGenotype) {
            context.regionVariants = NO_CALLS; // user requested we not proceed
            return;
        }

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (regionForGenotyping.size() == 0 && !hcArgs.disableOptimizations) {
            // no reads remain after filtering so nothing else to do!
            context.regionVariants = referenceModelForNoVariation(context.region, false, context.VCpriors);
            return;
        }

        // evaluate each sample's reads against all haplotypes
        final Map<String, List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, regionForGenotyping.getReads());

        // Calculate the likelihoods: CPU intensive part.
        context.readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(context.assemblyResult, samplesList, reads, true);

        alleleLikelihoodWriter.ifPresent(
                writer -> writer.writeAlleleLikelihoods(context.readLikelihoods));
    }

    private void uncollapse(final CallRegionContext context) {

        // no need for this step?
        if ( context.regionVariants != null ) {
            return;
        }

        // Realign reads to their best haplotype.
        final Map<GATKRead, GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(context.readLikelihoods, context.assemblyResult.getReferenceHaplotype(), context.assemblyResult.getPaddedReferenceLoc(), aligner, hcArgs.getHaplotypeToReferenceSWParameters());
        context.readLikelihoods.changeEvidence(readRealignments);
        List<Haplotype> haplotypes = context.readLikelihoods.alleles();

        // regenerate refview? use it to uncollapse haplotypes if needed
        final LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsing = (context.haplotypeCollapsing != null) ? context.haplotypeCollapsing : buildHaplotypeCollapsing(context);

        if ( haplotypeCollapsing != null ) {

            haplotypes = haplotypeCollapsing.uncollapseHmersInHaplotypes(haplotypes, false, null);
            logger.debug(String.format("%d haplotypes before uncollapsing", haplotypes.size()));
            Map<Haplotype, List<Haplotype>> identicalHaplotypesMap = LongHomopolymerHaplotypeCollapsingEngine.identicalBySequence(haplotypes);
            context.readLikelihoods.changeAlleles(haplotypes);
            context.readLikelihoods = context.readLikelihoods.marginalize(identicalHaplotypesMap);
            logger.debug(String.format("%d haplotypes after uncollapsing",  context.readLikelihoods.numberOfAlleles()));
        } else {
            logger.debug(String.format("Not performing uncollapsing with %d haplotypes", context.readLikelihoods.numberOfAlleles()));
        }
    }

    private void filter(final CallRegionContext context) {
        final HaplotypeCallerGenotypingEngine localGenotypingEngine = getLocalGenotypingEngine(context.region);

        try {
            // no need for this step?
            if (context.regionVariants != null) {
                return;
            }

            // ramp point: pre-filter off ramp
            if (preFilterOffRamp != null) {

                // if we are here, we are getting off before the filter
                preFilterOffRamp.add(context.region, "readLikelihoods", context.readLikelihoods, context.assemblyResult.getRegionForGenotyping(), context.region);
                context.regionVariants = NO_CALLS;
            } else if (postFilterOnRamp == null) {

                // if we are here, we are running the filter
                context.suspiciousLocations = new HashSet<>();
                if (hcArgs.filterAlleles) {
                    logger.debug("Filtering alleles");
                    AlleleFilteringHC alleleFilter = new AlleleFilteringHC(hcArgs, assemblyDebugOutStream, localGenotypingEngine);
                    //need to update haplotypes to find the alleles
                    EventMap.buildEventMapsForHaplotypes(context.readLikelihoods.alleles(),
                            context.assemblyResult.getFullReferenceWithPadding(),
                            context.assemblyResult.getPaddedReferenceLoc(),
                            hcArgs.assemblerArgs.debugAssembly,
                            hcArgs.maxMnpDistance);
                    context.readLikelihoods = alleleFilter.filterAlleles(context.readLikelihoods, context.assemblyResult.getPaddedReferenceLoc().getStart(), context.suspiciousLocations);

                } else {
                    logger.debug("Not filtering alleles");
                }

                // Reproscess with the HMM if we are in stepwise filtering mode
                if (hcArgs.stepwiseFiltering) {
                    final Map<String, List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, context.assemblyResult.getRegionForGenotyping().getReads());
                    context.readLikelihoods = (likelihoodCalculationEngine).computeReadLikelihoods(
                            context.readLikelihoods.alleles(),
                            context.assemblyResult.getRegionForGenotyping().getHeader(),
                            samplesList, reads, true);
                }


            } else {

                // if here, we are getting on after the filter
                final AssemblyRegion regionForGenotyping = context.assemblyResult.getRegionForGenotyping();
                final Map<String, List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, regionForGenotyping.getReads());
                final Map<String, List<GATKRead>> reads2 = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, getAllRegionReads(context));

                context.readLikelihoods = postFilterOnRamp.getAlleleLikelihoods(context.region, "readLikelihoods", samplesList, reads, reads2);
            }

            for (int sampleIndex = 0; sampleIndex < context.readLikelihoods.samples().size(); sampleIndex++)
                RampUtils.logReads(rpArgs.rampsDebugReads, "filter (END): " + sampleIndex, context.readLikelihoods.sampleEvidence(sampleIndex));
        } catch (IOException e) {
            throw new GATKException("filter failed", e);
        }
    }

    private void genotype(final CallRegionContext context) {
        final HaplotypeCallerGenotypingEngine localGenotypingEngine = getLocalGenotypingEngine(context.region);

        // no need for this step?
        if ( context.regionVariants != null ) {
            return;
        }

        final AssemblyRegion regionForGenotyping = context.assemblyResult.getRegionForGenotyping();

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);
        final Map<String, List<GATKRead>> perSampleFilteredReadList = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, filteredReads);

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]
        List<Haplotype>     haplotypes = context.readLikelihoods.alleles();
        final CalledHaplotypes calledHaplotypes = localGenotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                context.readLikelihoods,
                perSampleFilteredReadList,
                context.assemblyResult.getFullReferenceWithPadding(),
                context.assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getSpan(),
                context.features,
                context.givenAlleles,
                emitReferenceConfidence(),
                hcArgs.maxMnpDistance,
                readsHeader,
                haplotypeBAMWriter.isPresent(),
                context.suspiciousLocations != null ? context.suspiciousLocations : new HashSet<>(),
                context.readLikelihoods);

        if (haplotypeBAMWriter.isPresent()) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (hcArgs.disableOptimizations) {
                calledHaplotypeSet.add(context.assemblyResult.getReferenceHaplotype());
            }
            haplotypeBAMWriter.get().writeReadsAlignedToHaplotypes(haplotypes, context.assemblyResult.getPaddedReferenceLoc(), haplotypes,
                    calledHaplotypeSet, context.readLikelihoods, regionForGenotyping.getSpan());
        }

        if (hcArgs.assemblerArgs.debugAssembly) {
            logger.info("----------------------------------------------------------------------------------");
        }

        if (emitReferenceConfidence()) {
            if (!containsCalls(calledHaplotypes)) {
                // no called all of the potential haplotypes
                context.regionVariants = referenceModelForNoVariation(context.region, false, context.VCpriors);
            } else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section, then variant-containing section, then right flank
                context.nonVariantLeftFlankRegion.ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank, false, context.VCpriors)));

                result.addAll(referenceConfidenceModel.calculateRefConfidence(context.assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), context.assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        context.readLikelihoods, localGenotypingEngine.getPloidyModel(), calledHaplotypes.getCalls(), hcArgs.standardArgs.genotypeArgs.supportVariants != null,
                        context.VCpriors));

                context.nonVariantRightFlankRegion.ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank, false, context.VCpriors)));

                context.regionVariants = result;
            }
        } else {
            //TODO this should be updated once reducible annotations are handled properly.
            context.regionVariants = calledHaplotypes.getCalls()
                    .stream()
                    .map(RMSMappingQuality.getInstance()::finalizeRawMQ)
                    .collect(Collectors.toList());
        }
    }

    private Collection<GATKRead> getAllRegionReads(final CallRegionContext context) {
        return context.region.getReads();
    }

    private LongHomopolymerHaplotypeCollapsingEngine buildHaplotypeCollapsing(final CallRegionContext context) {

        // estblish reference mapper, if needed
        final byte[] fullReferenceWithPadding = context.region.getAssemblyRegionReference(referenceReader, AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY);
        final SimpleInterval paddedReferenceLoc = AssemblyBasedCallerUtils.getPaddedReferenceLoc(context.region, AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY, referenceReader);
        final Haplotype refHaplotype = AssemblyBasedCallerUtils.createReferenceHaplotype(context.region, paddedReferenceLoc, referenceReader);

        final LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsing =
                (hcArgs.flowAssemblyCollapseHKerSize > 0
                        && LongHomopolymerHaplotypeCollapsingEngine.needsCollapsing(refHaplotype.getBases(), hcArgs.flowAssemblyCollapseHKerSize, logger))
                        ? new LongHomopolymerHaplotypeCollapsingEngine(hcArgs.flowAssemblyCollapseHKerSize, hcArgs.flowAssemblyCollapsePartialMode, fullReferenceWithPadding,
                        paddedReferenceLoc, logger, hcArgs.assemblerArgs.debugAssembly, aligner, hcArgs.getHaplotypeToReferenceSWParameters())
                        : null;
        if ( haplotypeCollapsing != null ) {
            logger.debug("deploying refView on " + paddedReferenceLoc + ", region: " + context.region);
        }

        return haplotypeCollapsing;
    }
}
