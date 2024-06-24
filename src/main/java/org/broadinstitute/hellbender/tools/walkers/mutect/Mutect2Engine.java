package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import org.apache.commons.lang3.mutable.MutableLong;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.HomogeneousPloidyModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2CountsCollector;
import org.broadinstitute.hellbender.transformers.PalindromeArtifactClipReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 9/15/16.
 */
public final class Mutect2Engine implements AssemblyRegionEvaluator, AutoCloseable {

    private static final List<String> STANDARD_MUTECT_INFO_FIELDS = Arrays.asList(GATKVCFConstants.NORMAL_LOG_10_ODDS_KEY, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, GATKVCFConstants.NORMAL_ARTIFACT_LOG_10_ODDS_KEY,
            GATKVCFConstants.EVENT_COUNT_IN_REGION_KEY, GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, GATKVCFConstants.IN_PON_KEY, GATKVCFConstants.POPULATION_AF_KEY,
            GATKVCFConstants.GERMLINE_QUAL_KEY, GATKVCFConstants.CONTAMINATION_QUAL_KEY, GATKVCFConstants.SEQUENCING_QUAL_KEY,
            GATKVCFConstants.POLYMERASE_SLIPPAGE_QUAL_KEY, GATKVCFConstants.READ_ORIENTATION_QUAL_KEY,
            GATKVCFConstants.STRAND_QUAL_KEY, GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, GATKVCFConstants.N_COUNT_KEY, GATKVCFConstants.AS_UNIQUE_ALT_READ_SET_COUNT_KEY);
    private static final String MUTECT_VERSION = "2.2";

    public static final String TUMOR_SAMPLE_KEY_IN_VCF_HEADER = "tumor_sample";
    public static final String NORMAL_SAMPLE_KEY_IN_VCF_HEADER = "normal_sample";

    private static final Logger logger = LogManager.getLogger(Mutect2Engine.class);
    private final static List<VariantContext> NO_CALLS = Collections.emptyList();
    public static final int INDEL_START_QUAL = 30;
    public static final int INDEL_CONTINUATION_QUAL = 10;
    public static final double MAX_ALT_FRACTION_IN_NORMAL = 0.3;
    public static final int MAX_NORMAL_QUAL_SUM = 100;
    public static final int MIN_PALINDROME_SIZE = 5;

    public static final int HUGE_FRAGMENT_LENGTH = 1_000_000;
    public static final int MIN_TAIL_QUALITY = 9;

    final private M2ArgumentCollection MTAC;
    private SAMFileHeader header;
    private SAMSequenceDictionary sequenceDictionary;
    private final int minCallableDepth;
    public static final String CALLABLE_SITES_NAME = "callable";

    private static final int MIN_READ_LENGTH = 30;
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;
    public static final int MINIMUM_BASE_QUALITY = 6;   // for active region determination

    private final SampleList samplesList;
    private final Set<String> normalSamples;

    private final boolean forceCallingAllelesPresent;

    private CachingIndexedFastaSequenceFile referenceReader;
    private ReadThreadingAssembler assemblyEngine;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private SomaticGenotypingEngine genotypingEngine;
    private Optional<HaplotypeBAMWriter> haplotypeBAMWriter;
    private Optional<VariantContextWriter> assembledEventMapVcfOutputWriter;
    private Optional<PriorityQueue<VariantContext>> assembledEventMapVariants;
    private VariantAnnotatorEngine annotationEngine;
    private final SmithWatermanAligner aligner;
    private final AssemblyRegionTrimmer trimmer;
    private SomaticReferenceConfidenceModel referenceConfidenceModel = null;

    private final MutableLong callableSites = new MutableLong(0);

    private final Optional<F1R2CountsCollector> f1R2CountsCollector;

    private PileupQualBuffer tumorPileupQualBuffer;
    private PileupQualBuffer normalPileupQualBuffer;

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     *  @param MTAC command-line arguments for the HaplotypeCaller
     * @param assemblyRegionArgs
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param header header for the reads
     * @param referenceSpec reference specifier for the reference
     * @param annotatorEngine annotator engine built with desired annotations
     */
    public Mutect2Engine(final M2ArgumentCollection MTAC, AssemblyRegionArgumentCollection assemblyRegionArgs,
                         final boolean createBamOutIndex, final boolean createBamOutMD5, final SAMFileHeader header,
                         final SAMSequenceDictionary sequenceDictionary, final GATKPath referenceSpec, final VariantAnnotatorEngine annotatorEngine) {
        this.MTAC = Utils.nonNull(MTAC);
        this.header = Utils.nonNull(header);
        this.sequenceDictionary = sequenceDictionary;
        minCallableDepth = MTAC.callableDepth;
        referenceReader = ReferenceUtils.createReferenceReader(Utils.nonNull(referenceSpec));
        aligner = SmithWatermanAligner.getAligner(MTAC.smithWatermanImplementation);
        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(header)));

        tumorPileupQualBuffer = new PileupQualBuffer(MTAC.activeRegionMultipleSubstitutionBaseQualCorrection);
        normalPileupQualBuffer = new PileupQualBuffer(MTAC.activeRegionMultipleSubstitutionBaseQualCorrection);

        // optimize set operations for the common cases of no normal and one normal
        if (MTAC.normalSamples.isEmpty()) {
            normalSamples = Collections.emptySet();
        } else if (MTAC.normalSamples.size() == 1) {
            normalSamples = Collections.singleton(decodeSampleNameIfNecessary(MTAC.normalSamples.iterator().next()));
        } else {
            normalSamples = MTAC.normalSamples.stream().map(this::decodeSampleNameIfNecessary).collect(Collectors.toSet());
        }
        normalSamples.forEach(this::checkSampleInBamHeader);

        forceCallingAllelesPresent = MTAC.alleles != null;

        annotationEngine = Utils.nonNull(annotatorEngine);
        assemblyEngine = MTAC.createReadThreadingAssembler();
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs, MTAC.fbargs, true, MTAC.likelihoodArgs.likelihoodEngineImplementation);
        genotypingEngine = new SomaticGenotypingEngine(MTAC, normalSamples, annotationEngine, header, sequenceDictionary);
        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(MTAC, createBamOutIndex, createBamOutMD5, header);
        trimmer = new AssemblyRegionTrimmer(assemblyRegionArgs, header.getSequenceDictionary());

        boolean isFlowBased = (MTAC.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBased);

        referenceConfidenceModel = new SomaticReferenceConfidenceModel(samplesList, header, 0,
                MTAC.minAF, MTAC.refModelDelQual, !MTAC.overrideSoftclipFragmentCheck, isFlowBased);  //TODO: do something classier with the indel size arg
        final List<String> tumorSamples = ReadUtils.getSamplesFromHeader(header).stream().filter(this::isTumorSample).collect(Collectors.toList());
        f1R2CountsCollector = MTAC.f1r2TarGz == null ? Optional.empty() : Optional.of(new F1R2CountsCollector(MTAC.f1r2Args, header, MTAC.f1r2TarGz, tumorSamples));
        assembledEventMapVcfOutputWriter = Optional.ofNullable(MTAC.assemblerArgs.debugAssemblyVariantsOut != null ?
                GATKVariantContextUtils.createVCFWriter(
                        new GATKPath(MTAC.assemblerArgs.debugAssemblyVariantsOut).toPath(),
                        header.getSequenceDictionary(),
                        false,
                        Options.DO_NOT_WRITE_GENOTYPES, Options.INDEX_ON_THE_FLY)
                : null);
        assembledEventMapVariants = Optional.ofNullable(MTAC.assemblerArgs.debugAssemblyVariantsOut != null ?
                new PriorityQueue<>(200, new VariantContextComparator(header.getSequenceDictionary())) : null);
        assembledEventMapVcfOutputWriter.ifPresent(writer -> {VCFHeader head = new VCFHeader(); head.getSequenceDictionary(); writer.writeHeader(head);});
    }

    //default M2 read filters.  Cheap ones come first in order to fail fast.
    public static List<ReadFilter> makeStandardMutect2ReadFilters() {
        return Arrays.asList(new MappingQualityReadFilter(READ_QUALITY_FILTER_THRESHOLD),
                ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE,
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO,
                ReadFilterLibrary.MAPPED,
                ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT,
                ReadFilterLibrary.NOT_DUPLICATE,
                ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK,
                ReadFilterLibrary.NON_CHIMERIC_ORIGINAL_ALIGNMENT_READ_FILTER,
                ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT,
                new ReadLengthReadFilter(MIN_READ_LENGTH, Integer.MAX_VALUE),
                ReadFilterLibrary.GOOD_CIGAR,
                new WellformedReadFilter());
    }

    public static ReadTransformer makeStandardMutect2PostFilterReadTransformer(final Path referencePath, boolean clipITRArtifacts) {
        return !clipITRArtifacts ? ReadTransformer.identity() :
                new PalindromeArtifactClipReadTransformer(new ReferenceFileSource(referencePath), MIN_PALINDROME_SIZE);
    }

    /**
     * @return the default set of variant annotations for Mutect2
     */
    public static List<Class<? extends Annotation>> getStandardMutect2AnnotationGroups() {
        return Collections.singletonList(StandardMutectAnnotation.class);
    }

    public void writeHeader(final VariantContextWriter vcfWriter, final Set<VCFHeaderLine> defaultToolHeaderLines) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(new VCFHeaderLine("MutectVersion", MUTECT_VERSION));
        headerInfo.add(new VCFHeaderLine(FilterMutectCalls.FILTERING_STATUS_VCF_KEY, "Warning: unfiltered Mutect 2 calls.  Please run " + FilterMutectCalls.class.getSimpleName() + " to remove false positives."));
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerInfo.addAll(defaultToolHeaderLines);
        STANDARD_MUTECT_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);

        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ALLELE_FRACTION_KEY));
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        headerInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.PHASE_SET_KEY));

        ReadUtils.getSamplesFromHeader(header).stream().filter(this::isTumorSample)
                .forEach(s -> headerInfo.add(new VCFHeaderLine(TUMOR_SAMPLE_KEY_IN_VCF_HEADER, s)));

        normalSamples.forEach(sample -> headerInfo.add(new VCFHeaderLine(NORMAL_SAMPLE_KEY_IN_VCF_HEADER, sample)));
        if (emitReferenceConfidence()) {
            headerInfo.addAll(referenceConfidenceModel.getVCFHeaderLines());
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY));
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo, samplesList.asListOfSamples());
        vcfHeader.setSequenceDictionary(header.getSequenceDictionary());
        vcfWriter.writeHeader(vcfHeader);
    }

    public List<VariantContext> callRegion(final AssemblyRegion originalAssemblyRegion, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        // divide PCR qual by two in order to get the correct total qual when treating paired reads as independent
        AssemblyBasedCallerUtils.cleanOverlappingReadPairs(originalAssemblyRegion.getReads(), samplesList, header,
                false, OptionalInt.of(MTAC.pcrSnvQual /2), OptionalInt.of(MTAC.pcrIndelQual /2));

        if ( !originalAssemblyRegion.isActive() || originalAssemblyRegion.size() == 0 ) {
            return emitReferenceConfidence() ? referenceModelForNoVariation(originalAssemblyRegion) : NO_CALLS;  //TODD: does this need to be finalized?
        }

        removeUnmarkedDuplicates(originalAssemblyRegion);

        final List<Event> givenAlleles = featureContext.getValues(MTAC.alleles).stream()
                .filter(vc -> MTAC.forceCallFiltered || vc.isNotFiltered())
                .flatMap(vc -> GATKVariantContextUtils.splitVariantContextToEvents(vc, false, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false).stream())
                .collect(Collectors.toList());

        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(originalAssemblyRegion, MTAC, header, samplesList, logger, referenceReader, assemblyEngine, aligner, false, MTAC.fbargs, false);
        ReadThreadingAssembler.addAssembledVariantsToEventMapOutput(untrimmedAssemblyResult, assembledEventMapVariants, MTAC.maxMnpDistance, assembledEventMapVcfOutputWriter);
        final LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsing = untrimmedAssemblyResult.getHaplotypeCollapsingEngine();

        final SortedSet<Event> allVariationEvents = untrimmedAssemblyResult.getVariationEvents(MTAC.maxMnpDistance);

        Pair<Set<Event>, Set<Event>> goodAndBadPileupEvents =
                PileupBasedAlleles.goodAndBadPileupEvents(originalAssemblyRegion.getAlignmentData(), MTAC.pileupDetectionArgs, header, MTAC.minBaseQualityScore);
        final Set<Event> goodPileupEvents = goodAndBadPileupEvents.getLeft();
        final Set<Event> badPileupEvents = goodAndBadPileupEvents.getRight();

        goodPileupEvents.forEach(allVariationEvents::add);
        givenAlleles.forEach(allVariationEvents::add);

        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(originalAssemblyRegion, allVariationEvents, referenceContext);
        if (!trimmingResult.isVariationPresent()) {
            return emitReferenceConfidence() ? referenceModelForNoVariation(originalAssemblyRegion) : NO_CALLS;
        }

        AssemblyResultSet assemblyResult = untrimmedAssemblyResult.trimTo(trimmingResult.getVariantRegion());
        assemblyResult.removeHaplotypesWithBadAlleles(MTAC, badPileupEvents);
        assemblyResult.addGivenAlleles(givenAlleles, MTAC.maxMnpDistance, aligner, MTAC.getHaplotypeToReferenceSWParameters());
        assemblyResult.injectPileupEvents(originalAssemblyRegion, MTAC, aligner, goodPileupEvents);

        // we might find out after assembly that the "active" region actually has no variants
        if( ! assemblyResult.isVariationPresent() ) {
            return emitReferenceConfidence() ? referenceModelForNoVariation(originalAssemblyRegion) : NO_CALLS;
        }

        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        removeReadStubs(regionForGenotyping);

        final Map<String,List<GATKRead>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads, true);
        readLikelihoods.switchToNaturalLog();
        final SWParameters readToHaplotypeSWParameters = MTAC.getReadToHaplotypeSWParameters();
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner, readToHaplotypeSWParameters);
        readLikelihoods.changeEvidence(readRealignments);
        List<Haplotype> haplotypes = readLikelihoods.alleles();

        final AlleleLikelihoods<GATKRead, Haplotype> uncollapsedReadLikelihoods;
        if ( haplotypeCollapsing != null ) {

            haplotypes = haplotypeCollapsing.uncollapseHmersInHaplotypes(haplotypes, false, null);
            logger.debug(String.format("%d haplotypes before uncollapsing", haplotypes.size()));
            Map<Haplotype, List<Haplotype>> identicalHaplotypesMap = LongHomopolymerHaplotypeCollapsingEngine.identicalBySequence(haplotypes);
            readLikelihoods.changeAlleles(haplotypes);
            uncollapsedReadLikelihoods = readLikelihoods.marginalize(identicalHaplotypesMap);
            logger.debug(String.format("%d haplotypes after uncollapsing", uncollapsedReadLikelihoods.numberOfAlleles()));
        } else {
            logger.debug(String.format("Not performing uncollapsing with %d haplotypes", readLikelihoods.numberOfAlleles()));
            uncollapsedReadLikelihoods = readLikelihoods;
        }


        final AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal;
        // Optional: pre-filter haplotypes by removing weak/noisy alleles
        // this is important for the genotyping, as weak allele close to the real allele often decreases its quality

        Set<Integer> suspiciousLocations = new HashSet<>();
        if (MTAC.filterAlleles) {
            logger.debug("Filtering alleles");
            AlleleFilteringMutect alleleFilter = new AlleleFilteringMutect(MTAC, null, genotypingEngine);
            EventMap.buildEventMapsForHaplotypes(uncollapsedReadLikelihoods.alleles(),
                    assemblyResult.getFullReferenceWithPadding(),
                    assemblyResult.getPaddedReferenceLoc(),
                    MTAC.assemblerArgs.debugAssembly,
                    MTAC.maxMnpDistance);

            subsettedReadLikelihoodsFinal = alleleFilter.filterAlleles(uncollapsedReadLikelihoods,
                    assemblyResult.getPaddedReferenceLoc().getStart(), suspiciousLocations);
        } else {
            logger.debug("Not filtering alleles");
            subsettedReadLikelihoodsFinal = uncollapsedReadLikelihoods;
        }



        final CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                subsettedReadLikelihoodsFinal, assemblyResult, referenceContext,
                regionForGenotyping.getSpan(), featureContext, givenAlleles,
                header, haplotypeBAMWriter.isPresent(), emitReferenceConfidence(),
                suspiciousLocations);
        writeBamOutput(uncollapsedReadLikelihoods.alleles(),
                assemblyResult.getPaddedReferenceLoc(),
                subsettedReadLikelihoodsFinal,
                calledHaplotypes,
                regionForGenotyping.getSpan());

        if (emitReferenceConfidence()) {
            if ( !containsCalls(calledHaplotypes) ) {
                // no called all of the potential haplotypes
                return referenceModelForNoVariation(originalAssemblyRegion);
            } else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section, then variant-containing section, then right flank
                trimmingResult.nonVariantLeftFlankRegion().ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank)));

                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        subsettedReadLikelihoodsFinal, new HomogeneousPloidyModel(samplesList, 2), calledHaplotypes.getCalls()));

                trimmingResult.nonVariantRightFlankRegion().ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank)));

                return result;
            }
        }
        else {
            return calledHaplotypes.getCalls();
        }
    }

    private void removeReadStubs(final AssemblyRegion assemblyRegion) {
        final List<GATKRead> readStubs = assemblyRegion.getReads().stream()
                .filter(r -> r.getLength() < AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING).collect(Collectors.toList());
        assemblyRegion.removeAll(readStubs);
    }

    private void removeUnmarkedDuplicates(final AssemblyRegion assemblyRegion) {
        // PCR duplicates whose mates have MQ = 0 won't necessarily be marked as duplicates because
        // the mates map randomly to many different places
        final Map<ImmutablePair<String, Integer>, List<GATKRead>> possibleDuplicates = assemblyRegion.getReads().stream()
                .filter(read -> read.isPaired() && !read.mateIsUnmapped() &&
                        (!read.getMateContig().equals(read.getContig()) || Math.abs(read.getFragmentLength()) > HUGE_FRAGMENT_LENGTH))
                .collect(Collectors.groupingBy(
                        read -> ImmutablePair.of(ReadUtils.getSampleName(read, header), (read.isFirstOfPair() ? 1 : -1) * read.getUnclippedStart())));

        final List<GATKRead> duplicates = possibleDuplicates.values().stream().flatMap(list -> {
            final Map<String, List<GATKRead>> readsByContig = list.stream().collect(Collectors.groupingBy(GATKRead::getMateContig));
            // if no clear best contig, they're almost certainly all mapping errors (in addition to all but one being duplicates)
            // so we toss them all.  Otherwise we toss all but one
            return readsByContig.values().stream().flatMap(contigReads -> contigReads.stream().skip(contigReads.size() > list.size() / 2 ? 1 : 0));
        }).collect(Collectors.toList());
        assemblyRegion.removeAll(duplicates);
    }

    //TODO: refactor this
    private boolean containsCalls(final CalledHaplotypes calledHaplotypes) {
        return calledHaplotypes.getCalls().stream()
                .flatMap(call -> call.getGenotypes().stream())
                .anyMatch(Genotype::isCalled);
    }

    private void writeBamOutput(final Collection<Haplotype> initialHaplotypeList, Locatable paddedReferenceLoc, final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods, final CalledHaplotypes calledHaplotypes, final Locatable callableRegion) {
        if ( haplotypeBAMWriter.isPresent() ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            haplotypeBAMWriter.get().writeReadsAlignedToHaplotypes(
                    initialHaplotypeList,
                    paddedReferenceLoc,
                    initialHaplotypeList,
                    calledHaplotypeSet,
                    readLikelihoods,
                    callableRegion);
        }
    }

    //TODO: should be a variable, not a function
    private boolean hasNormal() {
        return !normalSamples.isEmpty();
    }

    private boolean isNormalSample(final String sample) {
        return normalSamples.contains(sample);
    }

    private boolean isTumorSample(final String sample) {
        return !isNormalSample(sample);
    }

    protected Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return AssemblyBasedCallerUtils.splitReadsBySample(samplesList, header, reads);
    }

    public void writeExtraOutputs(final File statsTable) {
        final List<MutectStats> stats = Arrays.asList(new MutectStats(CALLABLE_SITES_NAME, callableSites.getValue()));
        MutectStats.writeToFile(stats, statsTable);
        f1R2CountsCollector.ifPresent(collector -> {
            collector.writeHistograms();
            collector.closeAndArchiveFiles();
        });
    }

    @Override
    public void close() {
        likelihoodCalculationEngine.close();
        aligner.close();
        haplotypeBAMWriter.ifPresent(HaplotypeBAMWriter::close);
        assembledEventMapVcfOutputWriter.ifPresent(writer -> {assembledEventMapVariants.get().forEach(writer::add); writer.close();});
        referenceReader.close();
        genotypingEngine.close();
    }

    @Override
    public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext features) {
        if ( forceCallingAllelesPresent && features.getValues(MTAC.alleles, ref).stream().anyMatch(vc -> MTAC.forceCallFiltered || vc.isNotFiltered())) {
            return new ActivityProfileState(ref.getInterval(), 1.0);
        }

        final byte refBase = ref.getBase();
        final SimpleInterval refInterval = ref.getInterval();

        if( context == null || context.getBasePileup().isEmpty() ) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        final ReadPileup pileup = context.getBasePileup();
        if (MTAC.pileupDetectionArgs.usePileupDetection) {
            pileup.forEach(p -> PileupBasedAlleles.addMismatchPercentageToRead(p.getRead(), header, ref));
        }
        if (pileup.size() >= minCallableDepth) {
            callableSites.increment();
        }
        final ReadPileup tumorPileup = pileup.makeFilteredPileup(pe -> isTumorSample(ReadUtils.getSampleName(pe.getRead(), header)));
        f1R2CountsCollector.ifPresent(collector -> collector.process(tumorPileup, ref));
        tumorPileupQualBuffer.accumulateQuals(tumorPileup, refBase, MTAC.pcrSnvQual);
        final Pair<Integer, ByteArrayList> bestTumorAltAllele = tumorPileupQualBuffer.likeliestIndexAndQuals();
        final double tumorLogOdds = logLikelihoodRatio(tumorPileup.size() - bestTumorAltAllele.getRight().size(), bestTumorAltAllele.getRight());

        if (tumorLogOdds < MTAC.getInitialLogOdds()) {
            return new ActivityProfileState(refInterval, 0.0);
        } else if (MTAC.mutect3TrainingDataMode) {
            return new ActivityProfileState(ref.getInterval(), 1.0);
        } else if (hasNormal() && !MTAC.genotypeGermlineSites) {
            final ReadPileup normalPileup = pileup.makeFilteredPileup(pe -> isNormalSample(ReadUtils.getSampleName(pe.getRead(), header)));
            normalPileupQualBuffer.accumulateQuals(normalPileup, refBase, MTAC.pcrSnvQual);
            final Pair<Integer, ByteArrayList> bestNormalAltAllele = normalPileupQualBuffer.likeliestIndexAndQuals();
            if (bestNormalAltAllele.getLeft() == bestTumorAltAllele.getLeft()) {
                final int normalAltCount = bestNormalAltAllele.getRight().size();
                final double normalQualSum = normalPileupQualBuffer.qualSum(bestNormalAltAllele.getLeft());
                if (normalAltCount > normalPileup.size() * MAX_ALT_FRACTION_IN_NORMAL && normalQualSum > MAX_NORMAL_QUAL_SUM) {
                    return new ActivityProfileState(refInterval, 0.0);
                }
            }
        } else if (!MTAC.genotypeGermlineSites) {
            final List<VariantContext> germline = features.getValues(MTAC.germlineResource, refInterval);

            for (final VariantContext germlineVC : germline) {
                final List<Double> germlineAlleleFrequencies = getAttributeAsDoubleList(germlineVC, VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
                final Allele germlineRef = germlineVC.getReference();

                for (int germlineAltIdx = 0; germlineAltIdx < germlineAlleleFrequencies.size(); germlineAltIdx++) {
                    if (germlineAlleleFrequencies.get(germlineAltIdx) < MTAC.maxPopulationAlleleFrequency) {
                        continue;
                    }

                    final Allele germlineAlt = germlineVC.getAlternateAllele(germlineAltIdx);

                    // if it's a substitution that shares its first base with the dominant tumor allele, or if it's an
                    // indel and the dominant tumor allele is an indel, skip
                    if (PileupQualBuffer.likeliestIndexIsIndel(bestTumorAltAllele.getLeft()) && germlineAlt.length() != germlineRef.length()) {
                            return new ActivityProfileState(refInterval, 0.0);
                    } else if (PileupQualBuffer.likeliestIndexIsSubstitution(bestTumorAltAllele.getLeft()) && germlineAlt.length() == germlineRef.length()
                            && PileupQualBuffer.getSubstitutionBase(bestTumorAltAllele.getLeft()) == germlineRef.getBases()[0]) {
                        return new ActivityProfileState(refInterval, 0.0);
                    }

                }
            }
        }

        if (!MTAC.genotypePonSites && !features.getValues(MTAC.pon, new SimpleInterval(context.getContig(), (int) context.getPosition(), (int) context.getPosition())).isEmpty()) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        // if a site is active, count it toward the total of callable sites even if its depth is below the threshold
        if (pileup.size() < minCallableDepth) {
            callableSites.increment();
        }
        return new ActivityProfileState( refInterval, 1.0, ActivityProfileState.Type.NONE, null);
    }

    // NOTE: this is a hack to get around an htsjdk bug: https://github.com/samtools/htsjdk/issues/1228
    // htsjdk doesn't correctly detect the missing value string '.', so we have copied and fixed the htsjdk code
    public static List<Double> getAttributeAsDoubleList(final VariantContext vc, final String key, final double defaultValue) {
        return vc.getCommonInfo().getAttributeAsList(key).stream()
                .map(x -> {
                    if (x == null) {
                        return defaultValue;
                    } else if (x instanceof Number) {
                        return ((Number) x).doubleValue();
                    } else {
                        String string = (String) x;
                        return string.equals(VCFConstants.MISSING_VALUE_v4) ? defaultValue : Double.valueOf(string); // throws an exception if this isn't a string
                    }
                }).collect(Collectors.toList());
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return MTAC.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    private List<VariantContext> referenceModelForNoVariation(final AssemblyRegion region) {
        // don't correct overlapping base qualities because we did that upstream
        AssemblyBasedCallerUtils.finalizeRegion(region,
                false,
                true,
                (byte) MIN_TAIL_QUALITY,
                header,
                samplesList,
                false,
                false,
                false,
                MTAC.pileupDetectionArgs.usePileupDetection);  //take off soft clips and low Q tails before we calculate likelihoods


        final SimpleInterval paddedLoc = region.getPaddedSpan();
        final Haplotype refHaplotype = AssemblyBasedCallerUtils.createReferenceHaplotype(region, paddedLoc, referenceReader);
        final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
        return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                paddedLoc, region, AssemblyBasedCallerUtils.createDummyStratifiedReadMap(refHaplotype, samplesList, header, region),
                new HomogeneousPloidyModel(samplesList, 2), Collections.emptyList(), false, Collections.emptyList()); //TODO: clean up args
    }

    private static int getCurrentOrFollowingIndelLength(final PileupElement pe) {
        return pe.isDeletion() ? pe.getCurrentCigarElement().getLength() : pe.getLengthOfImmediatelyFollowingIndel();
    }

    private static byte indelQual(final int indelLength) {
        return (byte) Math.min(INDEL_START_QUAL + (indelLength - 1) * INDEL_CONTINUATION_QUAL, Byte.MAX_VALUE);
    }

    public static double logLikelihoodRatio(final int refCount, final List<Byte> altQuals) {
        return logLikelihoodRatio(refCount, altQuals, 1);
    }

    /**
     *  this implements the isActive() algorithm described in docs/mutect/mutect.pdf
     *  the multiplicative factor is for the special case where we pass a singleton list
     *  of alt quals and want to duplicate that alt qual over multiple reads
     * @param nRef          ref read count
     * @param altQuals      Phred-scaled qualities of alt-supporting reads
     * @param repeatFactor  Number of times each alt qual is duplicated
     * @param afPrior       Beta prior on alt allele fraction
     * @return
     */
    public static double logLikelihoodRatio(final int nRef, final List<Byte> altQuals, final int repeatFactor, final Optional<BetaDistributionShape> afPrior) {
        final int nAlt = repeatFactor * altQuals.size();
        final int n = nRef + nAlt;

        final double fTildeRatio = FastMath.exp(MathUtils.digamma(nRef + 1) - MathUtils.digamma(nAlt + 1));


        double readSum = 0;
        for (final byte qual : altQuals) {
            final double epsilon = QualityUtils.qualToErrorProb(qual);
            final double zBarAlt = (1 - epsilon) / (1 - epsilon + epsilon * fTildeRatio);
            final double logEpsilon = NaturalLogUtils.qualToLogErrorProb(qual);
            final double logOneMinusEpsilon = NaturalLogUtils.qualToLogProb(qual);
            readSum += zBarAlt * (logOneMinusEpsilon - logEpsilon) + MathUtils.fastBernoulliEntropy(zBarAlt);
        }

        final double betaEntropy;
        if (afPrior.isPresent()) {
            final double alpha = afPrior.get().getAlpha();
            final double beta = afPrior.get().getBeta();
            betaEntropy = Gamma.logGamma(alpha + beta) - Gamma.logGamma(alpha) - Gamma.logGamma(beta)
                    - Gamma.logGamma(alpha + beta + n) + Gamma.logGamma(alpha + nAlt) + Gamma.logGamma(beta + nRef);
        } else {
            betaEntropy = -Math.log(n + 1) - CombinatoricsUtils.binomialCoefficientLog(n, nAlt);
        }
        return betaEntropy + readSum * repeatFactor;
    }

    // the default case of a flat Beta(1,1) prior on allele fraction
    public static double logLikelihoodRatio(final int nRef, final List<Byte> altQuals, final int repeatFactor) {
        return logLikelihoodRatio(nRef, altQuals, repeatFactor, Optional.empty());
    }

    // same as above but with a constant error probability for several alts
    public static double logLikelihoodRatio(final int refCount, final int altCount, final double errorProbability) {
        final byte qual = QualityUtils.errorProbToQual(errorProbability);
        return logLikelihoodRatio(refCount, Collections.singletonList(qual), altCount);
    }

    // check that we're next to a soft clip that is not due to a read that got out of sync and ended in a bunch of BQ2's
    // we only need to check the next base's quality
    private static boolean isNextToUsefulSoftClip(final PileupElement pe) {
        final int offset = pe.getOffset();
        return pe.getQual() > MINIMUM_BASE_QUALITY &&
                ((pe.isBeforeSoftClip() && pe.getRead().getBaseQuality(offset + 1) > MINIMUM_BASE_QUALITY)
                        || (pe.isAfterSoftClip() && pe.getRead().getBaseQuality(offset - 1) > MINIMUM_BASE_QUALITY));
    }

    private void checkSampleInBamHeader(final String sample) {
        if (sample != null && !samplesList.asListOfSamples().contains(sample)) {
            throw new UserException.BadInput("Sample " + sample + " is not in BAM header: " + samplesList.asListOfSamples());
        }
    }

    private String decodeSampleNameIfNecessary(final String name) {
        return samplesList.asListOfSamples().contains(name) ? name : IOUtils.urlDecode(name);
    }

    /**
     * A resuable container class to accumulate qualities for each type of SNV and indels (all indels combined)
     */
    private static class PileupQualBuffer {
        private static final int OTHER_SUBSTITUTION = 4;
        private static final int INDEL = 5;

        private static double MULTIPLE_SUBSTITUTION_BASE_QUAL_CORRECTION;

        // indices 0-3 are A,C,G,T; 4 is other substitution (just in case it's some exotic protocol); 5 is indel
        private List<ByteArrayList> buffers = IntStream.range(0,6).mapToObj(n -> new ByteArrayList()).collect(Collectors.toList());

        public PileupQualBuffer(final double multipleSubstitutionQualCorrection) {
            MULTIPLE_SUBSTITUTION_BASE_QUAL_CORRECTION = multipleSubstitutionQualCorrection;
        }

        public void accumulateQuals(final ReadPileup pileup, final byte refBase, final int pcrErrorQual) {
            clear();
            final int position = pileup.getLocation().getStart();

            for (final PileupElement pe : pileup) {
                final int indelLength = getCurrentOrFollowingIndelLength(pe);
                if (indelLength > 0) {
                    accumulateIndel(indelQual(indelLength));
                } else if (isNextToUsefulSoftClip(pe)) {
                    accumulateIndel(indelQual(1));
                } else if (pe.getBase() != refBase && pe.getQual() > MINIMUM_BASE_QUALITY) {
                    final GATKRead read = pe.getRead();
                    final int mateStart = (!read.isProperlyPaired() || read.mateIsUnmapped()) ? Integer.MAX_VALUE : read.getMateStart();
                    final boolean overlapsMate = mateStart <= position && position < mateStart + read.getLength();
                    accumulateSubstitution(pe.getBase(), overlapsMate ? (byte) FastMath.min(pe.getQual(), pcrErrorQual/2) : pe.getQual());
                }
            }
        }

        public Pair<Integer, ByteArrayList> likeliestIndexAndQuals() {
            int bestIndex = 0;
            long bestSum = 0;
            for (int n = 0; n < buffers.size(); n++) {
                final long sum = qualSum(n);
                if (sum > bestSum) {
                    bestSum = sum;
                    bestIndex = n;
                }
            }
            return ImmutablePair.of(bestIndex, buffers.get(bestIndex));
        }

        public static boolean likeliestIndexIsIndel(final int index) {
            return index == INDEL;
        }

        public static boolean likeliestIndexIsSubstitution(final int index) {
            return index < OTHER_SUBSTITUTION;
        }

        public static byte getSubstitutionBase(final int index) {
            return BaseUtils.baseIndexToSimpleBase(index);
        }

        private void accumulateSubstitution(final byte base, final byte qual) {
            final int index = BaseUtils.simpleBaseToBaseIndex(base);
            if (index == -1) {  // -1 is the hard-coded value for non-simple bases in BaseUtils
                buffers.get(OTHER_SUBSTITUTION).add(qual);
            } else {
                buffers.get(index).add((byte) FastMath.min(qual + MULTIPLE_SUBSTITUTION_BASE_QUAL_CORRECTION, QualityUtils.MAX_QUAL));
            }
        }

        private void accumulateIndel(final byte qual) {
            buffers.get(INDEL).add(qual);
        }

        private void clear() {
            buffers.forEach(ByteArrayList::clear);
        }

        public long qualSum(final int index) {
            final ByteArrayList list = buffers.get(index);
            long result = 0;
            for (int n = 0; n < list.size(); n++) {
                result += list.getByte(n);
            }
            return result;
        }
    }
}
