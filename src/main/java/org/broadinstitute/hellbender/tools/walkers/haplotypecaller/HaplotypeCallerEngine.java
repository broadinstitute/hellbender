package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.NamedFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.OutputMode;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.transformers.IUPACReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.stream.Collectors;

/**
 * The core engine for the HaplotypeCaller that does all of the actual work of the tool.
 * Usage:
 * -Pass the HaplotypeCaller args into the constructor, which will initialize the HC engine completely.
 * -Get the appropriate VCF or GVCF writer (depending on our arguments) from {@link #makeVCFWriter}
 * -Write the appropriate VCF header via {@link #writeHeader}
 * -Repeatedly call {@link #isActive} to identify active vs. inactive regions
 * -Repeatedly call {@link #callRegion} to call variants in each region, and add them to your writer
 * -When done, call {@link #shutdown}. Close the writer you got from {@link #makeVCFWriter} yourself.
 */
public class HaplotypeCallerEngine implements AssemblyRegionEvaluator {

    private static final Logger logger = LogManager.getLogger(HaplotypeCallerEngine.class);

    protected final HaplotypeCallerArgumentCollection hcArgs;

    protected final SAMFileHeader readsHeader;

    protected ReferenceConfidenceModel referenceConfidenceModel = null;

    protected AssemblyRegionTrimmer trimmer;

    protected final OutputStreamWriter assemblyDebugOutStream;

    /**
     * List of interval file entries including regions and custom ploidy values to apply in that region.
     */
    protected final List<SimpleCount> ploidyRegions = new ArrayList<>();

    /**
     * An OverlapDetector object for checking whether region overlaps given ploidyRegions.
     */
    protected final OverlapDetector<SimpleCount> ploidyRegionsOverlapDetector;

    /**
     * List of all custom ploidies provided by user
     */
    private final LinkedHashSet<Integer> allCustomPloidies;

    /**
     * The default genotyping engine for the isActive() determination
     */
    private MinimalGenotypingEngine defaultActiveRegionEvaluationGenotyperEngine = null;

    /**
     * Map of user-provided ploidy values to corresponding active region genotyper. Values are checked as valid Integers during
     * initialization, but Strings are used as keys to avoid parsing repeatedly during runtime.
     */
    private final HashMap<Integer, MinimalGenotypingEngine> ploidyToActiveEvaluationGenotyper = new HashMap<>();

    protected ReadThreadingAssembler assemblyEngine = null;

    protected ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;
    private ReadLikelihoodCalculationEngine filterStepLikelihoodCalculationEngine = null;
    // If we are in PDHMM mode we need to hold onto two likelihoods engines for the fallback
    private ReadLikelihoodCalculationEngine pdhmmLikelihoodCalculationEngine = null;

    /**
     * The default genotyping engine to use for actual variant calling and genotyping in an active region.
     */
    protected HaplotypeCallerGenotypingEngine defaultGenotypingEngine = null;

    /**
     * Map of user-provided ploidy values to corresponding genotyper. Values are checked as valid Integers during
     * initialization, but Strings are used as keys to avoid parsing repeatedly during runtime.
     */
    protected final HashMap<Integer, HaplotypeCallerGenotypingEngine> ploidyToGenotyperMap = new HashMap<>();

    private VariantAnnotatorEngine annotationEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    protected final ReferenceSequenceFile referenceReader;

    // writes Haplotypes to a bam file when the -bamout option is specified
    protected Optional<HaplotypeBAMWriter> haplotypeBAMWriter;

    // writes Variants from assembly graph
    private Optional<VariantContextWriter> assembledEventMapVcfOutputWriter;
    private Optional<PriorityQueue<VariantContext>> assembledEventMapVariants;

    protected Optional<AlleleLikelihoodWriter> alleleLikelihoodWriter;

    private Set<String> sampleSet;
    protected SampleList samplesList;

    private final boolean forceCallingAllelesPresent;

    private byte minTailQuality;

    protected SmithWatermanAligner aligner;

    private final DragstrParams dragstrParams;

    public static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    /**
     * Minimum (exclusive) average number of high quality bases per soft-clip to consider that a set of soft-clips is a
     * high quality set.
     */
    private static final double AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD = 6.0;

    /**
     * Maximum-mininum confidence on a variant to exist to consider the position as a potential variant harbouring locus
     * when looking for active regions.
     */
    private static final double MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY = 4.0;

    /**
     * Minimum ploidy assumed when looking for loci that may harbour variation to identify active regions.
     * <p>
     * By default we take the putative ploidy provided by the user, but this turned out to be too insensitive
     * for low ploidy, notoriously with haploid samples. Therefore we impose this minimum.
     * </p>
     */
    private static final int MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY = 2;

    /**
     * Reads with length lower than this number, after clipping off overhangs outside the active region,
     * won't be considered for genotyping.
     */
    private static final int READ_LENGTH_FILTER_THRESHOLD = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    public static final int DEFAULT_READ_QUALITY_FILTER_THRESHOLD = 20;

    protected static final List<VariantContext> NO_CALLS = Collections.emptyList();

    private static final Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private static final Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     * @param hcArgs command-line arguments for the HaplotypeCaller
     * @param assemblyRegionArgs
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param readsHeader header for the reads
     * @param referenceReader reader to provide reference data, this reference reader will be closed when {@link #shutdown()} is called
     * @param annotationEngine variantAnnotatorEngine with annotations to process already added
     */
    public HaplotypeCallerEngine(final HaplotypeCallerArgumentCollection hcArgs, AssemblyRegionArgumentCollection assemblyRegionArgs, boolean createBamOutIndex,
                                 boolean createBamOutMD5, final SAMFileHeader readsHeader,
                                 CachingIndexedFastaSequenceFile referenceReader, VariantAnnotatorEngine annotationEngine) {
        this.dragstrParams = DragstrParamUtils.parse(hcArgs.likelihoodArgs.dragstrParams);
        this.hcArgs = Utils.nonNull(hcArgs);
        this.readsHeader = Utils.nonNull(readsHeader);
        this.referenceReader = Utils.nonNull(referenceReader);
        this.annotationEngine = Utils.nonNull(annotationEngine);
        this.aligner = SmithWatermanAligner.getAligner(hcArgs.smithWatermanImplementation);
        forceCallingAllelesPresent = hcArgs.alleles != null;

        // Add necessary debug streams to the output
        if (hcArgs.assemblyStateOutput != null) {
            try {
                assemblyDebugOutStream = new OutputStreamWriter(hcArgs.assemblyStateOutput.getOutputStream());
            } catch (final Exception e) {
                throw new UserException.CouldNotCreateOutputFile(hcArgs.assemblyStateOutput, "Provided argument for assembly debug graph location could not be created", e);
            }
        } else {
            assemblyDebugOutStream = null;
        }
        if (hcArgs.genotyperDebugOutStream != null) {
            HaplotypeCallerGenotypingDebugger.initialize(hcArgs.genotyperDebugOutStream);
        }

        // Parse the user provided custom ploidy regions into ploidyRegions object containing SimpleCounts
        if (this.hcArgs.ploidyRegions != null) {
            FeatureDataSource<NamedFeature> ploidyDataSource = new FeatureDataSource<>(this.hcArgs.ploidyRegions, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, NamedFeature.class);
            ploidyDataSource.forEach(r -> this.ploidyRegions.add(new SimpleCount(r)));
        }

        for (SimpleCount region : this.ploidyRegions) {
            if (!IntervalUtils.intervalIsOnDictionaryContig(region.getInterval(), readsHeader.getSequenceDictionary())) {
                throw new UserException("Invalid region provided for --ploidy-regions at " + region.getContig() + ":" + region.getStart() + "-" + region.getEnd() + ". Contig name or endpoint doesn't match read sequence dictionary.");
            }
        }

        this.ploidyRegionsOverlapDetector = OverlapDetector.create(this.ploidyRegions);
        this.allCustomPloidies = this.ploidyRegions.stream().map(SimpleCount::getCount).collect(Collectors.toCollection(LinkedHashSet::new));

        trimmer = new AssemblyRegionTrimmer(assemblyRegionArgs, readsHeader.getSequenceDictionary());
        initialize(createBamOutIndex, createBamOutMD5);
    }

    /**
     * Common method to use in order to remove unwanted annotations from the list returned by the plugin specifically
     * for reference confidence mode. Will also ensure StrandBiasBySample is present regardless of user requests.
     *
     * @param annotations a list of annotations to change
     * @return a list of annotations with non GVCF annotations removed
     */
    public static Collection<Annotation> filterReferenceConfidenceAnnotations(Collection<Annotation> annotations) {
        logger.info("Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled");

        // Override user preferences and add StrandBiasBySample
        if (!annotations.contains(new StrandBiasBySample())) {
            annotations.add(new StrandBiasBySample());
        }
        // Override user preferences and remove ChromosomeCounts, FisherStrand, StrandOddsRatio, and QualByDepth Annotations
        return annotations.stream()
                .filter(c -> (
                        c.getClass() != (ChromosomeCounts.class) &&
                        c.getClass() != (FisherStrand.class) &&
                        c.getClass() != (StrandOddsRatio.class) &&
                        c.getClass() != (QualByDepth.class))
                ).collect(Collectors.toList());
    }

    /**
     * @return the default set of variant annotations for use with HaplotypeCaller
     */
    public static List<Class<? extends Annotation>> getStandardHaplotypeCallerAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, StandardHCAnnotation.class);
    }

    private void initialize(boolean createBamOutIndex, final boolean createBamOutMD5) {
        // Note: order of operations matters here!

        initializeSamples();

        // Must be called after initializeSamples()
        validateAndInitializeArgs();
        minTailQuality = (byte)(hcArgs.minBaseQualityScore - 1);

        initializeActiveRegionEvaluationGenotyperEngine();

        defaultGenotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, ! hcArgs.doNotRunPhysicalPhasing, hcArgs.applyBQD);
        defaultGenotypingEngine.setAnnotationEngine(annotationEngine);

        // Create other custom genotyping engines if user provided ploidyRegions
        for (final int ploidy : this.allCustomPloidies) {
            HaplotypeCallerArgumentCollection newPloidyHcArgs = hcArgs.copyWithNewPloidy(ploidy);
            HaplotypeCallerGenotypingEngine newGenotypingEngine = new HaplotypeCallerGenotypingEngine(newPloidyHcArgs, samplesList, ! hcArgs.doNotRunPhysicalPhasing, hcArgs.applyBQD);
            newGenotypingEngine.setAnnotationEngine(annotationEngine);
            this.ploidyToGenotyperMap.put(ploidy, newGenotypingEngine);
        }

        boolean isFlowBased = (hcArgs.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBased)
                || (hcArgs.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBasedHMM);
        referenceConfidenceModel = new ReferenceConfidenceModel(samplesList, readsHeader,
                hcArgs.indelSizeToEliminateInRefModel,
                hcArgs.standardArgs.genotypeArgs.numRefIfMissing,
                hcArgs.refModelDelQual,
                !hcArgs.overrideSoftclipFragmentCheck,
                isFlowBased);

        //Allele-specific annotations are not yet supported in the VCF mode
        if (isAlleleSpecificExceptHmerLengthOrStrandBiasMode(annotationEngine) && isVCFMode()){
           throw new UserException("Allele-specific annotations are not yet supported in the VCF mode");
        }

        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(hcArgs, createBamOutIndex, createBamOutMD5, readsHeader);
        alleleLikelihoodWriter = AssemblyBasedCallerUtils.createAlleleLikelihoodWriter(hcArgs);
        assemblyEngine = hcArgs.createReadThreadingAssembler();
        assembledEventMapVcfOutputWriter = Optional.ofNullable(hcArgs.assemblerArgs.debugAssemblyVariantsOut != null ?
                GATKVariantContextUtils.createVCFWriter(
                        new GATKPath(hcArgs.assemblerArgs.debugAssemblyVariantsOut).toPath(),
                        readsHeader.getSequenceDictionary(),
                        false,
                        Options.DO_NOT_WRITE_GENOTYPES, Options.INDEX_ON_THE_FLY)
                : null);
        assembledEventMapVariants = Optional.ofNullable(hcArgs.assemblerArgs.debugAssemblyVariantsOut != null ?
                new PriorityQueue<>(200, new VariantContextComparator(readsHeader.getSequenceDictionary())) : null);
        assembledEventMapVcfOutputWriter.ifPresent(writer -> writeHeader(writer, readsHeader.getSequenceDictionary(), new HashSet<>()));
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs, !hcArgs.softClipLowQualityEnds, hcArgs.likelihoodArgs.likelihoodEngineImplementation);


        //Some sanity checking of the stepwise filtering approach
        if (hcArgs.stepwiseFiltering) {
            if (hcArgs.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBased) {
                throw new UserException("'--"+ HaplotypeCallerArgumentCollection.STEPWISE_FITLERING_ARGUMENT + "' is specified but the selected likelihoods engine is FlowBased. This is invalid.");
            }
            if (!hcArgs.filterAlleles) {
                throw new UserException("'--"+ HaplotypeCallerArgumentCollection.STEPWISE_FITLERING_ARGUMENT + "' is specified, but allele filtering is disabled. Please enable allele filtering with '--"+AssemblyBasedCallerArgumentCollection.FILTER_ALLELES+"'");
            }
        }

        filterStepLikelihoodCalculationEngine = (hcArgs.stepwiseFiltering?
                AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs, !hcArgs.softClipLowQualityEnds, ReadLikelihoodCalculationEngine.Implementation.FlowBased)
                : null);
        pdhmmLikelihoodCalculationEngine = (hcArgs.pileupDetectionArgs.generatePDHaplotypes?
                new PDPairHMMLikelihoodCalculationEngine((byte) hcArgs.likelihoodArgs.gcpHMM, hcArgs.likelihoodArgs.dontUseDragstrPairHMMScores ? null : DragstrParamUtils.parse(hcArgs.likelihoodArgs.dragstrParams),
                        hcArgs.likelihoodArgs.pairHMMNativeArgs.getPairHMMArgs(), hcArgs.likelihoodArgs.pairHMM, hcArgs.pileupDetectionArgs.pdhmmDebugOutputResults, AssemblyBasedCallerUtils.getGlobalMismatchingRateFromArgs(hcArgs.likelihoodArgs), hcArgs.likelihoodArgs.pcrErrorModel,
                        hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD, hcArgs.likelihoodArgs.enableDynamicReadDisqualification, hcArgs.likelihoodArgs.readDisqualificationThresholdConstant,
                        hcArgs.likelihoodArgs.expectedErrorRatePerBase, !hcArgs.likelihoodArgs.disableSymmetricallyNormalizeAllelesToReference, hcArgs.likelihoodArgs.disableCapReadQualitiesToMapQ, !hcArgs.softClipLowQualityEnds,
                        hcArgs.pileupDetectionArgs.pdhmmOptimization? hcArgs.informativeReadOverlapMargin : -1) //Logic to control whether we apply the pdhmm optimization
                : null);
    }

    private boolean isVCFMode() {
        return hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.NONE;
    }

    private boolean isAlleleSpecificExceptHmerLengthOrStrandBiasMode(final VariantAnnotatorEngine annotationEngine) {
        //HACK. Note: we can't use subclass information from ReducibleAnnotation (which would be the obvious choice)
        // because RMSMappingQuality is both a reducible annotation and a standard annotation.

        return annotationEngine.getInfoAnnotations().stream()
                .filter(infoFieldAnnotation -> !infoFieldAnnotation.getClass().getSimpleName().equals(AS_StrandBiasMutectAnnotation.class.getSimpleName()))
                .anyMatch(infoFieldAnnotation -> infoFieldAnnotation.getClass().getSimpleName().startsWith("AS_")) ||
                annotationEngine.getGenotypeAnnotations().stream()
                .anyMatch(genotypeAnnotation -> genotypeAnnotation.getClass().getSimpleName().startsWith("AS_"));
    }

    private void validateAndInitializeArgs() {
        if ( hcArgs.standardArgs.genotypeArgs.samplePloidy != HomoSapiensConstants.DEFAULT_PLOIDY && ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Currently, physical phasing is only available for diploid samples.");
        }

        if ( hcArgs.dontGenotype && emitReferenceConfidence() ) {
            throw new UserException("You cannot request gVCF output and 'do not genotype' at the same time");
        }

        if ( emitReferenceConfidence() ) {
            hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling = -0.0;
            logger.info("Standard Emitting and Calling confidence set to " + hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling + " for reference-model confidence output");
            if ( ! hcArgs.standardArgs.annotateAllSitesWithPLs ) {
                logger.info("All sites annotated with PLs forced to true for reference-model confidence output");
            }
            hcArgs.standardArgs.annotateAllSitesWithPLs = true;
        } else if ( ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Disabling physical phasing, which is supported only for reference-model confidence output");
        }

        if ( hcArgs.standardArgs.CONTAMINATION_FRACTION_FILE != null ) {
            hcArgs.standardArgs.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(hcArgs.standardArgs.CONTAMINATION_FRACTION_FILE, hcArgs.standardArgs.CONTAMINATION_FRACTION, sampleSet, logger));
        }

        Utils.validateArg(hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD >= QualityUtils.MIN_USABLE_Q_SCORE, "BASE_QUALITY_SCORE_THRESHOLD must be greater than or equal to " + QualityUtils.MIN_USABLE_Q_SCORE + " (QualityUtils.MIN_USABLE_Q_SCORE)");

        if ( emitReferenceConfidence() && samplesList.numberOfSamples() != 1 ) {
            throw new CommandLineException.BadArgumentValue(AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, "Can only be used in single sample mode currently. Use the --sample-name argument to run on a single sample out of a multi-sample BAM file.");
        }

        if (hcArgs.floorBlocks && !emitReferenceConfidence()) {
            throw new UserException(HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS + " refers to GVCF blocks," +
                    " so reference confidence mode (" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME +
                    ") must be specified.");
        }
    }

    private void initializeSamples() {
        sampleSet = ReadUtils.getSamplesFromHeader(readsHeader);
        samplesList = new IndexedSampleList(sampleSet);

        if ( hcArgs.sampleNameToUse != null ) {
            if ( ! sampleSet.contains(hcArgs.sampleNameToUse) ) {
                throw new CommandLineException.BadArgumentValue("--sample_name", "Specified name does not exist in input bam files");
            }
            if (sampleSet.size() == 1) {
                //No reason to incur performance penalty associated with filtering if they specified the name of the only sample
                hcArgs.sampleNameToUse = null;
            } else {
                samplesList = new IndexedSampleList(hcArgs.sampleNameToUse);
                sampleSet.clear();
                sampleSet.add(hcArgs.sampleNameToUse);
            }
        }
    }

    private void initializeActiveRegionEvaluationGenotyperEngine() {
        final StandardCallerArgumentCollection activeRegionArgs = new StandardCallerArgumentCollection();
        activeRegionArgs.copyStandardCallerArgsFrom(hcArgs.standardArgs);

        activeRegionArgs.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
        activeRegionArgs.genotypeArgs.standardConfidenceForCalling = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling); // low values used for isActive determination only, default/user-specified values used for actual calling
        activeRegionArgs.CONTAMINATION_FRACTION = 0.0;
        activeRegionArgs.CONTAMINATION_FRACTION_FILE = null;
        // Seems that at least with some test data we can lose genuine haploid variation if we use ploidy == 1
        activeRegionArgs.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, hcArgs.standardArgs.genotypeArgs.samplePloidy);

        defaultActiveRegionEvaluationGenotyperEngine = new MinimalGenotypingEngine(activeRegionArgs, samplesList);
        defaultActiveRegionEvaluationGenotyperEngine.setLogger(logger);

        // If custom ploidyRegions provided, create corresponding map for active region determination genotyper
        for (final int ploidy : this.allCustomPloidies) {
            StandardCallerArgumentCollection newPloidyActiveArgs = new StandardCallerArgumentCollection();
            newPloidyActiveArgs.copyStandardCallerArgsFrom(activeRegionArgs);
            newPloidyActiveArgs.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, ploidy);
            MinimalGenotypingEngine newActiveGenotypingEngine = new MinimalGenotypingEngine(newPloidyActiveArgs, samplesList);
            newActiveGenotypingEngine.setLogger(logger);
            this.ploidyToActiveEvaluationGenotyper.put(ploidy, newActiveGenotypingEngine);
        }
    }

    /**
     * @return the default set of read filters for use with the HaplotypeCaller
     */
    public static List<ReadFilter> makeStandardHCReadFilters() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(new MappingQualityReadFilter(DEFAULT_READ_QUALITY_FILTER_THRESHOLD));
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());

        return filters;
    }

    /**
     * @return the default read transformer to apply to HaplotypeCaller input
     */
    public static ReadTransformer makeStandardHCReadTransformer() {
        return new IUPACReadTransformer(true);
    }

    /**
     * Create a VCF or GVCF writer as appropriate, given our arguments
     *
     * @param outputVCF location to which the vcf should be written
     * @param readsDictionary sequence dictionary for the reads
     * @return a VCF or GVCF writer as appropriate, ready to use
     */
    public VariantContextWriter makeVCFWriter( final GATKPath outputVCF, final SAMSequenceDictionary readsDictionary,
                                               final boolean createOutputVariantIndex, final boolean  createOutputVariantMD5,
                                               final boolean sitesOnlyMode ) {
        Utils.nonNull(outputVCF);
        Utils.nonNull(readsDictionary);

        final List<Options> options = new ArrayList<>(2);
        if (createOutputVariantIndex) {options.add(Options.INDEX_ON_THE_FLY);}
        if (sitesOnlyMode) {options.add(Options.DO_NOT_WRITE_GENOTYPES);}

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                outputVCF.toPath(),
                readsDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()])
        );

        if ( hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
            try {
                writer = new GVCFWriter(writer, hcArgs.GVCFGQBands, hcArgs.floorBlocks);
            } catch ( IllegalArgumentException e ) {
                throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }
        }

        return writer;
    }

    /**
     * Create a VCF header.
     *
     * @param sequenceDictionary sequence dictionary for the reads
     * @return a VCF header
     */
    public VCFHeader makeVCFHeader( final SAMSequenceDictionary sequenceDictionary, final Set<VCFHeaderLine>  defaultToolHeaderLines ) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.addAll(defaultToolHeaderLines);

        headerInfo.addAll(defaultGenotypingEngine.getAppropriateVCFInfoHeaders());
        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions(emitReferenceConfidence()));
        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        if (hcArgs.standardArgs.genotypeArgs.supportVariants != null) {
            headerInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
            headerInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
            headerInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        }

        if ( ! hcArgs.doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
            headerInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.PHASE_SET_KEY));
        }

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        if (( hcArgs.flowAssemblyCollapseHKerSize > 0 ) || (hcArgs.flowAssemblyCollapseHKerSize==AssemblyBasedCallerUtils.DETERMINE_COLLAPSE_THRESHOLD)){
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXT_COLLAPSED_KEY));
        }

        if (hcArgs.filterAlleles) {
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSSIBLE_FP_ADJACENT_TP_KEY));
        }

        if ( emitReferenceConfidence() ) {
            headerInfo.addAll(referenceConfidenceModel.getVCFHeaderLines());
        }

        if (hcArgs.likelihoodArgs.dragstrParams != null) {
            headerInfo.addAll(DragstrVariantContextAnnotations.vcfHeaderLines());
        }

        if (hcArgs.standardArgs.genotypeArgs.genotypeAssignmentMethod == GenotypeAssignmentMethod.USE_POSTERIOR_PROBABILITIES) {
            headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_POSTERIORS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Float, "genotype posterior in Phred Scale" ));
            headerInfo.add(new VCFFormatHeaderLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Float, "genotype priors in Phred Scale"));
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo, sampleSet);
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        return vcfHeader;
    }

    /**
     * Writes an appropriate VCF header, given our arguments, to the provided writer
     *
     * @param vcfWriter writer to which the header should be written
     */
    public void writeHeader( final VariantContextWriter vcfWriter, final SAMSequenceDictionary sequenceDictionary,
                             final Set<VCFHeaderLine>  defaultToolHeaderLines) {
        Utils.nonNull(vcfWriter);
        vcfWriter.writeHeader(makeVCFHeader(sequenceDictionary, defaultToolHeaderLines));
    }

    /**
     * Determines the appropriate ploidy to use at given site for different genotyping engines.
     * @param region Current region of interest
     * @return Ploidy value to use here given user inputs, or -1 if fall back to default
     */
    private int getPloidyToUseAtThisSite(Locatable region) {
        Set<SimpleCount> overlaps = this.ploidyRegionsOverlapDetector.getOverlaps(region);
        // Return first engine for interval overlapping this region
        if (!overlaps.isEmpty()) {
            Iterator<SimpleCount> intervals = overlaps.iterator();
            int max = intervals.next().getCount();
            while (intervals.hasNext()) {
                int next = intervals.next().getCount();
                if (next > max) {
                    max = next;
                }
            }
            return max;
        } else {
            return -1; // Sentinel value to fall back to default genotyper
        }
    }

    /**
     * Selects appropriate active region genotyping engine for given region
     * @param region Current region of interest
     * @return Active genotyping engine with correct ploidy setting for given region
     */
    private MinimalGenotypingEngine getLocalActiveGenotyper(Locatable region) {
        int currentPloidy = getPloidyToUseAtThisSite(region);
        if (currentPloidy == -1) {
            return this.defaultActiveRegionEvaluationGenotyperEngine;
        } else {
            return this.ploidyToActiveEvaluationGenotyper.get(currentPloidy);
        }
    }

    /**
     * Selects appropriate genotyping engine for given region.
     * @param region Current region of interest, e.g. AssemblyRegion
     * @return Genotyping engine with correct ploidy setting for given region
     */
    protected HaplotypeCallerGenotypingEngine getLocalGenotypingEngine(Locatable region) {
        int currentPloidy = getPloidyToUseAtThisSite(region);
        if (currentPloidy == -1) {
            return this.defaultGenotypingEngine;
        } else {
            return this.ploidyToGenotyperMap.get(currentPloidy);
        }
    }

    /**
     * Given a pileup, returns an ActivityProfileState containing the probability (0.0 to 1.0) that it's an "active" site.
     *
     * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
     * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
     * of probabilities. This is consistent with GATK 3.
     *
     * @param context reads pileup to examine
     * @param ref reference base overlapping the pileup locus
     * @param features features overlapping the pileup locus
     * @return probability between 0.0 and 1.0 that the site is active (in practice with this implementation: either 0.0 or 1.0)
     */
    @Override
    public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext features) {
        MinimalGenotypingEngine localActiveGenotypingEngine = getLocalActiveGenotyper(ref);

        if (forceCallingAllelesPresent && features.getValues(hcArgs.alleles, ref).stream().anyMatch(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())) {
            return new ActivityProfileState(ref.getInterval(), 1.0);
        }

        if (context == null || context.getBasePileup().isEmpty()) {
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getInterval(), 0.0);
        }

        final int ploidy = localActiveGenotypingEngine.getConfiguration().genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied

        final Map<String, AlignmentContext> splitContexts;
        if (samplesList.numberOfSamples() == 1) {
            //If we know a priori that there's just one sample, take a shortcut and dont examine each read in the pileup
            splitContexts = context.splitContextBySampleName(samplesList.getSample(0), readsHeader);
        } else {
            splitContexts = context.splitContextBySampleName(readsHeader);
        }

        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        for (final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet()) {
            // We summarize the reads mismatches vs the ref at this stage to avoid issues due to reads being hardclipped down the road.
            if (hcArgs.pileupDetectionArgs.usePileupDetection) {
                sample.getValue().getBasePileup().forEach(p -> PileupBasedAlleles.addMismatchPercentageToRead(p.getRead(), readsHeader, ref));
            }
            // The ploidy here is not dictated by the sample but by the simple genotyping-engine used to determine whether regions are active or not.
            final int activeRegionDetectionHackishSamplePloidy = localActiveGenotypingEngine.getConfiguration().genotypeArgs.samplePloidy;
            final double[] genotypeLikelihoods = ((RefVsAnyResult) referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(
                    activeRegionDetectionHackishSamplePloidy,
                    sample.getValue().getBasePileup(), ref.getBase(),
                    hcArgs.minBaseQualityScore,
                    averageHQSoftClips, false)).genotypeLikelihoods;
            genotypes.add(new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make());
        }

        final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE, FAKE_ALT_ALLELE);
        final double isActiveProb;

        if (genotypes.size() == 1) {
            // Faster implementation using exact marginalization instead of iteration
            isActiveProb = localActiveGenotypingEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(genotypes.get(0).getLikelihoods().getAsVector());
        } else {
            final VariantContext vcOut = localActiveGenotypingEngine.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getEnd(), alleles).genotypes(genotypes).make());
            isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb(vcOut.getPhredScaledQual());
        }

        // Here we store the crude fast-likelihood score used to determine positional activity for the ActivityProfileState.
        // We need to use this score later for DRAGEN-GATK, as one of the Heuristics applied by the new Pileup-Caller is that
        // now we optionally will ignore pileup events from positions on the genome that did not ping as being "active," thus
        // we must take care to preserve this information for later.
        //TODO If you dare, the entire ##PileupBasedAlleles.getPileupVariantContexts() step could be moved here to cut down on storing these intermediate scores for later
        //TODO This equivalent operation does not currently exist in the #Mutect2.isActive() method
        ActivityProfileState output = new ActivityProfileState(ref.getInterval(), isActiveProb, averageHQSoftClips.mean() > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean());
        double[] gts = genotypes.get(0).getLikelihoods().getAsVector();
        int maxElementIndex = MathUtils.maxElementIndex(gts);
        output.setOriginalActiveProb(gts[maxElementIndex] - gts[0]);
        return output;

    }

    /**
     * Generate variant calls for an assembly region
     *
     * @param region region to assemble and perform variant calling on
     * @param features Features overlapping the assembly region
     * @return List of variants discovered in the region (may be empty)
     */
    public List<VariantContext> callRegion(final AssemblyRegion region, final FeatureContext features, final ReferenceContext referenceContext) {
        final HaplotypeCallerGenotypingEngine localGenotypingEngine = getLocalGenotypingEngine(region);

        if ( hcArgs.justDetermineActiveRegions ) {
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;
        }
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("calling for region: " +region.getSpan());
        }

        final List<VariantContext> VCpriors = new ArrayList<>();
        if (hcArgs.standardArgs.genotypeArgs.supportVariants != null) {
            features.getValues(hcArgs.standardArgs.genotypeArgs.supportVariants).stream().forEach(VCpriors::add);
        }

        if ( hcArgs.sampleNameToUse != null ) {
            removeReadsFromAllSamplesExcept(hcArgs.sampleNameToUse, region);
        }

        if( ! region.isActive() ) {
            // Not active so nothing to do!
            return referenceModelForNoVariation(region, true, VCpriors);
        }

        final List<Event> givenAlleles = features.getValues(hcArgs.alleles).stream()
                .filter(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())
                .flatMap(vc -> GATKVariantContextUtils.splitVariantContextToEvents(vc, false, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false).stream())
                .filter(event -> event.getStart() >= region.getSpan().getStart()) // filter out events that do not start within the region. This approach works because events that begin upstream of the calling window cannot be called by this region calling code in the frist place. 
                .collect(Collectors.toList());

        if( givenAlleles.isEmpty() && region.size() == 0 ) {
            // No reads here so nothing to do!
            return referenceModelForNoVariation(region, true, VCpriors);
        }

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.write("\n\n\n\n" + region.getSpan() + "\nNumber of reads in region: " + region.getReads().size() + "     they are:\n");
                for (GATKRead read : region.getReads()) {
                    assemblyDebugOutStream.write(read.getName() + "   " + read.convertToSAMRecord(region.getHeader()).getFlags() + "\n");
                }
            } catch (IOException e) {
                throw new UserException("Error writing to debug output stream", e);
            }
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(region, hcArgs, readsHeader, samplesList, logger, referenceReader, assemblyEngine, aligner, !hcArgs.doNotCorrectOverlappingBaseQualities, hcArgs.fbargs, false);
        ReadThreadingAssembler.addAssembledVariantsToEventMapOutput(untrimmedAssemblyResult, assembledEventMapVariants, hcArgs.maxMnpDistance, assembledEventMapVcfOutputWriter);

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.write("\nThere were " + untrimmedAssemblyResult.getHaplotypeList().size() + " haplotypes found. Here they are:\n");
                for (final String haplotype : untrimmedAssemblyResult.getHaplotypeList().stream().map(Haplotype::toString).sorted().collect(Collectors.toList())) {
                    assemblyDebugOutStream.write(haplotype);
                    assemblyDebugOutStream.append('\n');
                }
            } catch (IOException e) {
                throw new UserException("Error writing to debug output stream", e);
            }
        }

        Pair<Set<Event>, Set<Event>> goodAndBadPileupEvents =
                PileupBasedAlleles.goodAndBadPileupEvents(region.getAlignmentData(), hcArgs.pileupDetectionArgs, readsHeader, hcArgs.minBaseQualityScore);
        final Set<Event> goodPileupEvents = goodAndBadPileupEvents.getLeft();
        final Set<Event> badPileupEvents = goodAndBadPileupEvents.getRight();

        // Regenerate the list of AllVariationEvents, filtering out assembled variants that must be filtered according to the pileupcaller code.
        final SortedSet<Event> allVariationEvents = untrimmedAssemblyResult.getVariationEvents(hcArgs.maxMnpDistance).stream()
                .filter(event -> !badPileupEvents.contains(event))
                .collect(Collectors.toCollection(() -> new TreeSet<>(AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR)));

        goodPileupEvents.forEach(allVariationEvents::add);
        givenAlleles.forEach(allVariationEvents::add);

        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(region, allVariationEvents, referenceContext);

        if (!trimmingResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            return referenceModelForNoVariation(region, false, VCpriors);
        }

        AssemblyResultSet assemblyResult = untrimmedAssemblyResult.trimTo(trimmingResult.getVariantRegion());
        assemblyResult.addGivenAlleles(givenAlleles, hcArgs.maxMnpDistance, aligner, hcArgs.getHaplotypeToReferenceSWParameters());

        // Pre-work for the PDHMM, if we are in PDHMM mode then here is where we re-compute the haplotypes as PD haplotypes.
        if (hcArgs.pileupDetectionArgs.generatePDHaplotypes) {
            // Note: we ignore maxMNPDistance here because the PartiallyDeterminedHaplotypeComputationEngine does not currently handle MNPs
            SortedSet<Event> assemblyVariants = assemblyResult.getVariationEvents(0);
            if (hcArgs.pileupDetectionArgs.debugPileupStdout) {
                System.out.println("Generating PDHaplotypes for PDHMM");
                System.out.println("Assembled Variants to use:");
                assemblyVariants.forEach(System.out::println);
                System.out.println("Good Pileup Variants to use:");
                goodPileupEvents.forEach(System.out::println);
                System.out.println("Bad Pileup Variants to filter:");
                badPileupEvents.forEach(System.out::println);
                System.out.println("Adding Variants To Reference Haplotype:");
                System.out.println(assemblyResult.getReferenceHaplotype());
                System.out.println("FinalSpan: " + assemblyResult.getReferenceHaplotype().getGenomeLocation());
                System.out.println("CallingSpan: " + region.getSpan());
            }
            assemblyResult = PartiallyDeterminedHaplotypeComputationEngine.generatePDHaplotypes(assemblyResult,
                    badPileupEvents,
                    goodPileupEvents,
                    aligner,
                    hcArgs);
        }

        // Legacy Pileupcaller code. Supplement the assembly haps with artifical haps constructed from the discovered pileupcaller
        // alleles based off of the GenotypeGivenAlleles code approach.
        // NOTE: We might also call this if hcArgs.pileupDetectionArgs.useGGAFallback is set and we are making PD haplotypes. This
        //       fallback handles edge cases where the PartiallyDeterminedHaplotypeComputationEngine generates too-many haplotypes
        //       from a very complex assebmly region and we want to fall back to assembly.
        if (!hcArgs.pileupDetectionArgs.generatePDHaplotypes ||
                        (hcArgs.pileupDetectionArgs.useGGAFallback && !assemblyResult.hasOverwrittenHaps())) { // If we are generating PDHaps assert that it failed before calling this
            assemblyResult.removeHaplotypesWithBadAlleles(hcArgs, badPileupEvents);
            assemblyResult.injectPileupEvents(region, hcArgs, aligner, goodPileupEvents);
        }
        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        final List<GATKRead> readStubs = regionForGenotyping.getReads().stream()
                .filter(r -> AlignmentUtils.unclippedReadLength(r)  < AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING).collect(Collectors.toList());
        regionForGenotyping.removeAll(readStubs);

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);
        final Map<String, List<GATKRead>> perSampleFilteredReadList = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, filteredReads);

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (!assemblyResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            return referenceModelForNoVariation(region, false, VCpriors);
        }

        // For sure this is not true if gVCF is on.
        if (hcArgs.dontGenotype) {
            return NO_CALLS; // user requested we not proceed
        }

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (regionForGenotyping.size() == 0 && !hcArgs.disableOptimizations) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(region, false, VCpriors);
        }

        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("\n=================================================");
            HaplotypeCallerGenotypingDebugger.println("assemblyRegion: "+new SimpleInterval(region));
            HaplotypeCallerGenotypingDebugger.println("=================================================");
        }

        // evaluate each sample's reads against all haplotypes
        final Map<String,List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, regionForGenotyping.getReads());
        final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods;
        List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();


        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("\nUnclipped Haplotypes("+haplotypes.size()+"):");
            for (Haplotype haplotype : untrimmedAssemblyResult.getHaplotypeList()) {
                HaplotypeCallerGenotypingDebugger.println("["+haplotype.getStart()+"-"+haplotype.getEnd()+"] k="+haplotype.getKmerSize()+" len: "+haplotype.length()+" "+haplotype.getCigar()+(haplotype.isReference()?"ref":""));
                HaplotypeCallerGenotypingDebugger.println(haplotype.toString());
            }

            HaplotypeCallerGenotypingDebugger.println("\nClipped Haplotyes("+haplotypes.size()+"):");
            for (Haplotype haplotype : haplotypes) {
                HaplotypeCallerGenotypingDebugger.println("["+haplotype.getStart()+"-"+haplotype.getEnd()+"] k="+haplotype.getKmerSize()+" len: "+haplotype.length()+" "+haplotype.getCigar()+(haplotype.isReference()?"ref":""));
                HaplotypeCallerGenotypingDebugger.println(haplotype.toString());
            }
            HaplotypeCallerGenotypingDebugger.println("");
        }

        // Calculate the likelihoods: CPU intensive part.
        // flow based alignment might add an extra step of uncollapsing - implemented by possiblyUncollapseHaplotypesInReadLikelihoods
        // non-flow based alignment will not be affected.
        readLikelihoods = possiblyUncollapseHaplotypesInReadLikelihoods(untrimmedAssemblyResult,
                hcArgs.stepwiseFiltering
                        ? filterStepLikelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads, true)
                        : ( assemblyResult.isPartiallyDeterminedList() ?
                        pdhmmLikelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads, true) :
                        likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads, true)));

        alleleLikelihoodWriter.ifPresent(
                writer -> writer.writeAlleleLikelihoods(readLikelihoods));

        // Optional: pre-filter haplotypes by removing weak/noisy alleles
        // this is important for the genotyping, as weak allele close to the real allele often decreases its quality
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal;
        Set<Integer> suspiciousLocations = new HashSet<>();
        if (hcArgs.filterAlleles) {
            logger.debug("Filtering alleles");

            AlleleFilteringHC alleleFilter = new AlleleFilteringHC(hcArgs, assemblyDebugOutStream, localGenotypingEngine);
            //need to update haplotypes to find the alleles
            EventMap.buildEventMapsForHaplotypes(readLikelihoods.alleles(),
                    assemblyResult.getFullReferenceWithPadding(),
                    assemblyResult.getPaddedReferenceLoc(),
                    hcArgs.assemblerArgs.debugAssembly,
                    hcArgs.maxMnpDistance);
            subsettedReadLikelihoodsFinal = alleleFilter.filterAlleles(readLikelihoods, assemblyResult.getPaddedReferenceLoc().getStart(), suspiciousLocations);

        } else {
            subsettedReadLikelihoodsFinal = readLikelihoods;
        }


        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]
        haplotypes = subsettedReadLikelihoodsFinal.alleles();

        // Reproscess with the HMM if we are in stepwise filtering mode
        if (hcArgs.stepwiseFiltering) {
            subsettedReadLikelihoodsFinal = (likelihoodCalculationEngine).computeReadLikelihoods(
                    subsettedReadLikelihoodsFinal.alleles(),
                    assemblyResult.getRegionForGenotyping().getHeader(),
                    samplesList, reads, true);
        }

        //Realign reads to their best haplotype.
        final SWParameters readToHaplotypeSWParameters = hcArgs.getReadToHaplotypeSWParameters();
        // TODO Yes we skip realignment entirely when we are in DRAGEN-GATK PDHMM mode. Realignment of the reads makes no sense when
        // TODO the bases of the haplotypes used for calling no longer reflect specified variants present.
        if (!(hcArgs.pileupDetectionArgs.generatePDHaplotypes && !hcArgs.pileupDetectionArgs.useDeterminedHaplotypesDespitePdhmmMode)) {
            final Map<GATKRead, GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(subsettedReadLikelihoodsFinal, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner, readToHaplotypeSWParameters);
            subsettedReadLikelihoodsFinal.changeEvidence(readRealignments);
        }

        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            for (int counter = 0; counter < readLikelihoods.sampleEvidence(0).size(); counter++) {
                GATKRead read = readLikelihoods.sampleEvidence(0).get(counter);
                HaplotypeCallerGenotypingDebugger.println("read " + counter + ": " + read.getName() + " cigar: " + read.getCigar() + " mapQ: " + read.getMappingQuality() + " loc: [" + read.getStart() + "-" + read.getEnd() + "] unclippedloc: [" + read.getUnclippedStart() + "-" + read.getUnclippedEnd() + "] length:" + read.getLength());
                String hmmScores = "";
                for (int h = 0; h < readLikelihoods.numberOfAlleles(); h++) {
                    hmmScores = hmmScores + "," + readLikelihoods.sampleMatrix(0).get(h, counter);
                }
                HaplotypeCallerGenotypingDebugger.println(hmmScores);
            }
        }

        final CalledHaplotypes calledHaplotypes = localGenotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                subsettedReadLikelihoodsFinal,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getSpan(),
                features,
                givenAlleles,
                emitReferenceConfidence(),
                hcArgs.maxMnpDistance,
                readsHeader,
                haplotypeBAMWriter.isPresent(),
                suspiciousLocations,
                readLikelihoods);

        if ( haplotypeBAMWriter.isPresent() ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if ( hcArgs.disableOptimizations ) {
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            }
            haplotypeBAMWriter.get().writeReadsAlignedToHaplotypes(haplotypes, assemblyResult.getPaddedReferenceLoc(), haplotypes,
                                                             calledHaplotypeSet, subsettedReadLikelihoodsFinal,regionForGenotyping.getSpan());
        }

        if( hcArgs.assemblerArgs.debugAssembly) {
            logger.info("----------------------------------------------------------------------------------");
        }

        if ( emitReferenceConfidence() ) {
            if ( !containsCalls(calledHaplotypes) ) {
                // no called all of the potential haplotypes
                return referenceModelForNoVariation(region, false, VCpriors);
            }
            else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section, then variant-containing section, then right flank
                trimmingResult.nonVariantLeftFlankRegion().ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank, false, VCpriors)));

                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        subsettedReadLikelihoodsFinal, localGenotypingEngine.getPloidyModel(), calledHaplotypes.getCalls(), hcArgs.standardArgs.genotypeArgs.supportVariants != null,
                        VCpriors));

                trimmingResult.nonVariantRightFlankRegion().ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank, false, VCpriors)));

                return result;
            }
        }
        else {
            //TODO this should be updated once reducible annotations are handled properly.
            return calledHaplotypes.getCalls()
                    .stream()
                    .map(RMSMappingQuality.getInstance()::finalizeRawMQ)
                    .collect(Collectors.toList());
        }
    }

    /*
     * possibly uncollapse the haplotypes inside a read likelihood matrix
     *
     * At this stage, the Haplotypes in the likelihod matrix are derived from the reads. If the reads
     * are flow reads, then they are essentially collapsed, i.e. their maximal hmer size is limited
     * by the flow format's maxHmer. The method uncollapses the haplotypes to be consistant, as much
     * as possible (or applicable) according to the reference.
     */
    private AlleleLikelihoods<GATKRead, Haplotype> possiblyUncollapseHaplotypesInReadLikelihoods(final AssemblyResultSet assemblyResultSet, final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods) {

        if ( assemblyResultSet.getHaplotypeCollapsingEngine() != null ) {

            final List<Haplotype> haplotypes = assemblyResultSet.getHaplotypeCollapsingEngine().uncollapseHmersInHaplotypes(readLikelihoods.alleles(), false, null);
            logger.debug(() -> String.format("%d haplotypes before uncollapsing", haplotypes.size()));
            Map<Haplotype, List<Haplotype>> identicalHaplotypesMap = LongHomopolymerHaplotypeCollapsingEngine.identicalBySequence(haplotypes);
            readLikelihoods.changeAlleles(haplotypes);
            final AlleleLikelihoods<GATKRead, Haplotype>  uncollapsedReadLikelihoods = readLikelihoods.marginalize(identicalHaplotypesMap);
            logger.debug(() -> String.format("%d haplotypes after uncollapsing", uncollapsedReadLikelihoods.numberOfAlleles()));

            return uncollapsedReadLikelihoods;
        } else {
            return readLikelihoods;
        }
    }


    protected boolean containsCalls(final CalledHaplotypes calledHaplotypes) {
        return calledHaplotypes.getCalls().stream()
                .flatMap(call -> call.getGenotypes().stream())
                .anyMatch(Genotype::isCalled);
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    protected List<VariantContext> referenceModelForNoVariation(final AssemblyRegion region, final boolean needsToBeFinalized, final List<VariantContext> VCpriors) {
        final HaplotypeCallerGenotypingEngine localGenotypingEngine = getLocalGenotypingEngine(region);
        if ( emitReferenceConfidence() ) {
            if ( needsToBeFinalized ) {
                AssemblyBasedCallerUtils.finalizeRegion(region,
                        hcArgs.assemblerArgs.errorCorrectReads,
                        hcArgs.dontUseSoftClippedBases,
                        minTailQuality, readsHeader, samplesList,
                        ! hcArgs.doNotCorrectOverlappingBaseQualities,
                        hcArgs.softClipLowQualityEnds,
                        hcArgs.overrideSoftclipFragmentCheck,
                        hcArgs.pileupDetectionArgs.usePileupDetection);
            }

            filterNonPassingReads(region);

            final SimpleInterval paddedLoc = region.getPaddedSpan();
            final Haplotype refHaplotype = AssemblyBasedCallerUtils.createReferenceHaplotype(region, paddedLoc, referenceReader);
            final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
            return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                    paddedLoc, region, AssemblyBasedCallerUtils.createDummyStratifiedReadMap(refHaplotype, samplesList, readsHeader, region),
                    localGenotypingEngine.getPloidyModel(), Collections.emptyList(), hcArgs.standardArgs.genotypeArgs.supportVariants != null, VCpriors);
        }
        else {
            return NO_CALLS;
        }
    }

    /**
     * Shutdown this HC engine, closing resources as appropriate
     */
    public void shutdown() {
        likelihoodCalculationEngine.close();
        if (pdhmmLikelihoodCalculationEngine != null) pdhmmLikelihoodCalculationEngine.close();
        aligner.close();
        haplotypeBAMWriter.ifPresent(HaplotypeBAMWriter::close);
        assembledEventMapVcfOutputWriter.ifPresent(writer -> {assembledEventMapVariants.get().forEach(event -> writer.add(event)); writer.close();});
        if ( referenceReader != null) {
            try {
                referenceReader.close();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.close();
            } catch (IOException e) {
                throw new UserException("Error closing debug output stream", e);
            }
        }
        HaplotypeCallerGenotypingDebugger.close();
        // Write assembly region debug output if present
        assemblyEngine.printDebugHistograms();

    }

    protected Set<GATKRead> filterNonPassingReads( final AssemblyRegion activeRegion ) {
        // TODO: can we unify this additional filtering with makeStandardHCReadFilter()?

        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : activeRegion.getReads() ) {
            if( AlignmentUtils.unclippedReadLength(rec) < READ_LENGTH_FILTER_THRESHOLD ||
                    rec.getMappingQuality() < hcArgs.mappingQualityThreshold ||
                    !ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(rec) ||
                    (hcArgs.keepRG != null && !rec.getReadGroup().equals(hcArgs.keepRG)) ) {
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("Filtered before assembly the read: " + rec);
                }
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    protected void removeReadsFromAllSamplesExcept(final String targetSample, final AssemblyRegion activeRegion) {
        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : activeRegion.getReads() ) {
            if( ! ReadUtils.getSampleName(rec, readsHeader).equals(targetSample) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
    }
}
