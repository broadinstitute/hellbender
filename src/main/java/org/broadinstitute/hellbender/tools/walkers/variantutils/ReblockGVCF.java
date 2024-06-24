package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.collections.Permutation;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.broadinstitute.hellbender.utils.variant.writers.ReblockingGVCFBlockCombiner;
import org.broadinstitute.hellbender.utils.variant.writers.ReblockingGVCFWriter;
import org.broadinstitute.hellbender.utils.variant.writers.ReblockingOptions;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Condense homRef blocks in a single-sample GVCF
 *
 * <p>
 * ReblockGVCF compresses a GVCF by merging hom-ref blocks that were produced using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller according to new GQ band parameters.  Uncalled alleles and associated data will also be dropped unless --keep-all-alts is specified.
 * A joint callset produced with GVCFs reprocessed by ReblockGVCF will have
 * lower precision for hom-ref genotype qualities at variant sites, but the input data footprint can be greatly reduced
 * if the default GQ band parameters are used.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A HaplotypeCaller-produced GVCF to reblock
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A smaller GVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk ReblockGVCF \
 *   -GQB 20 -GQB 30 -GQB 40 --floor-blocks \
 *   -R reference.fasta \
 *   -V sample1.g.vcf \
 *   -O sample1.rb.g.vcf
 * </pre>
 *
 * Invocation as for smallest GVCFs to use with GnarlyGenotyper
 * <pre>
 *  gatk ReblockGVCF \
 *    -R reference.fasta \
 *    -V sample1.g.vcf \
 *    -drop-low-quals \
 *    -rgq-threshold 10 \
 *    -do-qual-approx \
 *    -O sample1.reblocked.g.vcf
 *  * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only single-sample GVCF files produced by HaplotypeCaller can be used as input for this tool.</p>
 * <p>Annotations and header lines that are uninformative for single-sample will be dropped: 
 *       MLEAC, MLEAF, DS, ExcessHet, HaplotypeScore, InbreedingCoeff, AS_InbreedingCoeff
 * <p>Note that when uncalled alleles are dropped, the original GQ may increase.  Use --keep-all-alts if GQ accuracy is a concern.</p>
 *
 */
@CommandLineProgramProperties(summary = "Compress a single-sample GVCF from HaplotypeCaller by merging homRef blocks using new GQ band parameters",
        oneLineSummary = "Condenses homRef blocks in a single-sample GVCF",
        programGroup = OtherProgramGroup.class)
@DocumentedFeature
public final class ReblockGVCF extends MultiVariantWalker {

    private static final Logger logger = LogManager.getLogger(ReblockGVCF.class);
    private static final OneShotLogger genotypeLogger = new OneShotLogger(ReblockGVCF.class);

    public static final String DROP_LOW_QUALS_ARG_NAME = "drop-low-quals";
    public static final String RGQ_THRESHOLD_LONG_NAME = "rgq-threshold-to-no-call";
    public static final String RGQ_THRESHOLD_SHORT_NAME = "rgq-threshold";
    public static final String TREE_SCORE_THRESHOLD_LONG_NAME = "tree-score-threshold-to-no-call";
    public static final String ANNOTATIONS_TO_KEEP_LONG_NAME = "annotations-to-keep";
    public static final String ANNOTATIONS_TO_REMOVE_LONG_NAME = "format-annotations-to-remove";
    public static final String KEEP_ALL_ALTS_ARG_NAME = "keep-all-alts";
    public static final String QUAL_APPROX_LONG_NAME = "do-qual-score-approximation";
    public static final String QUAL_APPROX_SHORT_NAME = "do-qual-approx";
    public static final String ALLOW_MISSING_LONG_NAME = "allow-missing-hom-ref-data";
    public static final String KEEP_SITE_FILTERS_LONG_NAME = "keep-site-filters";
    public static final String KEEP_SITE_FILTERS_SHORT_NAME = "keep-filters";
    public static final String ADD_FILTERS_TO_GENOTYPE = "add-site-filters-to-genotype";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    protected GATKPath outputFile;

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Output the band lower bound for each GQ block instead of the min GQ -- for better compression
     * Note that this argument also drops PLS for more efficient storage
     */
    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS, doc = "Output the band lower bound for each GQ block regardless of the data it represents", optional = true)
    private boolean floorBlocks = false;

    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.GQ_BAND_LONG_NAME, shortName=HaplotypeCallerArgumentCollection.GQ_BAND_SHORT_NAME,
            doc="Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)", optional = true)
    public List<Integer> GVCFGQBands = new ArrayList<>();
    {
        GVCFGQBands.add(20); GVCFGQBands.add(100);
    }

    @Advanced
    //@DeprecatedFeature(detail="This argument introduces 'holes', resulting in an invalid GVCF")
    @Argument(fullName=DROP_LOW_QUALS_ARG_NAME, shortName=DROP_LOW_QUALS_ARG_NAME, doc="Exclude variants and homRef blocks that are GQ0 from the reblocked GVCF to save space; drop low quality/uncalled alleles", optional = true)
    protected boolean dropLowQuals = false;

    @Advanced
    @Argument(fullName=RGQ_THRESHOLD_LONG_NAME, shortName=RGQ_THRESHOLD_SHORT_NAME, doc="Reference genotype quality (PL[0]) value below which variant sites will be converted to GQ0 homRef calls", optional = true)
    protected double rgqThreshold = 0.0;

    @Advanced
    @Argument(fullName=TREE_SCORE_THRESHOLD_LONG_NAME, doc="Tree score value below which variant sites will be converted to GQ0 homRef calls. Disabled when set to 0 (which is the default).", optional = true)
    private double treeScoreThreshold = 0.0;

    @Advanced
    @Argument(fullName=ANNOTATIONS_TO_KEEP_LONG_NAME, doc="Annotations that are not recognized by GATK to keep, that should be kept in final GVCF at variant sites.", optional = true)
    private List<String> annotationsToKeep = new ArrayList<>();

    @Advanced
    @Argument(fullName=ANNOTATIONS_TO_REMOVE_LONG_NAME, doc="FORMAT level annotations to remove from all genotypes in final GVCF.", optional = true)
    private List<String> annotationsToRemove = new ArrayList<>();

    @Advanced
    @Argument(fullName=QUAL_APPROX_LONG_NAME, shortName=QUAL_APPROX_SHORT_NAME, doc="Add necessary INFO field annotation to perform QUAL approximation downstream; required for GnarlyGenotyper", optional = true)
    protected boolean doQualApprox = false;

    @Advanced
    @Argument(fullName=ALLOW_MISSING_LONG_NAME, doc="Fill in homozygous reference genotypes with no PLs and no GQ with PL=[0,0,0].  Necessary for input from Regeneron's WeCall variant caller.", optional = true)
    protected boolean allowMissingHomRefData = false;

    @Advanced
    @Argument(fullName=KEEP_ALL_ALTS_ARG_NAME, doc="Keep all ALT alleles and full PL array for most accurate GQs", optional = true)
    protected boolean keepAllAlts = false;

    @Advanced
    @Argument(fullName=KEEP_SITE_FILTERS_LONG_NAME, shortName = KEEP_SITE_FILTERS_SHORT_NAME, doc="Keep site level filters for variants (not ref blocks).")
    private boolean keepFilters = false;

    @Advanced
    @Argument(fullName= ADD_FILTERS_TO_GENOTYPE, doc="Add site level filters to genotype level. Site level filters removed by default, if they should be kept, use --" + KEEP_SITE_FILTERS_LONG_NAME)
    private boolean addFiltersToFormatField = false;

    //TODO: this will be an argument when posteriors handling is fully implemented in AlleleSubsettingUtils
    protected String posteriorsKey = null;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;
    // the INFO field annotation key names to remove
    private static final List<String> infoFieldAnnotationKeyNamesToRemove = Arrays.asList(GVCFWriter.GVCF_BLOCK, GATKVCFConstants.HAPLOTYPE_SCORE_KEY,
            GATKVCFConstants.INBREEDING_COEFFICIENT_KEY, GATKVCFConstants.MLE_ALLELE_COUNT_KEY,
            GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, GATKVCFConstants.EXCESS_HET_KEY, GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY,
            GATKVCFConstants.DOWNSAMPLED_KEY);

    private CachingIndexedFastaSequenceFile referenceReader;

    private static final List<String> alleleBasedLengthAnnots = new ArrayList<>();

    public static class AlleleLengthComparator implements Comparator<Allele> {
        public int compare(Allele a1, Allele a2) {
            return a1.getBaseString().length() - a2.getBaseString().length();
        }
    }

    @VisibleForTesting
    ReblockingGVCFWriter vcfWriter;

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, AS_StandardAnnotation.class);
    }

    @Override
    public boolean requiresReference() {return true;}

    @Override
    public void onTraversalStart() {
        if (getSamplesForVariants().size() != 1) {
            throw new UserException.BadInput("ReblockGVCF can take multiple input GVCFs, but they must be "
                    + "non-overlapping shards from the same sample.  Found samples " + getSamplesForVariants());
        }

        final VCFHeader inputHeader = getHeaderForVariants();

        if (treeScoreThreshold > 0 && inputHeader.getInfoHeaderLine(GATKVCFConstants.TREE_SCORE) == null) {
            throw new UserException("-" + TREE_SCORE_THRESHOLD_LONG_NAME + " is set to value greater than 0: " + treeScoreThreshold
                + ", but the " + GATKVCFConstants.TREE_SCORE + " annotation is not present in the input GVCF.");
        }

        List<String> missingAnnotationsToRemove = annotationsToRemove.stream().filter(a -> inputHeader.getFormatHeaderLine(a)==null).toList();
        missingAnnotationsToRemove.forEach(a -> logger.warn("FORMAT level annotation " + a + ", which was requested to be removed by --" + ANNOTATIONS_TO_REMOVE_LONG_NAME + ", not found in input GVCF header."));

        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeaders);
        // Remove GCVFBlocks, legacy headers, and annotations that aren't informative for single samples
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCFWriter.GVCF_BLOCK) ||
                (vcfHeaderLine.getKey().equals("INFO")) && ((VCFInfoHeaderLine)vcfHeaderLine).getID().equals(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED) ||  //remove old (maybe wrong type) and add new with deprecated note
                (vcfHeaderLine.getKey().equals("INFO")) && infoFieldAnnotationKeyNamesToRemove.contains(((VCFInfoHeaderLine)vcfHeaderLine).getID()));

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        genotypingEngine = createGenotypingEngine(new IndexedSampleList(getSamplesForVariants()));
        createAnnotationEngine();

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED));  //NOTE: this is deprecated, but keep until we reprocess all GVCFs
        if (inputHeader.hasInfoLine(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED));
        }

        for(String annotation : annotationsToKeep) {
            if(inputHeader.hasInfoLine(annotation)) {
                VCFInfoHeaderLine headerLine = inputHeader.getInfoHeaderLine(annotation);
                if (headerLine.getCountType().equals(VCFHeaderLineCount.A)) {
                    throw new UserException(String.format("%s is an allele specific annotation which is currently unsupported. Only fixed length annotations may be kept in reblocking.", annotation));
                }
                headerLines.add(inputHeader.getInfoHeaderLine(annotation));
            } else {
                throw new UserException(String.format("%s is not in header of input GVCF but was requested to be kept by %s argument.", annotation, ANNOTATIONS_TO_KEEP_LONG_NAME));
            }
        }

        if (addFiltersToFormatField) {
            headerLines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
        }

        //Allele length and Genotype length annotations need to be subset or removed if alleles are dropped so we need to parse the header for annotation count types
        for(final VCFFormatHeaderLine formatHeaderLine : inputHeader.getFormatHeaderLines()) {
            if(formatHeaderLine.getCountType().equals(VCFHeaderLineCount.A) ||
                    formatHeaderLine.getCountType().equals(VCFHeaderLineCount.R) ||
                    formatHeaderLine.getCountType().equals(VCFHeaderLineCount.G)) {
                alleleBasedLengthAnnots.add(formatHeaderLine.getID());
            }
        }

        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        referenceReader = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());

        createVcfWriter(headerLines);
    }

    @VisibleForTesting
    public void createVcfWriter(Set<VCFHeaderLine> headerLines) {
        final VariantContextWriter writer = createVCFWriter(outputFile);

        final ReblockingOptions reblockingOptions = new ReblockingOptions(dropLowQuals, allowMissingHomRefData, rgqThreshold);

        try {
            vcfWriter = new ReblockingGVCFWriter(writer, new ArrayList<>(GVCFGQBands), floorBlocks, referenceReader, reblockingOptions);
        } catch ( final IllegalArgumentException e ) {
            throw new UserException.BadInput("GQBands are malformed: " + e.getMessage(), e);
        }
        vcfWriter.writeHeader(new VCFHeader(headerLines, getSamplesForVariants()));  //don't get samples from header -- multi-variant inputHeader doens't have sample names
    }

    private HaplotypeCallerGenotypingEngine createGenotypingEngine(final SampleList samples) {
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        // create the genotyping engine
        hcArgs.standardArgs.outputMode = OutputMode.EMIT_ALL_CONFIDENT_SITES;  //use confident vs. active mode so we can drop low quality variants
        hcArgs.standardArgs.annotateAllSitesWithPLs = true;
        hcArgs.standardArgs.genotypeArgs = genotypeArgs.clone();
        hcArgs.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;   //this is important to force emission of all alleles at a multiallelic site
        hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling = dropLowQuals ? genotypeArgs.standardConfidenceForCalling : 0.0;
        return new HaplotypeCallerGenotypingEngine(hcArgs, samples, true, false);

    }

    @VisibleForTesting
    protected void createAnnotationEngine() {
        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, false);
    }

    // get VariantContexts from input gVCFs and regenotype
    @Override
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        if (!variant.hasAllele(Allele.NON_REF_ALLELE)) {
            throw new UserException("Variant Context at " + variant.getContig() + ":" + variant.getStart() + " does not contain a <NON-REF> allele. This tool is only intended for use with GVCFs.");
        }
        VariantContext newVC = annotationsToRemove.size() > 0 ? removeVCFFormatAnnotations(variant) : variant;
        regenotypeVC(newVC);
    }

    /**
     * Remove format level annotations from genotype in variant context.
     *
     * @param vc variant context to remove format annotations from
     * @return variant context with format annotations removed from genotype
     */
    private VariantContext removeVCFFormatAnnotations(final VariantContext vc) {
        final Genotype genotype = vc.getGenotype(0);
        Map<String, Object> extendedAttributes = genotype.getExtendedAttributes();
        // If extendedAttributes is empty, we get an unmodifiable empty map, so skip removing any annotations
        if (extendedAttributes.isEmpty()) {
            return vc;
        }

        for (String annotation : annotationsToRemove) {
            extendedAttributes.remove(annotation);
        }
        final Genotype newGenotype = new GenotypeBuilder(genotype).noAttributes().attributes(extendedAttributes).make();
        return new VariantContextBuilder(vc).genotypes(newGenotype).make();
    }

    /**
     * Re-genotype (and re-annotate) a VariantContext, adding it to the VCF writer if necessary
     * Note that the GVCF write takes care of the actual homRef block merging based on {@code GVCFGQBands}
     *
     * @param originalVC     the combined genomic VC
     */
    private void regenotypeVC(final VariantContext originalVC) {

        //Pass back ref-conf homRef sites/blocks to be combined by the GVCFWriter
        if (isHomRefBlock(originalVC)) {
            //if this hom ref block is entirely overlapped by previous VCF output, then drop it
            if (originalVC.contigsMatch(vcfWriter.getVcfOutputEnd()) && originalVC.getEnd() <= vcfWriter.getVcfOutputEnd().getStart()) {
                return;
            }
            final Genotype genotype = originalVC.getGenotype(0);
            if (dropLowQuals && (!genotype.hasGQ() || genotype.getGQ() < rgqThreshold || genotype.getGQ() == 0)) {
                return;
            }
            if (!genotype.hasPL()) {
                if (genotype.hasGQ()) {
                    logger.warn("PL is missing for hom ref genotype at at least one position for sample " + genotype.getSampleName() + ": " + originalVC.getContig() + ":" + originalVC.getStart() +
                            ".  Using GQ to determine quality.");
                    vcfWriter.add(originalVC);
                } else {
                    final String message = "Homozygous reference genotypes must contain GQ or PL. Both are missing for hom ref genotype at "
                            + originalVC.getContig() + ":" + originalVC.getStart();
                    if (allowMissingHomRefData) {
                        logger.warn(message);
                        final VariantContextBuilder vcBuilder = new VariantContextBuilder(originalVC);
                        final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                        vcBuilder.genotypes(gBuilder.GQ(0).PL(new int[]{0,0,0}).make());
                        vcfWriter.add(vcBuilder.make());
                    } else {
                        throw new UserException.BadInput(message);
                    }
                }
            }
            vcfWriter.add(originalVC);
            return;
        }

        VariantContext result = originalVC;

        //Use the genotyping engine to do the QUAL thresholding if we're dropping low qual sites
        //don't need to calculate quals for sites with no data whatsoever or sites already genotyped homRef,
        //but if STAND_CALL_CONF > 0 we need to drop low quality alleles and regenotype
        //Note that spanning deletion star alleles will be considered low quality
        if (dropLowQuals && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 0 && !isMonomorphicCallWithAlts(originalVC)) {
            final VariantContext regenotyped = genotypingEngine.calculateGenotypes(originalVC);
            if (regenotyped == null) {
                return;
            }
            //make sure result has annotations so we don't have to keep originalVC around
            result = new VariantContextBuilder(regenotyped).attributes(subsetAnnotationsIfNecessary(annotationEngine, doQualApprox, posteriorsKey, originalVC, regenotyped, annotationsToKeep)).make();
        }


        //variants with PL[0] less than threshold get turned to homRef with PL=[0,0,0], shouldn't get INFO attributes
        //make sure we can call het variants with GQ >= rgqThreshold in joint calling downstream
        if(shouldBeReblocked(result)) {
            if (vcfWriter.getVcfOutputEnd() != null && result.contigsMatch(vcfWriter.getVcfOutputEnd()) && result.getEnd() <= vcfWriter.getVcfOutputEnd().getEnd()) {
                //variant is entirely overlapped by variants already output to the VCF, so drop it
                return;
            }
            final VariantContextBuilder newHomRefBuilder = lowQualVariantToGQ0HomRef(result);
            if (newHomRefBuilder != null) {  //can be null if we're dropping low quals
                vcfWriter.add(newHomRefBuilder.make());
            }
        }
        //high quality variant
        else {
            final VariantContext trimmedVariant = cleanUpHighQualityVariant(result);
            vcfWriter.add(trimmedVariant);
        }
    }

    /**
     *  Make sure annotations are updated if alleles are dropped
     * @param doQualApprox  output a QUALapprox INFO attribute?
     * @param posteriorsKey if not null, the key for genotype posteriors in this VCF version
     * @param originalVC    the original variant context with the full set of alleles
     * @param regenotyped   should have a strict subset of alleles in originalVC, untrimmed
     * @return  a set of annotations that correspond to the (possibly subset) alleles in the regenotyped VC
     */
    @VisibleForTesting
    static Map<String, Object> subsetAnnotationsIfNecessary(final VariantAnnotatorEngine annotationEngine,
                                                            final boolean doQualApprox, final String posteriorsKey,
                                                            final VariantContext originalVC, final VariantContext regenotyped,
                                                            final List<String> annotationsToKeep) {
        final Map<String, Object> newAnnotations;
        if (regenotyped.getNAlleles() != originalVC.getNAlleles()) {
            final Permutation<Allele> allelePermutation = new IndexedAlleleList<>(originalVC.getAlleles()).
                    permutation(new IndexedAlleleList<>(regenotyped.getAlleles()));
            final int[] relevantIndices = IntStream.range(0, regenotyped.getAlleles().size())
                    .map(n -> allelePermutation.fromIndex(n)).toArray();
            newAnnotations = new LinkedHashMap<>();
            composeUpdatedAnnotations(newAnnotations, doQualApprox, posteriorsKey, originalVC, annotationEngine, relevantIndices, regenotyped, annotationsToKeep);
        } else {
            newAnnotations = originalVC.getAttributes();
        }
        return newAnnotations;
    }

    /**
     * determine if a VC is a homRef block, i.e. has only a NON_REF alt allele
     * @param result VariantContext to process
     * @return true if VC is a homRef block and not a "call" with alt alleles and annotations
     */
    public static boolean isHomRefBlock(final VariantContext result) {
        return (result.getAlternateAlleles().size() == 1) && result.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);
    }

    /**
     * determine if VC is a homRef "call", i.e. an annotated variant with non-symbolic alt alleles and homRef genotypes
     * we treat these differently from het/homVar calls or homRef blocks
     * @param result VariantContext to process
     * @return true if VC is a 0/0 call and not a homRef block
     */
    private boolean isMonomorphicCallWithAlts(final VariantContext result) {
        final Genotype genotype = result.getGenotype(0);
        return (hasGenotypeValuesArray(posteriorsKey, genotype)
                && (genotype.isHomRef() || isNoCalledHomRef(posteriorsKey, genotype) || (genotype.hasPL() && MathUtils.minElementIndex(genotype.getPL()) == 0))
                && result.getAlternateAlleles().stream().anyMatch(this::isConcreteAlt));
    }

    /**
     *
     * @param posteriorsKey FORMAT field tag for genotype posterior probabilities, may be null
     * @param genotype  a genotype that may or may not have likelihood/probability data
     * @return true if we can quantify genotype call
     */
    private static boolean hasGenotypeValuesArray(final String posteriorsKey, final Genotype genotype) {
        return genotype.hasPL() || (posteriorsKey != null && genotype.hasExtendedAttribute(posteriorsKey));
    }

    private static boolean isNoCalledHomRef(final String posteriorsKey, final Genotype genotype) {
        return genotype.isNoCall()
                && hasGenotypeValuesArray(posteriorsKey, genotype)
                && getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey)[0] == 0;
    }

    /**
     * Should this variant context be turned into a reference block?
     * @param vc    a low quality variant or variant called homozygous reference
     * @return  true if this VariantContext is eligible to be combined with adjacent reference blocks
     */
    @VisibleForTesting
    boolean shouldBeReblocked(final VariantContext vc) {
        if (!vc.hasGenotypes()) {
            throw new IllegalStateException("Variant contexts must contain genotypes to be reblocked.");
        }
        final Genotype genotype = vc.getGenotype(0);
        final int[] pls = getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey);
        if (pls == null) {
            return true;
        }
        final int minLikelihoodIndex = MathUtils.minElementIndex(pls);
        final GenotypeAlleleCounts alleleCounts = GenotypesCache.get(genotype.getPloidy(), minLikelihoodIndex);

        final List<Allele> finalAlleles = alleleCounts.asAlleleList(vc.getAlleles());
        return (pls != null && pls[0] < rgqThreshold)
                || !genotypeHasConcreteAlt(finalAlleles)
                || finalAlleles.stream().anyMatch(a -> a.equals(Allele.NON_REF_ALLELE))
                || (!genotype.hasPL() && !genotype.hasGQ()
                // If TREE_SCORE threshold is set to be > 0 then sites without TREE_SCORE should be reblocked.
                || vc.getAttributeAsDouble(GATKVCFConstants.TREE_SCORE, 0) < treeScoreThreshold);
    }

    /**
     * Is there a real ALT allele that's not <NON_REF> or * ?
     * @param alleles   alleles in the called genotype
     * @return true if the genotype has a called allele that is a "real" alternate
     */
    private boolean genotypeHasConcreteAlt(final List<Allele> alleles) {
        return alleles.stream().anyMatch(this::isConcreteAlt);
    }

    private boolean isConcreteAlt(final Allele a) {
        return !a.isReference() && !a.isSymbolic() && !a.equals(Allele.SPAN_DEL);
    }

    /**
     * "reblock" a variant by converting its genotype to homRef, changing PLs, adding reblock END tags and other attributes
     * @param lowQualityVariant  a variant already determined to be low quality
     * @return a Builder that can be modified later, may be null
     */
    @VisibleForTesting
    public VariantContextBuilder lowQualVariantToGQ0HomRef(final VariantContext lowQualityVariant) {
        if(dropLowQuals && (!isMonomorphicCallWithAlts(lowQualityVariant) || !lowQualityVariant.getGenotype(0).isCalled())) {
            return null;
        }

        final Map<String, Object> attrMap = new LinkedHashMap<>();

        //this method does a lot of things, including fixing alleles and adding the END key
        final GenotypeBuilder gb = changeCallToHomRefVersusNonRef(lowQualityVariant, attrMap);  //note that gb has all zero PLs

        final VariantContextBuilder builder = new VariantContextBuilder(lowQualityVariant);

        final Genotype newG = gb.make();
        builder.alleles(Arrays.asList(newG.getAlleles().get(0), Allele.NON_REF_ALLELE)).genotypes(newG);
        if (vcfWriter.getVcfOutputEnd() != null && lowQualityVariant.contigsMatch(vcfWriter.getVcfOutputEnd()) && lowQualityVariant.getStart() <= vcfWriter.getVcfOutputEnd().getStart()) {
            final int newStart = vcfWriter.getVcfOutputEnd().getEnd() + 1;
            if (newStart > lowQualityVariant.getEnd()) {
                return null;
            }
            ReblockingGVCFBlockCombiner.moveBuilderStart(builder, newStart, referenceReader);
        }
        return builder.unfiltered()  //genotyping engine will add lowQual filter, so strip it off
                .log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap);
    }

    /**
     * Produce a GenotypeBuilder for a hom-ref call suitable to be merged into a reference block, i.e. set with PLs as
     * appropriate for a variant with only reference and <NON_REF> as alleles and END applied as appropriate
     * Note that this may modify {@code attrMap} as a side effect for END key
     * @param lowQualVariant a VC containing a genotype to be converted to a GQ0 homRef call; needed for alleles that correspond to PLs and other things
     * @param attrMap the new VC attribute map, to update the END tag for deletions
     * @return a GenotypeBuilder to make a 0/0 call with updated PLs and GQ
     */
    @VisibleForTesting
    protected GenotypeBuilder changeCallToHomRefVersusNonRef(final VariantContext lowQualVariant, final Map<String, Object> attrMap) {
        final Genotype genotype = lowQualVariant.getGenotype(0);
        final Allele inputRefAllele = lowQualVariant.getReference();
        GenotypeBuilder gb = new GenotypeBuilder(genotype);
        //if GT is not homRef, correct it and set GQ=0
        if (posteriorsKey == null && (!genotype.hasPL() || (genotype.hasPL() && genotype.getPL()[0] != 0))) {
            gb.PL(new int[GenotypeLikelihoods.numLikelihoods(2, genotype.getPloidy())]);  //2 alleles for ref and non-ref
            gb.GQ(0).noAD().alleles(Collections.nCopies(genotype.getPloidy(), inputRefAllele)).noAttributes();
        //for hom-ref variants, drop other ALTs and subset PLs, GQ is recalculated (may be > 0)
        } else {
            if (posteriorsKey != null && genotype.hasExtendedAttribute(posteriorsKey)) {
                subsetHomRefPosteriorsToRefVersusNonRef(lowQualVariant, gb);
            } else {
                //find best ALT so we can use its likelihood for NON_REF
                final List<Allele> bestAlleles = AlleleSubsettingUtils.calculateMostLikelyAlleles(lowQualVariant, genotype.getPloidy(), 1, true);
                final Allele bestAlt = bestAlleles.stream().filter(a -> !a.isReference()).findFirst().orElse(Allele.NON_REF_ALLELE);  //allow span dels
                //we care about the best alt even though it's getting removed because NON_REF should get the best likelihoods
                //it shouldn't matter that we're passing in different alt alleles since the GenotypesContext only knows
                // the called alleles and this is a reference genotype that will stay hom-ref
                final GenotypesContext context = AlleleSubsettingUtils.subsetAlleles(lowQualVariant.getGenotypes(),
                        genotype.getPloidy(), lowQualVariant.getAlleles(), Arrays.asList(inputRefAllele, bestAlt),
                        null, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);  //BEST_MATCH to avoid no-calling low qual genotypes
                final Genotype subsetG = context.get(0);
                gb = new GenotypeBuilder(subsetG).noAttributes();  //remove attributes because hom ref blocks shouldn't have posteriors
                //subsetting may strip GQ and PLs for low qual genotypes
                if (!subsetG.hasGQ()) {
                    gb.GQ(0);
                }
                if (!subsetG.hasPL()) {
                    gb.PL(new int[GenotypeLikelihoods.numLikelihoods(2, genotype.getPloidy())]);  //2 alleles for ref and non-ref
                }
            }
        }
        if (lowQualVariant.hasAttribute(VCFConstants.DEPTH_KEY)) {
            final int depth = lowQualVariant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            gb.DP(depth);
            gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, depth);
        } else if (genotype.hasAD()) {
            final int depth = (int) MathUtils.sum(genotype.getAD());
            gb.DP(depth);
            gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, depth);
        }
        //NOTE: If we're dropping a deletion allele, then we need to trim the reference and add an END tag with the vc stop position
        final Allele outputRefAllele;
        if (inputRefAllele.length() > 1 || genotype.getAlleles().contains(Allele.SPAN_DEL) || genotype.getAlleles().contains(Allele.NO_CALL)) {
            outputRefAllele = Allele.create(inputRefAllele.getBases()[0], true);
        } else {
            outputRefAllele = inputRefAllele;
        }
        attrMap.put(VCFConstants.END_KEY, lowQualVariant.getEnd());
        gb.alleles(Collections.nCopies(genotype.getPloidy(), outputRefAllele));
        return gb;
    }

    /**
     * Subset alleles as necessary and apply annotations
     * @param variant    VariantContext with full set of annotations (e.g. DP)
     * @return  an annotated VariantContext with data only for ref, non-ref and called alts
     */
    @VisibleForTesting
    VariantContext cleanUpHighQualityVariant(final VariantContext variant) {
        final Map<String, Object> attrMap = new LinkedHashMap<>();

        final Genotype genotype = getCalledGenotype(variant);
        VariantContextBuilder builder = new VariantContextBuilder(variant);  //QUAL from result is carried through
        builder.attributes(attrMap).genotypes(genotype);  //clear attributes

        final List<Allele> allelesToDrop = getAllelesToDrop(variant, genotype);

        final boolean allelesNeedSubsetting = !allelesToDrop.isEmpty();
        int[] relevantIndices = new int[variant.getNAlleles()];  //called alleles plus ref and non-ref
        final List<Allele> newAlleleSetUntrimmed = new ArrayList<>(variant.getAlleles());
        if(allelesNeedSubsetting && !keepAllAlts) {
            newAlleleSetUntrimmed.removeAll(allelesToDrop);
            final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(variant.getGenotypes(), genotype.getPloidy(), variant.getAlleles(),
                    newAlleleSetUntrimmed, null, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, alleleBasedLengthAnnots);
            if (gc.get(0).isHomRef() || !gc.get(0).hasGQ() || gc.get(0).getAlleles().contains(Allele.NO_CALL)) {  //could be low quality or no-call after subsetting
                if (dropLowQuals) {
                    return null;
                }
                final VariantContextBuilder newHomRefBuilder = lowQualVariantToGQ0HomRef(variant);
                if (newHomRefBuilder != null) {  //there's a chance the low quality variant may be entirely overlapped by a variant already output
                    vcfWriter.add(newHomRefBuilder.make());
                }
                return null;
            }
            //note that subsetting alleles can increase GQ, e.g. with one confident reference allele and a deletion allele that's either 4 or 5 bases long
            builder.genotypes(gc).alleles(newAlleleSetUntrimmed);
            //if deletions are dropped, alleles may need trimming
            final VariantContext newTrimmedAllelesVC = GATKVariantContextUtils.trimAlleles(builder.make(), false, true);
            builder = new VariantContextBuilder(newTrimmedAllelesVC);
            //save indices of new alleles for annotation processing
            relevantIndices = newAlleleSetUntrimmed.stream().mapToInt(a -> variant.getAlleles().indexOf(a)).toArray();
            final int refBlockDepth;
            if (variant.hasAttribute(VCFConstants.DEPTH_KEY)) {  //prefer INFO depth because HaplotypeCaller GVCF block code uses all reads, not just informative
                refBlockDepth = variant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            } else if (genotype.hasDP()) {
                refBlockDepth = genotype.getDP();
            } else {
                refBlockDepth = 0;
            }
            addRefBlockIfNecessary(variant, allelesToDrop, newTrimmedAllelesVC, refBlockDepth);
        }

        final VariantContext updatedAllelesVC = builder.make();
        final Genotype updatedAllelesGenotype = updatedAllelesVC.getGenotype(0);

        //remove any AD reads for the non-ref
        final Genotype g = removeNonRefADs(updatedAllelesGenotype, updatedAllelesVC.getAlleleIndex(Allele.NON_REF_ALLELE));
        if (addFiltersToFormatField) {
            final List<String> filters = updatedAllelesVC.getFilters().stream().toList();
            builder.genotypes(new GenotypeBuilder(g).filters(filters).make());
        } else {
            builder.genotypes(g);
        }

        composeUpdatedAnnotations(attrMap, doQualApprox, posteriorsKey, variant, annotationEngine, relevantIndices, updatedAllelesVC, annotationsToKeep);

        return keepFilters ? builder.attributes(attrMap).make() : builder.attributes(attrMap).unfiltered().make();
    }

    /**
     * If the ref allele is trimmed after alt deletions are dropped, add a reference block to account for the space covered before trimming
     * @param originalVC    VC with full set of alleles that may need to be trimmed
     * @param allelesToDrop alleles eligible to become the new non-ref likelihood
     * @param newTrimmedAllelesVC   VC with called alleles that may have been trimmed, i.e. differ from originalVC alleles
     */
    private void addRefBlockIfNecessary(final VariantContext originalVC, final List<Allele> allelesToDrop, final VariantContext newTrimmedAllelesVC, final int refBlockDepth) {
        //if deletion needs trimming, fill in the gap with a ref block
        final int oldLongestAlleleLength = originalVC.getReference().length();
        final int newLongestAlleleLength = newTrimmedAllelesVC.getReference().length();
        final Genotype genotype = originalVC.getGenotype(0);
        final int vcfOutputEnd = vcfWriter.getVcfOutputEnd() == null ? -1 : vcfWriter.getVcfOutputEnd().getEnd();
        if (newLongestAlleleLength < oldLongestAlleleLength) {
            //need to add a ref block to make up for the allele trimming or there will be a hole in the GVCF
            final int[] originalLikelihoods = getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey);
            if (originalLikelihoods != null) {
                final Allele oldShortestAltAllele;
                try {
                    oldShortestAltAllele = allelesToDrop.stream().filter(a -> !a.equals(Allele.SPAN_DEL))
                            .min(new AlleleLengthComparator()).orElseThrow(NoSuchElementException::new);
                } catch (final Exception e) {
                    throw new GATKException("No shortest ALT at " + originalVC.getStart() + " across alleles: " + allelesToDrop);
                }

                //subset PLs to ref and longest dropped allele (longest may not be most likely, but we'll approximate so we don't have to make more than one ref block)
                final int[] longestVersusRefPLIndices = AlleleSubsettingUtils.subsettedPLIndices(originalVC.getGenotype(0).getPloidy(),
                        originalVC.getAlleles(), Arrays.asList(originalVC.getReference(), oldShortestAltAllele));
                final int[] newRefBlockLikelihoods = MathUtils.normalizePLs(Arrays.stream(longestVersusRefPLIndices)
                        .map(idx -> originalLikelihoods[idx]).toArray());
                if (newRefBlockLikelihoods[0] != 0) {
                    for (int i = 0; i < newRefBlockLikelihoods.length; i++) {
                        newRefBlockLikelihoods[i] = Math.max(newRefBlockLikelihoods[i] - newRefBlockLikelihoods[0], 0);
                    }
                }

                //build the new reference block with updated likelihoods
                final GenotypeBuilder refBlockGenotypeBuilder = new GenotypeBuilder();
                final int refStart = Math.max(originalVC.getEnd() - (oldLongestAlleleLength - newLongestAlleleLength), vcfOutputEnd) + 1;
                final Allele newRef = Allele.create(ReferenceUtils.getRefBaseAtPosition(referenceReader, originalVC.getContig(), refStart), true);
                refBlockGenotypeBuilder.PL(newRefBlockLikelihoods)
                        .GQ(MathUtils.secondSmallestMinusSmallest(newRefBlockLikelihoods, 0))
                        .alleles(Arrays.asList(newRef, newRef)).DP(refBlockDepth);

                //add the new block to the buffer if it isn't covered by positions already output
                if (refStart > vcfOutputEnd && originalVC.getEnd() > vcfOutputEnd) {
                    final VariantContextBuilder trimBlockBuilder = new VariantContextBuilder();
                    trimBlockBuilder.chr(originalVC.getContig()).start(Math.max(refStart, vcfOutputEnd + 1)).stop(originalVC.getEnd()).
                            alleles(Arrays.asList(newRef, Allele.NON_REF_ALLELE)).attribute(VCFConstants.END_KEY, originalVC.getEnd())
                            .genotypes(refBlockGenotypeBuilder.make());
                    vcfWriter.add(trimBlockBuilder.make());
                }
            }
        }
    }

    /**
     * Update annotations if alleles were subset and add annotations specific to ReblockGVCF, like QUALapprox and RAW_GT_COUNT
     * @param destination   annotation map to modify with new annotations
     * @param doQualApprox  output a QUALapprox INFO attribute?
     * @param posteriorsKey if not null, the key for genotype posteriors in this VCF version
     * @param variant   variant context with full annotation data
     * @param annotationEngine  used to get the list of annotations
     * @param relevantIndices   indexes of alleles in updatedAllelesVC with respect to variant
     * @param updatedAllelesVC  variant context with final set of alleles
     */
    private static void composeUpdatedAnnotations(final Map<String, Object> destination, final boolean doQualApprox, final String posteriorsKey, final VariantContext variant, final VariantAnnotatorEngine annotationEngine,
                                                  final int[] relevantIndices, final VariantContext updatedAllelesVC, final List<String> annotationsToKeep) {
        updateMQAnnotations(destination, variant);

        final boolean allelesNeedSubsetting = relevantIndices.length < variant.getNAlleles();
        copyInfoAnnotations(destination, variant, infoFieldAnnotationKeyNamesToRemove, annotationEngine, allelesNeedSubsetting, relevantIndices);

        //generate qual annotations after we potentially drop alleles
        final Genotype updatedAllelesGenotype = updatedAllelesVC.getGenotype(0);
        if (doQualApprox) {
            if (hasGenotypeValuesArray(posteriorsKey, updatedAllelesGenotype)) {
                addQualAnnotations(destination, posteriorsKey, annotationEngine, updatedAllelesVC);
            }
        } else {  //manually copy annotations that might be from reblocking and aren't part of AnnotationEngine
            if (variant.hasAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY)) {
                destination.put(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, variant.getAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
            }
            if (variant.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY)) {
                destination.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, variant.getAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
            }
        }
        destination.put(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, updatedAllelesGenotype.getAlleles().stream().anyMatch(Allele::isReference) ?
                Arrays.asList(0,1,0) : Arrays.asList(0,0,1)); //ExcessHet currently uses rounded/integer genotype counts, so do the same here
        for (String annotation : annotationsToKeep) {
            if (variant.hasAttribute(annotation)) {
                destination.put(annotation, variant.getAttribute(annotation));
            }
        }
    }

    /**
     * Return the genotype from variant, call it if data is present and GT is no-call
     * @param variant   VariantContext with genotype data
     * @return  a called genotype (if possible)
     */
    private Genotype getCalledGenotype(final VariantContext variant) {
        final boolean hasPLAndPosteriorMismatch;
        final Genotype origG = variant.getGenotype(0);
        final int[] pls = getGenotypePosteriorsOtherwiseLikelihoods(origG, posteriorsKey);
        if (pls == null) {
            throw new IllegalStateException("Cannot verify called genotype without likelihoods or posteriors.  Error at "
                    + variant.getContig() + ":" + variant.getStart());
        }
        final int minLikelihoodIndex = MathUtils.minElementIndex(pls);
        final GenotypeAlleleCounts alleleCounts = GenotypesCache.get(origG.getPloidy(), minLikelihoodIndex);

        final List<Allele> finalAlleles = alleleCounts.asAlleleList(variant.getAlleles());
        hasPLAndPosteriorMismatch = !finalAlleles.containsAll(origG.getAlleles());

        if (origG.isNoCall() || hasPLAndPosteriorMismatch) {  //might be homRef if posteriors and PLs don't agree
            final Genotype noCallGT = variant.getGenotype(0);
            final GenotypeBuilder builderToCallAlleles = new GenotypeBuilder(noCallGT);
            //TODO: update to support DRAGEN posteriors
            GATKVariantContextUtils.makeGenotypeCall(noCallGT.getPloidy(), builderToCallAlleles, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN,
                    noCallGT.getLikelihoods().getAsVector(), variant.getAlleles(), noCallGT, null);
            return builderToCallAlleles.make();
        } else {
            return variant.getGenotype(0);
        }
    }

    /**
     * Get the list of concrete alternate alleles in variant that were not called in genotype
     * @param variant   has full set of alleles, not all of which may be called
     * @param calledGenotype    should not be no-call
     * @return  a list of (concrete) alt alleles that are not called in calledGenotype
     */
    private List<Allele> getAllelesToDrop(final VariantContext variant, final Genotype calledGenotype) {
        //always drop alleles that aren't called to reduce PL size
        final List<Allele> allelesToDrop = variant.getAlternateAlleles().stream()
                .filter(a -> !a.isSymbolic() && ! calledGenotype.getAlleles().contains(a))
                .collect(Collectors.toList());

        //if the position of variant overlaps the ref block buffer, then it means that its deletion has been converted to a ref block
        //(if there was a subsequent high quality deletion, it would have trimmed the buffer)
        if (calledGenotype.getAlleles().contains(Allele.SPAN_DEL) && (vcfWriter.siteOverlapsBuffer(variant)
                || vcfWriter.getVcfOutputEnd() == null
                || vcfWriter.getVcfOutputEnd().getEnd() < variant.getStart())) {
            allelesToDrop.add(Allele.SPAN_DEL);
        }
        return allelesToDrop;
    }

    /**
     * Add the "raw" annotations necessary for calculating QD and AS_QD
     * @param destination   has qual-related annotations added to it, but also potentially supplies DP value
     * @param posteriorsKey if not null, the key for genotype posteriors in this VCF version
     * @param updatedAllelesVC  variant context without uncalled alts
     */
    private static void addQualAnnotations(final Map<String, Object> destination, final String posteriorsKey, final VariantAnnotatorEngine annotationEngine, final VariantContext updatedAllelesVC) {
        final Genotype updatedAllelesGenotype = updatedAllelesVC.getGenotype(0);
        destination.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, updatedAllelesGenotype.getPL()[0]);
        int varDP = QualByDepth.getDepth(updatedAllelesVC.getGenotypes(), null);
        if (varDP == 0) {  //prevent QD=Infinity case
            //attrMap should have DP already from copyInfoAnnotations call above
            varDP = Integer.parseInt(destination.getOrDefault(VCFConstants.DEPTH_KEY, 1).toString()); //if there's no VarDP and no DP, just prevent Infs/NaNs and QD will get capped later
        }
        destination.put(GATKVCFConstants.VARIANT_DEPTH_KEY, varDP);
        if (annotationEngine.hasInfoAnnotation("AS_QualByDepth")) {
            final List<String> quals = new ArrayList<>();
            //get allele-specific QUAL approximation by subsetting PLs for each alt
            for (final Allele alt : updatedAllelesVC.getAlternateAlleles()) {
                if (alt.equals(Allele.NON_REF_ALLELE) || alt.equals(Allele.SPAN_DEL)) {
                    quals.add("0");
                    continue;
                }
                //TODO: this isn't going to work for DRAGEN's genotype posteriors
                final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(updatedAllelesVC.getGenotypes(),
                        updatedAllelesGenotype.getPloidy(), updatedAllelesVC.getAlleles(), Arrays.asList(updatedAllelesVC.getReference(), alt), null,
                        GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
                //assignment method doesn't really matter as long as we don't zero out PLs; don't need depth to get PLs for quals

                final Genotype subsettedGenotype = gc.get(0);
                final int[] likelihoods = getGenotypePosteriorsOtherwiseLikelihoods(subsettedGenotype, posteriorsKey);
                if (likelihoods != null) {
                    quals.add(Integer.toString(likelihoods[0]));
                }  else {  //AlleleSubsettingUtils can no-call genotypes with super duper low GQs
                    quals.add("0");
                }
            }
            destination.put(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM + String.join(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, quals));
            final List<Integer> as_varDP = AS_QualByDepth.getAlleleDepths(updatedAllelesVC.getGenotypes());
            if (as_varDP != null) {
                destination.put(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, as_varDP.stream().map(n -> Integer.toString(n)).collect(Collectors.joining(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM)));
            }
        }
    }

    /**
     * Any AD counts for the non-ref allele get propagated to every new allele when GVCFs are merged, so zero them out
     * @param g a genotype that may or may not contain AD
     * @param nonRefInd allele index of the non-ref, -1 if missing
     * @return  a Genotype with the non-ref AD zeroed out and new DP, or the original genotype if it doesn't have AD
     */
    private Genotype removeNonRefADs(final Genotype g, final int nonRefInd) {
        if (g.hasAD() && nonRefInd != -1) {
            final int[] ad = g.getAD();
            if (ad.length >= nonRefInd && ad[nonRefInd] > 0) { //only initialize a builder if we have to
                final GenotypeBuilder gb = new GenotypeBuilder(g);
                final int nonRefAdCount = ad[nonRefInd];
                ad[nonRefInd] = 0;
                gb.AD(ad);
                if (g.hasDP()) {
                    gb.DP(g.getDP() - nonRefAdCount);
                } else {
                    gb.DP((int) MathUtils.sum(ad));
                }
                return gb.make();
            } else {
                return g;
            }
        } else {
            return g;
        }
    }

    /**
     * Add the original annotations to the map for the new VC, subsetting AS annotations as necessary
     * Note that it is expected that both variant contexts contain a NON_REF allele
     * @param destinationAttrMap   map of new annotations, to be modified
     * @param sourceVC    VC with full set of alleles and INFO annotations
     * @param infoFieldAnnotationKeyNamesToRemove   annotations to exclude from the final output
     * @param annotationEngine  used to get the list of annotations
     * @param allelesNeedSubsetting do we have to subset allele-specific annotations?
     * @param relevantIndices   indexes for called alleles within the full alleles set from original VC
     */
    private static void copyInfoAnnotations(final Map<String, Object> destinationAttrMap, final VariantContext sourceVC,
                                            final List<String> infoFieldAnnotationKeyNamesToRemove, final VariantAnnotatorEngine annotationEngine,
                                            final boolean allelesNeedSubsetting, final int[] relevantIndices) {
        //copy over info annotations
        final Map<String, Object> origMap = sourceVC.getAttributes();
        final List<VariantAnnotation> engineAnnotations = new ArrayList<>();
        engineAnnotations.addAll(annotationEngine.getInfoAnnotations());
        engineAnnotations.addAll(annotationEngine.getJumboInfoAnnotations());
        for(final VariantAnnotation annotation : engineAnnotations) {
            for (final String key : annotation.getKeyNames()) {
                if (infoFieldAnnotationKeyNamesToRemove.contains(key)) {
                    continue;
                }
                if (origMap.containsKey(key)) {
                    destinationAttrMap.put(key, origMap.get(key));
                }
            }
            if (annotation instanceof AlleleSpecificAnnotation) {
                final boolean isReducible = annotation  instanceof ReducibleAnnotation;
                final List<String> keyNames = isReducible ? ((ReducibleAnnotation)annotation).getRawKeyNames() :
                        annotation.getKeyNames();
                for (final String rawKey : keyNames) {
                    if (infoFieldAnnotationKeyNamesToRemove.contains(rawKey)) {
                        continue;
                    }
                    if (origMap.containsKey(rawKey)) {
                        if (allelesNeedSubsetting && AnnotationUtils.isAlleleSpecific(annotation) && AnnotationUtils.isAlleleSpecificGatkKey(rawKey)) {
                            final List<String> alleleSpecificValues;
                            if (isReducible) {
                                alleleSpecificValues = AnnotationUtils.getAlleleLengthListOfStringFromRawData(sourceVC.getAttributeAsString(rawKey, null));
                            } else {
                                String value = sourceVC.getAttributeAsString(rawKey, null);
                                if (value == null) {
                                    alleleSpecificValues = Collections.emptyList();
                                } else {
                                    alleleSpecificValues = AnnotationUtils.getAlleleLengthListOfString(value);
                                }
                            }
                            final List<String> subsetList;
                            if (alleleSpecificValues.size() > 0) {
                                if (isReducible) {
                                    subsetList = AlleleSubsettingUtils.remapRLengthList(alleleSpecificValues, relevantIndices, "");
                                } else {
                                    subsetList = AlleleSubsettingUtils.remapALengthList(alleleSpecificValues, relevantIndices, "");
                                }
                                if (sourceVC.getAlleles().get(relevantIndices[relevantIndices.length - 1]).equals(Allele.NON_REF_ALLELE)) {
                                    //zero out non-ref value, just in case
                                    subsetList.set(subsetList.size() - 1, ((AlleleSpecificAnnotation) annotation).getEmptyRawValue());
                                }
                            } else {
                                subsetList = Collections.nCopies(relevantIndices.length, "");
                            }

                            if (isReducible) {
                                destinationAttrMap.put(rawKey, AnnotationUtils.encodeAnyASListWithRawDelim(subsetList));
                            } else {
                                destinationAttrMap.put(rawKey, AnnotationUtils.encodeStringList(subsetList));
                            }
                        } else {
                            destinationAttrMap.put(rawKey, origMap.get(rawKey));
                        }
                    }
                }
            }
        }
    }

    /**
     * Add the newest raw mapping quality annotations to the annotation map
     * @param destinationAttrMap   modified to add mapping quality raw annotations
     * @param sourceVC    VariantContext that may contain deprecated mapping quality annotations
     */
    private static void updateMQAnnotations(final Map<String, Object> destinationAttrMap, final VariantContext sourceVC) {
        //all VCs should get new RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, but preserve deprecated keys if present
        if (!sourceVC.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY)) {
            //we're going to approximate depth for MQ calculation with the site-level DP (should be informative and uninformative reads),
            //which is pretty safe because it will only differ if reads are missing MQ
            final int rawMqValue = sourceVC.hasAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED) ?
                    (int)Math.round(sourceVC.getAttributeAsDouble(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0.0)) :
                    (int)Math.round(sourceVC.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, 60.0) *
                            sourceVC.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, 60.0) *
                            sourceVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            destinationAttrMap.put(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,
                    rawMqValue + "," + sourceVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));

            //NOTE: this annotation is deprecated, but keep it here so we don't have to reprocess older GVCFs
            if (sourceVC.hasAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
                destinationAttrMap.put(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED,
                        sourceVC.getAttributeAsDouble(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0));
                destinationAttrMap.put(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED,
                        sourceVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            }
        }
    }

    /**
     * If genotype posterior probabilities are present, return those; otherwise use likelihoods
     * @param genotype  should have PLs
     * @param posteriorsKey if not null, the key for genotype posteriors in this VCF version
     * @return may be null
     */
    private static int[] getGenotypePosteriorsOtherwiseLikelihoods(final Genotype genotype, final String posteriorsKey) {
        if ((posteriorsKey != null && genotype.hasExtendedAttribute(posteriorsKey))) {
            final double[] posteriors = VariantContextGetters.getAttributeAsDoubleArray(genotype, posteriorsKey, () -> null, 0);
            return Arrays.stream(posteriors).mapToInt(x -> (int)Math.round(x)).toArray();
        } else {
            return genotype.getPL();
        }
    }

    /**
     * Given a variant with multi-allelic PLs, modify the GenotypeBuilder to have annotations as for just ref and non-ref
     * @param result    a VariantContext containing alternate alleles in addition to non-ref and a genotype that should be hom-ref
     * @param gb    a reference block GenotypeBuilder to be modified
     */
    private void subsetHomRefPosteriorsToRefVersusNonRef(final VariantContext result, final GenotypeBuilder gb) {
        //TODO: bestAlleles needs to be modified for posteriors
        final Genotype genotype = result.getGenotype(0);
        final List<Allele> bestAlleles = AlleleSubsettingUtils.calculateMostLikelyAlleles(result, genotype.getPloidy(), 1, false);
        final Allele bestAlt = bestAlleles.stream().filter(a -> !a.isReference()).findFirst().orElse(Allele.NON_REF_ALLELE);  //allow span dels
        final int[] idxVector = result.getGLIndicesOfAlternateAllele(bestAlt);
        final int[] multiallelicPLs = getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey);
        if (multiallelicPLs != null) {
            int[] newPLs = new int[GenotypeLikelihoods.numLikelihoods(2, genotype.getPloidy())];
            for (int i = 0; i < idxVector.length; i++) {
                newPLs[i] = multiallelicPLs[idxVector[i]];
            }
            //in the case of *, we need to renormalize to homref
            if (newPLs[0] != 0) {
                final int[] output = new int[newPLs.length];
                for (int i = 0; i < newPLs.length; i++) {
                    output[i] = Math.max(newPLs[i] - newPLs[0], 0);
                }
                newPLs = output;
            }
            gb.PL(newPLs);
            gb.GQ(MathUtils.secondSmallestMinusSmallest(newPLs, 0));
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
