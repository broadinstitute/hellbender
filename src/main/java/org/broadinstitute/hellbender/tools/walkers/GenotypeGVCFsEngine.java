package org.broadinstitute.hellbender.tools.walkers;


import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Engine class to allow for other classes to replicate the behavior of GenotypeGVCFs. See {@link GenotypeGVCFs} for details
 *
 * Usage:
 * -Pass the genotype args into the constructor, which will the initialize the engine completely
 * Get the appropriate writer and write the appropriate header via {@link #setupVCFWriter}
 * Repeatedly call {@link #callRegion} to call variants in each region, and add them to your writer
 */

public class GenotypeGVCFsEngine
{
    private static final Logger logger = LogManager.getLogger(GenotypeGVCFsEngine.class);
    private static final String GVCF_BLOCK = "GVCFBlock";

    //the annotation engine
    private VariantAnnotatorEngine annotationEngine = null;

    //the genotyping engine
    private GenotypingEngine<?> forceOutputGenotypingEngine = null;
    private MinimalGenotypingEngine genotypingEngine = null;

    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();

    private GenotypeCalculationArgumentCollection genotypeArgs;

    // INFO Header names that require alt alleles
    final LinkedHashSet<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();

    private boolean includeNonVariants;

    private VCFHeader outputHeader;

    private SampleList samples;

    private DragstrParams dragStrParams;

    final VCFHeader inputVCFHeader;

    final boolean   keepSB;

    /**
     * Create and initialize a new GenotypeGVCFsEngine given a collection of GenotypeGVCF arguments and a VCF header
     *  @param annotationEngine variantAnnotatorEngine with annotations to process already added
     * @param genotypeArgs command-line arguments for the GenotypeGVCFs caller
     * @param includeNonVariants true to save INFO header names that require alt alleles
     * @param inputVCFHeader header for the VCF
     * @param keepSB keep SB attribute (STRAND_BIAS_BY_SAMPLE)
     */
    public GenotypeGVCFsEngine(final VariantAnnotatorEngine annotationEngine, final GenotypeCalculationArgumentCollection genotypeArgs,
                               final boolean includeNonVariants, final VCFHeader inputVCFHeader, final boolean keepSB)
    {
        this.annotationEngine = annotationEngine;
        this.genotypeArgs = genotypeArgs;
        this.includeNonVariants = includeNonVariants;
        this.inputVCFHeader = inputVCFHeader;
        this.keepSB = keepSB;
        initialize();
    }

    private void initialize()
    {
        samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples()); //todo should this be getSampleNamesInOrder?

        // Request INFO field annotations inheriting from RankSumTest and RMSAnnotation added to remove list
        for ( final InfoFieldAnnotation annotation :  annotationEngine.getInfoAnnotations() ) {
            if ( annotation instanceof RankSumTest ||
                    annotation instanceof AS_RMSMappingQuality ||
                    annotation instanceof RMSMappingQuality) {
                final List<String> keyNames = annotation.getKeyNames();
                if ( !keyNames.isEmpty() ) {
                    infoFieldAnnotationKeyNamesToRemove.add(keyNames.get(0));
                }
            }
        }

        // We only want the engine to generate the AS_QUAL key if we are using AlleleSpecific annotations.
        genotypingEngine = new MinimalGenotypingEngine(createMinimalArgs(false), samples,
                annotationEngine.getInfoAnnotations().stream().anyMatch(AnnotationUtils::isAlleleSpecific), this.dragStrParams);
        forceOutputGenotypingEngine = new MinimalGenotypingEngine(createMinimalArgs(true), samples,
                annotationEngine.getInfoAnnotations().stream().anyMatch(AnnotationUtils::isAlleleSpecific));

        if ( includeNonVariants ) {
            // Save INFO header names that require alt alleles
            for ( final VCFHeaderLine headerLine : inputVCFHeader.getMetaDataInInputOrder() ) {
                if (headerLine instanceof VCFInfoHeaderLine ) {
                    if (((VCFInfoHeaderLine) headerLine).getCountType() == VCFHeaderLineCount.A) {
                        infoHeaderAltAllelesLineNames.add(((VCFInfoHeaderLine) headerLine).getID());
                    }
                }
            }
        }
    }

    public VariantContext callRegion(Locatable loc, List<VariantContext> variants, ReferenceContext ref, FeatureContext features,
                                     ReferenceConfidenceVariantContextMerger merger, boolean somaticInput, double tlodThreshold,
                                     double afTolerance, final boolean outputNonVariants) //do work for apply
    {
        final List<VariantContext> variantsToProcess = getVariantSubsetToProcess(loc, variants);

        if (dragStrParams == null || genotypeArgs.dontUseDragstrPriors) {
            ref.setWindow(10, 10); //TODO this matches the gatk3 behavior but may be unnecessary
        } else {
            ref.setWindow(dragStrParams.maximumLengthInBasePairs(), dragStrParams.maximumLengthInBasePairs());
        }
        genotypingEngine.setReferenceContext(ref);
        final VariantContext mergedVC = merger.merge(variantsToProcess, loc, ref.getBase(), true, false);
        final VariantContext regenotypedVC = somaticInput ? regenotypeSomaticVC(mergedVC, ref, features, outputNonVariants, tlodThreshold, afTolerance) :
                regenotypeVC(mergedVC, ref, features, outputNonVariants);

        return regenotypedVC;
    }


    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    private VariantContext regenotypeVC(final VariantContext originalVC, final ReferenceContext ref, final FeatureContext features, boolean includeNonVariants) {
        Utils.nonNull(originalVC);

        final VariantContext result;

        // only re-genotype polymorphic sites
        if ( originalVC.isVariant()  && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
            // note that the calculateGenotypes method also calculates the QUAL score
            final VariantContext regenotypedVC = calculateGenotypes(originalVC, includeNonVariants);
            if (regenotypedVC == null) {
                return null;
            }
            if (GATKVariantContextUtils.isProperlyPolymorphic(regenotypedVC) || includeNonVariants) {
                // Note that reversetrimAlleles must be performed after the annotations are finalized because the reducible annotation data maps
                // were generated and keyed on the un reverseTrimmed alleles from the starting VariantContexts. Thus reversing the order will make
                // it difficult to recover the data mapping due to the keyed alleles no longer being present in the variant context.
                final VariantContext withGenotypingAnnotations = addGenotypingAnnotations(originalVC.getAttributes(), regenotypedVC);
                final VariantContext withAnnotations = annotationEngine.finalizeAnnotations(withGenotypingAnnotations, originalVC);
                final int[] relevantIndices = regenotypedVC.getAlleles().stream().mapToInt(a -> originalVC.getAlleles().indexOf(a)).toArray();
                final VariantContext trimmed = GATKVariantContextUtils.reverseTrimAlleles(withAnnotations);
                final GenotypesContext updatedGTs = subsetAlleleSpecificFormatFields(outputHeader, trimmed.getGenotypes(), relevantIndices);
                result = new VariantContextBuilder(trimmed).genotypes(updatedGTs).make();
            } else {
                return null;
            }
        } else {
            result = originalVC;
        }


        // if it turned monomorphic then we either need to ignore or fix such sites
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        // Note that we also check depth because the previous conditions were also contingent on depth
        if (result.isPolymorphicInSamples()  && result.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0) {
            // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
            final VariantContextBuilder vcBuilder = new VariantContextBuilder(result);
            //don't count sites with no depth and no confidence towards things like AN and InbreedingCoeff
            vcBuilder.genotypes(assignNoCallsAnnotationExcludedGenotypes(result.getGenotypes()));
            VariantContext annotated = annotationEngine.annotateContext(vcBuilder.make(), features, ref, null, a -> true);
            return new VariantContextBuilder(annotated).genotypes(cleanupGenotypeAnnotations(annotated, false, keepSB)).make();
        } else if (includeNonVariants) {
            // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
            VariantContext preannotated = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true, false)).make();
            return annotationEngine.annotateContext(preannotated, features, ref, null, GenotypeGVCFsEngine::annotationShouldBeSkippedForHomRefSites);
        } else {
            return null;
        }
    }

    private static boolean annotationShouldBeSkippedForHomRefSites(VariantAnnotation annotation) {
        return annotation instanceof RankSumTest || annotation instanceof RMSMappingQuality || annotation instanceof AS_RMSMappingQuality;
    }

    private GenotypesContext subsetAlleleSpecificFormatFields(final VCFHeader outputHeader, final GenotypesContext originalGs, final int[] relevantIndices) {
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());
        for (final Genotype g : originalGs) {
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            final Set<String> keys = g.getExtendedAttributes().keySet();
            for (final String key : keys) {
                final VCFFormatHeaderLine headerLine = outputHeader.getFormatHeaderLine(key);
                final Object attribute;
                if (headerLine.getCountType().equals(VCFHeaderLineCount.INTEGER) && headerLine.getCount() == 1) {
                    attribute = g.getAnyAttribute(key);
                }
                else {
                    attribute = ReferenceConfidenceVariantContextMerger.generateAnnotationValueVector(headerLine.getCountType(),
                            VariantContextGetters.attributeToList(g.getAnyAttribute(key)), relevantIndices);
                }
                gb.attribute(key, attribute);
            }
            newGTs.add(gb.make());
        }
        return newGTs;
    }

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations and any other annotations generated by the genotyper.
        final Map<String, Object> attrs = new LinkedHashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY)) {
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
        }
        if (newVC.hasAttribute(GATKVCFConstants.AS_QUAL_KEY)) {
            attrs.put(GATKVCFConstants.AS_QUAL_KEY, newVC.getAttribute(GATKVCFConstants.AS_QUAL_KEY));
        }
        return new VariantContextBuilder(newVC).attributes(attrs).make();
    }

    private VariantContext calculateGenotypes(VariantContext vc, final boolean forceOutput) {
        return (forceOutput ? forceOutputGenotypingEngine : genotypingEngine).calculateGenotypes(vc, null, Collections.emptyList());
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    private VariantContext regenotypeSomaticVC(final VariantContext originalVC, final ReferenceContext ref, final FeatureContext features, boolean includeNonVariants, double tlodThreshold, double afTolerance) {
        Utils.nonNull(originalVC);

        final VariantContext result;
        if ( originalVC.isVariant()  && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
            result = callSomaticGenotypes(originalVC, tlodThreshold, afTolerance);
        } else if (includeNonVariants) {
            result = originalVC;
        } else {
            result = null;
        }
        return result;
    }

    /**
     * Drop low quality alleles and call genotypes
     * CombineGVCFs will convert calls to no-call (of varying ploidy, as is the case in somatic)
     *
     * @param vc input VariantContext with no-called genotypes
     * @return a VC with called genotypes and low quality alleles removed, may be null
     */
    private VariantContext callSomaticGenotypes(final VariantContext vc, double tlodThreshold, double afTolerance) {
        final List<Genotype> newGenotypes = new ArrayList<>();
        final GenotypesContext genotypes = vc.getGenotypes();
        final double[] perAlleleLikelihoodSums = new double[vc.getAlleles().size()];  //needs the ref for the subsetting utils

        for(final Genotype g : genotypes) {
            GenotypeBuilder gb = new GenotypeBuilder(g);
            final double[] tlodArray = VariantContextGetters.getAttributeAsDoubleArray(g, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
            final double[] variantAFArray = VariantContextGetters.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 0.0);
            double variantAFtotal = 0;
            final List<Allele> calledAlleles = new ArrayList<>();
            for(int i = 0; i < vc.getAlleles().size()-1; i++) {
                variantAFtotal += variantAFArray[i];
                if (tlodArray[i] > tlodThreshold) {
                    calledAlleles.add(vc.getAlternateAllele(i));
                    perAlleleLikelihoodSums[i+1] += tlodArray[i];
                }
            }
            //hack for weird Mutect2 ploidy -- if the variant is non-homoplasmic, call the reference allele too
            if(variantAFtotal < 1-afTolerance && (!g.hasAD() || g.getAD()[0] > 0)) {

                calledAlleles.add(0, vc.getReference());
            }
            //"ploidy" gets set according to the size of the alleles List in the Genotype
            gb.alleles(calledAlleles);
            newGenotypes.add(gb.make());
        }

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final VariantContext regenotypedVC = builder.genotypes(newGenotypes).make();

        final int maxAltAlleles = genotypingEngine.getConfiguration().genotypeArgs.maxAlternateAlleles;
        List<Allele> allelesToKeep;

        //we need to make sure all alleles pass the tlodThreshold
        allelesToKeep = new ArrayList<>(perAlleleLikelihoodSums.length-1);
        allelesToKeep.add(vc.getReference());
        for (int i = 1; i < perAlleleLikelihoodSums.length; i++) {
            if (perAlleleLikelihoodSums[i] > tlodThreshold) {
                allelesToKeep.add(vc.getAlternateAllele(i-1));
            }
        }

        if (regenotypedVC.getAlternateAlleles().size() > maxAltAlleles) {
            allelesToKeep = AlleleSubsettingUtils.filterToMaxNumberOfAltAllelesBasedOnScores(maxAltAlleles, allelesToKeep, perAlleleLikelihoodSums);
        }

        if (allelesToKeep.size() == 1) {
            return null;
        }

        //if we didn't drop alleles then we're done!
        if (allelesToKeep.size() == regenotypedVC.getAlleles().size()) {
            return regenotypedVC;
        }

        final int[] relevantIndices = allelesToKeep.stream().mapToInt(a -> regenotypedVC.getAlleles().indexOf(a)).toArray();

        //do another pass over genotypes to drop the alleles that aren't called
        final GenotypesContext reducedGenotypes = AlleleSubsettingUtils.subsetSomaticAlleles(outputHeader, regenotypedVC.getGenotypes(), allelesToKeep, relevantIndices);
        final VariantContext subsetVC = builder.alleles(allelesToKeep).genotypes(reducedGenotypes).make();
        final VariantContext trimmedVC = GATKVariantContextUtils.trimAlleles(subsetVC, true, true);
        if (GATKVariantContextUtils.isProperlyPolymorphic(trimmedVC)) {
            return trimmedVC;
        }
        else {
            return null;
        }
    }

    // If includeNonVariants is set, we're using group-by-locus traversal. To match GATK3 GenotypeGVCFs,
    // see if there is a variant in the overlapping group that starts exactly at the locus start position, and if so
    // prioritize and process only that variant. Otherwise process all of the overlapping variants.
    private List<VariantContext> getVariantSubsetToProcess(final Locatable loc, List<VariantContext> preProcessedVariants) {
        if (includeNonVariants) {
            final List<VariantContext> matchingStart =
                    preProcessedVariants.stream().filter(vc -> vc.getStart() == loc.getStart()).collect(Collectors.toList());
            if (matchingStart.size() == 0) {
                return preProcessedVariants;
            }
            else if (matchingStart.size() == 1) {
                return matchingStart;
            }
            // since this tool only accepts a single input source, there should never be
            // more than one variant at a given starting locus
            throw new IllegalStateException(
                    String.format(
                            "Variant input contains more than one variant starting at location: %s",
                            new SimpleInterval(matchingStart.get(0))));
        } else {
            return preProcessedVariants;
        }
    }

    /**
     * Creates a StandardCallerArgumentCollection with appropriate values filled in from the arguments in this walker
     */
    private StandardCallerArgumentCollection createMinimalArgs(final boolean forceOutput) {
        final StandardCallerArgumentCollection args = new StandardCallerArgumentCollection();
        args.genotypeArgs = genotypeArgs.clone();

        //keep hom ref calls even if no PLs
        args.genotypeArgs.genotypeAssignmentMethod = GenotypeAssignmentMethod.PREFER_PLS;

        //whether to emit non-variant sites is not contained in genotypeArgs and must be passed to args separately
        //Note: GATK3 uses OutputMode.EMIT_ALL_CONFIDENT_SITES when includeNonVariants is requested
        //GATK4 uses EMIT_ALL_ACTIVE_SITES to ensure LowQual sites are emitted.
        args.outputMode = forceOutput ? OutputMode.EMIT_ALL_ACTIVE_SITES : OutputMode.EMIT_VARIANTS_ONLY;
        return args;
    }

    /**
     * Create a VCF header in the writer
     *
     * @param vcfWriter
     * @return a VCF writer

     */
    public VariantContextWriter setupVCFWriter(Set<VCFHeaderLine> defaultToolVCFHeaderLines, boolean keepCombined, DbsnpArgumentCollection dbsnp, VariantContextWriter vcfWriter) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(defaultToolVCFHeaderLines);

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if (keepCombined) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_QUAL_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        }
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        outputHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        vcfWriter.writeHeader(outputHeader);

        return vcfWriter;
    }


    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     * 5. move GQ to RGQ if the site is monomorphic
     *
     * @param vc            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @param keepSB        keep value of SB attribute
     * @return a new set of Genotypes
     */
    @VisibleForTesting
    static List<Genotype> cleanupGenotypeAnnotations(final VariantContext vc, final boolean createRefGTs, final boolean keepSB) {
        final GenotypesContext oldGTs = vc.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) ) {
                depth = parseInt(oldGT.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY));
                builder.DP(depth);
                attrs.remove(GATKVCFConstants.MIN_DP_FORMAT_KEY);
            }

            if ( !keepSB ) {
                attrs.remove(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
            }

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, GenotypeGVCFs.PHASED_HOM_VAR_STRING);
            }

            // create AD if it's not there, but only if there's data
            if ( !oldGT.hasAD() && vc.isVariant() && depth > 0) {
                final int[] AD = new int[vc.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                 //keep 0 depth samples and 0 GQ samples as no-call
                if (depth > 0 && oldGT.hasGQ()) {
                    if (oldGT.getGQ() > 0) {
                        final List<Allele> refAlleles = Collections.nCopies(oldGT.getPloidy(), vc.getReference());
                        builder.alleles(refAlleles);
                    } else {
                        builder.alleles(Collections.nCopies(oldGT.getPloidy(),Allele.NO_CALL));
                    }

                    // move the GQ to RGQ
                    builder.noGQ();
                    attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
                } else {
                    builder.alleles(Collections.nCopies(oldGT.getPloidy(),Allele.NO_CALL));
                    builder.noGQ().noDP();
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    /**
     *  Does this genotype look like it has no reads and should be excluded from annotations?
     */
    static boolean excludeFromAnnotations(Genotype oldGT) {
        return (oldGT.isHomRef() || oldGT.isNoCall())
                && (!oldGT.hasDP() || oldGT.getDP() == 0)
                && oldGT.hasGQ() && oldGT.getGQ() == 0;
    }

    /**
     *  for GQ0 hom-refs that have no data, turn them back into no-calls so they don't affect the INFO annotations
     */
    private List<Genotype> assignNoCallsAnnotationExcludedGenotypes(final GenotypesContext genotypes) {
        final List<Genotype> returnList = new ArrayList<>();
        for (final Genotype oldGT : genotypes) {
            //convert GQ0s that were reblocked back to no-calls for better AN and InbreedingCoeff annotations
            if (excludeFromAnnotations(oldGT)) {
                final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
                builder.alleles(Collections.nCopies(oldGT.getPloidy(), Allele.NO_CALL));
                returnList.add(builder.make());
            } else {
                returnList.add(oldGT);
            }

        }
        return returnList;
    }


    private static int parseInt(Object attribute){
        if( attribute instanceof String) {
            return Integer.parseInt((String)attribute);
        } else if ( attribute instanceof Number){
            return ((Number) attribute).intValue();
        } else {
            throw new IllegalArgumentException("Expected a Number or a String but found something else.");
        }
    }
}
