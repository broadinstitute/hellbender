package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Variant context utilities related to merging variant-context instances.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@SuppressWarnings({"rawtypes","unchecked"}) //TODO fix uses of untyped Comparable.
public final class ReferenceConfidenceVariantContextMerger {

    private static VCFHeader vcfInputHeader = null;
    protected final VariantAnnotatorEngine annotatorEngine;
    private final boolean doSomaticMerge;
    protected boolean dropSomaticFilteringAnnotations;
    protected boolean callGTAlleles;
    protected final OneShotLogger oneShotAnnotationLogger = new OneShotLogger(this.getClass());
    protected final OneShotLogger oneShotHeaderLineLogger = new OneShotLogger(this.getClass());
    protected final OneShotLogger AS_Warning = new OneShotLogger(this.getClass());
    List<String> SOMATIC_INFO_ANNOTATIONS_TO_MOVE = Arrays.asList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY);

    private static final List<String> SOMATIC_FORMAT_ANNOTATIONS_TO_KEEP = Arrays.asList(
            GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY,
            GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY,
            GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY,
            VCFConstants.PHASE_SET_KEY);
    private static final List<String> SOMATIC_INFO_ANNOTATIONS_TO_DROP = Arrays.asList(
            GATKVCFConstants.POPULATION_AF_KEY);

    public ReferenceConfidenceVariantContextMerger(VariantAnnotatorEngine engine, final VCFHeader inputHeader) {
        this(engine, inputHeader, false, false);
    }

    public ReferenceConfidenceVariantContextMerger(VariantAnnotatorEngine engine, final VCFHeader inputHeader, boolean somaticInput) {
        this(engine, inputHeader, somaticInput, false);
    }

    public ReferenceConfidenceVariantContextMerger(VariantAnnotatorEngine engine, final VCFHeader inputHeader, boolean somaticInput, boolean dropSomaticFilteringAnnotations) {
        this(engine, inputHeader, somaticInput, dropSomaticFilteringAnnotations, false);
    }

    public ReferenceConfidenceVariantContextMerger(VariantAnnotatorEngine engine, final VCFHeader inputHeader, boolean somaticInput, boolean dropSomaticFilteringAnnotations, boolean callGTAlleles) {
        Utils.nonNull(inputHeader, "A VCF header must be provided");

        annotatorEngine = engine;
        vcfInputHeader = inputHeader;
        doSomaticMerge = somaticInput;
        this.dropSomaticFilteringAnnotations = dropSomaticFilteringAnnotations;
        this.callGTAlleles = callGTAlleles;
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * This version applies allele remapping prior to genotyping.
     *
     * TODO: THIS BEHAVIOR SHOULD BE MADE THE DEFAULT, AND THIS METHOD RETIRED IN FAVOR OF MERGE() BELOW
     * TODO: SEE https://github.com/broadinstitute/gatk/issues/8317
     *
     * @param vcs     collection of unsorted genomic vcs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @param samplesAreUniquified  if true, sample names have been uniquified
     * @return new VariantContext representing the merge of all vcs or null if it not relevant
     */
    public VariantContext mergeWithRemapping(final List<VariantContext> vcs, final Locatable loc, final Byte refBase,
                                             final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified) {
        // Call merge() with useRemappedAllelesForGenotyping == true
        return merge(vcs, loc, refBase, removeNonRefSymbolicAllele, samplesAreUniquified, true);
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * This overload does not apply allele remapping prior to genotyping.
     *
     * @param vcs     collection of unsorted genomic vcs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @param samplesAreUniquified  if true, sample names have been uniquified
     * @return new VariantContext representing the merge of all vcs or null if it not relevant
     */
    public VariantContext merge(final List<VariantContext> vcs, final Locatable loc, final Byte refBase,
                                final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified) {
        // Call merge() with useRemappedAllelesForGenotyping == false
        return merge(vcs, loc, refBase, removeNonRefSymbolicAllele, samplesAreUniquified, false);
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * @param vcs     collection of unsorted genomic vcs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @param samplesAreUniquified  if true, sample names have been uniquified
     * @param useRemappedAllelesForGenotyping  if true, remap alleles prior to making the genotype call
     * @return new VariantContext representing the merge of all vcs or null if it not relevant
     */
    public VariantContext merge(final List<VariantContext> vcs, final Locatable loc, final Byte refBase,
                                final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified,
                                final boolean useRemappedAllelesForGenotyping) {
        Utils.nonEmpty(vcs);

        // establish the baseline info (sometimes from the first VC)
        final String name = vcs.get(0).getSource();

        // ref allele
        final Allele refAllele = determineReferenceAlleleGivenReferenceBase(vcs, loc, refBase);
        if ( refAllele == null ) {
            return null;
        }

        // In this list we hold the mapping of each variant context alleles.
        final List<VCWithNewAlleles> vcAndNewAllelePairs = new ArrayList<>(vcs.size());
        final Set<String> aggregatedFilters = new HashSet<>();
        boolean sawPassSample = false;

        // cycle through and add info from the other vcs
        for ( final VariantContext vc : vcs ) {
            // if this context doesn't start at the current location then it must be a spanning event (deletion or ref block)
            final boolean isSpanningEvent = loc.getStart() != vc.getStart();
            vcAndNewAllelePairs.add(new VCWithNewAlleles(vc, isSpanningEvent ? replaceWithNoCallsAndDels(vc, doSomaticMerge) : remapAlleles(vc, refAllele), isSpanningEvent));
            if (doSomaticMerge && vc.filtersWereApplied()) {
                if (vc.isFiltered()) {
                    aggregatedFilters.addAll(vc.getFilters());
                }
                else {
                    sawPassSample = true;
                }
            }
        }

        final List<Allele> allelesList = collectTargetAlleles(vcAndNewAllelePairs, refAllele, removeNonRefSymbolicAllele);

        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id
        int depth = 0;
        final Map<String, List<?>> annotationMap = new LinkedHashMap<>();

        final GenotypesContext genotypes = GenotypesContext.create();

        for ( final VCWithNewAlleles vcWithNewAlleles : vcAndNewAllelePairs ) {
            final VariantContext vc = vcWithNewAlleles.getVc();
            final List<Allele> remappedAlleles = vcWithNewAlleles.getNewAlleles();

            genotypes.addAll(mergeRefConfidenceGenotypes(vc, remappedAlleles, allelesList, samplesAreUniquified, useRemappedAllelesForGenotyping));
            depth += calculateVCDepth(vc);

            if ( loc.getStart() != vc.getStart() ) {
                continue;
            }

            // special case ID (just preserve it)
            if ( vc.hasID() ) {
                rsIDs.add(vc.getID());
            }

            // add attributes
            addReferenceConfidenceAttributes(vcWithNewAlleles, annotationMap);
        }

        final Map<String, Object> attributes = mergeAttributes(depth, allelesList, annotationMap);

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : String.join(",", rsIDs);

        // note that in order to calculate the end position, we need a list of alleles that doesn't include anything symbolic
        final VariantContextBuilder builder = new VariantContextBuilder()
                .source(name)
                .id(ID)
                .alleles(allelesList)
                .chr(loc.getContig())
                .start(loc.getStart())
                .computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
                .genotypes(genotypes).unfiltered()
                .attributes(new TreeMap<>(attributes)).log10PError(CommonInfo.NO_LOG10_PERROR);  // we will need to re-genotype later
        if (doSomaticMerge) {
            //if all samples are filtered, this will apply all those filters to the VCF
            if (aggregatedFilters.isEmpty() || sawPassSample) {
                    builder.passFilters();
            } else {
                builder.filters(aggregatedFilters);
            }
        }
        return builder.make();
    }

    /**
     * Replaces any alleles in the VariantContext with NO CALLS or the symbolic deletion allele as appropriate, except for the generic ALT allele
     *
     * @param vc   VariantContext with the alleles to replace
     * @param doSomaticMerge    if we're merging somatic data, don't use spanning deletions
     * @return non-null list of alleles
     */
    private static List<Allele> replaceWithNoCallsAndDels(final VariantContext vc, final boolean doSomaticMerge) {
        Utils.nonNull(vc);

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());

        // no-call the reference allele
        result.add(Allele.NO_CALL);

        // handle the alternate alleles
        for ( final Allele allele : vc.getAlternateAlleles() ) {
            final Allele replacement;
            if ( allele.equals(Allele.NON_REF_ALLELE) ) {
                replacement = allele;
            } else if ( !doSomaticMerge && allele.length() < vc.getReference().length() ) {
                replacement = Allele.SPAN_DEL;
            } else {
                replacement = Allele.NO_CALL;
            }

            result.add(replacement);
        }
        return result;
    }

    /**
     * This method does a couple of things:
     * <ul><li>
     *     remaps the vc alleles considering the differences between the final reference allele and its own reference,</li>
     * <li>
     *     collects alternative alleles present in variant context and add them to the {@code finalAlleles} set.
     * </li></ul>
     *
     * @param vc           the variant context.
     * @param refAllele    final reference allele.
     * @return never {@code null}
     */
    //TODO as part of a larger refactoring effort {@link #remapAlleles} can be merged with {@link GATKVariantContextUtils#remapAlleles}.
    public static List<Allele> remapAlleles(final VariantContext vc, final Allele refAllele) {
        final Allele vcRef = vc.getReference();
        final byte[] refBases = refAllele.getBases();
        final int extraBaseCount = refBases.length - vcRef.getBases().length;
        if (extraBaseCount < 0) {
            throw new IllegalStateException("the wrong reference was selected");
        }

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());
        result.add(refAllele);

        for (final Allele a : vc.getAlternateAlleles()) {
            if (a.isSymbolic()) {
                result.add(a);
            } else if ( a == Allele.SPAN_DEL ) {
                // add SPAN_DEL directly so we don't try to extend the bases
                result.add(a);
            } else if (a.isCalled()) {
                result.add(extendAllele(a, extraBaseCount, refBases));
            } else { // NO_CALL and strange miscellanea
                result.add(a);
            }
        }
        return result;
    }

    @VisibleForTesting
    protected static class VCWithNewAlleles {
        private final VariantContext vc;
        private final List<Allele> newAlleles;
        private final boolean isSpanningEvent;

        VCWithNewAlleles(final VariantContext vc, final List<Allele> newAlleles, final boolean isSpanningEvent) {
            this.vc = vc;
            this.newAlleles = newAlleles;
            this.isSpanningEvent = isSpanningEvent;
        }

        private Stream<Allele> filterAllelesForFinalSet() {
            return newAlleles.stream().filter(a -> !a.equals(Allele.NON_REF_ALLELE))
                    .filter(a -> !a.isReference())
                    .filter(a -> !(a.isSymbolic() && vc.isSymbolic())) // skip <*DEL> if there isn't a real alternate allele.
                    .filter(Allele::isCalled) ; // skip NO_CALL
        }

        // record whether it's also a spanning deletion/event (we know this because the VariantContext type is no
        // longer "symbolic" but "mixed" because there are real alleles mixed in with the symbolic non-ref allele)
        boolean isSpanningDeletion() {
            return  (isSpanningEvent && vc.isMixed())
                    || vc.getAlleles().contains(Allele.SPAN_DEL)
                    || vc.getAlleles().contains(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
        }

        boolean isNonSpanningEvent() {
            return !isSpanningEvent && vc.isMixed();

        }
        public VariantContext getVc() {
            return vc;
        }

        List<Allele> getNewAlleles() {
            return newAlleles;
        }
    }

    private static List<Allele> collectTargetAlleles(final List<VCWithNewAlleles> vcAndNewAllelePairs, final Allele refAllele, final boolean removeNonRefSymbolicAllele) {
        // FinalAlleleSet contains the alleles of the new resulting VC
        // Using linked set in order to guarantee a stable order
        final Set<Allele> finalAlleleSet = new LinkedHashSet<>(10);
        // Reference goes first
        finalAlleleSet.add(refAllele);

        vcAndNewAllelePairs.stream()
                .flatMap(VCWithNewAlleles::filterAllelesForFinalSet)
                .forEachOrdered(finalAlleleSet::add);

        final boolean sawSpanningDeletion = vcAndNewAllelePairs.stream().anyMatch(VCWithNewAlleles::isSpanningDeletion);
        final boolean sawNonSpanningEvent = vcAndNewAllelePairs.stream().anyMatch(VCWithNewAlleles::isNonSpanningEvent);

        // Add <DEL> and <NON_REF> to the end if at all required in in the output.
        if ( sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele) ) {
            finalAlleleSet.add(Allele.SPAN_DEL);
        }
        if (!removeNonRefSymbolicAllele) {
            finalAlleleSet.add(Allele.NON_REF_ALLELE);
        }

        return new ArrayList<>(finalAlleleSet);
    }

    /**
     * lookup the depth from the VC DP field or calculate by summing the depths of the genotypes
     */
    protected static int calculateVCDepth(VariantContext vc) {
        if ( vc.hasAttribute(VCFConstants.DEPTH_KEY) ) {
            return vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
        } else { // handle the gVCF case from the HaplotypeCaller
            return vc.getGenotypes().stream()
                    .mapToInt(ReferenceConfidenceVariantContextMerger::getBestDepthValue)
                    .sum();
        }
    }

    private Map<String, Object> mergeAttributes(int depth, List<Allele> alleleList, Map<String, List<?>> annotationMap) {
        final Map<String, Object> attributes = new LinkedHashMap<>();

        attributes.putAll(annotatorEngine.combineAnnotations(alleleList, annotationMap));

        // Afterwards, we simply add the median value for all the annotations that weren't recognized as reducible
        for (String key : annotationMap.keySet()) {
            List<Comparable<?> > values = (List<Comparable<?>>) annotationMap.get(key);
            if (values!= null && values.size() > 0 ) {
                final int size = values.size();
                if (size == 1) {
                    attributes.put(key, values.get(0));
                } else {
                    attributes.put(key, Utils.getMedianValue(values));
                }
            }
        }

        if ( depth > 0 ) {
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
        }

        // remove stale AC and AF based attributes
        removeStaleAttributesAfterMerge(attributes);
        return attributes;
    }

    @VisibleForTesting
    static int getBestDepthValue(final Genotype gt) {
        if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
            return Integer.parseInt(gt.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
        } else {
            return gt.hasDP() ? gt.getDP() : 0;
        }
    }

    /**
     * @param alleles  the original alleles alleles
     * @return a non-null alleles of non-symbolic alleles
     */
    private static List<Allele> nonSymbolicAlleles(final List<Allele> alleles) {
        return alleles.stream()
                .filter(a -> !a.isSymbolic())
                .collect(Collectors.toList());
    }

    /**
     * Determines the ref allele given the provided reference base at this position
     *
     * @param VCs     collection of unsorted genomic VariantContexts
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning
     * @return new Allele or null if no reference allele/base is available
     */
    private static Allele determineReferenceAlleleGivenReferenceBase(final List<VariantContext> VCs, final Locatable loc, final Byte refBase) {
        final Allele refAllele = GATKVariantContextUtils.determineReferenceAllele(VCs, loc);
        if ( refAllele == null ) {
            if (refBase == null) {
                return null;
            } else {
                return Allele.create(refBase, true);
            }
        } else {
            return refAllele;
        }
    }

    /**
     * Remove the stale attributes from the merged set
     *
     * @param attributes the attribute map
     */
    protected static void removeStaleAttributesAfterMerge(final Map<String, Object> attributes) {
        attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
        attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.END_KEY);
        attributes.remove(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY);
        attributes.remove(GATKVCFConstants.EVENT_COUNT_IN_REGION_KEY); //median doesn't make sense here so drop it; used for ClusteredEventFilter, which doesn't apply to MT
    }

    /**
     * Adds attributes to the global map from the new context in a sophisticated manner
     *
     * @param vcPair                     variant context with attributes to add from
     * @param annotationMap              map of annotations for combining later
     */
    @VisibleForTesting
    private void addReferenceConfidenceAttributes(final VCWithNewAlleles vcPair, final Map<String, List<?>> annotationMap) {
        for ( final Map.Entry<String, Object> p : vcPair.getVc().getAttributes().entrySet() ) {
            final String key = p.getKey();
            if (SOMATIC_INFO_ANNOTATIONS_TO_MOVE.contains(key) ||
                    SOMATIC_INFO_ANNOTATIONS_TO_DROP.contains(key) || Mutect2FilteringEngine.STANDARD_MUTECT_INFO_FIELDS_FOR_FILTERING.contains(key)) {
                continue;
            }

            // If the key corresponds to a requested reducible key, store the data as AlleleSpecificAnnotationData
            if (annotatorEngine.isRequestedReducibleRawKey(key)) {
                final List<Object> valueList = vcPair.getVc().getAttributeAsList(key);

                List<ReducibleAnnotationData<?>> values = (List<ReducibleAnnotationData<?>>) annotationMap.get(key);
                if (values == null) {
                    values = new ArrayList<>();
                    annotationMap.put(key, values);
                }
                StringJoiner joiner = new StringJoiner(",");
                for (Object s : valueList) joiner.add(s.toString());

                ReducibleAnnotationData<Object> pairData = new AlleleSpecificAnnotationData<>(vcPair.getNewAlleles(), joiner.toString());
                values.add(pairData);

            // Otherwise simply treat it as a number
            } else{
                final Object value = p.getValue();

                // add the annotation values to a list for combining later
                List<Comparable<?>> values = (List<Comparable<?>>) annotationMap.get(key);
                if (values == null) {
                    values = new ArrayList<>();
                    annotationMap.put(key, values);
                }
                try {
                    if (!value.toString().contains(",")) {
                        values.add(parseNumericInfoAttributeValue(vcfInputHeader, key, value.toString()));
                    }
                    else {
                        String[] valueArray = value.toString().split("\\[|" + AnnotationUtils.LIST_DELIMITER +  "|\\]");
                        for (final Object val : valueArray) {
                            if (!val.equals("")) {
                                values.add(parseNumericInfoAttributeValue(vcfInputHeader, key, val.toString()));
                            }
                        }
                    }
                } catch (final NumberFormatException e) {
                    oneShotAnnotationLogger.warn(String.format("Detected invalid annotations: When trying to merge variant contexts at location %s:%d the annotation %s was not a numerical value and was ignored",vcPair.getVc().getContig(),vcPair.getVc().getStart(),p.toString()));
                    if (key.startsWith(GATKVCFConstants.ALLELE_SPECIFIC_ANNOTATION_PREFIX)) {
                        AS_Warning.warn(String.format("Reducible annotation '%s' detected, add -G StandardAnnotation -G AS_StandardAnnotation to the command to annotate in the final VC with this annotation.",key));
                    }
                }
            }
        }
    }

    // Use the VCF header's declared type for the given attribute to ensure that all the values for that attribute
    // across all the VCs being merged have the same boxed representation. Some VCs have a serialized value of "0"
    // for FLOAT attributes, with no embedded decimal point, but we still need to box those into Doubles, or the
    // subsequent sorting required to obtain the median will fail due to the list having a mix of Comparable<Integer>
    // and Comparable<Double>.
    private Comparable<?> parseNumericInfoAttributeValue(final VCFHeader vcfHeader, final String key, final String stringValue) {
        final VCFInfoHeaderLine infoLine = vcfHeader.getInfoHeaderLine(key);
        if (infoLine == null) {
            oneShotHeaderLineLogger.warn(String.format("At least one attribute was found (%s) for which there is no corresponding header line", key));
            if (stringValue.contains(".")) {
                return Double.parseDouble(stringValue);
            } else {
                return Integer.parseInt(stringValue);
            }
        }
        switch (infoLine.getType()) {
            case Integer:
                return Integer.parseInt(stringValue);
            case Float:
                return Double.parseDouble(stringValue);
            default:
                throw new NumberFormatException(
                        String.format(
                                "The VCF header specifies type %s type for INFO attribute key %s, but a numeric value is required",
                                infoLine.getType().name(),
                                key)
                );
        }
    }

    /**
     * prefix an allele with additional reference bases if extraBaseCount > 0
     */
    private static Allele extendAllele(Allele allele, int extraBaseCount, byte[] refBases) {
        if (extraBaseCount > 0) {
            final byte[] oldBases = allele.getBases();
            final byte[] newBases = Arrays.copyOf(oldBases, oldBases.length + extraBaseCount);
            System.arraycopy(refBases, refBases.length - extraBaseCount, newBases, oldBases.length, extraBaseCount);
            return Allele.create(newBases,false);
        } else {
            return allele;
        }
    }


    /**
     * Merge into the context a new genotype represented by the given VariantContext for the provided list of target alleles.
     * This method assumes that none of the alleles in the VC overlaps with any of the alleles in the set.
     *  @param vc                    the Variant Context for the sample
     * @param remappedAlleles       the list of remapped alleles for the sample
     * @param targetAlleles         the list of target alleles
     * @param samplesAreUniquified  true if sample names have been uniquified
     * @param useRemappedAllelesForGenotyping  if true, remap alleles prior to making the genotype call
     */
    private GenotypesContext mergeRefConfidenceGenotypes(final VariantContext vc,
                                                           final List<Allele> remappedAlleles,
                                                           final List<Allele> targetAlleles,
                                                           final boolean samplesAreUniquified,
                                                           final boolean useRemappedAllelesForGenotyping) {
        final GenotypesContext mergedGenotypes = GenotypesContext.create();
        final int maximumPloidy = vc.getMaxPloidy(GATKVariantContextUtils.DEFAULT_PLOIDY);
        // the map is different depending on the ploidy, so in order to keep this method flexible (mixed ploidies)
        // we need to get a map done (lazily inside the loop) for each ploidy, up to the maximum possible.
        final int[][] genotypeIndexMapsByPloidy = new int[maximumPloidy + 1][];

        for ( final Genotype g : vc.getGenotypes() ) {
            final String name;
            if (samplesAreUniquified) {
                name = g.getSampleName() + "." + vc.getSource();
            } else {
                name = g.getSampleName();
            }
            final int ploidy = g.getPloidy();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g);
            if (!doSomaticMerge) {
                //do attribute subsetting
                if (g.hasPL() || g.hasAD()) {
                    int[] perSampleIndexesOfRelevantAlleles = AlleleSubsettingUtils.getIndexesOfRelevantAllelesForGVCF(remappedAlleles, targetAlleles, vc.getStart(), g, false);
                    if (g.hasPL()) {
                        final int[] genotypeIndexMapByPloidy = genotypeIndexMapsByPloidy[ploidy] == null
                                ? GenotypeIndexCalculator.newToOldGenotypeMap(ploidy, perSampleIndexesOfRelevantAlleles) //probably horribly slow
                                : genotypeIndexMapsByPloidy[ploidy];
                        genotypeBuilder.PL(generatePL(g, genotypeIndexMapByPloidy));
                    }
                    if (g.hasAD()) {
                        genotypeBuilder.AD(AlleleSubsettingUtils.generateAD(g.getAD(), perSampleIndexesOfRelevantAlleles));
                    }

               }
                //clean up low confidence hom refs for better annotations later
                if (GenotypeGVCFsEngine.excludeFromAnnotations(g)) {
                    genotypeBuilder.alleles(Collections.nCopies(ploidy, Allele.NO_CALL));
                }
            }
            else {  // doSomaticMerge
                genotypeBuilder.noAttributes();
                if (g.hasDP()) {
                    genotypeBuilder.DP(g.getDP());
                }

                for (final String key : SOMATIC_FORMAT_ANNOTATIONS_TO_KEEP) {
                    if(g.hasExtendedAttribute(key)) {
                        genotypeBuilder.attribute(key, g.getExtendedAttribute(key));
                    }
                }

                // lazy initialization of the genotype index map by ploidy.
                int[] perSampleIndexesOfRelevantAlleles = AlleleSubsettingUtils.getIndexesOfRelevantAllelesForGVCF(remappedAlleles, targetAlleles, vc.getStart(), g, false);
                final int nonRefIndex = remappedAlleles.indexOf(Allele.NON_REF_ALLELE);
                final int[] AD;
                if (g.hasAD()) {
                    AD = AlleleSubsettingUtils.generateAD(g.getAD(), perSampleIndexesOfRelevantAlleles);
                    genotypeBuilder.AD(AD);
                } else if (g.hasDP()) {
                    AD = new int[targetAlleles.size()];
                    AD[0] = g.getDP();
                    genotypeBuilder.AD(AD);
                }
                if (g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY)) {  //homRef calls don't have AF
                    final double[] AF = AlleleSubsettingUtils.generateAF(VariantContextGetters.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> new double[]{0.0}, 0.0), perSampleIndexesOfRelevantAlleles);
                    genotypeBuilder.attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, AF);
                }
                else if ((g.isHomRef() || g.isNoCall()) && vc.getAlternateAlleles().size() == 1) {  //homRef blocks don't get an AF so assign it here; multi-sample GVCFs will have no-call GTs for ref blocks
                    genotypeBuilder.attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, new double[targetAlleles.size()-1]);
                }

                for (final String key : SOMATIC_INFO_ANNOTATIONS_TO_MOVE) {
                    setPerSampleSomaticAttributes(vc, perSampleIndexesOfRelevantAlleles, g, genotypeBuilder, key);
                }

                if (!dropSomaticFilteringAnnotations) {
                    for (final String key : Mutect2FilteringEngine.STANDARD_MUTECT_INFO_FIELDS_FOR_FILTERING) {
                        setPerSampleSomaticAttributes(vc, perSampleIndexesOfRelevantAlleles, g, genotypeBuilder, key);
                    }
                }
                //only copy filter status for single-sample VCs -- multi-sample VCs should already have GF updated
                if (vc.filtersWereApplied() && vc.getSampleNames().size() == 1 && !g.isHomRef()) {
                    //PASS has to have null filters, so we can't add an empty list
                    if (!vc.getFilters().isEmpty()) {
                        genotypeBuilder.filters(new ArrayList(vc.getFilters()));
                    }
                }
            }
            genotypeBuilder.name(name);
            final List<Allele> originalGTAlleles = useRemappedAllelesForGenotyping ?
                    g.getAlleles().stream().map(a -> (vc.getAlleleIndex(a) > -1 ? remappedAlleles.get(vc.getAlleleIndex(a)) : Allele.NO_CALL)).collect(Collectors.toList()) :
                    g.getAlleles();

            final GenotypeAssignmentMethod assignmentMethod;
            if (callGTAlleles && GenotypeUtils.shouldBeCalled(g)) {
                assignmentMethod = GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL;
            } else {
                assignmentMethod = GenotypeAssignmentMethod.SET_TO_NO_CALL;
            }
            GATKVariantContextUtils.makeGenotypeCall(g.getPloidy(),
                    genotypeBuilder, assignmentMethod,
                    g.hasLikelihoods() ? g.getLikelihoods().getAsVector() : null,
                    targetAlleles, new GenotypeBuilder(g).alleles(originalGTAlleles).make(), null);
            mergedGenotypes.add(genotypeBuilder.make());
        }

        return mergedGenotypes;
    }

    /**
     * Subset per-allele annotations based on the indices of the relevant alleles to keep
     * Modifies the input GenotypeBuilder with the revised annotations
     * If supplied annotation is for INFO field, it will be moved to FORMAT
     * @param vc
     * @param perSampleIndexesOfRelevantAlleles
     * @param g
     * @param genotypeBuilder
     * @param somaticKey    annotation key can be for INFO or FORMAT
     */
    private void setPerSampleSomaticAttributes(VariantContext vc, int[] perSampleIndexesOfRelevantAlleles, Genotype g, GenotypeBuilder genotypeBuilder, String somaticKey) {
        final VCFHeaderLineCount attributeCount;
        if (g.hasExtendedAttribute(somaticKey)) {  //if this is a multi-sample GVCF input, there may already be a FORMAT annotation
            attributeCount = vcfInputHeader.getFormatHeaderLine(somaticKey).getCountType();
            List<Object> originalTLODs = VariantContextGetters.attributeToList(g.getAnyAttribute(somaticKey));
            final Object newTLODs = generateAnnotationValueVector(attributeCount, originalTLODs, perSampleIndexesOfRelevantAlleles);
            if (newTLODs != null) {
                genotypeBuilder.attribute(somaticKey, newTLODs);
            }
        }
        else if (vc.hasAttribute(somaticKey)) {
            attributeCount = vcfInputHeader.getInfoHeaderLine(somaticKey).getCountType();
            final List<Object> originalTLODs = VariantContextGetters.attributeToList(vc.getAttribute(somaticKey));
            final Object newTLODs = generateAnnotationValueVector(attributeCount, originalTLODs, perSampleIndexesOfRelevantAlleles);
            if (newTLODs != null) {
                genotypeBuilder.attribute(somaticKey, newTLODs);
            }
        }
    }

    /**
     * Composes a new likelihood array given the original genotype and the genotype index map.
     *
     * @param g the original genotype.
     * @param genotypeIndexMapByPloidy genotype index map. The ith element indicates what genotype in {@code g} corresponds
     *                                 to the ith genotype in the return likelihoods array.
     *
     * @throws NullPointerException if {@code g} or {@code genotypeIndexMapByPloidy} is {@code null}, or if {@code g}
     *    does not contain likelihoods.
     * @throws IndexOutOfBoundsException if {@code genotypeIndexMapByPloidy} contain non valid
     *  genotype indices given the likelihood array in {@code g}.
     *
     * @return never {@code null} but an array of exactly {@code genotypeIndexMapByPloidy.length} positions.
     */
    private static int[] generatePL(final Genotype g, final int[] genotypeIndexMapByPloidy) {
        final int[] PLs = new int[genotypeIndexMapByPloidy.length];
        final int[] oldPLs = g.getPL();
        for (int i = 0; i < PLs.length; i++) {
            PLs[i] = oldPLs[genotypeIndexMapByPloidy[i]];
        }
        return PLs;
    }

    /**
     * Generates a new annotation value array by adding zeros for missing alleles given the set of indexes of the Genotype's current
     * alleles from the original annotation value array.
     *
     * @param originalList    the array containing the pre-parsed, original TLOD(s)
     * @param indexesOfRelevantAlleles the indexes of the original alleles corresponding to the new alleles
     * @return a List of new annotation values, may be null
     */
    @VisibleForTesting
    public static <T> List<T> generateAnnotationValueVector(VCFHeaderLineCount alleleCount,
                                                final List<T> originalList, final int[] indexesOfRelevantAlleles) {
        List<T> newLODs = null;
        if (alleleCount.equals(VCFHeaderLineCount.A)) {
            newLODs = AlleleSubsettingUtils.remapALengthList(originalList, indexesOfRelevantAlleles, null);
        } else if (alleleCount.equals(VCFHeaderLineCount.R)) {
            newLODs = AlleleSubsettingUtils.remapRLengthList(originalList, indexesOfRelevantAlleles, null);
        } else {  //count doesn't depend on alleles
            newLODs = originalList;
        }
        return newLODs;
    }

}