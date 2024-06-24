package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
 *
 * <p>
 * CombineGVCFs is meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs.
 * One could use this tool to genotype multiple individual GVCFs instead of GenomicsDBImport; one would first use
 * CombineGVCFs to combine them into a single GVCF and pass the results into GenotypeGVCFs. The main advantage of using CombineGVCFs
 * over GenomicsDBImport is the ability to combine multiple intervals at once without building a GenomicsDB. CombineGVCFs
 * is slower than GenomicsDBImport though, so it is recommended CombineGVCFs only be used when there are few samples to merge.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more HaplotypeCaller GVCFs to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined multi-sample gVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CombineGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf.gz \
 *   --variant sample2.g.vcf.gz \
 *   -O cohort.g.vcf.gz
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only GVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call GVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 * <p>If the GVCF files contain allele specific annotations, add `-G StandardAnnotation -G AS_StandardAnnotation` to the command line.</p>
 *
 * <p>Users generating large callsets (1000+ samples) may prefer GenomicsDBImport, which uses Intel's GenomicsDB and is capable of scaling to much larger sample sizes than CombineGVCFs.
 * This tool provides a pure java reference implementation of the combine operation which is available on all architectures.<p/>
 *
 */
@CommandLineProgramProperties(summary = "Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations", oneLineSummary = "Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations", programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
@WorkflowProperties
public final class CombineGVCFs extends MultiVariantWalkerGroupedOnStart {

    private VariantAnnotatorEngine annotationEngine;
    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary sequenceDictionary;
    private ReferenceConfidenceVariantContextMerger referenceConfidenceVariantContextMerger;

    public static final String BP_RES_LONG_NAME = "convert-to-base-pair-resolution";
    public static final String BREAK_BANDS_LONG_NAME = "break-bands-at-multiples-of";
    public static final String SOMATIC_INPUT_LONG_NAME = "input-is-somatic";
    public static final String DROP_SOMATIC_FILTERING_ANNOTATIONS_LONG_NAME = "drop-somatic-filtering-annotations";
    public static final String ALLELE_FRACTION_DELTA_LONG_NAME = "allele-fraction-error";

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined GVCF output file", optional=false)
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    private GATKPath outputFile;

    @Argument(fullName=BP_RES_LONG_NAME, doc = "If specified, convert banded gVCFs to all-sites gVCFs", optional=true)
    protected boolean useBpResolution = false;

    /**
     * To reduce file sizes our gVCFs group similar reference positions into bands.  However, there are cases when users will want to know that no bands
     * span across a given genomic position (e.g. when scatter-gathering jobs across a compute farm).  The option below enables users to break bands at
     * pre-defined positions.  For example, a value of 10,000 would mean that we would ensure that no bands span across chr1:10000, chr1:20000, etc.
     *
     * Note that the --convert-to-base-pair-resolution argument is just a special case of this argument with a value of 1.
     */
    @Argument(fullName=BREAK_BANDS_LONG_NAME, doc = "If > 0, reference bands will be broken up at genomic positions that are multiples of this number", optional=true)
    protected int multipleAtWhichToBreakBands = 0;

    /**
     * Merge somatic GVCFs, retaining LOD and haplotype event count information in FORMAT field
     * Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.
     */
    @Argument(fullName= SOMATIC_INPUT_LONG_NAME, doc = "Merge input GVCFs according to somatic (i.e. Mutect2) annotations (BETA)")
    protected boolean somaticInput = false;

    /**
     * Rather than move the per-sample INFO annotations used for filtering to the FORMAT field, drop them entirely.
     * This makes the FORMAT field more readable and reduces file sizes for large cohorts.
     */
    @Argument(fullName=DROP_SOMATIC_FILTERING_ANNOTATIONS_LONG_NAME, doc = "For input somatic GVCFs (i.e. from Mutect2) drop filtering annotations")
    protected boolean dropSomaticFilteringAnnotations = false;

    /**
     * By default CombineGVCFs reverts all genotypes to no-calls, but calls can be made if specified
     */
    @Argument(fullName = GenomicsDBArgumentCollection.CALL_GENOTYPES_LONG_NAME, doc = "Output called genotypes?", optional = true)
    protected boolean makeGenotypeCalls = false;

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Collections.singletonList(StandardAnnotation.class);
    }

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // State that gets accumulated between calls of apply()
    private final LinkedList<VariantContext> variantContextsOverlappingCurrentMerge = new LinkedList<>();
    private final Set<String> samples = new HashSet<>();
    private SimpleInterval prevPos = null;
    private byte refAfterPrevPos;
    private ReferenceContext storedReferenceContext;

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
        // Check that the input variant contexts do not contain MNPs as these may not be properly merged
        for (final VariantContext ctx : variantContexts) {
            if (GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(ctx)) {
                throw new UserException.BadInput(String.format(
                        "Combining gVCFs containing MNPs is not supported. %1s contained a MNP at %2s:%3d",
                        ctx.getSource(), ctx.getContig(), ctx.getStart()));
            }
        }

        // If we need to stop at an intermediate site since the last apply, do so (caused by gvcfBlocks, contexts ending, etc...)
        if (!variantContextsOverlappingCurrentMerge.isEmpty()) {
            Locatable last = prevPos!=null && prevPos.getContig().equals(variantContextsOverlappingCurrentMerge.get(0).getContig()) ?  prevPos : variantContextsOverlappingCurrentMerge.get(0);
            // If on a different contig, close out all the queued states on the current contig
            int end = last.getContig().equals(referenceContext.getWindow().getContig())
                    ? referenceContext.getInterval().getStart() - 1
                    : variantContextsOverlappingCurrentMerge.stream().mapToInt(VariantContext::getEnd).max().getAsInt();

            createIntermediateVariants( new SimpleInterval(last.getContig(), last.getStart(), end));
        }

        mergeWithNewVCs(variantContexts, referenceContext);

        // Update the stored reference if it has a later stop position than the current stored reference
        if ( (storedReferenceContext == null) ||
                (!referenceContext.getWindow().contigsMatch(storedReferenceContext.getWindow()) ) ||
                (storedReferenceContext.getWindow().getEnd() < referenceContext.getWindow().getEnd())) {
            storedReferenceContext = referenceContext;
        }
    }

    /**
     * calculates if there are any sites in the provided interval where we should expect the tool to create
     * a new variant context object by calling endPreviousStates() and closes them by providing appropriate reference
     * information and an empty list of new variant contexts.
     *
     */
    @VisibleForTesting
    void createIntermediateVariants(SimpleInterval intervalToClose) {
        resizeReferenceIfNeeded(intervalToClose);

        // Break up the GVCF according to the provided reference blocking scheme
        // The values returned from getIntermediateStopSites represent a proposed set of stop sites that may include
        // intervals that are outside the actual interval being closed. These sites are filtered out below.
        // Note: Precomputing these is really inefficient when large reference blocks are closed with
        // fine band resolution because it results in very large collections of stop sites (tens or hundreds of millions)
        // that must subsequently be sorted.
        final Set<Integer> sitesToStop = getIntermediateStopSites(intervalToClose, multipleAtWhichToBreakBands);

        // If any variant contexts ended (or were spanning deletions) the last context compute where we should stop them
        for (VariantContext vc : variantContextsOverlappingCurrentMerge) {

            // Asking if the number of alleles > 2 is a shorthand for a variant being present, as we expect <non-ref>
            // symbolic alleles to be present in all VariantContext. This might also be the case if we saw a spanning
            // deletion that reads into the current site, as we would expect ReferenceConfidenceVariantContextMerger to
            // insert symbolic alleles for those spanning variants.
            if (vc.getNAlleles() > 2) {
                for (int i = vc.getStart(); i <= vc.getEnd(); i++ ) {
                    sitesToStop.add(i);
                }
            } else if (vc.getEnd() <= intervalToClose.getEnd()) {
                sitesToStop.add(vc.getEnd());
            }
        }

        List<Integer> stoppedLocs = new ArrayList<>(sitesToStop);
        stoppedLocs.sort(Comparator.naturalOrder());

        // For each stopped loc that is within the interval being closed, create a fake QueuedContextState and pass it to endPreviousStats
        for (int stoppedLoc : stoppedLocs) {
            SimpleInterval loc = new SimpleInterval(intervalToClose.getContig(), stoppedLoc, stoppedLoc);
            if (( stoppedLoc <= intervalToClose.getEnd() && stoppedLoc>= intervalToClose.getStart()) && isWithinInterval(loc)) {
                byte[] refBases = Arrays.copyOfRange(storedReferenceContext.getBases(), stoppedLoc - storedReferenceContext.getWindow().getStart(), stoppedLoc - storedReferenceContext.getWindow().getStart() + 2);
                endPreviousStates(loc, refBases, Collections.emptyList(), true);
            }
        }

    }

    // Get any intermediate stop sites based on the break band multiple.
    @VisibleForTesting
    protected final static Set<Integer> getIntermediateStopSites(final SimpleInterval intervalToClose, final int breakBandMultiple) {
        final Set<Integer> sitesToStop = new HashSet<>();

        if ( breakBandMultiple > 0) {
            // if the intermediate interval to close starts before the end of the first band multiple,
            // create the first stop position at the end of the band multiple
            for (int blockEndPosition = intervalToClose.getStart() < (breakBandMultiple + 1) ?
                    Math.max(2, breakBandMultiple) :
                    (intervalToClose.getStart() / breakBandMultiple) * breakBandMultiple;
                 blockEndPosition <= intervalToClose.getEnd();
                 blockEndPosition += breakBandMultiple) {
                sitesToStop.add(blockEndPosition - 1); // Subtract 1 here because we want to split before this base
            }
        }
        return sitesToStop;
    }

    /**
     * Resize {@link #storedReferenceContext} to cover at least as much as intervalToClose
     * @param intervalToClose
     */
    private void resizeReferenceIfNeeded(SimpleInterval intervalToClose) {
        final int leftEdge = storedReferenceContext.getInterval().getStart() - intervalToClose.getStart();
        final int rightEdge = intervalToClose.getEnd() - storedReferenceContext.getInterval().getEnd();

        storedReferenceContext.setWindow(Math.max(1, leftEdge), Math.max(1, rightEdge));
    }


    @Override
    public void onTraversalStart() {
        if (somaticInput) {
            logger.warn("Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.");
        }

        // create the annotation engine
        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, false);

        vcfWriter = getVCFWriter();

        referenceConfidenceVariantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, getHeaderForVariants(), somaticInput, dropSomaticFilteringAnnotations, makeGenotypeCalls);

        //now that we have all the VCF headers, initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)'
        sequenceDictionary = getBestAvailableSequenceDictionary();

        // optimization to prevent mods when we always just want to break bands
        if ( multipleAtWhichToBreakBands == 1 || useBpResolution) {
            useBpResolution = true;
            multipleAtWhichToBreakBands = 1;
        }
    }

    private VariantContextWriter getVCFWriter() {
        final SortedSet<String> samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // add headers for annotations added by this tool
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        if (somaticInput) {
            //single-sample M2 variant filter status will get moved to genotype filter
            headerLines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));

            if (!dropSomaticFilteringAnnotations) {
                //standard M2 INFO annotations for filtering will get moved to FORMAT field
                for (final String key : Mutect2FilteringEngine.STANDARD_MUTECT_INFO_FIELDS_FOR_FILTERING) {
                    headerLines.add(GATKVCFHeaderLines.getEquivalentFormatHeaderLine(key));
                }
            }
        }

        final VariantContextWriter writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /**
     * Method which calls endPreviousStates at the appropriate places on the given a new startingStates object
     * and an OverallState object corresponding to the currently accumulated reads.
     *
     * @param variantContexts list of variant contexts with the same start position to be reduced
     * @param referenceContext ReferenceContext object overlapping the provided VariantContexts
     */
    public void mergeWithNewVCs(final List<VariantContext> variantContexts, final ReferenceContext referenceContext) {
        if ( !variantContexts.isEmpty() ) {
            if ( ! okayToSkipThisSite(variantContexts, referenceContext) ) {
                SimpleInterval loc = referenceContext.getInterval();
                if (loc.getStart()-1 > 0) {
                    endPreviousStates(new SimpleInterval(loc.getContig(), loc.getStart() - 1, loc.getStart() - 1),
                            Arrays.copyOfRange(referenceContext.getBases(), 1, referenceContext.getWindow().getLengthOnReference()),
                            variantContexts,
                            false);
                }
            }
            variantContextsOverlappingCurrentMerge.addAll(variantContexts);
            for(final VariantContext vc : variantContextsOverlappingCurrentMerge){
                samples.addAll(vc.getSampleNames());
            }
        }
    }

    /**
     * Is it okay to skip the given position?
     *
     * @param variantContexts  the query variant contexts representing the current position
     * @param referenceContext  Reference context object overlapping the variant contexts
     * @return true if it is okay to skip this position, false otherwise
     */
    private boolean okayToSkipThisSite(List<VariantContext> variantContexts, ReferenceContext referenceContext) {
        Set<String> intersection = new HashSet<>(getSamples(variantContexts));
        intersection.retainAll(samples);

        //if there's a starting VC with a sample that's already in a current VC, don't skip this position
        return prevPos != null && referenceContext.getInterval().getStart() == prevPos.getStart() + 1 && intersection.isEmpty();
    }

    private Set<String> getSamples(List<VariantContext> variantContexts) {
        Set<String> output = new HashSet<>();
        for (final VariantContext vc : variantContexts) {
            output.addAll(vc.getSampleNames());
        }
        return output;
    }

    /**
     * Disrupt the VariantContexts so that they all stop at the given pos, write them out, and put the remainder back in the list.
     * @param pos   the position for the starting variantContexts
     * @param variantContexts the current variant contexts with the same starting position
     * @param forceOutputAtCurrentPosition  indicates whether we output a variant at the current position, independent of VCF start/end, i.e. in BP resolution mode
     */
    private void endPreviousStates(final SimpleInterval pos, final byte[] refBases, final List<VariantContext> variantContexts, boolean forceOutputAtCurrentPosition) {
        Set<String> newSamples = getSamples(variantContexts);

        final byte refBase = refBases[0];
        //if we're in BP resolution mode or a VC ends at the current position then the reference for the next output VC (refNextBase)
        // will be advanced one base
        final byte refNextBase = (forceOutputAtCurrentPosition) ? (refBases.length > 1 ? refBases[1] : (byte)'N' ): refBase;

        final List<VariantContext> stoppedVCs = new ArrayList<>(variantContextsOverlappingCurrentMerge.size());

        for (int i = variantContextsOverlappingCurrentMerge.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = variantContextsOverlappingCurrentMerge.get(i);
            //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
            if ( vc.getStart() <= pos.getStart() || !vc.contigsMatch(pos)) {

                stoppedVCs.add(vc);

                // if it was ending anyways, then remove it from the future state
                // or if ending vc is the same sample as a starting VC, then remove it from the future state
                if((vc.getEnd() == pos.getStart()) || (variantContexts.size() > 0 && !forceOutputAtCurrentPosition && newSamples.containsAll(vc.getSampleNames()))) {
                    samples.removeAll(vc.getSampleNames());
                    variantContextsOverlappingCurrentMerge.remove(i);
                }
            }
        }

        //output the stopped variantContexts if there is no previous output (state.prevPos == null) or our current position is past
        // the last write position (state.prevPos)
        //NOTE: BP resolution with have current position == state.prevPos because it gets output via a different control flow
        if ( !stoppedVCs.isEmpty() &&  (prevPos == null || IntervalUtils.isAfter(pos,prevPos,sequenceDictionary) )) {
            final SimpleInterval closingSpot = new SimpleInterval(stoppedVCs.get(0).getContig(), pos.getStart(), pos.getStart());

            // we need the specialized merge if the site contains anything other than ref blocks
            final VariantContext mergedVC;
            if ( containsTrueAltAllele(stoppedVCs) ) {
                mergedVC = referenceConfidenceVariantContextMerger.merge(stoppedVCs, closingSpot, refBase, false, false);
            } else {
                mergedVC = referenceBlockMerge(stoppedVCs, pos.getStart());
            }

            vcfWriter.add(mergedVC);
            prevPos = closingSpot;
            refAfterPrevPos = refNextBase;
        }
    }

    /**
     * Combine a list of reference block VariantContexts.
     * We can't use GATKVariantContextUtils.simpleMerge() because it is just too slow for this sort of thing.
     *
     * @param vcs   the variant contexts to merge
     * @param end   the end of this block (inclusive)
     * @return a new merged VariantContext
     */
    private VariantContext referenceBlockMerge(final List<VariantContext> vcs, final int end) {

        final VariantContext first = vcs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( prevPos == null || !prevPos.getContig().equals(first.getContig()) || first.getStart() >= prevPos.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = prevPos.getStart() + 1;
            refAllele = Allele.create(refAfterPrevPos, true);
        }
        final List<Allele> allelesToUse = Arrays.asList(refAllele, Allele.NON_REF_ALLELE);

        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);
        if ( !useBpResolution && end != start ) {
            attrs.put(VCFConstants.END_KEY, Integer.toString(end));
        }

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        final GenotypeAssignmentMethod assignmentMethod = makeGenotypeCalls ? GenotypeAssignmentMethod.PREFER_PLS : GenotypeAssignmentMethod.SET_TO_NO_CALL;
        for (final VariantContext vc : vcs) {
            for (final Genotype g : vc.getGenotypes()) {
                final GenotypeBuilder gBuilder = new GenotypeBuilder(g);
                GATKVariantContextUtils.makeGenotypeCall(g.getPloidy(),
                            gBuilder, assignmentMethod,
                            g.hasLikelihoods() ? g.getLikelihoods().getAsVector() : null, allelesToUse, g, null);
                genotypes.add(gBuilder.make());
            }
        }
        return new VariantContextBuilder("", first.getContig(), start, end, allelesToUse).attributes(attrs).genotypes(genotypes).make();
    }

    /**
     * Does the given list of VariantContexts contain any with an alternate allele other than <NON_REF>?
     *
     * @param VCs  list of VariantContexts
     * @return true if there are one or more variantContexts that contain a true alternate allele, false otherwise
     */
    private static boolean containsTrueAltAllele(final List<VariantContext> VCs) {

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 ) {
                return true;
            }
        }
        return false;
    }

    @Override
    public Object onTraversalSuccess() {

        if (storedReferenceContext == null) {
            logger.warn("Error: The requested interval contained no data in source VCF files");
            return null;
        }

        if ( !variantContextsOverlappingCurrentMerge.isEmpty() ) {
            // finish off the last blocks
            final SimpleInterval lastInterval = new SimpleInterval(
                    variantContextsOverlappingCurrentMerge.get(0).getContig(),
                    variantContextsOverlappingCurrentMerge.get(0).getStart(),
                    variantContextsOverlappingCurrentMerge.stream().map(VariantContext::getEnd).max(Comparator.naturalOrder()).get());
                createIntermediateVariants(lastInterval);
            // there shouldn't be any state left unless the user cut in the middle of a gVCF block
            if ( !variantContextsOverlappingCurrentMerge.isEmpty() ) {
                logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
            }
        }

        return null;
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
