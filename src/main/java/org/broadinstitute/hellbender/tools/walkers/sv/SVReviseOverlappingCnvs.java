package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.*;

/**
 * Completes a series of cleaning steps for a VCF produced by the GATK-SV pipeline.
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>
 *         TODO
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>
 *         TODO
 *     </li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <pre>
 *     TODO
 * </pre>
 *
 * <h3>Processing Steps</h3>
 * <ol>
 *     <li>
 *         TODO
 *     </li>
 * </ol>
 */
@CommandLineProgramProperties(
        summary = "Clean and format SV VCF",
        oneLineSummary = "Clean and format SV VCF",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVReviseOverlappingCnvs extends MultiplePassVariantWalker {
    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;

    // Data structure for overlap detection
    private static final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();

    // Data structure for mutliallelic CNV detection
    private static final Set<String> multiCnvs = new HashSet<>();

    // Data structures for revising genotypes
    private static final Map<String, Map<String, Pair<String, String>>> revisedEventsAll = new HashMap<>();
    private static final Map<String, Set<String>> revisedEventsFiltered = new HashMap<>();
    private static final Map<String, Map<String, Integer>> currentCopyNumbers = new HashMap<>();

    // Data structures for revising copy numbers
    private static final Map<String, Set<String>> abnormalRdCn = new HashMap<>();
    private static final Map<String, Map<String, Integer>> revisedCopyNumbers = new HashMap<>();
    private static final Set<String> revisedComplete = new HashSet<>();

    // Data structures for cached data
    private final Map<String, Set<String>> nonRefSamplesCache = new HashMap<>();
    private final Map<String, Map<String, Set<String>>> supportCache = new HashMap<>();
    private final Map<String, Map<String, Integer>> rdCnCache = new HashMap<>();

    private static final int MIN_VARIANT_SIZE = 5000;

    @Override
    protected int numberOfPasses() { return 3; }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            processCollectedVariants();
            clearAllCaches();
        }
    }

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outputVcf);
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MULTI_CNV, 0, VCFHeaderLineType.Flag, "Variant is a multiallelic CNV"));
        vcfWriter.writeHeader(header);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant, final ReadsContext readsContext,
                                final ReferenceContext referenceContext, final FeatureContext featureContext, final int n) {
        switch (n) {
            case 0:
                firstPassApply(variant);
                break;
            case 1:
                secondPassApply(variant);
                break;
            case 2:
                thirdPassApply(variant);
                break;
        }
    }

    public void firstPassApply(final VariantContext variant) {
        // Skip processing if not CNV
        if (!isDelDup(variant)) {
            return;
        }

        // Flag variant as being a multiallelic CNV
        final boolean isDel = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        for (final Genotype genotype : variant.getGenotypes()) {
            if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

            final int rdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            if ((isDel && rdCn > 3) || (!isDel && (rdCn < 1 || rdCn > 4))) {
                multiCnvs.add(variant.getID());
                break;
            }
        }

        // Skip processing if below size threshold
        if (!isLarge(variant, MIN_VARIANT_SIZE)) {
            return;
        }

        // Flag sample as having an abnormal copy number
        for (final Genotype genotype : variant.getGenotypes()) {
            if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

            final int rdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && rdCn < 2) || (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && rdCn > 2)) {
                abnormalRdCn.computeIfAbsent(variant.getID(), k -> new HashSet<>()).add(genotype.getSampleName());
            }
        }

        // Remove variants not in current context window
        overlappingVariantsBuffer.removeIf(overlappingVariant -> {
            boolean shouldRemove = !overlappingVariant.getContig().equals(variant.getContig())
                    || (overlappingVariant.getStart() + overlappingVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)) < variant.getStart();
            if (shouldRemove) removeVariantFromCaches(overlappingVariant.getID());
            return shouldRemove;
        });

        // Process overlaps with variants in the buffer
        for (final VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(bufferedVariant, variant)) {
                processGt(bufferedVariant, variant);
                processCn(bufferedVariant, variant);
            }
        }
        overlappingVariantsBuffer.add(variant);
    }

    public void secondPassApply(final VariantContext variant) {
        // Skip processing if not in revised events map
        if (!revisedEventsFiltered.containsKey(variant.getID())) {
            return;
        }

        // Initialize data structures
        final String variantId = variant.getID();
        final Set<String> samples = revisedEventsFiltered.get(variantId);
        final Map<String, Integer> variantRdCn = new HashMap<>();

        // Initialize revisedRdCn values for each variant
        for (final String sampleName : samples) {
            final Genotype genotype = variant.getGenotype(sampleName);
            if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

            final int rdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            variantRdCn.put(sampleName, rdCn);
        }
        currentCopyNumbers.put(variantId, variantRdCn);
    }

    public void thirdPassApply(final VariantContext variant) {
        final VariantContextBuilder builder = new VariantContextBuilder(variant);

        // Revise genotypes
        if (revisedEventsAll.containsKey(variant.getID())) {
            processRevisedGt(builder, variant);
        }

        // Revise copy numbers
        if (revisedCopyNumbers.containsKey(variant.getID())) {
            processRevisedCn(builder, variant);
        }

        // Tag multiallelic CNVs
        if (multiCnvs.contains((variant.getID()))) {
            builder.attribute(GATKSVVCFConstants.MULTI_CNV, true);
        }

        vcfWriter.add(builder.make());
    }

    private void processGt(final VariantContext v1, final VariantContext v2) {
        // Determine larger variant, swapping if necessary
        VariantContext largerVariant = v1;
        VariantContext smallerVariant = v2;
        final int length1 = v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        final int length2 = v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        if (length2 > length1) {
            largerVariant = v2;
            smallerVariant = v1;
        }

        // Skip if coverage below expected
        final double coverage = getCoverage(largerVariant, smallerVariant);
        if (coverage < 0.5) {
            return;
        }

        // Skip if same variant ID, SV type or sample sets
        final String largerId = largerVariant.getID();
        final String smallerId = smallerVariant.getID();
        final String largerSvType = largerVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final String smallerSvType = smallerVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final Set<String> largerSamples = getNonReferenceSamples(largerVariant);
        final Set<String> smallerSamples = getNonReferenceSamples(smallerVariant);
        if (largerId.equals(smallerId) || largerSvType.equals(smallerSvType) || largerSamples.equals(smallerSamples)) {
            return;
        }

        // Skip if no non-overlapping samples
        final Set<String> nonCommonSamples = new HashSet<>(largerSamples);
        nonCommonSamples.removeAll(smallerSamples);
        if (nonCommonSamples.isEmpty()) {
            return;
        }

        // Add variant pair to data structure
        for (final String sample : nonCommonSamples) {
            revisedEventsAll.computeIfAbsent(smallerId, k -> new HashMap<>())
                    .put(sample, new ImmutablePair<>(largerId, largerSvType));
        }
    }

    private void processCn(final VariantContext v1, final VariantContext v2) {
        // Determine larger variant, swapping if necessary
        VariantContext largerVariant = v1;
        VariantContext smallerVariant = v2;
        int length1 = v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        int length2 = v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        if (length2 > length1) {
            largerVariant = v2;
            smallerVariant = v1;
            length1 = length2;
        }

        // Calculate overlap
        final double coverage = getCoverage(largerVariant, smallerVariant);
        if (coverage < 0.5) {
            return;
        }

        // Skip if no common abnormal samples
        final Set<String> samples = new HashSet<>(abnormalRdCn.getOrDefault(largerVariant.getID(), Collections.emptySet()));
        samples.retainAll(abnormalRdCn.getOrDefault(smallerVariant.getID(), Collections.emptySet()));
        if (samples.isEmpty()) {
            return;
        }

        // Cached non-boolean fields
        final String largerId = largerVariant.getID();
        final String smallerId = smallerVariant.getID();
        final Map<String, Integer> largerRdCn = getRdCn(largerVariant);
        final Map<String, Integer> smallerRdCn = getRdCn(smallerVariant);
        final Map<String, Set<String>> largerSupport = getSupport(largerVariant);
        final Map<String, Set<String>> smallerSupport = getSupport(smallerVariant);
        final String largerSvType = largerVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final String smallerSvType = smallerVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");

        // Cached length fields
        final int minEnd = Math.min(
                largerVariant.getStart() + largerVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0),
                smallerVariant.getStart() + smallerVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)
        );
        final int maxStart = Math.max(largerVariant.getStart(), smallerVariant.getStart());
        final int lengthOverlap = minEnd - maxStart + 1;
        final double largerOverlap = (double) lengthOverlap / (double) length1;

        // Cached boolean fields
        final boolean largerIsMultiCnv = multiCnvs.contains(largerId);
        final boolean smallerIsMultiCnv = multiCnvs.contains(smallerId);
        final boolean isMatchingSvType = largerSvType.equals(smallerSvType);
        final boolean isOverlapping = (largerOverlap > 0.5);
        final boolean isLargerThanMin = length1 > MIN_VARIANT_SIZE;

        // Iterate through samples to test against conditions
        for (final String sample : samples) {
            final String largerFullId = largerId + "@" + sample;
            final String smallerFullId = smallerId + "@" + sample;
            if (revisedComplete.contains(largerFullId)) {
                continue;
            }

            // Initialize variables for evaluation
            final int largerSampleRdCn = revisedCopyNumbers.getOrDefault(largerId, Collections.emptyMap()).getOrDefault(sample, largerRdCn.get(sample));
            final int smallerSampleRdCn = revisedCopyNumbers.getOrDefault(smallerId, Collections.emptyMap()).getOrDefault(sample, smallerRdCn.get(sample));
            final Set<String> largerSampleSupport = largerSupport.get(sample);
            final Set<String> smallerSampleSupport = smallerSupport.get(sample);
            final Genotype genotype2 = smallerVariant.getGenotype(sample);

            // Condition 1: Smaller depth call is being driven by larger call
            if (largerSampleSupport.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && largerSampleSupport.size() > 1
                    && smallerSampleSupport.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1))) && !largerIsMultiCnv) {
                if (largerSampleRdCn == 0) {
                    makeRevision(smallerFullId, smallerSampleRdCn + 2);
                } else if (largerSampleRdCn == 1) {
                    makeRevision(smallerFullId, smallerSampleRdCn + largerSampleRdCn);
                } else if (largerSampleRdCn > 1) {
                    int newCN = smallerSampleRdCn - largerSampleRdCn + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(smallerFullId, newCN);
                }
            }

            // Condition 2: Smaller CNV is driven by larger CNV genotype
            else if (smallerSampleSupport.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && smallerSampleSupport.size() > 1
                    && largerSampleSupport.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && !genotype2.isHomRef() && !smallerIsMultiCnv && isOverlapping) {
                if (smallerSampleRdCn == 0) {
                    makeRevision(largerFullId, largerSampleRdCn + 2);
                } else if (smallerSampleRdCn == 1) {
                    makeRevision(largerFullId, largerSampleRdCn + smallerSampleRdCn);
                } else if (smallerSampleRdCn > 1) {
                    int newCN = largerSampleRdCn - smallerSampleRdCn + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(largerFullId, newCN);
                }
            }

            // Condition 3: Depth-only calls where smaller call is driven by larger call
            else if (largerSampleSupport.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && smallerSampleSupport.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && !largerIsMultiCnv && isMatchingSvType) {
                if (largerSampleRdCn == 0 && largerSampleRdCn != smallerSampleRdCn) {
                    makeRevision(smallerFullId, smallerSampleRdCn + 2);
                } else if (largerSampleRdCn == 1 && largerSampleRdCn > smallerSampleRdCn) {
                    makeRevision(smallerFullId, 1);
                } else if (largerSampleRdCn > 1 && largerSampleRdCn < smallerSampleRdCn) {
                    makeRevision(smallerFullId, Math.max(smallerSampleRdCn - largerSampleRdCn + 2, 0));
                } else {
                    makeRevision(smallerFullId, 2);
                }
            }

            // Condition 4: Any other time a larger call drives a smaller call
            else if (largerSampleSupport.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && !largerIsMultiCnv && isLargerThanMin) {
                if (largerSampleRdCn == 0) {
                    makeRevision(smallerFullId, smallerSampleRdCn + 2);
                } else if (largerSampleRdCn == 1) {
                    makeRevision(smallerFullId, smallerSampleRdCn + largerSampleRdCn);
                } else if (largerSampleRdCn > 1) {
                    int newCN = smallerSampleRdCn - largerSampleRdCn + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(smallerFullId, newCN);
                }
            }
        }
    }

    private void processRevisedGt(final VariantContextBuilder builder, final VariantContext variant) {
        // Initialize data structures
        final String variantId = variant.getID();
        final Map<String, Pair<String, String>> variantEvents = revisedEventsAll.get(variantId);
        final Map<String, Genotype> revisedGenotypes = new HashMap<>();
        final List<Genotype> oldGenotypes = variant.getGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(oldGenotypes.size());

        // Populate genotypes that need revising
        for (final Map.Entry<String, Pair<String, String>> entry : variantEvents.entrySet()) {
            final Pair<String, String> event = entry.getValue();
            if (event == null) {
                continue;
            }

            final String sampleName = entry.getKey();
            final Genotype genotype = variant.getGenotype(sampleName);
            final String largerId = event.getLeft();
            final String largerSvType = event.getRight();
            final int currentRdCn = currentCopyNumbers.get(variantId).getOrDefault(sampleName, 0);
            final int largerRdCn = currentCopyNumbers.getOrDefault(largerId, new HashMap<>()).getOrDefault(sampleName, 0);

            int newVal = -1;
            if (largerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && currentRdCn == 2 && largerRdCn == 3) {
                newVal = 1;
            } else if (largerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && currentRdCn == 2 && largerRdCn == 1) {
                newVal = 3;
            }

            if (newVal != -1) {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

                final int rdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
                gb.GQ(rdCn);
                revisedGenotypes.put(sampleName, gb.make());
            }
        }

        // Populate genotypes that don't need revising
        for (final Genotype genotype : oldGenotypes) {
            final String sampleName = genotype.getSampleName();
            newGenotypes.add(revisedGenotypes.getOrDefault(sampleName, genotype));
        }

        builder.genotypes(newGenotypes);
    }

    private void processRevisedCn(final VariantContextBuilder builder, final VariantContext variant) {
        // Initialize data structures
        final String variantId = variant.getID();
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());

        // Replace revised alleles and copy numbers
        for (final Genotype genotype : genotypes) {
            final String sampleName = genotype.getSampleName();
            if (revisedCopyNumbers.get(variantId).containsKey(sampleName)) {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                gb.attribute(GATKSVVCFConstants.RD_CN, revisedCopyNumbers.get(variantId).get(sampleName));
                updatedGenotypes.add(gb.make());
            } else {
                updatedGenotypes.add(genotype);
            }
        }
        builder.genotypes(updatedGenotypes);
    }

    private void processCollectedVariants() {
        // Prune variant-sample pairs we need RD_CN values for
        for (final Map.Entry<String, Map<String, Pair<String, String>>> entry : revisedEventsAll.entrySet()) {
            for (final Map.Entry<String, Pair<String, String>> innerEntry : entry.getValue().entrySet()) {
                final String sampleName = innerEntry.getKey();
                final String variantId = entry.getKey();
                final String largerVariantId = innerEntry.getValue().getLeft();
                final String largerSvType = innerEntry.getValue().getRight();
                if (largerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) || largerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL)) {
                    revisedEventsFiltered.computeIfAbsent(variantId, k -> new HashSet<>()).add(sampleName);
                    revisedEventsFiltered.computeIfAbsent(largerVariantId, k -> new HashSet<>()).add(sampleName);
                }
            }
        }
    }

    private Set<String> getNonReferenceSamples(final VariantContext variant) {
        final String variantId = variant.getID();
        if (nonRefSamplesCache.containsKey(variantId)) {
            return nonRefSamplesCache.get(variantId);
        }

        final Set<String> samples = new HashSet<>();
        for (final Genotype genotype : variant.getGenotypes()) {
            if (genotype.isCalled() && !genotype.isHomRef()) {
                samples.add(genotype.getSampleName());
            }
        }

        nonRefSamplesCache.put(variantId, samples);
        return samples;
    }

    private Map<String, Set<String>> getSupport(final VariantContext variant) {
        final String variantId = variant.getID();
        if (supportCache.containsKey(variantId)) {
            return supportCache.get(variantId);
        }

        Map<String, Set<String>> supportMap = new HashMap<>();
        for (final Genotype genotype : variant.getGenotypes()) {
            final String supportStr = genotype.hasExtendedAttribute(GATKSVVCFConstants.EV)
                    ? genotype.getExtendedAttribute(GATKSVVCFConstants.EV).toString()
                    : "";
            final Set<String> supportSet = new HashSet<>();
            if (!supportStr.isEmpty()) {
                supportSet.addAll(Arrays.asList(supportStr.split(",")));
            }
            supportMap.put(genotype.getSampleName(), supportSet);
        }

        supportCache.put(variantId, supportMap);
        return supportMap;
    }

    private Map<String, Integer> getRdCn(final VariantContext variant) {
        final String variantId = variant.getID();
        if (rdCnCache.containsKey(variantId)) {
            return rdCnCache.get(variantId);
        }

        final Map<String, Integer> rdCnMap = new HashMap<>();
        for (final Genotype genotype : variant.getGenotypes()) {
            if (genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                rdCnMap.put(genotype.getSampleName(), Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
            }
        }

        rdCnCache.put(variantId, rdCnMap);
        return rdCnMap;
    }

    private void clearAllCaches() {
        nonRefSamplesCache.clear();
        supportCache.clear();
        rdCnCache.clear();
    }

    private void removeVariantFromCaches(final String variantID) {
        nonRefSamplesCache.remove(variantID);
        supportCache.remove(variantID);
        rdCnCache.remove(variantID);
    }

    private void makeRevision(final String id, final int val) {
        final String[] tokens = id.split("@");
        final String variantId = tokens[0];
        final String sample = tokens[1];
        revisedCopyNumbers.computeIfAbsent(variantId, k -> new HashMap<>()).put(sample, val);
        if (val == 2) {
            revisedComplete.add(id);
        }
    }

    private boolean isDelDup(final VariantContext variant) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        return svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP);
    }

    private boolean isLarge(final VariantContext variant, final int minSize) {
        final int variantLength = Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
        return variantLength >= minSize;
    }

    private boolean overlaps(final VariantContext v1, final VariantContext v2) {
        return v1.getContig().equals(v2.getContig())
                && v1.getStart() <= (v2.getStart() + v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0))
                && v2.getStart() <= (v1.getStart() + v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
    }

    private double getCoverage(final VariantContext larger, final VariantContext smaller) {
        final int largerStart = larger.getStart();
        final int smallerStart = smaller.getStart();
        final int largerStop = larger.getEnd();
        final int smallerStop = smaller.getEnd();

        if (largerStart <= smallerStop && smallerStart <= largerStop) {
            final int intersectionSize = Math.min(smallerStop, largerStop) - Math.max(smallerStart, largerStart) + 1;
            return (double) intersectionSize / (smallerStop - smallerStart + 1);
        }
        return 0.0;
    }
}