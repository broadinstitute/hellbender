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

    // Data structures to hold accumulated data across variants
    private static final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();

    private static final Map<String, Map<String, Pair<String, String>>> revisedEventsAll = new HashMap<>();
    private static final Map<String, Set<String>> revisedEventsFiltered = new HashMap<>();
    private static final Map<String, Map<String, Integer>> currentCopyNumbers = new HashMap<>();

    private static final Map<String, Set<String>> abnormalRdCn = new HashMap<>();
    private static final Map<String, Map<String, Integer>> revisedCopyNumbers = new HashMap<>();
    private static final Set<String> revisedComplete = new HashSet<>();

    private static final Set<String> multiCnvs = new HashSet<>();

    private static final int MIN_VARIANT_SIZE = 5000;

    @Override
    protected int numberOfPasses() { return 3; }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            processCollectedVariants();
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
        if (!isDelDup(variant)) {
            return;
        }

        // Flag variant as being a multiallelic CNV if it passes certain conditions
        final boolean isDel = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        for (String sample : variant.getSampleNamesOrderedByName()) {
            final Genotype genotype = variant.getGenotype(sample);
            if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

            final int rdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            if ((isDel && rdCn > 3) || (!isDel && (rdCn < 1 || rdCn > 4))) {
                multiCnvs.add(variant.getID());
                break;
            }
        }

        if (!isLarge(variant, MIN_VARIANT_SIZE)) {
            return;
        }

        // Flag sample as having an abnormal copy number if it passes certain conditions
        for (final String sample : variant.getSampleNames()) {
            final Genotype genotype = variant.getGenotype(sample);
            if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

            final int rdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && rdCn < 2) || (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && rdCn > 2)) {
                abnormalRdCn.computeIfAbsent(variant.getID(), k -> new HashSet<>()).add(sample);
            }
        }

        // Process overlaps with variants in the buffer
        overlappingVariantsBuffer.removeIf(vc -> !vc.getContig().equals(variant.getContig())
                || (vc.getStart() + vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)) < variant.getStart());
        for (final VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(bufferedVariant, variant)) {
                processOverlap(bufferedVariant, variant);
                adjustCopyNumber(bufferedVariant, variant);
            }
        }
        overlappingVariantsBuffer.add(variant);
    }

    public void secondPassApply(final VariantContext variant) {
        if (!revisedEventsFiltered.containsKey(variant.getID())) {
            return;
        }

        // Initialize data structures
        final String variantId = variant.getID();
        final Set<String> samples = revisedEventsFiltered.get(variantId);
        final Map<String, Integer> variantRdCn = new HashMap<>();

        // Initialize revisedRdCn value for each variant
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
        if (revisedEventsAll.containsKey(variant.getID())) {
            processRevisedEvent(builder, variant);
        }
        if (revisedCopyNumbers.containsKey(variant.getID())) {
            processRevisedCn(builder, variant);
        }
        if (multiCnvs.contains((variant.getID()))) {
            builder.attribute(GATKSVVCFConstants.MULTI_CNV, true);
        }
        vcfWriter.add(builder.make());
    }

    private void processCollectedVariants() {
        // Prune variant-sample pairs we need RD_CN values for
        for (final Map.Entry<String, Map<String, Pair<String, String>>> entry : revisedEventsAll.entrySet()) {
            for (final Map.Entry<String, Pair<String, String>> innerEntry : entry.getValue().entrySet()) {
                final String sampleName = innerEntry.getKey();
                final String variantId = entry.getKey();
                final String widerVariantId = innerEntry.getValue().getLeft();
                final String svType = innerEntry.getValue().getRight();
                if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL)) {
                    revisedEventsFiltered.computeIfAbsent(variantId, k -> new HashSet<>()).add(sampleName);
                    revisedEventsFiltered.computeIfAbsent(widerVariantId, k -> new HashSet<>()).add(sampleName);
                }
            }
        }
    }

    private void processOverlap(final VariantContext v1, final VariantContext v2) {
        // Get overlap data
        VariantContext wider;
        VariantContext narrower;
        final int length1 = v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        final int length2 = v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        if (length1 > length2) {
            wider = v1;
            narrower = v2;
        } else if (length2 > length1) {
            wider = v2;
            narrower = v1;
        } else {
            return;
        }
        final String widerID = wider.getID();
        final String narrowerID = narrower.getID();

        // Skip processing if same variant ID, SV type or samples
        final String widerSvType = wider.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final String narrowerSvType = narrower.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final Set<String> widerSamples = getNonReferenceSamples(wider);
        final Set<String> narrowerSamples = getNonReferenceSamples(narrower);
        if (widerID.equals(narrowerID) || widerSvType.equals(narrowerSvType) || widerSamples.equals(narrowerSamples)) {
            return;
        }

        // Get samples present in wider but not in narrower
        final Set<String> nonCommonSamples = new HashSet<>(widerSamples);
        nonCommonSamples.removeAll(narrowerSamples);
        if (nonCommonSamples.isEmpty()) {
            return;
        }

        // Revise variant if coverage exceeds threshold
        final double coverage = getCoverage(wider, narrower);
        if (coverage >= 0.5) {
            for (final String sample : nonCommonSamples) {
                revisedEventsAll.computeIfAbsent(narrowerID, k -> new HashMap<>())
                        .put(sample, new ImmutablePair<>(widerID, widerSvType));
            }
        }
    }

    private void processRevisedEvent(final VariantContextBuilder builder, final VariantContext variant) {
        // Initialize data structures
        final String variantId = variant.getID();
        final Map<String, Pair<String, String>> variantEvents = revisedEventsAll.get(variantId);
        final List<Genotype> newGenotypes = new ArrayList<>();

        // Create updated genotypes
        for (String sample : variant.getSampleNamesOrderedByName()) {
            final Genotype oldGenotype = variant.getGenotype(sample);
            final Pair<String, String> event = variantEvents.get(sample);

            if (event != null) {
                final String widerVariantId = event.getLeft();
                final String widerSvType = event.getRight();
                final int currentRdCn = currentCopyNumbers.get(variantId).getOrDefault(sample, 0);
                final int widerRdCn = currentCopyNumbers.getOrDefault(widerVariantId, new HashMap<>()).getOrDefault(sample, 0);

                int newVal = -1;
                if (widerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && currentRdCn == 2 && widerRdCn == 3) {
                    newVal = 1;
                } else if (widerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && currentRdCn == 2 && widerRdCn == 1) {
                    newVal = 3;
                }

                if (newVal != -1) {
                    final GenotypeBuilder gb = new GenotypeBuilder(oldGenotype);
                    gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                    if (!oldGenotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) continue;

                    final int rdCn = Integer.parseInt(oldGenotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
                    gb.GQ(rdCn);
                    newGenotypes.add(gb.make());
                } else {
                    newGenotypes.add(oldGenotype);
                }
            } else {
                newGenotypes.add(oldGenotype);
            }
        }
        builder.genotypes(newGenotypes);
    }

    private void adjustCopyNumber(final VariantContext v1, final VariantContext v2) {
        // Determine larger variant
        final int length1 = v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        final int length2 = v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        VariantContext largerVariant = v1;
        VariantContext smallerVariant = v2;
        int largerLength = length1;
        int smallerLength = length2;

        // Swap variants if necessary
        if (length2 > length1) {
            largerVariant = v2;
            smallerVariant = v1;
            largerLength = length2;
            smallerLength = length1;
        }

        // Get variant attributes
        final String variantId1 = largerVariant.getID();
        final String variantId2 = smallerVariant.getID();
        final Map<String, Integer> variantRdCn1 = getRdCn(largerVariant);
        final Map<String, Integer> variantRdCn2 = getRdCn(smallerVariant);
        final Map<String, Set<String>> variantSupport1 = getSupport(largerVariant);
        final Map<String, Set<String>> variantSupport2 = getSupport(smallerVariant);
        final String svType1 = largerVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final String svType2 = smallerVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");

        // Calculate overlap
        final int minEnd = Math.min(
                largerVariant.getStart() + largerVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0),
                smallerVariant.getStart() + smallerVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)
        );
        final int maxStart = Math.max(largerVariant.getStart(), smallerVariant.getStart());
        final int lengthOverlap = minEnd - maxStart + 1;
        final double overlap1 = (double) lengthOverlap / (double) largerLength;
        final double overlap2 = (double) lengthOverlap / (double) smallerLength;

        // Get samples with abnormal CN across both variants
        final Set<String> samples = new HashSet<>(abnormalRdCn.getOrDefault(variantId1, Collections.emptySet()));
        samples.retainAll(abnormalRdCn.getOrDefault(variantId2, Collections.emptySet()));

        // Iterate through samples to test against conditions
        for (String sample : samples) {
            final String id1 = variantId1 + "@" + sample;
            final String id2 = variantId2 + "@" + sample;
            if (revisedComplete.contains(id1)) {
                continue;
            }

            // Initialize variables for evaluation
            final int rdCn1 = revisedCopyNumbers.getOrDefault(variantId1, Collections.emptyMap()).getOrDefault(sample, variantRdCn1.get(sample));
            final int rdCn2 = revisedCopyNumbers.getOrDefault(variantId2, Collections.emptyMap()).getOrDefault(sample, variantRdCn2.get(sample));
            final Set<String> support1 = variantSupport1.get(sample);
            final Set<String> support2 = variantSupport2.get(sample);
            final Genotype genotype2 = smallerVariant.getGenotype(sample);

            // Condition 1: Smaller depth call is being driven by larger call
            if (support1.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && support1.size() > 1
                    && support2.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && overlap2 > 0.5 && !largerVariant.hasAttribute(GATKSVVCFConstants.MULTI_CNV)) {
                if (rdCn1 == 0) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1) {
                    makeRevision(id2, rdCn2 + rdCn1);
                } else if (rdCn1 > 1) {
                    int newCN = rdCn2 - rdCn1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                }
            }

            // Condition 2: Smaller CNV is driven by larger CNV genotype
            else if (support1.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && support2.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && support2.size() > 1
                    && overlap1 > 0.5 && overlap2 > 0.5 && !smallerVariant.hasAttribute(GATKSVVCFConstants.MULTI_CNV)
                    && !genotype2.isHomRef()) {
                if (rdCn2 == 0) {
                    makeRevision(id1, rdCn1 + 2);
                } else if (rdCn2 == 1) {
                    makeRevision(id1, rdCn1 + rdCn2);
                } else if (rdCn2 > 1) {
                    int newCN = rdCn1 - rdCn2 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id1, newCN);
                }
            }

            // Condition 3: Depth-only calls where smaller call is driven by larger call
            else if (support1.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && support2.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && overlap2 > 0.5 && !largerVariant.hasAttribute(GATKSVVCFConstants.MULTI_CNV) && svType1.equals(svType2)) {
                if (rdCn1 == 0 && rdCn1 != rdCn2) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1 && rdCn1 > rdCn2) {
                    makeRevision(id2, 1);
                } else if (rdCn1 > 1 && rdCn1 < rdCn2) {
                    makeRevision(id2, Math.max(rdCn2 - rdCn1 + 2, 0));
                } else {
                    makeRevision(id2, 2);
                }
            }

            // Condition 4: Any other time a larger call drives a smaller call
            else if (support1.contains(GATKSVVCFConstants.EV_VALUES.get(1))
                    && overlap2 > 0.5 && !largerVariant.hasAttribute(GATKSVVCFConstants.MULTI_CNV) && largerLength > MIN_VARIANT_SIZE) {
                if (rdCn1 == 0) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1) {
                    makeRevision(id2, rdCn2 + rdCn1);
                } else if (rdCn1 > 1) {
                    int newCN = rdCn2 - rdCn1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                }
            }
        }
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

    private void processRevisedCn(final VariantContextBuilder builder, final VariantContext variant) {
        // Initialize data structures
        final String variantID = variant.getID();
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());

        // Replace revised alleles and copy numbers
        for (final Genotype genotype : genotypes) {
            final String sampleName = genotype.getSampleName();
            if (revisedCopyNumbers.get(variantID).containsKey(sampleName)) {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                gb.attribute(GATKSVVCFConstants.RD_CN, revisedCopyNumbers.get(variantID).get(sampleName));
                updatedGenotypes.add(gb.make());
            } else {
                updatedGenotypes.add(genotype);
            }
        }
        builder.genotypes(updatedGenotypes);
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

    private Set<String> getNonReferenceSamples(final VariantContext variant) {
        final Set<String> samples = new HashSet<>();
        for (final String sampleName : variant.getSampleNames()) {
            final Genotype genotype = variant.getGenotype(sampleName);
            if (genotype.isCalled() && !genotype.isHomRef()) {
                samples.add(sampleName);
            }
        }
        return samples;
    }

    private double getCoverage(final VariantContext wider, final VariantContext narrower) {
        final int nStart = narrower.getStart();
        final int nStop = narrower.getEnd();
        final int wStart = wider.getStart();
        final int wStop = wider.getEnd();

        if (wStart <= nStop && nStart <= wStop) {
            final int intersectionSize = Math.min(nStop, wStop) - Math.max(nStart, wStart) + 1;
            return (double) intersectionSize / (nStop - nStart + 1);
        }
        return 0.0;
    }

    private Map<String, Set<String>> getSupport(final VariantContext variant) {
        Map<String, Set<String>> supportMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            final Genotype genotype = variant.getGenotype(sample);
            final String supportStr = genotype.hasExtendedAttribute(GATKSVVCFConstants.EV) ? genotype.getExtendedAttribute(GATKSVVCFConstants.EV).toString() : "";
            final Set<String> supportSet = new HashSet<>();
            if (!supportStr.isEmpty()) {
                supportSet.addAll(Arrays.asList(supportStr.split(",")));
            }
            supportMap.put(sample, supportSet);
        }
        return supportMap;
    }

    private Map<String, Integer> getRdCn(final VariantContext variant) {
        final Map<String, Integer> rdCnMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            final Genotype genotype = variant.getGenotype(sample);
            if (genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                rdCnMap.put(sample, Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
            }
        }
        return rdCnMap;
    }
}