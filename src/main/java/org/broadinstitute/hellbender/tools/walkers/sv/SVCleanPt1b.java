package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;
import java.util.HashMap;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

/**
 * Completes an initial series of cleaning steps for a VCF produced by the GATK-SV pipeline.
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>
 *         VCF containing structural variant (SV) records from the GATK-SV pipeline.
 *     </li>
 *     <li>
 *         TODO
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>
 *         Cleansed VCF.
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
public class SVCleanPt1b extends MultiplePassVariantWalker {
    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;

    private final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();
    final private Map<String, Map<String, Pair<String, String>>> revisedEventsAll = new HashMap<>();
    final private Map<String, Set<String>> revisedEventsFiltered = new HashMap<>();
    final private Map<String, Map<String, Integer>> revisedRdCn = new HashMap<>();

    private static final int MIN_VARIANT_SIZE_CNV = 1000;
    private static final int MIN_VARIANT_SIZE = 5000;

    @Override
    protected int numberOfPasses() {
        return 3;
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext, final int n) {
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
            default:
                throw new IllegalArgumentException("Invalid pass number: " + n);
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        switch (n) {
            case 0:
                processCollectedVariants();
                break;
            case 1:
                vcfWriter = createVCFWriter(outputVcf);
                final VCFHeader header = getHeaderForVariants();
                header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MULTI_CNV, 0, VCFHeaderLineType.Flag, "Variant is a multiallelic CNV"));
                vcfWriter.writeHeader(header);
                break;
        }
    }

    public void firstPassApply(final VariantContext variant) {
        if (!isDelDup(variant) || !isLarge(variant, MIN_VARIANT_SIZE)) {
            return;
        }

        // Process overlaps with variants in the buffer
        overlappingVariantsBuffer.removeIf(vc -> !vc.getContig().equals(variant.getContig()) || vc.getEnd() < variant.getStart());
        for (VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(bufferedVariant, variant)) {
                processOverlap(bufferedVariant, variant);
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
            final String rdCn = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN);
            variantRdCn.put(sampleName, Integer.parseInt(rdCn));
        }
        revisedRdCn.put(variantId, variantRdCn);
    }

    public void thirdPassApply(final VariantContext variant) {
        VariantContextBuilder builder = new VariantContextBuilder(variant);
        if (revisedEventsAll.containsKey(variant.getID())) {
            processVariant(builder, variant);
        }
        if (isDelDup(variant) && isLarge(variant, MIN_VARIANT_SIZE_CNV)) {
            processCnvs(builder, variant);
        }
        vcfWriter.add(builder.make());
    }

    private void processOverlap(final VariantContext v1, final VariantContext v2) {
        // Get overlap data
        VariantContext wider;
        VariantContext narrower;
        if (v1.getLengthOnReference() > v2.getLengthOnReference()) {
            wider = v1;
            narrower = v2;
        } else if (v2.getLengthOnReference() > v1.getLengthOnReference()) {
            wider = v2;
            narrower = v1;
        } else {
            return;
        }
        String widerID = wider.getID();
        String narrowerID = narrower.getID();

        // Skip processing if same variant ID, SV type or samples
        String widerSvType = wider.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        String narrowerSvType = narrower.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        Set<String> widerSamples = getNonReferenceSamples(wider);
        Set<String> narrowerSamples = getNonReferenceSamples(narrower);
        if (widerID.equals(narrowerID) || widerSvType.equals(narrowerSvType) || widerSamples.equals(narrowerSamples)) {
            return;
        }

        // Get samples present in wider but not in narrower
        Set<String> nonCommonSamples = new HashSet<>(widerSamples);
        nonCommonSamples.removeAll(narrowerSamples);
        if (nonCommonSamples.isEmpty()) {
            return;
        }

        // Revise variant if coverage exceeds threshold
        double coverage = getCoverage(wider, narrower);
        if (coverage >= 0.5) {
            for (String sample : nonCommonSamples) {
                revisedEventsAll.computeIfAbsent(narrowerID, k -> new HashMap<>())
                        .put(sample, new ImmutablePair<>(widerID, widerSvType));
            }
        }
    }

    private void processCollectedVariants() {
        for (final Map.Entry<String, Map<String, Pair<String, String>>> entry : revisedEventsAll.entrySet()) {
            for (final Map.Entry<String, Pair<String, String>> innerEntry : entry.getValue().entrySet()) {
                // Identifies variant-sample pairs we need RD_CN values for to improve speed
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

    private void processVariant(final VariantContextBuilder builder, final VariantContext variant) {
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
                final int currentRdCn = revisedRdCn.get(variantId).getOrDefault(sample, 0);
                final int widerRdCn = revisedRdCn.getOrDefault(widerVariantId, new HashMap<>()).getOrDefault(sample, 0);
                if (!revisedEventsFiltered.getOrDefault(widerVariantId, new HashSet<>()).contains(sample)) {
                    System.err.println(sample + " " + widerVariantId);
                }

                int newVal = -1;
                if (widerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && currentRdCn == 2 && widerRdCn == 3) {
                    newVal = 1;
                } else if (widerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && currentRdCn == 2 && widerRdCn == 1) {
                    newVal = 3;
                }

                if (newVal != -1) {
                    System.out.println(variantId);
                    final GenotypeBuilder gb = new GenotypeBuilder(oldGenotype);
                    gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                    gb.GQ(Integer.parseInt((String) oldGenotype.getExtendedAttribute(GATKSVVCFConstants.RD_GQ)));
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

    private void processCnvs(final VariantContextBuilder builder, final VariantContext variant) {
        final boolean isDel = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        for (String sample : variant.getSampleNamesOrderedByName()) {
            final Genotype genotype = variant.getGenotype(sample);
            final String rdCnString = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN);
            final int rdCn = Integer.parseInt(rdCnString);
            if ((isDel && rdCn > 3) || (!isDel && (rdCn < 1 || rdCn > 4))) {
                builder.attribute(GATKSVVCFConstants.MULTI_CNV, true);
                break;
            }
        }
    }

    private boolean isDelDup(final VariantContext variant) {
        String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        return svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP);
    }

    private boolean isLarge(final VariantContext variant, final int minSize) {
        int variantLength = Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
        return variantLength >= minSize;
    }

    private boolean overlaps(final VariantContext v1, final VariantContext v2) {
        return v1.getContig().equals(v2.getContig()) && v1.getStart() <= v2.getEnd() && v2.getStart() <= v1.getEnd();
    }

    private Set<String> getNonReferenceSamples(final VariantContext variant) {
        Set<String> samples = new HashSet<>();
        for (String sampleName : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sampleName);
            if (genotype.isCalled() && !genotype.isHomRef()) {
                samples.add(sampleName);
            }
        }
        return samples;
    }

    private double getCoverage(final VariantContext wider, final VariantContext narrower) {
        int nStart = narrower.getStart();
        int nStop = narrower.getEnd();
        int wStart = wider.getStart();
        int wStop = wider.getEnd();

        if (wStart <= nStop && nStart <= wStop) {
            int intersectionSize = Math.min(nStop, wStop) - Math.max(nStart, wStart) + 1;
            return (double) intersectionSize / (nStop - nStart + 1);
        }
        return 0.0;
    }
}
