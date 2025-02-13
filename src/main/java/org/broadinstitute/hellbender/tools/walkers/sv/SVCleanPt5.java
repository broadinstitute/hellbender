package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;

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
public class SVCleanPt5 extends MultiplePassVariantWalker {
    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;

    private final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();
    private final Set<String> filteredVariantIds = new HashSet<>();

    @Override
    protected int numberOfPasses() { return 2; }

    @Override
    public void onTraversalStart() {
        // Remove unnecessary header lines
        final Set<VCFHeaderLine> newHeaderLines = new HashSet<>();
        final VCFHeader header = getHeaderForVariants();
        for (final VCFHeaderLine line : header.getMetaDataInInputOrder()) {
            if (line instanceof VCFInfoHeaderLine) {
                if (GATKSVVCFConstants.FILTER_VCF_INFO_LINES.contains(((VCFInfoHeaderLine) line).getID())) {
                    continue;
                }
            } else if (line instanceof VCFFormatHeaderLine) {
                if (((VCFFormatHeaderLine) line).getID().equals(GATKSVVCFConstants.EV)) {
                    continue;
                }
            } else if (line instanceof VCFAltHeaderLine) {
                if (((VCFAltHeaderLine) line).getID().equals(GATKSVVCFConstants.UNR)) {
                    continue;
                }
            }
            if (GATKSVVCFConstants.FILTER_VCF_LINES.stream().anyMatch(line.toString()::contains)) {
                continue;
            }
            newHeaderLines.add(line);
        }

        // Add new header lines
        VCFHeader newHeader = new VCFHeader(newHeaderLines, header.getGenotypeSamples());
        newHeader.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.EV, 1, VCFHeaderLineType.String, "Classes of evidence supporting final genotype"));

        // Write header
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(newHeader);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext, int n) {
        switch (n) {
            case 0:
                firstPassApply(variant);
                break;
            case 1:
                secondPassApply(variant);
                break;
        }
    }

    @Override
    protected void afterNthPass(int n) {}

    public void firstPassApply(final VariantContext variant) {
        if (!variant.getFilters().contains(GATKSVVCFConstants.MULTIALLELIC)) {
            return;
        }

        overlappingVariantsBuffer.removeIf(vc -> !vc.getContig().equals(variant.getContig())
                || (vc.getStart() + vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)) < variant.getStart());
        for (VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(bufferedVariant, variant)) {
                processVariantPair(bufferedVariant, variant);
            }
        }
        overlappingVariantsBuffer.add(variant);
    }

    public void secondPassApply(final VariantContext variant) {
        if (filteredVariantIds.contains(variant.getID())) {
            return;
        }

        VariantContextBuilder builder = new VariantContextBuilder(variant);
        processSvType(variant, builder);
        cleanseInfoFields(builder);
        vcfWriter.add(builder.make());
    }

    private void processVariantPair(VariantContext v1, VariantContext v2) {
        // Determine larger variant
        VariantContext largerVariant = v1;
        VariantContext smallerVariant = v2;
        int length1 = v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        int length2 = v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        int smallerLength = length2;

        // Swap variants if necessary
        if (length2 > length1) {
            largerVariant = v2;
            smallerVariant = v1;
            smallerLength = length1;
        }

        // Calculate overlap
        int minEnd = Math.min(
                largerVariant.getStart() + largerVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0),
                smallerVariant.getStart() + smallerVariant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)
        );
        int maxStart = Math.max(largerVariant.getStart(), smallerVariant.getStart());
        int overlapLength = minEnd - maxStart + 1;
        if (overlapLength <= 0) {
            return;
        }

        // Filter variant based on conditions
        double coverage = (double) overlapLength / smallerLength;
        if (coverage > 0.5 && !filteredVariantIds.contains(largerVariant.getID())) {
            filteredVariantIds.add(smallerVariant.getID());
        }
    }

    private void processSvType(final VariantContext variant, final VariantContextBuilder builder) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
        boolean hasMobileElement = variant.getAlleles().stream()
                .map(GATKSVVariantContextUtils::getSymbolicAlleleSymbols)
                .flatMap(Arrays::stream)
                .anyMatch(symbol -> symbol.equals(GATKSVVCFConstants.ME));
        if (svType == null || hasMobileElement) {
            return;
        }

        final Allele refAllele = variant.getReference();
        final Allele altAllele = Allele.create("<" + svType + ">", false);
        List<Allele> newAlleles = Arrays.asList(refAllele, altAllele);

        List<Genotype> genotypes = variant.getGenotypes();
        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        for (Genotype genotype : genotypes) {
            GenotypeBuilder gb = new GenotypeBuilder(genotype);
            long altCount = genotype.getAlleles().stream().filter(allele -> allele.isCalled() && !allele.isReference()).count();
            if (altCount == 1) { // Heterozygous (0/1)
                gb.alleles(Arrays.asList(refAllele, altAllele));
            } else if (altCount == 2) { // Homozygous Alternate (1/1)
                gb.alleles(Arrays.asList(altAllele, altAllele));
            }
            updatedGenotypes.add(gb.make());
        }

        builder.alleles(newAlleles);
        builder.genotypes(updatedGenotypes);
    }

    private void cleanseInfoFields(final VariantContextBuilder builder) {
        Map<String, Object> attributes = builder.getAttributes();
        for (String field : GATKSVVCFConstants.FILTER_VCF_INFO_LINES) {
            if (attributes.containsKey(field)) {
                builder.rmAttribute(field);
            }
        }
    }

    private boolean overlaps(final VariantContext v1, final VariantContext v2) {
        return v1.getContig().equals(v2.getContig())
                && v1.getStart() <= (v2.getStart() + v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0))
                && v2.getStart() <= (v1.getStart() + v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
    }
}
