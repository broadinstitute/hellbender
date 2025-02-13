package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.io.IOException;
import java.nio.file.Files;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;
import java.util.HashMap;

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
public class SVReviseLargeCnvs extends VariantWalker {
    public static final String OUTLIERS_LIST_LONG_NAME = "outliers-list";

    @Argument(
            fullName = OUTLIERS_LIST_LONG_NAME,
            doc = "File with outlier samples",
            optional = true
    )
    private GATKPath outliersListPath;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;

    private Set<String> outlierSamples;

    private double maxVF;

    private static final int MIN_LARGE_EVENT_SIZE = 1000;
    private static final int MIN_MULTIALLELIC_EVENT_SIZE = 5000;

    @Override
    public void onTraversalStart() {
        // Read and parse input files
        try {
            outlierSamples = new HashSet<>();
            if (outliersListPath != null) {
                outlierSamples = new HashSet<>(Files.readAllLines(outliersListPath.toPath()));
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading input file", e);
        }

        // Populate maxVf based on sample information
        maxVF = Math.max((getHeaderForVariants().getGenotypeSamples().size() - outlierSamples.size()) * 0.01, 2);

        // Filter specific header lines
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PESR_GT_OVERDISPERSION, 0, VCFHeaderLineType.Flag, "High PESR dispersion count"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "Predicted copy state"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, 1, VCFHeaderLineType.Integer, "Read-depth genotype quality"));
        header.addMetaDataLine(new VCFFilterHeaderLine(GATKSVVCFConstants.MULTIALLELIC, "Multiallelic site"));

        // Write header
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(header);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Initialize data structures
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        List<Genotype> genotypes = variant.getGenotypes();

        // Process variants
        processMultiallelic(builder, genotypes);
        genotypes = processLargeDeletions(variant, builder, genotypes);
        genotypes = processLargeDuplications(variant, builder, genotypes);

        // Build genotypes
        if (isCalled(builder, genotypes)) {
            builder.genotypes(genotypes);
            vcfWriter.add(builder.make());
        }
    }

    private void processMultiallelic(final VariantContextBuilder builder, final List<Genotype> genotypes) {
        int numGtOver2 = 0;
        for (final Genotype genotype : genotypes) {
            final Integer peGt = genotype.hasExtendedAttribute(GATKSVVCFConstants.PE_GT) ?
                    Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.PE_GT).toString()) : null;
            final Integer srGt = genotype.hasExtendedAttribute(GATKSVVCFConstants.SR_GT) ?
                    Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.SR_GT).toString()) : null;
            int gt;
            if (peGt == null) {
                continue;
            } else if (srGt == null) {
                gt = peGt;
            } else if (peGt > 0 && srGt == 0) {
                gt = peGt;
            } else if (peGt == 0) {
                gt = srGt;
            } else {
                final Integer peGq = genotype.hasExtendedAttribute(GATKSVVCFConstants.PE_GQ) ?
                        Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.PE_GQ).toString()) : null;
                final Integer srGq = genotype.hasExtendedAttribute(GATKSVVCFConstants.SR_GQ) ?
                        Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.SR_GQ).toString()) : null;
                if (peGq != null && srGq != null && peGq >= srGq) {
                    gt = peGt;
                } else {
                    gt = srGt;
                }
            }
            if (gt > 2) {
                numGtOver2 += 1;
            }
        }
        if (numGtOver2 > maxVF) {
            builder.attribute(GATKSVVCFConstants.PESR_GT_OVERDISPERSION, true);
        }
    }

    private List<Genotype> processLargeDeletions(final VariantContext variant, final VariantContextBuilder builder, List<Genotype> genotypes) {
        if (!variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL)) {
            return genotypes;
        }

        boolean multiallelicFilter = false;
        if (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_LARGE_EVENT_SIZE) {
            Map<String, Integer> sampleRdCn = new HashMap<>();
            for (final Genotype genotype : genotypes) {
                final String sample = genotype.getSampleName();
                if (!outlierSamples.contains(sample)  && genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                    sampleRdCn.put(sample, Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
                }
            }
            if (sampleRdCn.values().stream().filter(value -> value > 3).count() > maxVF) {
                multiallelicFilter = true;
            }
        }

        boolean gt5kbFilter = false;
        final List<Integer> allowedAlleleIndices = Arrays.asList(-1, 0, 1);
        if (genotypes.stream().anyMatch(g -> g.getAlleles().stream().anyMatch(a -> !allowedAlleleIndices.contains(variant.getAlleleIndex(a))))) {
            gt5kbFilter = true;
        } else if (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_MULTIALLELIC_EVENT_SIZE && !multiallelicFilter) {
            gt5kbFilter = true;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        if (gt5kbFilter) {
            for (final Genotype genotype : genotypes) {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                if (!genotype.isNoCall()) {
                    if (genotype.hasGQ() && Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, 0).toString()) >= 2) {
                        gb.alleles(Arrays.asList(variant.getReference(), variant.getReference()));
                    } else if (genotype.hasGQ() && Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, 0).toString()) == 1) {
                        gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                    } else if (genotype.hasGQ()) {
                        gb.alleles(Arrays.asList(variant.getAlternateAllele(0), variant.getAlternateAllele(0)));
                    }
                }
                updatedGenotypes.add(gb.make());
            }
            genotypes = updatedGenotypes;
        }

        updatedGenotypes = new ArrayList<>(genotypes.size());
        if (multiallelicFilter) {
            for (final Genotype genotype : genotypes) {
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.noGQ();
                gb.alleles(Arrays.asList(Allele.NO_CALL));
                gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN));
                gb.attribute(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.RD_GQ));
                updatedGenotypes.add(gb.make());
            }
            genotypes = updatedGenotypes;

            builder.filter(GATKSVVCFConstants.MULTIALLELIC);
            builder.attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CNV);
            builder.alleles(Arrays.asList(variant.getReference(), Allele.create("<" + GATKSVVCFConstants.CNV + ">", false)));
        }

        return genotypes;
    }

    private List<Genotype> processLargeDuplications(final VariantContext variant, final VariantContextBuilder builder, List<Genotype> genotypes) {
        if (!variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP)) {
            return genotypes;
        }

        boolean multiallelicFilter = false;
        if (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_LARGE_EVENT_SIZE) {
            Map<String, Integer> sampleRdCn = new HashMap<>();
            for (final Genotype genotype : genotypes) {
                final String sample = genotype.getSampleName();
                if (!outlierSamples.contains(sample) && genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                    sampleRdCn.put(sample, Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
                }
            }
            if (sampleRdCn.values().stream().filter(value -> value > 4).count() > maxVF) {
                multiallelicFilter = true;
            }
            if (sampleRdCn.values().stream().filter(value -> (value < 1 || value > 4)).count() > 4) {
                if (sampleRdCn.values().stream().filter(value -> (value < 1 || value > 4)).distinct().count() > maxVF) {
                    multiallelicFilter = true;
                }
            }
        }

        boolean gt5kbFilter = false;
        final List<Integer> allowedAlleleIndices = Arrays.asList(-1, 0, 1);
        if (genotypes.stream().anyMatch(g -> g.getAlleles().stream().anyMatch(a -> !allowedAlleleIndices.contains(variant.getAlleleIndex(a))))) {
            gt5kbFilter = true;
        } else if (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_MULTIALLELIC_EVENT_SIZE && !multiallelicFilter) {
            gt5kbFilter = true;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        if (gt5kbFilter) {
            for (final Genotype genotype : genotypes) {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                if (!genotype.isNoCall()) {
                    if (genotype.hasGQ() && Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, 3).toString()) <= 2) {
                        gb.alleles(Arrays.asList(variant.getReference(), variant.getReference()));
                    } else if (genotype.hasGQ() && Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, 0).toString()) == 3) {
                        gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                    } else if (genotype.hasGQ()) {
                        gb.alleles(Arrays.asList(variant.getAlternateAllele(0), variant.getAlternateAllele(0)));
                    }
                }
                updatedGenotypes.add(gb.make());
            }
            genotypes = updatedGenotypes;
        }

        updatedGenotypes = new ArrayList<>(genotypes.size());
        if (multiallelicFilter) {
            for (final Genotype genotype : genotypes) {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.noGQ();
                gb.alleles(Arrays.asList(Allele.NO_CALL));
                gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN));
                gb.attribute(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.RD_GQ));
                updatedGenotypes.add(gb.make());
            }
            genotypes = updatedGenotypes;

            builder.filter(GATKSVVCFConstants.MULTIALLELIC);
            builder.attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CNV);
            builder.alleles(Arrays.asList(variant.getReference(), Allele.create("<" + GATKSVVCFConstants.CNV + ">", false)));
        }

        return genotypes;
    }

    public boolean isCalled(final VariantContextBuilder builder, final List<Genotype> genotypes) {
        for (final Genotype genotype : genotypes) {
            if (!isNoCallGt(genotype.getAlleles())) {
                return true;
            }
        }

        if (builder.getAttributes().getOrDefault(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.CNV)) {
            for (final Genotype genotype : genotypes) {
                if (Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2).toString()) != 2) {
                    return true;
                }
            }
        }

        return false;
    }

    private boolean isNoCallGt(final List<Allele> alleles) {
        if (alleles.size() == 1 && alleles.get(0).isReference()) return true;
        else if (alleles.size() == 2 && alleles.get(0).isReference() && alleles.get(1).isReference()) return true;
        else if (alleles.size() == 1 && alleles.get(0).isNoCall()) return true;
        return false;
    }
}
