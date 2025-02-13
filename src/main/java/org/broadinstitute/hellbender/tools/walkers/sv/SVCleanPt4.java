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

import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.io.BufferedReader;
import java.io.FileReader;
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
        summary = "SClean and format SV VCF",
        oneLineSummary = "Clean and format SV VCF",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVCleanPt4 extends VariantWalker {
    public static final String REVISED_CN_LIST_LONG_NAME = "revised-cn-list";
    public static final String OUTLIERS_LIST_LONG_NAME = "outliers-list";

    @Argument(
            fullName = REVISED_CN_LIST_LONG_NAME,
            doc = "File with variant IDs, sample IDs, and RD_CN values"
    )
    private GATKPath cnReviseList;

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

    private Map<String, Map<String, Integer>> revisedCopyNumbers;
    private Set<String> outlierSamples;

    private double recordStart;
    private double recordEnd;
    private long recordIdx;
    private double maxVF;

    private static final int MIN_LARGE_EVENT_SIZE = 1000;
    private static final int MIN_MULTIALLELIC_EVENT_SIZE = 5000;

    @Override
    public void onTraversalStart() {
        // Read and parse input files
        try {
            revisedCopyNumbers = readRevisedEvents(cnReviseList);
            outlierSamples = new HashSet<>();
            if (outliersListPath != null) {
                outlierSamples = new HashSet<>(Files.readAllLines(outliersListPath.toPath()));
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading input file", e);
        }

        // Parse batch-level metadata
        String cnReviseListFileName = cnReviseList.toPath().getFileName().toString();
        String[] regenoFileNameTokens = cnReviseListFileName.split("\\.");
        String[] batchTokens = regenoFileNameTokens[1].split("_");
        int batchNum = Math.max(Integer.parseInt(batchTokens[0]), 1);
        int totalBatch = Math.max(Integer.parseInt(batchTokens[1]), 1);
        long totalNumVariants = 0;
        String inputVcfPath = getDrivingVariantsFeatureInput().getFeaturePath();
        try (FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(inputVcfPath)) {
            for (VariantContext ignored : dataSource) {
                totalNumVariants++;
            }
        }
        double segments = totalNumVariants / (double) totalBatch;
        recordStart = (batchNum - 1) * segments;
        recordEnd = batchNum * segments;
        maxVF = Math.max((getHeaderForVariants().getGenotypeSamples().size() - outlierSamples.size()) * 0.01, 2);
        recordIdx = 0;

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
        // Exit if outside batch range
        if (recordIdx < recordStart || recordIdx >= recordEnd) {
            recordIdx++;
            return;
        }
        recordIdx++;

        // Initialize data structures
        VariantContextBuilder builder = new VariantContextBuilder(variant);
        List<Genotype> genotypes = variant.getGenotypes();

        // Process variants
        genotypes = processRevisedCn(variant, genotypes);
        processMultiallelic(builder, genotypes);
        genotypes = processLargeDeletions(variant, builder, genotypes);
        genotypes = processLargeDuplications(variant, builder, genotypes);
        genotypes = processRevisedSex(variant, genotypes);

        // Build genotypes
        if (isCalled(builder, genotypes)) {
            builder.genotypes(genotypes);
            vcfWriter.add(builder.make());
        }
    }

    private List<Genotype> processRevisedCn(final VariantContext variant, final List<Genotype> genotypes) {
        final String variantID = variant.getID();
        if (!revisedCopyNumbers.containsKey(variantID)) {
            return genotypes;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        for (Genotype genotype : genotypes) {
            String sampleName = genotype.getSampleName();
            if (revisedCopyNumbers.get(variantID).containsKey(sampleName)) {
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                gb.attribute(GATKSVVCFConstants.RD_CN, revisedCopyNumbers.get(variantID).get(sampleName));
                updatedGenotypes.add(gb.make());
            } else {
                updatedGenotypes.add(genotype);
            }
        }
        return updatedGenotypes;
    }

    private void processMultiallelic(final VariantContextBuilder builder, final List<Genotype> genotypes) {
        int numGtOver2 = 0;
        for (Genotype genotype : genotypes) {
            Integer peGt = genotype.hasExtendedAttribute(GATKSVVCFConstants.PE_GT) ?
                    Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.PE_GT).toString()) : null;
            Integer srGt = genotype.hasExtendedAttribute(GATKSVVCFConstants.SR_GT) ?
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
                Integer peGq = genotype.hasExtendedAttribute(GATKSVVCFConstants.PE_GQ) ?
                        Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.PE_GQ).toString()) : null;
                Integer srGq = genotype.hasExtendedAttribute(GATKSVVCFConstants.SR_GQ) ?
                        Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.SR_GQ).toString()) : null;
                if (peGq != null && srGq != null && peGq >= srGq) {
                    gt = peGt;
                } else {
                    gt = srGt;
                }
            }
            if (gt > 2) {
                numGtOver2++;
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
            for (Genotype genotype : genotypes) {
                if (!outlierSamples.contains(genotype.getSampleName())  && genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                    sampleRdCn.put(genotype.getSampleName(), Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
                }
            }
            if (sampleRdCn.values().stream().filter(value -> value > 3).count() > maxVF) {
                multiallelicFilter = true;
            }
        }

        boolean gt5kbFilter = false;
        List<Integer> allowedAlleleIndices = Arrays.asList(-1, 0, 1);
        if (genotypes.stream().anyMatch(g -> g.getAlleles().stream().anyMatch(a -> !allowedAlleleIndices.contains(variant.getAlleleIndex(a))))) {
            gt5kbFilter = true;
        } else if (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_MULTIALLELIC_EVENT_SIZE && !multiallelicFilter) {
            gt5kbFilter = true;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        if (gt5kbFilter) {
            for (Genotype genotype : genotypes) {
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
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
            for (Genotype genotype : genotypes) {
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
            for (Genotype genotype : genotypes) {
                if (!outlierSamples.contains(genotype.getSampleName()) && genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                    sampleRdCn.put(genotype.getSampleName(), Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
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
        List<Integer> allowedAlleleIndices = Arrays.asList(-1, 0, 1);
        if (genotypes.stream().anyMatch(g -> g.getAlleles().stream().anyMatch(a -> !allowedAlleleIndices.contains(variant.getAlleleIndex(a))))) {
            gt5kbFilter = true;
        } else if (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_MULTIALLELIC_EVENT_SIZE && !multiallelicFilter) {
            gt5kbFilter = true;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        if (gt5kbFilter) {
            for (Genotype genotype : genotypes) {
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
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
            for (Genotype genotype : genotypes) {
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

    private List<Genotype> processRevisedSex(final VariantContext variant, List<Genotype> genotypes) {
        if (!variant.getAttributeAsBoolean(GATKSVVCFConstants.REVISED_EVENT, false)) {
            return genotypes;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        for (Genotype genotype : genotypes) {
            if (Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, 0).toString()) > 0) {
                int newRdCn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()) - 1;
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.attribute(GATKSVVCFConstants.RD_CN, newRdCn);
                if (genotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                    gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, newRdCn);
                }
                updatedGenotypes.add(gb.make());
            } else {
                updatedGenotypes.add(genotype);
            }
        }
        return updatedGenotypes;
    }

    public boolean isCalled(final VariantContextBuilder builder, final List<Genotype> genotypes) {
        for (Genotype genotype : genotypes) {
            if (!isNoCallGt(genotype.getAlleles())) {
                return true;
            }
        }

        if (builder.getAttributes().getOrDefault(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.CNV)) {
            for (Genotype genotype : genotypes) {
                final int cn = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2).toString());
                if (cn != 2) {
                    return true;
                }
            }
        }

        return false;
    }

    private boolean isNoCallGt(List<Allele> alleles) {
        if (alleles.size() == 1 && alleles.get(0).isReference()) return true;
        else if (alleles.size() == 2 && alleles.get(0).isReference() && alleles.get(1).isReference()) return true;
        else if (alleles.size() == 1 && alleles.get(0).isNoCall()) return true;
        return false;
    }

    private Map<String, Map<String, Integer>> readRevisedEvents(final GATKPath filePath) {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toPath().toFile()))) {
            final Map<String, Map<String, Integer>> result = new HashMap<>();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields.length < 3) continue;

                String variantId = fields[0];
                String sampleId = fields[1];
                int rdCn = Integer.parseInt(fields[2]);

                result.computeIfAbsent(variantId, k -> new HashMap<>()).put(sampleId, rdCn);
            }
            return result;
        } catch (IOException e) {
            throw new RuntimeException("Error reading input file", e);
        }
    }
}
