package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.stream.Collectors;

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
public final class SVCleanPt1a extends VariantWalker {
    public static final String CHRX_LONG_NAME = "chr-X";
    public static final String CHRY_LONG_NAME = "chr-Y";
    public static final String FAIL_LIST_LONG_NAME = "fail-list";
    public static final String PASS_LIST_LONG_NAME = "pass-list";
    public static final String OUTPUT_SAMPLES_LIST_LONG_NAME = "output-samples-list";

    @Argument(
            fullName = CHRX_LONG_NAME,
            doc = "chrX column name",
            optional = true
    )
    private final String chrX = "chrX";

    @Argument(
            fullName = CHRY_LONG_NAME,
            doc = "chrY column name",
            optional = true
    )
    private final String chrY = "chrY";

    @Argument(
            fullName = FAIL_LIST_LONG_NAME,
            doc = "File with variants failing the background test"
    )
    private GATKPath failList;

    @Argument(
            fullName = PASS_LIST_LONG_NAME,
            doc = "File with variants passing both sides"
    )
    private GATKPath passList;

    @Argument(
            fullName = OUTPUT_SAMPLES_LIST_LONG_NAME,
            doc = "Output file with samples"
    )
    private GATKPath outputSamplesList;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;
    private BufferedWriter samplesWriter = null;

    private Set<String> failSet;
    private Set<String> passSet;

    private static final int MIN_ALLOSOME_EVENT_SIZE = 5000;

    @Override
    public void onTraversalStart() {
        // Read supporting files
        failSet = readLastColumn(failList);
        passSet = readLastColumn(passList);

        // Add new header lines
        VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(new VCFFilterHeaderLine(GATKSVVCFConstants.UNRESOLVED, "Variant is unresolved"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HIGH_SR_BACKGROUND, 0, VCFHeaderLineType.Flag, "High number of SR splits in background samples indicating messy region"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BOTHSIDES_SUPPORT, 0, VCFHeaderLineType.Flag, "Variant has read-level support for both sides of breakpoint"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.REVISED_EVENT, 0, VCFHeaderLineType.Flag, "Variant has been revised due to a copy number mismatch"));

        // Write header
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(header);

        // Write samples list
        try {
            samplesWriter = new BufferedWriter(new FileWriter(outputSamplesList.toPath().toFile()));
            for (String sample : header.getGenotypeSamples()) {
                samplesWriter.write(sample);
                samplesWriter.newLine();
            }
            samplesWriter.flush();
        } catch (IOException e) {
            throw new RuntimeException("Can't create output file", e);
        }
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Create core data
        VariantContextBuilder builder = new VariantContextBuilder(variant);
        List<Genotype> genotypes = variant.getGenotypes();

        // Process variant
        genotypes = processEV(genotypes);
        processVarGQ(variant, builder);
        processMultiallelic(variant, builder);
        processUnresolved(variant, builder);
        processNoisyEvents(variant, builder);
        processBothsidesSupportEvents(variant, builder);
        genotypes = processAllosomes(variant, builder, genotypes);

        builder.genotypes(genotypes);
        vcfWriter.add(builder.make());
    }

    private List<Genotype> processEV(final List<Genotype> genotypes) {
        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        for (Genotype genotype : genotypes) {
            if (genotype.hasExtendedAttribute(GATKSVVCFConstants.EV)) {
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
                String evAttribute = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.EV);
                final int evIndex = Integer.parseInt(evAttribute);
                if (evIndex >= 0 && evIndex < GATKSVVCFConstants.EV_VALUES.size()) {
                    gb.attribute(GATKSVVCFConstants.EV, GATKSVVCFConstants.EV_VALUES.get(evIndex));
                }
                updatedGenotypes.add(gb.make());
            } else {
                updatedGenotypes.add(genotype);
            }
        }
        return updatedGenotypes;
    }

    private void processVarGQ(final VariantContext variant, final VariantContextBuilder builder) {
        if (variant.hasAttribute(GATKSVVCFConstants.VAR_GQ)) {
            final double varGQ = variant.getAttributeAsDouble(GATKSVVCFConstants.VAR_GQ, 0);
            builder.rmAttribute(GATKSVVCFConstants.VAR_GQ);
            builder.log10PError(varGQ / -10.0);
        }
    }

    private void processMultiallelic(final VariantContext variant, final VariantContextBuilder builder) {
        if (variant.hasAttribute(GATKSVVCFConstants.MULTIALLELIC)) {
            builder.rmAttribute(GATKSVVCFConstants.MULTIALLELIC);
        }
    }

    private void processUnresolved(final VariantContext variant, final VariantContextBuilder builder) {
        if (variant.hasAttribute(GATKSVVCFConstants.UNRESOLVED)) {
            builder.rmAttribute(GATKSVVCFConstants.UNRESOLVED);
            builder.filter(GATKSVVCFConstants.UNRESOLVED);
        }
    }

    private void processNoisyEvents(final VariantContext variant, final VariantContextBuilder builder) {
        if (failSet.contains(variant.getID())) {
            builder.attribute(GATKSVVCFConstants.HIGH_SR_BACKGROUND, true);
        }
    }

    private void processBothsidesSupportEvents(final VariantContext variant, final VariantContextBuilder builder) {
        if (passSet.contains(variant.getID())) {
            builder.attribute(GATKSVVCFConstants.BOTHSIDES_SUPPORT, true);
        }
    }

    private List<Genotype> processAllosomes(final VariantContext variant, final VariantContextBuilder builder, final List<Genotype> genotypes) {
        final String chromosome = variant.getContig();
        if (!chromosome.equals(chrX) && !chromosome.equals(chrY)) {
            return genotypes;
        }

        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        for (Genotype genotype : genotypes) {
            final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP)) &&
                    (variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0) >= MIN_ALLOSOME_EVENT_SIZE)) {
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
                final boolean isY = chromosome.equals(chrY);
                final int sex = (int) genotype.getExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT);
                if (sex == 1 && isRevisableEvent(variant, isY, sex)) { // Male
                    builder.attribute(GATKSVVCFConstants.REVISED_EVENT, true);
                    adjustMaleGenotype(genotype, gb, svType);
                } else if (sex == 2 && isY) { // Female
                    gb.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
                } else if (sex == 0) { // Unknown
                    gb.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
                }
                updatedGenotypes.add(gb.make());
            }
            else {
                updatedGenotypes.add(genotype);
            }
        }
        return updatedGenotypes;
    }

    private void adjustMaleGenotype(final Genotype genotype, final GenotypeBuilder gb, final String svType) {
        if (genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
            final int rdCN = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            gb.attribute(GATKSVVCFConstants.RD_CN, rdCN + 1);

            final Allele refAllele = genotype.getAllele(0);
            final Allele altAllele = genotype.getAllele(1);
            if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL)) {
                if (rdCN >= 1) gb.alleles(Arrays.asList(refAllele, refAllele));
                else if (rdCN == 0) gb.alleles(Arrays.asList(refAllele, altAllele));
            } else if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP)) {
                if (rdCN <= 1) gb.alleles(Arrays.asList(refAllele, refAllele));
                else if (rdCN == 2) gb.alleles(Arrays.asList(refAllele, altAllele));
                else gb.alleles(Arrays.asList(altAllele, altAllele));
            }
        }
    }

    private boolean isRevisableEvent(final VariantContext variant, final boolean isY, final int sex) {
        final List<Genotype> genotypes = variant.getGenotypes();
        final int[] maleCounts = new int[4];
        final int[] femaleCounts = new int[4];
        for (final Genotype genotype : genotypes) {
            final int rdCN = (int) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, -1);
            final int rdCNVal = Math.min(rdCN, 3);
            if (rdCNVal == -1) continue;

            if (sex == 1) {
                maleCounts[rdCNVal]++;
            } else if (sex == 2) {
                femaleCounts[rdCNVal]++;
            }
        }

        final int maleMedian = calcMedianDistribution(maleCounts);
        final int femaleMedian = calcMedianDistribution(femaleCounts);
        return maleMedian == 2 && (isY ? femaleMedian == 0 : femaleMedian == 4);
    }

    private int calcMedianDistribution(final int[] counts) {
        final int total = Arrays.stream(counts).sum();
        if (total == 0) return -1;

        final int target = total / 2;
        int runningTotal = 0;
        for (int i = 0; i < 4; i++) {
            runningTotal += counts[i];
            if (runningTotal == target) {
                return i * 2 + 1;
            } else if (runningTotal > target) {
                return i * 2;
            }
        }
        throw new RuntimeException("Error calculating median");
    }

    private Set<String> readLastColumn(final GATKPath filePath) {
        try {
            final Path path = filePath.toPath();
            final TableReader<String> reader = TableUtils.reader(path, (columns, exceptionFactory) ->
                    (dataline) -> dataline.get(columns.columnCount() - 1)
            );

            Set<String> result = reader.stream().collect(Collectors.toSet());
            reader.close();
            return result;
        } catch (IOException e) {
            throw new RuntimeException("Error reading variant list file: " + filePath, e);
        }
    }
}