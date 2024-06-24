package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.mutable.MutableLong;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.broadinstitute.hellbender.engine.*;

import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Given a VCF of known variants from multiple samples, calculate how much each sample
 * contributes to a pooled BAM.  We first estimate each sample's mixing fraction as the average alt fraction
 * in the bam of singleton hets belonging to that sample, then normalize these initial estimates to sum to one.
 *
 * <p>
 * In the CRSP sensitivity validation, we have a bam derived from a pool of 5, 10, or 20 samples
 * and a vcf of all known variants in those samples.  The pooled bam is a simulated tumor and
 * the vcf of individual variants is our truth data.  We annotate the truth data with the estimated
 * allele fractions of each variant in the pooled bam in order to bin our results by "tumor" allele fraction.
 * To estimate allele fractions we need to know the sample mixing fractions as an intermediate step.
 * </p>
 *
 * <p>
 * One could use this tool for generating precise truth data any time multiple samples with known ground
 * truth are combined to simulate tumors or any situation in which low allele fractions occur.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk --java-options "-Xmx4g" CalculateMixingFractions \
 *   -V input.vcf \
 *   -O output.table
 * </pre>
 *
 * Created by David Benjamin on 1/30/17.
 */
@CommandLineProgramProperties(
        summary = "Calculate proportions of different samples in a pooled bam by first estimating each sample's" +
                " mixing fraction as the average alt fraction in the bam of singleton hets belonging to that sample, then" +
                " normalizing these initial estimates to sum to one.",
        oneLineSummary = "(Internal) Calculate proportions of different samples in a pooled bam",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class CalculateMixingFractions extends VariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table of samplemixing fractions",
            optional=false)
    private final File outputFile = null;

    private final Map<String, AltAndTotalReadCounts> sampleCounts = new HashMap<>();

    @Override
    public void onTraversalStart() {
        getHeaderForVariants().getGenotypeSamples().forEach(s -> sampleCounts.put(s, new AltAndTotalReadCounts()));
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        if (!isBiallelicSingletonHetSnp(vc)) {
            return;
        }

        final Optional<String> variantSample = StreamSupport.stream(vc.getGenotypes().spliterator(), false)
                .filter(genotype -> genotype.isHet())
                .map(genotype -> genotype.getSampleName())
                .findFirst();

        if (!variantSample.isPresent()) {
            return;
        }

        final List<GATKRead> reads = new ArrayList<>();
        final List<Integer> offsets = new ArrayList<>();

        for (final GATKRead read : readsContext) {
            if (read.failsVendorQualityCheck()) {
                continue;
            }
            final AlignmentStateMachine asm = new AlignmentStateMachine(read);
            while ( asm.stepForwardOnGenome() != null && asm.getGenomePosition() < vc.getStart()) {

            }

            if (asm.getGenomePosition() == vc.getStart()) {
                reads.add(read);
                offsets.add(asm.getReadOffset());
            }
        }

        final ReadPileup pileup = new ReadPileup(vc, reads, offsets);

        final byte altBase = vc.getAlternateAllele(0).getBases()[0];
        final long altCount = StreamSupport.stream(pileup.spliterator(), false)
                .filter(pe -> pe.getBase() == altBase)
                .count();
        final long totalCount = pileup.size();
        sampleCounts.get(variantSample.get()).addCounts(altCount, totalCount);
    }

    @Override
    public Object onTraversalSuccess() {
        final double mixingFractionNormalizer = sampleCounts.values().stream()
                .mapToDouble(AltAndTotalReadCounts::getAltFraction).sum();

        final List<MixingFraction> mixingFractions = sampleCounts.entrySet().stream()
                .map(e -> new MixingFraction(e.getKey(), e.getValue().getAltFraction() / mixingFractionNormalizer))
                .collect(Collectors.toList());

        MixingFraction.writeMixingFractions(mixingFractions, outputFile);

        return "SUCCESS";
    }

    private class AltAndTotalReadCounts {
        private final MutableLong altCount = new MutableLong(0);
        private final MutableLong totalCount = new MutableLong(0);

        public void addCounts(final long numAlts, final long numTotal) {
            altCount.add(numAlts);
            totalCount.add(numTotal);
        }

        public double getAltFraction() { return altCount.doubleValue() / totalCount.doubleValue(); }
    }

    private boolean isBiallelicSingletonHetSnp(final VariantContext vc) {
        return vc.isBiallelic() && vc.isSNP()
                && ( (vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)
                && getArrayAttribute(vc, VCFConstants.ALLELE_COUNT_KEY)[0] == 1)
                || vc.getGenotypes().stream().filter(genotype -> genotype.isHet()).count() == 1);
    }

    private static double[] getArrayAttribute(final VariantContext vc, final String attribute) {
        return VariantContextGetters.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }
}
