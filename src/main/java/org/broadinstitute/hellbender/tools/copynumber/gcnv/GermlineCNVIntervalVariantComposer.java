package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntervalCopyNumberGenotypingData;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of
 * {@link GermlineCNVCaller} calls into genotyped intervals.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVIntervalVariantComposer extends GermlineCNVVariantComposer<IntervalCopyNumberGenotypingData> {

    /* VCF FORMAT header keys */

    /**
     * Copy number maximum a posteriori value
     */
    static final String CN = "CN";

    /**
     * Copy number log posterior (in Phred-scale)
     */
    static final String CNLP = "CNLP";

    /**
     * Genotype call quality
     */
    static final String CNQ = "CNQ";

    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;

    /**
     * Constructor for {@link PostprocessGermlineCNVCalls} Postprocessor
     *
     * @param outputWriter variant context writer
     * @param sampleName sample name
     * @param refAutosomalCopyNumberState ref copy-number state on autosomal contigs
     * @param allosomalContigSet set of allosomal contigs (ref copy-number allele be chosen according to
     *                           given contig baseline copy-number states)
     */
    public GermlineCNVIntervalVariantComposer(final VariantContextWriter outputWriter,
                                              final String sampleName,
                                              final IntegerCopyNumberState refAutosomalCopyNumberState,
                                              final Set<String> allosomalContigSet) {
        super(outputWriter, sampleName);
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
    }

    @Override
    public void composeVariantContextHeader(final SAMSequenceDictionary sequenceDictionary,
                                            final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        result.setSequenceDictionary(sequenceDictionary);

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(result::addMetaDataLine);

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.String, "Genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori value"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNLP, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Integer, "Copy number log posterior (in Phred-scale) rounded down"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNQ, 1,
                VCFHeaderLineType.Integer, "Genotype call quality as the difference between" +
                " the best and second best phred-scaled log posterior scores"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));
        outputWriter.writeHeader(result);
    }

    @VisibleForTesting
    VariantContext composeVariantContext(final IntervalCopyNumberGenotypingData intervalCopyNumberGenotypingData) {
        final CopyNumberPosteriorDistribution copyNumberPosteriorDistribution =
                intervalCopyNumberGenotypingData.getCopyNumberPosteriorDistribution();
        final List<Integer> copyNumberPLVector = calculateCopyNumberPLVector(copyNumberPosteriorDistribution);
        final int copyNumberMAP = calculateMAPCopyNumberState(copyNumberPosteriorDistribution).getCopyNumber();
        final int GQ = calculateGenotypeQuality(copyNumberPosteriorDistribution);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(ALL_ALLELES);
        variantContextBuilder.chr(intervalCopyNumberGenotypingData.getContig());
        variantContextBuilder.start(intervalCopyNumberGenotypingData.getStart());
        variantContextBuilder.stop(intervalCopyNumberGenotypingData.getEnd());
        variantContextBuilder.id(String.format(VARIANT_PREFIX + "_%s_%d_%d",
                intervalCopyNumberGenotypingData.getContig(),
                intervalCopyNumberGenotypingData.getStart(),
                intervalCopyNumberGenotypingData.getEnd()));

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final String contig = intervalCopyNumberGenotypingData.getContig();
        final IntegerCopyNumberState refCopyNumber = allosomalContigSet.contains(contig)
                ? intervalCopyNumberGenotypingData.getBaselineIntegerCopyNumberState()
                : refAutosomalCopyNumberState;

        final Allele allele;
        if (copyNumberMAP > refCopyNumber.getCopyNumber()) {
            allele = GATKSVVCFConstants.DUP_ALLELE;
        } else if (copyNumberMAP < refCopyNumber.getCopyNumber()) {
            allele = GATKSVVCFConstants.DEL_ALLELE;
        } else {
            allele = REF_ALLELE;
        }
        genotypeBuilder.alleles(Collections.singletonList(allele));
        genotypeBuilder.attribute(CN, copyNumberMAP);
        genotypeBuilder.attribute(CNLP, copyNumberPLVector);
        genotypeBuilder.attribute(CNQ, GQ);
        final Genotype genotype = genotypeBuilder.make();

        variantContextBuilder.attribute(VCFConstants.END_KEY, intervalCopyNumberGenotypingData.getEnd());
        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }

    /**
     * Calculate copy-number log posterior scores
     *
     * @param copyNumberPosteriorDistribution copy-number log posterior distribution
     * @return a list of copy number log posterior scores
     */
    @VisibleForTesting
    static List<Integer> calculateCopyNumberPLVector(final CopyNumberPosteriorDistribution copyNumberPosteriorDistribution) {
        final Double largestLogProb = copyNumberPosteriorDistribution.getIntegerCopyNumberStateList().stream()
                .mapToDouble(copyNumberPosteriorDistribution::getCopyNumberPosterior)
                .max().orElse(Double.NaN);
        return copyNumberPosteriorDistribution.getIntegerCopyNumberStateList().stream()
                .mapToDouble(copyNumberPosteriorDistribution::getCopyNumberPosterior)
                .map(value -> value - largestLogProb)
                .mapToInt(GermlineCNVIntervalVariantComposer::convertLogProbabilityToPhredScore)
                .boxed()
                .collect(Collectors.toList());
    }

    /**
     * Find MAP copy-number state
     *
     * @param copyNumberPosteriorDistribution copy-number log posterior distribution
     * @return MAP copy number state
     */
    @VisibleForTesting
    static IntegerCopyNumberState calculateMAPCopyNumberState(
            final CopyNumberPosteriorDistribution copyNumberPosteriorDistribution) {
        final Optional<IntegerCopyNumberState> copyNumberStateMAP = copyNumberPosteriorDistribution
                .getIntegerCopyNumberStateList()
                .stream()
                .max(Comparator.comparingDouble(copyNumberPosteriorDistribution::getCopyNumberPosterior));
        return copyNumberStateMAP.get();
    }

    /**
     * Calculate genotype quality score defined as a difference of second smallest and smallest copy number log
     * posterior scores in phred-scale
     *
     * @param copyNumberPosteriorDistribution copy-number log posterior distribution
     * @return genotype quality score
     */
    @VisibleForTesting
    static int calculateGenotypeQuality(final CopyNumberPosteriorDistribution copyNumberPosteriorDistribution) {
        final List<Integer> sortedPosteriors = calculateCopyNumberPLVector(copyNumberPosteriorDistribution).stream()
                .sorted().collect(Collectors.toList());
        /* sanity check */
        Utils.validate(sortedPosteriors.get(0) == 0, "Something went wrong. Smallest copy" +
                " number posterior score must be 0");
        return sortedPosteriors.get(1);
    }

    private static int convertLogProbabilityToPhredScore(final double posteriorProbInLogSpace) {
        return (int) FastMath.floor(-10.0 * posteriorProbInLogSpace * MathUtils.LOG10_E);
    }
}
