package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import picard.sam.util.Pair;

import java.io.OutputStreamWriter;

/**
 * Evaluate GVCF reference block concordance of an input GVCF against a truth GVCF.
 *
 * <p>This tool evaluates reference blocks of two GVCF files against each other and produces three histograms:</p>
 *
 * <ul>
 *     <li>Truth block histogram: Indicates the number of occurrences of reference blocks with a given confidence score and length in the truth GVCF</li>
 *     <li>Eval block histogram: Indicates the number of occurrences of reference blocks with a given confidence score and length in the eval GVCF</li>
 *     <li>Confidence concordance histogram: Reflects the confidence scores of bases in reference blocks in the truth and eval VCF, respectively. An entry of 10 at bin "80,90" means that there are 10 bases which simultaneously have a reference confidence of 80 in the truth GVCF and a reference confidence of 90 in the eval GVCF.</li>
 * </ul>
 *
 * <p>This tool only considers bases in reference blocks and, in contrast to the {@link Concordance} tool, regardless of passing or failing filters.</p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk ReferenceBlockConcordance \
 *   -R reference.fa \
 *   -eval eval.vcf \
 *   --truth truth.vcf \
 *   --truth-block-histogram truth_block_histogram.tsv \
 *   --eval-block-histogram eval_block_histogram.tsv \
 *   --confidence-concordance-histogram confidence_concordance_histogram.tsv
 * </pre>
 */

@CommandLineProgramProperties(
        summary = ReferenceBlockConcordance.USAGE_SUMMARY,
        oneLineSummary = ReferenceBlockConcordance.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
public class ReferenceBlockConcordance extends AbstractConcordanceWalker {
    public static final String TRUTH_BLOCK_HISTOGRAM_LONG_NAME = "truth-block-histogram";
    public static final String TRUTH_BLOCK_HISTOGRAM_SHORT_NAME = "tbh";
    public static final String EVAL_BLOCK_HISTOGRAM_LONG_NAME = "eval-block-histogram";
    public static final String EVAL_BLOCK_HISTOGRAM_SHORT_NAME = "ebh";
    public static final String CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME = "confidence-concordance-histogram";
    public static final String CONFIDENCE_CONCORDANCE_HISTOGRAM_SHORT_NAME = "cch";

    static final String USAGE_ONE_LINE_SUMMARY = "Evaluate GVCF reference block concordance of an input GVCF against a truth GVCF";
    static final String USAGE_SUMMARY = "This tool evaluates two GVCF files against each other and produces three histograms:\n" +
            "Truth block histogram: Indicates the number of occurrence of reference blocks with a given confidence score and length in the truth GVCF\n" +
            "Eval block histogram: Indicates the number of occurrence of reference blocks with a given confidence score and length in the eval GVCF\n" +
            "Confidence concordance histogram: Reflects the confidence scores of bases in reference blocks in the truth and eval VCF, respectively. An entry of 10 at bin \"80,90\" means that 10 bases in the truth GVCF have a confidence score of 80 while those same bases have a score of 90 in the eval GVCF.\n" +
            "In contrast to the Concordance tool, this tool considers all variants, regardless of passing or failing filters.";

    @Argument(doc = "A histogram of block lengths and their associated confidence scores for the truth sample",
            fullName = TRUTH_BLOCK_HISTOGRAM_LONG_NAME,
            shortName = TRUTH_BLOCK_HISTOGRAM_SHORT_NAME)
    protected GATKPath truthBlockHistogramFile;
    @Argument(doc = "A histogram of block lengths and their associated confidence scores for the eval sample",
            fullName = EVAL_BLOCK_HISTOGRAM_LONG_NAME,
            shortName = EVAL_BLOCK_HISTOGRAM_SHORT_NAME)
    protected GATKPath evalBlockHistogramFile;
    @Argument(doc = "Reflects the confidence scores of bases in reference blocks in the truth and eval VCF, respectively. An entry of 10 at bin \"80,90\" means that 10 bases in the truth GVCF have a confidence score of 80 while those same bases have a score of 90 in the eval GVCF.",
            fullName = CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME,
            shortName = CONFIDENCE_CONCORDANCE_HISTOGRAM_SHORT_NAME)
    protected GATKPath confidenceConcordanceHistogramFile;

    // TODO this should be a Histogram<Pair<Integer, Integer>>, however, the MetricsFile class cannot read
    // arbitrary types, therefore, it must be converted to a String, which may be slower
    private final Histogram<String> truthBlockHistogram = new Histogram<>();
    private final Histogram<String> evalBlockHistogram = new Histogram<>();
    private final Histogram<String> confidenceConcordanceHistogram = new Histogram<>();

    private VariantContext currentTruthVariantContext = null;
    private VariantContext currentEvalVariantContext = null;

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        // Only allow HomRef sites
        return this::isHomRef;
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() {
        // Only allow HomRef sites
        return this::isHomRef;
    }

    private boolean isHomRef(VariantContext variantContext) {
        return variantContext.getGenotypes().get(0).isHomRef();
    }

    private Pair<Integer, Integer> extractLengthAndGQ(VariantContext variant) {
        if (variant.getGenotypes().size() != 1) {
            throw new IllegalStateException(String.format("A multisample GVCF file was provided, however, only single sample GVCFs are currently supported. This occurred when reading \"%s\".", variant.toStringDecodeGenotypes()));
        }
        final Genotype genotype = variant.getGenotype(0);
        return new Pair<>(variant.getLengthOnReference(), genotype.getGQ());
    }

    @Override
    protected void apply(TruthVersusEval truthVersusEval, ReadsContext readsContext, ReferenceContext refContext) {
        // Truth
        if (truthVersusEval.hasTruth()) {
            currentTruthVariantContext = truthVersusEval.getTruth();

            truthBlockHistogram.increment(extractLengthAndGQ(truthVersusEval.getTruth()).toString());
        }

        // Eval
        if (truthVersusEval.hasEval()) {
            currentEvalVariantContext = truthVersusEval.getEval();

            evalBlockHistogram.increment(extractLengthAndGQ(truthVersusEval.getEval()).toString());
        }

        if (currentTruthVariantContext != null && !currentTruthVariantContext.overlaps(truthVersusEval)) {
            currentTruthVariantContext = null;
        }
        if (currentEvalVariantContext != null && !currentEvalVariantContext.overlaps(truthVersusEval)) {
            currentEvalVariantContext = null;
        }

        // Evaluate only when currently seeing two NON_REF blocks
        if (currentTruthVariantContext != null && currentEvalVariantContext != null) {
            final SimpleInterval truthInterval = new SimpleInterval(currentTruthVariantContext);
            final SimpleInterval evalInterval = new SimpleInterval(currentEvalVariantContext);
            if (truthInterval.overlaps(evalInterval)) {
                confidenceConcordanceHistogram.increment(new Pair<>(currentTruthVariantContext.getGenotype(0).getGQ(), currentEvalVariantContext.getGenotype(0).getGQ()).toString(), truthInterval.intersect(evalInterval).getLengthOnReference());
            }
        }
    }

    @Override
    protected boolean shouldVariantsBeMatched(VariantContext truth, VariantContext eval) {
        return true;
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();

        // Truth histogram
        final MetricsFile<?, String> truthBlockMetricsFile = getMetricsFile();
        truthBlockMetricsFile.addHistogram(truthBlockHistogram);
        truthBlockMetricsFile.write(new OutputStreamWriter(truthBlockHistogramFile.getOutputStream()));

        // Eval histogram
        final MetricsFile<?, String> evalBlockMetricsFile = getMetricsFile();
        evalBlockMetricsFile.addHistogram(evalBlockHistogram);
        evalBlockMetricsFile.write(new OutputStreamWriter(evalBlockHistogramFile.getOutputStream()));

        // Confidence concordance
        final MetricsFile<?, String> confidenceConcordanceMetricsFile = getMetricsFile();
        confidenceConcordanceMetricsFile.addHistogram(confidenceConcordanceHistogram);
        confidenceConcordanceMetricsFile.write(new OutputStreamWriter(confidenceConcordanceHistogramFile.getOutputStream()));

        return "SUCCESS";
    }
}
