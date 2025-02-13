package org.broadinstitute.hellbender.tools.walkers.validation;

import java.nio.file.Paths;
import java.io.IOException;

import org.apache.commons.collections4.Predicate;

import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.AlleleConcordanceWalker;
import org.broadinstitute.hellbender.exceptions.UserException;

import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

/**
 * Compare INFO field values between two VCFs or compare two different INFO fields from one VCF.
 * We only evaluate sites that are in both VCFs.
 * Although we use the arguments eval and truth, we only compare the scores, we do not determine the correct score.
 * Either VCF can be used as eval or truth, or the same VCF can be used for both.
 * Differences greater than the epsilon argument will trigger a warning.
 *
 * <h3>Compare the CNN_2D info fields for the same sites from two different VCFs:</h3>
 *
 * <pre>
 * gatk EvaluateInfoFieldConcordance \
 *  -eval a.vcf \
 *  -truth another.vcf \
 *  -S summary.txt \
 *  -eval-info-key CNN_2D \
 *  -truth-info-key CNN_2D \
 *  -epsilon 0.01
 * </pre>
 *
 * <h3>Compare the CNN_2D info field with the CNN_1D field from the same sites in one VCF:</h3>
 *
 * <pre>
 * gatk EvaluateInfoFieldConcordance \
 *  -eval my.vcf \
 *  -truth my.vcf \
 *  -S summary.txt \
 *  -eval-info-key CNN_2D \
 *  -truth-info-key CNN_1D \
 *  -epsilon 0.01
 * </pre>
 */
@CommandLineProgramProperties(
        summary=EvaluateInfoFieldConcordance.USAGE_SUMMARY,
        oneLineSummary=EvaluateInfoFieldConcordance.USAGE_ONE_LINE_SUMMARY,
        programGroup=VariantEvaluationProgramGroup.class)
@DocumentedFeature
@BetaFeature
public class EvaluateInfoFieldConcordance extends AlleleConcordanceWalker {
    static final String USAGE_ONE_LINE_SUMMARY = "Evaluate concordance of info fields in an input VCF against a validated truth VCF";
    static final String USAGE_SUMMARY = "This tool evaluates info fields from an input VCF against a VCF that has been validated and is considered to represent ground truth.\n";
    public static final String SUMMARY_LONG_NAME = "summary";
    public static final String SUMMARY_SHORT_NAME = "S";

    @Argument(doc="A table of summary statistics (true positives, sensitivity, etc.)", fullName=SUMMARY_LONG_NAME, shortName=SUMMARY_SHORT_NAME)
    protected String summary;

    @Argument(fullName="eval-info-key", shortName="eval-info-key", doc="Info key from eval vcf")
    protected String evalInfoKey;

    @Argument(fullName="truth-info-key", shortName="truth-info-key", doc="Info key from truth vcf")
    protected String truthInfoKey;

    @Advanced
    @Argument(fullName="warn-big-differences",
            shortName="warn-big-differences",
            doc="If set differences in the info key values greater than epsilon will trigger warnings.",
            optional=true)
    protected boolean warnBigDifferences = false;

    @Advanced
    @Argument(fullName="epsilon", shortName="epsilon", doc="Difference tolerance", optional=true)
    protected double epsilon = 0.1;

    private int snpCount = 0;
    private int indelCount = 0;

    private double snpSumDelta = 0.0;
    private double snpSumDeltaSquared = 0.0;
    private double indelSumDelta = 0.0;
    private double indelSumDeltaSquared = 0.0;

    @Override
    public void onTraversalStart() {
        if(getEvalHeader().getInfoHeaderLine(evalInfoKey) == null){
            throw new UserException("Missing key:"+evalInfoKey+" in Eval VCF:"+evalVariantsFile);
        }

        if(getTruthHeader().getInfoHeaderLine(truthInfoKey) == null){
            throw new UserException("Missing key:"+truthInfoKey+" in Truth VCF:"+truthVariantsFile);
        }
    }

    @Override
    protected void apply(AlleleTruthVersusEval truthVersusEval, ReadsContext readsContext, ReferenceContext refContext) {
        ConcordanceState concordanceState = truthVersusEval.getConcordance();
        switch (concordanceState) {
            case TRUE_POSITIVE: {
                if(truthVersusEval.getEval().isSNP()){
                    snpCount++;
                } else if (truthVersusEval.getEval().isIndel()) {
                    indelCount++;
                }
                this.infoDifference(truthVersusEval.getEval(), truthVersusEval.getTruth());
                break;
            }
            case FALSE_POSITIVE:
            case FALSE_NEGATIVE:
            case FILTERED_TRUE_NEGATIVE:
            case FILTERED_FALSE_NEGATIVE: {
                break;
            }
            default: {
                throw new IllegalStateException("Unexpected ConcordanceState: " + concordanceState.toString());
            }
        }
    }

    private void infoDifference(final VariantContext eval, final VariantContext truth) {
        if(eval.hasAttribute(this.evalInfoKey) && truth.hasAttribute(truthInfoKey)) {
            final double evalVal = Double.valueOf((String) eval.getAttribute(this.evalInfoKey));
            final double truthVal = Double.valueOf((String) truth.getAttribute(this.truthInfoKey));
            final double delta = evalVal - truthVal;
            final double deltaSquared = delta * delta;
            if (eval.isSNP()) {
                this.snpSumDelta += Math.sqrt(deltaSquared);
                this.snpSumDeltaSquared += deltaSquared;
            } else if (eval.isIndel()) {
                this.indelSumDelta += Math.sqrt(deltaSquared);
                this.indelSumDeltaSquared += deltaSquared;
            }
            if (warnBigDifferences && Math.abs(delta) > this.epsilon) {
                this.logger.warn(String.format("Difference (%f) greater than epsilon (%f) at %s:%d %s:", delta, this.epsilon, eval.getContig(), eval.getStart(), eval.getAlleles().toString()));
                this.logger.warn(String.format("\t\tTruth info: " + truth.getAttributes().toString()));
                this.logger.warn(String.format("\t\tEval info: " + eval.getAttributes().toString()));
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final double snpMean = this.snpSumDelta / snpCount;
        final double snpVariance = (this.snpSumDeltaSquared - this.snpSumDelta * this.snpSumDelta / snpCount) / snpCount;
        final double snpStd = Math.sqrt(snpVariance);
        final double indelMean = this.indelSumDelta / indelCount;
        final double indelVariance = (this.indelSumDeltaSquared - this.indelSumDelta * this.indelSumDelta / indelCount) / indelCount;
        final double indelStd = Math.sqrt(indelVariance);

        this.logger.info(String.format("SNP average delta %f and standard deviation: %f", snpMean, snpStd));
        this.logger.info(String.format("INDEL average delta %f and standard deviation: %f", indelMean, indelStd));

        try (final InfoConcordanceRecord.InfoConcordanceWriter
                     concordanceWriter = InfoConcordanceRecord.getWriter(Paths.get(this.summary))){
            concordanceWriter.writeRecord(new InfoConcordanceRecord(VariantContext.Type.SNP, this.evalInfoKey, this.truthInfoKey, snpMean, snpStd));
            concordanceWriter.writeRecord(new InfoConcordanceRecord(VariantContext.Type.INDEL, this.evalInfoKey, this.truthInfoKey, indelMean, indelStd));
        } catch (IOException e) {
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        return "SUCCESS";
    }

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> !vc.isFiltered() && !vc.isSymbolicOrSV();
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() {
        return vc -> !vc.isFiltered() && !vc.isSymbolicOrSV();
    }

}
