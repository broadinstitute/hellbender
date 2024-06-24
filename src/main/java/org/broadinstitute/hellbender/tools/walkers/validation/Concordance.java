package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.lang3.mutable.MutableLong;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Evaluate site-level concordance of an input VCF against a truth VCF.
 *
 * <p>This tool evaluates two variant callsets against each other and produces a six-column summary metrics table. The summary:</p>
 *
 * <ul>
 *     <li>stratifies SNP and INDEL calls,</li>
 *     <li>tallies true-positive, false-positive and false-negative calls,</li>
 *     <li>and calculates sensitivity and precision.</li>
 * </ul>
 *
 * <p>The tool assumes all records in the --truth VCF are passing truth variants. For the -eval VCF, the tool uses only unfiltered passing calls.</p>
 *
 * <p>Optionally, the tool can be set to produce VCFs of the following variant records, annotated with each variant's concordance status:</p>
 * <ul>
 *     <li>True positives and false negatives (i.e. all variants in the truth VCF): useful for calculating sensitivity</li>
 *     <li>True positives and false positives (i.e. all variants in the eval VCF): useful for obtaining a training data
 *     set for machine learning classifiers of artifacts</li>
 * </ul>
 *
 * <p>These output VCFs can be passed to {@link VariantsToTable} to produce a TSV file for statistical analysis in R
 * or Python.</p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk Concordance \
 *   -R reference.fa \
 *   -eval eval.vcf \
 *   --truth truth.vcf \
 *   --summary summary.tsv    
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = Concordance.USAGE_SUMMARY,
        oneLineSummary = Concordance.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class Concordance extends AbstractConcordanceWalker {

    static final String USAGE_ONE_LINE_SUMMARY = "Evaluate concordance of an input VCF against a validated truth VCF";
    static final String USAGE_SUMMARY = "This tool evaluates an input VCF against a VCF that has been validated" +
            " and is considered to represent ground truth.\n" +
            " The summary statistics (# true positives, # false positives, # false negatives, sensitivity, precision) are reported \n" +
            " in a TSV file (--summary). Note that this tool assumes that the truth VCF only contains PASS variants.";

    public static final String SUMMARY_LONG_NAME = "summary";
    public static final String SUMMARY_SHORT_NAME = "S";

    public static final String FILTER_ANALYSIS_LONG_NAME = "filter-analysis";

    public static final String TRUE_POSITIVES_AND_FALSE_NEGATIVES_LONG_NAME = "true-positives-and-false-negatives";
    public static final String TRUE_POSITIVES_AND_FALSE_NEGATIVES_SHORT_NAME = "tpfn";
    public static final String TRUE_POSITIVES_AND_FALSE_POSITIVES_LONG_NAME = "true-positives-and-false-positives";
    public static final String TRUE_POSITIVES_AND_FALSE_POSITIVES_SHORT_NAME = "tpfp";
    public static final String FILTERED_TRUE_NEGATIVES_AND_FALSE_NEGATIVES_LONG_NAME = "filtered-true-negatives-and-false-negatives";
    public static final String FILTERED_TRUE_NEGATIVES_AND_FALSE_NEGATIVES_SHORT_NAME = "ftnfn";
    
    public static final String TRUTH_STATUS_VCF_ATTRIBUTE = "STATUS";
    public static VCFInfoHeaderLine TRUTH_STATUS_HEADER_LINE =
            new VCFInfoHeaderLine(TRUTH_STATUS_VCF_ATTRIBUTE, 1,VCFHeaderLineType.String, "Truth status: TP/FP/FN for true positive/false positive/false negative.");

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)",
            fullName = SUMMARY_LONG_NAME,
            shortName = SUMMARY_SHORT_NAME)
    protected GATKPath summary;

    @Argument(doc = "A table of the contribution of each filter to true and false negatives",
            fullName = FILTER_ANALYSIS_LONG_NAME,
            optional = true)
    protected GATKPath filterAnalysis;

    @Argument(doc = "A vcf to write true positives and false negatives",
            fullName = TRUE_POSITIVES_AND_FALSE_NEGATIVES_LONG_NAME,
            shortName = TRUE_POSITIVES_AND_FALSE_NEGATIVES_SHORT_NAME,
            optional = true)
    protected GATKPath truePositivesAndFalseNegativesVcf = null;

    @Argument(doc = "A vcf to write true positives and false positives",
            fullName = TRUE_POSITIVES_AND_FALSE_POSITIVES_LONG_NAME,
            shortName = TRUE_POSITIVES_AND_FALSE_POSITIVES_SHORT_NAME,
            optional = true)
    protected GATKPath truePositivesAndFalsePositivesVcf = null;

    @Argument(doc = "A vcf to write filtered true negatives and false negatives",
            fullName = FILTERED_TRUE_NEGATIVES_AND_FALSE_NEGATIVES_LONG_NAME,
            shortName = FILTERED_TRUE_NEGATIVES_AND_FALSE_NEGATIVES_SHORT_NAME,
            optional = true)
    protected GATKPath filteredTrueNegativesAndFalseNegativesVcf = null;

    // we count true positives, false positives, false negatives for snps and indels
    private final EnumMap<ConcordanceState, MutableLong> snpCounts = new EnumMap<>(ConcordanceState.class);
    private final EnumMap<ConcordanceState, MutableLong> indelCounts = new EnumMap<>(ConcordanceState.class);
    private VariantContextWriter truePositivesAndFalseNegativesVcfWriter;
    private VariantContextWriter truePositivesAndFalsePositivesVcfWriter;
    private VariantContextWriter filteredTrueNegativesAndFalseNegativesVcfWriter;

    private final Map<String, FilterAnalysisRecord> filterAnalysisRecords = new HashMap<>();

    @Override
    public void onTraversalStart() {
        Set<VCFHeaderLine> defaultToolHeaderLines = getDefaultToolVCFHeaderLines();
        for (final ConcordanceState state : ConcordanceState.values()) {
            snpCounts.put(state, new MutableLong(0));
            indelCounts.put(state, new MutableLong(0));
        }

        final VCFHeader evalHeader = getEvalHeader();

        if (truePositivesAndFalseNegativesVcf != null) {
            truePositivesAndFalseNegativesVcfWriter = createVCFWriter(truePositivesAndFalseNegativesVcf);
            final VCFHeader truthHeader = getTruthHeader();
            truthHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);
            defaultToolHeaderLines.forEach(truthHeader::addMetaDataLine);
            truePositivesAndFalseNegativesVcfWriter.writeHeader(truthHeader);
        }

        if (truePositivesAndFalsePositivesVcf != null) {
            truePositivesAndFalsePositivesVcfWriter = createVCFWriter(truePositivesAndFalsePositivesVcf);

            defaultToolHeaderLines.forEach(evalHeader::addMetaDataLine);
            evalHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);
            truePositivesAndFalsePositivesVcfWriter.writeHeader(evalHeader);
        }

        if (filteredTrueNegativesAndFalseNegativesVcf != null) {
            filteredTrueNegativesAndFalseNegativesVcfWriter = createVCFWriter(filteredTrueNegativesAndFalseNegativesVcf);
            evalHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);
            defaultToolHeaderLines.forEach(evalHeader::addMetaDataLine);
            filteredTrueNegativesAndFalseNegativesVcfWriter.writeHeader(evalHeader);
        }

        final List<String> filtersInVcf = evalHeader.getFilterLines().stream().map(VCFFilterHeaderLine::getID).collect(Collectors.toList());
        filtersInVcf.forEach(filter -> filterAnalysisRecords.put(filter, new FilterAnalysisRecord(filter, 0,0,0,0)));

    }

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        final ConcordanceState concordanceState = truthVersusEval.getConcordance();
        if (truthVersusEval.getTruthIfPresentElseEval().isSNP()) {
            snpCounts.get(concordanceState).increment();
        } else {
            indelCounts.get(concordanceState).increment();
        }

        switch (concordanceState) {
            case TRUE_POSITIVE:
                writeTruePositive(truthVersusEval);
                break;
            case FALSE_POSITIVE:
                writeFalsePositive(truthVersusEval);
                break;
            case FALSE_NEGATIVE:
                writeFalseNegative(truthVersusEval);
                break;
            case FILTERED_TRUE_NEGATIVE:
                writeFilteredTrueNegative(truthVersusEval);
                break;
            case FILTERED_FALSE_NEGATIVE:
                writeFilteredFalseNegative(truthVersusEval);
                break;
            default:
                throw new IllegalStateException("Unexpected ConcordanceState: " + concordanceState.toString());
        }

        if (filterAnalysis != null && concordanceState == ConcordanceState.FILTERED_TRUE_NEGATIVE || concordanceState == ConcordanceState.FILTERED_FALSE_NEGATIVE) {
            final Set<String> filters = truthVersusEval.getEval().getFilters();
            final boolean unique = filters.size() == 1;
            filters.stream().map(filterAnalysisRecords::get).forEach(record -> updateFilterAnalysisRecord(record, concordanceState, unique));
        }
    }

    private void updateFilterAnalysisRecord(final FilterAnalysisRecord record, final ConcordanceState state, final boolean isOnlyFilter) {
        if (state == ConcordanceState.FILTERED_TRUE_NEGATIVE) {
            record.incrementTrueNegative();
            if (isOnlyFilter) {
                record.incrementUniqueTrueNegative();
            }
        } else if (state == ConcordanceState.FILTERED_FALSE_NEGATIVE) {
            record.incrementFalseNegative();
            if (isOnlyFilter) {
                record.incrementUniqueFalseNegative();
            }
        } else {
            throw new IllegalStateException("This method should only be called on a filtered ConcordanceState.");
        }
    }

    private void writeTruePositive(final TruthVersusEval truthVersusEval) {
        final ConcordanceState state = truthVersusEval.getConcordance();
        Utils.validateArg(state == ConcordanceState.TRUE_POSITIVE, "This is not a true positive.");
        tryToWrite(truePositivesAndFalseNegativesVcfWriter, annotateWithConcordanceState(truthVersusEval.getTruth(), state));
        tryToWrite(truePositivesAndFalsePositivesVcfWriter, annotateWithConcordanceState(truthVersusEval.getEval(), state));
    }

    private void writeFalsePositive(final TruthVersusEval truthVersusEval) {
        final ConcordanceState state = truthVersusEval.getConcordance();
        Utils.validateArg(state == ConcordanceState.FALSE_POSITIVE, "This is not a false positive.");
        tryToWrite(truePositivesAndFalsePositivesVcfWriter, annotateWithConcordanceState(truthVersusEval.getEval(), state));
    }

    private void writeFalseNegative(final TruthVersusEval truthVersusEval) {
        final ConcordanceState state = truthVersusEval.getConcordance();
        Utils.validateArg(state == ConcordanceState.FALSE_NEGATIVE, "This is not a false negative.");
        tryToWrite(truePositivesAndFalseNegativesVcfWriter, annotateWithConcordanceState(truthVersusEval.getTruth(), state));
    }

    private void writeFilteredFalseNegative(final TruthVersusEval truthVersusEval) {
        final ConcordanceState state = truthVersusEval.getConcordance();
        Utils.validateArg(state == ConcordanceState.FILTERED_FALSE_NEGATIVE, "This is not a filtered false negative.");
        tryToWrite(truePositivesAndFalseNegativesVcfWriter, annotateWithConcordanceState(truthVersusEval.getTruth(), state));
        tryToWrite(filteredTrueNegativesAndFalseNegativesVcfWriter, annotateWithConcordanceState(truthVersusEval.getEval(), state));
    }

    private void writeFilteredTrueNegative(final TruthVersusEval truthVersusEval) {
        final ConcordanceState state = truthVersusEval.getConcordance();
        Utils.validateArg(state == ConcordanceState.FILTERED_TRUE_NEGATIVE, "This is not a filtered true negative.");
        tryToWrite(filteredTrueNegativesAndFalseNegativesVcfWriter, annotateWithConcordanceState(truthVersusEval.getEval(), state));
    }

    private static void tryToWrite(final VariantContextWriter writer, final VariantContext vc) {
        if (writer != null) {
            writer.add(vc);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try ( ConcordanceSummaryRecord.Writer concordanceSummaryWriter = ConcordanceSummaryRecord.getWriter(summary) ){
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.SNP,
                    snpCounts.get(ConcordanceState.TRUE_POSITIVE).longValue(),
                    snpCounts.get(ConcordanceState.FALSE_POSITIVE).longValue(),
                    snpCounts.get(ConcordanceState.FALSE_NEGATIVE).longValue() + snpCounts.get(ConcordanceState.FILTERED_FALSE_NEGATIVE).longValue()));
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.INDEL,
                    indelCounts.get(ConcordanceState.TRUE_POSITIVE).longValue(),
                    indelCounts.get(ConcordanceState.FALSE_POSITIVE).longValue(),
                    indelCounts.get(ConcordanceState.FALSE_NEGATIVE).longValue() + indelCounts.get(ConcordanceState.FILTERED_FALSE_NEGATIVE).longValue()));
        } catch (IOException e){
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        if (filterAnalysis != null) {
            FilterAnalysisRecord.writeToFile(filterAnalysisRecords.values(), filterAnalysis);
        }

        if (truePositivesAndFalsePositivesVcfWriter != null) {
            truePositivesAndFalsePositivesVcfWriter.close();
        }

        if (truePositivesAndFalseNegativesVcfWriter != null) {
            truePositivesAndFalseNegativesVcfWriter.close();
        }

        if (filteredTrueNegativesAndFalseNegativesVcfWriter != null) {
            filteredTrueNegativesAndFalseNegativesVcfWriter.close();
        }

        return "SUCCESS";
    }

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());

        // We make sure that the truth has at least one alt allele.
        // If it does, we pick the first for comparison:
        final boolean containsAltAllele =
                (truth.getAlternateAlleles().size() == eval.getAlternateAlleles().size()) &&
                ((truth.getAlternateAlleles().size() > 0) &&
                        eval.getAlternateAlleles().contains(truth.getAlternateAllele(0)));

        return sameRefAllele && containsAltAllele;
    }

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> !vc.isFiltered() && ! vc.isSymbolicOrSV();
    }

    private VariantContext annotateWithConcordanceState(final VariantContext vc, final ConcordanceState state) {
        return new VariantContextBuilder(vc).attribute(TRUTH_STATUS_VCF_ATTRIBUTE, state.getAbbreviation()).make();
    }
}
