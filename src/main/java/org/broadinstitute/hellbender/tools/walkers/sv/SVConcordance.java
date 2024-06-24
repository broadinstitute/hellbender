package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.concordance.ClosestSVFinder;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceAnnotator;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceLinkage;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import picard.vcf.GenotypeConcordance;

import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>This tool calculates SV genotype concordance between an "evaluation" VCF and a "truth" VCF. For each evaluation
 * variant, a single truth variant is matched based on the following order of criteria:</p>
 *
 * <ol>
 *     <li>Total breakend distance</li>
 *     <li>Min breakend distance (among the two sides)</li>
 *     <li>Genotype concordance</li>
 * </ol>
 *
 * after meeting minimum overlap criteria. Evaluation VCF variants that are sucessfully matched are annotated with
 * genotype concordance metrics, including allele frequency of the truth variant. Concordance metrics are computed
 * on the intersection of sample sets of the two VCFs, but all other annotations including variant truth status
 * and allele frequency use all records and samples available. See output header for descriptions
 * of the specific fields. For multi-allelic CNVs, only a copy state concordance metric is
 * annotated. Allele frequencies will be recalculated automatically if unavailable in the provided VCFs.
 *
 * For large inputs, users may enable the --do-not-sort flag to reduce memory usage.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Evaluation VCF
 *     </li>
 *     <li>
 *         Truth VCF
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         The evaluation VCF annotated with genotype concordance metrics
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVConcordance \
 *       --sequence-dictionary ref.dict \
 *       --eval evaluation.vcf.gz \
 *       --truth truth.vcf.gz \
 *       -O output.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Annotates structural variant genotype concordance",
        oneLineSummary = "Annotates structural variant genotype concordance",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVConcordance extends AbstractConcordanceWalker {

    public static final String UNSORTED_OUTPUT_LONG_NAME = "do-not-sort";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "Skip output sorting. This can substantially reduce memory usage.",
            fullName = UNSORTED_OUTPUT_LONG_NAME
    )
    private boolean doNotSort;

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private SVConcordanceLinkage linkage;
    private ClosestSVFinder engine;
    private String currentContig = null;

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> true;
    }

    @Override
    public void onTraversalStart() {
        // Use master sequence dictionary i.e. hg38 .dict file since the "best" dictionary is grabbed
        // from the VCF, which is sometimes out of order
        dictionary = getMasterSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        // Check that vcfs are sorted the same
        SequenceDictionaryUtils.validateDictionaries("eval", getEvalHeader().getSequenceDictionary(),
                "truth", getTruthHeader().getSequenceDictionary(), false, true);

        linkage = new SVConcordanceLinkage(dictionary);
        linkage.setDepthOnlyParams(clusterParameterArgs.getDepthParameters());
        linkage.setMixedParams(clusterParameterArgs.getMixedParameters());
        linkage.setEvidenceParams(clusterParameterArgs.getPESRParameters());

        // Concordance computations should be done on common samples only
        final Set<String> commonSamples = Sets.intersection(
                new HashSet<>(getEvalHeader().getGenotypeSamples()),
                new HashSet<>(getTruthHeader().getGenotypeSamples()));
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(commonSamples);
        engine = new ClosestSVFinder(linkage, collapser::annotate, !doNotSort, dictionary);

        if (doNotSort) {
            createOutputVariantIndex = false;
        }
        writer = createVCFWriter(outputFile);
        writer.writeHeader(createHeader(getEvalHeader()));
    }

    @Override
    public Object onTraversalSuccess() {
        flushClusters(true);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    @Override
    public void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        if (truthVersusEval.hasTruth()) {
            add(truthVersusEval.getTruth(), true);
        }
        if (truthVersusEval.hasEval()) {
            add(truthVersusEval.getEval(), false);
        }
    }

    private void add(final VariantContext variant, final boolean isTruth) {
        SVCallRecord record = SVCallRecordUtils.create(variant);
        if (!record.getContigA().equals(currentContig)) {
            flushClusters(true);
            currentContig = record.getContigA();
        }
        if (isTruth) {
            record = minimizeTruthFootprint(record);
        }
        engine.add(record, isTruth);
        flushClusters(false);
    }

    /**
     * Strips unneeded attributes from a truth variant to save memory.
     */
    protected SVCallRecord minimizeTruthFootprint(final SVCallRecord item) {
        final List<Genotype> genotypes = item.getGenotypes().stream().map(SVConcordance::stripTruthGenotype).collect(Collectors.toList());
        return new SVCallRecord(item.getId(), item.getContigA(), item.getPositionA(),
                item.getStrandA(), item.getContigB(), item.getPositionB(), item.getStrandB(), item.getType(),
                item.getCpxSubtype(), item.getLength(), item.getAlgorithms(), item.getAlleles(), genotypes,
                item.getAttributes(), item.getFilters(), item.getLog10PError(), dictionary);
    }

    private static Genotype stripTruthGenotype(final Genotype genotype) {
        final GenotypeBuilder builder = new GenotypeBuilder(genotype.getSampleName()).alleles(genotype.getAlleles());
        if (genotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
            builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT));
        }
        return builder.make();
    }

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        return true;
    }

    private void flushClusters(final boolean force) {
        engine.flush(force).stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .forEach(writer::add);
    }

    private VCFHeader createHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFFormatHeaderLine(GenotypeConcordance.CONTINGENCY_STATE_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "The genotype concordance contingency state"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.TRUTH_CN_EQUAL_FORMAT, 1, VCFHeaderLineType.Integer, "Truth CNV copy state is equal (1=True, 0=False)"));
        header.addMetaDataLine(Concordance.TRUTH_STATUS_HEADER_LINE);
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "Genotype concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "CNV copy number concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NON_REF_GENOTYPE_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_PPV_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_PPV_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_PPV_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_SPECIFICITY_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype specificity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, 1, VCFHeaderLineType.String, "Matching truth set variant id"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Truth set allele count"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Truth set allele number"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Truth set allele frequency"));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        return header;
    }


}
