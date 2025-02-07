package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
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
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.ClusteringParameters;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.tools.sv.cluster.StratifiedClusteringTableParser;
import org.broadinstitute.hellbender.tools.sv.concordance.ClosestSVFinder;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceAnnotator;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceLinkage;
import org.broadinstitute.hellbender.tools.sv.concordance.StratifiedConcordanceEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.OptionalSVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import picard.vcf.GenotypeConcordance;

import java.io.IOException;
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
 * This tool also allows supports stratification of SVs into groups with specified matching criteria including SV type,
 * size range, and interval overlap. Please see the {@link GroupedSVCluster} tool documentation for further details
 * on how to specify stratification groups.
 *
 * Note that unlike {@link GroupedSVCluster}, this tool allows any variant to
 * match more than one stratification group. If this occurs, groups will be prioritized by their ordering in the input
 * stratification table, with groups appearing first receiving higher priority. While all matching groups will
 * be listed in the STRAT INFO field, the variant ID pertaining to the highest-priority group will be populated in
 * the TRUTH_VID field (groups with no matching variant are ignored). It is therefore recommended that the groups with
 * the most specific clustering criteria be listed as higher priority.
 *
 * The "default" stratification group, with clustering parameters specified directly through the clustering program
 * arguments (e.g. --depth-breakend-window, --pesr-interval-overlap, etc.), is always present and given lowest priority.
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

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected GATKPath outputFile;

    /**
     * Expected format is tab-delimited and contains columns NAME, RECIPROCAL_OVERLAP, SIZE_SIMILARITY, BREAKEND_WINDOW,
     * SAMPLE_OVERLAP. First line must be a header with column names. Comment lines starting with
     * {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    @Argument(
            doc = "Configuration file (.tsv) containing the clustering parameters for each group",
            fullName = GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME,
            optional = true
    )
    public GATKPath strataClusteringConfigFile;

    @Argument(fullName = SVClusterWalker.MAX_RECORDS_IN_RAM_LONG_NAME,
            doc = "When writing VCF files that need to be sorted, this will specify the number of records stored in " +
                    "RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a " +
                    "VCF file, and increases the amount of RAM needed.",
            optional=true)
    public int maxRecordsInRam = 10000;

    @ArgumentCollection
    protected final SVClusterEngineArgumentsCollection defaultClusteringArgs = new SVClusterEngineArgumentsCollection();
    @ArgumentCollection
    private final OptionalSVStratificationEngineArgumentsCollection stratArgs = new OptionalSVStratificationEngineArgumentsCollection();

    protected StratifiedConcordanceEngine engine;
    protected SAMSequenceDictionary dictionary;
    protected SortingCollection<VariantContext> sortingBuffer;
    protected VariantContextWriter writer;


    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        // Use master sequence dictionary i.e. hg38 .dict file since the "best" dictionary is grabbed
        // from the VCF, which is sometimes out of order
        dictionary = getMasterSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        // Check that vcfs are sorted the same
        SequenceDictionaryUtils.validateDictionaries("eval", getEvalHeader().getSequenceDictionary(),
                "truth", getTruthHeader().getSequenceDictionary(), false, true);
        writer = createVCFWriter(outputFile);
        final VCFHeader header = createHeader(getEvalHeader());
        writer.writeHeader(header);
        sortingBuffer = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(header, true),
                header.getVCFRecordComparator(),
                maxRecordsInRam,
                tmpDir.toPath());

        // Concordance computations should be done on common samples only
        final Set<String> commonSamples = Sets.intersection(
                new HashSet<>(getEvalHeader().getGenotypeSamples()),
                new HashSet<>(getTruthHeader().getGenotypeSamples()));
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(commonSamples);

        final Map<String, ClosestSVFinder> clusterEngineMap = new HashMap<>();
        if (strataClusteringConfigFile != null) {
            try (final TableReader<StratifiedClusteringTableParser.StratumParameters> tableReader = TableUtils.reader(strataClusteringConfigFile.toPath(), StratifiedClusteringTableParser::tableParser)) {
                for (final StratifiedClusteringTableParser.StratumParameters parameters : tableReader) {
                    // Identical parameters for each linkage type
                    final ClusteringParameters pesrParams = ClusteringParameters.createPesrParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                    final ClusteringParameters mixedParams = ClusteringParameters.createMixedParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                    final ClusteringParameters depthParams = ClusteringParameters.createDepthParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                    final SVConcordanceLinkage linkage = new SVConcordanceLinkage(dictionary);
                    linkage.setDepthOnlyParams(depthParams);
                    linkage.setMixedParams(mixedParams);
                    linkage.setEvidenceParams(pesrParams);
                    final ClosestSVFinder engine = new ClosestSVFinder(linkage, collapser::annotate, false, dictionary);
                    clusterEngineMap.put(parameters.name(), engine);
                }
            } catch (final IOException e) {
                throw new GATKException("IO error while reading config table", e);
            }
        }

        if ((stratArgs.configFile == null) ^ (strataClusteringConfigFile == null)) {
            throw new UserException.BadInput("Both --" + OptionalSVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME
                    + " and --" + GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME + " must be used together, but only one was specified.");
        }
        final SVStratificationEngine stratEngine = SVStratify.loadStratificationConfig(stratArgs.configFile, stratArgs, dictionary);
        engine = new StratifiedConcordanceEngine(clusterEngineMap, stratEngine, stratArgs, defaultClusteringArgs, collapser, dictionary);
    }


    @Override
    public Object onTraversalSuccess() {
        for (final VariantContext variant : engine.flush(true)) {
            sortingBuffer.add(variant);
        }
        if (!engine.isEmpty()) {
            throw new GATKException("Concordance engine is not empty, but it should be");
        }
        for (final VariantContext variant : sortingBuffer) {
            writer.add(variant);
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        if (sortingBuffer != null) {
            sortingBuffer.cleanup();
        }
        if (writer != null) {
            writer.close();
        }
        super.closeTool();
    }

    @Override
    public void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        if (truthVersusEval.hasTruth()) {
            final SVCallRecord record = minimizeTruthFootprint(SVCallRecordUtils.create(truthVersusEval.getTruth(), dictionary));
            engine.addTruthVariant(record);
        }
        if (truthVersusEval.hasEval()) {
            final SVCallRecord record = SVCallRecordUtils.create(truthVersusEval.getEval(), dictionary);
            engine.addEvalVariant(record);
        }
        for (final VariantContext variant : engine.flush(false)) {
            sortingBuffer.add(variant);
        }
    }

    /**
     * Strips unneeded attributes from a truth variant to save memory.
     */
    private SVCallRecord minimizeTruthFootprint(final SVCallRecord item) {
        final List<Genotype> genotypes = item.getGenotypes().stream().map(SVConcordance::stripTruthGenotype).collect(Collectors.toList());
        return new SVCallRecord(item.getId(), item.getContigA(), item.getPositionA(),
                item.getStrandA(), item.getContigB(), item.getPositionB(), item.getStrandB(), item.getType(),
                item.getComplexSubtype(), item.getComplexEventIntervals(), item.getLength(), item.getEvidence(), item.getAlgorithms(),
                item.getAlleles(), genotypes, item.getAttributes(), item.getFilters(), item.getLog10PError(), dictionary);
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

    protected VCFHeader createHeader(final VCFHeader header) {
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
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_RECIPROCAL_OVERLAP_INFO, 1, VCFHeaderLineType.Float, "Reciprocal overlap with truth variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_SIZE_SIMILARITY_INFO, 1, VCFHeaderLineType.Float, "Size similarity with truth variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_DISTANCE_START_INFO, 1, VCFHeaderLineType.Integer, "Start coordinate distance in bp to truth variant's start"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_DISTANCE_END_INFO, 1, VCFHeaderLineType.Integer, "End coordinate distance in bp to truth variant's end"));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        SVStratify.addStratifyMetadata(header);
        return header;
    }
}
