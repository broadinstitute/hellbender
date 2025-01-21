package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.tools.sv.stratify.RequiredSVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * <p>Clusters structural variants using the same base algorithms as {@link SVCluster}. In addition, variants are
 * grouped according to customizable stratification criteria including:
 * <ul>
 *     <li>SV type</li>
 *     <li>Size range</li>
 *     <li>Reference track overlap</li>
 * </ul>
 * The first step is to define these groups in a stratification configuration TSV file. Please see the
 * {@link SVStratify} tool for a description of the stratification method and expected table format.
 *
 * <p>Each SV is only clustered with other SVs in its own group. Each group must be mutually exclusive, meaning that
 * any given SV should only belong to one group. Furthermore, SVs that do not fall into any of the groups will not be
 * clustered.</p>
 *
 * <p>The second step is to define the clustering configuration for each group. This is again done by creating a TSV
 * file with the following columns defined on the first line:
 * <ol>
 *     <li>NAME</li>
 *     <li>RECIPROCAL_OVERLAP</li>
 *     <li>SIZE_SIMILARITY</li>
 *     <li>BREAKEND_WINDOW</li>
 *     <li>SAMPLE_OVERLAP</li>
 * </ol>
 * where NAME corresponds to the same name given in the stratification configuration. Every group needs to be given
 * a configuration here. That is, there should be a 1:1 correspondence of the rows in the two configuration files
 * (order does not matter).
 * </p>
 *
 * <p>The remaining columns define the clustering parameters for the group. See {@link SVCluster} for more information
 * on the different parameters. Note that, unlike {@link SVCluster}, distinct parameter sets for depth-only,
 * PESR, and "mixed" clustering cannot be defined for this tool. Instead, the same parameters are applied to
 * all three cases.</p>
 *
 * <p>For example,</p>
 * <table border="1">
 *   <tr>
 *     <td>NAME</td><td>RECIPROCAL_OVERLAP</td><td>SIZE_SIMILARITY</td><td>BREAKEND_WINDOW</td><td>SAMPLE_OVERLAP</td>
 *   </tr>
 *   <tr>
 *     <td>DEL_large_SD</td><td>0.3</td><td>0.5</td><td>1000000</td><td>0.1</td>
 *   </tr>
 *   <tr>
 *     <td>DUP_large_SD</td><td>0.3</td><td>0.5</td><td>1000000</td><td>0.1</td>
 *   </tr>
 *   <tr>
 *     <td>DEL_small_SR_RM</td><td>0.5</td><td>0.5</td><td>100</td><td>0.1</td>
 *   </tr>
 *   <tr>
 *     <td>DUP_small_SR_RM</td><td>0.5</td><td>0.5</td><td>100</td><td>0.1</td>
 *   </tr>
 *   <tr>
 *     <td>INS_SR</td><td>0.5</td><td>0.5</td><td>100</td><td>0</td>
 *   </tr>
 * </table>
 *
 * <p>This tool accepts multiple VCF inputs with no restrictions on site or sample overlap.</p>
 *
 * <p>This tool does not support CNV defragmentation via the {@link #algorithm} parameter.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SV VCFs
 *     </li>
 *     <li>
 *         Stratification configuration TSV file
 *     </li>
 *     <li>
 *         Clustering configuration TSV file
 *     </li>
 *     <li>
 *         Reference fasta
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Clustered VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk GroupedSVCluster \
 *       -R reference.fasta \
 *       -V variants.vcf.gz \
 *       -O clustered.vcf.gz \
 *       --track-name repeatmasker \
 *       --track-intervals repeatmasker.bed \
 *       --stratify-config strata.tsv \
 *       --clustering-config cluster.tsv
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Clusters structural variants within independent stratification groups",
        oneLineSummary = "Clusters structural variants grouping by type, size, and track overlap",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class GroupedSVCluster extends SVClusterWalker {
    public static final String CLUSTERING_CONFIG_FILE_LONG_NAME = "clustering-config";

    @ArgumentCollection
    private final RequiredSVStratificationEngineArgumentsCollection stratArgs = new RequiredSVStratificationEngineArgumentsCollection();

    /**
     * Expected format is tab-delimited and contains columns NAME, RECIPROCAL_OVERLAP, SIZE_SIMILARITY, BREAKEND_WINDOW,
     * SAMPLE_OVERLAP. First line must be a header with column names. Comment lines starting with
     * {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    @Argument(
            doc = "Configuration file (.tsv) containing the clustering parameters for each group",
            fullName = CLUSTERING_CONFIG_FILE_LONG_NAME
    )
    public GATKPath strataClusteringConfigFile;

    private SVStratificationEngine stratEngine;
    private final Map<String, SVClusterEngine> clusterEngineMap = new HashMap<>();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        // sorting not guaranteed
        createOutputVariantIndex = false;
        stratEngine = SVStratify.loadStratificationConfig(stratArgs.configFile, stratArgs, dictionary);
        Utils.validate(!stratEngine.getStrata().isEmpty(),
                "No strata defined with --" + SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME);
        readStrataClusteringConfig();
        Utils.validate(stratEngine.getStrata().size() == clusterEngineMap.size(),
                "Stratification and clustering configurations have a different number of groups.");
        for (final SVStratificationEngine.Stratum stratum : stratEngine.getStrata()) {
            Utils.validate(clusterEngineMap.containsKey(stratum.getName()),
                    "Could not find group " + stratum.getName() + " in clustering configuration.");
        }
    }

    @Override
    protected VCFHeader createHeader() {
        final VCFHeader header = super.createHeader();
        SVStratify.addStratifyMetadata(header);
        return header;
    }

    private void readStrataClusteringConfig() {
        try (final TableReader<StratifiedClusteringTableParser.StratumParameters> tableReader = TableUtils.reader(strataClusteringConfigFile.toPath(), StratifiedClusteringTableParser::tableParser)) {
            for (final StratifiedClusteringTableParser.StratumParameters parameters : tableReader) {
                // Identical parameters for each linkage type
                final ClusteringParameters pesrParams = ClusteringParameters.createPesrParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                final ClusteringParameters mixedParams = ClusteringParameters.createMixedParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                final ClusteringParameters depthParams = ClusteringParameters.createDepthParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                final SVClusterEngine clusterEngine = createClusteringEngine(pesrParams, mixedParams, depthParams);
                clusterEngineMap.put(parameters.name(), clusterEngine);
            }
        } catch (final IOException e) {
            throw new GATKException("IO error while reading config table", e);
        }
    }

    private SVClusterEngine createClusteringEngine(final ClusteringParameters pesrParams, final ClusteringParameters mixedParams, final ClusteringParameters depthParams) {
        if (algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE || algorithm == CLUSTER_ALGORITHM.MAX_CLIQUE) {
            final SVClusterEngine.CLUSTERING_TYPE type = algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE ?
                    SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE : SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
            return SVClusterEngineFactory.createCanonical(type, breakpointSummaryStrategy,
                    altAlleleSummaryStrategy, dictionary, reference, enableCnv,
                    depthParams, mixedParams, pesrParams);
        } else {
            throw new IllegalArgumentException("Unsupported algorithm: " + algorithm.name());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        for (final SVClusterEngine engine : clusterEngineMap.values()) {
            engine.flush().stream().forEach(this::write);
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
    }

    @Override
    public void applyRecord(final SVCallRecord record) {
        final Collection<SVStratificationEngine.Stratum> stratifications = stratEngine.getMatches(record,
                stratArgs.overlapFraction, stratArgs.numBreakpointOverlaps, stratArgs.numBreakpointOverlapsInterchrom);
        if (stratifications.size() > 1) {
            // don't allow more than one match since it would proliferate variants
            final String matchesString = String.join(", ", stratifications.stream().map(SVStratificationEngine.Stratum::getName).collect(Collectors.toList()));
            throw new GATKException("Record " + record.getId() + " matched multiple groups: " + matchesString +
                    ". Groups must be mutually exclusive. Please modify the group configurations and/or tracks so that " +
                    "no variant can match more than one group.");
        } else if (stratifications.isEmpty()) {
            // no match, don't cluster
            record.getAttributes().put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, Collections.singletonList(record.getId()));
            record.getAttributes().put(GATKSVVCFConstants.STRATUM_INFO_KEY, Collections.singletonList(SVStratify.DEFAULT_STRATUM));
            write(record);
        } else {
            // exactly one match
            final SVStratificationEngine.Stratum stratum = stratifications.iterator().next();
            Utils.validate(clusterEngineMap.containsKey(stratum.getName()), "Group undefined: " + stratum.getName());
            record.getAttributes().put(GATKSVVCFConstants.STRATUM_INFO_KEY, Collections.singletonList(stratum.getName()));
            clusterAndWrite(record, clusterEngineMap.get(stratum.getName()));
        }
    }

    private void clusterAndWrite(final SVCallRecord record, final SVClusterEngine clusterEngine) {
        clusterEngine.addAndFlush(record).stream().forEach(this::write);
    }
}
