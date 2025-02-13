package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.collect.Lists;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.io.Serializable;
import java.util.List;

/**
 * Set of arguments related to the {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler}
 */
public abstract class ReadThreadingAssemblerArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final double DEFAULT_PRUNING_LOG_ODDS_THRESHOLD = MathUtils.log10ToLog(1.0);
    public static final double DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD = MathUtils.log10ToLog(4.0);

    public static final String ERROR_CORRECT_READS_LONG_NAME = "error-correct-reads";
    public static final String PILEUP_ERROR_CORRECTION_LOG_ODDS_LONG_NAME = "error-correction-log-odds";
    public static final String CAPTURE_ASSEMBLY_FAILURE_BAM_LONG_NAME = "capture-assembly-failure-bam";
    public static final String KMER_SIZE_LONG_NAME = "kmer-size";
    public static final String DONT_INCREASE_KMER_SIZE_LONG_NAME = "dont-increase-kmer-sizes-for-cycles";
    public static final String LINKED_DE_BRUIJN_GRAPH_LONG_NAME = "linked-de-bruijn-graph";
    public static final String ALLOW_NON_UNIQUE_KMERS_IN_REF_LONG_NAME = "allow-non-unique-kmers-in-ref";
    public static final String NUM_PRUNING_SAMPLES_LONG_NAME = "num-pruning-samples";
    public static final String MIN_DANGLING_BRANCH_LENGTH_LONG_NAME = "min-dangling-branch-length";
    public static final String RECOVER_ALL_DANGLING_BRANCHES_LONG_NAME = "recover-all-dangling-branches";
    public static final String MAX_NUM_HAPLOTYPES_IN_POPULATION_LONG_NAME = "max-num-haplotypes-in-population";
    public static final String MIN_PRUNING_LONG_NAME = "min-pruning";
    public static final String ADAPTIVE_PRUNING_INITIAL_ERROR_RATE_LONG_NAME = "adaptive-pruning-initial-error-rate";
    public static final String PRUNING_LOD_THRESHOLD_LONG_NAME = "pruning-lod-threshold";
    public static final String PRUNING_SEEDING_LOD_THRESHOLD_LONG_NAME = "pruning-seeding-lod-threshold";
    public static final String MAX_UNPRUNED_VARIANTS_LONG_NAME = "max-unpruned-variants";
    public static final String DISABLE_ARTIFICIAL_HAPLOTYPE_RECOVERY_LONG_NAME = "disable-artificial-haplotype-recovery";
    public static final String ENABLE_LEGACY_GRAPH_CYCLE_DETECTION_LONG_NAME = "enable-legacy-graph-cycle-detection";
    public static final String DEBUG_ASSEMBLY_LONG_NAME = "debug-assembly";
    public static final String DEBUG_GRAPH_TRANSFORMATIONS_LONG_NAME = "debug-graph-transformations";
    public static final String DEBUG_ASSEMBLY_VARIANTS_OUT_LONG_NAME = "debug-assembly-variants-out";
    public static final String GRAPH_OUTPUT_LONG_NAME = "graph-output";
    public static final String HAPLOTYPE_DEBUG_HISTOGRAM_OUTPUT_LONG_NAME = "haplotype-debug-histogram-output";
    public static final String NUM_MATCHING_BASES_IN_DANGLING_END_TO_RECOVER_LONG_NAME = "num-matching-bases-in-dangling-end-to-recover";
    public static final String KMER_LENGTH_FOR_READ_ERROR_CORRECTION_LONG_NAME = "kmer-length-for-read-error-correction";
    public static final String MIN_OBSERVATIONS_FOR_KMER_TO_BE_SOLID_LONG_NAME = "min-observations-for-kmer-to-be-solid";



    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the read threading assembler
    // -----------------------------------------------------------------------------------------------

    /**
     * Multiple kmer sizes can be specified, using e.g. `--kmer-size 10 --kmer-size 25`.
     */
    @Advanced
    @Argument(fullName= KMER_SIZE_LONG_NAME, doc="Kmer size to use in the read threading assembler", optional = true)
    public List<Integer> kmerSizes = Lists.newArrayList(45);

    /**
     * When graph cycles are detected, the normal behavior is to increase kmer sizes iteratively until the cycles are
     * resolved. Disabling this behavior may cause the program to give up on assembling the ActiveRegion.
     */
    @Advanced
    @Argument(fullName= DONT_INCREASE_KMER_SIZE_LONG_NAME, doc="Disable iterating over kmer sizes when graph cycles are detected", optional = true)
    public boolean dontIncreaseKmerSizesForCycles = false;

    /**
     * By default, the program does not allow processing of reference sections that contain non-unique kmers. Disabling
     * this check may cause problems in the assembly graph.
     */
    @Advanced
    @Argument(fullName= ALLOW_NON_UNIQUE_KMERS_IN_REF_LONG_NAME, doc="Allow graphs that have non-unique kmers in the reference", optional = true)
    public boolean allowNonUniqueKmersInRef = false;

    /**
     * If fewer samples than the specified number pass the minPruning threshold for a given path, that path will be eliminated from the graph.
     */
    @Advanced
    @Argument(fullName= NUM_PRUNING_SAMPLES_LONG_NAME, doc="Number of samples that must pass the minPruning threshold", optional = true)
    public int numPruningSamples = 1;

    /**
     * When constructing the assembly graph we are often left with "dangling" branches.  The assembly engine attempts to rescue these branches
     * by merging them back into the main graph.  This argument describes the minimum length of a dangling branch needed for the engine to
     * try to rescue it.  A smaller number here will lead to higher sensitivity to real variation but also to a higher number of false positives.
     */
    @Advanced
    @Argument(fullName= MIN_DANGLING_BRANCH_LENGTH_LONG_NAME, doc="Minimum length of a dangling branch to attempt recovery", optional = true)
    public int minDanglingBranchLength = 4;

    /**
     * By default, the read threading assembler does not recover dangling branches that fork after splitting from the reference.  This argument
     * tells the assembly engine to recover all dangling branches.
     */
    @Advanced
    @Argument(fullName= RECOVER_ALL_DANGLING_BRANCHES_LONG_NAME, doc="Recover all dangling branches", optional = true)
    public boolean recoverAllDanglingBranches = false;

    /**
     * The assembly graph can be quite complex, and could imply a very large number of possible haplotypes.  Each haplotype
     * considered requires N PairHMM evaluations if there are N reads across all samples.  In order to control the
     * run of the haplotype caller we only take maxNumHaplotypesInPopulation paths from the graph, in order of their
     * weights, no matter how many paths are possible to generate from the graph.  Putting this number too low
     * will result in dropping true variation because paths that include the real variant are not even considered.
     * You can consider increasing this number when calling organisms with high heterozygosity.
     */
    @Advanced
    @Argument(fullName= MAX_NUM_HAPLOTYPES_IN_POPULATION_LONG_NAME, doc="Maximum number of haplotypes to consider for your population", optional = true)
    public int maxNumHaplotypesInPopulation = 128;

    /**
     * Paths with fewer supporting kmers than the specified threshold will be pruned from the graph.
     *
     * Be aware that this argument can dramatically affect the results of variant calling and should only be used with great caution.
     * Using a prune factor of 1 (or below) will prevent any pruning from the graph, which is generally not ideal; it can make the
     * calling much slower and even less accurate (because it can prevent effective merging of "tails" in the graph).  Higher values
     * tend to make the calling much faster, but also lowers the sensitivity of the results (because it ultimately requires higher
     * depth to produce calls).
     */
    @Advanced
    @Argument(fullName= MIN_PRUNING_LONG_NAME, doc = "Minimum support to not prune paths in the graph", optional = true)
    public int minPruneFactor = 2;

    /**
     * Initial base error rate guess for the probabilistic adaptive pruning model.  Results are not very sensitive to this
     * parameter because it is only a starting point from which the algorithm discovers the true error rate.
     */
    @Advanced
    @Argument(fullName= ADAPTIVE_PRUNING_INITIAL_ERROR_RATE_LONG_NAME, doc = "Initial base error rate estimate for adaptive pruning", optional = true)
    public double initialErrorRateForPruning = 0.001;

    /**
     * Log-10 likelihood ratio threshold for adaptive pruning algorithm.
     */
    @Advanced
    @Argument(fullName= PRUNING_LOD_THRESHOLD_LONG_NAME, doc = "Ln likelihood ratio threshold for adaptive pruning algorithm", optional = true)
    public double pruningLogOddsThreshold = DEFAULT_PRUNING_LOG_ODDS_THRESHOLD;

    /**
     * Log-10 likelihood ratio threshold for adaptive pruning algorithm.
     */
    @Advanced
    @Argument(fullName= PRUNING_SEEDING_LOD_THRESHOLD_LONG_NAME, doc = "Ln likelihood ratio threshold for seeding subgraph of good variation in adaptive pruning algorithm", optional = true)
    public double pruningSeedingLogOddsThreshold = DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD;

    /**
     * The maximum number of variants in graph the adaptive pruner will allow
     */
    @Advanced
    @Argument(fullName= MAX_UNPRUNED_VARIANTS_LONG_NAME, doc = "Maximum number of variants in graph the adaptive pruner will allow", optional = true)
    public int maxUnprunedVariants = 100;

    /**
     * Disables graph simplification into a seq graph, opts to construct a proper De Bruijn graph with potential loops
     *
     * NOTE: --linked-de-bruijn-graph is currently an experimental feature that does not directly match with
     *        the regular HaplotypeCaller. Specifically the haplotype finding code does not perform correctly at complicated
     *        sites. Use this mode at your own risk.
     */
    @Advanced
    @Argument(fullName= LINKED_DE_BRUIJN_GRAPH_LONG_NAME, doc = "If enabled, the Assembly Engine will construct a Linked De Bruijn graph to recover better haplotypes", optional = true)
    public boolean useLinkedDeBruijnGraph = false;


    /**
     * This is used to disable the recovery of paths that were dropped in the graph based on the junction trees. Disabling this
     * will affect sensitivity but improve phasing and runtime somewhat.
     */
    @Hidden
    @Argument(fullName= DISABLE_ARTIFICIAL_HAPLOTYPE_RECOVERY_LONG_NAME, doc = "If in 'linked-de-bruijn-graph' mode, disable recovery of haplotypes based on graph edges that are not in junction trees", optional = true)
    public boolean disableArtificialHaplotypeRecovery = false;

    /**
     * This option exists purely to support concordance with DRAGEN-GATK, it is not recommended to enable this in most use cases.
     * Use this toggle to re-enable the GATK3 behavior of checking for graph cycles (and consequently throwing away the graph) before
     * purning low weight chains.
     */
    @Hidden
    @Argument(fullName= ENABLE_LEGACY_GRAPH_CYCLE_DETECTION_LONG_NAME, doc = "Use to revert the change to assembly graph code that moved pruning to before cycle detection")
    public boolean enableLegacyGraphCycleDetection = false;

    @Advanced
    @Argument(fullName= DEBUG_ASSEMBLY_LONG_NAME, shortName="debug", doc="Print out verbose debug information about each assembly region", optional = true)
    public boolean debugAssembly;

    @Hidden
    @Argument(fullName= DEBUG_GRAPH_TRANSFORMATIONS_LONG_NAME, doc="Write DOT formatted graph files out of the assembler for only this graph size", optional = true)
    public boolean debugGraphTransformations = false;
    @Hidden
    @Argument(fullName = DEBUG_ASSEMBLY_VARIANTS_OUT_LONG_NAME, doc="Write a VCF file containing all of the variants from the event map (pre-filter) for each active region", optional = true)
    public String debugAssemblyVariantsOut = null;
    /**
     * This argument is meant for debugging and is not immediately useful for normal analysis use.
     */
    @Argument(fullName= GRAPH_OUTPUT_LONG_NAME, shortName="graph", doc="Write debug assembly graph information to this file", optional = true)
    public String graphOutput = null;

    /**
     * This argument is meant for debugging and is not immediately useful for normal analysis use.
     */
    @Hidden
    @Argument(fullName= HAPLOTYPE_DEBUG_HISTOGRAM_OUTPUT_LONG_NAME, doc="Write debug assembly graph information to this file", optional = true)
    public String haplotypeHistogramOutput = null;

    @Hidden
    @Argument(fullName = CAPTURE_ASSEMBLY_FAILURE_BAM_LONG_NAME, doc = "Write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", optional = true)
    public boolean captureAssemblyFailureBAM = false;

    @Hidden
    @Argument(fullName = NUM_MATCHING_BASES_IN_DANGLING_END_TO_RECOVER_LONG_NAME, doc = "Sets the number of exactly matching bases in the suffix of a dangling tail and the prefix for a dangling head necessary in order to recover the path. (-1 indicates legacy behavior)", optional = true)
    public int minMatchingBasesToDanglingEndRecovery = -1;

    //---------------------------------------------------------------------------------------------------------------
    //
    // Read Error Corrector Related Parameters
    //
    // ---------------------------------------------------------------------------------------------------------------

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName = PILEUP_ERROR_CORRECTION_LOG_ODDS_LONG_NAME, doc = "Log odds threshold for pileup error correction.  Off by default", optional = true)
    public double pileupErrorCorrectionLogOdds = Double.NEGATIVE_INFINITY;


    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName = ERROR_CORRECT_READS_LONG_NAME, doc = "Use an exploratory algorithm to error correct the kmers used during assembly", optional = true)
    public boolean errorCorrectReads = false;

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName= KMER_LENGTH_FOR_READ_ERROR_CORRECTION_LONG_NAME, doc = "Use an exploratory algorithm to error correct the kmers used during assembly", optional = true)
    public int kmerLengthForReadErrorCorrection = 25;

    @Hidden
    @Argument(fullName= MIN_OBSERVATIONS_FOR_KMER_TO_BE_SOLID_LONG_NAME, doc = "A k-mer must be seen at least these times for it considered to be solid", optional = true)
    public int minObservationsForKmerToBeSolid = 20;

    public abstract ReadThreadingAssembler makeReadThreadingAssembler();
}
