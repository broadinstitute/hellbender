package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.vcf.*;

import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.*;

/**
 * This class contains the {@link VCFHeaderLine} definitions for the annotation keys in {@link GATKVCFConstants}.
 * VCF-standard header lines are in {@link VCFStandardHeaderLines}, in htsjdk
 */
public class GATKVCFHeaderLines {

    public static VCFInfoHeaderLine getInfoLine(final String id, final boolean lookInVCFStandardLines) {
        if (infoLines.containsKey(id)) {
            return infoLines.get(id);
        } else if (lookInVCFStandardLines) {
            return VCFStandardHeaderLines.getInfoLine(id);
        } else {
            throw new IllegalStateException("No VCF INFO header line found for key " + id);
        }
    }

    public static VCFInfoHeaderLine getInfoLine(final String id) {
        return getInfoLine(id, false);
    }

    public static VCFFormatHeaderLine getFormatLine(final String id, final boolean lookInVCFStandardLines) {
        if (formatLines.containsKey(id)) {
            return formatLines.get(id);
        } else if (lookInVCFStandardLines) {
            return VCFStandardHeaderLines.getFormatLine(id);
        } else {
            throw new IllegalStateException("No VCF INFO header line found for key " + id);
        }
    }

    public static VCFFormatHeaderLine getFormatLine(final String id) {
        return getFormatLine(id, false);
    }

    public static VCFFilterHeaderLine getFilterLine(final String id) {
        if (!filterLines.containsKey(id)) {
            throw new IllegalStateException("No VCF FILTER header line found for key " + id);
        }
        return filterLines.get(id);
    }

    public static Set<VCFInfoHeaderLine> getAllInfoLines() { return Collections.unmodifiableSet(new HashSet<>(infoLines.values())); }
    public static Set<VCFFormatHeaderLine> getAllFormatLines() { return Collections.unmodifiableSet(new HashSet<>(formatLines.values())); }
    public static Set<VCFFilterHeaderLine> getAllFilterLines() { return Collections.unmodifiableSet(new HashSet<>(filterLines.values())); }

    private static final Map<String, VCFInfoHeaderLine> infoLines = new LinkedHashMap<>(60);
    private static final Map<String, VCFFormatHeaderLine> formatLines = new LinkedHashMap<>(25);
    private static final Map<String, VCFFilterHeaderLine> filterLines = new LinkedHashMap<>(2);

    private static void addFormatLine(final VCFFormatHeaderLine line) {
        Utils.nonNull(line);
        formatLines.put(line.getID(), line);
    }

    private static void addInfoLine(final VCFInfoHeaderLine line) {
        Utils.nonNull(line);
        infoLines.put(line.getID(), line);
    }

    private static void addFilterLine(final VCFFilterHeaderLine line) {
        Utils.nonNull(line);
        filterLines.put(line.getID(), line);
    }

    public static VCFFormatHeaderLine getEquivalentFormatHeaderLine(final String infoFieldKey) {
        final VCFInfoHeaderLine infoLine = getInfoLine(infoFieldKey);
        if (infoLine.isFixedCount()) {
            return new VCFFormatHeaderLine(infoLine.getID(), infoLine.getCount(), infoLine.getType(), infoLine.getDescription());
        } else {
            return new VCFFormatHeaderLine(infoLine.getID(), infoLine.getCountType(), infoLine.getType(), infoLine.getDescription());
        }
    }

    static {
        addInfoLine(new VCFInfoHeaderLine(SB_TABLE_KEY, 4, VCFHeaderLineType.Integer, "Forward/reverse read counts for strand bias tests"));
        addInfoLine(new VCFInfoHeaderLine(GATKVCFConstants.AC_ADJUSTED_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Allele count for each ALT, adjusted to represent high quality genotypes only"));

        addFilterLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER_NAME, "Low quality"));
        addFilterLine(new VCFFilterHeaderLine(VCFConstants.PASSES_FILTERS_v4, "Site contains at least one allele that passes filters"));

        addFilterLine(new VCFFilterHeaderLine(NAY_FROM_YNG, "Considered a NAY in the Yay, Nay, Grey table"));
        addFilterLine(new VCFFilterHeaderLine(EXCESS_ALLELES,"Site has an excess of alternate alleles based on the input threshold"));
        addFilterLine(new VCFFilterHeaderLine(NO_HQ_GENOTYPES, "Site has no high quality variant genotypes"));
        addFilterLine(new VCFFilterHeaderLine(EXCESS_HET_KEY, "Site has excess het value larger than the threshold"));

        // M2-related filters
        addFilterLine(new VCFFilterHeaderLine(ALIGNMENT_ARTIFACT_FILTER_NAME, "Alignment artifact"));
        addFilterLine(new VCFFilterHeaderLine(CLUSTERED_EVENTS_FILTER_NAME, "Clustered events observed in the tumor"));
        addFilterLine(new VCFFilterHeaderLine(GERMLINE_RISK_FILTER_NAME, "Evidence indicates this site is germline, not somatic"));
        addFilterLine(new VCFFilterHeaderLine(PON_FILTER_NAME, "Blacklisted site in panel of normals"));
        addFilterLine(new VCFFilterHeaderLine(TUMOR_EVIDENCE_FILTER_NAME, "Mutation does not meet likelihood threshold"));
        addFilterLine(new VCFFilterHeaderLine(POLYMERASE_SLIPPAGE, "Site filtered due to contraction of short tandem repeat region"));
        addFilterLine(new VCFFilterHeaderLine(MULTIALLELIC_FILTER_NAME, "Site filtered because too many alt alleles pass tumor LOD"));
        addFilterLine(new VCFFilterHeaderLine(STRAND_ARTIFACT_FILTER_NAME, "Evidence for alt allele comes from one read direction only"));
        addFilterLine(new VCFFilterHeaderLine(ARTIFACT_IN_NORMAL_FILTER_NAME, "artifact_in_normal"));
        addFilterLine(new VCFFilterHeaderLine(MEDIAN_BASE_QUALITY_FILTER_NAME, "alt median base quality"));
        addFilterLine(new VCFFilterHeaderLine(MEDIAN_MAPPING_QUALITY_FILTER_NAME, "ref - alt median mapping quality"));
        addFilterLine(new VCFFilterHeaderLine(MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME, "abs(ref - alt) median fragment length"));
        addFilterLine(new VCFFilterHeaderLine(READ_POSITION_FILTER_NAME, "median distance of alt variants from end of reads"));
        addFilterLine(new VCFFilterHeaderLine(CONTAMINATION_FILTER_NAME, "contamination"));
        addFilterLine(new VCFFilterHeaderLine(DUPLICATED_EVIDENCE_FILTER_NAME, "evidence for alt allele is overrepresented by apparent duplicates"));
        addFilterLine(new VCFFilterHeaderLine(READ_ORIENTATION_ARTIFACT_FILTER_NAME, "orientation bias detected by the orientation bias mixture model"));
        addFilterLine(new VCFFilterHeaderLine(BAD_HAPLOTYPE_FILTER_NAME, "Variant near filtered variant on same haplotype."));
        addFilterLine(new VCFFilterHeaderLine(STRICT_STRAND_BIAS_FILTER_NAME, "Evidence for alt allele is not represented in both directions"));
        addFilterLine(new VCFFilterHeaderLine(N_RATIO_FILTER_NAME, "Ratio of N to alt exceeds specified ratio"));
        addFilterLine(new VCFFilterHeaderLine(ALLELE_FRACTION_FILTER_NAME, "Allele fraction is below specified threshold"));

        //Mitochondrial M2-related filters
        addFilterLine(new VCFFilterHeaderLine(POSSIBLE_NUMT_FILTER_NAME, "Allele depth is below expected coverage of NuMT in autosome"));
        addFilterLine(new VCFFilterHeaderLine(LOW_HET_FILTER_NAME, "All low heteroplasmy sites are filtered when at least x low het sites pass all other filters"));
        addFilterLine(new VCFFilterHeaderLine(FAIL, "Fail the site if all alleles fail but for different reasons."));
        addFilterLine(new VCFFilterHeaderLine(SITE_LEVEL_FILTERS, "There are no allele specific filters that apply to this allele. Only site level filters apply."));
        addFilterLine(new VCFFilterHeaderLine(LOW_HET_FILTER_NAME, "All low heteroplasmy sites are filtered when at least x low het sites pass all other filters"));

        addFormatLine(new VCFFormatHeaderLine(ALLELE_BALANCE_KEY, 1, VCFHeaderLineType.Float, "Allele balance for each het genotype"));
        addFormatLine(new VCFFormatHeaderLine(MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY, 1, VCFHeaderLineType.Integer, "Number of Mapping Quality Zero Reads per sample"));
        addFormatLine(new VCFFormatHeaderLine(STRAND_COUNT_BY_SAMPLE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of reads on the forward and reverse strand supporting each allele (including reference)"));
        addFormatLine(new VCFFormatHeaderLine(STRAND_BIAS_BY_SAMPLE_KEY, 4, VCFHeaderLineType.Integer, "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias."));
        addFormatLine(new VCFFormatHeaderLine(HAPLOTYPE_CALLER_PHASING_ID_KEY, 1, VCFHeaderLineType.String, "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"));
        addFormatLine(new VCFFormatHeaderLine(HAPLOTYPE_CALLER_PHASING_GT_KEY, 1, VCFHeaderLineType.String, "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles"));
        addFormatLine(new VCFFormatHeaderLine(PHASE_SET_KEY, 1, VCFHeaderLineType.Integer, "Phasing set (typically the position of the first variant in the set)"));
        addFormatLine(new VCFFormatHeaderLine(MIN_DP_FORMAT_KEY, 1, VCFHeaderLineType.Integer, "Minimum DP observed within the GVCF block"));
        addFormatLine(new VCFFormatHeaderLine(REFERENCE_GENOTYPE_QUALITY, 1, VCFHeaderLineType.Integer, "Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)"));
        addFormatLine(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_KEY, 1, VCFHeaderLineType.Integer, "Phred score of the genotype combination and phase given that the genotypes are correct"));
        addFormatLine(new VCFFormatHeaderLine(PHRED_SCALED_POSTERIORS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Phred-scaled Posterior Genotype Probabilities"));
        addFormatLine(new VCFFormatHeaderLine(JOINT_LIKELIHOOD_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred-scaled joint likelihood of the genotype combination (before applying family priors)"));
        addFormatLine(new VCFFormatHeaderLine(JOINT_POSTERIOR_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred-scaled joint posterior probability of the genotype combination (after applying family priors)"));

        // M2-related format lines
        addFormatLine(new VCFFormatHeaderLine(ALLELE_FRACTION_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele fractions of alternate alleles in the tumor"));
        addFormatLine(new VCFFormatHeaderLine(F1R2_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Count of reads in F1R2 pair orientation supporting each allele"));
        addFormatLine(new VCFFormatHeaderLine(F2R1_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Count of reads in F2R1 pair orientation supporting each allele"));
        addFormatLine(new VCFFormatHeaderLine(FEATURIZED_READ_SETS_KEY, 1, VCFHeaderLineType.Integer, "Featurized read sets by allele for Mutect3 training"));
        addFormatLine(new VCFFormatHeaderLine(FRAGMENT_ALLELE_DEPTHS, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Count of fragments supporting each allele."));

        addFormatLine(new VCFFormatHeaderLine(FEATURIZED_READ_SETS_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Featurized read sets by allele for Mutect3 training"));
        addFormatLine(new VCFFormatHeaderLine(FEATURIZED_READ_SETS_COUNTS_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Counts by allele of reads contributing to the featurized read sets annotation"));
        addInfoLine(new VCFInfoHeaderLine(MLE_ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed"));
        addInfoLine(new VCFInfoHeaderLine(MLE_ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"));
        addInfoLine(new VCFInfoHeaderLine(DOWNSAMPLED_KEY, 0, VCFHeaderLineType.Flag, "Were any of the samples downsampled?"));
        addInfoLine(new VCFInfoHeaderLine(BASE_QUAL_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities"));
        addInfoLine(new VCFInfoHeaderLine(AS_BASE_QUAL_RANK_SUM_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "allele specific Z-score from Wilcoxon rank sum test of each Alt Vs. Ref base qualities"));
        addInfoLine(new VCFInfoHeaderLine(AS_RAW_BASE_QUAL_RANK_SUM_KEY, 1, VCFHeaderLineType.String, "raw data for allele specific rank sum test of base qualities"));
        addInfoLine(new VCFInfoHeaderLine(CLIPPING_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases"));
        addInfoLine(new VCFInfoHeaderLine(FISHER_STRAND_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled p-value using Fisher's exact test to detect strand bias"));
        addInfoLine(new VCFInfoHeaderLine(AS_FISHER_STRAND_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "allele specific phred-scaled p-value using Fisher's exact test to detect strand bias of each alt allele"));
        addInfoLine(new VCFInfoHeaderLine(AS_SB_TABLE_KEY, 1, VCFHeaderLineType.String, "Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |."));
        addInfoLine(new VCFInfoHeaderLine(AS_SBP_TABLE_KEY, 1, VCFHeaderLineType.String, "Allele-specific forward/reverse read counts for strand bias tests (with probability). Includes the reference and alleles separated by |."));
        addInfoLine(new VCFInfoHeaderLine(NOCALL_CHROM_KEY, 1, VCFHeaderLineType.Integer, "Number of no-called samples"));
        addInfoLine(new VCFInfoHeaderLine(GQ_MEAN_KEY, 1, VCFHeaderLineType.Float, "Mean of all GQ values"));
        addInfoLine(new VCFInfoHeaderLine(GQ_STDEV_KEY, 1, VCFHeaderLineType.Float, "Standard deviation of all GQ values"));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPE_SCORE_KEY, 1, VCFHeaderLineType.Float, "Consistency of the site with at most two segregating haplotypes"));
        addInfoLine(new VCFInfoHeaderLine(INBREEDING_COEFFICIENT_KEY, 1, VCFHeaderLineType.Float, "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"));
        addInfoLine(new VCFInfoHeaderLine(AS_INBREEDING_COEFFICIENT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele-specific inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"));
        addInfoLine(new VCFInfoHeaderLine(EXCESS_HET_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled p-value for exact test of excess heterozygosity"));
        addInfoLine(new VCFInfoHeaderLine(RAW_GENOTYPE_COUNT_KEY, 3, VCFHeaderLineType.Integer, "Counts of genotypes w.r.t. the reference allele in the following order: 0/0, 0/*, */*, i.e. all alts lumped together; for use in calculating excess heterozygosity"));
        addInfoLine(new VCFInfoHeaderLine(LIKELIHOOD_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref haplotype likelihoods"));
        addInfoLine(new VCFInfoHeaderLine(MAP_QUAL_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"));
        addInfoLine(new VCFInfoHeaderLine(AS_MAP_QUAL_RANK_SUM_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "allele specific Z-score From Wilcoxon rank sum test of each Alt vs. Ref read mapping qualities"));
        addInfoLine(new VCFInfoHeaderLine(MAPPING_QUALITY_DEPTH_DEPRECATED, 1, VCFHeaderLineType.Integer, "Depth over variant samples for better MQ calculation (deprecated -- use " + RAW_MAPPING_QUALITY_WITH_DEPTH_KEY + " instead.)"));
        addInfoLine(new VCFInfoHeaderLine(RAW_RMS_MAPPING_QUALITY_DEPRECATED, 1, VCFHeaderLineType.Float, "Raw data for RMS Mapping Quality (deprecated -- use " + RAW_MAPPING_QUALITY_WITH_DEPTH_KEY + " instead.)"));
        addInfoLine(new VCFInfoHeaderLine(RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, 2, VCFHeaderLineType.Integer, "Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated " + RMSMappingQuality.getDeprecatedRawKeyName() + " formulation."));
        addInfoLine(new VCFInfoHeaderLine(AS_RAW_RMS_MAPPING_QUALITY_KEY, 1, VCFHeaderLineType.String, "Allele-specfic raw data for RMS Mapping Quality"));
        addInfoLine(new VCFInfoHeaderLine(AS_RMS_MAPPING_QUALITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele-specific RMS Mapping Quality"));
        addInfoLine(new VCFInfoHeaderLine(AS_RAW_MAP_QUAL_RANK_SUM_KEY, 1, VCFHeaderLineType.String, "Allele-specfic raw data for Mapping Quality Rank Sum"));
        addInfoLine(new VCFInfoHeaderLine(AS_MAP_QUAL_RANK_SUM_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele-specific Mapping Quality Rank Sum"));
        addInfoLine(new VCFInfoHeaderLine(AS_FILTER_STATUS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Filter status for each allele, as assessed by ApplyVQSR. Note that the VCF filter field will reflect the most lenient/sensitive status across all alleles."));
        addInfoLine(new VCFInfoHeaderLine(AS_CULPRIT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each alt allele, the annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out"));
        addInfoLine(new VCFInfoHeaderLine(AS_VQS_LOD_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each alt allele, the log odds of being a true variant versus being false under the trained gaussian mixture model"));
        addInfoLine(new VCFInfoHeaderLine(AS_YNG_STATUS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each alt allele, the yay/nay/grey status (yay are known good alleles, nay are known false positives, grey are unknown)"));
        addInfoLine(new VCFInfoHeaderLine(SCORE_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Score according to the model applied by ScoreVariantAnnotations"));
        addInfoLine(new VCFInfoHeaderLine(CALIBRATION_SENSITIVITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, String.format("Calibration sensitivity corresponding to the value of %s", SCORE_KEY)));

        addInfoLine(new VCFInfoHeaderLine(TRANSMITTED_SINGLETON, 1, VCFHeaderLineType.String, "Possible transmitted singleton (site with AC=2 from parent and child). Parent ID is listed."));
        addInfoLine(new VCFInfoHeaderLine(NON_TRANSMITTED_SINGLETON, 1, VCFHeaderLineType.String, "Possible non transmitted singleton (site with AC=1 in just one parent). Parent ID is listed."));
        addInfoLine(new VCFInfoHeaderLine(HI_CONF_DENOVO_KEY, 1, VCFHeaderLineType.String, "High confidence possible de novo mutation (GQ >= 20 for all trio members)=[comma-delimited list of child samples]"));
        addInfoLine(new VCFInfoHeaderLine(LO_CONF_DENOVO_KEY, 1, VCFHeaderLineType.String, "Low confidence possible de novo mutation (GQ >= 10 for child, GQ > 0 for parents)=[comma-delimited list of child samples]"));
        addInfoLine(new VCFInfoHeaderLine(QUAL_BY_DEPTH_KEY, 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth"));
        addInfoLine(new VCFInfoHeaderLine(AS_QUAL_BY_DEPTH_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele-specific Variant Confidence/Quality by Depth"));
        addInfoLine(new VCFInfoHeaderLine(AS_QUAL_KEY, 1, VCFHeaderLineType.Float, "Allele-specific Variant Qual Score"));
        addInfoLine(new VCFInfoHeaderLine(AS_RAW_QUAL_APPROX_KEY, 1, VCFHeaderLineType.String, "Allele-specific QUAL approximations"));
        addInfoLine(new VCFInfoHeaderLine(RAW_QUAL_APPROX_KEY, 1, VCFHeaderLineType.Integer, "Sum of PL[0] values; used to approximate the QUAL score"));
        addInfoLine(new VCFInfoHeaderLine(VARIANT_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "(informative) depth over variant genotypes"));
        addInfoLine(new VCFInfoHeaderLine(AS_VARIANT_DEPTH_KEY, 1, VCFHeaderLineType.String, "Allele-specific (informative) depth over variant genotypes -- including ref, RAW format"));
        addInfoLine(new VCFInfoHeaderLine(AS_ALT_ALLELE_DEPTH_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Allele-specific (informative) depth for alt alleles over variant genotypes; effectively sum of ADs"));
        addInfoLine(new VCFInfoHeaderLine(READ_POS_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"));
        addInfoLine(new VCFInfoHeaderLine(AS_READ_POS_RANK_SUM_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "allele specific Z-score from Wilcoxon rank sum test of each Alt vs. Ref read position bias"));
        addInfoLine(new VCFInfoHeaderLine(AS_RAW_READ_POS_RANK_SUM_KEY, 1, VCFHeaderLineType.String, "allele specific raw data for rank sum test of read position bias"));
        addInfoLine(new VCFInfoHeaderLine(SAMPLE_LIST_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "List of polymorphic samples"));
        addInfoLine(new VCFInfoHeaderLine(STRAND_ODDS_RATIO_KEY, 1, VCFHeaderLineType.Float, "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias"));
        addInfoLine(new VCFInfoHeaderLine(AS_STRAND_ODDS_RATIO_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele specific strand Odds Ratio of 2x|Alts| contingency table to detect allele specific strand bias"));
        addInfoLine(new VCFInfoHeaderLine(AS_STRAND_ODDS_RATIO_PROB_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele specific strand Odds Ratio of 2x|Alts| contingency table to detect allele specific strand bias (using probability)"));
        addInfoLine(new VCFInfoHeaderLine(STR_PRESENT_KEY, 0, VCFHeaderLineType.Flag, "Variant is a short tandem repeat"));
        addInfoLine(new VCFInfoHeaderLine(REPEAT_UNIT_KEY, 1, VCFHeaderLineType.String, "Tandem repeat unit (bases)"));
        addInfoLine(new VCFInfoHeaderLine(REPEATS_PER_ALLELE_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Number of times tandem repeat unit is repeated, for each allele (including reference)"));
        addInfoLine(new VCFInfoHeaderLine(NUMBER_OF_DISCOVERED_ALLELES_KEY, 1, VCFHeaderLineType.Integer, "Number of alternate alleles discovered (but not necessarily genotyped) at this site"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_AC_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Original AC"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_AF_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Original AF"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_AN_KEY, 1, VCFHeaderLineType.Integer, "Original AN"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_DP_KEY, 1, VCFHeaderLineType.Integer, "Original DP"));
        addInfoLine(new VCFInfoHeaderLine(VQS_LOD_KEY, 1, VCFHeaderLineType.Float, "Log odds of being a true variant versus being false under the trained gaussian mixture model"));
        addInfoLine(new VCFInfoHeaderLine(CNN_1D_KEY, 1, VCFHeaderLineType.Float, "Log odds of being a true variant versus being false under the trained 1D Convolutional Neural Network"));
        addInfoLine(new VCFInfoHeaderLine(CNN_2D_KEY, 1, VCFHeaderLineType.Float, "Log odds of being a true variant versus being false under the trained 2D Convolutional Neural Network"));
        addInfoLine(new VCFInfoHeaderLine(CULPRIT_KEY, 1, VCFHeaderLineType.String, "The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out"));
        addInfoLine(new VCFInfoHeaderLine(POSITIVE_LABEL_KEY, 1, VCFHeaderLineType.Flag, "This variant was used to build the positive training set of good variants"));
        addInfoLine(new VCFInfoHeaderLine(NEGATIVE_LABEL_KEY, 1, VCFHeaderLineType.Flag, "This variant was used to build the negative training set of bad variants"));
        addInfoLine(new VCFInfoHeaderLine(GENOTYPE_AND_VALIDATE_STATUS_KEY, 1, VCFHeaderLineType.String, "Value from the validation VCF"));
        addInfoLine(new VCFInfoHeaderLine(INTERVAL_GC_CONTENT_KEY, 1, VCFHeaderLineType.Float, "GC Content of the interval"));
        addInfoLine(new VCFInfoHeaderLine(GENOTYPE_PRIOR_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Genotype Likelihood Prior"));
        addInfoLine(new VCFInfoHeaderLine(BASE_QUAL_HISTOGRAM_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer,
                "Base quality counts for each allele represented sparsely as alternating entries of qualities and counts for each allele." +
                        "For example [10,1,0,20,0,1] means one ref base with quality 10 and one alt base with quality 20."));

        addInfoLine(new VCFInfoHeaderLine(TREE_SCORE, 1, VCFHeaderLineType.Float, "Score from single sample filtering with random forest model."));

        // M2-related info lines
        addInfoLine(new VCFInfoHeaderLine(EVENT_COUNT_IN_HAPLOTYPE_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of somatic events in best supporting haplotype for each alt allele"));
        addInfoLine(new VCFInfoHeaderLine(EVENT_COUNT_IN_REGION_KEY, 1, VCFHeaderLineType.Integer, "Number of potential somatic events in the assembly region"));
        addInfoLine(new VCFInfoHeaderLine(NORMAL_LOG_10_ODDS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Normal log 10 likelihood ratio of diploid het or hom alt genotypes"));
        addInfoLine(new VCFInfoHeaderLine(TUMOR_LOG_10_ODDS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Log 10 likelihood ratio score of variant existing versus not existing"));
        addFormatLine(new VCFFormatHeaderLine(TUMOR_LOG_10_ODDS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Log 10 likelihood ratio score of variant existing versus not existing"));
        addInfoLine(new VCFInfoHeaderLine(IN_PON_KEY, 0, VCFHeaderLineType.Flag, "site found in panel of normals"));
        addInfoLine(new VCFInfoHeaderLine(POPULATION_AF_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "negative log 10 population allele frequencies of alt alleles"));
        addInfoLine(new VCFInfoHeaderLine(GERMLINE_QUAL_KEY, 1, VCFHeaderLineType.Integer, "Phred-scaled quality that alt alleles are not germline variants"));
        addInfoLine(new VCFInfoHeaderLine(SEQUENCING_QUAL_KEY, 1, VCFHeaderLineType.Integer, "Phred-scaled quality that alt alleles are not sequencing errors"));
        addInfoLine(new VCFInfoHeaderLine(POLYMERASE_SLIPPAGE_QUAL_KEY, 1, VCFHeaderLineType.Integer, "Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors"));
        addInfoLine(new VCFInfoHeaderLine(STRAND_QUAL_KEY, 1, VCFHeaderLineType.Integer, "Phred-scaled quality of strand bias artifact"));
        addInfoLine(new VCFInfoHeaderLine(CONTAMINATION_QUAL_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled qualities that alt allele are not due to contamination"));
        addInfoLine(new VCFInfoHeaderLine(READ_ORIENTATION_QUAL_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled qualities that alt allele are not due to read orientation artifact"));
        addInfoLine(new VCFInfoHeaderLine(NORMAL_ARTIFACT_LOG_10_ODDS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Log 10 odds of artifact in normal with same allele fraction as tumor"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_CONTIG_MISMATCH_KEY, 1, VCFHeaderLineType.Integer, "Number of alt reads whose original alignment doesn't match the current contig."));
        addInfoLine(new VCFInfoHeaderLine(N_COUNT_KEY, 1, VCFHeaderLineType.Integer, "Count of N bases in the pileup"));
        addInfoLine(new VCFInfoHeaderLine(AS_UNIQUE_ALT_READ_SET_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of reads with unique start and mate end positions for each alt at a variant site"));
        addInfoLine(new VCFInfoHeaderLine(MEDIAN_BASE_QUALITY_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "median base quality by allele"));
        addInfoLine(new VCFInfoHeaderLine(MEDIAN_FRAGMENT_LENGTH_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "median fragment length by allele"));
        addInfoLine(new VCFInfoHeaderLine(MEDIAN_MAPPING_QUALITY_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "median mapping quality by allele"));
        addInfoLine(new VCFInfoHeaderLine(MEDIAN_READ_POSITON_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "median distance from end of read"));
        addInfoLine(new VCFInfoHeaderLine(AS_SB_TABLE_KEY, 1, VCFHeaderLineType.String, "Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |."));
        addInfoLine(new VCFInfoHeaderLine(UNITIG_SIZES_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Sizes of reassembled unitigs"));
        addInfoLine(new VCFInfoHeaderLine(JOINT_ALIGNMENT_COUNT_KEY, 1, VCFHeaderLineType.Integer, "Number of joint alignments"));
        addInfoLine(new VCFInfoHeaderLine(ALIGNMENT_SCORE_DIFFERENCE_KEY, 1, VCFHeaderLineType.Integer, "Difference in alignment score between best and next-best alignment"));
        addInfoLine(new VCFInfoHeaderLine(EXT_COLLAPSED_KEY,1, VCFHeaderLineType.Integer, "Indicates longer hmer collapsing took place (this is a flow-based specific tag)"));
        addInfoLine(new VCFInfoHeaderLine(POSSIBLE_FP_ADJACENT_TP_KEY,1, VCFHeaderLineType.Flag, "Indicates a locus where false positive allele might be affecting a true positive allele"));

        // flow annotations
        addInfoLine(new VCFInfoHeaderLine(FLOW_INDEL_CLASSIFY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Flow: indel class: ins, del, NA"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_INDEL_LENGTH, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Flow: length of indel"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_HMER_INDEL_LENGTH, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Flow: length of the hmer indel, if so"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_HMER_INDEL_NUC, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Flow: nucleotide of the hmer indel, if so"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_LEFT_MOTIF, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Flow: motif to the left of the indel"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_RIGHT_MOTIF, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Flow: motif to the right of the indel"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_GC_CONTENT, 1, VCFHeaderLineType.Float, "Flow: percentage of G or C in the window around hmer"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_CYCLESKIP_STATUS, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Flow: cycle skip status: cycle-skip, possible-cycle-skip, non-skip"));
        addInfoLine(new VCFInfoHeaderLine(FLOW_VARIANT_TYPE, 1, VCFHeaderLineType.String, "Flow: type of variant: SNP/NON-H-INDEL/H-INDEL"));
        addInfoLine(new VCFInfoHeaderLine(REFERENCE_BASES_KEY, 1, VCFHeaderLineType.String, "local reference bases."));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPE_EQUIVALENCE_COUNTS_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Counts of support for haplotype groups excluding difference at the site in question."));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPE_COMPLEXITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Edit distances of each alt allele's most common supporting haplotype from closest germline haplotype, excluding differences at the site in question."));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPE_DOMINANCE_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "For each alt allele, fraction of read support that best fits the most-supported haplotype containing the allele"));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPES_BEFORE_FILTERING_KEY, 1, VCFHeaderLineType.Integer, "Haplotypes detected by the assembly region before haplotype filtering is applied"));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPES_FILTERED_KEY, 1, VCFHeaderLineType.Integer, "Haplotypes filtered out by the haplotype filtering code"));
    }
}
