package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFConstants;

import java.util.Arrays;
import java.util.List;

/**
 * This class contains any constants (primarily FORMAT/INFO keys) in VCF files used by the GATK.
 * Note that VCF-standard constants are in VCFConstants, in htsjdk.  Keys in header lines should
 * have matching entries in GATKVCFHeaderLines
 */
public final class GATKVCFConstants {

    public static final String CONTIG_ID_KEY =                      "ID";
    public static final String CONTIG_LENGTH_KEY =                  "length";
    public static final String ASSEMBLY_NAME_KEY =                  "assembly";

    //INFO keys
    public static final String ALLELE_SPECIFIC_PREFIX =             "AS_";
    public static final String AS_FILTER_STATUS_KEY =               "AS_FilterStatus";
    public static final String RAW_RMS_MAPPING_QUALITY_DEPRECATED =        "RAW_MQ";  //NOTE: this is deprecated in favor of the new RAW_MQandDP below
    public static final String MAPPING_QUALITY_DEPTH_DEPRECATED =   "MQ_DP";  //NOTE: this is deprecated in favor of the new RAW_MQandDP below
    public static final String RAW_MAPPING_QUALITY_WITH_DEPTH_KEY = "RAW_MQandDP";
    public static final String AS_RMS_MAPPING_QUALITY_KEY =         "AS_MQ";
    public static final String AS_RAW_RMS_MAPPING_QUALITY_KEY =     "AS_RAW_MQ";
    public static final String AS_CULPRIT_KEY =                     "AS_culprit";
    public static final String AS_VQS_LOD_KEY =                     "AS_VQSLOD";
    public static final String AS_YNG_STATUS_KEY =                  "AS_YNG";
    public static final String SCORE_KEY =                          "SCORE";
    public static final String CALIBRATION_SENSITIVITY_KEY =        "CALIBRATION_SENSITIVITY";
    public static final String ORIGINAL_AC_KEY =                    "AC_Orig"; //SelectVariants
    public static final String ORIGINAL_AF_KEY =                    "AF_Orig"; //SelectVariants
    public static final String ORIGINAL_AN_KEY =                    "AN_Orig"; //SelectVariants
    public static final String AC_ADJUSTED_KEY =                    "AC_adj"; //GnarlyGenotyper
    public static final String BASE_QUAL_RANK_SUM_KEY =             "BaseQRankSum";
    public static final String BASE_QUAL_HISTOGRAM_KEY =            "BQHIST";
    public static final String AS_BASE_QUAL_RANK_SUM_KEY =          "AS_BaseQRankSum";
    public static final String AS_RAW_BASE_QUAL_RANK_SUM_KEY =      "AS_RAW_BaseQRankSum";
    public static final String GENOTYPE_AND_VALIDATE_STATUS_KEY =   "callStatus";
    public static final String CLIPPING_RANK_SUM_KEY =              "ClippingRankSum";
    public static final String CULPRIT_KEY =                        "culprit";
    public static final String ORIGINAL_DP_KEY =                    "DP_Orig"; //SelectVariants
    public static final String DOWNSAMPLED_KEY =                    "DS";
    public static final String EVENT_COUNT_IN_REGION_KEY =          "ECNT"; //M2
    public static final String EVENT_COUNT_IN_HAPLOTYPE_KEY =       "ECNTH"; //M2
    public static final String FISHER_STRAND_KEY =                  "FS";
    public static final String AS_FISHER_STRAND_KEY =               "AS_FS";
    public static final String AS_SB_TABLE_KEY =                    "AS_SB_TABLE";
    public static final String AS_SBP_TABLE_KEY =                   "AS_SBP_TABLE";
    public static final String POSSIBLE_FP_ADJACENT_TP_KEY =        "SUSP_NOISY_ADJACENT_TP_VARIANT";
    public static final String SB_TABLE_KEY =                       "SB_TABLE";
    public static final String GQ_MEAN_KEY =                        "GQ_MEAN";
    public static final String GQ_STDEV_KEY =                       "GQ_STDDEV";
    public static final String HAPLOTYPE_SCORE_KEY =                "HaplotypeScore";
    public static final String HI_CONF_DENOVO_KEY =                 "hiConfDeNovo";
    public static final String TRANSMITTED_SINGLETON =              "transmittedSingleton";
    public static final String NON_TRANSMITTED_SINGLETON =          "nonTransmittedSingleton";
    public static final String INTERVAL_GC_CONTENT_KEY =            "IGC";
    public static final String INBREEDING_COEFFICIENT_KEY =         "InbreedingCoeff";
    public static final String AS_INBREEDING_COEFFICIENT_KEY =      "AS_InbreedingCoeff";
    public static final String EXCESS_HET_KEY =                     "ExcessHet";
    public static final String RAW_GENOTYPE_COUNT_KEY =             "RAW_GT_COUNT";
    public static final String LIKELIHOOD_RANK_SUM_KEY =            "LikelihoodRankSum";
    public static final String LO_CONF_DENOVO_KEY =                 "loConfDeNovo";
    public static final String MLE_ALLELE_COUNT_KEY =               "MLEAC";
    public static final String MLE_ALLELE_FREQUENCY_KEY =           "MLEAF";
    public static final String MAP_QUAL_RANK_SUM_KEY =              "MQRankSum";
    public static final String RAW_MAP_QUAL_RANK_SUM_KEY =          "RAW_MQRankSum";
    public static final String AS_MAP_QUAL_RANK_SUM_KEY =           "AS_MQRankSum";
    public static final String AS_RAW_MAP_QUAL_RANK_SUM_KEY =       "AS_RAW_MQRankSum";
    public static final String NOCALL_CHROM_KEY =                   "NCC";
    public static final String NUMBER_OF_DISCOVERED_ALLELES_KEY =   "NDA";
    public static final String NEGATIVE_LABEL_KEY =                 "NEGATIVE_TRAIN_SITE";
    public static final String GENOTYPE_PRIOR_KEY =                 "PG";
    public static final String POSITIVE_LABEL_KEY =                 "POSITIVE_TRAIN_SITE";
    public static final String QUAL_BY_DEPTH_KEY =                  "QD";
    public static final String AS_QUAL_BY_DEPTH_KEY =               "AS_QD";
    public static final String AS_QUAL_KEY =                        "AS_QUAL";
    public static final String RAW_QUAL_APPROX_KEY =                "QUALapprox";
    public static final String AS_RAW_QUAL_APPROX_KEY =             "AS_QUALapprox";
    public static final String VARIANT_DEPTH_KEY =                  "VarDP";
    public static final String AS_VARIANT_DEPTH_KEY =               "AS_VarDP";
    public static final String AS_ALT_ALLELE_DEPTH_KEY =            "AS_AltDP";
    public static final String READ_POS_RANK_SUM_KEY =              "ReadPosRankSum";
    public static final String AS_READ_POS_RANK_SUM_KEY =           "AS_ReadPosRankSum";
    public static final String AS_RAW_READ_POS_RANK_SUM_KEY =       "AS_RAW_ReadPosRankSum";
    public static final String REPEATS_PER_ALLELE_KEY =             "RPA";
    public static final String REPEAT_UNIT_KEY =                    "RU";
    public static final String SAMPLE_LIST_KEY =                    "Samples";
    public static final String STRAND_ODDS_RATIO_KEY =              "SOR";
    public static final String AS_STRAND_ODDS_RATIO_KEY =           "AS_SOR";
    public static final String AS_STRAND_ODDS_RATIO_PROB_KEY =      "AS_SORP";
    public static final String STR_PRESENT_KEY =                    "STR";
    public static final String VQS_LOD_KEY =                        "VQSLOD";
    public static final String CNN_1D_KEY =                         "CNN_1D";
    public static final String CNN_2D_KEY =                         "CNN_2D";
    public static final String F1R2_KEY =                           "F1R2";
    public static final String F2R1_KEY =                           "F2R1";

    public static final String TREE_SCORE =                         "TREE_SCORE";

    // Mutect2-specific INFO keys
    public static final String TUMOR_LOG_10_ODDS_KEY =              "TLOD";
    public static final String NORMAL_LOG_10_ODDS_KEY =             "NLOD";
    public static final String IN_PON_KEY =                         "PON";
    public static final String NORMAL_ARTIFACT_LOG_10_ODDS_KEY =    "NALOD";
    public static final String POPULATION_AF_KEY =                  "POPAF";
    public static final String GERMLINE_QUAL_KEY =                  "GERMQ";
    public static final String SEQUENCING_QUAL_KEY =                "SEQQ";
    public static final String POLYMERASE_SLIPPAGE_QUAL_KEY =       "STRQ";
    public static final String STRAND_QUAL_KEY =                    "STRANDQ";
    public static final String CONTAMINATION_QUAL_KEY =             "CONTQ";
    public static final String READ_ORIENTATION_QUAL_KEY =          "ROQ";
    public static final String ORIGINAL_CONTIG_MISMATCH_KEY =       "OCM";
    public static final String N_COUNT_KEY =                        "NCount";
    public static final String AS_UNIQUE_ALT_READ_SET_COUNT_KEY =   "AS_UNIQ_ALT_READ_COUNT";
    public static final String MEDIAN_BASE_QUALITY_KEY =            "MBQ";
    public static final String MEDIAN_MAPPING_QUALITY_KEY =         "MMQ";
    public static final String MEDIAN_FRAGMENT_LENGTH_KEY =         "MFRL";
    public static final String MEDIAN_READ_POSITON_KEY =            "MPOS";
    public static final String UNITIG_SIZES_KEY =                   "UNITIGS";
    public static final String ALIGNMENT_SCORE_DIFFERENCE_KEY =     "ALIGN_DIFF";
    public static final String JOINT_ALIGNMENT_COUNT_KEY =          "NALIGNS";
    public static final String REFERENCE_BASES_KEY =                "REF_BASES";

    // Methylation-specific INFO Keys
    public static final String UNCONVERTED_BASE_COVERAGE_KEY =      "UNCONVERTED_BASE_COV";
    public static final String CONVERTED_BASE_COVERAGE_KEY =        "CONVERTED_BASE_COV";
    public static final String METHYLATION_REFERENCE_CONTEXT_KEY =  "REFERENCE_CONTEXT";


    // FORMAT keys
    public static final String ALLELE_BALANCE_KEY =                 "AB";
    public static final String JOINT_LIKELIHOOD_TAG_NAME =          "JL"; //FamilyLikelihoodsUtils
    public static final String JOINT_POSTERIOR_TAG_NAME =           "JP"; //FamilyLikelihoodsUtils
    public final static String MIN_DP_FORMAT_KEY =                  "MIN_DP";
    public static final String MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY = "MQ0";
    public static final String HAPLOTYPE_CALLER_PHASING_GT_KEY =    "PGT";
    public static final String HAPLOTYPE_CALLER_PHASING_ID_KEY =    "PID";
    public static final String PHASE_SET_KEY =                      "PS";
    public static final String PHRED_SCALED_POSTERIORS_KEY =        "PP"; //FamilyLikelihoodsUtils / PosteriorLikelihoodsUtils
    public static final String REFERENCE_GENOTYPE_QUALITY =         "RGQ";
    public static final String STRAND_COUNT_BY_SAMPLE_KEY =         "SAC";
    public static final String STRAND_BIAS_BY_SAMPLE_KEY =          "SB";
    public static final String FEATURIZED_READ_SETS_KEY =           "FRS";
    public static final String FEATURIZED_READ_SETS_COUNTS_KEY =    "FRSCNT";
    public static final String HAPLOTYPE_EQUIVALENCE_COUNTS_KEY =   "HEC";
    public static final String HAPLOTYPE_COMPLEXITY_KEY =           "HAPCOMP";
    public static final String HAPLOTYPE_DOMINANCE_KEY =            "HAPDOM";
    public static final String HAPLOTYPES_BEFORE_FILTERING_KEY =    "ASSEMBLED_HAPS";
    public static final String HAPLOTYPES_FILTERED_KEY =            "FILTERED_HAPS";
    public final static String TRANSMISSION_PROBABILITY_KEY =       "TP"; //PhaseByTransmission
    public static final String FRAGMENT_ALLELE_DEPTHS =             "FAD";

    // flow annotations
    public final static String FLOW_INDEL_CLASSIFY =                "X_IC";
    public final static String FLOW_INDEL_LENGTH =                  "X_IL";
    public final static String FLOW_HMER_INDEL_LENGTH =             "X_HIL";
    public final static String FLOW_HMER_INDEL_NUC =                "X_HIN";
    public final static String FLOW_LEFT_MOTIF =                    "X_LM";
    public final static String FLOW_RIGHT_MOTIF =                   "X_RM";
    public final static String FLOW_GC_CONTENT =                    "X_GCC";
    public final static String FLOW_CYCLESKIP_STATUS =              "X_CSS";
    public final static String FLOW_VARIANT_TYPE =                  "VARIANT_TYPE";


    // M2-specific FORMAT keys
    public static final String ALLELE_FRACTION_KEY =                "AF";

    // Collapse keys
    public static final String EXT_COLLAPSED_KEY =                  "XC";

    //FILTERS
    /* Note that many filters used throughout GATK (most notably in VariantRecalibration) are dynamic,
       their names (or descriptions) depend on some threshold.  Those filters are not included here
     */
    public static final String CLUSTERED_EVENTS_FILTER_NAME =                 "clustered_events"; //M2
    public static final String GERMLINE_RISK_FILTER_NAME =                    "germline"; //M2
    public static final String LOW_QUAL_FILTER_NAME =                         "LowQual";
    public static final String ALIGNMENT_ARTIFACT_FILTER_NAME =               "alignment";
    public static final String PON_FILTER_NAME =                              "panel_of_normals"; //M2
    public static final String POLYMERASE_SLIPPAGE =                          "slippage"; //M2
    public static final String TUMOR_EVIDENCE_FILTER_NAME =                   "weak_evidence"; //M2
    public static final String MULTIALLELIC_FILTER_NAME =                     "multiallelic"; //M2
    public static final String STRAND_ARTIFACT_FILTER_NAME =                  "strand_bias"; // M2
    public static final String DUPLICATED_EVIDENCE_FILTER_NAME =              "duplicate";
    public final static String ARTIFACT_IN_NORMAL_FILTER_NAME =               "normal_artifact";
    public final static String MEDIAN_BASE_QUALITY_FILTER_NAME =              "base_qual";
    public final static String MEDIAN_MAPPING_QUALITY_FILTER_NAME =           "map_qual";
    public final static String MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME = "fragment";
    public final static String READ_POSITION_FILTER_NAME =                    "position";
    public final static String CONTAMINATION_FILTER_NAME =                    "contamination";
    public final static String READ_ORIENTATION_ARTIFACT_FILTER_NAME =        "orientation";
    public final static String BAD_HAPLOTYPE_FILTER_NAME =                    "haplotype";
    public final static String STRICT_STRAND_BIAS_FILTER_NAME =               "strict_strand";
    public final static String N_RATIO_FILTER_NAME =                           "n_ratio";
    public final static String ALLELE_FRACTION_FILTER_NAME =                   "low_allele_frac";
    public static final String POSSIBLE_NUMT_FILTER_NAME =                     "possible_numt";
    public static final String LOW_HET_FILTER_NAME =                           "mt_many_low_hets";
    public static final String FAIL =                                           "FAIL";
    public static final String SITE_LEVEL_FILTERS =                             "SITE";

    public static final String SNP = "SNP";
    public static final String INDEL = "INDEL";
    public static final String EXCESS_ALLELES = "EXCESS_ALLELES";
    public static final String NO_HQ_GENOTYPES = "NO_HQ_GENOTYPES";

    public static final String NAY_FROM_YNG = "NAY";
    public static final String VQSR_FAILURE_PREFIX = "low_VQSLOD_";
    public static final String VQSR_FAILURE_SNP = VQSR_FAILURE_PREFIX + SNP;
    public static final String VQSR_FAILURE_INDEL = VQSR_FAILURE_PREFIX + INDEL;

    public static final String CALIBRATION_SENSITIVITY_FAILURE_PREFIX = "high_CALIBRATION_SENSITIVITY_";
    public static final String CALIBRATION_SENSITIVITY_FAILURE_SNP = CALIBRATION_SENSITIVITY_FAILURE_PREFIX + SNP;
    public static final String CALIBRATION_SENSITIVITY_FAILURE_INDEL = CALIBRATION_SENSITIVITY_FAILURE_PREFIX + INDEL;



    public static final List<String> MUTECT_FILTER_NAMES = Arrays.asList(VCFConstants.PASSES_FILTERS_v4, POLYMERASE_SLIPPAGE,
            PON_FILTER_NAME, CLUSTERED_EVENTS_FILTER_NAME, TUMOR_EVIDENCE_FILTER_NAME, GERMLINE_RISK_FILTER_NAME,
            MULTIALLELIC_FILTER_NAME, STRAND_ARTIFACT_FILTER_NAME, ARTIFACT_IN_NORMAL_FILTER_NAME,
            MEDIAN_BASE_QUALITY_FILTER_NAME, MEDIAN_MAPPING_QUALITY_FILTER_NAME,
            MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME,
            READ_POSITION_FILTER_NAME, CONTAMINATION_FILTER_NAME, DUPLICATED_EVIDENCE_FILTER_NAME,
            READ_ORIENTATION_ARTIFACT_FILTER_NAME, BAD_HAPLOTYPE_FILTER_NAME,
            STRICT_STRAND_BIAS_FILTER_NAME, N_RATIO_FILTER_NAME, ALLELE_FRACTION_FILTER_NAME, POSSIBLE_NUMT_FILTER_NAME, FAIL);

    public static final List<String> MUTECT_AS_FILTER_NAMES = Arrays.asList(AS_FILTER_STATUS_KEY);

    // Symbolic alleles
    public final static String SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG = "ALT";
    public final static String NON_REF_SYMBOLIC_ALLELE_NAME = "NON_REF";
    public final static String SPANNING_DELETION_SYMBOLIC_ALLELE_NAME_DEPRECATED = "*:DEL";
    public final static Allele SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED = Allele.create("<" + SPANNING_DELETION_SYMBOLIC_ALLELE_NAME_DEPRECATED + ">", false); // represents any possible spanning deletion allele at this si
    public static final String ALLELE_SPECIFIC_ANNOTATION_PREFIX = "AS";
    public static final String PSEUDO_DEPTH_KEY = "DD";
    public static final String PSEUDO_FRACTION_KEY = "DF";

    public static boolean isSpanningDeletion(final Allele allele){
        return allele.equals(Allele.SPAN_DEL) || allele.equals(SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
    }
}
