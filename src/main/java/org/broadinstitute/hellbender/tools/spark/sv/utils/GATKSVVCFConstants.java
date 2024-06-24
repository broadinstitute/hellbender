package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.variant.variantcontext.Allele;

public final class GATKSVVCFConstants {

    // todo: add these and the other standard SV info fields from the VCF spec to htsjdk VCFStandardHeaderLines
    // VCF standard keys reserved for sv
    public static final String SVTYPE = "SVTYPE";
    public static final String SVLEN = "SVLEN";
    public static final String IMPRECISE = "IMPRECISE";
    public static final String CIPOS = "CIPOS";
    public static final String CIEND = "CIEND";

    public static final String BREAKEND_STR = "BND";
    public static final String BND_MATEID_STR = "MATEID";

    // symbolic alt allele names
    public static final String SYMB_ALT_STRING_INV = "INV";
    public static final String SYMB_ALT_STRING_DEL = "DEL";
    public static final String SYMB_ALT_STRING_INS = "INS";
    public static final String SYMB_ALT_STRING_DUP = "DUP";
    public static final String SYMB_ALT_STRING_INVDUP = "DUP:INV";

    // symbolic alt alleles
    public static final Allele DEL_ALLELE = Allele.create("<DEL>", false);
    public static final Allele DUP_ALLELE = Allele.create("<DUP>", false);

    // GATK-SV specific header lines
    // TODO: 10/3/17 the following comment is a goal we are trying to achieve
    // applicable to all records all the time
    public static final String READ_PAIR_SUPPORT = "READ_PAIR_SUPPORT";
    public static final String SPLIT_READ_SUPPORT = "SPLIT_READ_SUPPORT";

    // applicable to all precise records all the time
    public static final String CONTIG_NAMES = "CTG_NAMES";
    public static final String TOTAL_MAPPINGS = "TOTAL_MAPPINGS";
    public static final String MAPPING_QUALITIES = "MAPPING_QUALITIES";
    public static final String HQ_MAPPINGS = "HQ_MAPPINGS";
    public static final String ALIGN_LENGTHS = "ALIGN_LENGTHS";
    public static final String MAX_ALIGN_LENGTH = "MAX_ALIGN_LENGTH";

    // applicable to all precise records when available
    public static final String SEQ_ALT_HAPLOTYPE = "SEQ_ALT_HAPLOTYPE";
    public static final String INSERTED_SEQUENCE = "INSSEQ";
    public static final String INSERTED_SEQUENCE_LENGTH = "INSLEN";
    public static final String INSERTED_SEQUENCE_MAPPINGS = "INSSEQ_MAP";
    public static final String HOMOLOGY = "HOMSEQ";
    public static final String HOMOLOGY_LENGTH = "HOMLEN";
    public static final String LINK = "LINK";
    public static final String EXTERNAL_CNV_CALLS = "EXTERNAL_CNV_CALLS";
    public static final String CTG_GOOD_NONCANONICAL_MAPPING = "CTG_GOOD_NONCANONICAL_MAPPING";

    // type specific: tandem duplication
    public static final String DUP_REPEAT_UNIT_REF_SPAN = "DUP_REPEAT_UNIT_REF_SPAN";
    public static final String DUP_SEQ_CIGARS = "DUP_SEQ_CIGARS";
    public static final String DUPLICATION_NUMBERS = "DUP_NUM";
    public static final String DUP_ANNOTATIONS_IMPRECISE = "DUP_ANNOTATIONS_IMPRECISE";
    public static final String DUP_IMPRECISE_AFFECTED_RANGE = "DUP_IMPRECISE_AFFECTED_RANGE";

    public static final String DUP_TAN_CONTRACTION_STRING = "CONTRACTION";
    public static final String DUP_TAN_EXPANSION_STRING = "EXPANSION";

    // type specific: inverted duplication
    public static final String DUP_ORIENTATIONS = "DUP_ORIENTATIONS";

    // type specific: inversion
    public static final String INV33 = "INV33";
    public static final String INV55 = "INV55";

    // type specific: CPX
    public static final String CPX_SV_SYB_ALT_ALLELE_STR = "CPX";
    public static final String CPX_EVENT_ALT_ARRANGEMENTS = "ALT_ARRANGEMENT";
    public static final String CPX_SV_REF_SEGMENTS = "SEGMENTS";
    public static final String CPX_EVENT_KEY = "CPX_EVENT";
    public static final String CPX_INTERVALS = "CPX_INTERVALS";
    public static final String CPX_TYPE = "CPX_TYPE";

    public enum ComplexVariantSubtype {
        delINV,
        INVdel,
        dupINV,
        INVdup,
        delINVdel,
        dupINVdup,
        delINVdup,
        dupINVdel,
        piDUP_FR,
        piDUP_RF,
        dDUP,
        dDUP_iDEL,
        INS_iDEL,
        CTX_PP_QQ,
        CTX_PQ_QP,
        CTX_INV
    }

    // not defined in output vcf header but used in internal id that is currently output in the ID column
    public static final String INTERVAL_VARIANT_ID_FIELD_SEPARATOR = "_";
    public static final String DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING = "DEL-DUPLICATION-TANDEM-CONTRACTION";
    public static final String DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING = "INS-DUPLICATION-TANDEM-EXPANSION";
    public static final String DUP_INV_INTERNAL_ID_START_STRING = "INS-DUPLICATION-INVERTED-EXPANSION";

    // for breakpoint segmentation
    public static final String ALGORITHMS_ATTRIBUTE = "ALGORITHMS";
    public static final String STRANDS_ATTRIBUTE = "STRANDS";
    public static final String DEPTH_ALGORITHM = "depth";
    public static final String CONTIG2_ATTRIBUTE = "CHR2";
    public static final String END2_ATTRIBUTE = "END2";
    public static final String BND_DELETION_STRANDS = "+-";
    public static final String BND_DUPLICATION_STRANDS = "-+";

    // format block
    public static final String COPY_NUMBER_FORMAT = "CN";
    public static final String EXPECTED_COPY_NUMBER_FORMAT = "ECN";
    public static final String COPY_NUMBER_QUALITY_FORMAT = "CNQ";

    // filter block
    public static final String ASSEMBLY_BASED_VARIANT_MQ_FILTER_KEY = "LOW_MQ";
    public static final String ASSEMBLY_BASED_VARIANT_ALN_LENGTH_FILTER_KEY = "SHORT_ALN";
    public static final String LOW_QS_SCORE_FILTER_KEY = "LOW_QS";
    public static final String FREQUENCY_FILTER_KEY = "FREQ";

    // Clustering
    public static final String CLUSTER_MEMBER_IDS_KEY = "MEMBERS";

    // Concordance
    public static final String GENOTYPE_CONCORDANCE_INFO = "GENOTYPE_CONCORDANCE";
    public static final String NON_REF_GENOTYPE_CONCORDANCE_INFO = "NON_REF_GENOTYPE_CONCORDANCE";

    public static final String HET_PPV_INFO = "HET_PPV";
    public static final String HET_SENSITIVITY_INFO = "HET_SENSITIVITY";
    public static final String HET_SPECIFICITY_INFO = "HET_SPECIFICITY";

    public static final String HOMVAR_PPV_INFO = "HOMVAR_PPV";
    public static final String HOMVAR_SENSITIVITY_INFO = "HOMVAR_SENSITIVITY";
    public static final String HOMVAR_SPECIFICITY_INFO = "HOMVAR_SPECIFICITY";

    public static final String VAR_PPV_INFO = "VAR_PPV";
    public static final String VAR_SENSITIVITY_INFO = "VAR_SENSITIVITY";
    public static final String VAR_SPECIFICITY_INFO = "VAR_SPECIFICITY";

    public static final String TRUTH_CN_EQUAL_FORMAT = "TRUTH_CN_EQUAL";
    public static final String COPY_NUMBER_CONCORDANCE_INFO = "CNV_CONCORDANCE";

    public static final String TRUTH_VARIANT_ID_INFO = "TRUTH_VID";

    public static final String TRUTH_ALLELE_COUNT_INFO = "TRUTH_AC";
    public static final String TRUTH_ALLELE_NUMBER_INFO = "TRUTH_AN";
    public static final String TRUTH_ALLELE_FREQUENCY_INFO = "TRUTH_AF";

    // functional annotations
    public static final String LOF = "PREDICTED_LOF";
    public static final String INT_EXON_DUP = "PREDICTED_INTRAGENIC_EXON_DUP";
    public static final String COPY_GAIN = "PREDICTED_COPY_GAIN";
    public static final String DUP_PARTIAL = "PREDICTED_DUP_PARTIAL";
    public static final String PARTIAL_EXON_DUP = "PREDICTED_PARTIAL_EXON_DUP";
    public static final String INTRONIC = "PREDICTED_INTRONIC";
    public static final String INV_SPAN = "PREDICTED_INV_SPAN";
    public static final String UTR = "PREDICTED_UTR";
    public static final String MSV_EXON_OVERLAP = "PREDICTED_MSV_EXON_OVERLAP";
    public static final String PROMOTER = "PREDICTED_PROMOTER";
    public static final String BREAKEND_EXON = "PREDICTED_BREAKEND_EXONIC";
    public static final String INTERGENIC = "PREDICTED_INTERGENIC";
    public static final String NONCODING_SPAN = "PREDICTED_NONCODING_SPAN";
    public static final String NONCODING_BREAKPOINT = "PREDICTED_NONCODING_BREAKPOINT";
    public static final String NEAREST_TSS = "PREDICTED_NEAREST_TSS";
    public static final String TSS_DUP = "PREDICTED_TSS_DUP";
    public static final String PARTIAL_DISPERSED_DUP = "PREDICTED_PARTIAL_DISPERSED_DUP";

    // SVTYPE classes
    public enum StructuralVariantAnnotationType {
        DEL,
        DUP,
        INS,
        INV,
        CPX,
        BND,
        CTX,
        CNV
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
