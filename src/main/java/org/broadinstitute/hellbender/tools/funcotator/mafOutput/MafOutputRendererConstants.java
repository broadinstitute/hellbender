package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.funcotator.mafOutput.CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_FIELD;

/**
 * Class to hold all the constants required for the {@link MafOutputRenderer}.
 * Designed to be a simple container class with no methods.
 */
public class MafOutputRendererConstants {

    //==================================================================================================================
    // Static initializers:

    static {
        final Map<String, String> variantClassMap = new HashMap<>();

        variantClassMap.put(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),     "In_Frame_Del");
        variantClassMap.put(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),     "In_Frame_Ins");
        variantClassMap.put(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),  "Frame_Shift_Ins");
        variantClassMap.put(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),  "Frame_Shift_Del");
        variantClassMap.put(GencodeFuncotation.VariantClassification.MISSENSE.toString(),         "Missense_Mutation");
        variantClassMap.put(GencodeFuncotation.VariantClassification.NONSENSE.toString(),         "Nonsense_Mutation");
        variantClassMap.put(GencodeFuncotation.VariantClassification.SILENT.toString(),           "Silent");
        variantClassMap.put(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),      "Splice_Site");
        variantClassMap.put(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),  "Translation_Start_Site");
        variantClassMap.put(GencodeFuncotation.VariantClassification.NONSTOP.toString(),          "Nonstop_Mutation");
        variantClassMap.put(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),   "5'UTR");
        variantClassMap.put(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),  "3'UTR");
        variantClassMap.put(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(), "5'Flank");
        variantClassMap.put(GencodeFuncotation.VariantClassification.INTRON.toString(),           "Intron");
        variantClassMap.put(GencodeFuncotation.VariantClassification.LINCRNA.toString(),          "RNA");

        VariantClassificationMap = variantClassMap;
        VariantClassificationMapInverse = variantClassMap.entrySet().stream().collect(Collectors.toMap(Map.Entry::getValue, Map.Entry::getKey));
    }

    //==================================================================================================================
    // High-Level Constants:

    /**
     * Value to insert into unused annotation columns.
     */
    static final String UNUSED_STRING = "NA";

    /**
     * The string representing a comment in a MAF file.
     */
    static final String COMMENT_STRING = "#";

    /**
     * Delimiter for fields in the output MAF file.
     */
    static final String FIELD_DELIMITER = "\t";

    /**
     * Used for creating funcotations while rendering the MAF.
     */
    static final String MAF_COUNT_RENDERING_DATASOURCE_DUMMY_NAME = "MAF_COUNT_OUTPUT";
    static final String MAF_DBSNP_RENDERING_DATASOURCE_DUMMY_NAME = "MAF_DBSNP_OUTPUT";

    //==================================================================================================================
    // Specific Field Values:

    // Field Names:
    public static final String FieldName_Hugo_Symbol                            = "Hugo_Symbol";
    public static final String FieldName_Entrez_Gene_Id                         = "Entrez_Gene_Id";
    public static final String FieldName_Center                                 = "Center";
    public static final String FieldName_NCBI_Build                             = "NCBI_Build";
    public static final String FieldName_Chromosome                             = "Chromosome";
    public static final String FieldName_Start_Position                         = "Start_Position";
    public static final String FieldName_End_Position                           = "End_Position";
    public static final String FieldName_Strand                                 = "Strand";
    public static final String FieldName_Variant_Classification                 = "Variant_Classification";
    public static final String FieldName_Variant_Type                           = "Variant_Type";
    public static final String FieldName_Reference_Allele                       = "Reference_Allele";
    public static final String FieldName_Tumor_Seq_Allele1                      = "Tumor_Seq_Allele1";
    public static final String FieldName_Tumor_Seq_Allele2                      = "Tumor_Seq_Allele2";
    public static final String FieldName_dbSNP_RS                               = "dbSNP_RS";
    public static final String FieldName_dbSNP_Val_Status                       = "dbSNP_Val_Status";
    public static final String FieldName_Tumor_Sample_Barcode                   = "Tumor_Sample_Barcode";
    public static final String FieldName_Matched_Norm_Sample_Barcode            = "Matched_Norm_Sample_Barcode";
    public static final String FieldName_Match_Norm_Seq_Allele1                 = "Match_Norm_Seq_Allele1";
    public static final String FieldName_Match_Norm_Seq_Allele2                 = "Match_Norm_Seq_Allele2";
    public static final String FieldName_Tumor_Validation_Allele1               = "Tumor_Validation_Allele1";
    public static final String FieldName_Tumor_Validation_Allele2               = "Tumor_Validation_Allele2";
    public static final String FieldName_Match_Norm_Validation_Allele1          = "Match_Norm_Validation_Allele1";
    public static final String FieldName_Match_Norm_Validation_Allele2          = "Match_Norm_Validation_Allele2";
    public static final String FieldName_Verification_Status                    = "Verification_Status";
    public static final String FieldName_Validation_Status                      = "Validation_Status";
    public static final String FieldName_Mutation_Status                        = "Mutation_Status";
    public static final String FieldName_Sequencing_Phase                       = "Sequencing_Phase";
    public static final String FieldName_Sequence_Source                        = "Sequence_Source";
    public static final String FieldName_Validation_Method                      = "Validation_Method";
    public static final String FieldName_Score                                  = "Score";
    public static final String FieldName_BAM_File                               = "BAM_File";
    public static final String FieldName_Sequencer                              = "Sequencer";
    public static final String FieldName_Tumor_Sample_UUID                      = "Tumor_Sample_UUID";
    public static final String FieldName_Matched_Norm_Sample_UUID               = "Matched_Norm_Sample_UUID";
    public static final String FieldName_Genome_Change                          = "Genome_Change";
    public static final String FieldName_Annotation_Transcript                  = "Annotation_Transcript";
    public static final String FieldName_Transcript_Strand                      = "Transcript_Strand";
    public static final String FieldName_Transcript_Exon                        = "Transcript_Exon";
    public static final String FieldName_Transcript_Position                    = "Transcript_Position";
    public static final String FieldName_cDNA_Change                            = "cDNA_Change";
    public static final String FieldName_Codon_Change                           = "Codon_Change";
    public static final String FieldName_Protein_Change                         = "Protein_Change";
    public static final String FieldName_Other_Transcripts                      = "Other_Transcripts";
    public static final String FieldName_Refseq_mRNA_Id                         = "Refseq_mRNA_Id";
    public static final String FieldName_Refseq_prot_Id                         = "Refseq_prot_Id";
    public static final String FieldName_SwissProt_acc_Id                       = "SwissProt_acc_Id";
    public static final String FieldName_SwissProt_entry_Id                     = "SwissProt_entry_Id";
    public static final String FieldName_Description                            = "Description";
    public static final String FieldName_UniProt_AApos                          = "UniProt_AApos";
    public static final String FieldName_UniProt_Region                         = "UniProt_Region";
    public static final String FieldName_UniProt_Site                           = "UniProt_Site";
    public static final String FieldName_UniProt_Natural_Variations             = "UniProt_Natural_Variations";
    public static final String FieldName_UniProt_Experimental_Info              = "UniProt_Experimental_Info";
    public static final String FieldName_GO_Biological_Process                  = "GO_Biological_Process";
    public static final String FieldName_GO_Cellular_Component                  = "GO_Cellular_Component";
    public static final String FieldName_GO_Molecular_Function                  = "GO_Molecular_Function";
    public static final String FieldName_COSMIC_overlapping_mutations           = "COSMIC_overlapping_mutations";
    public static final String FieldName_COSMIC_fusion_genes                    = "COSMIC_fusion_genes";
    public static final String FieldName_COSMIC_tissue_types_affected           = "COSMIC_tissue_types_affected";
    public static final String FieldName_COSMIC_total_alterations_in_gene       = "COSMIC_total_alterations_in_gene";
    public static final String FieldName_Tumorscape_Amplification_Peaks         = "Tumorscape_Amplification_Peaks";
    public static final String FieldName_Tumorscape_Deletion_Peaks              = "Tumorscape_Deletion_Peaks";
    public static final String FieldName_TCGAscape_Amplification_Peaks          = "TCGAscape_Amplification_Peaks";
    public static final String FieldName_TCGAscape_Deletion_Peaks               = "TCGAscape_Deletion_Peaks";
    public static final String FieldName_DrugBank                               = "DrugBank";
    public static final String FieldName_ref_context                            = "ref_context";
    public static final String FieldName_gc_content                             = "gc_content";
    public static final String FieldName_CCLE_ONCOMAP_overlapping_mutations     = "CCLE_ONCOMAP_overlapping_mutations";
    public static final String FieldName_CCLE_ONCOMAP_total_mutations_in_gene   = "CCLE_ONCOMAP_total_mutations_in_gene";
    public static final String FieldName_CGC_Mutation_Type                      = "CGC_Mutation_Type";
    public static final String FieldName_CGC_Translocation_Partner              = "CGC_Translocation_Partner";
    public static final String FieldName_CGC_Tumor_Types_Somatic                = "CGC_Tumor_Types_Somatic";
    public static final String FieldName_CGC_Tumor_Types_Germline               = "CGC_Tumor_Types_Germline";
    public static final String FieldName_CGC_Other_Diseases                     = "CGC_Other_Diseases";
    public static final String FieldName_DNARepairGenes_Activity_linked_to_OMIM = "DNARepairGenes_Activity_linked_to_OMIM";
    public static final String FieldName_FamilialCancerDatabase_Syndromes       = "FamilialCancerDatabase_Syndromes";
    public static final String FieldName_MUTSIG_Published_Results               = "MUTSIG_Published_Results";
    public static final String FieldName_OREGANNO_ID                            = "OREGANNO_ID";
    public static final String FieldName_OREGANNO_Values                        = "OREGANNO_Values";
    public static final String FieldName_tumor_f                                = "tumor_f";
    public static final String FieldName_t_alt_count                            = "t_alt_count";
    public static final String FieldName_t_ref_count                            = "t_ref_count";
    public static final String FieldName_n_alt_count                            = "n_alt_count";
    public static final String FieldName_n_ref_count                            = "n_ref_count";


    // db SNP specific
    /**
     * The funcotation factory name that must be specified in the config for the custom MAF renderer to recognize it as a dbSNP annotation.
     */
    static final String DBSNP_DS_NAME = "dbSNP";

    /**
     * The funcotation factory result we should expect for the dbsnp validation flag.
     */
    static final String DBSNP_VLD_NAME = DBSNP_DS_NAME + "_VLD";

    // Field Values:
    static final String FieldValue_Strand                  = "+";
    static final String OTHER_TRANSCRIPT_DELIMITER         = "|";
    static final String FieldValue_Variant_Type_Insertion  = "INS";
    static final String FieldValue_Variant_Type_Deletion   = "DEL";
    static final String FieldValue_Gencode_Chromosome_Mito = "chrM";
    static final String FieldValue_Chromosome_Mito         = "MT";
    static final String EmptyAllele                        = "-";

    // Variant Classification Map:
    static final Map<String, String> VariantClassificationMap;
    static final Map<String, String> VariantClassificationMapInverse;

    // Output Field Name Map Defaults:
    //TODO these are hardcoded to gencode versions and should be updated to generalize to any version of gencode (see https://github.com/broadinstitute/gatk/issues/8482)
    static final List<String> OutputFieldNameMap_Hugo_Symbol                            = Arrays.asList(FieldName_Hugo_Symbol, "Gencode_19_hugoSymbol", "Gencode_27_hugoSymbol", "Gencode_28_hugoSymbol", "Gencode_34_hugoSymbol", "Gencode_43_hugoSymbol", "gene", "Gene");
    static final List<String> OutputFieldNameMap_Entrez_Gene_Id                         = Arrays.asList(FieldName_Entrez_Gene_Id, "HGNC_Entrez_Gene_ID", "HGNC_Entrez Gene ID", "HGNC_Entrez_Gene_ID(supplied_by_NCBI)", "HGNC_Entrez Gene ID(supplied by NCBI)", "entrez_id", "gene_id");
    static final List<String> OutputFieldNameMap_Center                                 = Arrays.asList(FieldName_Center, "center");
    static final List<String> OutputFieldNameMap_NCBI_Build                             = Arrays.asList(FieldName_NCBI_Build, "Gencode_19_ncbiBuild", "Gencode_27_ncbiBuild", "Gencode_28_ncbiBuild", "Gencode_34_ncbiBuild", "Gencode_43_ncbiBuild", "ncbi_build");
    static final List<String> OutputFieldNameMap_Chromosome                             = Arrays.asList(FieldName_Chromosome, "Gencode_19_chromosome", "Gencode_27_chromosome", "Gencode_28_chromosome", "Gencode_34_chromosome", "Gencode_43_chromosome", "chr", "contig", "chromosome", "chrom", "Chrom");
    static final List<String> OutputFieldNameMap_Start_Position                         = Arrays.asList(FieldName_Start_Position, "Start_position", "Gencode_19_start", "Gencode_27_start", "Gencode_28_start", "Gencode_34_start", "Gencode_43_start", "start", "Start", "start_pos", "pos");
    static final List<String> OutputFieldNameMap_End_Position                           = Arrays.asList(FieldName_End_Position, "End_position", "Gencode_19_end", "Gencode_27_end", "Gencode_28_end", "Gencode_34_end", "Gencode_43_end", "end", "End", "end_pos");
    static final List<String> OutputFieldNameMap_Strand                                 = Collections.singletonList(FieldName_Strand);
    static final List<String> OutputFieldNameMap_Variant_Classification                 = Arrays.asList(FieldName_Variant_Classification, "Gencode_19_variantClassification", "Gencode_27_variantClassification", "Gencode_28_variantClassification", "Gencode_34_variantClassification", "Gencode_43_variantClassification", "variant_classification");
    static final List<String> OutputFieldNameMap_Variant_Type                           = Arrays.asList(FieldName_Variant_Type, "Gencode_19_variantType", "Gencode_27_variantType", "Gencode_28_variantType", "Gencode_34_variantType", "Gencode_43_variantType", "variant_type");
    static final List<String> OutputFieldNameMap_Reference_Allele                       = Arrays.asList(FieldName_Reference_Allele, "Gencode_19_refAllele", "Gencode_27_refAllele", "Gencode_28_refAllele", "Gencode_34_refAllele", "Gencode_43_refAllele", "ref", "ref_allele", "reference_allele");
    static final List<String> OutputFieldNameMap_Tumor_Seq_Allele1                      = Arrays.asList(FieldName_Tumor_Seq_Allele1, "Gencode_19_tumorSeqAllele1", "Gencode_27_tumorSeqAllele1", "Gencode_28_tumorSeqAllele1", "Gencode_34_tumorSeqAllele1", "Gencode_43_tumorSeqAllele1", "ref", "ref_allele", "reference_allele");
    static final List<String> OutputFieldNameMap_Tumor_Seq_Allele2                      = Arrays.asList(FieldName_Tumor_Seq_Allele2, "Gencode_19_tumorSeqAllele2", "Gencode_27_tumorSeqAllele2", "Gencode_28_tumorSeqAllele2", "Gencode_34_tumorSeqAllele2", "Gencode_43_tumorSeqAllele2", "alt", "alt_allele", "alt2", "alt_allele2", "alternate_allele2", "observed_allele2", "alternate_allele", "observed_allele", "alt1", "alt_allele1", "alternate_allele1", "observed_allele1");
    static final List<String> OutputFieldNameMap_dbSNP_RS                               = Arrays.asList(FieldName_dbSNP_RS, "dbsnp_rs", "dbSNP_RSPOS");
    static final List<String> OutputFieldNameMap_dbSNP_Val_Status                       = Arrays.asList(FieldName_dbSNP_Val_Status, MAF_DBSNP_VAL_STATUS_FIELD, "dbsnp_val_status", DBSNP_VLD_NAME);
    static final List<String> OutputFieldNameMap_Tumor_Sample_Barcode                   = Arrays.asList(FieldName_Tumor_Sample_Barcode, "tumor_barcode", "tumor_id", "case_barcode", "case_id", "tumor_name");
    static final List<String> OutputFieldNameMap_Matched_Norm_Sample_Barcode            = Arrays.asList(FieldName_Matched_Norm_Sample_Barcode, "normal_barcode", "normal_id", "control_barcode", "control_id", "normal_name", "sample_name");
    static final List<String> OutputFieldNameMap_Match_Norm_Seq_Allele1                 = Arrays.asList(FieldName_Match_Norm_Seq_Allele1, "Match_Norm_Seq_Allele1");
    static final List<String> OutputFieldNameMap_Match_Norm_Seq_Allele2                 = Arrays.asList(FieldName_Match_Norm_Seq_Allele2, "Match_Norm_Seq_Allele2");
    static final List<String> OutputFieldNameMap_Tumor_Validation_Allele1               = Arrays.asList(FieldName_Tumor_Validation_Allele1, "Tumor_Validation_Allele1");
    static final List<String> OutputFieldNameMap_Tumor_Validation_Allele2               = Arrays.asList(FieldName_Tumor_Validation_Allele2, "Tumor_Validation_Allele2");
    static final List<String> OutputFieldNameMap_Match_Norm_Validation_Allele1          = Arrays.asList(FieldName_Match_Norm_Validation_Allele1, "Match_Norm_Validation_Allele1");
    static final List<String> OutputFieldNameMap_Match_Norm_Validation_Allele2          = Arrays.asList(FieldName_Match_Norm_Validation_Allele2, "Match_Norm_Validation_Allele2");
    static final List<String> OutputFieldNameMap_Verification_Status                    = Arrays.asList(FieldName_Verification_Status, "Verification_Status");
    static final List<String> OutputFieldNameMap_Validation_Status                      = Arrays.asList(FieldName_Validation_Status, "validation_status");
    static final List<String> OutputFieldNameMap_Mutation_Status                        = Arrays.asList(FieldName_Mutation_Status, "status");
    static final List<String> OutputFieldNameMap_Sequencing_Phase                       = Arrays.asList(FieldName_Sequencing_Phase, "phase");
    static final List<String> OutputFieldNameMap_Sequence_Source                        = Arrays.asList(FieldName_Sequence_Source, "source");
    static final List<String> OutputFieldNameMap_Validation_Method                      = Arrays.asList(FieldName_Validation_Method, "Validation_Method");
    static final List<String> OutputFieldNameMap_Score                                  = Arrays.asList(FieldName_Score, "Score");
    static final List<String> OutputFieldNameMap_BAM_File                               = Arrays.asList(FieldName_BAM_File, "BAM_file", "bam", "bam_file");
    static final List<String> OutputFieldNameMap_Sequencer                              = Arrays.asList(FieldName_Sequencer, "sequencer", "platform");
    static final List<String> OutputFieldNameMap_Tumor_Sample_UUID                      = Arrays.asList(FieldName_Tumor_Sample_UUID, "tumor_uuid", "case_uuid", "tumor_barcode", "tumor_id", "case_barcode", "case_id", "tumor_name", "Tumor_Sample_Barcode");
    static final List<String> OutputFieldNameMap_Matched_Norm_Sample_UUID               = Arrays.asList(FieldName_Matched_Norm_Sample_UUID, "normal_uuid", "control_uuid", "normal_barcode", "normal_id", "control_barcode", "control_id", "normal_name", "sample_name", "Matched_Norm_Sample_Barcode");
    static final List<String> OutputFieldNameMap_Genome_Change                          = Arrays.asList(FieldName_Genome_Change, "Gencode_19_genomeChange", "Gencode_27_genomeChange", "Gencode_28_genomeChange", "Gencode_34_genomeChange", "Gencode_43_genomeChange", "genome_change");
    static final List<String> OutputFieldNameMap_Annotation_Transcript                  = Arrays.asList(FieldName_Annotation_Transcript, "Gencode_19_annotationTranscript", "Gencode_27_annotationTranscript", "Gencode_28_annotationTranscript", "Gencode_34_annotationTranscript", "Gencode_43_annotationTranscript", "annotation_transcript", "transcript_id");
    static final List<String> OutputFieldNameMap_Transcript_Strand                      = Arrays.asList(FieldName_Transcript_Strand, "Gencode_19_transcriptStrand", "Gencode_27_transcriptStrand", "Gencode_28_transcriptStrand", "Gencode_34_transcriptStrand", "Gencode_43_transcriptStrand", "transcript_strand");
    static final List<String> OutputFieldNameMap_Transcript_Exon                        = Arrays.asList(FieldName_Transcript_Exon, "Gencode_19_transcriptExon", "Gencode_27_transcriptExon", "Gencode_28_transcriptExon", "Gencode_34_transcriptExon", "Gencode_43_transcriptExon", "transcript_exon");
    static final List<String> OutputFieldNameMap_Transcript_Position                    = Arrays.asList(FieldName_Transcript_Position, "Gencode_19_transcriptPos", "Gencode_27_transcriptPos", "Gencode_28_transcriptPos", "Gencode_34_transcriptPos",  "Gencode_43_transcriptPos", "transcript_position");
    static final List<String> OutputFieldNameMap_cDNA_Change                            = Arrays.asList(FieldName_cDNA_Change, "Gencode_19_cDnaChange", "Gencode_27_cDnaChange", "Gencode_28_cDnaChange", "Gencode_34_cDnaChange", "Gencode_43_cDnaChange", "transcript_change");
    static final List<String> OutputFieldNameMap_Codon_Change                           = Arrays.asList(FieldName_Codon_Change, "Gencode_19_codonChange", "Gencode_27_codonChange", "Gencode_28_codonChange", "Gencode_34_codonChange", "Gencode_43_codonChange", "codon_change");
    static final List<String> OutputFieldNameMap_Protein_Change                         = Arrays.asList(FieldName_Protein_Change, "Gencode_19_proteinChange", "Gencode_27_proteinChange", "Gencode_28_proteinChange", "Gencode_34_proteinChange", "Gencode_43_proteinChange", "protein_change");
    static final List<String> OutputFieldNameMap_Other_Transcripts                      = Arrays.asList(FieldName_Other_Transcripts, "Gencode_19_otherTranscripts", "Gencode_27_otherTranscripts", "Gencode_28_otherTranscripts", "Gencode_34_otherTranscripts", "Gencode_43_otherTranscripts", "other_transcripts");
    static final List<String> OutputFieldNameMap_Refseq_mRNA_Id                         = Arrays.asList(FieldName_Refseq_mRNA_Id, "Gencode_XRefSeq_mRNA_id", "gencode_xref_refseq_mRNA_id", "ENSEMBL_RefSeq_mRNA_accession", "RefSeq_mRNA_Id", "HGNC_RefSeq IDs");
    static final List<String> OutputFieldNameMap_Refseq_prot_Id                         = Arrays.asList(FieldName_Refseq_prot_Id, "Gencode_XRefSeq_prot_acc", "gencode_xref_refseq_prot_acc", "ENSEMBL_RefSeq_protein_accession", "RefSeq_prot_Id");
    static final List<String> OutputFieldNameMap_SwissProt_acc_Id                       = Arrays.asList(FieldName_SwissProt_acc_Id, "Simple_Uniprot_uniprot_accession", "uniprot_accession", "UniProt_uniprot_accession");
    static final List<String> OutputFieldNameMap_SwissProt_entry_Id                     = Arrays.asList(FieldName_SwissProt_entry_Id, "Simple_Uniprot_uniprot_entry_name", "uniprot_entry_name", "UniProt_uniprot_entry_name");
    static final List<String> OutputFieldNameMap_Description                            = Arrays.asList(FieldName_Description, "RefSeq_Description", "HGNC_Approved_Name", "HGNC_Approved Name");
    static final List<String> OutputFieldNameMap_UniProt_AApos                          = Arrays.asList(FieldName_UniProt_AApos, "UniProt_AAxform_aapos", "uniprot_AA_pos");
    static final List<String> OutputFieldNameMap_UniProt_Region                         = Arrays.asList(FieldName_UniProt_Region, "UniProt_AA_region");
    static final List<String> OutputFieldNameMap_UniProt_Site                           = Arrays.asList(FieldName_UniProt_Site, "UniProt_AA_site");
    static final List<String> OutputFieldNameMap_UniProt_Natural_Variations             = Arrays.asList(FieldName_UniProt_Natural_Variations, "UniProt_AA_natural_variation");
    static final List<String> OutputFieldNameMap_UniProt_Experimental_Info              = Arrays.asList(FieldName_UniProt_Experimental_Info, "UniProt_AA_experimental_info");
    static final List<String> OutputFieldNameMap_GO_Biological_Process                  = Arrays.asList(FieldName_GO_Biological_Process, "Simple_Uniprot_GO_Biological_Process", "UniProt_GO_Biological_Process");
    static final List<String> OutputFieldNameMap_GO_Cellular_Component                  = Arrays.asList(FieldName_GO_Cellular_Component, "Simple_Uniprot_GO_Cellular_Component", "UniProt_GO_Cellular_Component");
    static final List<String> OutputFieldNameMap_GO_Molecular_Function                  = Arrays.asList(FieldName_GO_Molecular_Function, "Simple_Uniprot_GO_Molecular_Function", "UniProt_GO_Molecular_Function");
    static final List<String> OutputFieldNameMap_COSMIC_overlapping_mutations           = Arrays.asList(FieldName_COSMIC_overlapping_mutations, "Cosmic_overlapping_mutations", "COSMIC_overlapping_mutations", "COSMIC_overlapping_mutation_AAs");
    static final List<String> OutputFieldNameMap_COSMIC_fusion_genes                    = Arrays.asList(FieldName_COSMIC_fusion_genes, "CosmicFusion_fusion_genes", "COSMIC_FusionGenes_fusion_genes");
    static final List<String> OutputFieldNameMap_COSMIC_tissue_types_affected           = Arrays.asList(FieldName_COSMIC_tissue_types_affected, "CosmicTissue_tissue_types_affected", "COSMIC_tissue_types_affected", "COSMIC_Tissue_tissue_types_affected");
    static final List<String> OutputFieldNameMap_COSMIC_total_alterations_in_gene       = Arrays.asList(FieldName_COSMIC_total_alterations_in_gene, "CosmicTissue_total_alterations_in_gene", "COSMIC_total_alterations_in_gene", "COSMIC_Tissue_total_alterations_in_gene");
    static final List<String> OutputFieldNameMap_Tumorscape_Amplification_Peaks         = Arrays.asList(FieldName_Tumorscape_Amplification_Peaks, "TUMORScape_Amplification_Peaks");
    static final List<String> OutputFieldNameMap_Tumorscape_Deletion_Peaks              = Arrays.asList(FieldName_Tumorscape_Deletion_Peaks, "TUMORScape_Deletion_Peaks");
    static final List<String> OutputFieldNameMap_TCGAscape_Amplification_Peaks          = Arrays.asList(FieldName_TCGAscape_Amplification_Peaks, "TCGAScape_Amplification_Peaks");
    static final List<String> OutputFieldNameMap_TCGAscape_Deletion_Peaks               = Arrays.asList(FieldName_TCGAscape_Deletion_Peaks, "TCGAScape_Deletion_Peaks");
    static final List<String> OutputFieldNameMap_DrugBank                               = Arrays.asList(FieldName_DrugBank, "Simple_Uniprot_DrugBank", "UniProt_DrugBank");
    static final List<String> OutputFieldNameMap_ref_context                            = Arrays.asList(FieldName_ref_context, "Gencode_19_referenceContext", "Gencode_27_referenceContext", "Gencode_28_referenceContext", "Gencode_34_referenceContext", "Gencode_43_referenceContext", "ref_context");
    static final List<String> OutputFieldNameMap_gc_content                             = Arrays.asList(FieldName_gc_content, "Gencode_19_gcContent", "Gencode_27_gcContent", "Gencode_28_gcContent", "Gencode_34_gcContent", "Gencode_43_gcContent", "gc_content");
    static final List<String> OutputFieldNameMap_CCLE_ONCOMAP_overlapping_mutations     = Arrays.asList(FieldName_CCLE_ONCOMAP_overlapping_mutations, "CCLE_By_GP_overlapping_mutations");
    static final List<String> OutputFieldNameMap_CCLE_ONCOMAP_total_mutations_in_gene   = Arrays.asList(FieldName_CCLE_ONCOMAP_total_mutations_in_gene, "CCLE_By_Gene_total_mutations_in_gene");
    static final List<String> OutputFieldNameMap_CGC_Mutation_Type                      = Arrays.asList(FieldName_CGC_Mutation_Type, "CGC_Mutation Type");
    static final List<String> OutputFieldNameMap_CGC_Translocation_Partner              = Arrays.asList(FieldName_CGC_Translocation_Partner, "CGC_Translocation Partner");
    static final List<String> OutputFieldNameMap_CGC_Tumor_Types_Somatic                = Arrays.asList(FieldName_CGC_Tumor_Types_Somatic, "CGC_Tumour Types  (Somatic Mutations)", "CGC_Tumour_Types__(Somatic_Mutations)");
    static final List<String> OutputFieldNameMap_CGC_Tumor_Types_Germline               = Arrays.asList(FieldName_CGC_Tumor_Types_Germline, "CGC_Tumour Types (Germline Mutations)", "CGC_Tumour_Types_(Germline_Mutations)");
    static final List<String> OutputFieldNameMap_CGC_Other_Diseases                     = Arrays.asList(FieldName_CGC_Other_Diseases, "CGC_Other Syndrome/Disease", "CGC_Other_Syndrome/Disease");
    static final List<String> OutputFieldNameMap_DNARepairGenes_Activity_linked_to_OMIM = Collections.singletonList(FieldName_DNARepairGenes_Activity_linked_to_OMIM);
    static final List<String> OutputFieldNameMap_FamilialCancerDatabase_Syndromes       = Arrays.asList(FieldName_FamilialCancerDatabase_Syndromes, "Familial_Cancer_Genes_Syndrome");
    static final List<String> OutputFieldNameMap_MUTSIG_Published_Results               = Arrays.asList(FieldName_MUTSIG_Published_Results, "MutSig Published Results_Published_Results");
    static final List<String> OutputFieldNameMap_OREGANNO_ID                            = Arrays.asList(FieldName_OREGANNO_ID, "Oreganno_ID", "ORegAnno_ID");
    static final List<String> OutputFieldNameMap_OREGANNO_Values                        = Arrays.asList(FieldName_OREGANNO_Values, "Oreganno_Values", "ORegAnno_Values");
    static final List<String> OutputFieldNameMap_tumor_f                                = Arrays.asList(FieldName_tumor_f, "sample_allelic_fraction");
    static final List<String> OutputFieldNameMap_t_alt_count                            = Collections.singletonList(FieldName_t_alt_count);
    static final List<String> OutputFieldNameMap_t_ref_count                            = Collections.singletonList(FieldName_t_ref_count);
    static final List<String> OutputFieldNameMap_n_alt_count                            = Collections.singletonList(FieldName_n_alt_count);
    static final List<String> OutputFieldNameMap_n_ref_count                            = Collections.singletonList(FieldName_n_ref_count);
}
