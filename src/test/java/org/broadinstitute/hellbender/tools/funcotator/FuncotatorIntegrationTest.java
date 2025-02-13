package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.CustomMafFuncotationCreator;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRendererConstants;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.validation.constraints.AssertFalse;
import javax.validation.constraints.AssertTrue;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * An integration test for the {@link Funcotator} tool.
 * Created by jonn on 8/29/17.
 */
public class FuncotatorIntegrationTest extends CommandLineProgramTest {

    // Temp directory in which to place output files.
    private static final File tmpOutDir;

    // Whether to do debug output (i.e. leave output around).
    // This should always be false when checked in.
    // These tests would take ~30 minutes to complete each.
    private static final boolean enableFullScaleValidationTest = false;
    private static final String  LARGE_DATASOURCES_FOLDER      = "funcotator_dataSources_latest";
    private static final String  GERMLINE_DATASOURCES_FOLDER   = "funcotator_dataSources_germline_latest";

    private static final String XSV_CLINVAR_MULTIHIT_TEST_VCF = toolsTestDir + "funcotator" + File.separator + "clinvar_hg19_multihit_test.vcf";
    private static final String FILTER_TEST_VCF               = toolsTestDir + "funcotator" + File.separator + "FILTER_test.vcf";
    private static final String VCF_FIELD_ORDER_SWAP_TEST_VCF = toolsTestDir + "funcotator" + File.separator + "vcfBugRepro.vcf";
    private static final String DS_XSV_CLINVAR_TESTS          = largeFileTestDir + "funcotator" + File.separator + "small_ds_clinvar_hg19" + File.separator;
    private static final String DS_FILTER_PARSE_TESTS         = largeFileTestDir + "funcotator" + File.separator + "small_ds_FILTER_test" + File.separator;
    public static final String VCF_FIELD_ORDER_TEST_DATA_SOURCES = largeFileTestDir + "funcotator" + File.separator + "vcfFuncotationOrderingBugRepro" + File.separator;

    private static final String NOT_M2_TEST_HG19 = toolsTestDir + "funcotator/NotM2_test_custom_maf_fields.vcf";
    private static final String M2_TEST_HG19 = toolsTestDir + "funcotator/M2_test_custom_maf_fields.vcf";
    private static final String NOT_M2_TEST_HG19_TUMOR_ONLY = toolsTestDir + "funcotator/NotM2_test_custom_maf_fields_tumor_only.vcf";
    private static final String M2_TEST_HG19_TUMOR_ONLY = toolsTestDir + "funcotator/M2_test_custom_maf_fields_tumor_only.vcf";
    private static final String THREE_SAMPLE_SOMATIC = toolsTestDir + "funcotator/Three_sample_somatic.vcf";
    private static final String PIK3CA_VCF_HG19          = toolsTestDir + "funcotator" + File.separator + "0816201804HC0_R01C01.pik3ca.vcf";
    private static final String PIK3CA_VCF_HG38          = toolsTestDir + "funcotator" + File.separator + "hg38_trio.pik3ca.vcf";
    private static final String PIK3CA_VCF_HG19_SNPS     = toolsTestDir + "funcotator" + File.separator + "PIK3CA_SNPS_3.vcf";
    private static final String PIK3CA_VCF_HG19_INDELS   = toolsTestDir + "funcotator" + File.separator + "PIK3CA_INDELS_3.vcf";
    private static final String MUC16_VCF_HG19           = toolsTestDir + "funcotator" + File.separator + "MUC16_MNP.vcf";
    private static final String PIK3CA_VCF_HG19_ALTS     = toolsTestDir + "funcotator" + File.separator + "PIK3CA_3_miss_clinvar_alt_only.vcf";
    private static final String SPANNING_DEL_VCF         = toolsTestDir + "funcotator" + File.separator + "spanning_del.vcf";

    // TODO: Get rid of this variable and use the general data sources path (issue #5350 - https://github.com/broadinstitute/gatk/issues/5350):
    private static final String DS_PIK3CA_DIR            = largeFileTestDir + "funcotator" + File.separator + "small_ds_pik3ca" + File.separator;
    private static final String DS_PIK3CA_DIR_v43_MANE   = largeFileTestDir + "funcotator" + File.separator + "ds_pik3ca_v43_WITH_MANE_SELECT_PLUS_CLINICAL_ANNOTATIONS_FOR_TEST" + File.separator;

    private static final String  MEDIUM_DATASOURCES_DIR     = largeFileTestDir + "funcotator" + File.separator + "funcotator_dataSources" + File.separator;

    private static final String MAF_TEST_CONFIG          = toolsTestDir + "funcotator" + File.separator + "maf.config";
    private static final String XSV_CLINVAR_COL_TEST_VCF = toolsTestDir + "funcotator" + File.separator + "clinvar_hg19_column_test.vcf";
    private static final String DS_XSV_CLINVAR_COL_TEST  = largeFileTestDir + "funcotator" + File.separator + "small_ds_clinvar_hg19" + File.separator;
    private static final String EMPTY_VCF  = publicTestDir + File.separator + "empty.vcf";
    private static final String PIK3CA_DBSNP_DS          = toolsTestDir + "funcotator" + File.separator + "small_pik3ca_dbsnp_ds";
    private static final String MAF_DBSNP_TEST           = toolsTestDir + "funcotator" + File.separator + "maf_dbsnp_test_input.vcf";

    // E coli data sources:
    private static final String DS_ECOLI_DIR             = largeFileTestDir + "funcotator" + File.separator + "ecoli_ds" + File.separator;
    private static final String E_COLI_EXPECTED_OUT      = largeFileTestDir + "funcotator" + File.separator + "e_coli.MG1655.expected_output.vcf";

    // Non-locatable funcotation file:
    private static final String NON_LOCATABLE_FUNCOTATED_INPUT_VCF = toolsTestDir + "funcotator" + File.separator + "non_locatable_proof_input.vcf";

    // Symbollic allele funcotation file:
    private static final String SYMBOLLIC_ALLELE_FUNCOTATED_INPUT_VCF = toolsTestDir + "funcotator" + File.separator + "symbollic_allele_proof_input.vcf";

    private static final List<String> VCF_FIELDS_GENCODE_19_DS = Arrays.asList("Gencode_19_hugoSymbol","Gencode_19_ncbiBuild","Gencode_19_chromosome","Gencode_19_start","Gencode_19_end","Gencode_19_variantClassification","Gencode_19_variantType","Gencode_19_refAllele","Gencode_19_tumorSeqAllele1","Gencode_19_tumorSeqAllele2","Gencode_19_genomeChange","Gencode_19_annotationTranscript","Gencode_19_transcriptStrand","Gencode_19_transcriptExon","Gencode_19_transcriptPos","Gencode_19_cDnaChange","Gencode_19_codonChange","Gencode_19_proteinChange","Gencode_19_gcContent","Gencode_19_referenceContext","Gencode_19_otherTranscripts");//,"Achilles_Top_Genes","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_HGMD_ID","ClinVar_SYM","ClinVar_TYPE","ClinVar_ASSEMBLY","ClinVar_rs","Cosmic_overlapping_mutations","CosmicFusion_fusion_genes","CosmicFusion_fusion_id","CosmicTissue_total_alterations_in_gene","CosmicTissue_tissue_types_affected","DNARepairGenes_Activity_linked_to_OMIM","DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView","DNARepairGenes_Accession_number_linked_to_NCBI_Entrez","Familial_Cancer_Genes_Syndrome","Familial_Cancer_Genes_Synonym","Familial_Cancer_Genes_Reference","Gencode_XHGNC_hgnc_id","Gencode_XRefSeq_mRNA_id","Gencode_XRefSeq_prot_acc","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","Oreganno_Build","Oreganno_ID","Oreganno_Values","Simple_Uniprot_uniprot_entry_name","Simple_Uniprot_DrugBank","Simple_Uniprot_alt_uniprot_accessions","Simple_Uniprot_uniprot_accession","Simple_Uniprot_GO_Biological_Process","Simple_Uniprot_GO_Cellular_Component","Simple_Uniprot_GO_Molecular_Function","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID");
    private static final List<String> VCF_FIELDS_GENCODE_28_DS = Arrays.asList("Gencode_28_hugoSymbol","Gencode_28_ncbiBuild","Gencode_28_chromosome","Gencode_28_start","Gencode_28_end","Gencode_28_variantClassification","Gencode_28_variantType","Gencode_28_refAllele","Gencode_28_tumorSeqAllele1","Gencode_28_tumorSeqAllele2","Gencode_28_genomeChange","Gencode_28_annotationTranscript","Gencode_28_transcriptStrand","Gencode_28_transcriptExon","Gencode_28_transcriptPos","Gencode_28_cDnaChange","Gencode_28_codonChange","Gencode_28_proteinChange","Gencode_28_gcContent","Gencode_28_referenceContext","Gencode_28_otherTranscripts");//,"Achilles_Top_Genes","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_HGMD_ID","ClinVar_SYM","ClinVar_TYPE","ClinVar_ASSEMBLY","ClinVar_rs","Cosmic_overlapping_mutations","CosmicFusion_fusion_genes","CosmicFusion_fusion_id","CosmicTissue_total_alterations_in_gene","CosmicTissue_tissue_types_affected","DNARepairGenes_Activity_linked_to_OMIM","DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView","DNARepairGenes_Accession_number_linked_to_NCBI_Entrez","Familial_Cancer_Genes_Syndrome","Familial_Cancer_Genes_Synonym","Familial_Cancer_Genes_Reference","Gencode_XHGNC_hgnc_id","Gencode_XRefSeq_mRNA_id","Gencode_XRefSeq_prot_acc","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","Oreganno_Build","Oreganno_ID","Oreganno_Values","Simple_Uniprot_uniprot_entry_name","Simple_Uniprot_DrugBank","Simple_Uniprot_alt_uniprot_accessions","Simple_Uniprot_uniprot_accession","Simple_Uniprot_GO_Biological_Process","Simple_Uniprot_GO_Cellular_Component","Simple_Uniprot_GO_Molecular_Function","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID");
    private static final List<String> MAF_FIELDS_GENCODE_DS = Arrays.asList(MafOutputRendererConstants.FieldName_Hugo_Symbol, MafOutputRendererConstants.FieldName_NCBI_Build, MafOutputRendererConstants.FieldName_Chromosome,
            MafOutputRendererConstants.FieldName_Start_Position, MafOutputRendererConstants.FieldName_End_Position, MafOutputRendererConstants.FieldName_Variant_Classification,  MafOutputRendererConstants.FieldName_Variant_Type,
            MafOutputRendererConstants.FieldName_Reference_Allele, MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, MafOutputRendererConstants.FieldName_Genome_Change,
            MafOutputRendererConstants.FieldName_Annotation_Transcript, MafOutputRendererConstants.FieldName_Transcript_Strand, MafOutputRendererConstants.FieldName_Transcript_Exon, MafOutputRendererConstants.FieldName_Transcript_Position,
            MafOutputRendererConstants.FieldName_cDNA_Change, MafOutputRendererConstants.FieldName_Codon_Change, MafOutputRendererConstants.FieldName_Protein_Change, MafOutputRendererConstants.FieldName_gc_content,
            MafOutputRendererConstants.FieldName_ref_context, MafOutputRendererConstants.FieldName_Other_Transcripts);

    private static String hg38Chr3Ref;
    private static String hg19Chr3Ref;
    private static String b37Chr2Ref;
    private static String hg19Chr19Ref;

    private static String eColiRef;

    static {
        // This is intentionally set here so that output can be examined in the case of running full scale tests.
        if (!enableFullScaleValidationTest ) {
            tmpOutDir = createTempDir("funcotatorTmpFolder");
        } else {
            tmpOutDir = new File("funcotatorTmpFolder" + File.separator);
            if (!tmpOutDir.mkdirs() && !tmpOutDir.exists()) {
                throw new GATKException("Error making output folder for test: " + FuncotatorIntegrationTest.class.getName());
            }
        }

        hg38Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg38Chr3Ref();
        hg19Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref();
        b37Chr2Ref = FuncotatorReferenceTestUtils.retrieveB37Chr2Ref();
        hg19Chr19Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref();
        eColiRef = FuncotatorReferenceTestUtils.retrieveEcoliReference();
    }

    //==================================================================================================================
    // Disabled tests to regenerate expected outputs for integration tests:
    @Test(dataProvider = "provideForNonTrivialLargeDataValidationTest",
               enabled = false)
    public void regenerateExpectedOutputsForNonTrivialLargeDataValidationTest(
                                                  final String inputVcfName,
                                                  final String referencePath,
                                                  final String referenceVersion,
                                                  final String dataSourcesPath,
                                                  final String expectedOutputPath,
                                                  final List<String> excludedFields,
                                                  final boolean reannotate) {

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values() ) {
            if (outputFormatType.equals(FuncotatorArgumentDefinitions.OutputFormatType.SEG)) {
                continue;
            }
            final File outputFile = createTempFile(tmpOutDir + File.separator + inputVcfName + ".funcotator", "." + outputFormatType.toString().toLowerCase());

            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    inputVcfName,
                    outputFile,
                    referencePath,
                    dataSourcesPath,
                    referenceVersion,
                    outputFormatType,
                    true);

            // Run the tool with our args:
            runCommandLine(arguments);

            // ---------------------------------

            final String typeCorrectedExpectedOutputPath = FilenameUtils.removeExtension(expectedOutputPath) + "." + outputFormatType.toString().toLowerCase();

            // Copy the output file to our expected output area:
            try {
                Files.copy(outputFile.toPath(), IOUtils.getPath(typeCorrectedExpectedOutputPath), StandardCopyOption.REPLACE_EXISTING);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Unable to copy over generated data to expected output path!", ex);
            }
        }
    }

    //==================================================================================================================
    // Helper methods to create output files and maybe leave them around to debug the test output.

    private static File getOutputFile(final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType) {
        return getOutputFile("funcotator_tmp_out", outputFormatType.toString().toLowerCase());
    }

    private static File getOutputFile(final String outfileBaseName,
                                      final String outFileExtension) {
        final File outputFile;
        if (!enableFullScaleValidationTest ) {
            outputFile = createTempFile(tmpOutDir + File.separator + outfileBaseName, "." + outFileExtension);
        } else {
            outputFile = new File(tmpOutDir, outfileBaseName + "." + outFileExtension);
        }
        return outputFile;
    }

    /**
     * @return A standard list of manual annotations and their values that can be added to tests in order to test the manual annotations piping.
     */
    private static List<String> getManualAnnotations() {

        final List<String> annotationList = new ArrayList<>();

        annotationList.add("dbSNP_RS:0");
        annotationList.add("dbSNP_Val_Status:No_Value");
        annotationList.add("Center:broad.mit.edu");
        annotationList.add("source:WES");
        annotationList.add("normal_barcode:normal_sample");
        annotationList.add("tumor_barcode:tumor_sample");
        annotationList.add("NCBI_Build:37");
        annotationList.add("Strand:+");
        annotationList.add("status:Somatic");
        annotationList.add("phase:Phase_I");
        annotationList.add("sequencer:Illumina");
        annotationList.add("Tumor_Validation_Allele1:");
        annotationList.add("Tumor_Validation_Allele2:");
        annotationList.add("Match_Norm_Validation_Allele1:");
        annotationList.add("Match_Norm_Validation_Allele2:");
        annotationList.add("Verification_Status:");
        annotationList.add("Validation_Status:");
        annotationList.add("Validation_Method:");
        annotationList.add("Score:");
        annotationList.add("BAM_file:");
        annotationList.add("Match_Norm_Seq_Allele1:");
        annotationList.add("Match_Norm_Seq_Allele2:");

        return annotationList;
    }


    private static List<String> getAnnotationOverrides() {
        final List<String> annotationList = new ArrayList<>();
        annotationList.add("Oreganno_Build:BUILDED_GOOD_REAL_BIG");
        return annotationList;
    }

    private static void addManualAnnotationsToArguments(final ArgumentsBuilder arguments) {

        // ================================================================================
        // Annotation Defaults:
        for ( final String manualAnnotation : getManualAnnotations() ) {
            arguments.add(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, manualAnnotation);
        }
        // ================================================================================
        // Annotation Overrides:
        for ( final String annotationOverride : getAnnotationOverrides() ) {
            arguments.add(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME, annotationOverride);
        }
    }

    private ArgumentsBuilder createBaselineArgumentsForFuncotator(final String variantFileName,
                                                                  final File outputFile,
                                                                  final String referenceFileName,
                                                                  final String dataSourcesPath,
                                                                  final String refVer,
                                                                  final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                                                  final boolean shouldValidateSeqDicts) {
        return createBaselineArgumentsForFuncotator(variantFileName,
            outputFile,
            referenceFileName,
            dataSourcesPath,
            refVer,
            outputFormatType,
            shouldValidateSeqDicts,
            Collections.emptyList(),
            false);
    }

    private ArgumentsBuilder createBaselineArgumentsForFuncotator(final String variantFileName,
                                                                  final File outputFile,
                                                                  final String referenceFileName,
                                                                  final String dataSourcesPath,
                                                                  final String refVer,
                                                                  final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                                                  final boolean shouldValidateSeqDicts,
                                                                  final List<String> excludedFields,
                                                                  final boolean reannotate) {

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(variantFileName));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(referenceFileName));
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, dataSourcesPath);
        arguments.add(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, refVer);
        arguments.add(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.add("verbosity", "INFO");
        excludedFields.forEach(ef -> arguments.add(FuncotatorArgumentDefinitions.EXCLUSION_FIELDS_LONG_NAME, ef));

        if ( !shouldValidateSeqDicts ) {
            // Disable the sequence dictionary check for the tests:
            arguments.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);
        }

        if ( reannotate ) {
            arguments.add(FuncotatorArgumentDefinitions.REANNOTATE_VCF_LONG_NAME, true);
        }

        return arguments;
    }

    //==================================================================================================================

    @DataProvider
    private Object[][] provideForNonTrivialLargeDataValidationTest() {

        return new Object[][] {
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1_EXPECTED_OUTPUT,
                        Collections.emptyList(),
                        false
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_2,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_2_EXPECTED_OUTPUT,
                        Collections.emptyList(),
                        false
                },
                {
                    //This tests https://github.com/broadinstitute/gatk/issues/6173
                        FuncotatorTestConstants.SINGLE_LINE,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.SINGLE_LINE_EXPECTED,
                        Collections.emptyList(),
                        false
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG38,
                        hg38Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_EXPECTED_OUTPUT,
                        Collections.emptyList(),
                        false
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_LARGE_DATA_SET,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_LARGE_DATA_SET_EXPECTED_OUTPUT,
                        Collections.emptyList(),
                        true
                },
        };
    }

    @DataProvider
    Iterator<Object[]> provideForIntegrationTest() {

        final ArrayList<Object[]> testCases = new ArrayList<>();

        // TODO: Fix this set of tests!  THEY DON'T ALL PASS!  (Issue: https://github.com/broadinstitute/gatk/issues/4344)
        // VCF SNP / MNP / INDEL tests
//        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values() ) {
//            for ( final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values() ) {
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//            }
//        }

        // Basic Test Cases:
        for (final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            for (final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values()) {

                // Funcotator, the command line tool, does not support SEG outputs.  Users must use FuncotateSegments.
                if (outputFormatType.equals(FuncotatorArgumentDefinitions.OutputFormatType.SEG)) {
                    continue;
                }

                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                                hg19Chr3Ref,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3,
                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                                hg19Chr19Ref,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19,
                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
            }
        }

        return testCases.iterator();
    }

    @DataProvider
    public Object[][] provideForLargeDataValidationTest() {
        return new Object[][]{
                {
                        "0816201804HC0_R01C01.vcf",
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GERMLINE_DATASOURCES_FOLDER,
                        true
                },
                {
                        "0816201804HC0_R01C01.vcf",
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_LOCAL_CLOUD_FOLDER,
                        false
                },
                {
                        "0816201804HC0_R01C01.vcf",
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_REMOTE_CLOUD_FOLDER,
                        false
                },
                {
                        "hg38_test_variants.vcf",
                        hg38Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        LARGE_DATASOURCES_FOLDER,
                        true
                },
                {
                        "hg38_trio.vcf",
                        hg38Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        LARGE_DATASOURCES_FOLDER,
                        true
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        false
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_2,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        false
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG38,
                        hg38Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        false
                },
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_LARGE_DATA_SET,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        false
                },
        };
    }

    //==================================================================================================================

    // This test is to make sure we don't create a bunch of temp files anywhere.
    // It will force anyone who changes the outputToTmpDir flag to make it true when they check in this test file.
    //
    // DO NOT ADD THIS TO ANY TEST GROUPS!
    @Test
    public void metaTestEnsureTempDirs() {
        Assert.assertFalse(enableFullScaleValidationTest);
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void testFuncotatorWithoutValidatingResults(final String dataSourcesPath,
                                                       final String refVer,
                                                       final String referenceFileName,
                                                       final String variantFileName,
                                                       final String transcriptName,
                                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                                       final String expectedOutputFilePath) throws IOException {

        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(variantFileName, outputFile, referenceFileName, dataSourcesPath, refVer, outputFormatType, false);

        runCommandLine(arguments);
    }

    private void validateFuncotationsOnVcf(final Iterable<VariantContext> vcfIterable, final String[] funcotationFieldNames) {
        for (final VariantContext vc : vcfIterable ) {
            final String funcotation = vc.getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, "");

            Assert.assertNotEquals(funcotation, "");

            final String rawFuncotations = funcotation.substring(1,funcotation.length()-1);

            Assert.assertEquals(StringUtils.countMatches(rawFuncotations, VcfOutputRenderer.FIELD_DELIMITER), funcotationFieldNames.length - 1);

            // This is here to make sure we can create the FuncotationMap object without exploding.
            // It serves as a secondary check.
            final FuncotationMap funkyMap = FuncotationMap.createAsAllTableFuncotationsFromVcf(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, funcotationFieldNames,
                    funcotation, vc.getAlternateAllele(0), "VCF");
        }
    }

    // This test will take no less than 2 hours to run, hence it is disabled by default.
    @Test(enabled = enableFullScaleValidationTest,
          groups = {"funcotatorValidation"},
          dataProvider = "provideForLargeDataValidationTest")
    public void largeDataValidationTest(final String inputVcfName,
                                        final String referencePath,
                                        final String referenceVersion,
                                        final String dataSourcesPath,
                                        final boolean isDsEnvironmentPath) throws IOException {

        // Get our main test folder path from our environment:
        final String testFolderInputPath = getFuncotatorLargeDataValidationTestInputPath();

        long startTime = 0, endTime = 0;
        final long overallStartTime = System.nanoTime();

        final String outFileBaseName = inputVcfName + ".funcotator";

        final String dataSourcesPathString;
        if (isDsEnvironmentPath) {
            dataSourcesPathString = getFuncotatorLargeDataValidationTestInputPath() + dataSourcesPath;
        }
        else {
            dataSourcesPathString = dataSourcesPath;
        }

        for (final FuncotatorArgumentDefinitions.OutputFormatType outFormat : FuncotatorArgumentDefinitions.OutputFormatType.values()) {

            startTime = System.nanoTime();

            final File outputFile;
            if (outFormat == FuncotatorArgumentDefinitions.OutputFormatType.VCF) {
                outputFile = getOutputFile(outFileBaseName, outFormat.toString().toLowerCase());
            } else {
                outputFile = getOutputFile(outFileBaseName + ".maf", "tsv");
            }

            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    testFolderInputPath + inputVcfName,
                    outputFile,
                    referencePath,
                    dataSourcesPathString,
                    referenceVersion,
                    outFormat,
                    true);

            // Add our manual annotations to the arguments:
            addManualAnnotationsToArguments(arguments);

            // Run the tool with our args:
            runCommandLine(arguments);

            endTime = System.nanoTime();

            logger.warn("  Elapsed Time (" + outFormat.toString() + "): " + (endTime - startTime) / 1e9 + "s");

            // Perform a simple validation if the file was a VCF:
            if ( outFormat == FuncotatorArgumentDefinitions.OutputFormatType.VCF) {

                try ( final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(outputFile.getAbsolutePath()) ) {
                    Assert.assertTrue(vcfReader.getHeader() instanceof VCFHeader, "Header is not a VCFHeader!");
                    final VCFHeader vcfHeader = (VCFHeader)vcfReader.getHeader();

                    final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
                    final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

                    validateFuncotationsOnVcf(vcfReader, funcotationFieldNames);

                }
            }
        }

        logger.warn("Total Elapsed Time: " + (endTime - overallStartTime) / 1e9 + "s");
    }

    @Test(dataProvider = "provideForNonTrivialLargeDataValidationTest")
    public void nonTrivialLargeDataValidationTest(final String inputVcfName,
                                                  final String referencePath,
                                                  final String referenceVersion,
                                                  final String dataSourcesPath,
                                                  final String expectedOutputPath,
                                                  final List<String> excludedFields,
                                                  final boolean reannotate) {

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            // The CLI for Funcotator does not support SEG output.  Must use FuncotateSegments for that.
            if (outputFormatType.equals(FuncotatorArgumentDefinitions.OutputFormatType.SEG)) {
                continue;
            }
            final String typeCorrectedExpectedOutPath = FilenameUtils.removeExtension(expectedOutputPath) + "." + outputFormatType.toString().toLowerCase();

            final File outputFile = createTempFile(tmpOutDir + File.separator + inputVcfName + ".funcotator", "." + outputFormatType.toString().toLowerCase());

            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    inputVcfName,
                    outputFile,
                    referencePath,
                    dataSourcesPath,
                    referenceVersion,
                    outputFormatType,
                    true,
                    excludedFields,
                    reannotate);

            // Run the tool with our args:
            long startTime = 0, endTime = 0;
            startTime = System.nanoTime();
            runCommandLine(arguments);
            endTime = System.nanoTime();

            logger.warn("  " + outputFormatType.toString() + " Elapsed Time: " + (endTime - startTime) / 1e9 + "s");

            // ========================================================
            // Validate our output:

            if ( outputFormatType == FuncotatorArgumentDefinitions.OutputFormatType.VCF ) {
                // Get the actual data:
                assertEqualVariantFiles(outputFile, typeCorrectedExpectedOutPath);
            }
            else {
                try {
                    IntegrationTestSpec.assertEqualTextFiles(outputFile, new File(typeCorrectedExpectedOutPath), "#");
                }
                catch ( final IOException ex ) {
                    throw new GATKException("Error opening expected file: " + expectedOutputPath, ex);
                }
            }
        }
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void exhaustiveArgumentTest(final String dataSourcesPath,
                                       final String refVer,
                                       final String referenceFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                       final String expectedOutputFilePath) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + transcriptName + "." + outputFormatType.toString().toLowerCase();

        final File outputFile = getOutputFile(outFileName.substring(0, outFileName.length() - 4), outFileName.substring(outFileName.length() - 3));

        // Set up our transcript of interest in a file:
        final File transcriptIdFile = getSafeNonExistentFile("TranscriptIdFile.txt");
        try (final BufferedWriter bw = new BufferedWriter(new FileWriter(transcriptIdFile))) {
            bw.write(transcriptName);
        }

        // Required Args:
        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                variantFileName,
                outputFile,
                referenceFileName,
                dataSourcesPath,
                refVer,
                outputFormatType,
                false);

        // Transcript selection:
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.BEST_EFFECT.toString());
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME, transcriptIdFile.toString());

        // Annotation Defaults and Overrides:
        arguments.add(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "GARBAGEDAY:SPUMONI");
        arguments.add(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME, "Oreganno_Build:BUILDED_GOOD_REAL_BIG");

        // Run the beast:
        runCommandLine(arguments);

        // Only test for content-correctness if the output file was specified:
        if (expectedOutputFilePath != null) {
            // Get the expected output file:
            final File expectedOutputFile = new File(expectedOutputFilePath);

            // Make sure that the actual and expected output files are the same:
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedOutputFile, "#");
        }
    }

    /**
     * Test that we can annotate a b37 (here it is clinvar) datasource when GENCODE is using "chr*" and the datasource is not.
     */
    @Test
    public void testCanAnnotateMixedContigHg19Clinvar() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        final int NUM_VARIANTS = 21;
        final int NUM_CLINVAR_HITS = 4;
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS, "Found unexpected number of variants!");

        // Look for "MedGen" to know that we have a clinvar hit.
        Assert.assertEquals(variantContexts.stream()
                .filter(vc -> StringUtils.contains(vc.getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, ""), "MedGen"))
                .count(), NUM_CLINVAR_HITS, "Found unexpected number of ClinVar hits!");
    }

    /**
     * Test that the manual annotations and overrides will be correctly rendered on output, and will occur only once each.
     */
    @Test
    public void testManualAnnotationsCorrectness() {

        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        // Add our manual annotations to the arguments:
        addManualAnnotationsToArguments(arguments);

        // Run the tool:
        runCommandLine(arguments);

        // ===========================================
        // Now let's validate that everything's there:

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();

        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        // Ensure that the field names are all unique:
        final Set<String> fieldNameSet = new HashSet<>(Arrays.asList(funcotationFieldNames));
        Assert.assertEquals(fieldNameSet.size(), funcotationFieldNames.length);

        for ( int i = 0; i < variantContexts.size(); ++i ) {

            final VariantContext vc = variantContexts.get(i);

            final String tx = FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY;
            final Allele altAllele = vc.getAlternateAllele(0);

            // Get our Funcotations in a map:
            final String funcotation = vc.getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, "");
            Assert.assertNotEquals(funcotation, "");
            final FuncotationMap funkyMap =
                    FuncotationMap.createAsAllTableFuncotationsFromVcf(
                            tx,
                            funcotationFieldNames,
                            funcotation,
                            altAllele,
                            "VCF");

            // Verify our manual annotations are in the output:
            for ( final String manualAnnotation : getManualAnnotations() ) {

                // Get our expected field name and value:
                final String expectedFieldName = manualAnnotation.substring(0, manualAnnotation.indexOf(':'));
                final String expectedFieldValue = manualAnnotation.substring(manualAnnotation.indexOf(':')+1, manualAnnotation.length());

                // Now verify the value of the funcotation:
                Assert.assertEquals(funkyMap.getFieldValue(tx, expectedFieldName, altAllele), expectedFieldValue);
            }

            // Verify our annotation overrides are in the output:
            for ( final String annotationOverride : getAnnotationOverrides() ) {

                // Get our expected field name and value:
                final String expectedFieldName = annotationOverride.substring(0, annotationOverride.indexOf(':'));
                final String expectedFieldValue = annotationOverride.substring(annotationOverride.indexOf(':')+1, annotationOverride.length());

                // Now verify the value of the funcotation:
                Assert.assertEquals(funkyMap.getFieldValue(tx, expectedFieldName, altAllele), expectedFieldValue);
            }
        }
    }

    @Test
    public void testXsvLocatableAnnotationsHaveOnlyOneEntryForMultiHitLocations() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                XSV_CLINVAR_MULTIHIT_TEST_VCF,
                outputFile,
                b37Reference,
                DS_XSV_CLINVAR_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final VCFInfoHeaderLine funcotationHeaderLine = vcfInfo.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final int EXPECTED_NUM_VARIANTS = 1;
        Assert.assertEquals(vcfInfo.getRight().size(), EXPECTED_NUM_VARIANTS, "Found more than " + EXPECTED_NUM_VARIANTS + " variants!");

        validateFuncotationsOnVcf(vcfInfo.getRight(), funcotationFieldNames);
    }

    @Test
    public void testXsvLocatableAnnotationsHaveCorrectColsForOnlyOnePositionSpecified() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                XSV_CLINVAR_COL_TEST_VCF,
                outputFile,
                b37Reference,
                DS_XSV_CLINVAR_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final VCFInfoHeaderLine funcotationHeaderLine = vcfInfo.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final int EXPECTED_NUM_VARIANTS = 10;
        Assert.assertEquals(vcfInfo.getRight().size(), EXPECTED_NUM_VARIANTS);

        validateFuncotationsOnVcf(vcfInfo.getRight(), funcotationFieldNames);
    }

    @Test
    public void testCanAnnotateHg38ClinvarAndGencodeV28() {
        // Clinvar datasource did  go through one round of preprocessing to make contig names "1" --> "chr1" (for example).  This is an issue with ClinVar, not GATK.
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG38,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);

        runCommandLine(arguments);

        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        final int NUM_VARIANTS = 100;
        final int NUM_CLINVAR_HITS = 1; // This was verified with SelectVariants
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS);

        // Look for "MedGen" to know that we have a clinvar hit.
        Assert.assertEquals(variantContexts.stream()
                .filter(vc -> StringUtils.contains(vc.getAttributeAsString("FUNCOTATION", ""), "MedGen"))
                .count(), NUM_CLINVAR_HITS);
    }

    @Test
    public void testMANESelectAnnotationDifferencesAndGencodeV43() {
        // Clinvar datasource did  go through one round of preprocessing to make contig names "1" --> "chr1" (for example).  This is an issue with ClinVar, not GATK.
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFileMANESelectMode = getOutputFile(outputFormatType);
        final File outputFileDefaultMode = getOutputFile(outputFormatType);

        // Run the tool twice with FuncotatorArgumentDefinitions.PREFER_MANE_TRANSCRIPT_MODE on and off and compare the results
        final ArgumentsBuilder arguments1 = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG38,
                outputFileMANESelectMode,
                hg38Chr3Ref,
                DS_PIK3CA_DIR_v43_MANE,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);
        arguments1.add(FuncotatorArgumentDefinitions.PREFER_MANE_TRANSCRIPT_MODE, "true");
        arguments1.add("intervals", "chr3:179148357-179235107");
        runCommandLine(arguments1);

        final ArgumentsBuilder arguments2 = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG38,
                outputFileDefaultMode,
                hg38Chr3Ref,
                DS_PIK3CA_DIR_v43_MANE,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);
        arguments2.add("intervals", "chr3:179148357-179235107");
        runCommandLine(arguments2);

        // Assert that we NEVER seen the lower priority transcript when we are in MANE select mode
        try (final FeatureDataSource<VariantContext> featureDataSourceMANE = new FeatureDataSource<>(outputFileMANESelectMode);
             final FeatureDataSource<VariantContext> featureDataSourceDefault = new FeatureDataSource<>(outputFileDefaultMode)) {
            boolean hasSeenENST00000435560 = false;
            for (final VariantContext vc : featureDataSourceMANE) {
                if (vc.getAttributeAsString("FUNCOTATION", null) != null) {
                    final String funcotationMANE = vc.getAttributeAsString("FUNCOTATION", "");
                    Assert.assertFalse(funcotationMANE.contains("ENST00000263967.4")); // this transcript is MANE_SELECT not MANE_PLUS_CLINICAL so should never be selected for these variants
                    Assert.assertTrue(funcotationMANE.contains("ENST00000643187.1") || funcotationMANE.contains("ENST00000435560.1"));
                    if (funcotationMANE.contains("ENST00000435560.1")) {
                        hasSeenENST00000435560 = true;
                    }
                }
            }
            Assert.assertTrue(hasSeenENST00000435560); // ENST00000435560.1 is a basic annotation with no MANE_SELECT transcripts overlapping it, so it should still be emitted
            hasSeenENST00000435560 = false;
            // Assert that we see a mix of transcripts when we are in MANE select mode
            boolean hasSeenMANESelect = false;
            boolean hasSeenMANEPlusClinical = false;
            for (final VariantContext vc : featureDataSourceDefault) {
                if (vc.getAttributeAsString("FUNCOTATION", null) != null) {
                    final String funcotationDefault = vc.getAttributeAsString("FUNCOTATION", "");
                    if (funcotationDefault.contains("ENST00000643187.1")) {
                        hasSeenMANESelect = true;
                    }
                    if (funcotationDefault.contains("ENST00000643187.1")) {
                        hasSeenMANEPlusClinical = true;
                    }
                    if (funcotationDefault.contains("ENST00000435560.1")) {
                        hasSeenENST00000435560 = true;
                    }
                }
            }
            // assert that we have seen both versions of the transcript in the default case
            Assert.assertTrue(hasSeenENST00000435560); // ENST00000435560.1 is a basic annotation with no MANE_SELECT transcripts overlapping it, so it should still be emitted
            Assert.assertTrue(hasSeenMANESelect);
            Assert.assertTrue(hasSeenMANEPlusClinical);
        }
    }

    //Test for https://github.com/broadinstitute/gatk/issues/6173
    @Test
    public void testVCFColumnsArentShuffled() {
        final File outputFile = createTempFile("tmpTestFilterParsing", "vcf");

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                VCF_FIELD_ORDER_SWAP_TEST_VCF,
                outputFile,
                b37Reference,
                VCF_FIELD_ORDER_TEST_DATA_SOURCES,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                FuncotatorArgumentDefinitions.OutputFormatType.VCF,
                false);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> tempVcf =  VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        Assert.assertEquals( tempVcf.getRight().size(), 1 );

        final String[] funcotatorKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(tempVcf.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME).getDescription());

        final VariantContext variantContext = tempVcf.getRight().get(0);
        final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotatorKeys, variantContext, "Gencode_19_annotationTranscript", "FAKE_SOURCE");

        Allele allele = variantContext.getAlternateAllele(0);
        final String txId = funcs.get(allele).getTranscriptList().get(0);
        Assert.assertEquals( funcs.get(allele).get(txId).size(), 1 );

        final Funcotation funcotation = funcs.get(allele).get(txId).get(0);

        //Assert that the value of the field F# is F#|F# encoded in Funcotator's percent encoding scheme
        for(int i = 1; i <= 9; i++){
            Assert.assertEquals(funcotation.getField("dbSnp_F"+i), "F"+i+"_%7C_"+"F"+i);
        }
    }


    @Test
    public void testFilterParsing() {

        final File outputFile = createTempFile("tmpTestFilterParsing", "vcf");

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                FILTER_TEST_VCF,
                outputFile,
                b37Reference,
                DS_FILTER_PARSE_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                FuncotatorArgumentDefinitions.OutputFormatType.VCF,
                false);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> tempVcf =  VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        Assert.assertEquals( tempVcf.getRight().size(), 1 );

        final String[] funcotatorKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(tempVcf.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME).getDescription());

        final VariantContext variantContext = tempVcf.getRight().get(0);
        final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotatorKeys, variantContext, "Gencode_19_annotationTranscript", "FAKE_SOURCE");

        final String txId = funcs.get(variantContext.getAlternateAllele(0)).getTranscriptList().get(0);
        Assert.assertEquals( funcs.get(variantContext.getAlternateAllele(0)).get(txId).size(), 1 );

        final Funcotation funcotation = funcs.get(variantContext.getAlternateAllele(0)).get(txId).get(0);

        Assert.assertEquals(funcotation.getField("dbSnp_FILTER"), "FILTER_8");
    }

    @Test
    public void testExclusionFromDatasourceVcfToVcf() {
        // Clinvar datasource did  go through one round of preprocessing to make contig names "1" --> "chr1" (for example).  This is an issue with ClinVar, not GATK.
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final List<String> excludedFields = Arrays.asList("dummy_ClinVar_VCF_DBVARID", "dummy_ClinVar_VCF_CLNVI");
        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG38,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false, excludedFields, false);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> tempVcf =  VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());

        final String[] funcotatorKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(tempVcf.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME).getDescription());

        // Ensure that the header does not contain the excluded fields
        Stream.of(funcotatorKeys).forEach(k -> Assert.assertFalse(excludedFields.contains(k)));

        final List<VariantContext> variantContexts = tempVcf.getRight();
        for (final VariantContext vc : variantContexts) {
            final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                    funcotatorKeys, vc, "Gencode_28_annotationTranscript", "FAKE_SOURCE");
            for (final String txId: funcs.get(vc.getAlternateAllele(0)).getTranscriptList()) {
                final List<Funcotation> funcotations = funcs.get(vc.getAlternateAllele(0)).get(txId);
                for (final Funcotation funcotation : funcotations) {
                    funcotation.getFieldNames().forEach(f -> Assert.assertFalse(excludedFields.contains(f)));
                }
            }
        }
    }

    @DataProvider(name = "provideForMafVcfConcordance")
    final Object[][] provideForMafVcfConcordance() {
        return new Object[][]{
                {PIK3CA_VCF_HG19_SNPS, b37Reference, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList("Gencode_19_proteinChange"), Collections.singletonList(MafOutputRendererConstants.FieldName_Protein_Change), DS_PIK3CA_DIR, true, 15},
                {PIK3CA_VCF_HG19_INDELS, b37Reference, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList("Gencode_19_proteinChange"), Collections.singletonList(MafOutputRendererConstants.FieldName_Protein_Change), DS_PIK3CA_DIR, true, 57},
                {MUC16_VCF_HG19, hg19Chr19Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList("Gencode_19_proteinChange"), Collections.singletonList(MafOutputRendererConstants.FieldName_Protein_Change), FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, false, 2057},
                {
                        PIK3CA_VCF_HG38,
                        hg38Chr3Ref,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        VCF_FIELDS_GENCODE_28_DS,
                        MAF_FIELDS_GENCODE_DS,
                        DS_PIK3CA_DIR,
                        false,
                        104,
                },
                {PIK3CA_VCF_HG19_INDELS, b37Reference, FuncotatorTestConstants.REFERENCE_VERSION_HG19, VCF_FIELDS_GENCODE_19_DS, MAF_FIELDS_GENCODE_DS, DS_PIK3CA_DIR, true, 57},
        };
    }

    private void createConfigFileForMAF(final File mafConfigFile) {
        try ( final PrintWriter printWriter = new PrintWriter(mafConfigFile) ) {
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_CONTIG_COLUMN + " = " + MafOutputRendererConstants.FieldName_Chromosome);
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_START_COLUMN + " = " + MafOutputRendererConstants.FieldName_Start_Position);
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_END_COLUMN + " = " + MafOutputRendererConstants.FieldName_End_Position);
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_DELIMITER + " = \\t");
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_NAME + " = ");
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_SRC_FILE + " = ");
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_VERSION + " = ");
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_ORIGIN_LOCATION + " = ");
            printWriter.println(DataSourceUtils.CONFIG_FILE_FIELD_NAME_PREPROCESSING_SCRIPT + " = ");

        }
        catch (final FileNotFoundException ex) {
            throw new GATKException("Could not create the tmp config file to test maf/vcf concorance: " + mafConfigFile.toURI().toString(), ex);
        }
    }

    /**
     * Make sure that VCFs and MAFs have exactly the same annotation strings.  This test does not look for
     *  multiallelics.
     */
    @Test(dataProvider = "provideForMafVcfConcordance")
    public void testVcfMafConcordance(final String inputVcf,
                                      final String inputRef,
                                      final String funcotatorRef,
                                      final List<String> annotationsToCheckVcf,
                                      final List<String> annotationsToCheckMaf,
                                      final String datasourceDir,
                                      final boolean forceB37Hg19Conversion,
                                      final int gtNumVariants) {

        // ===========================================================================================
        // Run in VCF Mode:
        // ===========================================
        final FuncotatorArgumentDefinitions.OutputFormatType vcfOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File vcfOutputFile = getOutputFile(vcfOutputFormatType);

        final ArgumentsBuilder argumentsVcf = createBaselineArgumentsForFuncotator(
                inputVcf,
                vcfOutputFile,
                inputRef,
                datasourceDir,
                funcotatorRef,
                vcfOutputFormatType,
                false);

        argumentsVcf.add(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);
        argumentsVcf.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        if ( forceB37Hg19Conversion ) {
            // We need this argument since we are testing on a subset of b37
            argumentsVcf.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        }

        runCommandLine(argumentsVcf);

        // ===========================================================================================
        // Run in MAF Mode:
        // ===========================================

        final FuncotatorArgumentDefinitions.OutputFormatType mafOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File mafOutputFile = getOutputFile(mafOutputFormatType);

        final ArgumentsBuilder argumentsMaf = createBaselineArgumentsForFuncotator(
                inputVcf,
                mafOutputFile,
                inputRef,
                datasourceDir,
                funcotatorRef,
                mafOutputFormatType,
                false);

        argumentsMaf.add(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);
        argumentsMaf.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        if ( forceB37Hg19Conversion ) {
            // We need this argument since we are testing on a subset of b37
            argumentsMaf.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        }

        runCommandLine(argumentsMaf);

        // TODO: Create another config file for this MAF that has good column headers in it.
        final File mafConfigFile = createTempFile(mafOutputFile.getName(), ".config");
        createConfigFileForMAF(mafConfigFile);

        // ===========================================================================================
        // Read in and validate VCF:
        // ===========================================

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(vcfOutputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        Assert.assertEquals(variantContexts.stream().map(vc -> vc.getAlleles().size() - 1).mapToInt(Integer::intValue).sum(), gtNumVariants);
        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.hasAttribute(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME)));

        // ===========================================================================================
        // Read in and validate MAF:
        // ===========================================

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(mafOutputFile.toPath(), mafConfigFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(), gtNumVariants);

        // Some errors manifest as all of the variant classifications being IGR.  Check to make sure that is not the case.
        Assert.assertTrue(maf.getRecords().stream()
                .anyMatch(v -> !v.getAnnotationValue(MafOutputRendererConstants.FieldName_Variant_Classification).equals("IGR")), "Output produced only IGR annotations!");

        // ===========================================================================================
        // Compare VCF and MAF Annotations:
        // ===========================================

        // Get the annotation fields:
        final String[] funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final Map<String, Pair<Pair<String, String>, Pair<String, String>>> unequalFieldValuesMap = new HashMap<>();

        for (int i = 0; i < annotationsToCheckMaf.size(); i++) {
            final String annotationToCheckVcf = annotationsToCheckVcf.get(i);
            final String annotationToCheckMaf = annotationsToCheckMaf.get(i);

            // Have to get the contig / start / end from the interval itself for MAF:
            List<String> mafFieldValues;
            switch ( annotationsToCheckMaf.get(i) ) {
                case MafOutputRendererConstants.FieldName_Chromosome:
                    mafFieldValues = maf.getRecords().stream().map(AnnotatedInterval::getContig).collect(Collectors.toList());
                    break;
                case MafOutputRendererConstants.FieldName_Start_Position:
                    mafFieldValues = maf.getRecords().stream().map(AnnotatedInterval::getStart).map(Object::toString).collect(Collectors.toList());
                    break;
                case MafOutputRendererConstants.FieldName_End_Position:
                    mafFieldValues = maf.getRecords().stream().map(AnnotatedInterval::getEnd).map(Object::toString).collect(Collectors.toList());
                    break;
                case MafOutputRendererConstants.FieldName_Hugo_Symbol:
                    mafFieldValues = maf.getRecords().stream().map(x -> x.getAnnotationValue(MafOutputRendererConstants.FieldName_Hugo_Symbol)).map(a -> a.isEmpty() ? "Unknown" : a).collect(Collectors.toList());
                    break;
                default:
                    mafFieldValues = maf.getRecords().stream().map(v -> v.getAnnotationValue(annotationToCheckMaf)).collect(Collectors.toList());
                    break;
            }
            mafFieldValues = mafFieldValues.stream().map(val -> MafOutputRenderer.mafTransformInvert(annotationToCheckMaf, val, funcotatorRef)).collect(Collectors.toList());

            // Note that we assume that each variant context has one allele and one transcript.  This is true due to the
            // datasources and input VCF.
            // Don't try to refactor this for-loop to a stream here.
            final List<String> vcfFieldValues = new ArrayList<>();
            for (final VariantContext v: variantContexts) {

                final String transcriptFuncotationName = "Gencode_" + (funcotatorRef.equals("hg19") ? "19" : "28") + "_annotationTranscript";

                final Map<Allele, FuncotationMap> alleleFuncotationMapMap = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                        funcotationKeys, v, transcriptFuncotationName, "TEST");

                for (final Allele alternateAllele : v.getAlternateAlleles() ) {
                    final FuncotationMap funcotationMap = alleleFuncotationMapMap.get(alternateAllele);
                    vcfFieldValues.add(funcotationMap.getFieldValue(funcotationMap.getTranscriptList().get(0), annotationToCheckVcf, alternateAllele));
                }
            }

            // "Fuzzy" matching between VCF and MAF by allowing blank VCF fields to be Unknown in the MAF.
            // This is a temporary hold-over until we can do full aliased comparisons of fields.
            // TODO: Update with better comparison using aliases when they are fully implemented.
            if ( !Objects.equals(mafFieldValues, vcfFieldValues) ) {
                Assert.assertEquals(mafFieldValues.size(), vcfFieldValues.size());

                int numDifferentValues = 0;
                String diffVcfValueString = "";
                String diffMafValueString = "";

                for ( int valueIndex = 0 ; valueIndex < mafFieldValues.size(); ++valueIndex ) {

                    final String vcfValue = vcfFieldValues.get(valueIndex);
                    final String mafValue = mafFieldValues.get(valueIndex);

                    if ( mafValue.isEmpty() ) {
                        if (!vcfValue.equals("Unknown")) {
                            diffVcfValueString = diffVcfValueString + "\t" + vcfValue + "[" + valueIndex + "]";
                            diffMafValueString = diffMafValueString + "\t" + mafValue + "[" + valueIndex + "]";
                            ++numDifferentValues;
                        }
                    }
                    else {
                        if (!vcfValue.equals(mafValue)) {
                            diffVcfValueString = diffVcfValueString + "\t" + vcfValue + "[" + valueIndex + "]";
                            diffMafValueString = diffMafValueString + "\t" + mafValue + "[" + valueIndex + "]";
                            ++numDifferentValues;
                        }
                    }
                }

                // Special case for end coordinates:
                if ( annotationToCheckMaf.equals("End_Position") ) {
                    // End position is post-processed in MAF for indels,
                    // and therefore is not expected to always be equal to the same value as in the VCF:
                    logger.info("End positions are not the same (this is OK): \nVCF:\t" + diffVcfValueString + "\nMAF:\t" + diffMafValueString);
                }
                else {
                    // Add our info to our map:
                    unequalFieldValuesMap.put(annotationToCheckVcf,
                            Pair.of(
                                    Pair.of(annotationToCheckVcf, diffVcfValueString),
                                    Pair.of(annotationToCheckMaf, diffMafValueString)
                            ));
                }
            }
        }

        // We have had a problem.  Alert the user in as helpful a manner as possible:
        if ( !unequalFieldValuesMap.isEmpty() ) {
            final StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.append("Failed Matching VCF and MAF fields:\n");
            for ( final Map.Entry<String, Pair<Pair<String, String>, Pair<String, String>>> entry : unequalFieldValuesMap.entrySet() ) {

                final int formatLength = 5 + Math.max(entry.getValue().getLeft().getLeft().length(), entry.getValue().getRight().getLeft().length());
                final String formatString = "\t%s (%-" + formatLength + "s";

                stringBuilder.append( String.format(formatString, "VCF", entry.getValue().getLeft().getLeft() + "):") );
                stringBuilder.append( entry.getValue().getLeft().getRight() );
                stringBuilder.append('\n');

                stringBuilder.append( String.format(formatString, "MAF", entry.getValue().getRight().getLeft() + "):") );
                stringBuilder.append( entry.getValue().getRight().getRight() );
                stringBuilder.append('\n');

                stringBuilder.append("----\n");
            }

            throw new AssertionError(stringBuilder.toString());
        }
    }

    @Test
    public void testVcfDatasourceAccountsForAltAlleles() {
        final FuncotatorArgumentDefinitions.OutputFormatType vcfOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File vcfOutputFile = getOutputFile(vcfOutputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19_ALTS,
                vcfOutputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                vcfOutputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(vcfOutputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        Assert.assertTrue(variantContexts.size() > 0);
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        // The first variant context should have clinvar annotations, since it hit on the alt allele.  None of the rest.
        // This test assumes that each test variant context has only one alt allele.
        // The rest should not have any clinvar hits.
        for (int i = 0; i < variantContexts.size(); i++) {
            final String gtString = (i == 0) ? "MedGen:C0027672,SNOMED_CT:699346009" : "";
            final Map<Allele, FuncotationMap> alleleToFuncotationMap = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(funcotationKeys,
                    variantContexts.get(i), "Gencode_19_annotationTranscript", "TEST");
            Assert.assertEquals(alleleToFuncotationMap.entrySet().size(), 1, "Found more than 1 alternate allele!");

            final FuncotationMap funcotationMap = alleleToFuncotationMap.values().iterator().next();
            Assert.assertEquals(funcotationMap.getTranscriptList().size(), 1, "Found more than 1 funcotation!");
            Assert.assertTrue(funcotationMap.getTranscriptList().stream().noneMatch(k -> k.equals(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY)), "Found a funcotation with an unknown transcript name: " + FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY);
            Assert.assertTrue(funcotationMap.getTranscriptList().stream().noneMatch(StringUtils::isEmpty), "Found a funcotation with an empty transcript!");
            final List<Funcotation> funcotations = funcotationMap.get(funcotationMap.getTranscriptList().get(0));
            Assert.assertEquals(funcotations.size(), 1, "Found more than one funcotation in the funcotation map!");
            final Funcotation funcotation = funcotations.get(0);

            Assert.assertEquals(funcotation.getField("dummy_ClinVar_VCF_CLNDISDB"), FuncotatorUtils.sanitizeFuncotationFieldForVcf(gtString), "Field (dummy_ClinVar_VCF_CLNDISDB) was unsanititzed: " + funcotation.getField("dummy_ClinVar_VCF_CLNDISDB"));
        }
    }

    @Test
    public void testCanAnnotateSpanningDeletions() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                SPANNING_DEL_VCF,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        Assert.assertTrue(variantContexts.size() > 0);
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final int NUM_VARIANTS = 10;
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS);
        Assert.assertEquals(variantContexts.stream().filter(v -> v.hasAllele(Allele.SPAN_DEL)).count(), NUM_VARIANTS);

        for (final VariantContext vc : variantContexts) {
            final Map<Allele, FuncotationMap> alleleToFuncotationMap = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(funcotationKeys,
                    vc, "Gencode_28_annotationTranscript", "TEST");

            // Make sure that all spanning deletions are funcotated with could not determine and that the others are something else.
            for (final Allele allele : alleleToFuncotationMap.keySet()) {
                final List<String> txIds = alleleToFuncotationMap.get(allele).getTranscriptList();
                for (final String txId : txIds) {
                    if (allele.equals(Allele.SPAN_DEL)) {
                        Assert.assertEquals(alleleToFuncotationMap.get(allele).get(txId).get(0).getField("Gencode_28_variantClassification"), GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE.toString());
                        Assert.assertEquals(alleleToFuncotationMap.get(allele).get(txId).get(0).getField("Gencode_28_variantType"), GencodeFuncotation.VariantType.NA.toString());
                    } else {
                        Assert.assertNotEquals(alleleToFuncotationMap.get(allele).get(txId).get(0).getField("Gencode_28_variantClassification"), GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE.toString());
                    }
                }
            }
        }
    }

    @Test
    public void testNoSpanningDeletionWriteWithMAF() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                SPANNING_DEL_VCF,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);

        runCommandLine(arguments);

        // There should only be one variant in the MAF, not two.
        //  TODO:  This input VCF has an "END" field which (currently) throws off the reading of an AnnotatedIntervalCollection.  The MAF_TEST_CONFIG is a workaround/hack.  See issue https://github.com/broadinstitute/gatk/issues/4897
        final AnnotatedIntervalCollection annotatedIntervalCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), Paths.get(MAF_TEST_CONFIG), null);
        Assert.assertEquals(annotatedIntervalCollection.getRecords().size(), 10);
    }

    @Test
    public void testVCFToMAFPreservesFields() {

        final AnnotatedIntervalCollection maf = runPik3caHg19VcfToMaf(new HashSet<>());
        Assert.assertTrue(maf.getRecords().size() > 0);
        Assert.assertTrue(maf.getRecords().stream().allMatch(r -> r.hasAnnotation("ILLUMINA_BUILD")));
        Assert.assertTrue(maf.getRecords().stream().allMatch(r -> r.getAnnotationValue("ILLUMINA_BUILD").startsWith("37")));

        // Needs to get aliases from the MAF, since AF (and maybe more) has its name changed.  So create a dummy
        //  MafOutputRenderer that mimics the one that is used in the command line invocation above and get the aliases.
        final File dummyOutputFile = getOutputFile(FuncotatorArgumentDefinitions.OutputFormatType.MAF);
        final MafOutputRenderer dummyMafOutputRenderer = new MafOutputRenderer(dummyOutputFile.toPath(), Collections.emptyList(), new VCFHeader(), new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), "b37", new HashSet<String>(), "Unknown");
        final Map<String, Set<String>> mafAliasMap = dummyMafOutputRenderer.getReverseOutputFieldNameMap();

        // Get all of the alias lists
        final Pair<VCFHeader, List<VariantContext>> vcf  = VariantContextTestUtils.readEntireVCFIntoMemory(PIK3CA_VCF_HG19);
        final Set<String> vcfHeaderInfoSet = vcf.getLeft().getInfoHeaderLines().stream()
                .map(VCFCompoundHeaderLine::getID)
                .map(s -> mafAliasMap.getOrDefault(s, Collections.singleton(s)))
                .flatMap(Set::stream)
                .collect(Collectors.toSet());
        Assert.assertTrue(vcfHeaderInfoSet.size() > 0);
        Assert.assertTrue(maf.getAnnotations().containsAll(vcfHeaderInfoSet));
    }

    private AnnotatedIntervalCollection runPik3caHg19VcfToMaf(final Set<String> excludedFields) {
        final File outputFile = getOutputFile(FuncotatorArgumentDefinitions.OutputFormatType.MAF);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                FuncotatorArgumentDefinitions.OutputFormatType.MAF,
                false);

        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        excludedFields.forEach(f -> arguments.add(FuncotatorArgumentDefinitions.EXCLUSION_FIELDS_LONG_NAME, f));

        // Disable the sequence dictionary check for the tests:
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        return AnnotatedIntervalCollection.create(outputFile.toPath(), null);
    }

    @Test
    public void testVcfToMafHonorsExcludedFields() {
        final String fieldToEnsureIsIncluded = "dummy_ClinVar_VCF_CLNVC";
        final HashSet<String> excludedFields = com.google.common.collect.Sets.newHashSet("dummy_ClinVar_VCF_AF_EXAC", "dummy_ClinVar_VCF_CLNSIGCONF");
        final AnnotatedIntervalCollection maf = runPik3caHg19VcfToMaf(excludedFields);
        Assert.assertTrue(maf.getRecords().size() > 0);
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation(fieldToEnsureIsIncluded)));
        maf.getRecords().forEach(r -> Assert.assertEquals(Sets.intersection(r.getAnnotations().keySet(), excludedFields).size(), 0));
    }

    @Test
    public void testVCFToVCFPreservesFields() {

        final Set<String> PIK3CA_VCF_HG19_INPUT_FIELDS = new HashSet<>(Arrays.asList("FUNCOTATION", "AC", "AF", "ALLELE_A", "ALLELE_B", "AN", "BEADSET_ID", "DP", "GC_SCORE", "ILLUMINA_BUILD", "ILLUMINA_CHR", "ILLUMINA_POS", "ILLUMINA_STRAND", "N_AA", "N_AB", "N_BB", "PROBE_A", "PROBE_B", "SOURCE", "devR_AA", "devR_AB", "devR_BB", "devTHETA_AA", "devTHETA_AB", "devTHETA_BB", "devX_AA", "devX_AB", "devX_BB", "devY_AA", "devY_AB", "devY_BB", "meanR_AA", "meanR_AB", "meanR_BB", "meanTHETA_AA", "meanTHETA_AB", "meanTHETA_BB", "meanX_AA", "meanX_AB", "meanX_BB", "meanY_AA", "meanY_AB", "meanY_BB", "refSNP", "zthresh_X", "zthresh_Y"));

        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        // We need this argument since we are testing on a subset of b37
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcf  = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());

        final Set<String> vcfFields = vcf.getLeft().getInfoHeaderLines().stream()
                .map(VCFCompoundHeaderLine::getID)
                .collect(Collectors.toSet());

        final Set<String> missingFields = Sets.difference(PIK3CA_VCF_HG19_INPUT_FIELDS, vcfFields);
        final Set<String> additionalFields = Sets.difference(vcfFields, PIK3CA_VCF_HG19_INPUT_FIELDS);

        Assert.assertTrue(missingFields.size() == 0, "Fields were missing in the output: " + missingFields.stream().collect(Collectors.joining(", ")));
        Assert.assertTrue(additionalFields.size() == 0, "Fields were added in the output: " + additionalFields.stream().collect(Collectors.joining(", ")));

        final List<VariantContext> variantContexts = vcf.getRight();

        Assert.assertTrue(variantContexts.size() > 0);
        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.hasAttribute("ILLUMINA_BUILD")));
        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.getAttributeAsString("ILLUMINA_BUILD", "").startsWith("37")));
        final VCFInfoHeaderLine funcotationHeaderLine = vcf.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        // Make sure that none of the input fields appear in the funcotation header.
        final Set<String> funcotationKeys = new HashSet<>(Arrays.asList(FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription())));
        Assert.assertEquals(Sets.intersection(funcotationKeys, PIK3CA_VCF_HG19_INPUT_FIELDS).size(), 0);
    }

    @DataProvider
    public Object[][] provideTNVcfs() {
        // These two VCFs are exactly the same, except how the sample names are handled.
        return new Object[][] {
                {NOT_M2_TEST_HG19},
                {M2_TEST_HG19}
        };
    }
    @Test(dataProvider = "provideTNVcfs")
    public void testMafCustomCountFields(final String tnVcf) {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                tnVcf,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Painstakingly lifted from the input VCF file
        final int[] tRefCounts = new int[]{42, 44, 117, 16, 24, 23};
        final int[] tAltCounts = new int[]{3, 2, 4, 6, 8, 7};
        final int[] nRefCounts = new int[]{43, 45, 145, 22, 45, 13};
        final int[] nAltCounts = new int[]{1, 1, 4, 0, 0, 0};
        final double[] tumorF = new double[]{0.075, 0.049, 0.034, 0.278, 0.259, 0.19};

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_alt_count)),
                   tAltCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_ref_count)),
                    tRefCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_ref_count)),
                    nRefCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_alt_count)),
                    nAltCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Double.parseDouble(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_tumor_f)),
                    tumorF[i], "Record did not match GT at entry " + i));
    }
    @DataProvider
    public Object[][] provideTOnlyVcfs() {
        // These two VCFs are exactly the same, except how the sample names are handled.
        return new Object[][] {
                {NOT_M2_TEST_HG19_TUMOR_ONLY},
                {M2_TEST_HG19_TUMOR_ONLY}
        };
    }
    @Test(dataProvider = "provideTOnlyVcfs")
    public void testMafCustomCountFieldsTumorOnly(final String tnVcf) {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                tnVcf,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Painstakingly lifted from the input VCF file
        final int[] tRefCounts = new int[]{42, 44, 117, 16, 24, 23};
        final int[] tAltCounts = new int[]{3, 2, 4, 6, 8, 7};
        final double[] tumorF = new double[]{0.075, 0.049, 0.034, 0.278, 0.259, 0.19};

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_alt_count)),
                   tAltCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_ref_count)),
                    tRefCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_ref_count),
                    "", "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_alt_count),
                        "","Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Double.parseDouble(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_tumor_f)),
                    tumorF[i], "Record did not match GT at entry " + i));
    }

    /**
     * Test that an input file with more than one possible tumor-normal pairing fails.
     */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testMoreThanOneTNPair() {
        runSomaticVcf(THREE_SAMPLE_SOMATIC);
    }

    private void runSomaticVcf(final String vcfFile) {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                vcfFile,
                outputFile,
                b37Reference,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testSequenceDictionaryCheck() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        // The input VCF and Reference file are incompatible because
        // the reference file dictionary has only chromosome 2 and the
        // input VCF has a dictionary that contains all contigs for HG19.
        // Therefore, the reference dictionary is NOT a superset of the input VCF dictionary.

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                XSV_CLINVAR_COL_TEST_VCF,
                outputFile,
                b37Chr2Ref,
                DS_XSV_CLINVAR_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                true);

        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);
    }

    @Test
    public void testEColiFuncotations() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(FuncotatorTestConstants.ECOLI_VCF_FILE_NAME));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(eColiRef));
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_ECOLI_DIR);
        arguments.add(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_ECOLI);
        arguments.add(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());
        runCommandLine(arguments);
        assertEqualVariantFiles(outputFile, E_COLI_EXPECTED_OUT);
    }

    private void assertEqualVariantFiles(final File outputFile, final String eColiExpectedOut) {
        // Get the actual data:
        final Pair<VCFHeader, List<VariantContext>> actualVcfInfo               = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final List<VariantContext>                  actualVariantContexts       = actualVcfInfo.getRight();
        final VCFHeader                             actualVcfHeader             = actualVcfInfo.getLeft();
        final VCFInfoHeaderLine                     actualFuncotationHeaderLine = actualVcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        // Get the expected data:
        final Pair<VCFHeader, List<VariantContext>> expectedVcfInfo               = VariantContextTestUtils.readEntireVCFIntoMemory(new File(eColiExpectedOut).getAbsolutePath());
        final List<VariantContext>                  expectedVariantContexts       = expectedVcfInfo.getRight();
        final VCFHeader                             expectedVcfHeader             = expectedVcfInfo.getLeft();
        final VCFInfoHeaderLine                     expectedFuncotationHeaderLine = expectedVcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        // Check that they're equal:
        Assert.assertEquals(actualFuncotationHeaderLine, expectedFuncotationHeaderLine);
        VariantContextTestUtils.assertEqualVariants(actualVariantContexts, expectedVariantContexts);
    }

    @Test
    public void testNoVariantsProduceMaf() {
        // Make sure that a MAF is actually produced (not a blank file).  Testing https://github.com/broadinstitute/gatk/issues/4937
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(EMPTY_VCF));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(b37Reference));
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);
        arguments.add(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG19);
        arguments.add(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());
        arguments.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);
        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(),  0);
        assertCustomFieldsArePresent(maf);
    }

    private void assertCustomFieldsArePresent(final AnnotatedIntervalCollection maf) {
        // Double-check that the custom MAF fields are present.
        Assert.assertTrue(CustomMafFuncotationCreator.COUNT_FIELD_NAMES.stream()
                .allMatch(f -> maf.getAnnotations().contains(f)));

        Assert.assertTrue(maf.getAnnotations().contains(MafOutputRendererConstants.FieldName_dbSNP_Val_Status));
    }

    @Test
    public void testEnsureDbSnpInMaf() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(MAF_DBSNP_TEST));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(b37Reference));
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, PIK3CA_DBSNP_DS);
        arguments.add(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG19);
        arguments.add(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());
        arguments.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);
        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(),  5);
        assertCustomFieldsArePresent(maf);

        final List<String> gtDbSnpValStatus = Arrays.asList("",
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG,
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG,
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG,
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG);
        final List<String> guessDbSnpValStatus = maf.getRecords().stream().map(r -> r.getAnnotationValue(MafOutputRendererConstants.FieldName_dbSNP_Val_Status))
                .collect(Collectors.toList());
        Assert.assertEquals(guessDbSnpValStatus, gtDbSnpValStatus);
    }

    @Test
    public void testCanCreateNonLocatableFuncotations() {

        final File outputFile = createTempFile(tmpOutDir + File.separator + NON_LOCATABLE_FUNCOTATED_INPUT_VCF + ".funcotator", ".vcf");

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                NON_LOCATABLE_FUNCOTATED_INPUT_VCF,
                outputFile,
                b37Reference,
                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                FuncotatorArgumentDefinitions.OutputFormatType.VCF,
                true);

        // Run the tool with our args:
        runCommandLine(arguments);

        // ===============================

        final Pair<VCFHeader, List<VariantContext>> vcfInfo
                = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());

        final VCFInfoHeaderLine funcotationHeaderLine = vcfInfo.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final VariantContext variant = vcfInfo.getRight().get(0);

        final Map<Allele, FuncotationMap> alleleToFuncotationMap =
                FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                        funcotationFieldNames,
                        variant,
                        "Gencode_19_annotationTranscript",
                        "TEST");

        // Make sure we get the correct transcript here:
        Assert.assertEquals(alleleToFuncotationMap.get(variant.getAlternateAllele(0)).getTranscriptList().size(), 1);

        // Now get the transcript annotations:
        final List<Funcotation> funcotations = alleleToFuncotationMap.get(variant.getAlternateAllele(0)).get("ENST00000379410.3");

        // Now assert that we got what we should have gotten from a few HGNC (non-locatable data source) fields:
        Assert.assertEquals(funcotations.get(0).getField("HGNC_HGNC_ID"), "HGNC:25284");
        Assert.assertEquals(funcotations.get(0).getField("HGNC_Status"), "Approved");
        Assert.assertEquals(funcotations.get(0).getField("HGNC_Locus_Type"), "gene_%20_with_%20_protein_%20_product");
        Assert.assertEquals(funcotations.get(0).getField("HGNC_Locus_Group"), "protein-coding_%20_gene");
        Assert.assertEquals(funcotations.get(0).getField("HGNC_Previous_Name"), "\"pleckstrin_%20_homology_%20_domain_%20_containing_%2C__%20_family_%20_N_%20_member_%20_1\"");
        Assert.assertEquals(funcotations.get(0).getField("HGNC_Synonyms"), "DKFZP434H2010");
    }

    @Test(expectedExceptions = UserException.class)
    public void testUserExceptionOnAlleleDepthFieldSizeOneForMafOutput() {
        final File outputFile = createTempFile(tmpOutDir + File.separator + SYMBOLLIC_ALLELE_FUNCOTATED_INPUT_VCF + ".funcotator", ".vcf");

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                FuncotatorTestConstants.BAD_DATA_ONLY_ONE_AD_FIELD_VALUE,
                outputFile,
                hg38Reference,
                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                FuncotatorArgumentDefinitions.OutputFormatType.MAF,
                true);

        // Run the tool with our args:
        runCommandLine(arguments);
    }

    @Test
    public void testCanHandleSymbollicAlleleFuncotations() {
        // NOTE: This is an integration test because of the plumbing to convert B37 variants into compatible intervals
        //       for the HG19 data sources.

        final File outputFile = createTempFile(tmpOutDir + File.separator + SYMBOLLIC_ALLELE_FUNCOTATED_INPUT_VCF + ".funcotator", ".vcf");

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                SYMBOLLIC_ALLELE_FUNCOTATED_INPUT_VCF,
                outputFile,
                b37Reference,
                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                FuncotatorArgumentDefinitions.OutputFormatType.VCF,
                true);

        // Run the tool with our args:
        runCommandLine(arguments);

        // ===============================

        final Pair<VCFHeader, List<VariantContext>> vcfInfo
                = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());

        final VCFInfoHeaderLine funcotationHeaderLine = vcfInfo.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());


        // Set up our expected values:
        final String[][] expectedValues = new String[][]{
                {
                    "ENST00000453464.2",
                    "RNF223",
                    "COULD_NOT_DETERMINE"
                },
                {
                    "ENST00000400809.3",
                    "CCNL2",
                    "FIVE_PRIME_FLANK"
                },
        };

        // Now we validate our output:

        Assert.assertEquals(vcfInfo.getRight().size(), 2);

        for ( int i = 0; i < vcfInfo.getRight().size() ; ++i ) {

            final VariantContext variant = vcfInfo.getRight().get(i);
            final String[] expected = expectedValues[i];

            final Map<Allele, FuncotationMap> alleleToFuncotationMap =
                    FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                            funcotationFieldNames,
                            variant,
                            "Gencode_19_annotationTranscript",
                            "TEST");

            // Make sure we get the correct transcript here:
            Assert.assertEquals(alleleToFuncotationMap.get(variant.getAlternateAllele(0)).getTranscriptList().size(), 1);

            // Now get the transcript annotations:
            final List<Funcotation> funcotations = alleleToFuncotationMap.get(variant.getAlternateAllele(0)).get(expected[0]);

            Assert.assertEquals(funcotations.size(), 1);

            // Now assert that we got what we should have gotten:
            Assert.assertEquals(funcotations.get(0).getField("Gencode_19_hugoSymbol"), expected[1]);
            Assert.assertEquals(funcotations.get(0).getField("Gencode_19_variantClassification"), expected[2]);
        }
    }

    @DataProvider
    public Object[][] provideForUnannotatedInputWithOverrideArgument() {
        return new Object[][] {
                {
                    FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1_EXPECTED_OUTPUT
                }
        };
    }
    /**
     * An integration test for the Funcotator tool, and specifically,
     * for the function checkIfAlreadyAnnotated which is included in this
     * tool. In this test, we are giving as input a vcf which has not already been
     * annotated, but we are also using the flag "true" to indicate that we want it
     * to be reannotated to make sure that this does not cause a problem.
     * Created by Hailey on 7/9/21.
     */
    @Test(dataProvider = "provideForUnannotatedInputWithOverrideArgument")
    public void testUnannotatedInputWithOverrideArgument(final String inputVcfName,
                                                      final String referencePath,
                                                      final String referenceVersion,
                                                      final String dataSourcesPath,
                                                      final String expectedOutputPath ) {

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType :
        FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            if (outputFormatType.equals(FuncotatorArgumentDefinitions.OutputFormatType.SEG)) {
                continue;
            }
            final String typeCorrectedExpectedOutPath = FilenameUtils.removeExtension(expectedOutputPath) +
                    "." + outputFormatType.toString().toLowerCase();
            final File outputFile = createTempFile(tmpOutDir + File.separator + inputVcfName + ".funcotator", "." + outputFormatType.toString().toLowerCase());

            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    inputVcfName,
                    outputFile,
                    referencePath,
                    dataSourcesPath,
                    referenceVersion,
                    outputFormatType,
                    true,
                    Collections.emptyList(),
                    true);
            runCommandLine(arguments);
        }

     }

    @DataProvider
    public Object[][] provideForAlreadyAnnotatedInput() {
        return new Object[][] {
                {
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1_EXPECTED_OUTPUT,
                        b37Reference,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                        FuncotatorTestConstants.NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1_EXPECTED_OUTPUT,

                }
        };
    }


    /** In this test, we are giving as input a vcf which has already been annotated.
     * Because we are not giving the flag argument which indicates that we should reannotate,
     * then the default argument will be used which is to no reannotate and instead to throw
     * an error.
     */
    @Test(dataProvider = "provideForAlreadyAnnotatedInput", expectedExceptions = UserException.BadInput.class)
    public void testAlreadyAnnotatedInputWithoutOverrideArgument(final String inputVcfName,
                                                           final String referencePath,
                                                           final String referenceVersion,
                                                           final String dataSourcesPath,
                                                           final String expectedOutputPath)
    {

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType :
                FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            if (outputFormatType.equals(FuncotatorArgumentDefinitions.OutputFormatType.SEG)) {
                continue;
            }
            final String typeCorrectedExpectedOutPath = FilenameUtils.removeExtension(expectedOutputPath) +
                    "." + outputFormatType.toString().toLowerCase();
            final File outputFile = createTempFile(tmpOutDir + File.separator + inputVcfName + ".funcotator", "." + outputFormatType.toString().toLowerCase());

            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    inputVcfName,
                    outputFile,
                    referencePath,
                    dataSourcesPath,
                    referenceVersion,
                    outputFormatType,
                    true);
            runCommandLine(arguments);
        }

    }

    /** In this instance, we are again giving an already annotated file, but
     * the --reannotate-vcf argument is being turned on using the flag "true"
     * to indicate that we want this file to be reannotated. There should be no errors.
     */
    @Test(dataProvider = "provideForAlreadyAnnotatedInput")
    public void testAlreadyAnnotatedInputWithOverrideArgument(final String inputVcfName,
                                                final String referencePath,
                                                final String referenceVersion,
                                                final String dataSourcesPath,
                                                final String expectedOutputPath)
    {

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType :
                FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            if (outputFormatType.equals(FuncotatorArgumentDefinitions.OutputFormatType.SEG)) {
                continue;
            }
            final String typeCorrectedExpectedOutPath = FilenameUtils.removeExtension(expectedOutputPath) +
                    "." + outputFormatType.toString().toLowerCase();
            final File outputFile = createTempFile(tmpOutDir + File.separator + inputVcfName + ".funcotator", "." + outputFormatType.toString().toLowerCase());

            //the --reannotate-vcf is being turned on here using the "true" flag
            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    inputVcfName,
                    outputFile,
                    referencePath,
                    dataSourcesPath,
                    referenceVersion,
                    outputFormatType,
                    true,
                    Collections.emptyList(),
                    true);
            runCommandLine(arguments);
        }

    }

    @Test
    public void testCustomVariantClassificationOrder() {

        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                largeFileTestDir + "funcotator" + File.separator + "custom_vc_order_files" + File.separator + "custom_vc_input_test.vcf",
                outputFile,
                b37Reference,
                MEDIUM_DATASOURCES_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.add(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        // It's best to make this run on BEST_EFFECT with this test
        arguments.add(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.BEST_EFFECT);

        // Add the new ordering of the variant classifications:
        final GATKPath custom_vc_path = new GATKPath(largeFileTestDir + "funcotator" + File.separator
                + "custom_vc_order_files" + File.separator + "custom_vc_order_for_int_test.tsv");
        arguments.add(FuncotatorArgumentDefinitions.CUSTOM_VARIANT_CLASS_ORDER_FILE, custom_vc_path.toPath().toUri().toString());

        // Run Funcotator:
        runCommandLine(arguments);

        // Validate results:
        final String expectedFilePath = largeFileTestDir + "funcotator" + File.separator + "custom_vc_order_files" + File.separator + "custom_vc_expected_out.vcf";

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        final Pair<VCFHeader, List<VariantContext>> expectedVcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(expectedFilePath);
        final List<VariantContext> expectedVariantContexts = expectedVcfInfo.getRight();

        Assert.assertEquals(variantContexts.size(), expectedVariantContexts.size());

        try {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, new File(expectedFilePath), "#");
        }
        catch ( final IOException ex ) {
            throw new GATKException("Error opening expected file: " + expectedFilePath, ex);
        }

        // Reset severity:
        GencodeFuncotation.VariantClassification.resetSeveritiesToDefault();
    }

}


