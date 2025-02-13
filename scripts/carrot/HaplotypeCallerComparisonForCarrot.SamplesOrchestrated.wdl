version 1.0

 import "https://raw.githubusercontent.com/broadinstitute/gatk/je_carrotHCTestsScripts/scripts/carrot/HaplotypeCallerComparisonForCarrot.monitoring.wdl" as RunSampleHeadToHead

 workflow VariantCallingCarrotOrchestrated {

  input {
    ## Samples:
    File CHM_input_bam
    File CHM_input_bam_index
    String CHM_base_file_name
    File CHM_calling_interval_list
    Float CHM_contamination

    File NIST_input_bam
    File NIST_input_bam_index
    String NIST_base_file_name
    File NIST_calling_interval_list
    Float NIST_contamination

    File exome1_input_bam
    File exome1_input_bam_index
    String exome1_base_file_name
    File exome1_calling_interval_list


    Boolean run_dragen_mode_variant_calling = false
    Boolean use_spanning_event_genotyping = true
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Float? contamination = 0
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_str
    Int agg_preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gatk_control_docker = "broadinstitute/gatk-nightly:latest"
    Boolean make_gvcf = true
    Boolean make_bamout = false
    Boolean use_gatk3_haplotype_caller = false
    Boolean skip_reblocking = false
    Boolean use_dragen_hard_filtering = false
    File monitoring_script = "gs://emeryj-testing/cromwell_monitoring_script.sh"
  }

  call RunSampleHeadToHead.VariantCallingCarrot as CHMSampleHeadToHead {
       input:
            input_bam            = CHM_input_bam,
            input_bam_index      = CHM_input_bam_index,
            base_file_name       = CHM_base_file_name,
            final_vcf_base_name  = CHM_base_file_name,
            calling_interval_list = CHM_calling_interval_list,

            haplotype_scatter_count = haplotype_scatter_count,
            break_bands_at_multiples_of = break_bands_at_multiples_of,
            contamination = CHM_contamination,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            agg_preemptible_tries = 0,
            gatk_docker = gatk_docker,
            gatk_control_docker = gatk_control_docker,
            monitoring_script = monitoring_script
  }

    call RunSampleHeadToHead.VariantCallingCarrot as NISTSampleHeadToHead {
         input:
              input_bam            = NIST_input_bam,
              input_bam_index      = NIST_input_bam_index,
              base_file_name       = NIST_base_file_name,
              final_vcf_base_name  = NIST_base_file_name,
              calling_interval_list = NIST_calling_interval_list,

              haplotype_scatter_count = haplotype_scatter_count,
              break_bands_at_multiples_of = break_bands_at_multiples_of,
              contamination = NIST_contamination,
              ref_fasta = ref_fasta,
              ref_fasta_index = ref_fasta_index,
              ref_dict = ref_dict,
              agg_preemptible_tries = 0,
              gatk_docker = gatk_docker,
              gatk_control_docker = gatk_control_docker,
              monitoring_script = monitoring_script
    }

    call RunSampleHeadToHead.VariantCallingCarrot as EXOME1SampleHeadToHead {
         input:
              input_bam            = exome1_input_bam,
              input_bam_index      = exome1_input_bam_index,
              base_file_name       = exome1_base_file_name,
              final_vcf_base_name  = exome1_base_file_name,
              calling_interval_list = exome1_calling_interval_list,

              haplotype_scatter_count = haplotype_scatter_count,
              break_bands_at_multiples_of = break_bands_at_multiples_of,
              contamination = 0,
              ref_fasta = ref_fasta,
              ref_fasta_index = ref_fasta_index,
              ref_dict = ref_dict,
              agg_preemptible_tries = 0,
              gatk_docker = gatk_docker,
              gatk_control_docker = gatk_control_docker,
              monitoring_script = monitoring_script
    }


    output {
            File CHM_output_vcf = CHMSampleHeadToHead.output_vcf
            File CHM_output_vcf_index = CHMSampleHeadToHead.output_vcf_index
            Array[File] CHM_output_runtimes = CHMSampleHeadToHead.output_runtimes
            File CHM_representative_benchmarking = CHMSampleHeadToHead.representative_benchmarking
            File CHM_control_vcf = CHMSampleHeadToHead.control_vcf
            File CHM_control_vcf_index = CHMSampleHeadToHead.control_vcf_index
            Array[File] CHM_control_output_runtimes = CHMSampleHeadToHead.control_output_runtimes
            File CHM_control_representative_benchmarking = CHMSampleHeadToHead.control_representative_benchmarking

            File NIST_output_vcf = NISTSampleHeadToHead.output_vcf
            File NIST_output_vcf_index = NISTSampleHeadToHead.output_vcf_index
            Array[File] NIST_output_runtimes = NISTSampleHeadToHead.output_runtimes
            File NIST_representative_benchmarking = NISTSampleHeadToHead.representative_benchmarking
            File NIST_control_vcf = NISTSampleHeadToHead.control_vcf
            File NIST_control_vcf_index = NISTSampleHeadToHead.control_vcf_index
            Array[File] NIST_control_output_runtimes = NISTSampleHeadToHead.control_output_runtimes
            File NIST_control_representative_benchmarking = NISTSampleHeadToHead.control_representative_benchmarking

            File EXOME1_output_vcf = EXOME1SampleHeadToHead.output_vcf
            File EXOME1_output_vcf_index = EXOME1SampleHeadToHead.output_vcf_index
            Array[File] EXOME1_output_runtimes = EXOME1SampleHeadToHead.output_runtimes
            File EXOME1_representative_benchmarking = EXOME1SampleHeadToHead.representative_benchmarking
            File EXOME1_control_vcf = EXOME1SampleHeadToHead.control_vcf
            File EXOME1_control_vcf_index = EXOME1SampleHeadToHead.control_vcf_index
            Array[File] EXOME1_control_output_runtimes = EXOME1SampleHeadToHead.control_output_runtimes
            File EXOME1_control_representative_benchmarking = EXOME1SampleHeadToHead.control_representative_benchmarking
  }
}
 
