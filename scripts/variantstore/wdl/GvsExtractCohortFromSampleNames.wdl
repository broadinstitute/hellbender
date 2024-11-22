version 1.0

import "GvsPrepareRangesCallset.wdl" as GvsPrepareCallset
import "GvsExtractCallset.wdl" as GvsExtractCallset
import "GvsUtils.wdl" as Utils

# Workflow used by AoU to extract variants for a given cohort of sample_names

workflow GvsExtractCohortFromSampleNames {

  input {
    # cohort_sample_names_array will take precedence over cohort_sample_names if both are set
    Array[String]? cohort_sample_names_array
    File? cohort_sample_names
    Boolean is_wgs = true

    String gvs_project
    String gvs_dataset
    String call_set_identifier
    String cohort_table_prefix = call_set_identifier
    String query_project = gvs_project

    # not using the defaults in GvsPrepareCallset because we might be using pre created datasets defined by the caller
    String destination_dataset_name = gvs_dataset
    String destination_project_id = gvs_project
    String? fq_gvs_extraction_temp_tables_dataset
    String extraction_uuid = call_set_identifier
    String filter_set_name
    String output_file_base_name = call_set_identifier
    Boolean? control_samples

    String? output_gcs_dir
    # set to "NONE" if all the reference data was loaded into GVS in GvsImportGenomes
    String drop_state = "NONE"
    Boolean bgzip_output_vcfs = false
    String ploidy_table_name = "sample_chromosome_ploidy"
    Boolean collect_variant_calling_metrics = false

    File? interval_list
    Int? extract_preemptible_override
    Int? extract_maxretries_override
    Int? extract_scatter_count_override
    Int? extract_memory_override
    Int? extract_disk_override
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override

    File? target_interval_list

    String? git_branch_or_tag
    String? gatk_docker
    File? gatk_override
    String? cloud_sdk_docker
    String? variants_docker
  }

  Boolean write_cost_to_db = if ((gvs_project != destination_project_id) || (gvs_project != query_project)) then false else true

  if (!defined(git_branch_or_tag) || !defined(gatk_docker)  || !defined(cloud_sdk_docker) || !defined(variants_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_git_hash = select_first([git_branch_or_tag, GetToolVersions.git_hash])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])


  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      project_id = query_project,
      fq_table = "~{gvs_project}.~{gvs_dataset}.sample_info",
      cloud_sdk_docker = effective_cloud_sdk_docker
  }

  call Utils.GetNumSamplesLoaded {
    input:
      fq_sample_table = "~{gvs_project}.~{gvs_dataset}.sample_info",
      project_id = query_project,
      sample_table_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp,
      cloud_sdk_docker = effective_cloud_sdk_docker
  }

  # writing the array to a file has to be done in a task
  # https://support.terra.bio/hc/en-us/community/posts/360071465631-write-lines-write-map-write-tsv-write-json-fail-when-run-in-a-workflow-rather-than-in-a-task
  if (defined(cohort_sample_names_array)) {
    call write_array_task {
      input:
        input_array = select_first([cohort_sample_names_array]),
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }
  }

  File cohort_sample_names_file = select_first([write_array_task.output_file, cohort_sample_names])

  # scatter for WGS and exome samples based on past successful runs and NOT optimized
  Int effective_scatter_count = if defined(extract_scatter_count_override) then select_first([extract_scatter_count_override])
                                else if is_wgs then
                                     if GetNumSamplesLoaded.num_samples < 5000 then 1 # This results in 1 VCF per chromosome.
                                     else if GetNumSamplesLoaded.num_samples < 20000 then 2000 # Stroke Anderson
                                          else if GetNumSamplesLoaded.num_samples < 50000 then 10000
                                               else 20000
                                     else
                                     if GetNumSamplesLoaded.num_samples < 5000 then 1 # This results in 1 VCF per chromosome.
                                     else if GetNumSamplesLoaded.num_samples < 20000 then 1000
                                          else if GetNumSamplesLoaded.num_samples < 50000 then 2500
                                               else 7500

  # allow an interval list to be passed in, but default it to our standard one if no args are here
  File effective_interval_list = select_first([interval_list, "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"])


  call GvsPrepareCallset.GvsPrepareCallset {
    input:
      call_set_identifier = call_set_identifier,
      extract_table_prefix = cohort_table_prefix,
      sample_names_to_extract = cohort_sample_names_file,
      project_id = gvs_project,
      query_labels = ["extraction_uuid=~{extraction_uuid}"],
      query_project = query_project,
      dataset_name = gvs_dataset, # unused if fq_* args are given
      destination_project = destination_project_id,
      destination_dataset = destination_dataset_name,
      fq_temp_table_dataset = fq_gvs_extraction_temp_tables_dataset,
      write_cost_to_db = write_cost_to_db,
      enable_extract_table_ttl = true,
      interval_list = effective_interval_list,
      control_samples = control_samples,
      cloud_sdk_docker = effective_cloud_sdk_docker,
      git_hash = effective_git_hash,
      variants_docker = effective_variants_docker,
  }

  call GvsExtractCallset.GvsExtractCallset {
    input:
      go = GvsPrepareCallset.done,
      project_id = gvs_project,
      query_project = query_project,
      dataset_name = gvs_dataset,
      call_set_identifier = call_set_identifier,
      cohort_project_id = destination_project_id,
      cohort_dataset_name = destination_dataset_name,
      extract_table_prefix = cohort_table_prefix,

      scatter_count = effective_scatter_count,
      filter_set_name = filter_set_name,
      output_file_base_name = output_file_base_name,
      output_gcs_dir = output_gcs_dir,

      drop_state = drop_state,
      bgzip_output_vcfs = bgzip_output_vcfs,
      collect_variant_calling_metrics = collect_variant_calling_metrics,
      ploidy_table_name = ploidy_table_name,
      extract_preemptible_override = extract_preemptible_override,
      extract_maxretries_override = extract_maxretries_override,
      split_intervals_disk_size_override = split_intervals_disk_size_override,
      split_intervals_mem_override = split_intervals_mem_override,
      extract_memory_override_gib = extract_memory_override,
      disk_override = extract_disk_override,
      interval_list = effective_interval_list,
      control_samples = control_samples,

      cloud_sdk_docker = effective_cloud_sdk_docker,
      gatk_override = gatk_override,
      gatk_docker = effective_gatk_docker,
      git_hash = effective_git_hash,
      variants_docker = effective_variants_docker,
      write_cost_to_db = write_cost_to_db,
      target_interval_list = target_interval_list,
  }

  output {
    Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
    Array[File] output_vcfs = GvsExtractCallset.output_vcfs
    Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
    String recorded_git_hash = effective_git_hash
  }

}

task write_array_task {
  input {
    Array[String] input_array
    String cloud_sdk_docker
  }

  command <<<
  >>>

  output {
    File output_file = write_lines(input_array)
  }

  runtime {
    docker: cloud_sdk_docker
  }
}

