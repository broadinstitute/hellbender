version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsExtractCallset {
  input {
    Boolean go = true
    String dataset_name
    String project_id
    String call_set_identifier

    String cohort_project_id = project_id
    String cohort_dataset_name = dataset_name
    Boolean do_not_filter_override = false
    Boolean control_samples = false
    String extract_table_prefix
    String filter_set_name
    String query_project = project_id
    # This is optional now since the workflow will choose an appropriate value below if this is unspecified.
    Int? scatter_count
    Int? extract_memory_override_gib
    # The amount of memory on the extract VM *not* allocated to the GATK process.
    Int extract_overhead_memory_override_gib = 3
    Int? disk_override
    Boolean bgzip_output_vcfs = false
    Boolean zero_pad_output_vcf_filenames = true
    Boolean collect_variant_calling_metrics = false

    # set to "NONE" if all the reference data was loaded into GVS in GvsImportGenomes
    String drop_state = "NONE"

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File interval_weights_bed = "gs://gvs_quickstart_storage/weights/gvs_full_vet_weights_1kb_padded_orig.bed"

    File? target_interval_list

    String? variants_docker
    String? cloud_sdk_docker
    String? gatk_docker
    String? git_branch_or_tag
    String? git_hash
    File? gatk_override

    String output_file_base_name = filter_set_name

    String? ploidy_table_name
    Int? extract_maxretries_override
    Int? extract_preemptible_override
    String? output_gcs_dir
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override
    Float x_bed_weight_scaling = 4
    Float y_bed_weight_scaling = 4
    Boolean is_wgs = true
    Boolean convert_filtered_genotypes_to_nocalls = false
    Boolean write_cost_to_db = true
    Int maximum_alternate_alleles = 1000
  }

  File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  File dbsnp_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  File dbsnp_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

  String fq_gvs_dataset = "~{project_id}.~{dataset_name}"
  String fq_cohort_dataset = "~{cohort_project_id}.~{cohort_dataset_name}"

  String full_extract_prefix = if (control_samples) then "~{extract_table_prefix}_controls" else extract_table_prefix

  String fq_filter_set_info_table = "~{fq_gvs_dataset}.filter_set_info"
  String fq_filter_set_site_table = "~{fq_gvs_dataset}.filter_set_sites"
  String fq_filter_set_tranches_table = "~{fq_gvs_dataset}.filter_set_tranches"
  String fq_sample_table = "~{fq_gvs_dataset}.sample_info"
  String fq_cohort_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__DATA"
  String fq_ranges_cohort_ref_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__REF_DATA"
  String fq_ranges_cohort_vet_extract_table_name = "~{full_extract_prefix}__VET_DATA"
  String fq_ranges_cohort_vet_extract_table = "~{fq_cohort_dataset}.~{fq_ranges_cohort_vet_extract_table_name}"

  # make the fully qualified version of the ploidy table if present, otherwise leave it undefined
  String? fq_ploidy_mapping_table = if (defined(ploidy_table_name)) then "~{fq_gvs_dataset}.~{ploidy_table_name}" else ploidy_table_name

  String fq_samples_to_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__SAMPLES"
  Array[String] tables_patterns_for_datetime_check = ["~{full_extract_prefix}__%"]

  Boolean emit_pls = false
  Boolean emit_ads = true

  String intervals_file_extension = if (zero_pad_output_vcf_filenames) then '-~{output_file_base_name}.interval_list' else '-scattered.interval_list'
  String vcf_extension = if (bgzip_output_vcfs) then '.vcf.bgz' else '.vcf.gz'

  if (!defined(git_hash) || !defined(gatk_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  call Utils.ScaleXYBedValues {
    input:
      interval_weights_bed = interval_weights_bed,
      x_bed_weight_scaling = x_bed_weight_scaling,
      y_bed_weight_scaling = y_bed_weight_scaling,
      variants_docker = effective_variants_docker,
  }

  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      project_id = query_project,
      fq_table = fq_sample_table,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call Utils.GetNumSamplesLoaded {
    input:
      fq_sample_table = fq_sample_table,
      project_id = query_project,
      control_samples = control_samples,
      sample_table_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  # scatter for WGS and exome samples based on past successful runs and NOT optimized
  Int effective_scatter_count = if defined(scatter_count) then select_first([scatter_count])
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

  Int effective_split_intervals_disk_size_override = select_first([split_intervals_disk_size_override,
                                                                  if GetNumSamplesLoaded.num_samples < 100 then 50 # Quickstart
                                                                  else 500])

  Int effective_extract_memory_gib = if defined(extract_memory_override_gib) then select_first([extract_memory_override_gib])
                                     else if effective_scatter_count <= 100 then 37 + extract_overhead_memory_override_gib
                                          else if effective_scatter_count <= 500 then 17 + extract_overhead_memory_override_gib
                                               else 9 + extract_overhead_memory_override_gib

  # WDL 1.0 trick to set a variable ('none') to be undefined.
  if (false) {
    File? none = ""
  }

  call Utils.SplitIntervals {
    input:
      intervals = interval_list,
      ref_fasta = reference,
      ref_fai = reference_index,
      ref_dict = reference_dict,
      interval_weights_bed = ScaleXYBedValues.xy_scaled_bed,
      intervals_file_extension = intervals_file_extension,
      scatter_count = effective_scatter_count,
      output_gcs_dir = output_gcs_dir,
      split_intervals_disk_size_override = effective_split_intervals_disk_size_override,
      split_intervals_mem_override = split_intervals_mem_override,
      gatk_docker = effective_gatk_docker,
      gatk_override = gatk_override,
  }

  call Utils.GetBQTableLastModifiedDatetime as FilterSetInfoTimestamp {
    input:
      project_id = project_id,
      fq_table = "~{fq_filter_set_info_table}",
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  if ( !do_not_filter_override ) {
    call Utils.ValidateFilterSetName {
      input:
        project_id = query_project,
        fq_filter_set_info_table = "~{fq_filter_set_info_table}",
        filter_set_name = filter_set_name,
        filter_set_info_timestamp = FilterSetInfoTimestamp.last_modified_timestamp,
        cloud_sdk_docker = effective_cloud_sdk_docker,
      }

    call Utils.IsVETS {
      input:
        project_id = query_project,
        fq_filter_set_info_table = "~{fq_filter_set_info_table}",
        filter_set_name = filter_set_name,
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }
  }

  # If we're not using the VQSR filters, set it to VETS (really shouldn't matter one way or the other)
  # Otherwise use the auto-derived flag.
  Boolean use_VETS = select_first([IsVETS.is_vets, true])

  call Utils.GetBQTablesMaxLastModifiedTimestamp {
    input:
      query_project = query_project,
      data_project = project_id,
      dataset_name = dataset_name,
      table_patterns = tables_patterns_for_datetime_check,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call Utils.GetExtractVetTableVersion {
    input:
      query_project = query_project,
      data_project = project_id,
      dataset_name = dataset_name,
      table_name = fq_ranges_cohort_vet_extract_table_name,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  scatter(i in range(length(SplitIntervals.interval_files))) {
    String interval_filename = basename(SplitIntervals.interval_files[i])
    String vcf_filename = if (zero_pad_output_vcf_filenames) then sub(interval_filename, ".interval_list", "") else "~{output_file_base_name}_${i}"

    call ExtractTask {
      input:
        go                                    = select_first([ValidateFilterSetName.done, true]),
        dataset_name                          = dataset_name,
        call_set_identifier                   = call_set_identifier,
        use_VETS                              = use_VETS,
        gatk_docker                           = effective_gatk_docker,
        gatk_override                         = gatk_override,
        reference                             = reference,
        reference_index                       = reference_index,
        reference_dict                        = reference_dict,
        fq_samples_to_extract_table           = fq_samples_to_extract_table,
        interval_index                        = i,
        intervals                             = SplitIntervals.interval_files[i],
        fq_cohort_extract_table               = fq_cohort_extract_table,
        fq_ranges_cohort_ref_extract_table    = fq_ranges_cohort_ref_extract_table,
        fq_ranges_cohort_vet_extract_table    = fq_ranges_cohort_vet_extract_table,
        vet_extract_table_version             = GetExtractVetTableVersion.version,
        read_project_id                       = query_project,
        do_not_filter_override                = do_not_filter_override,
        fq_filter_set_info_table              = fq_filter_set_info_table,
        fq_filter_set_site_table              = fq_filter_set_site_table,
        fq_ploidy_mapping_table               = fq_ploidy_mapping_table,
        fq_filter_set_tranches_table          = if (use_VETS) then none else fq_filter_set_tranches_table,
        filter_set_name                       = filter_set_name,
        drop_state                            = drop_state,
        output_file                           = vcf_filename + vcf_extension,
        max_last_modified_timestamp           = GetBQTablesMaxLastModifiedTimestamp.max_last_modified_timestamp,
        extract_preemptible_override          = extract_preemptible_override,
        extract_maxretries_override           = extract_maxretries_override,
        disk_override                         = disk_override,
        memory_gib                            = effective_extract_memory_gib,
        overhead_memory_gib                   = extract_overhead_memory_override_gib,
        emit_pls                              = emit_pls,
        emit_ads                              = emit_ads,
        convert_filtered_genotypes_to_nocalls = convert_filtered_genotypes_to_nocalls,
        write_cost_to_db                      = write_cost_to_db,
        maximum_alternate_alleles             = maximum_alternate_alleles,
        target_interval_list                  = target_interval_list,
    }

    if (collect_variant_calling_metrics) {
      call CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ExtractTask.output_vcf,
          input_vcf_index = ExtractTask.output_vcf_index,
          metrics_filename_prefix = call_set_identifier + "." + i,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = SplitIntervals.interval_files[i],
          ref_dict = reference_dict,
          gatk_docker = effective_gatk_docker
      }
    }
  }

  if (collect_variant_calling_metrics) {
    call GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = call_set_identifier,
        output_gcs_dir = output_gcs_dir,
        gatk_docker = effective_gatk_docker
    }
  }

  call SumBytes {
    input:
      file_sizes_bytes = flatten([ExtractTask.output_vcf_bytes, ExtractTask.output_vcf_index_bytes]),
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call CreateManifestAndOptionallyCopyOutputs {
    input:
      interval_indices = ExtractTask.interval_number,
      output_vcfs = ExtractTask.output_vcf,
      output_vcf_indices = ExtractTask.output_vcf_index,
      output_vcf_bytes = ExtractTask.output_vcf_bytes,
      output_vcf_index_bytes = ExtractTask.output_vcf_index_bytes,
      output_gcs_dir = output_gcs_dir,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call Utils.GetBQTableLastModifiedDatetime {
    input:
      project_id = query_project,
      fq_table = fq_samples_to_extract_table,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GenerateSampleListFile {
    input:
      fq_samples_to_extract_table = fq_samples_to_extract_table,
      samples_to_extract_table_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp,
      output_gcs_dir = output_gcs_dir,
      query_project = query_project,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  output {
    Array[File] output_vcfs = ExtractTask.output_vcf
    Array[File] output_vcf_indexes = ExtractTask.output_vcf_index
    Array[File] output_vcf_interval_files = SplitIntervals.interval_files
    Float total_vcfs_size_mb = SumBytes.total_mb
    File manifest = CreateManifestAndOptionallyCopyOutputs.manifest
    File sample_name_list = GenerateSampleListFile.sample_name_list
    File? summary_metrics_file = GatherVariantCallingMetrics.summary_metrics_file
    File? detail_metrics_file = GatherVariantCallingMetrics.detail_metrics_file
    String recorded_git_hash = effective_git_hash
    Boolean done = true
  }
}

task ExtractTask {
  input {
    Boolean go

    String dataset_name
    String call_set_identifier

    Boolean use_VETS

    File reference
    File reference_index
    File reference_dict

    String fq_samples_to_extract_table

    Int interval_index
    File intervals
    String drop_state

    String fq_cohort_extract_table
    String fq_ranges_cohort_ref_extract_table
    String fq_ranges_cohort_vet_extract_table
    String? fq_ploidy_mapping_table
    String? vet_extract_table_version
    String read_project_id
    String output_file

    String cost_observability_tablename = "cost_observability"

    Boolean emit_pls
    Boolean emit_ads
    Boolean convert_filtered_genotypes_to_nocalls = false

    Boolean do_not_filter_override
    String fq_filter_set_info_table
    String fq_filter_set_site_table
    String? fq_filter_set_tranches_table
    String? filter_set_name
    Boolean write_cost_to_db

    # Runtime Options:
    String gatk_docker
    File? gatk_override
    Int? extract_preemptible_override
    Int? extract_maxretries_override
    Int? disk_override
    Int memory_gib
    Int overhead_memory_gib

    Int? local_sort_max_records_in_ram = 10000000
    Int? maximum_alternate_alleles

    File? target_interval_list

    # for call-caching -- check if DB tables haven't been updated since the last run
    String max_last_modified_timestamp
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  String intervals_name = basename(intervals)
  String cost_observability_line = if (write_cost_to_db == true) then "--cost-observability-tablename ~{cost_observability_tablename}" else ""

  String inferred_reference_state = if (drop_state == "NONE") then "ZERO" else drop_state

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    bash ~{monitoring_script} > monitoring.log &

    export GATK_LOCAL_JAR="~{default="/root/gatk.jar" gatk_override}"

    if [ ~{do_not_filter_override} = true ]; then
      FILTERING_ARGS=''
    elif [ ~{use_VETS} = false ]; then
      FILTERING_ARGS='--filter-set-info-table ~{fq_filter_set_info_table}
        --filter-set-site-table ~{fq_filter_set_site_table}
        --tranches-table ~{fq_filter_set_tranches_table}
        --filter-set-name ~{filter_set_name}'
    else
      FILTERING_ARGS='--filter-set-info-table ~{fq_filter_set_info_table}
        --filter-set-site-table ~{fq_filter_set_site_table}
        --filter-set-name ~{filter_set_name}'
    fi

    # This tool may get invoked with "Retry with more memory" with a different amount of memory than specified in
    # the input `memory_gib`. If so, use the memory-related environment variables rather than the `memory_gib` input.
    # But also be prepared if those memory-related variables are not set and fall back to using `memory_gib`.
    # https://support.terra.bio/hc/en-us/articles/4403215299355-Out-of-Memory-Retry
    if [[ -z "${MEM_UNIT:-}" ]]
    then
        memory_mb=$(python3 -c "from math import floor; print(floor((~{memory_gib} - ~{overhead_memory_gib}) * 1000))")
    elif [[ ${MEM_UNIT} == "GB" ]]
    then
        memory_mb=$(python3 -c "from math import floor; print(floor((${MEM_SIZE} - ~{overhead_memory_gib}) * 1000))")
    else
        echo "Unexpected memory unit: ${MEM_UNIT}" 1>&2
        exit 1
    fi

    gatk --java-options "-Xmx${memory_mb}m" \
      ExtractCohortToVcf \
        --vet-ranges-extract-fq-table ~{fq_ranges_cohort_vet_extract_table} \
        ~{"--vet-ranges-extract-table-version " + vet_extract_table_version} \
        --ref-ranges-extract-fq-table ~{fq_ranges_cohort_ref_extract_table} \
        --ref-version 38 \
        -R ~{reference} \
        -O ~{output_file} \
        --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
        --sample-table ~{fq_samples_to_extract_table} \
        ~{"--inferred-reference-state " + inferred_reference_state} \
        -L ~{intervals} \
        --project-id ~{read_project_id} \
        ~{true='--emit-pls' false='' emit_pls} \
        ~{true='--emit-ads' false='' emit_ads} \
        ~{true='' false='--use-vqsr-scoring' use_VETS} \
        ~{true='--convert-filtered-genotypes-to-no-calls' false='' convert_filtered_genotypes_to_nocalls} \
        ~{'--maximum-alternate-alleles ' + maximum_alternate_alleles} \
        ${FILTERING_ARGS} \
        ~{"--sample-ploidy-table " + fq_ploidy_mapping_table} \
        --dataset-id ~{dataset_name} \
        --call-set-identifier ~{call_set_identifier} \
        --wdl-step GvsExtractCallset \
        --wdl-call ExtractTask \
        --shard-identifier ~{intervals_name} \
        ~{cost_observability_line}

    if [[ -n "~{target_interval_list}" ]]
    then
      pre_off_target_vcf="pre_off_target_~{output_file}"
      mv ~{output_file} ${pre_off_target_vcf}
      mv ~{output_file}.tbi "${pre_off_target_vcf}.tbi"

      gatk --java-options "-Xmx${memory_mb}m" \
        IndexFeatureFile \
          -I ~{target_interval_list}

      gatk --java-options "-Xmx${memory_mb}m" \
        VariantFiltration \
          ~{"--filter-not-in-mask --mask-name OUTSIDE_OF_TARGETS --mask-description 'Outside of sequencing target intervals' --mask " + target_interval_list} \
          -O ~{output_file} \
          -V ${pre_off_target_vcf}
    fi

    du -b ~{output_file} | cut -f1 > vcf_bytes.txt
    du -b ~{output_file}.tbi | cut -f1 > vcf_index_bytes.txt

  >>>
  runtime {
    docker: gatk_docker
    memory: memory_gib + " GB"
    disks: "local-disk " + select_first([disk_override, 150]) + " HDD"
    bootDiskSizeGb: 15
    preemptible: select_first([extract_preemptible_override, "2"])
    maxRetries: select_first([extract_maxretries_override, "3"])
    cpu: 2
    noAddress: true
  }

  # files sizes are floats instead of ints because they can be larger
  output {
    Int interval_number = interval_index
    File output_vcf = "~{output_file}"
    Float output_vcf_bytes = read_float("vcf_bytes.txt")
    File output_vcf_index = "~{output_file}.tbi"
    Float output_vcf_index_bytes = read_float("vcf_index_bytes.txt")
    File monitoring_log = "monitoring.log"
  }
}

task CreateManifestAndOptionallyCopyOutputs {
  input {
    Array[Int] interval_indices
    Array[File] output_vcfs
    Array[File] output_vcf_indices
    Array[Float] output_vcf_bytes
    Array[Float] output_vcf_index_bytes
    String? output_gcs_dir
    String cloud_sdk_docker
  }
  parameter_meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    output_vcfs: {
      localization_optional: true
    }
    output_vcf_indices: {
      localization_optional: true
    }
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    declare -a interval_indices=(~{sep=' ' interval_indices})
    declare -a output_vcfs=(~{sep=' ' output_vcfs})
    declare -a output_vcf_indices=(~{sep=' ' output_vcf_indices})
    declare -a output_vcf_bytes=(~{sep=' ' output_vcf_bytes})
    declare -a output_vcf_index_bytes=(~{sep=' ' output_vcf_index_bytes})

    # (Possibly) create a manifest of VCFs and indexes to bulk copy with `gcloud storage cp`.
    echo -n > vcf_manifest.txt

    echo -n >> manifest_lines.txt
    for (( i=0; i<${#interval_indices[@]}; ++i));
      do
        echo "Interval " + $i

        OUTPUT_VCF=${output_vcfs[$i]}
        LOCAL_VCF=$(basename $OUTPUT_VCF)
        OUTPUT_VCF_INDEX=${output_vcf_indices[$i]}
        LOCAL_VCF_INDEX=$(basename $OUTPUT_VCF_INDEX)

        if [ -n "${OUTPUT_GCS_DIR}" ]; then
          echo $OUTPUT_VCF >> vcf_manifest.txt
          echo $OUTPUT_VCF_INDEX >> vcf_manifest.txt
          OUTPUT_FILE_DEST="${OUTPUT_GCS_DIR}/$LOCAL_VCF"
          OUTPUT_FILE_INDEX_DEST="${OUTPUT_GCS_DIR}/$LOCAL_VCF_INDEX"
        else
          OUTPUT_FILE_DEST=$LOCAL_VCF
          OUTPUT_FILE_INDEX_DEST=$LOCAL_VCF_INDEX
        fi

        echo ${interval_indices[$i]},${OUTPUT_FILE_DEST},${output_vcf_bytes[$i]},${OUTPUT_FILE_INDEX_DEST},${output_vcf_index_bytes[$i]} >> manifest_lines.txt

      done;

    echo "vcf_file_location, vcf_file_bytes, vcf_index_location, vcf_index_bytes" >> manifest.txt
    sort -n manifest_lines.txt | cut -d',' -f 2- >> manifest.txt

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      # Copy VCFs, indexes and the manifest to the output directory.
      echo manifest.txt >> vcf_manifest.txt
      cat vcf_manifest.txt | xargs -IFILE gsutil cp FILE ${OUTPUT_GCS_DIR}
    fi
  >>>
  output {
    File manifest_lines = "manifest_lines.txt"
    File manifest = "manifest.txt"
  }

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 500 HDD"
    preemptible: 3
    cpu: 1
  }
}

task SumBytes {
  input {
    Array[Float] file_sizes_bytes
    String cloud_sdk_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "~{sep=" " file_sizes_bytes}" | tr " " "\n" | python3 -c "
    import sys;
    total_bytes = sum(float(i.strip()) for i in sys.stdin);
    total_mb = total_bytes/10**6;
    print(total_mb);"
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 500 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Float total_mb = read_float(stdout())
  }
}

task GenerateSampleListFile {
  input {
    String fq_samples_to_extract_table
    String samples_to_extract_table_timestamp
    String query_project

    String? output_gcs_dir
    String cloud_sdk_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:extract_callset"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    # bq query --max_rows check: max rows set to at least the number of samples
    bq --apilog=false --project_id=~{query_project} --format=csv query --max_rows 1000000000 --use_legacy_sql=false ~{bq_labels} \
      'SELECT sample_name FROM `~{fq_samples_to_extract_table}`' | sed 1d > sample-name-list.txt

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil cp sample-name-list.txt ${OUTPUT_GCS_DIR}/
    fi
  >>>
  output {
    File sample_name_list = "sample-name-list.txt"
  }

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 500 HDD"
    preemptible: 3
    cpu: 1
  }
}

task CollectVariantCallingMetrics {
  input {
    File input_vcf
    File input_vcf_index
    File dbsnp_vcf
    File dbsnp_vcf_index
    File interval_list
    File ref_dict
    String metrics_filename_prefix

    Int memory_mb = 7500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 200
    String gatk_docker
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    bash ~{monitoring_script} > monitoring.log &

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      CollectVariantCallingMetrics \
      --INPUT ~{input_vcf} \
      --DBSNP ~{dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ~{ref_dict} \
      --OUTPUT ~{metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ~{interval_list}
  >>>

  output {
    File summary_metrics_file = "~{metrics_filename_prefix}.variant_calling_summary_metrics"
    File detail_metrics_file = "~{metrics_filename_prefix}.variant_calling_detail_metrics"
    File monitoring_log = "monitoring.log"
  }

  runtime {
    docker: gatk_docker
    cpu: 2
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 2
    noAddress: true
  }
}

task GatherVariantCallingMetrics {

  input {
    Array[File] input_details
    Array[File] input_summaries
    String output_prefix
    String? output_gcs_dir

    Int memory_mb = 3000
    Int disk_size_gb = 200
    String gatk_docker
  }

  parameter_meta {
    input_details: {
      localization_optional: true
    }
    input_summaries: {
      localization_optional: true
    }
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    bash ~{monitoring_script} > monitoring.log &

    input_details_fofn=~{write_lines(input_details)}
    input_summaries_fofn=~{write_lines(input_summaries)}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat $input_details_fofn | gsutil -m cp -L cp.log -c -I metrics/; do
      sleep 1
      ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
      echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat $input_summaries_fofn | gsutil -m cp -L cp.log -c -I metrics/; do
      sleep 1
      ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
      echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=$(cat $input_details_fofn | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("--INPUT metrics/%s ", $1)}')

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      AccumulateVariantCallingMetrics \
      $INPUT \
      --OUTPUT ~{output_prefix}

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil cp ~{output_prefix}.variant_calling_summary_metrics ${OUTPUT_GCS_DIR}/
      gsutil cp ~{output_prefix}.variant_calling_detail_metrics ${OUTPUT_GCS_DIR}/
    fi
  >>>

  runtime {
    docker: gatk_docker
    cpu: 1
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
    noAddress: true
  }

  output {
    File summary_metrics_file = "~{output_prefix}.variant_calling_summary_metrics"
    File detail_metrics_file = "~{output_prefix}.variant_calling_detail_metrics"
    File monitoring_log = "monitoring.log"
  }
}
