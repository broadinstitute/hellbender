version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsImportGenomes {

  input {
    Boolean go = true
    String? git_branch_or_tag
    String? git_hash
    String dataset_name
    String project_id

    Int num_samples
    File external_sample_names
    File input_vcfs
    File input_vcf_indexes

    Boolean skip_loading_vqsr_fields = false
    Boolean use_compressed_references = false

    # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
    String drop_state = "NONE"
    # beta customers will almost always have a naive GCP account, and as such will not be able to cross over their quotas
    # without Google shutting import down by throwing them API errors.  For them, we limit our scattering.
    Boolean is_rate_limited_beta_customer = false
    # This was determined to be the point at which we come close to but don't cross over the "AppendRows throughput per
    # project for small regions per minute per region" default quota of ~19G.  Uses up to ~90% of the quota at peaks
    # without going over
    Int beta_customer_max_scatter = 200
    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    Int? load_data_scatter_width
    Int? load_data_preemptible_override
    Int? load_data_maxretries_override
    # At least one of these "load" inputs must be true
    Boolean load_vet_and_ref_ranges = true
    Boolean load_vcf_headers = false
    String? basic_docker
    String? cloud_sdk_docker
    String? variants_docker
    String? gatk_docker
    File? load_data_gatk_override
    String? billing_project_id
    Boolean is_wgs = true
  }

  Int max_auto_batch_size = if is_wgs then 25000 else 100000
  String genome_type = if is_wgs then "WGS" else "exome"

  # Broad users enjoy higher quotas and can scatter more widely than beta users before BigQuery smacks them
  # We don't expect this to be changed at runtime, so we can keep this as a constant defined in here
  Int broad_user_max_scatter = 1000

  # figure out max scatter depending on whether they're a Broad internal user or a beta customer.
  Int max_scatter_for_user =  if is_rate_limited_beta_customer then beta_customer_max_scatter
                              else broad_user_max_scatter

  if (!defined(git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker) || !defined(gatk_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  if (!load_vcf_headers && !load_vet_and_ref_ranges) {
    call Utils.TerminateWorkflow as MustLoadAtLeastOneThing {
      input:
        message = "GvsImportGenomes called with both load_vcf_headers and load_vet_and_ref_ranges set to false",
        basic_docker = effective_basic_docker,
    }
  }

  if ((num_samples > max_auto_batch_size) && !(defined(load_data_scatter_width))) {
    call Utils.TerminateWorkflow as DieDueToTooManySamplesWithoutExplicitLoadDataBatchSize {
      input:
        message = "Importing " + num_samples + " samples but 'load_data_scatter_width' is not explicitly specified; the limit for auto batch-sizing is " + max_auto_batch_size + " for " + genome_type + " samples.",
        basic_docker = effective_basic_docker,
    }
  }

  # At least 1, per limits above not more than 20.
  # But if it's a beta customer, use the number computed above
  Int effective_load_data_batch_size = if (defined(load_data_scatter_width)) then num_samples / select_first([load_data_scatter_width])
                                       else if num_samples < max_scatter_for_user then 1
                                         else if is_wgs then num_samples / max_scatter_for_user
                                           else if num_samples < 5001 then 20
                                             else if num_samples < 20001 then 100
                                               else if num_samples < 50001 then 500
                                                 else 1000

  # Both preemptible and maxretries should be scaled up alongside import batch size since the likelihood of preemptions
  # and retryable random BQ import errors increases with import batch size / job run time.

  # At least 3, per limits above not more than 5.
  Int effective_load_data_preemptible = if (defined(load_data_preemptible_override)) then select_first([load_data_preemptible_override])
                                        else if effective_load_data_batch_size < 12 then 3
                                          else effective_load_data_batch_size / 4

  Int effective_load_data_maxretries = select_first([load_data_maxretries_override, 5])

  call GetUningestedSampleIds {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      external_sample_names = external_sample_names,
      num_samples = num_samples,
      table_name = "sample_info",
      load_vet_and_ref_ranges = load_vet_and_ref_ranges,
      load_vcf_headers = load_vcf_headers,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call CurateInputLists {
    input:
      input_vcf_index_list = input_vcf_indexes,
      input_vcf_list = input_vcfs,
      input_sample_name_list = external_sample_names,
      input_samples_to_be_loaded_map = GetUningestedSampleIds.sample_map,
      variants_docker = effective_variants_docker,
  }

  call CreateFOFNs {
    input:
      batch_size = effective_load_data_batch_size,
      input_vcf_index_list = CurateInputLists.input_vcf_indexes,
      input_vcf_list = CurateInputLists.input_vcfs,
      sample_name_list = CurateInputLists.sample_name_list,
      basic_docker = effective_basic_docker,
  }

  scatter (i in range(length(CreateFOFNs.vcf_sample_name_fofns))) {
    call LoadData {
      input:
        index = i,
        dataset_name = dataset_name,
        project_id = project_id,
        skip_loading_vqsr_fields = skip_loading_vqsr_fields,
        drop_state = drop_state,
        drop_state_includes_greater_than = false,
        input_vcf_indexes = read_lines(CreateFOFNs.vcf_batch_vcf_index_fofns[i]),
        input_vcfs = read_lines(CreateFOFNs.vcf_batch_vcf_fofns[i]),
        interval_list = interval_list,
        gatk_docker = effective_gatk_docker,
        gatk_override = load_data_gatk_override,
        load_data_preemptible = effective_load_data_preemptible,
        load_data_maxretries = effective_load_data_maxretries,
        sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
        sample_map = GetUningestedSampleIds.sample_map,
        load_vet_and_ref_ranges = load_vet_and_ref_ranges,
        load_vcf_headers = load_vcf_headers,
        billing_project_id = billing_project_id,
        use_compressed_references = use_compressed_references,
    }
  }

  if (load_vcf_headers) {
    call ProcessVCFHeaders {
      input:
        variants_docker = effective_variants_docker,
        load_done = LoadData.done,
        dataset_name = dataset_name,
        project_id = project_id,
    }
  }

  if (load_vet_and_ref_ranges) {
    call SetIsLoadedColumn {
      input:
        load_done = LoadData.done,
        project_id = project_id,
        dataset_name = dataset_name,
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }
  }

  output {
    Boolean done = true
    Boolean used_tighter_gcp_quotas = is_rate_limited_beta_customer
    String recorded_git_hash = effective_git_hash
    Array[File] load_data_stderrs = LoadData.stderr
  }
}

task CreateFOFNs {
  input {
    Int batch_size
    File input_vcf_index_list
    File input_vcf_list
    File sample_name_list
    String basic_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    split -a 5 -l ~{batch_size} ~{input_vcf_list} batched_vcfs.
    split -a 5 -l ~{batch_size} ~{input_vcf_index_list} batched_vcf_indexes.
    split -a 5 -l ~{batch_size} ~{sample_name_list} batched_sample_names.
  >>>
  runtime {
    docker: basic_docker
    bootDiskSizeGb: 15
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[File] vcf_batch_vcf_fofns = glob("batched_vcfs.*")
    Array[File] vcf_batch_vcf_index_fofns = glob("batched_vcf_indexes.*")
    Array[File] vcf_sample_name_fofns = glob("batched_sample_names.*")
  }
}

task LoadData {
  input {
    Int index
    String dataset_name
    String project_id
    String? billing_project_id

    Array[File] input_vcf_indexes
    Array[File] input_vcfs
    File interval_list
    File sample_map
    Array[String] sample_names

    String? drop_state
    Boolean? drop_state_includes_greater_than = false
    Boolean force_loading_from_non_allele_specific = false
    Boolean skip_loading_vqsr_fields = false
    Boolean use_compressed_references = false
    Boolean load_vet_and_ref_ranges
    Boolean load_vcf_headers

    String gatk_docker
    File? gatk_override
    Int load_data_preemptible
    Int load_data_maxretries
  }

  meta {
    description: "Load data into BigQuery using the Write Api"
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }

    input_vcf_indexes: {
      localization_optional: true
    }
  }

  Int samples_per_table = 4000
  Int num_samples = length(sample_names)
  String temp_table = "~{dataset_name}.sample_names_to_load_~{index}"
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"
  String table_name = "sample_info"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # workaround for https://github.com/broadinstitute/cromwell/issues/3647
    export TMPDIR=/tmp

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    ## check which samples still need loading by looking in the BQ database for the loaded status of these samples

    # Create temp table with the sample_names and load external sample names into temp table -- make sure it doesn't exist already
    set +o errexit
    bq --apilog=false show --project_id=~{project_id} ~{temp_table} > /dev/null
    BQ_SHOW_RC=$?
    set -o errexit

    # If there is already a table of sample names or something else is wrong, burn it down to start fresh.
    if [ $BQ_SHOW_RC -eq 0 ]; then
      bq --apilog=false rm -t -f --project_id=~{project_id} ~{temp_table}
    fi

    echo "Creating the external sample name list table ~{temp_table}"
    bq --apilog=false --project_id=~{project_id} mk ~{temp_table} "sample_name:STRING"
    NAMES_FILE=~{write_lines(sample_names)}
    bq --apilog=false load --project_id=~{project_id} ~{temp_table} $NAMES_FILE "sample_name:STRING"

    # Get the current min/max id, or 0 if there are none. Withdrawn samples still have IDs so don't filter them out.
    # bq query --max_rows check: ok one row
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} '
      SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM `~{dataset_name}.~{table_name}`
        AS samples JOIN `~{temp_table}` AS temp ON samples.sample_name = temp.sample_name' > results.csv

    # Get sample map of samples that haven't been loaded yet
    # Break out individual queries into "status buckets" for all of the statuses we care about.
    for status in ~{true="REFERENCES_LOADED VARIANTS_LOADED" false="" load_vet_and_ref_ranges} ~{true="HEADERS_LOADED" false="" load_vcf_headers}
    do
      echo "
        SELECT sample_id, samples.sample_name FROM \`~{dataset_name}.~{table_name}\` AS samples JOIN \`~{temp_table}\` AS temp ON
        samples.sample_name = temp.sample_name WHERE
        samples.sample_id NOT IN (SELECT sample_id FROM \`~{dataset_name}.sample_load_status\` WHERE status = '$status') AND
        samples.withdrawn is NULL" > query.txt

      # bq query --max_rows check: ok sets max rows explicitly
      cat query.txt |
        bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} -n ~{num_samples} > \
        $status.status_bucket.csv
    done

    ## delete the table that was only needed for this ingest test
    bq --apilog=false --project_id=~{project_id} rm -f=true ~{temp_table}

    #  If a given sample shows up in any status bucket it should appear in the final sample map exactly once.
    # Add a header manually:
    echo "sample_id,sample_name" > sample_map.csv
    # The real header sorts to the bottom of the file, delete that.
    cat *.status_bucket.csv | sort -u | sed '$d' >> sample_map.csv

    ## now we want to create a sub list of these samples (without the ones that have already been loaded)

    python3 /gatk/scripts/variantstore/scripts/curate_input_array_files.py \
      --sample_map_to_be_loaded_file_name sample_map.csv \
      --sample_name_list_file_name $NAMES_FILE \
      --vcf_list_file_name ~{write_lines(input_vcfs)} \
      --vcf_index_list_file_name  ~{write_lines(input_vcf_indexes)}

    # translate files created by the python script into BASH arrays---but only of the samples that aren't there already
    VCFS_ARRAY=($(cat output_vcf_list_file |tr "\n" " "))
    VCF_INDEXES_ARRAY=($(cat output_vcf_index_list_file |tr "\n" " "))
    SAMPLE_NAMES_ARRAY=($(cat output_sample_name_list_file |tr "\n" " "))


    # loop over the BASH arrays (See https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value)
    for i in "${!VCFS_ARRAY[@]}"; do
      gs_input_vcf="${VCFS_ARRAY[$i]}"
      gs_input_vcf_index="${VCF_INDEXES_ARRAY[$i]}"
      sample_name="${SAMPLE_NAMES_ARRAY[$i]}"

      # we always do our own localization
      gcloud storage ~{"--billing-project " + billing_project_id} cp $gs_input_vcf input_vcf_$i.vcf.gz
      gcloud storage ~{"--billing-project " + billing_project_id} cp $gs_input_vcf_index input_vcf_$i.vcf.gz.tbi
      updated_input_vcf=input_vcf_$i.vcf.gz

      gatk --java-options "-Xmx2g" CreateVariantIngestFiles \
        -V ${updated_input_vcf} \
        -L ~{interval_list} \
        ~{"--ref-block-gq-to-ignore " + drop_state} \
        --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
        --force-loading-from-non-allele-specific ~{force_loading_from_non_allele_specific} \
        --project-id ~{project_id} \
        --dataset-name ~{dataset_name} \
        --output-type BQ \
        --enable-reference-ranges ~{load_vet_and_ref_ranges} \
        --enable-vet ~{load_vet_and_ref_ranges} \
        -SN ${sample_name} \
        -SNM ~{sample_map} \
        --ref-version 38 \
        --skip-loading-vqsr-fields ~{skip_loading_vqsr_fields} \
        --enable-vcf-headers ~{load_vcf_headers} \
        --use-compressed-refs ~{use_compressed_references}

      rm input_vcf_$i.vcf.gz
      rm input_vcf_$i.vcf.gz.tbi

    done
  >>>

  runtime {
    docker: gatk_docker
    maxRetries: load_data_maxretries
    memory: "3.75 GB"
    disks: "local-disk 50 HDD"
    preemptible: load_data_preemptible
    cpu: 1
    noAddress: true
  }
  output {
    Boolean done = true
    File stderr = stderr()
  }
}

task ProcessVCFHeaders {
  input {
    String dataset_name
    String project_id
    Array[String] load_done
    String variants_docker
  }
  meta {
    volatile: true
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    python3 /app/process_sample_vcf_headers.py \
      --project_id=~{project_id} \
      --dataset_name ~{dataset_name}
  >>>

  runtime {
    docker: variants_docker
    disks: "local-disk 500 HDD"
  }
}


task SetIsLoadedColumn {
  input {
    String dataset_name
    String project_id

    Array[String] load_done
    String cloud_sdk_docker
  }
  meta {
    # This is doing some tricky stuff with `INFORMATION_SCHEMA` so just punt and let it be `volatile`.
    volatile: true
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # set is_loaded to true if there is a corresponding vet table partition with rows for that sample_id

    # Note that we tried modifying CreateVariantIngestFiles to UPDATE sample_info.is_loaded on a per-sample basis.
    # The major issue that was found is that BigQuery allows only 20 such concurrent DML statements. Considered using
    # an exponential backoff, but at the number of samples that are being loaded this would introduce significant delays
    # in workflow processing. So this method is used to set *all* of the saple_info.is_loaded flags at one time.

    # bq query --max_rows check: ok update
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
    'UPDATE `~{dataset_name}.sample_info` SET is_loaded = true
    WHERE sample_id IN (SELECT CAST(partition_id AS INT64)
    from `~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS`
    WHERE partition_id NOT LIKE "__%" AND total_logical_bytes > 0 AND table_name LIKE "vet_%") OR sample_id IN
    (SELECT references_status.sample_id FROM `~{dataset_name}.sample_load_status` AS references_status
                     INNER JOIN `~{dataset_name}.sample_load_status` AS variants_status
                     ON references_status.sample_id = variants_status.sample_id
                     AND references_status.status = "REFERENCES_LOADED"
                     AND variants_status.status = "VARIANTS_LOADED")'
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    String done = "done"
  }
}

task GetUningestedSampleIds {
  input {
    String dataset_name
    String project_id

    File external_sample_names
    Int num_samples
    String table_name
    String cloud_sdk_docker
    # At least one of these "load" inputs must be true
    Boolean load_vet_and_ref_ranges
    Boolean load_vcf_headers
  }
  meta {
    # Do not call cache this, we want to read the database state every time.
    volatile: true
  }

  Int samples_per_table = 4000
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"
  String temp_table="~{dataset_name}.sample_names_to_load"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # Create temp table with the sample_names and load external sample names into temp table -- make sure it doesn't exist already
    set +o errexit
    bq --apilog=false show --project_id=~{project_id} ~{temp_table} > /dev/null
    BQ_SHOW_RC=$?
    set -o errexit

    # If there is already a table of sample names or something else is wrong, bail.
    if [ $BQ_SHOW_RC -eq 0 ]; then
      echo "There is already a list of sample names. This may need manual cleanup. Exiting"
      exit 1
    fi

    echo "Creating the external sample name list table ~{temp_table}"
    bq --apilog=false --project_id=~{project_id} mk ~{temp_table} "sample_name:STRING"
    bq --apilog=false load --project_id=~{project_id} ~{temp_table} ~{external_sample_names} "sample_name:STRING"

    # Get the current min/max id, or 0 if there are none. Withdrawn samples still have IDs so don't filter them out.
    # bq query --max_rows check: ok one row
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} '

      SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM `~{dataset_name}.~{table_name}`
        AS samples JOIN `~{temp_table}` AS temp ON samples.sample_name = temp.sample_name

    ' > results.csv

    # prep for being able to return min table id
    min_sample_id=$(tail -1 results.csv | cut -d, -f1)
    max_sample_id=$(tail -1 results.csv | cut -d, -f2)

    # no samples have been loaded or we don't have the right external_sample_names or something else is wrong, bail
    if [ $max_sample_id -eq 0 ]; then
      echo "Max id is 0. Exiting"
      exit 1
    fi

    python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_sample_id
    python3 -c "from math import ceil; print(ceil($min_sample_id/~{samples_per_table}))" > min_sample_id

    # Get sample map of samples that haven't been loaded yet
    # Break out individual queries into "status buckets" for all of the statuses we care about.
    for status in ~{true="REFERENCES_LOADED VARIANTS_LOADED" false="" load_vet_and_ref_ranges} ~{true="HEADERS_LOADED" false="" load_vcf_headers}
    do
      echo "
        SELECT sample_id, samples.sample_name FROM \`~{dataset_name}.~{table_name}\` AS samples JOIN \`~{temp_table}\` AS temp ON
        samples.sample_name = temp.sample_name WHERE
        samples.sample_id NOT IN (SELECT sample_id FROM \`~{dataset_name}.sample_load_status\` WHERE status = '$status') AND
        samples.withdrawn is NULL" > query.txt

      # bq query --max_rows check: ok sets max rows explicitly
      cat query.txt |
        bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} -n ~{num_samples} > \
        $status.status_bucket.csv
    done

    ## delete the table that was only needed for this ingest test
    bq --apilog=false --project_id=~{project_id} rm -f=true ~{temp_table}

    #  If a given sample shows up in any status bucket it should appear in the final sample map exactly once.
    # Add a header manually:
    echo "sample_id,sample_name" > sample_map.csv
    # The real header sorts to the bottom of the file, delete that.
    cat *.status_bucket.csv | sort -u | sed '$d' >> sample_map.csv

    cut -d, -f1 sample_map.csv > gvs_ids.csv

    ## delete the table that was only needed for this ingest
    bq --apilog=false --project_id=~{project_id} rm -f=true ~{temp_table}
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 5
    cpu: 1
  }
  output {
    Int max_table_id = ceil(read_float("max_sample_id"))
    Int min_table_id = ceil(read_float("min_sample_id"))
    File sample_map = "sample_map.csv"
    File gvs_ids = "gvs_ids.csv"
    Array[File] status_buckets = glob("*.status_bucket.csv")
  }
}

task CurateInputLists {
  input {
    File input_vcf_index_list
    File input_vcf_list
    File input_samples_to_be_loaded_map
    File input_sample_name_list
    String variants_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    python3 /app/curate_input_array_files.py --sample_map_to_be_loaded_file_name ~{input_samples_to_be_loaded_map} \
                                             --sample_name_list_file_name ~{input_sample_name_list} \
                                             --vcf_list_file_name ~{input_vcf_list} \
                                             --vcf_index_list_file_name  ~{input_vcf_index_list}
  >>>
  runtime {
    docker: variants_docker
    memory: "3 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    cpu: 1
  }

  output {
    File input_vcf_indexes = "output_vcf_index_list_file"
    File input_vcfs = "output_vcf_list_file"
    File sample_name_list = "output_sample_name_list_file"
  }
}
