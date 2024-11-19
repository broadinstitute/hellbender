version 1.0

import "GvsCreateTables.wdl" as GvsCreateTables
import "GvsUtils.wdl" as Utils

workflow GvsAssignIds {

  input {
    Boolean go = true
    String? git_branch_or_tag
    String? git_hash
    String dataset_name
    String project_id

    File external_sample_names
    Boolean samples_are_controls = false

    Boolean load_vet_and_ref_ranges
    Boolean load_vcf_headers

    Boolean use_compressed_references = false
    String? cloud_sdk_docker
  }

  String sample_info_table = "sample_info"
  String sample_info_schema_json = '[{"name": "sample_name","type": "STRING","mode": "REQUIRED"},{"name": "sample_id","type": "INTEGER","mode": "NULLABLE"},{"name":"is_loaded","type":"BOOLEAN","mode":"NULLABLE"},{"name":"is_control","type":"BOOLEAN","mode":"REQUIRED"},{"name":"withdrawn","type":"TIMESTAMP","mode":"NULLABLE"}]'
  String vcf_header_lines_scratch_schema_json = '[{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name":"vcf_header_lines","type":"STRING","mode":"NULLABLE"},{"name":"vcf_header_lines_hash","type":"STRING","mode":"REQUIRED"},{"name":"is_expected_unique","type":"BOOLEAN","mode":"REQUIRED"}]'
  String vcf_header_lines_schema_json = '[{"name":"vcf_header_lines_hash","type":"STRING","mode":"REQUIRED"}, {"name":"vcf_header_lines","type":"STRING","mode":"REQUIRED"},{"name":"is_expected_unique","type":"BOOLEAN","mode":"REQUIRED"}]'
  String sample_vcf_header_schema_json = '[{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"}, {"name":"vcf_header_lines_hash","type":"STRING","mode":"REQUIRED"}]'
  String sample_load_status_schema_json = '[{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name":"status","type":"STRING","mode":"REQUIRED"}, {"name":"event_timestamp","type":"TIMESTAMP","mode":"REQUIRED"}]'
  # The schema below for the `sample_chromosome_ploidy` table matches that in the AoU `echo` dataset, which differs from
  # the standard `sample_chromosome_ploidy` schema by the addition of the `genotype` column. At the time of this writing
  # it is not expected that any additional samples will be loaded into the `echo` dataset, but if there are additional
  # samples then this code will match the actual `echo` schema.
  String sample_chromosome_ploidy_schema_json = '[{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name": "chromosome","type": "INTEGER","mode": "REQUIRED"},{"name":"genotype","type":"STRING","mode":"NULLABLE"},{"name": "ploidy","type": "INTEGER","mode": "REQUIRED"}]'

  if (!defined(git_hash) || !defined(cloud_sdk_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  call ValidateSamples {
    input:
      sample_names_file = external_sample_names,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GvsCreateTables.CreateTables as CreateSampleInfoTable {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      go = ValidateSamples.done,
      datatype = "sample_info",
      schema_json = sample_info_schema_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false",
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GvsCreateTables.CreateTables as CreateSampleLoadStatusTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      go = ValidateSamples.done,
      datatype = "sample_load_status",
      schema_json = sample_load_status_schema_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false",
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GvsCreateTables.CreateTables as CreateSamplePloidyMapTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      go = ValidateSamples.done,
      datatype = "sample_chromosome_ploidy",
      schema_json = sample_chromosome_ploidy_schema_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false",
      clustering_field = "sample_id",
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  if (load_vcf_headers) {
    call GvsCreateTables.CreateTables as CreateScratchVCFHeaderLinesTable {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        go = ValidateSamples.done,
        datatype = "vcf_header_lines_scratch",
        schema_json = vcf_header_lines_scratch_schema_json,
        max_table_id = 1,
        superpartitioned = "false",
        partitioned = "false",
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call GvsCreateTables.CreateTables as CreateVCFHeaderLinesTable {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        go = ValidateSamples.done,
        datatype = "vcf_header_lines",
        schema_json = vcf_header_lines_schema_json,
        max_table_id = 1,
        superpartitioned = "false",
        partitioned = "false",
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call GvsCreateTables.CreateTables as CreateSampleVCFHeaderTable {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        go = ValidateSamples.done,
        datatype = "sample_vcf_header",
        schema_json = sample_vcf_header_schema_json,
        max_table_id = 1,
        superpartitioned = "false",
        partitioned = "false",
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }
  }

  call CreateCostObservabilityTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      go = ValidateSamples.done,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call AssignIds {
    input:
      sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      sample_info_table = sample_info_table,
      samples_are_controls = samples_are_controls,
      table_creation_done = CreateSampleInfoTable.done,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  if (load_vet_and_ref_ranges) {
    call GvsCreateTables.CreateBQTables as CreateTablesForMaxId {
      input:
        project_id = project_id,
        git_branch_or_tag = git_branch_or_tag,
        git_hash = effective_git_hash,
        dataset_name = dataset_name,
        max_table_id = AssignIds.max_table_id,
        cloud_sdk_docker = effective_cloud_sdk_docker,
        use_compressed_references = use_compressed_references,
    }
  }

  output {
    Boolean done = true
    File gvs_ids_tsv = AssignIds.gvs_ids_tsv
    String recorded_git_hash = effective_git_hash
  }
}

task AssignIds {
  input {
    String dataset_name
    String project_id

    String sample_info_table
    File sample_names
    Boolean samples_are_controls
    Boolean table_creation_done
    String cloud_sdk_docker
  }
  meta {
    description: "Assigns Ids to samples"
    # Not `volatile: true` since successfully assigned IDs should not be assigned again.
  }

  Int samples_per_table = 4000
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:assign_ids"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    num_samples=$(wc -l < ~{sample_names})

    # make sure that sample names were actually passed, fail if empty
    if [ $num_samples -eq 0 ]; then
      echo "No sample names found in file ~{sample_names}. Exiting"
      exit 1
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # Check that lock table has not been created yet
    set +e
    bq --apilog=false --project_id=~{project_id} show ~{dataset_name}.sample_id_assignment_lock > /dev/null
    BQ_SHOW_RC=$?
    set -e
    if [ $BQ_SHOW_RC -eq 0 ]; then
      echo "lock table already exists. Exiting"
      exit 1
    fi

    # create the lock table
    bq --apilog=false --project_id=~{project_id} mk ~{dataset_name}.sample_id_assignment_lock "sample_name:STRING"

    # first load name into the lock table - will check for dupes when adding to sample_info table
    bq --apilog=false load --project_id=~{project_id} ~{dataset_name}.sample_id_assignment_lock ~{sample_names} "sample_name:STRING"

    # add sample_name to sample_info_table
    # bq query --max_rows check: ok insert
    bq --apilog=false --project_id=~{project_id} query --use_legacy_sql=false ~{bq_labels} \
      'INSERT into `~{dataset_name}.~{sample_info_table}` (sample_name, is_control) select sample_name, ~{samples_are_controls} from `~{dataset_name}.sample_id_assignment_lock` m where m.sample_name not in (SELECT sample_name FROM `~{dataset_name}.~{sample_info_table}`)'

    # get the current maximum id, or 0 if there are none
    # bq query --max_rows check: ok single row
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} 'SELECT IFNULL(MAX(sample_id),0) FROM `~{dataset_name}.~{sample_info_table}`' > maxid
    offset=$(tail -1 maxid)

    # perform actual id assignment
    # bq query --max_rows check: ok update
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} --parameter=offset:INTEGER:$offset \
      'UPDATE `~{dataset_name}.~{sample_info_table}` m SET m.sample_id = id_assign.id FROM (SELECT sample_name, @offset + ROW_NUMBER() OVER(order by sample_name) as id FROM `~{dataset_name}.~{sample_info_table}` WHERE sample_id IS NULL) id_assign WHERE m.sample_name = id_assign.sample_name;'

    # retrieve the list of assigned ids and samples to update the datamodel
    echo "entity:sample_id,gvs_id" > update.tsv
    # bq query --max_rows check: ok num samples explicit
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} -n $num_samples --parameter=offset:INTEGER:$offset \
      'SELECT sample_name, sample_id from `~{dataset_name}.~{sample_info_table}` WHERE sample_id >= @offset' > update.tsv
    cat update.tsv | sed -e 's/sample_id/gvs_id/' -e 's/sample_name/entity:sample_id/' -e 's/,/\t/g' > gvs_ids.tsv

    # Get the max id for which to create tables.
    # We can't safely pipe to `head -1` because while `head` will exit successfully after reading the first line, the
    # pipeline will continue trying to write data to the `head` process. If this happens we'll get a 141 exit code and
    # with `set -o pipefail` turned on this will fail our task. As a workaround use this `<(...)` temp file construct.
    # https://news.ycombinator.com/item?id=9255830
    max_sample_id=$(head -1 <(cat update.tsv | cut -d, -f2 | sort -r -n))
    python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_table_id

    # remove the lock table
    bq --apilog=false --project_id=~{project_id} rm -f -t ~{dataset_name}.sample_id_assignment_lock
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "3.75 GB"
    disks: "local-disk " + 10 + " HDD"
    cpu: 1
  }
  output {
    File gvs_ids_tsv = "gvs_ids.tsv"
    Int max_table_id = read_int("max_table_id")
  }
}

task CreateCostObservabilityTable {
  input {
    String project_id
    String dataset_name
    Boolean go
    String cloud_sdk_docker
  }

  String cost_observability_json =
                                 '[ { "name": "call_set_identifier", "type": "STRING", "mode": "REQUIRED", "description": "The name by which we refer to the callset" }, ' +
                                 '  { "name": "step", "type": "STRING", "mode": "REQUIRED", "description": "The name of the core GVS workflow to which this belongs" }, ' +
                                 '  { "name": "call", "type": "STRING", "mode": "NULLABLE", "description": "The WDL call to which this belongs" }, ' +
                                 '  { "name": "shard_identifier", "type": "STRING", "mode": "NULLABLE", "description": "A unique identifier for this shard, may or may not be its index" }, ' +
                                 '  { "name": "call_start_timestamp", "type": "TIMESTAMP", "mode": "REQUIRED", "description": "When the call logging this event was started" }, ' +
                                 '  { "name": "event_timestamp", "type": "TIMESTAMP", "mode": "REQUIRED", "description": "When the observability event was logged" }, ' +
                                 '  { "name": "event_key", "type": "STRING", "mode": "REQUIRED", "description": "The type of observability event being logged" }, ' +
                                 '  { "name": "event_bytes", "type": "INTEGER", "mode": "REQUIRED", "description": "Number of bytes reported for this observability event" } ] '

  meta {
    # not volatile: true, always run this when asked
  }
  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    TABLE="~{dataset_name}.cost_observability"

    # Check that the table has not been created yet
    set +o errexit
    bq --apilog=false show --project_id=~{project_id} $TABLE > /dev/null
    BQ_SHOW_RC=$?
    set -o errexit

    if [ $BQ_SHOW_RC -ne 0 ]; then
      PARTITION_STRING="--time_partitioning_field call_start_timestamp --time_partitioning_type DAY"
      echo "making table $TABLE"
      echo '~{cost_observability_json}' > schema.json
      bq --apilog=false mk ${PARTITION_STRING} --project_id=~{project_id} $TABLE schema.json
    fi
  >>>
  runtime {
    docker: cloud_sdk_docker
  }
  output {
    Boolean done = true
  }
}

task ValidateSamples {
  input {
    File sample_names_file
    String cloud_sdk_docker
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    if [[ ! -s ~{sample_names_file} ]]
    then
      echo "ERROR: The input file ~{sample_names_file} is empty"
      exit 1;
    fi

    sort ~{sample_names_file} | uniq -d > output.txt
    if [[ -s output.txt ]]
    then
      echo "ERROR: The input file ~{sample_names_file} contains the following duplicate entries:"
      cat output.txt
      exit 1;
    fi

  >>>

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    cpu: "1"
    preemptible: 1
    maxRetries: 0
    disks: "local-disk 100 HDD"
  }

  output {
    Boolean done = true
  }
}
