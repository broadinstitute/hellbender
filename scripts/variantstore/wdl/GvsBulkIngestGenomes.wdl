version 1.0

import "GvsUtils.wdl" as Utils
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes


workflow GvsBulkIngestGenomes {
    input {
        # Begin GenerateImportFofnFromDataTable
        # for now set the entity type names with a default
        String data_table_name = "sample" ## Note that it is possible an advanced user has a different name for the table. We could glean some information from the sample_set name if that has been defined, but this has not, and try to use that information instead of simply using the default "sample"
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
        # End GenerateImportFofnFromDataTable

        # Begin GvsAssignIds
        String dataset_name
        String project_id
        String? basic_docker
        String? cloud_sdk_docker
        String? variants_docker
        String? gatk_docker
        String? git_branch_or_tag
        String? git_hash

        File? gatk_override
        # End GvsAssignIds

        # Begin GvsImportGenomes
        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
        String drop_state = "NONE"

        # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable BigQuery errors,
        # so if specifying `load_data_batch_size`, adjust preemptible and maxretries accordingly. Or just take the defaults, as those should work fine in most cases.
        Int? load_data_batch_size
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override
        String? billing_project_id
        Boolean use_compressed_references = true
        # End GvsImportGenomes
    }

    parameter_meta {
        data_table_name: "The name of the data table; This table holds the GVCFs to be ingested; `sample` is the default."
        sample_id_column_name: "The column that will be used for the sample name / id in GVS; the <data_table_name>_id will be used as the default"
        vcf_files_column_name: "The column that supplies the path for the GVCFs to be ingested. If not specified, the workflow will attempt to derive the column name."
        vcf_index_files_column_name: "The column that supplies the path for the GVCF index files to be ingested. If not specified, the workflow will attempt to derive the column name."
        sample_set_name: "The recommended way to load samples; Sample sets must be created by the user. If no sample_set_name is specified, all samples will be loaded into GVS"
    }

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

    call GenerateImportFofnFromDataTable {
        input:
            variants_docker = effective_variants_docker,
            sample_set_name = sample_set_name,
            data_table_name = data_table_name,
            user_defined_sample_id_column_name = sample_id_column_name, ## NOTE: the user needs to define this, or it will default to the <entity>_id column
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
    }

    call SplitBulkImportFofn {
        input:
            import_fofn = GenerateImportFofnFromDataTable.output_fofn,
            basic_docker = effective_basic_docker,
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = SplitBulkImportFofn.sample_name_fofn,
            samples_are_controls = false,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call ImportGenomes.GvsImportGenomes as ImportGenomes {
        input:
            go = AssignIds.done,
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = SplitBulkImportFofn.sample_name_fofn,
            num_samples = SplitBulkImportFofn.sample_num,
            input_vcfs = SplitBulkImportFofn.vcf_file_name_fofn,
            input_vcf_indexes = SplitBulkImportFofn.vcf_index_file_name_fofn,
            interval_list = interval_list,

            # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
            # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
            # those should work fine in most cases.
            load_data_batch_size = load_data_batch_size,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            load_data_gatk_override = gatk_override,
            drop_state = drop_state,
            billing_project_id = billing_project_id,
    }

    output {
        Boolean done = true
        String workspace_bucket = GenerateImportFofnFromDataTable.workspace_bucket
        String workspace_id = GenerateImportFofnFromDataTable.workspace_id
        String submission_id = GenerateImportFofnFromDataTable.submission_id
        String workflow_id = GenerateImportFofnFromDataTable.workflow_id
        String recorded_git_hash = effective_git_hash
    }
}


task GenerateImportFofnFromDataTable {
    ## In order to get the names of the columns with the GVCF and GVCF Index file paths, without requiring that the user input it manually, we apply heuristics
    input {
        String data_table_name ## NOTE: if not specified by the user, this has been set to "sample"
        String? sample_set_name
        String? user_defined_sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String variants_docker
    }

    ## set some output files
    String vcf_files_column_name_output_file = "vcf_files_column_name.txt"
    String vcf_index_files_column_name_output_file = "vcf_index_files_column_name.txt"

    String entity_id = data_table_name + "_id"

    String output_fofn_name = "output.tsv"
    String error_file_name = "errors.txt"

    String workspace_id_output = "workspace_id.txt"
    String workspace_name_output = "workspace_name.txt"
    String workspace_namespace_output = "workspace_namespace.txt"
    String workspace_bucket_output = "workspace_bucket.txt"
    String submission_id_output = "submission_id.txt"
    String workflow_id_output = "workflow_id.txt"

    String sample_name_column = if (defined(user_defined_sample_id_column_name)) then select_first([user_defined_sample_id_column_name]) else entity_id

    command <<<

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Sniff the workspace bucket out of the delocalization script and extract the workspace id from that.
        sed -n -E 's!.*gs://fc-(secure-)?([^\/]+).*!\2!p' /cromwell_root/gcs_delocalization.sh | sort -u > ~{workspace_id_output}
        sed -n -E 's!.*gs://(fc-(secure-)?[^\/]+).*!\1!p' /cromwell_root/gcs_delocalization.sh | sort -u > ~{workspace_bucket_output}
        sed -n -E 's!.*gs://fc-(secure-)?([^\/]+)/submissions/([^\/]+).*!\3!p' /cromwell_root/gcs_delocalization.sh | sort -u > ~{submission_id_output}
        sed -n -E 's!.*gs://fc-(secure-)?([^\/]+)/submissions/([^\/]+)/([^\/]+)/([^\/]+).*!\5!p' /cromwell_root/gcs_delocalization.sh | sort -u > ~{workflow_id_output}

        export WORKSPACE_ID="$(cat '~{workspace_id_output}')"
        export WORKSPACE_BUCKET="$(cat '~{workspace_bucket_output}')"

        # Hit rawls with the workspace ID

        python3 /app/get_workspace_name_for_import.py \
            --workspace_id ${WORKSPACE_ID} \
            --workspace_name_output '~{workspace_name_output}' \
            --workspace_namespace_output '~{workspace_namespace_output}'

        export WORKSPACE_NAME="$(cat '~{workspace_name_output}')"
        export WORKSPACE_NAMESPACE="$(cat '~{workspace_namespace_output}')"

        # Get a list of all columns in the table. Apply basic heuristics to write the resulting vcf_files_column_name and vcf_index_files_column_name.

        python3 /app/get_columns_for_import.py \
             ~{"--user_defined_sample_id " + user_defined_sample_id_column_name} \
             ~{"--entity_set_name " + sample_set_name} \
             ~{"--user_defined_vcf " + vcf_files_column_name} \
             ~{"--user_defined_index " + vcf_index_files_column_name} \
            --entity_type ~{data_table_name} \
            --vcf_output ~{vcf_files_column_name_output_file} \
            --vcf_index_output ~{vcf_index_files_column_name_output_file}

        if [[ -z "~{vcf_files_column_name}" ]]
        then
            export VCF_COLUMN_NAME="$(cat ~{vcf_files_column_name_output_file})"
        else
            export VCF_COLUMN_NAME="~{vcf_files_column_name}"
        fi

        if [[ -z "~{vcf_index_files_column_name}" ]]
        then
            export VCF_INDEX_COLUMN_NAME="$(cat ~{vcf_index_files_column_name_output_file})"
        else
            export VCF_INDEX_COLUMN_NAME="~{vcf_index_files_column_name}"
        fi

        export GOOGLE_PROJECT="${WORKSPACE_NAMESPACE}"
        python3 /app/generate_fofn_for_import.py \
            --data-table-name ~{data_table_name} \
            --sample-id-column-name ~{sample_name_column} \
            --vcf-files-column-name "${VCF_COLUMN_NAME}" \
            --vcf-index-files-column-name "${VCF_INDEX_COLUMN_NAME}" \
            ~{"--sample-set-name " + sample_set_name} \
            --output-file-name ~{output_fofn_name} \
            --error-file-name ~{error_file_name}

        if [ -s ~{error_file_name} ]; then
            echo ""
            echo "-------- the following issues were found with the sample data, no samples were ingested in this run --------"
            cat ~{error_file_name}
            echo ""
            exit 1
        fi

    >>>

    runtime {
        docker: variants_docker
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        String workspace_bucket = read_string("~{workspace_bucket_output}")
        String workspace_id = read_string("~{workspace_id_output}")
        String submission_id = read_string("~{submission_id_output}")
        String workflow_id = read_string("~{workflow_id_output}")
        File output_fofn = output_fofn_name
    }
}


task SplitBulkImportFofn {
    input {
        File import_fofn
        String basic_docker
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        PS4='\D{+%F %T} \w $ '

        cut -f 1 ~{import_fofn} > sample_names.txt
        cut -f 2 ~{import_fofn} > vcf_file_names.txt
        cut -f 3 ~{import_fofn} > vcf_index_file_names.txt
        wc -l < ~{import_fofn} > sample_num.txt
    >>>

    runtime {
        docker: basic_docker
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File sample_name_fofn = "sample_names.txt"
        File vcf_file_name_fofn = "vcf_file_names.txt"
        File vcf_index_file_name_fofn = "vcf_index_file_names.txt"
        Int sample_num = read_int("sample_num.txt")
    }
}
