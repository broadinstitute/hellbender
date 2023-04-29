version 1.0

import "GvsUtils.wdl" as Utils
import "GvsPrepareBulkImport.wdl" as PrepareBulkImport
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes


workflow GvsBulkIngestGenomes {
    input {
        # Begin GvsPrepareBulkImport
        String terra_project_id # env vars
        String samples_table_name
        String sample_id_column_name
        String vcf_files_column_name
        String vcf_index_files_column_name
        # End GvsPrepareBulkImport

        # Begin GvsAssignIds
        String dataset_name
        String bq_project_id
        String call_set_identifier

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
        # End GvsImportGenomes
    }

    call GetWorkspaceId

    call GetWorkspaceName {
        input:
            workspace_id = GetWorkspaceId.workspace_bucket
    }

    call GetColumnNames

    call PrepareBulkImport.GvsPrepareBulkImport as PrepareBulkImport {
        input:
            project_id = terra_project_id,
            workspace_name = GetWorkspaceName.workspace_name,
            workspace_bucket = GetWorkspaceId.workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            dataset_name = dataset_name,
            project_id = bq_project_id,
            external_sample_names = read_lines(PrepareBulkImport.sampleFOFN),
            samples_are_controls = false
    }

    call ImportGenomes.GvsImportGenomes as ImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = bq_project_id,

            external_sample_names = read_lines(PrepareBulkImport.sampleFOFN),
            input_vcfs = read_lines(PrepareBulkImport.vcfFOFN),
            input_vcf_indexes = read_lines(PrepareBulkImport.vcfIndexFOFN),

            ## TODO the following will be kept short term for testing
            ## is_rate_limited_beta_customer = true,
            ## beta_customer_max_scatter = 25,

            interval_list = interval_list,

            # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
            # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
            # those should work fine in most cases.
            load_data_batch_size = load_data_batch_size,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_gatk_override = gatk_override,
    }

    output {
        Array[File] load_data_stderrs = ImportGenomes.load_data_stderrs
    }
}

    task GetColumnNames {
        command <<<
            # Get a list of all columns in the table

            # First we will check for the default named columns and make sure that each row has a value

            curl -X 'GET' \
            'https://rawls.dsde-prod.broadinstitute.org/api/workspaces/{workspace_id}/{workspace_name}/entities?useCache=true' \
            -H 'accept: application/json' \
            -H 'Authorization: Bearer XXX

            {
            "participant": {
            "attributeNames": [],
            "count": 249045,
            "idName": "participant_id"
            },
            "sample": {
            "attributeNames": [
            "blah blah these are the column names we want"
            ],
            "count": 249049,
            "idName": "sample_id"
            },
            "sample_set": {
            "attributeNames": [
            "samples"
            ],
            "count": 12,
            "idName": "sample_set_id"
            }
            }
        >>>

        runtime {
            docker: "ubuntu:latest"
        }

        output {
            Array[String] bag_of_column_names = []
        }
    }

    task GetWorkspaceName {
        input {
            String workspace_id
        }
        command <<<
            # Hit rawls with the workspace ID

            # python get_workspace_name()
            python3 -c "import requests; response=requests.get('https://rawls.dsde-prod.broadinstitute.org/api/workspaces/id/~{workspace_id}?fields=workspace.namespace,workspace.googleProject'); print(response.workspace.namespace)" > workspace_name.txt

    >>>

    runtime {
        docker: "ubuntu:latest"
    }

    output {
        #String terra_project_id = proj_id
        String workspace_name = read_string("workspace_name.txt")
        #String workspace_name = workspace_name
    }
}



    task GetWorkspaceId {
        command <<<
            # Prepend date, time and pwd to xtrace log entries.
            PS4='\D{+%F %T} \w $ '
            set -o errexit -o nounset -o pipefail -o xtrace

            # Sniff the workspace bucket out of the delocalization script and extract the workspace id from that.
            sed -n -E 's!.*gs://fc-(secure-)?([^\/]+).*!\2!p' /cromwell_root/gcs_delocalization.sh | sort -u > workspace_id.txt
            sed -n -E 's!.*gs://(fc-(secure-)?[^\/]+).*!\1!p' /cromwell_root/gcs_delocalization.sh | sort -u > workspace_bucket.txt
        >>>

        runtime {
            docker: "ubuntu:latest"
        }

        output {
            String workspace_id = read_string("workspace_id.txt")
            String workspace_bucket = read_string("workspace_bucket.txt")
        }
    }

