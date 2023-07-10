version 1.0

import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset

workflow GvsUnified {
    input {
        # Begin GvsAssignIds
        String dataset_name
        String project_id
        String call_set_identifier

        File? gatk_override
        # End GvsAssignIds

        # Begin GvsImportGenomes
        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
        String drop_state = "NONE"
        # for beta users, rate limit their ingest to stay below quotas
        Boolean is_beta_user = false

        # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
        # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
        # those should work fine in most cases.
        Int? load_data_batch_size
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override
        # End GvsImportGenomes

        # Begin GvsCreateFilterSet
        String filter_set_name = call_set_identifier
        Boolean use_VQSR_lite = true

        Int? INDEL_VQSR_CLASSIC_max_gaussians_override = 4
        Int? INDEL_VQSR_CLASSIC_mem_gb_override
        Int? SNP_VQSR_CLASSIC_max_gaussians_override = 6
        Int? SNP_VQSR_CLASSIC_mem_gb_override
        # End GvsCreateFilterSet

        # Begin GvsPrepareRangesCallset
        String extract_table_prefix

        String query_project = project_id
        String destination_project = project_id
        String destination_dataset = dataset_name
        String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}"

        Array[String]? query_labels
        File? sample_names_to_extract
        # End GvsPrepareRangesCallset

        # Begin GvsExtractCallset
        Int? extract_scatter_count

        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"

        String extract_output_file_base_name = sub(filter_set_name, " ", "-")

        Int? extract_maxretries_override
        Int? extract_preemptible_override
        String? extract_output_gcs_dir
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Boolean extract_do_not_filter_override = false
        # End GvsExtractCallset
        String branch_name
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
    }

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            gatk_override = gatk_override,
            branch_name = branch_name,
            interval_list = interval_list,
            drop_state = drop_state,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
    }

    call PopulateAltAllele.GvsPopulateAltAllele {
        input:
            call_set_identifier = call_set_identifier,
            go = BulkIngestGenomes.done,
            dataset_name = dataset_name,
            project_id = project_id,
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = GvsPopulateAltAllele.done,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            filter_set_name = filter_set_name,
            use_VQSR_lite = use_VQSR_lite,
            interval_list = interval_list,
            gatk_override = gatk_override,
            INDEL_VQSR_CLASSIC_max_gaussians_override = INDEL_VQSR_CLASSIC_max_gaussians_override,
            INDEL_VQSR_CLASSIC_mem_gb_override = INDEL_VQSR_CLASSIC_mem_gb_override,
            SNP_VQSR_CLASSIC_max_gaussians_override = SNP_VQSR_CLASSIC_max_gaussians_override,
            SNP_VQSR_CLASSIC_mem_gb_override = SNP_VQSR_CLASSIC_mem_gb_override,
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            call_set_identifier = call_set_identifier,
            go = GvsCreateFilterSet.done,
            dataset_name = dataset_name,
            project_id = project_id,
            extract_table_prefix = extract_table_prefix,
            query_project = query_project,
            destination_project = destination_project,
            destination_dataset = destination_dataset,
            fq_temp_table_dataset = fq_temp_table_dataset,
            query_labels = query_labels,
            sample_names_to_extract = sample_names_to_extract,
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            extract_table_prefix = extract_table_prefix,
            filter_set_name = filter_set_name,
            query_project = query_project,
            scatter_count = extract_scatter_count,
            interval_list = interval_list,
            interval_weights_bed = interval_weights_bed,
            gatk_override = gatk_override,
            output_file_base_name = extract_output_file_base_name,
            extract_maxretries_override = extract_maxretries_override,
            extract_preemptible_override = extract_preemptible_override,
            output_gcs_dir = extract_output_gcs_dir,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
            do_not_filter_override = extract_do_not_filter_override,
            drop_state = drop_state,
    }

    output {
        Array[File] output_vcfs = GvsExtractCallset.output_vcfs
        Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
        Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
        Array[File] output_vcf_interval_files = GvsExtractCallset.output_vcf_interval_files
        File? sample_name_list = GvsExtractCallset.sample_name_list
        File manifest = GvsExtractCallset.manifest
        Boolean done = true
    }
}
