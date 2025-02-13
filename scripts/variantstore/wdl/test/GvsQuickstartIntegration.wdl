version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration
import "GvsQuickstartVATIntegration.wdl" as QuickstartVATIntegration
import "../GvsJointVariantCalling.wdl" as JointVariantCalling
import "../GvsUtils.wdl" as Utils

workflow GvsQuickstartIntegration {
    input {
        String git_branch_or_tag
        Boolean use_default_dockers = false
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
        Boolean run_exome_integration = true
        Boolean run_beta_integration = true
        Boolean run_bge_integration = true
        Boolean run_vat_integration = true
        Boolean run_vat_integration_test_from_vds = true        # If false, will use sites-only VCF
        String sample_id_column_name = "sample_id"
        String vcf_files_column_name = "hg38_reblocked_gvcf"
        String vcf_index_files_column_name = "hg38_reblocked_gvcf_index"
        String wgs_sample_set_name = "wgs_integration_sample_set"
        String exome_sample_set_name = "exome_integration_sample_set"
        String bge_sample_set_name = "bge_integration_sample_set"
        File? target_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/bge_TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list"
        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? variants_nirvana_docker
        String? gatk_docker
        String? hail_version
        Boolean chr20_X_Y_only = true
        Int? maximum_alternate_alleles
        String ploidy_table_name = "sample_chromosome_ploidy"
    }

    File full_wgs_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File full_exome_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"
    String expected_subdir = if (!chr20_X_Y_only) then "all_chrs/"  else ""
    File expected_output_prefix = "gs://gvs-internal-quickstart/integration/2024-10-29/" + expected_subdir
    File truth_data_prefix = "gs://gvs-internal-quickstart/integration/test_data/2025-01-17/" + expected_subdir

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
        File? none = ""
    }

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_variants_nirvana_docker = select_first([variants_nirvana_docker, GetToolVersions.variants_nirvana_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    if (chr20_X_Y_only && (run_hail_integration || run_vcf_integration)) {
        call FilterIntervalListChromosomes {
            input:
                full_interval_list = full_wgs_interval_list,
                chromosomes = ["chrX", "chrY", "chr20"],
                variants_docker = effective_variants_docker,
        }
    }

    if (!use_default_dockers) {
        call Utils.BuildGATKJar {
            input:
                git_branch_or_tag = git_branch_or_tag,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
        }
    }

    String workspace_bucket = GetToolVersions.workspace_bucket
    String workspace_id = GetToolVersions.workspace_id
    String submission_id = GetToolVersions.submission_id

    # Note for `GvsQuickstartIntegration` we use the git_branch_or_tag *input* and its corresponding git hash. This is not
    # necessarily the same as the branch name selected in Terra for the integration `GvsQuickstartIntegration` workflow,
    # though in practice likely they are the same.
    if (run_hail_integration) {
        # This test workflow is probably best representative of the AoU workflow. Parameters used here should be those used for AoU callsets
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVETSIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                extract_do_not_filter_override = false,
                dataset_suffix = "vets_hail",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = select_first([FilterIntervalListChromosomes.out, full_wgs_interval_list]),
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = wgs_sample_set_name,
                bgzip_output_vcfs = true,
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = workspace_bucket,
                workspace_id = workspace_id,
                submission_id = submission_id,
                hail_version = effective_hail_version,
                maximum_alternate_alleles = maximum_alternate_alleles,
                ploidy_table_name = ploidy_table_name,
        }

        if (GvsQuickstartHailVETSIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as HailVETSQuotaFail {
                input:
                    message = "GvsQuickstartHailVETSIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVETSIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VETS = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "vets_vcf",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = wgs_sample_set_name,
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = workspace_bucket,
                workspace_id = workspace_id,
                submission_id = submission_id,
                maximum_alternate_alleles = maximum_alternate_alleles,
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VETS = false,
                extract_do_not_filter_override = true,
                dataset_suffix = "vqsr_vcf",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = wgs_sample_set_name,
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = workspace_bucket,
                workspace_id = workspace_id,
                submission_id = submission_id,
                maximum_alternate_alleles = maximum_alternate_alleles,
        }

        if (QuickstartVcfVQSRIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as VcfVQSRQuotaFail {
                input:
                    message = "QuickstartVcfVQSRIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }

        if (QuickstartVcfVETSIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as VcfVETSQuotaFail {
                input:
                    message = "QuickstartVcfVETSIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_exome_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfExomeIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VETS = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "exome",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = false,
                interval_list = full_exome_interval_list,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = exome_sample_set_name,
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = workspace_bucket,
                workspace_id = workspace_id,
                submission_id = submission_id,
                maximum_alternate_alleles = maximum_alternate_alleles,
                target_interval_list = target_interval_list,
        }

        if (QuickstartVcfExomeIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as ExomeQuotaFail {
                input:
                    message = "QuickstartVcfExomeIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_bge_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfBgeIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VETS = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "bge",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = false,
                interval_list = full_exome_interval_list,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = bge_sample_set_name,
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = workspace_bucket,
                workspace_id = workspace_id,
                submission_id = submission_id,
                maximum_alternate_alleles = maximum_alternate_alleles,
                target_interval_list = target_interval_list,
        }

        if (QuickstartVcfBgeIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as BgeQuotaFail {
                input:
                    message = "QuickstartVcfBgeIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_beta_integration) {
        String project_id = "gvs-internal"

        String extract_output_gcs_dir = "~{workspace_bucket}/output_vcfs/by_submission_id/~{submission_id}/beta"
        Boolean collect_variant_calling_metrics = true

        call Utils.CreateDatasetForTest {
            input:
                git_branch_or_tag = git_branch_or_tag,
                dataset_prefix = "quickit",
                dataset_suffix = "beta",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call JointVariantCalling.GvsJointVariantCalling as QuickstartBeta {
            input:
                go = CreateDatasetForTest.done,
                call_set_identifier = git_branch_or_tag,
                dataset_name = CreateDatasetForTest.dataset_name,
                project_id = project_id,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                extract_output_file_base_name = "quickit",
                filter_set_name = "quickit",
                extract_table_prefix = "quickit",
                sample_set_name = wgs_sample_set_name,
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = workspace_bucket,
                workspace_id = workspace_id,
                submission_id = submission_id,
                maximum_alternate_alleles = maximum_alternate_alleles,
                git_branch_or_tag = git_branch_or_tag,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                extract_output_gcs_dir = extract_output_gcs_dir,
                collect_variant_calling_metrics = collect_variant_calling_metrics,
        }

        if (!QuickstartBeta.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as QuickstartBetaQuotaFail {
                input:
                    message = "QuickstartBeta should have used tighter GCP quotas but did not!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_vat_integration) {
        String extract_vat_output_gcs_dir = "~{workspace_bucket}/output_vat/by_submission_id/~{submission_id}/vat"

        call QuickstartVATIntegration.GvsQuickstartVATIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_default_dockers = use_default_dockers,
                truth_data_prefix = truth_data_prefix,
                expected_output_prefix = expected_output_prefix,
                dataset_suffix = "vat",
                output_path = extract_vat_output_gcs_dir,
                use_vds_as_input = run_vat_integration_test_from_vds,
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                variants_nirvana_docker = effective_variants_nirvana_docker,
                gatk_docker = effective_gatk_docker,
        }
    }

    output {
        String recorded_git_hash = GetToolVersions.git_hash
    }
}

task FilterIntervalListChromosomes {
    input {
        File full_interval_list
        Array[String]+ chromosomes
        String variants_docker
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python3 /app/filter_interval_list_chromosomes.py --input-interval-list ~{full_interval_list} \
          --output-interval-list "filtered.interval_list" --chromosome ~{sep=' --chromosome ' chromosomes}
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File out = "filtered.interval_list"
    }
}


