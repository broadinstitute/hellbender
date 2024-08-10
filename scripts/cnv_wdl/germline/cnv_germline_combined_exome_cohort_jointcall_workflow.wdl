version 1.0

import "cnv_germline_cohort_workflow.wdl" as CohortWorkflow
import "../cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCombinedCohortJointcalling {

    input {
        #####
        ##### Cohort workflow
        #####

          ##################################
          #### required basic arguments ####
          ##################################
          File intervals
          File? blacklist_intervals
          Array[String]+ normal_bams
          Array[String]+ normal_bais
          String cohort_entity_id
          File contig_ploidy_priors
          Int num_intervals_per_scatter
          File ref_fasta_dict
          File ref_fasta_fai
          File ref_fasta
          String gatk_docker

          ##################################
          #### optional basic arguments ####
          ##################################
          # If true, AnnotateIntervals will be run to create GC annotations and explicit
          # GC correction will be performed by the model generated by
          Boolean? do_explicit_gc_correction
          File? gatk4_jar_override
          Int? preemptible_attempts

          # Required if BAM/CRAM is in a requester pays bucket
          String? gcs_project_for_requester_pays

          ####################################################
          #### optional arguments for PreprocessIntervals ####
          ####################################################
          Int? padding
          Int? bin_length

          ##################################################
          #### optional arguments for AnnotateIntervals ####
          ##################################################
          File? mappability_track_bed
          File? mappability_track_bed_idx
          File? segmental_duplication_track_bed
          File? segmental_duplication_track_bed_idx
          Int? feature_query_lookahead
          Int? mem_gb_for_annotate_intervals

          #################################################
          #### optional arguments for FilterIntervals ####
          ################################################
          File? blacklist_intervals_for_filter_intervals
          Float? minimum_gc_content
          Float? maximum_gc_content
          Float? minimum_mappability
          Float? maximum_mappability
          Float? minimum_segmental_duplication_content
          Float? maximum_segmental_duplication_content
          Int? low_count_filter_count_threshold
          Float? low_count_filter_percentage_of_samples
          Float? extreme_count_filter_minimum_percentile
          Float? extreme_count_filter_maximum_percentile
          Float? extreme_count_filter_percentage_of_samples
          Int? mem_gb_for_filter_intervals

          ##############################################
          #### optional arguments for CollectCounts ####
          ##############################################
          Array[String]? disabled_read_filters_for_collect_counts
          String? collect_counts_format
          Boolean? collect_counts_enable_indexing
          Int? mem_gb_for_collect_counts

          ########################################################################
          #### optional arguments for DetermineGermlineContigPloidyCohortMode ####
          ########################################################################
          Float? ploidy_mean_bias_standard_deviation
          Float? ploidy_mapping_error_rate
          Float? ploidy_global_psi_scale
          Float? ploidy_sample_psi_scale
          Int? mem_gb_for_determine_germline_contig_ploidy
          Int? cpu_for_determine_germline_contig_ploidy

          ############################################################
          #### optional arguments for GermlineCNVCallerCohortMode ####
          ############################################################
          Float? gcnv_p_alt
          Float? gcnv_p_active
          Float? gcnv_cnv_coherence_length
          Float? gcnv_class_coherence_length
          Int? gcnv_max_copy_number
          Int? mem_gb_for_germline_cnv_caller
          Int? cpu_for_germline_cnv_caller

          # optional arguments for germline CNV denoising model
          Int? gcnv_max_bias_factors
          Float? gcnv_mapping_error_rate
          Float? gcnv_interval_psi_scale
          Float? gcnv_sample_psi_scale
          Float? gcnv_depth_correction_tau
          Float? gcnv_log_mean_bias_standard_deviation
          Float? gcnv_init_ard_rel_unexplained_variance
          Int? gcnv_num_gc_bins
          Float? gcnv_gc_curve_standard_deviation
          String? gcnv_copy_number_posterior_expectation_mode
          Boolean? gcnv_enable_bias_factors
          Int? gcnv_active_class_padding_hybrid_mode

          # optional arguments for Hybrid ADVI
          Float? gcnv_learning_rate
          Float? gcnv_adamax_beta_1
          Float? gcnv_adamax_beta_2
          Int? gcnv_log_emission_samples_per_round
          Float? gcnv_log_emission_sampling_median_rel_error
          Int? gcnv_log_emission_sampling_rounds
          Int? gcnv_max_advi_iter_first_epoch
          Int? gcnv_max_advi_iter_subsequent_epochs
          Int? gcnv_min_training_epochs
          Int? gcnv_max_training_epochs
          Float? gcnv_initial_temperature
          Int? gcnv_num_thermal_advi_iters
          Int? gcnv_convergence_snr_averaging_window
          Float? gcnv_convergence_snr_trigger_threshold
          Int? gcnv_convergence_snr_countdown_window
          Int? gcnv_max_calling_iters
          Float? gcnv_caller_update_convergence_threshold
          Float? gcnv_caller_internal_admixing_rate
          Float? gcnv_caller_external_admixing_rate
          Boolean? gcnv_disable_annealing

          ###################################################
          #### arguments for PostprocessGermlineCNVCalls ####
          ###################################################
          Int ref_copy_number_autosomal_contigs
          Int? mem_gb_for_postprocess_germline_cnv_calls
          Int? disk_space_gb_for_postprocess_germline_cnv_calls
          Array[String]? allosomal_contigs

          ##########################
          #### arguments for QC ####
          ##########################
          Int maximum_number_events_per_sample
          Int maximum_number_pass_events_per_sample







        #####
        ##### inputs for joint calling
        #####

        ##################################
        #### required basic arguments ####
        ##################################
              Int num_samples_per_scatter_block
              File intervals
              File? blacklist_intervals

#              File contig_ploidy_calls_tar_path_list
#              File gcnv_calls_tars_path_list
#              File genotyped_intervals_vcf_indexes_path_list
#              File genotyped_intervals_vcfs_path_list
#              File genotyped_segments_vcf_indexes_path_list
#              File genotyped_segments_vcfs_path_list

              #qc arguments
              Int maximum_number_events
              Int maximum_number_pass_events

#              Array[File] gcnv_model_tars
#              Array[File] calling_configs
#              Array[File] denoising_configs
#              Array[File] gcnvkernel_version
#              Array[File] sharded_interval_lists
              Array[String]? allosomal_contigs
              Int ref_copy_number_autosomal_contigs
              File ref_fasta_dict
              File ref_fasta_fai
              File ref_fasta
              String x_contig_name
              File protein_coding_gtf
              File linc_rna_gtf
              File promoter_bed
              File noncoding_bed
              String gatk_docker
              String gatk_docker_clustering
              String gatk_docker_qual_calc
              String sv_pipeline_docker



    }

    call CohortWorkflow.CNVGermlineCohortWorkflow as CohortWF {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            normal_bams = normal_bams,
            normal_bais = normal_bais,
            cohort_entity_id = cohort_entity_id,
            contig_ploidy_priors = contig_ploidy_priors,
            num_intervals_per_scatter = num_intervals_per_scatter,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta = ref_fasta,
            gatk_docker = gatk_docker,
            do_explicit_gc_correction = do_explicit_gc_correction,
            gatk4_jar_override = gatk4_jar_override,
            preemptible_attempts = preemptible_attempts,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays,
            padding = padding,
            bin_length = bin_length,

            mappability_track_bed = mappability_track_bed,
            mappability_track_bed_idx = mappability_track_bed_idx,

            segmental_duplication_track_bed = segmental_duplication_track_bed,
            segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
            feature_query_lookahead = feature_query_lookahead,
            mem_gb_for_annotate_intervals = mem_gb_for_annotate_intervals,

            blacklist_intervals_for_filter_intervals = blacklist_intervals_for_filter_intervals,
            minimum_gc_content = minimum_gc_content,
            maximum_gc_content = maximum_gc_content,
            minimum_mappability = minimum_mappability,
            maximum_mappability = maximum_mappability,
            minimum_segmental_duplication_content = minimum_segmental_duplication_content,
            maximum_segmental_duplication_content = maximum_segmental_duplication_content,
            low_count_filter_count_threshold = low_count_filter_count_threshold,
            low_count_filter_percentage_of_samples = low_count_filter_percentage_of_samples,
            extreme_count_filter_minimum_percentile = extreme_count_filter_minimum_percentile,
            extreme_count_filter_maximum_percentile = extreme_count_filter_maximum_percentile,
            extreme_count_filter_percentage_of_samples = extreme_count_filter_percentage_of_samples,
            mem_gb_for_filter_intervals = mem_gb_for_filter_intervals,

            disabled_read_filters_for_collect_counts = disabled_read_filters_for_collect_counts,
            collect_counts_format = collect_counts_format,
            collect_counts_enable_indexing = collect_counts_enable_indexing,
            mem_gb_for_collect_counts = mem_gb_for_collect_counts,

            ploidy_mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
            ploidy_mapping_error_rate = ploidy_mapping_error_rate,
            ploidy_global_psi_scale = ploidy_global_psi_scale,
            ploidy_sample_psi_scale = ploidy_sample_psi_scale,
            mem_gb_for_determine_germline_contig_ploidy = mem_gb_for_determine_germline_contig_ploidy,
            cpu_for_determine_germline_contig_ploidy = cpu_for_determine_germline_contig_ploidy,

            gcnv_p_alt = gcnv_p_alt,
            gcnv_p_active = gcnv_p_active,
            gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
            gcnv_class_coherence_length = gcnv_class_coherence_length,
            gcnv_max_copy_number = gcnv_max_copy_number,
            mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
            cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,

            gcnv_max_bias_factors = gcnv_max_bias_factors,
            gcnv_mapping_error_rate = gcnv_mapping_error_rate,
            gcnv_interval_psi_scale = gcnv_interval_psi_scale,
            gcnv_sample_psi_scale = gcnv_sample_psi_scale,
            gcnv_depth_correction_tau = gcnv_depth_correction_tau,
            gcnv_log_mean_bias_standard_deviation = gcnv_log_mean_bias_standard_deviation,
            gcnv_init_ard_rel_unexplained_variance = gcnv_init_ard_rel_unexplained_variance,
            gcnv_num_gc_bins = gcnv_num_gc_bins,
            gcnv_gc_curve_standard_deviation = gcnv_gc_curve_standard_deviation,
            gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
            gcnv_enable_bias_factors = gcnv_enable_bias_factors,
            gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,

            gcnv_learning_rate = gcnv_learning_rate,
            gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
            gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
            gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
            gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
            gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
            gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
            gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
            gcnv_min_training_epochs = gcnv_min_training_epochs,
            gcnv_max_training_epochs = gcnv_max_training_epochs,
            gcnv_initial_temperature = gcnv_initial_temperature,
            gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
            gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
            gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
            gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
            gcnv_max_calling_iters = gcnv_max_calling_iters,
            gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
            gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
            gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
            gcnv_disable_annealing = gcnv_disable_annealing,

            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            mem_gb_for_postprocess_germline_cnv_calls = mem_gb_for_postprocess_germline_cnv_calls,
            disk_space_gb_for_postprocess_germline_cnv_calls = disk_space_gb_for_postprocess_germline_cnv_calls,
            allosomal_contigs = allosomal_contigs,

            maximum_number_events_per_sample = maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample,


    }

    File contig_ploidy_calls_tar_path_list = CohortWF.contig_ploidy_calls_tar_path_list
    File gcnv_calls_tars_path_list = CohortWF.gcnv_calls_tars_path_list
    File genotyped_intervals_vcf_indexes_path_list = CohortWF.genotyped_intervals_vcf_indexes_path_list
    File genotyped_intervals_vcfs_path_list = CohortWF.genotyped_intervals_vcfs_path_list
    File genotyped_segments_vcf_indexes_path_list = CohortWF.genotyped_segments_vcf_indexes_path_list
    File genotyped_segments_vcfs_path_list = CohortWF.genotyped_segments_vcfs_path_list

    Array[File] gcnv_model_tars = CohortWF.gcnv_model_tars
    Array[File] calling_configs = CohortWF.calling_configs
    Array[File] denoising_configs = CohortWF.denoising_configs
    Array[File] gcnvkernel_version = CohortWF.gcnvkernel_version
    Array[File] sharded_interval_lists = CohortWF.sharded_interval_lists

#    call JointCallWorkflow.JointCallExomeCNVs as JointCallWF {
#        input:
#             num_samples_per_scatter_block = num_samples_per_scatter_block,
#             intervals = intervals,
#             blacklist_intervals = blacklist_intervals,
#
#             contig_ploidy_calls_tar_path_list = CohortWF.contig_ploidy_calls_tar_path_list,
#             gcnv_calls_tars_path_list = CohortWF.gcnv_calls_tars_path_list,
#             genotyped_intervals_vcf_indexes_path_list = CohortWF.genotyped_intervals_vcf_indexes_path_list,
#             genotyped_intervals_vcfs_path_list = CohortWF.genotyped_intervals_vcfs_path_list,
#             genotyped_segments_vcf_indexes_path_list = CohortWF.genotyped_segments_vcf_indexes_path_list,
#             genotyped_segments_vcfs_path_list = CohortWF.genotyped_segments_vcfs_path_list,
#
#             maximum_number_events = maximum_number_events,
#             maximum_number_pass_events = maximum_number_pass_events,
#
#             gcnv_model_tars = CohortWF.gcnv_model_tars,
#             calling_configs = CohortWF.calling_configs,
#             denoising_configs = CohortWF.denoising_configs,
#             gcnvkernel_version = CohortWF.gcnvkernel_version,
#             sharded_interval_lists = CohortWF.sharded_interval_lists,
#             allosomal_contigs = allosomal_contigs,
#             ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
#             ref_fasta_dict = ref_fasta_dict,
#             ref_fasta_fai = ref_fasta_fai,
#             ref_fasta = ref_fasta,
#             x_contig_name = x_contig_name,
#             protein_coding_gtf = protein_coding_gtf,
#             linc_rna_gtf = linc_rna_gtf,
#             promoter_bed = promoter_bed,
#             noncoding_bed = noncoding_bed,
#             gatk_docker = gatk_docker,
#             gatk_docker_clustering = gatk_docker_clustering,
#             gatk_docker_qual_calc = gatk_docker_qual_calc,
#             sv_pipeline_docker = sv_pipeline_docker,
#
#    }





    #we do these as FoFNs for Terra compatibility
        Array[File] contig_ploidy_calls_tars = read_lines(contig_ploidy_calls_tar_path_list)
        Array[File] segments_vcfs = read_lines(genotyped_segments_vcfs_path_list)
        Array[File] segments_vcf_indexes = read_lines(genotyped_segments_vcf_indexes_path_list)
        Array[File] intervals_vcfs = read_lines(genotyped_intervals_vcfs_path_list)
        Array[File] intervals_vcf_indexes = read_lines(genotyped_intervals_vcf_indexes_path_list)
        Array[Array[File]] gcnv_calls_tars = read_tsv(gcnv_calls_tars_path_list)
        Array[Array[File]] call_tars_sample_by_shard = transpose(gcnv_calls_tars)

        #create a ped file to use for allosome copy number (e.g. XX, XY)
        call MakePedFile {
          input:
            contig_ploidy_calls_tar = read_lines(contig_ploidy_calls_tar_path_list),
            x_contig_name = x_contig_name
        }

        call CNVTasks.SplitInputArray as SplitSegmentsVcfsList {
            input:
                input_array = segments_vcfs,
                num_inputs_in_scatter_block = num_samples_per_scatter_block,
                gatk_docker = gatk_docker
        }

        call CNVTasks.SplitInputArray as SplitSegmentsIndexesList {
            input:
                input_array = segments_vcf_indexes,
                num_inputs_in_scatter_block = num_samples_per_scatter_block,
                gatk_docker = gatk_docker
        }

        Array[Array[String]] split_segments = SplitSegmentsVcfsList.split_array
        Array[Array[String]] split_segments_indexes = SplitSegmentsIndexesList.split_array

        #for more than num_samples_per_scatter_block, do an intermediate combine first
        if (length(split_segments) > 1) {
          scatter (subarray_index in range(length(split_segments))) {
            call JointSegmentation as ScatterJointSegmentation {
              input:
                segments_vcfs = split_segments[subarray_index],
                segments_vcf_indexes = split_segments_indexes[subarray_index],
                ped_file = MakePedFile.ped_file,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                gatk_docker = gatk_docker_clustering,
                model_intervals = intervals
            }
          }
        }

        #refine breakpoints over all samples
        call JointSegmentation as GatherJointSegmentation {
          input:
            segments_vcfs = select_first([ScatterJointSegmentation.clustered_vcf, segments_vcfs]),
            segments_vcf_indexes = select_first([ScatterJointSegmentation.clustered_vcf_index, segments_vcfs]),
            ped_file = MakePedFile.ped_file,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_docker = gatk_docker_clustering,
            model_intervals = intervals
        }

        #recalculate each sample's quality scores based on new breakpoints and filter low QS or high AF events;
        #exclude samples with too many events
        scatter (scatter_index in range(length(segments_vcfs))) {
          call CNVTasks.PostprocessGermlineCNVCalls as RecalcQual {
            input:
                  entity_id = sub(sub(basename(intervals_vcfs[scatter_index]), ".vcf.gz", ""), "intervals_output_", ""),
                  gcnv_calls_tars = call_tars_sample_by_shard[scatter_index],
                  gcnv_model_tars = gcnv_model_tars,
                  calling_configs = calling_configs,
                  denoising_configs = denoising_configs,
                  gcnvkernel_version = gcnvkernel_version,
                  sharded_interval_lists = sharded_interval_lists,
                  contig_ploidy_calls_tar = read_lines(contig_ploidy_calls_tar_path_list)[0],  #this is always a list of one tar
                  allosomal_contigs = allosomal_contigs,
                  ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                  sample_index = scatter_index,
                  maximum_number_events = maximum_number_events,
                  maximum_number_pass_events = maximum_number_pass_events,
                  intervals_vcf = intervals_vcfs[scatter_index],
                  intervals_vcf_index = intervals_vcf_indexes[scatter_index],
                  clustered_vcf = GatherJointSegmentation.clustered_vcf,
                  clustered_vcf_index = GatherJointSegmentation.clustered_vcf_index,
                  gatk_docker = gatk_docker_qual_calc
          }
        }

        #only put samples that passed QC into the combined VCF
        scatter(idx in range(length(RecalcQual.genotyped_segments_vcf))) {
          if (RecalcQual.qc_status_string[idx] == "PASS") {
            String subset = RecalcQual.genotyped_segments_vcf[idx]
            String subset_indexes = RecalcQual.genotyped_segments_vcf_index[idx]
          }
          if (RecalcQual.qc_status_string[idx] != "PASS") {
            String failed = sub(sub(basename(RecalcQual.genotyped_segments_vcf[idx]), ".vcf.gz", ""), "segments_output_", "")
          }
        }
        Array[String] subset_arr = select_all(subset)
        Array[String] subset_index_arr = select_all(subset_indexes)
        Array[String] failed_qc_samples = select_all(failed)

        call FastCombine {
          input:
            input_vcfs = subset_arr,
            input_vcf_indexes = subset_index_arr,
            sv_pipeline_docker = sv_pipeline_docker
        }






    output {
#        File contig_ploidy_model_tar = CohortWF.contig_ploidy_model_tar
#        File filtered_intervals = CohortWF.filtered_intervals
#        File clustered_vcf = JointSegmentation.clustered_vcf
#        File clustered_vcf_index = JointSegmentation.clustered_vcf_index
        File clustered_vcf = FastCombine.combined_vcf
        File clustered_vcf_index = FastCombine.combined_vcf_index
#        Array[File] gcnv_model_tars = CohortWF.gcnv_model_tar



#        File preprocessed_intervals = CohortWF.preprocessed_intervals
#        Array[File] read_counts_entity_ids = CohortWF.entity_id
#        Array[File] read_counts = CohortWF.counts
#        File? annotated_intervals = CohortWF.annotated_intervals
#        File filtered_intervals = CohortWF.filtered_intervals
#        File contig_ploidy_model_tar = CohortWF.contig_ploidy_model_tar
#        File contig_ploidy_calls_tar = CohortWF.contig_ploidy_calls_tar
#        File contig_ploidy_calls_tar_path_list = CohortWF.path_list
#        Array[File] sample_contig_ploidy_calls_tars = CohortWF.sample_contig_ploidy_calls_tar
#        Array[File] gcnv_model_tars = CohortWF.gcnv_model_tar
#        Array[Array[File]] gcnv_calls_tars = CohortWF.gcnv_call_tars
#        File gcnv_calls_tars_path_list = CohortWF.path_list
#        Array[File] gcnv_tracking_tars = CohortWF.gcnv_tracking_tar

#        Array[File] genotyped_intervals_vcfs = CohortWF.genotyped_intervals_vcf
#        File genotyped_intervals_vcfs_path_list = CohortWF.path_list
#        Array[File] genotyped_intervals_vcf_indexes = CohortWF.genotyped_intervals_vcf_index
#        File genotyped_intervals_vcf_indexes_path_list = CohortWF.path_list
#        Array[File] genotyped_segments_vcfs = CohortWF.genotyped_segments_vcf
#        File genotyped_segments_vcfs_path_list = CohortWF.path_list
#        Array[File] genotyped_segments_vcf_indexes = CohortWF.genotyped_segments_vcf_index
#        File genotyped_segments_vcf_indexes_path_list = CohortWF.path_list
#
#        Array[File] denoised_copy_ratios = CohortWF.denoised_copy_ratios
#        Array[File] sample_qc_status_files = CohortWF.qc_status_file
#        Array[String] sample_qc_status_strings = CohortWF.qc_status_string
#        File model_qc_status_file = CohortWF.qc_status_file
#        String model_qc_string = CohortWF.qc_status_string
#        Array[File] denoised_copy_ratios = CohortWF.denoised_copy_ratios
#
#        Array[File] gcnv_model_tars = CohortWF.gcnv_model_tar
#        Array[File] calling_configs = CohortWF.calling_config_json
#        Array[File] denoising_configs = CohortWF.denoising_config_json
#        Array[File] gcnvkernel_version = CohortWF.gcnvkernel_version_json
#        Array[File] sharded_interval_lists = CohortWF.sharded_interval_list

    }



}

task MakePedFile {
    input {
        Array[File] contig_ploidy_calls_tar
        String x_contig_name
    }

    command <<<
        set -e

        while read tar
        do
        mkdir callsDir
        tar -xf $tar -C callsDir

        for sample in $(ls -d -1 callsDir/SAMPLE*)
        do
        sample_name=$(cat $sample/sample_name.txt)
        x_ploidy=$(grep ^~{x_contig_name} $sample/contig_ploidy.tsv | cut -f 2)
        [[ -z "$x_ploidy" ]] && { echo "Chromosome ~{x_contig_name} ploidy call not found for sample " $sample_name; exit 1; }
        printf "%s\t%s\t0\t0\t%s\t0\n" $sample_name $sample_name $x_ploidy >> cohort.ped
        done
        rm -rf callsDir
        done < ~{write_lines(contig_ploidy_calls_tar)}
    >>>

    output {
        File ped_file = "cohort.ped"
    }

    runtime {
        docker: "gatksv/sv-base-mini:b3af2e3"
        memory: "3000 MB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 2
    }
}


task JointSegmentation {
    input {
        Array[File] segments_vcfs
        Array[File] segments_vcf_indexes
        File ped_file
        File ref_fasta_dict
        File ref_fasta_fai
        File ref_fasta
        File model_intervals

        # Runtime parameters
        String gatk_docker
        Int? mem_gb
        Int? disk_space_gb
        Boolean use_ssd = false
        Int? cpu
        Int? preemptible_attempts
    }

    parameter_meta {
        segments_vcfs: {localization_optional: true}
        segments_vcf_indexes: {localization_optional: true}
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    #NOTE: output has to be gzipped to be read in by pyvcf in the next step
    command <<<
        set -e
        gatk --java-options "-Xmx~{command_mem_mb}m" JointGermlineCNVSegmentation \
        -R ~{ref_fasta} -O clustered.vcf.gz -V ~{sep=' -V ' segments_vcfs} --model-call-intervals ~{model_intervals} -ped ~{ped_file}
    >>>

    output {
        File clustered_vcf = "clustered.vcf.gz"
        File clustered_vcf_index = "clustered.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 2])
    }
}

task FastCombine {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String sv_pipeline_docker
        Int? preemptible_tries
        Int? disk_size
        Float? mem_gb
    }

    command <<<
        #bcftools gets pissy if the indexes look older than their VCFs
        index_fofn=~{write_lines(input_vcf_indexes)}
        while read index; do touch $index; done < $index_fofn

        bcftools merge -l ~{write_lines(input_vcfs)} -o combined.vcf.gz -O z --threads 4 -m all -0

        tabix combined.vcf.gz
    >>>

    output {
        File combined_vcf = "combined.vcf.gz"
        File combined_vcf_index = "combined.vcf.gz.tbi"
    }

    runtime {
        docker: sv_pipeline_docker
        preemptible: select_first([preemptible_tries, 2])
        memory: select_first([mem_gb, 3.5]) + " GiB"
        cpu: "1"
        disks: "local-disk " + select_first([disk_size, 50]) + " HDD"
    }
}
