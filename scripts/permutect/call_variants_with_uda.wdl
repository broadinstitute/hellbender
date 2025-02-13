version 1.0

# run Mutect2 to get both training AND test datasets.  The training dataset is preprocessed and combined with
# high-quality labeled data to make a UDA dataset, then used to train an artifact model.  The test dataset is used
# for the posterior model and filtering.

# note that the artifact model can be trained before the Mutect2 workflow runs FilterMutectCalls

import "https://api.firecloud.org/ga4gh/v1/tools/davidben:mutect2/versions/18/plain-WDL/descriptor" as m2
import "permutect-uda-dataset.wdl" as uda
import "permutect-train-artifact-model.wdl" as training
import "permutect-call-variants.wdl" as calling

workflow CallVariantsWithUDA {
    input {
        # basic inputs for Mutect2
        File? intervals
        File? masked_intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File primary_bam
        File primary_bai
        File? control_bam
        File? control_bai
        File? gnomad
        File? gnomad_idx
        String? m2_extra_args
        File? dragstr_model
        Boolean make_bamout = false
        Boolean compress_vcfs = false

        # Mutect2 filtering
        Boolean skip_m2_filtering
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? realignment_index_bundle
        String? realignment_extra_args
        Boolean? run_orientation_bias_mixture_model_filter

        # preprocessing arguments
        Int chunk_size

        # training arguments for both artifact model and posterior model
        Int batch_size
        Int inference_batch_size
        Int num_workers
        Int? gpu_count
        Int? training_mem

        # UDA training arguments
        File base_model
        File source_train_tar
        String source_edit_type = "keep_everything"
        String target_edit_type = "unlabel_everything"
        Int num_epochs
        Int num_calibration_epochs
        Float dropout_p
        Array[Int] aggregation_layers
        Array[Int] calibration_layers
        String? training_extra_args
        Boolean learn_artifact_spectra
        Float? genomic_span

        # Permutect filtering / posterior model
        File? test_dataset_truth_vcf    # used for evaluation
        File? test_dataset_truth_vcf_idx
        Int? num_spectrum_iterations
        Float? spectrum_learning_rate
        String? permutect_filtering_extra_args
        String bcftools_docker = "us.gcr.io/broad-dsde-methods/davidben/bcftools"
        File? obscene_hack_leave_unset


        # runtime
        String gatk_docker
        String permutect_docker
        File? gatk_override
        String basic_bash_docker = "ubuntu:16.04"
        Int scatter_count
        Int preemptible = 2
        Int max_retries = 1
        Int small_task_cpu = 2
        Int small_task_mem = 4
        Int small_task_disk = 100
        Int boot_disk_size = 12
        Int learn_read_orientation_mem = 8000
        Int filter_alignment_artifacts_mem = 9000
        String? gcs_project_for_requester_pays

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_disk = 0
    }

    # note: we make both training and test datasets
    # note: for speed we may skip filtering in order to begin UDA artifact model training immediately
    # the only M2 filtering we may need is contamination, and that may be skipped
    call m2.Mutect2 {
        input:
            intervals = intervals,
            masked_intervals = masked_intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            tumor_reads = primary_bam,
            tumor_reads_index = primary_bai,
            normal_reads = control_bam,
            normal_reads_index = control_bai,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            realignment_index_bundle = realignment_index_bundle,
            realignment_extra_args = realignment_extra_args,
            run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
            m2_extra_args = m2_extra_args,
            dragstr_model = dragstr_model,
            make_bamout = make_bamout,
            make_permutect_training_dataset = true,
            make_permutect_test_dataset = true,
            permutect_test_dataset_truth_vcf = test_dataset_truth_vcf,
            permutect_test_dataset_truth_vcf_idx = test_dataset_truth_vcf_idx,
            skip_filtering = skip_m2_filtering,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            scatter_count = scatter_count,
            preemptible = preemptible,
            max_retries =  max_retries,
            small_task_cpu = small_task_cpu,
            small_task_mem = small_task_mem,
            small_task_disk = small_task_disk,
            boot_disk_size = boot_disk_size,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays,
            emergency_extra_disk = emergency_extra_disk
    }

    # preprocess the training data from Mutect2
    call Preprocess {
        input:
            training_dataset = select_first([Mutect2.permutect_training_dataset]),
            chunk_size = chunk_size,
            permutect_docker = permutect_docker
    }

    # combine the source_tar and preprocessed training data into a UDA dataset
    call uda.PermutectUDADataset {
        input:
            source_train_tar = source_train_tar,
            target_train_tar = Preprocess.train_tar,
            source_edit_type = source_edit_type,
            target_edit_type = target_edit_type,
            chunk_size = chunk_size,
            permutect_docker = permutect_docker,
            preemptible = 0,
            max_retries = 0
    }

    # train an artifact model on the UDA dataset
    call training.TrainPermutect {
        input:
            train_tar = PermutectUDADataset.uda_train_tar,
            base_model = base_model,
            num_epochs = num_epochs,
            num_calibration_epochs = num_calibration_epochs,
            batch_size = batch_size,
            inference_batch_size = inference_batch_size,
            num_workers = num_workers,
            mem = training_mem,
            gpu_count = gpu_count,
            dropout_p = dropout_p,
            aggregation_layers = aggregation_layers,
            calibration_layers = calibration_layers,
            extra_args = training_extra_args,
            learn_artifact_spectra = learn_artifact_spectra,
            genomic_span = genomic_span,
            permutect_docker = permutect_docker,
            preemptible = 0,
            max_retries = 0
    }

    # we already ran M2 so we don't need the entire calling workflow, just the post-M2 parts of it
    call calling.SplitMultiallelics {
        input:
            input_vcf = Mutect2.output_vcf,
            input_vcf_idx = Mutect2.output_vcf_idx,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            bcftools_docker = bcftools_docker
    }

    call calling.IndexVCF as IndexAfterSplitting {
        input:
            unindexed_vcf = SplitMultiallelics.output_vcf,
            gatk_docker = gatk_docker
    }

    call calling.PermutectFiltering {
        input:
            mutect2_vcf = IndexAfterSplitting.vcf,
            mutect2_vcf_idx = IndexAfterSplitting.vcf_index,
            permutect_model = TrainPermutect.artifact_model,
            test_dataset = select_first([Mutect2.permutect_test_dataset]),
            contigs_table = Mutect2.permutect_contigs_table,
            maf_segments = Mutect2.maf_segments,
            mutect_stats = Mutect2.mutect_stats,
            batch_size = batch_size,
            num_workers = num_workers,
            gpu_count = gpu_count,
            num_spectrum_iterations = num_spectrum_iterations,
            spectrum_learning_rate = spectrum_learning_rate,
            chunk_size = chunk_size,
            permutect_filtering_extra_args = permutect_filtering_extra_args,
            permutect_docker = permutect_docker,
    }


    call calling.IndexVCF as IndexAfterFiltering {
        input:
            unindexed_vcf = PermutectFiltering.output_vcf,
            gatk_docker = gatk_docker
    }

    output {
        File? bamout = Mutect2.bamout
        File? bamout_index = Mutect2.bamout_index
        File mutect_stats = Mutect2.mutect_stats
        File permutect_contigs_table = Mutect2.permutect_contigs_table
        File permutect_read_groups_table = Mutect2.permutect_read_groups_table
        File train_tar = Preprocess.train_tar
        File training_tensorboard_tar = TrainPermutect.training_tensorboard_tar
        File output_vcf = IndexAfterFiltering.vcf
        File output_vcf_idx = IndexAfterFiltering.vcf_index
        File calling_tensorboard_tar = PermutectFiltering.tensorboard_report
    }

}

task Preprocess {
    input {
        File training_dataset
        Int chunk_size
        Int? source_label

        String permutect_docker
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Int? mem
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 16000
    Int command_mem = machine_mem - 500

    command <<<
        set -e

        gatk PermutectPreprocessDataset --training-datasets ~{training_dataset} --chunk-size ~{chunk_size} ~{"--sources " + source_label} --output train.tar
    >>>

    runtime {
        docker: permutect_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 2])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File train_tar = "train.tar"
    }
}