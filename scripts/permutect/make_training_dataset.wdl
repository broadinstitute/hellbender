version 1.0

# run Mutect2 without filtering to get plain text training data, then run preprocess_dataset.

import "https://api.firecloud.org/ga4gh/v1/tools/davidben:mutect2/versions/18/plain-WDL/descriptor" as m2

workflow MakeTrainingDataset {
    input {
        # basic inputs
        File? intervals
        File? masked_intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File reads
        File reads_index

        File? gnomad
        File? gnomad_idx

        # extra arguments
        String? m2_extra_args

        # preprocessing arguments
        Int chunk_size

        # additional modes and outputs
        File? dragstr_model
        Boolean make_bamout = false
        Boolean compress_vcfs = false
        File? permutect_training_dataset_truth_vcf
        File? permutect_training_dataset_truth_vcf_idx


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

    call m2.Mutect2 {
        input:
            intervals = intervals,
            masked_intervals = masked_intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            tumor_reads = reads,
            tumor_reads_index = reads_index,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx,
            m2_extra_args = m2_extra_args,
            dragstr_model = dragstr_model,
            make_bamout = make_bamout,
            make_permutect_training_dataset = true,
            permutect_training_dataset_truth_vcf = permutect_training_dataset_truth_vcf,
            permutect_training_dataset_truth_vcf_idx = permutect_training_dataset_truth_vcf_idx,
            skip_filtering = true,
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

    call Preprocess {
        input:
            training_dataset = select_first([Mutect2.permutect_training_dataset]),
            chunk_size = chunk_size,
            permutect_docker = permutect_docker
    }

    output {
        File? bamout = Mutect2.bamout
        File? bamout_index = Mutect2.bamout_index
        File mutect_stats = Mutect2.mutect_stats
        File permutect_contigs_table = Mutect2.permutect_contigs_table
        File permutect_read_groups_table = Mutect2.permutect_read_groups_table
        File train_tar = Preprocess.train_tar
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
        Boolean use_ssd = true
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
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 2])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File train_tar = "train.tar"
    }
}