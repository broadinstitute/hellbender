version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/davidben:mutect2/versions/17/plain-WDL/descriptor" as m2

workflow Permutect {
    input {
        File permutect_model

        File? intervals
        File? masks
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count
        Int? num_spectrum_iterations
        Float? spectrum_learning_rate
        File primary_bam
        File primary_bai
        File? control_bam
        File? control_bai
        File? gnomad
        File? gnomad_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? realignment_index_bundle
        File? dragstr_model
        String? realignment_extra_args
        Boolean? run_orientation_bias_mixture_model_filter
        String? m2_extra_args
        String? split_intervals_extra_args
        Int batch_size
        Int num_workers
        Int? gpu_count
        Int chunk_size
        File? test_dataset_truth_vcf    # used for evaluation
        File? test_dataset_truth_vcf_idx

        String? permutect_filtering_extra_args
        String gatk_docker
        String bcftools_docker
        File? gatk_override
        String permutect_docker
        Int? preemptible
        Int? max_retries
        File? obscene_hack_leave_unset
    }

    call m2.Mutect2 {
        input:
            make_permutect_training_dataset = false,
            make_permutect_test_dataset = true,
            permutect_training_dataset_truth_vcf = test_dataset_truth_vcf,
            permutect_training_dataset_truth_vcf_idx = test_dataset_truth_vcf_idx,
            intervals = intervals,
            masked_intervals = masks,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            tumor_reads = primary_bam,
            tumor_reads_index = primary_bai,
            normal_reads = if control_bam == "" then obscene_hack_leave_unset else control_bam,
            normal_reads_index = if control_bam == "" then obscene_hack_leave_unset else control_bai,

            scatter_count = scatter_count,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            realignment_index_bundle = realignment_index_bundle,
            realignment_extra_args = realignment_extra_args,
            dragstr_model = dragstr_model,
            run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
            m2_extra_args = m2_extra_args,
            make_bamout = false,

            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries
    }

    call SplitMultiallelics {
        input:
            input_vcf = Mutect2.output_vcf,
            input_vcf_idx = Mutect2.output_vcf_idx,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            bcftools_docker = bcftools_docker
    }

    call IndexVCF as IndexAfterSplitting {
        input:
            unindexed_vcf = SplitMultiallelics.output_vcf,
            gatk_docker = gatk_docker
    }

    call PermutectFiltering {
        input:
            mutect2_vcf = IndexAfterSplitting.vcf,
            mutect2_vcf_idx = IndexAfterSplitting.vcf_index,
            permutect_model = permutect_model,
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

    call IndexVCF as IndexAfterFiltering {
        input:
            unindexed_vcf = PermutectFiltering.output_vcf,
            gatk_docker = gatk_docker
    }

    output {
        File output_vcf = IndexAfterFiltering.vcf
        File output_vcf_idx = IndexAfterFiltering.vcf_index
        File tensorboard_report = PermutectFiltering.tensorboard_report
        File test_dataset = select_first([Mutect2.permutect_test_dataset])
        File mutect2_vcf = Mutect2.output_vcf
        File mutect2_vcf_idx = Mutect2.output_vcf_idx
    }
}

 task PermutectFiltering {
    input {
        File permutect_model
        File test_dataset
        File contigs_table
        File mutect2_vcf
        File mutect2_vcf_idx
        File? maf_segments
        File? normal_maf_segments
        File mutect_stats
        Int? num_spectrum_iterations
        Float? spectrum_learning_rate
        Int batch_size
        Int num_workers
        Int? gpu_count
        Int chunk_size
        String? permutect_filtering_extra_args

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
        # set -e
        genomic_span=`grep "callable" ~{mutect_stats} | while read name value; do echo $value; done`

        gatk PermutectFilterVariants --input ~{mutect2_vcf} --test-dataset ~{test_dataset} \
            --permutect-model ~{permutect_model} \
            --contigs-table ~{contigs_table} \
            --output permutect-filtered.vcf \
            --tensorboard-dir tensorboard \
            --batch-size ~{batch_size} --num-workers ~{num_workers} --chunk-size ~{chunk_size} \
            ~{" --num-spectrum-iterations " + num_spectrum_iterations} \
            ~{" --spectrum-learning-rate " + spectrum_learning_rate} \
            ~{" --maf-segments " + maf_segments} ~{" --normal-maf-segments " + normal_maf_segments} \
            --genomic-span $genomic_span ~{permutect_filtering_extra_args}

        tar cvf tensorboard.tar tensorboard/
    >>>

    runtime {
        docker: permutect_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 0])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 2])
        gpuType: "nvidia-tesla-t4"
        gpuCount: select_first([gpu_count, 1])
        nvidiaDriverVersion: "535.183.01"
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
    }

    output {
        File output_vcf = "permutect-filtered.vcf"
        File tensorboard_report = "tensorboard.tar"
    }
}

task SplitMultiallelics {
    input {
        File input_vcf
        File input_vcf_idx
        File ref_fasta
        File ref_fai
        File ref_dict
        String bcftools_docker
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Int? mem
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 4000
    Int command_mem = machine_mem - 500

    command <<<

        bcftools norm -m -any -f ~{ref_fasta} ~{input_vcf} > output.vcf

    >>>

    runtime {
        docker: bcftools_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 0])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 2])
    }

    output {
        File output_vcf = "output.vcf"
    }
}

task IndexVCF {
    input {
        File unindexed_vcf
        String gatk_docker
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Int? mem
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 4000
    Int command_mem = machine_mem - 500

    command <<<

        cp ~{unindexed_vcf} indexed.vcf

        gatk --java-options "-Xmx~{command_mem}m" IndexFeatureFile -I indexed.vcf

        gatk --java-options "-Xmx~{command_mem}m" SelectVariants -V indexed.vcf -O output.vcf --lenient \
            -DGA DP -DGA AF -DGA F1R2 -DGA F2R1 -DGA FAD -DGA SB \
            -DA AS_FilterStatus -DA AS_SB_TABLE -DA ECNT -DA GERMQ -DA MBQ -DA MFRL -DA MMQ -DA MPOS

        set -e
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 0])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File vcf = "output.vcf"
        File vcf_index = "output.vcf.idx"
    }
}
