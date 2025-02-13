version 1.0

# source and target train tar have already been preprocessed
# source needs to be given the '0' source label
# target needs to be given the '1' source label
workflow PermutectUDADataset {
    input {
        File source_train_tar
        File target_train_tar

        # most likely KEEP_EVERYTHING for the source and UNLABEL_ARTIFACTS or UNLABEL_EVERYTHING for target
        String source_edit_type
        String target_edit_type

        Int chunk_size

        String permutect_docker
        Int? preemptible
        Int? max_retries
    }

    call EditDataset as EditSource {
        input:
            train_tar = [source_train_tar],
            new_source = 0,
            edit_type = source_edit_type,
            chunk_size = chunk_size,
            permutect_docker = permutect_docker
    }

    call EditDataset as EditTarget {
        input:
            train_tar = [target_train_tar],
            new_source = 1,
            edit_type = target_edit_type,
            chunk_size = chunk_size,
            permutect_docker = permutect_docker
    }

    call EditDataset as Merge {
        input:
            train_tar = [EditSource.output_tarfile, EditTarget.output_tarfile],
            edit_type = "keep_everything",
            chunk_size = chunk_size,
            permutect_docker = permutect_docker
    }

    output {
        File uda_train_tar = Merge.output_tarfile
    }
}

task EditDataset {
    input {
        Array[File] train_tar
        Int chunk_size
        Int? new_source
        String edit_type
        String? extra_args

        String permutect_docker
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Int? mem
        Boolean use_ssd = false
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 16000
    Int command_mem = machine_mem - 500

    command <<<
        set -e

        gatk PermutectEditDataset \
            --train-tar ~{sep=' --train-tar ' train_tar} \
            --chunk-size ~{chunk_size} \
            ~{" --source " + new_source} \
            --dataset-edit ~{edit_type} \
            --output edited.tar \
            ~{extra_args}
    >>>

    runtime {
        docker: permutect_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File output_tarfile = "edited.tar"
    }
}