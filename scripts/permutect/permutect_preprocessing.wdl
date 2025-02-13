version 1.0


workflow PreprocessPermutect {
    input {
        File training_dataset
        Int chunk_size
        Int? source_label

        String permutect_docker
        Int? preemptible
        Int? max_retries
    }

    call Preprocess {
        input:
            training_dataset = training_dataset,
            permutect_docker = permutect_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            chunk_size = chunk_size,
            source_label = source_label
    }

    output {
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
        Boolean use_ssd = false
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 16000
    Int command_mem = machine_mem - 500

    command <<<
        set -e

        preprocess_dataset --training_datasets ~{training_dataset} --chunk_size ~{chunk_size} ~{"--sources " + source_label} --output train.tar
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