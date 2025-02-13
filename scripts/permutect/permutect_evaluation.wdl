version 1.0

workflow PermutectEvaluation {
    input {
        File permutect_model
        File evaluation_tar
        Int batch_size
        Int? num_workers

        String permutect_docker
        Int? preemptible
        Int? max_retries
    }


   call Evaluate {
        input:
            evaluation_tar = evaluation_tar,
            permutect_model = permutect_model,
            batch_size = batch_size,
            num_workers = num_workers,
            permutect_docker = permutect_docker
    }


    output {
        File tensorboard_tar = Evaluate.tensorboard_tar
    }
}

task Evaluate {
    input {
        File evaluation_tar
        File permutect_model
        Int batch_size
        Int? num_workers
        String? extra_args

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

        evaluate_model \
            --evaluation_tar ~{evaluation_tar} \
            --permutect_model ~{permutect_model} \
            --batch_size ~{batch_size} \
            ~{"--num_workers " + num_workers} \
            --tensorboard_dir tensorboard \
            ~{extra_args}

        tar cvf tensorboard.tar tensorboard/
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
        File tensorboard_tar = "tensorboard.tar"
    }
}
