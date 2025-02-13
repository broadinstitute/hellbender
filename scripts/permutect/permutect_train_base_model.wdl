version 1.0


workflow TrainPermutectBaseModel {
    input {
        File train_tar
        File? pretrained_model
        Int num_epochs
        Int batch_size
        Int inference_batch_size
        Int? num_workers
        Float dropout_p
        Float reweighting_range
        Array[Int] read_layers
        Int self_attention_hidden_dimension
        Int num_self_attention_layers
        Array[Int] info_layers
        Array[Int] aggregation_layers
        Array[String] ref_seq_layer_strings
        String? extra_args
        Int? gpu_count

        String permutect_docker
        Int? preemptible
        Int? max_retries
    }

    call TrainPermutectBase {
        input:
            train_tar = train_tar,
            pretrained_model = pretrained_model,
            permutect_docker = permutect_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            num_epochs = num_epochs,
            batch_size = batch_size,
            inference_batch_size = inference_batch_size,
            num_workers = num_workers,
            gpu_count = gpu_count,
            dropout_p = dropout_p,
            reweighting_range = reweighting_range,
            read_layers = read_layers,
            self_attention_hidden_dimension = self_attention_hidden_dimension,
            num_self_attention_layers = num_self_attention_layers,
            info_layers = info_layers,
            aggregation_layers = aggregation_layers,
            ref_seq_layer_strings = ref_seq_layer_strings,
            extra_args = extra_args
    }

    output {
        File base_model = TrainPermutectBase.base_model
        File training_tensorboard_tar = TrainPermutectBase.tensorboard_tar
    }
}


task TrainPermutectBase {
    input {
        File train_tar
        File? pretrained_model

        Int num_epochs
        Int batch_size
        Int inference_batch_size
        Int? num_workers
        Int? gpu_count
        Float dropout_p
        Float reweighting_range
        Array[Int] read_layers
        Int self_attention_hidden_dimension
        Int num_self_attention_layers
        Array[Int] info_layers
        Array[Int] aggregation_layers
        Array[String] ref_seq_layer_strings

        String? extra_args

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

        gatk PermutectTrainBaseModel \
            --train-tar ~{train_tar} \
            ~{"--pretrained-model " + pretrained_model} \
            --read-layers ~{sep=' --read-layers ' read_layers} \
            --self-attention-hidden-dimension ~{self_attention_hidden_dimension} \
            --num-self-attention-layers ~{num_self_attention_layers} \
            --info-layers ~{sep=' --info-layers ' info_layers} \
            --aggregation-layers ~{sep=' --aggregation-layers ' aggregation_layers} \
            --ref-seq-layer-strings ~{sep=' --ref-seq-layer-strings ' ref_seq_layer_strings} \
            --dropout-p ~{dropout_p} \
            --reweighting-range ~{reweighting_range} \
            --batch-size ~{batch_size} \
            --inference-batch-size ~{inference_batch_size} \
            ~{"--num-workers " + num_workers} \
            --num-epochs ~{num_epochs} \
            --output base_model.pt \
            --tensorboard-dir tensorboard \
            ~{extra_args}

        tar cvf tensorboard.tar tensorboard/
    >>>

    runtime {
        docker: permutect_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 0])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
        gpuType: "nvidia-tesla-t4"
        gpuCount: select_first([gpu_count, 1])
        nvidiaDriverVersion: "535.183.01"
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
    }

    output {
        File base_model = "base_model.pt"
        File tensorboard_tar = "tensorboard.tar"
    }
}

