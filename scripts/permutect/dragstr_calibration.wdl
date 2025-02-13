version 1.0

## NOTE: this is essentially copied from https://github.com/broadinstitute/warp/blob/develop/tasks/broad/DragenTasks.wdl
## with minor modifications

workflow DragstrCalibration {
    input {
        File ref_fasta
        File ref_fai
        File ref_dict
        File reads
        File reads_index
    }

    call CalibrateDragstrModel {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            reads = reads,
            reads_index = reads_index
    }

    output {
        File dragstr_model = CalibrateDragstrModel.dragstr_model
    }
}

task CalibrateDragstrModel {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File reads
    File reads_index

    String docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int preemptible_tries = 1
    Int threads = 4
    Int? disk_space
    Int? mem
  }

  String parallel_args  = "--threads " + threads

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 4000
  Int command_mem = machine_mem - 500

  parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      reads: {localization_optional: true}
      reads_index: {localization_optional: true}
    }

  command <<<
    set -x

    gatk ComposeSTRTableFile \
      -R ~{ref_fasta} \
      -O str_table.tsv


    gatk --java-options "-Xmx~{command_mem}m" CalibrateDragstrModel \
        -R ~{ref_fasta} \
        -I ~{reads} \
        -str str_table.tsv \
        -O params.dragstr \
        ~{parallel_args}
  >>>

  runtime {
     docker: docker
     disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
     memory: machine_mem + " MB"
     preemptible: preemptible_tries
     cpu: threads
  }

  output {
    File dragstr_model = "params.dragstr"
  }
}