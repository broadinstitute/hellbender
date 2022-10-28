version 1.0

workflow GvsCalculatePrecisionAndSensitivity {
  input {
    Array[File] input_vcfs
    String output_basename

    String chromosome = "chr20"

    Array[String] sample_names
    Array[File] truth_vcfs
    Array[File] truth_vcf_indices
    Array[File] truth_beds

    File ref_fasta
  }

  parameter_meta {
    input_vcfs: "A collection of VCFS used for analysis (these are generated by `GvsExtractCallSet'). Note: these need not be subsetted down to the chromosome."
    output_basename: "The base name for the output files generated by the pipeline."
    chromosome: "The chromosome on which to run the analysis of Precision and Sensitivity. The default value for this is `chr20`. If it is set to `all` then the analysis will be run across *all* chromosomes."
    sample_names: "A list of the sample names that are controls and that will be used for the analysis. For every element on the list of sample names there must be a corresponding element on the list of `truth_vcfs`, `truth_vcf_indices`, and `truth_beds`."
    truth_vcfs: "A list of the VCFs that contain the truth data used for analyzing the samples in `sample_names`."
    truth_vcf_indices: "A list of the VCF indices for the truth data VCFs supplied above."
    truth_beds: "A list of the bed files for the truth data used for analyzing the samples in `sample_names`."
    ref_fasta: "The cloud path for the reference fasta sequence."
  }

  # WDL 1.0 trick to set a variable ('none') to be undefined.
  if (false) {
    String? none = "None"
  }

  # Couldn't figure out how to use a workflow input that is optional, but predefined (e.g. String? chromosome = "chr20"
  # So instead, the workflow takes as input the (non-optional) 'chromosome' argument, predefined to "chr20".
  # In order to have the analysis run on ALL chrosomses, set the input 'chromosome' to 'all', which then sets an internal variable 'contig' to be undefined.
  # This variable is then passed to tasks that use this value, where undefined indicates to analyze across all chromosomes.
  String? contig = if (chromosome == "all") then none else chromosome

  if ((length(sample_names) != length(truth_vcfs)) || (length(sample_names) != length(truth_vcf_indices)) || (length(sample_names) != length(truth_beds))) {
    call ErrorWithMessage {
      input:
        message = "The inputs 'sample_names', 'truth_vcfs', 'truth_vcf_indices', and 'truth_beds' must all contain the same number of elements"
    }
  }

  String output_chr_basename = if defined(contig) then output_basename + "." + contig else output_basename

  if (defined(contig)) {
    scatter(i in range(length(input_vcfs))) {
      call IsVcfOnChromosome {
        input:
          input_vcf = input_vcfs[i],
          chromosome = select_first([contig])
      }
    }
  }

  call GatherVcfs {
    input:
      input_vcfs = select_first([IsVcfOnChromosome.output_vcf, input_vcfs]),
      output_basename = output_chr_basename
  }

  scatter(i in range(length(sample_names))) {
    String sample_name = sample_names[i]
    String output_sample_basename = output_chr_basename + "." + sample_name

    call SelectVariants {
      input:
        input_vcf = GatherVcfs.output_vcf,
        input_vcf_index = GatherVcfs.output_vcf_index,
        contig = contig,
        sample_name = sample_name,
        output_basename = output_sample_basename
    }

    call Add_AS_MAX_VQSLOD_ToVcf {
      input:
        input_vcf = SelectVariants.output_vcf,
        output_basename = output_sample_basename + ".maxas"
    }

    call BgzipAndTabix {
      input:
        input_vcf = Add_AS_MAX_VQSLOD_ToVcf.output_vcf,
        output_basename = output_sample_basename + ".maxas"
    }

    call EvaluateVcf as EvaluateVcfFiltered {
      input:
        input_vcf = BgzipAndTabix.output_vcf,
        input_vcf_index = BgzipAndTabix.output_vcf_index,
        truth_vcf = truth_vcfs[i],
        truth_vcf_index = truth_vcf_indices[i],
        truth_bed = truth_beds[i],
        contig = contig,
        output_basename = sample_name + "-bq_roc_filtered",
        ref_fasta = ref_fasta
    }

    call EvaluateVcf as EvaluateVcfAll {
      input:
        input_vcf = BgzipAndTabix.output_vcf,
        input_vcf_index = BgzipAndTabix.output_vcf_index,
        truth_vcf = truth_vcfs[i],
        truth_vcf_index = truth_vcf_indices[i],
        truth_bed = truth_beds[i],
        contig = contig,
        all_records = true,
        output_basename = sample_name + "-bq_all",
        ref_fasta = ref_fasta
    }
  }

  call CollateReports {
    input:
      all_reports = EvaluateVcfAll.report,
      filtered_reports = EvaluateVcfFiltered.report
  }

  output {
    File report = CollateReports.report
    Array[Array[File]] filtered_eval_outputs = EvaluateVcfFiltered.outputs
    Array[Array[File]] all_eval_outputs = EvaluateVcfAll.outputs
  }
}

task IsVcfOnChromosome {
  input {
    File input_vcf
    String chromosome
  }

  String output_vcf_name = basename(input_vcf)

  command {
    set -e -o pipefail

    cat ~{input_vcf} | gunzip | grep -v '^#' | cut -f 1 | sort | uniq > chrom.txt
    NL=$(cat chrom.txt | wc -l)
    if [ $NL -ne 1 ]; then
      echo "~{input_vcf} has either no records or records on multiple chromosomes"
      exit 1
    fi

    mkdir output
    touch output/~{output_vcf_name}
    CHR=$(cat chrom.txt)
    if [ $CHR = ~{chromosome} ]; then
      cp ~{input_vcf} output/~{output_vcf_name}
      echo "Including ~{input_vcf} as it is on chromosome ~{chromosome}."
    else
      touch output/~{output_vcf_name}
      echo "Skipping ~{input_vcf} as it is not on chromosome ~{chromosome}."
    fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    File output_vcf = "output/~{output_vcf_name}"
  }
}

task GatherVcfs {
  input {
    Array[File] input_vcfs
    String output_basename

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(3*size(input_vcfs, "GiB"))
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    set -euo pipefail

    CHR_VCFS_ARG=""
    for file in ~{sep=' ' input_vcfs}
    do
      if [ -s $file ]; then
        CHR_VCFS_ARG+=" --INPUT $file "
      fi
    done
    echo $CHR_VCFS_ARG

    # --REORDER_INPUT_BY_FIRST_VARIANT means that the vcfs supplied here need not be ordered by location.
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      GatherVcfs \
        --REORDER_INPUT_BY_FIRST_VARIANT \
        $CHR_VCFS_ARG \
        --OUTPUT ~{output_basename}.vcf.gz

    tabix ~{output_basename}.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task SelectVariants {
  input {
    File input_vcf
    File input_vcf_index
    String? contig
    String sample_name

    String output_basename

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 50
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      SelectVariants \
        -V ~{input_vcf} \
        ~{"-L " + contig} \
        --sample-name ~{sample_name} \
        --select-type-to-exclude NO_VARIATION \
        -O ~{output_basename}.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task Add_AS_MAX_VQSLOD_ToVcf {
  input {
    File input_vcf
    String output_basename

    String docker = "us.gcr.io/broad-dsde-methods/variantstore:2022-10-28-alpine"
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 50
  }

  command <<<
    set -e

    python3 /app/add_max_as_vqslod.py ~{input_vcf} > ~{output_basename}.vcf
  >>>
  runtime {
    docker: docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File output_vcf = "~{output_basename}.vcf"
  }
}

task BgzipAndTabix {
  input {
    File input_vcf
    String output_basename

    String docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.5-1.10.2-0.1.16-1649948623"
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(3 * size(input_vcf, "GiB")) + 50
  }

  command {
    # note that bgzip has an option (-i) to index the bgzipped output, but this file is not a tabix file
    # note also that we use '-c' so that bgzip doesn't create the bgzipped file in place, rather it's in a location
    # where it's easy to output from the task.
    bgzip -c ~{input_vcf} > ~{output_basename}.vcf.gz
    tabix ~{output_basename}.vcf.gz
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"

  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task EvaluateVcf {
  input {
    File input_vcf
    File input_vcf_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed

    String? contig
    Boolean all_records = false

    File ref_fasta

    String output_basename

    String docker = "docker.io/realtimegenomics/rtg-tools:latest"
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(2 * size(ref_fasta, "GiB")) + 50
  }

  command <<<
    set -e -o pipefail

    rtg format --output human_REF_SDF ~{ref_fasta}

    rtg vcfeval \
      ~{"--region " + contig} \
      ~{if all_records then "--all-records" else ""} \
      --roc-subset snp,indel \
      --vcf-score-field=INFO.MAX_AS_VQSLOD \
      -t human_REF_SDF \
      -b ~{truth_vcf} \
      -e ~{truth_bed}\
      -c ~{input_vcf} \
      -o ~{output_basename}

    # Touch a file with the name of the sample in that directory, so that it's identifiable among the globbed outputs.
    touch ~{output_basename}/~{output_basename}

    touch report.txt
    for type in "snp" "indel"
      do
        d=$(cat ~{output_basename}/${type}_roc.tsv.gz | gunzip | tail -1 | cut -f3,5,6,7)
        echo -e "~{output_basename}\t$type\t$d" >> report.txt
    done
  >>>
  runtime {
    docker: docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"

  }
  output {
    File report = "report.txt"
    Array[File] outputs = glob("~{output_basename}/*")
  }
}

task CollateReports {
  input {
    Array[File] all_reports
    Array[File] filtered_reports
  }

  command {
    set -e -o pipefail

    echo "sample  type  FPs FNs precision sensitivity"
    while read -r a;
    do
      cat $a
    done < ~{write_lines(all_reports)}

    while read -r a;
    do
      cat $a
    done < ~{write_lines(filtered_reports)}
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    File report = stdout()
  }
}

# Print given message to stderr and return an error
task ErrorWithMessage{
  input {
    String message
  }
  command <<<
    >&2 echo "Error: ~{message}"
    exit 1
  >>>

  runtime {
    docker: "ubuntu:20.04"
  }
}
