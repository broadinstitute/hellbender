version 1.0
 
 workflow VariantCallingCarrot {

  String pipeline_version = "2.0.3"

  input {
    Boolean run_dragen_mode_variant_calling = false
    Boolean use_spanning_event_genotyping = true
    File calling_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Float? contamination = 0
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_str
    String base_file_name
    String final_vcf_base_name
    Int agg_preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gatk_control_docker = "broadinstitute/gatk-nightly:latest"
    Boolean make_gvcf = true
    Boolean make_bamout = false
    Boolean use_gatk3_haplotype_caller = false
    Boolean skip_reblocking = false
    Boolean use_dragen_hard_filtering = false
    File monitoring_script = "gs://emeryj-testing/cromwell_monitoring_script.sh"
  }
 
  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call ScatterIntervalList as ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }
  
  
    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (i in range(length(ScatterIntervalList.out))) {
  
      
      # Generate GVCF by interval
      call HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
        input:
          contamination = if run_dragen_mode_variant_calling then 0 else contamination,
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = ScatterIntervalList.out[i],
          vcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          hc_scatter = hc_divisor,
          make_gvcf = make_gvcf,
          make_bamout = make_bamout,
          run_dragen_mode_variant_calling = run_dragen_mode_variant_calling,
          use_dragen_hard_filtering = use_dragen_hard_filtering,
          use_spanning_event_genotyping = use_spanning_event_genotyping,
          preemptible_tries = agg_preemptible_tries,
          gatk_docker = gatk_docker,
          monitoring_script = monitoring_script,
          scatter_i = i
       }

      # Generate GVCF by interval (control)
      call HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4CONTROL {
        input:
          contamination = if run_dragen_mode_variant_calling then 0 else contamination,
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = ScatterIntervalList.out[i],
          vcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          hc_scatter = hc_divisor,
          make_gvcf = make_gvcf,
          make_bamout = make_bamout,
          run_dragen_mode_variant_calling = run_dragen_mode_variant_calling,
          use_dragen_hard_filtering = use_dragen_hard_filtering,
          use_spanning_event_genotyping = use_spanning_event_genotyping,
          preemptible_tries = agg_preemptible_tries,
          gatk_docker = gatk_control_docker,
          monitoring_script = monitoring_script,
          scatter_i = i
      }

      if (use_dragen_hard_filtering) {
        call DragenHardFilterVcf as DragenHardFilterVcf {
          input:
            input_vcf = HaplotypeCallerGATK4.output_vcf,
            input_vcf_index = HaplotypeCallerGATK4.output_vcf_index,
            make_gvcf = make_gvcf,
            vcf_basename = base_file_name,
            preemptible_tries = agg_preemptible_tries
        }
        call DragenHardFilterVcf as DragenHardFilterVcfCONTROL {
          input:
            input_vcf = HaplotypeCallerGATK4CONTROL.output_vcf,
            input_vcf_index = HaplotypeCallerGATK4CONTROL.output_vcf_index,
            make_gvcf = make_gvcf,
            vcf_basename = base_file_name,
            preemptible_tries = agg_preemptible_tries
        }
      }
    
    File vcfs_to_merge = select_first([DragenHardFilterVcf.output_vcf, HaplotypeCallerGATK4.output_vcf])
    File vcf_indices_to_merge = select_first([DragenHardFilterVcf.output_vcf_index, HaplotypeCallerGATK4.output_vcf_index])

    File vcfsCONTROL_to_merge = select_first([DragenHardFilterVcfCONTROL.output_vcf, HaplotypeCallerGATK4CONTROL.output_vcf])
    File vcfCONTROL_indices_to_merge = select_first([DragenHardFilterVcfCONTROL.output_vcf_index, HaplotypeCallerGATK4CONTROL.output_vcf_index])
  }


  # Combine by-interval (g)VCFs into a single sample (g)VCF file
  String hard_filter_suffix = if use_dragen_hard_filtering then ".hard-filtered" else ""
  String merge_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  call MergeVCFs as MergeVCFs {
    input:
      input_vcfs = vcfs_to_merge,
      input_vcfs_indexes = vcf_indices_to_merge,
      output_vcf_name = final_vcf_base_name + hard_filter_suffix + merge_suffix,
      preemptible_tries = agg_preemptible_tries
  }
  call MergeVCFs as MergeVCFsCONTROL {
    input:
      input_vcfs = vcfsCONTROL_to_merge,
      input_vcfs_indexes = vcfCONTROL_indices_to_merge,
      output_vcf_name = final_vcf_base_name + hard_filter_suffix + merge_suffix,
      preemptible_tries = agg_preemptible_tries
  }
  
    output {
    File output_vcf = select_first([MergeVCFs.output_vcf])
    File output_vcf_index = select_first([MergeVCFs.output_vcf_index])
    Array[File] output_runtimes = HaplotypeCallerGATK4.hc_time_out
    File representative_benchmarking = HaplotypeCallerGATK4.monitoring[39]

    File control_vcf = select_first([MergeVCFsCONTROL.output_vcf])
    File control_vcf_index = select_first([MergeVCFsCONTROL.output_vcf_index])
    Array[File] control_output_runtimes = HaplotypeCallerGATK4CONTROL.hc_time_out
    File control_representative_benchmarking = HaplotypeCallerGATK4CONTROL.monitoring[39]
  }
  meta {
    allowNestedInputs: true
  }
}


task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File input_bam_index
    File interval_list
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    Int preemptible_tries
    Int hc_scatter
    Boolean run_dragen_mode_variant_calling = false
    Boolean use_dragen_hard_filtering = false
    Boolean use_spanning_event_genotyping = true
    File? dragstr_model
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.2.0"
    Int memory_multiplier = 1
    File monitoring_script
    Int scatter_i
  }
  
  Int memory_size_mb = ceil(8000 * memory_multiplier)

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_file_name = vcf_basename + output_suffix

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 20

  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command <<<
    set -e

    bash ~{monitoring_script} > monitoring.log &

    # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
    # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
    # memory_size_gb because of Cromwell's retry with more memory feature.
    # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
    #       which do not rely on the output format of the `free` command.
    available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
    let java_memory_size_mb=available_memory_mb-2024
    echo Total available memory: ${available_memory_mb} MB >&2
    echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
    apt update && apt install time
    /usr/bin/time -f $'%e\n%U\n%S' -o runtime-~{scatter_i}.txt gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_file_name} \
      -contamination ~{default=0 contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      ~{true="--dragen-mode" false="" run_dragen_mode_variant_calling} \
      ~{false="--disable-spanning-event-genotyping" true="" use_spanning_event_genotyping} \
      ~{if defined(dragstr_model) then "--dragstr-params-path " + dragstr_model else ""} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam
  >>>

  runtime {
    maxRetries: 3
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "~{memory_size_mb} MiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{output_file_name}"
    File output_vcf_index = "~{output_file_name}.tbi"
    File bamout =  "~{vcf_basename}.bamout.bam"
    File monitoring = "monitoring.log"
    File hc_time_out = "runtime-~{scatter_i}.txt"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3000 MiB"
    disks: "local-disk ~{disk_size} HDD"
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}


# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  input {
    File interval_list
    Int scatter_count
    Int break_bands_at_multiples_of
  }

  command <<<
    set -e
    mkdir out
    java -Xms1000m -Xmx1500m -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    memory: "2000 MiB"
  }
}

# This hard filtering matches DRAGEN 3.4.12. For later DRAGEN versions, this needs to be updated.
task DragenHardFilterVcf {
  input {
    File input_vcf
    File input_vcf_index
    Boolean make_gvcf
    String vcf_basename
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.2.0"
  }

  Int disk_size = ceil(2 * size(input_vcf, "GiB")) + 20

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_vcf_name = vcf_basename + ".hard-filtered" + output_suffix

  command {
     gatk --java-options "-Xms2000m -Xmx2500m" \
      VariantFiltration \
      -V ~{input_vcf} \
      --filter-expression "QUAL < 10.4139" \
      --filter-name "DRAGENHardQUAL" \
      -O ~{output_vcf_name}
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3000 MiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
}



  
  
