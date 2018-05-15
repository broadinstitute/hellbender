# This workflow takes an input CRAM to call variants with HaplotypeCaller
# Then filters the calls with the CNNVariant neural net tool
# The site-level scores are added to the INFO field of the VCF.
# The architecture arguments, info_key and tensor type arguments MUST be in agreement
# (e.g. 2D models must have tensor_type of read_tensor and info_key CNN_2D, 1D models have tensor_type reference and info_key CNN_1D)
# The INFO field key will be "1D_CNN" or "2D_CNN" depending on the neural net architecture used for inference.
# The architecture arguments specify pre-trained networks.
# New networks can be trained by the GATK tools: CNNVariantWriteTensors and CNNVariantTrain
# The CRAM could be generated by the single-sample pipeline
# (https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl)
# Also accepts a BAM as the input file in which case a BAM index is required as well.

import "cnn_variant_common_tasks.wdl" as CNNTasks

workflow Cram2FilteredVcf {
    File input_file                  # Aligned CRAM file or Aligned BAM files
    File? input_file_index           # Index for an aligned BAM file if that is the input, unneeded if input is a CRAM
    File reference_fasta 
    File reference_dict
    File reference_fasta_index
    File resource_fofn               # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
    File resource_fofn_index         # File of VCF file indices of resources
    File? architecture_json          # Neural Net configuration for CNNScoreVariants
    File? architecture_hd5           # Pre-Trained weights and architecture for CNNScoreVariants
    Int? inference_batch_size        # Batch size for python in CNNScoreVariants
    Int? transfer_batch_size         # Batch size for java in CNNScoreVariants
    Int? intra_op_threads            # Tensorflow threading within nodes
    Int? inter_op_threads            # Tensorflow threading between nodes
    String output_prefix             # Identifying string for this run will be used to name all output files
    String? tensor_type              # What kind of tensors the Neural Net expects (e.g. reference, read_tensor)
    String info_key                  # The score key for the info field of the vcf (e.g. CNN_1D, CNN_2D)
    String snp_tranches              # Filtering threshold(s) for SNPs in terms of sensitivity to overlapping known variants in resources
    String indel_tranches            # Filtering threshold(s) for INDELs in terms of sensitivity to overlapping known variants in resources
    File? gatk_override
    String gatk_docker
    File calling_intervals
    Int scatter_count                # Number of shards for parallelization of HaplotypeCaller and CNNScoreVariants
    String extra_args                # Extra arguments for HaplotypeCaller

    # Runtime parameters
    Int? mem_gb
    Int? preemptible_attempts
    Float? disk_space_gb
    Int? cpu

    Int? increase_disk_size
    Int additional_disk = select_first([increase_disk_size, 20])
    Float ref_size = size(reference_fasta, "GB") + size(reference_fasta_index, "GB") + size(reference_dict, "GB")

    # Clunky check to see if the input is a BAM or a CRAM
    if (basename(input_file) == basename(input_file, ".bam")){
        call CNNTasks.CramToBam {
            input:
              reference_fasta = reference_fasta,
              reference_dict = reference_dict,
              reference_fasta_index = reference_fasta_index,
              cram_file = input_file,
              output_prefix = output_prefix,
              disk_space_gb = round(4*size(input_file, "GB") + ref_size + additional_disk),
              preemptible_attempts = preemptible_attempts
        }
    }

    call CNNTasks.SplitIntervals {
        input:
            gatk_override = gatk_override,
            scatter_count = scatter_count,
            intervals = calling_intervals,
            ref_fasta = reference_fasta,
            ref_dict = reference_dict,
            ref_fai = reference_fasta_index,
            gatk_docker = gatk_docker,
            disk_space = round(additional_disk + ref_size)
    }

    String input_bam = select_first([CramToBam.output_bam, input_file])
    Float bam_size = size(input_bam, "GB")

    scatter (calling_interval in SplitIntervals.interval_files) {
        call CNNTasks.RunHC4 {
            input:
                input_bam = input_bam,
                input_bam_index = select_first([CramToBam.output_bam_index, input_file_index]),
                reference_fasta = reference_fasta,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                extra_args = extra_args,
                disk_space_gb = round(bam_size + ref_size + additional_disk)
        }

        call CNNTasks.CNNScoreVariants {
            input:
                input_vcf = RunHC4.raw_vcf,
                input_vcf_index = RunHC4.raw_vcf_index,
                bam_file = RunHC4.bamout,
                bam_file_index = RunHC4.bamout_index,
                architecture_json = architecture_json,
                architecture_hd5 = architecture_hd5,
                reference_fasta = reference_fasta,
                tensor_type = tensor_type,
                inference_batch_size = inference_batch_size,
                transfer_batch_size = transfer_batch_size,
                intra_op_threads = intra_op_threads,
                inter_op_threads = inter_op_threads,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,               
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts,
                mem_gb = mem_gb,
                disk_space_gb = round((bam_size/scatter_count) + ref_size + additional_disk)
        }
    }

    call CNNTasks.MergeVCFs as MergeVCF_HC4 {
        input: 
            input_vcfs = CNNScoreVariants.cnn_annotated_vcf,
            output_prefix = output_prefix,
            gatk_override = gatk_override,
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker,
            disk_space_gb = additional_disk
    }

    call CNNTasks.FilterVariantTranches {
        input:
            input_vcf = MergeVCF_HC4.merged_vcf,
            input_vcf_index = MergeVCF_HC4.merged_vcf_index,
            resource_fofn = resource_fofn,
            resource_fofn_index = resource_fofn_index,
            output_prefix = output_prefix,
            snp_tranches = snp_tranches,
            indel_tranches = indel_tranches,
            info_key = info_key,
            gatk_override = gatk_override,
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker,
            disk_space_gb = additional_disk
    }

    call CNNTasks.SamtoolsMergeBAMs {
        input:
            input_bams = RunHC4.bamout,
            output_prefix = output_prefix,
            disk_space_gb = round(bam_size + ref_size + additional_disk)
    }

    output {
        FilterVariantTranches.*
    }
}
