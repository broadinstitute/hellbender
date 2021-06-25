version 1.0

import "GvsWarpTasks.wdl" as Tasks

workflow GvsCreateFilterSet {
   input {
        File wgs_intervals
        Int scatter_count

        File reference
        File reference_index
        File reference_dict

        String data_project
        String default_dataset
        String tmp_output_directory

        String query_project = data_project
        Array[String]? query_labels

        String fq_sample_table = "~{data_project}.~{default_dataset}.sample_info"
        String fq_alt_allele_table = "~{data_project}.~{default_dataset}.alt_allele"
        String fq_info_destination_table = "~{data_project}.~{default_dataset}.filter_set_info"
        String fq_tranches_destination_table = "~{data_project}.~{default_dataset}.filter_set_tranches"
        String fq_filter_sites_destination_table = "~{data_project}.~{default_dataset}.filter_set_sites"

        String filter_set_name

        File? excluded_intervals

        String output_file_base_name
        File? gatk_override

        File dbsnp_vcf
        File dbsnp_vcf_index

        File? snps_model
        File? indels_model

        Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]
        Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]
        Array[String] indel_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]
        Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]

        File? excluded_sites_bed

        File hapmap_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf
        File omni_resource_vcf_index
        File one_thousand_genomes_resource_vcf
        File one_thousand_genomes_resource_vcf_index
        File mills_resource_vcf
        File mills_resource_vcf_index
        File axiomPoly_resource_vcf
        File axiomPoly_resource_vcf_index
        File dbsnp_resource_vcf = dbsnp_vcf
        File dbsnp_resource_vcf_index = dbsnp_vcf_index

        # Effectively disable by default
        Int? excess_alleles_threshold = 1000000

        # Runtime attributes
        Int? small_disk_override
        Int small_disk = select_first([small_disk_override, "100"])
        Int? medium_disk_override
        Int medium_disk = select_first([medium_disk_override, "200"])
        Int? large_disk_override
        Int large_disk = select_first([large_disk_override, "300"])
        Int? huge_disk_override
        Int huge_disk = select_first([huge_disk_override, "400"])

        String? preemptible_tries_override
        Int preemptible_tries = select_first([preemptible_tries_override, "3"])

        File? service_account_json

        Int? SNP_VQSR_machine_mem_gb
        Int SNP_VQSR_downsampleFactor = 1

        Int? INDEL_VQSR_machine_mem_gb

        Int snps_variant_recalibration_threshold = 20000
    }

    call GetNumSamples {
        input:
            fq_sample_table = fq_sample_table,
            service_account_json = service_account_json,
            project_id = query_project,
            dataset_name = default_dataset
    }

    call SplitIntervals {
          input:
              intervals = wgs_intervals,
              ref_fasta = reference,
              ref_fai = reference_index,
              ref_dict = reference_dict,
              scatter_count = scatter_count
        }

    scatter(i in range(scatter_count)) {
        call ExtractFilterTask {
            input:
                gatk_override            = gatk_override,
                reference                = reference,
                reference_index          = reference_index,
                reference_dict           = reference_dict,
                fq_sample_table          = fq_sample_table,
                intervals                = SplitIntervals.interval_files[i],
                excluded_intervals       = excluded_intervals,
                fq_alt_allele_table      = fq_alt_allele_table,
                excess_alleles_threshold = excess_alleles_threshold,
                read_project_id          = query_project,
                output_file              = "${output_file_base_name}_${i}.vcf.gz",
                service_account_json     = service_account_json,
                query_project            = query_project,
                query_labels             = query_labels
        }
    }

    call MergeVCFs {
       input:
           input_vcfs = ExtractFilterTask.output_vcf,
           input_vcfs_indexes = ExtractFilterTask.output_vcf_index,
           output_vcf_name = "${output_file_base_name}.vcf.gz",
           preemptible_tries = 3,
           gatk_override = gatk_override
    }

    call IndelsVariantRecalibrator {
        input:
        sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
        sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
        model_report = indels_model,
        recalibration_filename = filter_set_name + ".indels.recal",
        tranches_filename = filter_set_name + ".indels.tranches",
        recalibration_tranche_values = indel_recalibration_tranche_values,
        recalibration_annotation_values = indel_recalibration_annotation_values,
        excluded_sites_bed = excluded_sites_bed,
        mills_resource_vcf = mills_resource_vcf,
        mills_resource_vcf_index = mills_resource_vcf_index,
        axiomPoly_resource_vcf = axiomPoly_resource_vcf,
        axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = true,
        disk_size = large_disk,
        machine_mem_gb = INDEL_VQSR_machine_mem_gb,
    }

    if (GetNumSamples.num_samples > snps_variant_recalibration_threshold) {
        call Tasks.SNPsVariantRecalibratorCreateModel {
            input:
                sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
                sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
                recalibration_filename = filter_set_name + ".snps.recal",
                tranches_filename = filter_set_name + ".snps.tranches",
                recalibration_tranche_values = snp_recalibration_tranche_values,
                recalibration_annotation_values = snp_recalibration_annotation_values,
                downsampleFactor = SNP_VQSR_downsampleFactor,
                model_report_filename = filter_set_name + ".snps.model.report",
                hapmap_resource_vcf = hapmap_resource_vcf,
                hapmap_resource_vcf_index = hapmap_resource_vcf_index,
                omni_resource_vcf = omni_resource_vcf,
                omni_resource_vcf_index = omni_resource_vcf_index,
                one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
                one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
                dbsnp_resource_vcf = dbsnp_resource_vcf,
                dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
                use_allele_specific_annotations = true,
                disk_size = small_disk
        }

        scatter (idx in range(length(ExtractFilterTask.output_vcf))) {
            call SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
                input:
                    sites_only_variant_filtered_vcf = ExtractFilterTask.output_vcf[idx],
                    sites_only_variant_filtered_vcf_index = ExtractFilterTask.output_vcf_index[idx],
                    recalibration_filename = filter_set_name + ".snps." + idx + ".recal",
                    tranches_filename = filter_set_name + ".snps." + idx + ".tranches",
                    recalibration_tranche_values = snp_recalibration_tranche_values,
                    recalibration_annotation_values = snp_recalibration_annotation_values,
                    model_report = SNPsVariantRecalibratorCreateModel.model_report,
                    hapmap_resource_vcf = hapmap_resource_vcf,
                    hapmap_resource_vcf_index = hapmap_resource_vcf_index,
                    omni_resource_vcf = omni_resource_vcf,
                    omni_resource_vcf_index = omni_resource_vcf_index,
                    one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
                    one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
                    dbsnp_resource_vcf = dbsnp_resource_vcf,
                    dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
                    use_allele_specific_annotations = true,
                    disk_size = small_disk
            }
        }

        call Tasks.GatherTranches as SNPGatherTranches {
            input:
                tranches = SNPsVariantRecalibratorScattered.tranches,
                output_filename = filter_set_name + ".snps.gathered.tranches",
                mode = "SNP",
                disk_size = small_disk,
                gatk_override = gatk_override
        }
    }

    if (GetNumSamples.num_samples <= snps_variant_recalibration_threshold) {
        call SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
            input:
                sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
                sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
                model_report = snps_model,
                recalibration_filename = filter_set_name + ".snps.recal",
                tranches_filename = filter_set_name + ".snps.tranches",
                recalibration_tranche_values = snp_recalibration_tranche_values,
                recalibration_annotation_values = snp_recalibration_annotation_values,
                excluded_sites_bed = excluded_sites_bed,
                hapmap_resource_vcf = hapmap_resource_vcf,
                hapmap_resource_vcf_index = hapmap_resource_vcf_index,
                omni_resource_vcf = omni_resource_vcf,
                omni_resource_vcf_index = omni_resource_vcf_index,
                one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
                one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
                dbsnp_resource_vcf = dbsnp_resource_vcf,
                dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
                use_allele_specific_annotations = true,
                disk_size = large_disk,
                machine_mem_gb = SNP_VQSR_machine_mem_gb,
                downsampleFactor= SNP_VQSR_downsampleFactor
        }
    }

    scatter (idx in range(length(ExtractFilterTask.output_vcf))) {
        call CreateFilterSetFiles {
            input:
                gatk_override = gatk_override,
                filter_set_name = filter_set_name,
                output_directory = tmp_output_directory,
                sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
                sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
                snp_recal_file = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
                snp_recal_file_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
                snp_recal_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
                indel_recal_file = IndelsVariantRecalibrator.recalibration,
                indel_recal_file_index = IndelsVariantRecalibrator.recalibration_index,
                indel_recal_tranches = IndelsVariantRecalibrator.tranches,
                index = idx,
                service_account_json = service_account_json,
                query_project = query_project
        }
    }

    call UploadFilterSetFilesToBQ {
        input:
            filter_set_name = filter_set_name,
            output_directory = tmp_output_directory,
            fq_info_destination_table = fq_info_destination_table,
            fq_tranches_destination_table = fq_tranches_destination_table,
            fq_filter_sites_destination_table = fq_filter_sites_destination_table,
            service_account_json = service_account_json,
            query_project = query_project
    }

    output {
        File output_vcf = MergeVCFs.output_vcf
        File output_vcf_idx = MergeVCFs.output_vcf_index
    }
}

################################################################################
task GetNumSamples {
    input {
        String fq_sample_table
        File? service_account_json
        String project_id
        String dataset_name
    }

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gcloud auth activate-service-account --key-file='~{service_account_json}'
            gcloud config set project ~{project_id}
            echo "project_id = ~{project_id}" > ~/.bigqueryrc
        fi

        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false \
        "SELECT COUNT(*) as num_rows FROM ~{fq_sample_table}" > num_rows.csv

        NUMROWS=$(python3 -c "csvObj=open('num_rows.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

        [[ $NUMROWS =~ ^[0-9]+$ ]] && echo $NUMROWS || exit 1

    >>>

    output {
        Int num_samples = read_int(stdout())
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        cpu: 1
    }
}

task ExtractFilterTask {
    # indicates that this task should NOT be call cached

    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}

    input {
        # ------------------------------------------------
        # Input args:
        File reference
        File reference_index
        File reference_dict

        String fq_sample_table

        File intervals
        File? excluded_intervals

        String fq_alt_allele_table
        String read_project_id
        String output_file
        Int? excess_alleles_threshold

        # Runtime Options:
        File? gatk_override
        File? service_account_json
        String query_project
        Array[String]? query_labels

        Int? local_sort_max_records_in_ram = 1000000
    }


    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'
    # Note the coercion of optional query_labels using select_first([expr, default])
    Array[String] query_label_args = if defined(query_labels) then prefix("--query-labels ", select_first([query_labels])) else []

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        if [ ~{has_service_account_file} = 'true' ]; then
            export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
            gcloud auth activate-service-account --key-file='~{service_account_json}'
            gcloud config set project ~{query_project}
        fi

        df -h

        gatk --java-options "-Xmx4g" \
            ExtractFeatures \
                --ref-version 38  \
                -R "~{reference}" \
                -O "~{output_file}" \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_sample_table} \
                --alt-allele-table ~{fq_alt_allele_table} \
                ~{"--excess-alleles-threshold " + excess_alleles_threshold} \
                ~{sep=" " query_label_args} \
                -L ~{intervals} \
                ~{"-XL " + excluded_intervals} \
                --project-id ~{read_project_id}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "7 GB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 15
        preemptible: 3
        cpu: 2
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf = "~{output_file}"
        File output_vcf_index = "~{output_file}.tbi"
    }
 }

 task SplitIntervals {
     input {
         File intervals
         File ref_fasta
         File ref_fai
         File ref_dict
         Int scatter_count
         String? split_intervals_extra_args

         File? gatk_override
     }

     parameter_meta {
         intervals: {
             localization_optional: true
         }
         ref_fasta: {
             localization_optional: true
         }
         ref_fai: {
             localization_optional: true
         }
         ref_dict: {
             localization_optional: true
         }
      }

      command {
          set -e
          export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

          mkdir interval-files
          gatk --java-options "-Xmx5g" SplitIntervals \
              -R ~{ref_fasta} \
              ~{"-L " + intervals} \
              -scatter ~{scatter_count} \
              -O interval-files \
              ~{split_intervals_extra_args}
          cp interval-files/*.interval_list .
      }

      runtime {
          docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
          bootDiskSizeGb: 15
          memory: "6 GB"
          disks: "local-disk 10 HDD"
          preemptible: 3
          cpu: 1
      }

      output {
          Array[File] interval_files = glob("*.interval_list")
      }
  }

task CreateFilterSetFiles {
    # indicates that this task should NOT be call cached

    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}

    input {
        # ------------------------------------------------
        # Input args:

        # Runtime Options:
        File? gatk_override

        String filter_set_name
        String output_directory

        File sites_only_variant_filtered_vcf
        File sites_only_variant_filtered_vcf_index

        File snp_recal_file
        File snp_recal_file_index
        File snp_recal_tranches

        File indel_recal_file
        File indel_recal_file_index
        File indel_recal_tranches

        String index

        File? service_account_json
        String query_project
    }

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -eo pipefail

        if [ ~{has_service_account_file} = 'true' ]; then
            export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
            gcloud auth activate-service-account --key-file='~{service_account_json}'
            gcloud config set project ~{query_project}
        fi

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx1g" \
            CreateFilteringFiles \
            --ref-version 38 \
            --filter-set-name ~{filter_set_name} \
            -mode SNP \
            -V ~{snp_recal_file} \
            -O ~{filter_set_name}.snps.recal.tsv

        gatk --java-options "-Xmx1g" \
            CreateFilteringFiles \
            --ref-version 38 \
            --filter-set-name ~{filter_set_name} \
            -mode INDEL \
            -V ~{indel_recal_file} \
            -O ~{filter_set_name}.indels.recal.tsv

        gatk --java-options "-Xmx1g" \
            CreateSiteFilteringFiles \
            --ref-version 38 \
            --filter-set-name ~{filter_set_name} \
            -V ~{sites_only_variant_filtered_vcf} \
            -O filter_sites_load.~{index}.tsv

        # merge into a single file
        echo "Merging SNP + INDELs"
        cat ~{filter_set_name}.snps.recal.tsv ~{filter_set_name}.indels.recal.tsv | grep -v filter_set_name | grep -v "#"  > filter_set_load.~{index}.tsv

        echo "Merging Tranches"
        cat ~{snp_recal_tranches} ~{indel_recal_tranches} | grep -v targetTruthSensitivity | grep -v "#" | awk -v CALLSET=~{filter_set_name} '{ print CALLSET "," $0 }' > tranches_load.~{index}.csv

        if [ -n "~{output_directory}" ]; then
            gsutil cp filter_sites_load.~{index}.tsv ~{output_directory}/
            gsutil cp filter_set_load.~{index}.tsv ~{output_directory}/
            gsutil cp tranches_load.~{index}.csv ~{output_directory}/
        fi
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "3500 MB"
        disks: "local-disk 200 HDD"
        bootDiskSizeGb: 15
        preemptible: 0
        cpu: 1
    }

    # ------------------------------------------------
    # Outputs:
    output {
        String filter_sites_load = "~{output_directory}/filter_sites_load.~{index}.tsv"
        String filter_set_load = "~{output_directory}/filter_set_load.~{index}.tsv"
        String tranches_load = "~{output_directory}/tranches_load.~{index}.tsv"
    }
}

task UploadFilterSetFilesToBQ {
    # indicates that this task should NOT be call cached

    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}

    input {
        # ------------------------------------------------
        # Input args:

        String filter_set_name
        String output_directory

        String fq_info_destination_table
        String fq_tranches_destination_table
        String fq_filter_sites_destination_table

        File? service_account_json
        String query_project
    }

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -eo pipefail

        if [ ~{has_service_account_file} = 'true' ]; then
            export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
            gcloud auth activate-service-account --key-file='~{service_account_json}'
            gcloud config set project ~{query_project}
        fi

        # BQ load likes a : instead of a . after the project
        bq_info_table=$(echo ~{fq_info_destination_table} | sed s/\\./:/)

        echo "Loading Filter Set into BQ"
        bq load --project_id=~{query_project} --skip_leading_rows 0 -F "tab" \
        ${bq_info_table} \
        ~{output_directory}/filter_set_load.*.tsv \
        "filter_set_name:string,type:string,location:integer,ref:string,alt:string,vqslod:float,culprit:string,training_label:string,yng_status:string" > status_bq_load_info_table

        # BQ load likes a : instead of a . after the project
        bq_tranches_table=$(echo ~{fq_tranches_destination_table} | sed s/\\./:/)

        echo "Loading Tranches into BQ"
        bq load --project_id=~{query_project} --skip_leading_rows 0 -F "," \
        ${bq_tranches_table} \
        ~{output_directory}/tranches_load.*.tsv \
        "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float"  > status_bq_load_tranches_table

        # BQ load likes a : instead of a . after the project
        bq_filter_sites_table=$(echo ~{fq_filter_sites_destination_table} | sed s/\\./:/)

        # Creating site
        echo "Loading Filter Set Sites into BQ"
        bq load --project_id=~{query_project} --skip_leading_rows 1 -F "tab" \
        ${bq_filter_sites_table} \
        ~{output_directory}/filter_sites_load.*.tsv \
        "filter_set_name:string,location:integer,filters:string"  > status_bq_load_filter_sites_table
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "3500 MB"
        disks: "local-disk 200 HDD"
        bootDiskSizeGb: 15
        preemptible: 0
        cpu: 1
    }

    # ------------------------------------------------
    # Outputs:
    output {
        String status_bq_load_info_table = read_string("status_bq_load_info_table")
        String status_bq_load_tranches_table = read_string("status_bq_load_tranches_table")
        String status_bq_load_filter_sites_table = read_string("status_bq_load_filter_sites_table")
    }
 }

 task MergeVCFs {
    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}

   input {
     Array[File] input_vcfs
     Array[File] input_vcfs_indexes
     String output_vcf_name

     File? gatk_override
     Int preemptible_tries
   }

   Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

    parameter_meta {
         input_vcfs: {
             localization_optional: true
         }
         input_vcfs_indexes: {
             localization_optional: true
         }
      }

   command {
       export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

       gatk --java-options -Xmx3g GatherVcfsCloud \
            --ignore-safety-checks --gather-type BLOCK \
            -I ~{sep=' -I ' input_vcfs} \
            --output ~{output_vcf_name}

       tabix ~{output_vcf_name}

   }
   runtime {
     docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
     preemptible: preemptible_tries
     memory: "3 GiB"
     disks: "local-disk ~{disk_size} HDD"
   }
   output {
     File output_vcf = "~{output_vcf_name}"
     File output_vcf_index = "~{output_vcf_name}.tbi"
   }
 }

task IndelsVariantRecalibrator {

  input {
    String recalibration_filename
    String tranches_filename
    File? model_report

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File? excluded_sites_bed
    File mills_resource_vcf
    File axiomPoly_resource_vcf
    File dbsnp_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 4

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    Int? machine_mem_gb
  }

  Int machine_mem = select_first([machine_mem_gb, 30])
  Int java_mem = machine_mem - 5


  command <<<
    set -euo pipefail

    gatk --java-options -Xmx~{java_mem}g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      ~{"-XL " + excluded_sites_bed} \
      -O ~{recalibration_filename} \
      --output-model indels.model \
      --rscript-file indels.Rscript \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode INDEL \
      ~{"--input-model " + model_report} \
      --max-gaussians ~{max_gaussians} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} GiB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
    File model = "indels.model"
    File rscript = "indels.Rscript"
  }
}

task SNPsVariantRecalibrator {

  input {
    String recalibration_filename
    String tranches_filename
    File? model_report

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File? excluded_sites_bed
    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 6
    Int? downsampleFactor = 1

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    Int? machine_mem_gb

  }

  Int auto_mem = ceil(2 * size([sites_only_variant_filtered_vcf,
                                hapmap_resource_vcf,
                                omni_resource_vcf,
                                one_thousand_genomes_resource_vcf,
                                dbsnp_resource_vcf],
                     "GiB"))
  Int machine_mem = select_first([machine_mem_gb, if auto_mem < 7 then 7 else auto_mem])
  Int java_mem = machine_mem - 5
  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  command <<<
    set -euo pipefail

    MODEL_REPORT=~{model_report}

    gatk --java-options -Xmx~{java_mem}g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      ~{"-XL " + excluded_sites_bed} \
      -O ~{recalibration_filename} \
      --output-model snps.model \
      --rscript-file snps.Rscript \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode SNP \
      --sample-every-Nth-variant ~{downsampleFactor} \
      ~{model_report_arg} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} GiB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 0
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
    File model = "snps.model"
    File rscript = "snps.Rscript"
  }
}





