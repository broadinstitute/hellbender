version 1.0

import "GvsUtils.wdl" as Utils
import "../variant_annotations_table/GvsCreateVATFilesFromBigQuery.wdl" as GvsCreateVATFilesFromBigQuery

workflow GvsCreateVATfromVDS {
    input {
        String project_id
        String dataset_name
        File ancestry_file
        String filter_set_name
        File? sites_only_vcf
        File? sites_only_vcf_index
        File? vds_path
        String output_path

        String? basic_docker
        String? git_branch_or_tag
        String? hail_version
        String? vat_version
        String? workspace_gcs_project

        Boolean leave_hail_cluster_running_at_end = false
        Int? merge_vcfs_disk_size_override
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Int? split_intervals_scatter_count
        Boolean use_classic_VQSR = false
        Boolean use_reference_disk = true

        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? gatk_docker
        String? variants_docker
        String? variants_nirvana_docker
    }

    parameter_meta {
        project_id: {
            help: "Google project ID for the GVS BigQuery dataset"
        }
        dataset_name: {
            help: "BigQuery dataset name for GVS"
        }
        ancestry_file: {
            help: "TSV file in GCS where the first column is the research ID and the last is the derived ancestry"
        }
        filter_set_name: {
            help: "name of the filter set used to generate the callset in GVS"
        }
        sites_only_vcf: {
            help: "Optional sites-only VCF file. If defined, generation of a sites-only VCF from the VDS will be skipped. If defined, then 'vds_path' must NOT be defined."
        }
        output_path: {
            help: "GCS location (with a trailing '/') to put temporary and output files for the VAT pipeline"
        }
        vds_path: {
            help: "Optional top-level directory of the GVS VDS to be used to create the VAT. If defined, then 'sites_only_vcf' must NOT be defined"
        }
    }


    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    String region = "us-central1"
    String gcs_subnetwork_name = "subnetwork"

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_variants_nirvana_docker = select_first([variants_nirvana_docker, GetToolVersions.variants_nirvana_docker])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])
    String effective_google_project = select_first([workspace_gcs_project, GetToolVersions.google_project])

    # If the vat version is undefined or v1 then the vat tables would be named like filter_vat, otherwise filter_vat_v2.
    String effective_vat_version = if (defined(vat_version) && select_first([vat_version]) != "v1") then "_" + select_first([vat_version]) else ""
    String vat_table_name = filter_set_name + "_vat" + effective_vat_version

    String output_path_without_a_trailing_slash = sub(output_path, "/$", "")
    String effective_output_path = if (output_path == output_path_without_a_trailing_slash) then output_path + "/" else output_path

    if ((defined(sites_only_vcf)) && (defined(vds_path))) {
        call Utils.TerminateWorkflow as IfSitesOnlyVcfSetDontSetCreateParameters {
            input:
                message = "Error: If 'sites_only_vcf' is set as an input, you may not set 'vds_path'",
                basic_docker = effective_basic_docker,
        }
    }

    if (!defined(sites_only_vcf) && (!defined(vds_path))) {
        call Utils.TerminateWorkflow as MustSetSitesOnlyVcfCreationParameters {
            input:
                message = "Error: If 'sites_only_vcf' is not set as an input, you MUST set 'vds_path'",
                basic_docker = effective_basic_docker,
        }
    }

    if (defined(sites_only_vcf) || (defined(vds_path))) {
        if (!defined(split_intervals_scatter_count)) {
            call Utils.GetBQTableLastModifiedDatetime as SampleDateTime {
                input:
                    project_id = project_id,
                    fq_table = "~{project_id}.~{dataset_name}.sample_info",
                    cloud_sdk_docker = effective_cloud_sdk_docker,
            }

            call Utils.GetNumSamplesLoaded {
                input:
                    fq_sample_table ="~{project_id}.~{dataset_name}.sample_info",
                    project_id = project_id,
                    sample_table_timestamp = SampleDateTime.last_modified_timestamp,
                    cloud_sdk_docker = effective_cloud_sdk_docker,
            }

            Int calculated_scatter_count = if (GetNumSamplesLoaded.num_samples < 11) then 10 else
                                               if (GetNumSamplesLoaded.num_samples < 250000) then 500 else
                                                   if (GetNumSamplesLoaded.num_samples < 450000) then 1000 else 2000
        }

        Int effective_scatter_count = select_first([split_intervals_scatter_count, calculated_scatter_count])

        call MakeSubpopulationFilesAndReadSchemaFiles {
            input:
                input_ancestry_file = ancestry_file,
                variants_docker = effective_variants_docker,
        }

        if (!defined(sites_only_vcf)) {
            call GenerateSitesOnlyVcf {
                input:
                    vds_path = select_first([vds_path]),
                    use_classic_VQSR = use_classic_VQSR,
                    workspace_project = effective_google_project,
                    hail_version = effective_hail_version,
                    ancestry_file_path = MakeSubpopulationFilesAndReadSchemaFiles.ancestry_file_path,
                    workspace_bucket = GetToolVersions.workspace_bucket,
                    region = region,
                    gcs_subnetwork_name = gcs_subnetwork_name,
                    leave_cluster_running_at_end = leave_hail_cluster_running_at_end,
                    variants_docker = effective_variants_docker,
            }
        }

        call Utils.CopyFile as CopySitesOnlyVcf {
            input:
                input_file = select_first([sites_only_vcf, GenerateSitesOnlyVcf.sites_only_vcf]),
                output_gcs_dir = effective_output_path + "sites_only_vcf",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        if (!defined(sites_only_vcf) || !defined(sites_only_vcf_index)) {
            call Utils.IndexVcf {
                input:
                    input_vcf = CopySitesOnlyVcf.output_file_path,
                    gatk_docker = effective_gatk_docker,
            }
        }

        call Utils.CopyFile as CopySitesOnlyVcfIndex {
            input:
                input_file = select_first([sites_only_vcf_index, IndexVcf.output_vcf_index]),
                output_gcs_dir = effective_output_path + "sites_only_vcf",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call Utils.SplitIntervals {
            input:
                intervals = interval_list,
                ref_fasta = reference,
                ref_fai = reference_index,
                ref_dict = reference_dict,
                scatter_count = effective_scatter_count,
                output_gcs_dir = effective_output_path + "intervals",
                split_intervals_disk_size_override = split_intervals_disk_size_override,
                split_intervals_mem_override = split_intervals_mem_override,
                gatk_docker = effective_gatk_docker,
        }

        String sites_only_vcf_basename = basename(CopySitesOnlyVcf.output_file_path, ".sites-only.vcf")

        scatter(i in range(length(SplitIntervals.interval_files))) {
            String interval_file_basename = basename(SplitIntervals.interval_files[i], ".interval_list")
            String vcf_filename = interval_file_basename + "." + sites_only_vcf_basename

            call Utils.SelectVariants {
                input:
                    input_vcf = CopySitesOnlyVcf.output_file_path,
                    input_vcf_index = CopySitesOnlyVcfIndex.output_file_path,
                    interval_list = SplitIntervals.interval_files[i],
                    output_basename = vcf_filename,
                    gatk_docker = effective_gatk_docker,
            }

            call RemoveDuplicatesFromSitesOnlyVCF {
                input:
                    sites_only_vcf = SelectVariants.output_vcf,
                    ref = reference,
                    variants_docker = effective_variants_docker,
            }

            call StripCustomAnnotationsFromSitesOnlyVCF {
                input:
                    input_vcf = RemoveDuplicatesFromSitesOnlyVCF.output_vcf,
                    custom_annotations_header = MakeSubpopulationFilesAndReadSchemaFiles.custom_annotations_template_file,
                    output_vcf_name = "${vcf_filename}.unannotated.sites_only.vcf",
                    output_custom_annotations_filename = "${vcf_filename}.custom_annotations.tsv",
                    variants_docker = effective_variants_docker,
            }

            ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
            call AnnotateVCF {
                input:
                    input_vcf = StripCustomAnnotationsFromSitesOnlyVCF.output_vcf,
                    output_annotated_file_name = "${vcf_filename}_annotated",
                    custom_annotations_file = StripCustomAnnotationsFromSitesOnlyVCF.output_custom_annotations_file,
                    variants_nirvana_docker = effective_variants_nirvana_docker,
                    use_reference_disk = use_reference_disk,
            }

            call PrepVtAnnotationJson {
                input:
                    positions_annotation_json = AnnotateVCF.positions_annotation_json,
                    output_file_suffix = "${vcf_filename}.json.gz",
                    output_path = effective_output_path,
                    variants_docker = effective_variants_docker,
            }

            call PrepGenesAnnotationJson {
                input:
                    genes_annotation_json = AnnotateVCF.genes_annotation_json,
                    output_file_suffix = "${vcf_filename}.json.gz",
                    output_path = effective_output_path,
                    variants_docker = effective_variants_docker,
            }

        }

        call Utils.MergeTsvs {
            input:
                input_files = RemoveDuplicatesFromSitesOnlyVCF.track_dropped,
                output_file_name = "${sites_only_vcf_basename}.dropped.tsv",
                basic_docker = effective_basic_docker,
        }

        call BigQueryLoadJson {
            input:
                nirvana_schema = MakeSubpopulationFilesAndReadSchemaFiles.vat_schema_json_file,
                vt_schema = MakeSubpopulationFilesAndReadSchemaFiles.variant_transcript_schema_json_file,
                genes_schema = MakeSubpopulationFilesAndReadSchemaFiles.genes_schema_json_file,
                project_id = project_id,
                dataset_name = dataset_name,
                output_path = effective_output_path,
                base_vat_table_name = vat_table_name,
                prep_vt_json_done = PrepVtAnnotationJson.done,
                prep_genes_json_done = PrepGenesAnnotationJson.done,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call DeduplicateVatInBigQuery {
            input:
                input_vat_table_name = BigQueryLoadJson.vat_table,
                output_vat_table_name = vat_table_name,
                nirvana_schema = MakeSubpopulationFilesAndReadSchemaFiles.vat_schema_json_file,
                project_id = project_id,
                dataset_name = dataset_name,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call GvsCreateVATFilesFromBigQuery.GvsCreateVATFilesFromBigQuery {
            input:
                project_id = project_id,
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                dataset_name = dataset_name,
                vat_table_name = DeduplicateVatInBigQuery.vat_table,
                output_path = effective_output_path,
                merge_vcfs_disk_size_override = merge_vcfs_disk_size_override,
                precondition_met = BigQueryLoadJson.done,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
       }
    }

    output {
        String? cluster_name = GenerateSitesOnlyVcf.cluster_name
        File? dropped_sites_file = MergeTsvs.output_file
        File? final_tsv_file = GvsCreateVATFilesFromBigQuery.final_tsv_file
        String recorded_git_hash = GetToolVersions.git_hash
    }
}


task GenerateSitesOnlyVcf {
    input {
        String vds_path
        Boolean use_classic_VQSR
        String workspace_project
        String workspace_bucket
        String region
        String gcs_subnetwork_name
        Boolean leave_cluster_running_at_end
        String hail_version
        String ancestry_file_path
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction

        String variants_docker
    }
    String prefix = "sites-only-vcf"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        pip3 install --upgrade pip
        pip3 install hail~{'==' + hail_version}

        pip3 install --upgrade google-cloud-dataproc ijson

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        # hex="$(head -c4 < /dev/urandom | xxd -p)-$(head -c2 < /dev/urandom | xxd -p)"
        # HACK hard code the hex so sites-only Hail table pieces have matching names between runs.
        hex="997d8faf-6cf5"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        sites_only_vcf_filename="~{workspace_bucket}/~{prefix}-${hex}.sites-only.vcf.bgz"
        echo ${sites_only_vcf_filename} > sites_only_vcf_filename.txt

        hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"

        # construct a JSON of arguments for python script to be run in the hail cluster
        cat > script-arguments.json <<FIN
        {
            "vds_input_path": "~{vds_path}",
            "temp_path": "${hail_temp_path}",
            "ancestry_input_path": "~{ancestry_file_path}",
            "sites_only_output_path" : "${sites_only_vcf_filename}"
        }
        FIN

        # Run the hail python script to make a sites-only VCF from a VDS
        # - The autoscaling policy gvs-autoscaling-policy will exist already from the VDS creation
        python3 /app/run_in_hail_cluster.py \
            --script-path /app/hail_create_vat_inputs.py \
            --secondary-script-path-list /app/create_vat_inputs.py \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            --autoscaling-policy gvs-autoscaling-policy \
            --region ~{region} \
            --gcs-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}
    >>>

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 0
        docker: variants_docker
        bootDiskSizeGb: 10
    }

    output {
        String cluster_name = read_string("cluster_name.txt")
        File sites_only_vcf = read_string("sites_only_vcf_filename.txt")
    }
}


task MakeSubpopulationFilesAndReadSchemaFiles {
    input {
        File input_ancestry_file

        String schema_filepath = "/data/variant_annotation_table/schema/"
        String vat_schema_json_filename = "vat_schema.json"
        String variant_transcript_schema_json_filename = "vt_schema.json"
        String genes_schema_json_filename = "genes_schema.json"
        String variants_docker
    }
    String output_ancestry_filename =  "ancestry_mapping.tsv"
    String custom_annotations_template_filename =  "custom_annotations_template.tsv"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cp ~{schema_filepath}* .

        ## the ancestry file is processed down to a simple mapping from sample to subpopulation
        python3 /app/extract_subpop.py \
            --input_path ~{input_ancestry_file} \
            --output_path ~{output_ancestry_filename} \
            --custom_annotations_template_path ~{custom_annotations_template_filename}
    >>>

    runtime {
        docker: variants_docker
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    output {
        File vat_schema_json_file = vat_schema_json_filename
        File variant_transcript_schema_json_file = variant_transcript_schema_json_filename
        File genes_schema_json_file = genes_schema_json_filename

        File ancestry_mapping_list = output_ancestry_filename
        File custom_annotations_template_file = custom_annotations_template_filename
        String ancestry_file_path = input_ancestry_file
    }
}


task StripCustomAnnotationsFromSitesOnlyVCF {
    input {
        File input_vcf
        File custom_annotations_header
        String output_vcf_name
        String output_custom_annotations_filename
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil((size(input_vcf, "GB") + size(custom_annotations_header, "GB")) * 4) + 100

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        python3 /app/strip_custom_annotations_from_sites_only_vcf.py \
        --input_vcf ~{input_vcf} \
        --input_custom_annotations_tsv ~{custom_annotations_header} \
        --output_vcf ~{output_vcf_name} \
        --output_custom_annotations_tsv ~{output_custom_annotations_filename}

    >>>

    runtime {
        docker: variants_docker
        memory: "7 GiB"
        cpu: "2"
        preemptible: 3
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf = output_vcf_name
        File output_custom_annotations_file = output_custom_annotations_filename
        File monitoring_log = "monitoring.log"
    }
}


task RemoveDuplicatesFromSitesOnlyVCF {
    input {
        File sites_only_vcf
        File ref
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil(size(sites_only_vcf, "GB") * 5) + 100

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # custom function to prepend the current datetime to an echo statement
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o sites_only.bcf ~{sites_only_vcf}

        echo_date "VAT: Calculating number of sites with Ns"

        ## track the dropped variants with N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view --threads 4 -i 'REF~"N"' -O u sites_only.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with N's in the reference AND sites with AC=0"
        ## NOTE: Sites that were filtered out because of AC=0 are not recorded in the 'track_dropped.tsv' file, but can be
        ##       determined by examining the sites-only VCF provided to this WDL.
        bcftools view --threads 4 -e 'REF~"N" || AC=0' -O b sites_only.bcf -o filtered_sites_only.bcf
        rm sites_only.bcf

        echo_date "VAT: normalize, left align and split multi allelic sites to new lines, remove duplicate lines"
        ## note that normalization may create sites with more than 50 alt alleles
        bcftools norm --threads 4 -m- --check-ref w -f ~{ref} filtered_sites_only.bcf -O b -o normalized.bcf
        rm filtered_sites_only.bcf

        echo_date "VAT: detecting and removing duplicate rows from sites-only VCF"

        ## During normalization, sometimes duplicate variants appear but with different calculations. This seems to be a bug in bcftools. For now we are dropping all duplicate variants
        ## to locate the duplicates, we first make a file of just the first 5 columns
        bcftools query normalized.bcf -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | sort | uniq -d > duplicates.tsv

        echo_date "VAT: done with duplicate detection"
        wc -l duplicates.tsv
        echo_date "VAT: Duplicates may have been found"

        # If there ARE dupes to remove
        if [ -s duplicates.tsv ]; then
            ## remove those rows (that match up to the first 5 cols)
            echo_date "VAT: Removing those rows"
            bcftools view --threads 4 normalized.bcf | grep -v -wFf duplicates.tsv > deduplicated.vcf
        else
            # There are no duplicates to remove
            echo_date "VAT: No duplicates found"
            bcftools view --threads 4 normalized.bcf -o deduplicated.vcf
        fi
        rm normalized.bcf

        ## add duplicates to the file that's tracking dropped variants
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        echo_date "VAT: finished"
    >>>

    runtime {
        docker: variants_docker
        maxRetries: 3
        memory: "16 GB"
        preemptible: 3
        cpu: "8"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File track_dropped = "track_dropped.tsv"
        File output_vcf = "deduplicated.vcf"
        File monitoring_log = "monitoring.log"
    }
}


task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File custom_annotations_file

        # Mentioning this path in the inputs section of the task combined with checking the 'Use reference disks' option
        # in Terra UI tells Cromwell to arrange for the Nirvana reference disk to be attached to this VM.
        File summon_reference_disk =
            "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/MITOMAP_20200819.nsa.idx"

        String variants_nirvana_docker

        File omim_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/OMIM_20220516.nga"
        File cosmic_gene_fusion_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/COSMIC_GeneFusions_94.gfj"
        File primate_ai_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/PrimateAI_0.2.nsa"
        File primate_ai_annotations_idx = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/PrimateAI_0.2.nsa.idx"
        File splice_ai_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/SpliceAi_1.3.nsa"
        File splice_ai_annotations_idx = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/SpliceAi_1.3.nsa.idx"
        Boolean use_reference_disk
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String gene_annotation_json_name = output_annotated_file_name + ".genes.json.gz"
    String positions_annotation_json_name = output_annotated_file_name + ".positions.json.gz"
    String nirvana_location = "/Nirvana/Nirvana.dll"
    String custom_creation_location = "/Nirvana/SAUtils.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        if [[ "~{use_reference_disk}" == "true" ]]
        then
            # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
            # disk image was built: while all the reference files do exist on the image they are not at the expected
            # locations. The following code works around this issue and should continue to work even after a corrected
            # version of the Nirvana reference image is deployed into Terra.

            # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
            # mount point" because we do not actually confirm this is a reference disk until the following code block.
            CANDIDATE_MOUNT_POINT=$(lsblk | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
            if [[ -z ${CANDIDATE_MOUNT_POINT} ]]; then
                >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
                exit 1
            fi

            # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
            # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
            REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
            REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
            if [[ -z ${REFERENCE_PATH} ]]; then
                >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
                exit 1
            fi

            # Take the parent of the parent directory of this file as root of the locally mounted references:
            DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"
        else
            DATA_SOURCES_FOLDER=/cromwell_root/nirvana_references
            mkdir ${DATA_SOURCES_FOLDER}

            # Download the references
            dotnet /Nirvana/Downloader.dll --ga GRCh38 --out ${DATA_SOURCES_FOLDER}

            # As of 2024-01-24 OMIM is no longer included among the bundle of annotation resources pulled down by the
            # Nirvana downloader. As this annotation set is currently central for our VAT logic, special-case link in
            # the OMIM .nsa bundle we downloaded back when we made the Delta reference disk:
            ln ~{omim_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            # Similarly, the following annotations were removed from the latest Nirvana annotations (3.18.1), but we
            # re-add them as desired by Lee
            ln ~{cosmic_gene_fusion_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{primate_ai_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{primate_ai_annotations_idx} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{splice_ai_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{splice_ai_annotations_idx} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
        fi

        # =======================================
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        ## use --skip-ref once you are on a version of nirvana later than 3.14 (once they have created a docker image for it)
        dotnet ~{custom_creation_location} customvar \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -i ~{custom_annotations_file} \
            -o $CUSTOM_ANNOTATIONS_FOLDER

        # =======================================
        # Create Nirvana annotations:

        dotnet ~{nirvana_location} \
            -i ~{input_vcf} \
            -c $DATA_SOURCES_FOLDER~{path} \
            --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
            --sd $CUSTOM_ANNOTATIONS_FOLDER \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -o ~{output_annotated_file_name}

        # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
        # Parse out the Genes section into a separate annotated json
        dotnet  ~{jasix_location} \
            --in ~{annotation_json_name} \
            --section genes \
            --out ~{gene_annotation_json_name}

        # Parse out the Positions section into a separate annotated json
        dotnet  ~{jasix_location} \
        --in ~{annotation_json_name} \
        --section positions \
        --out ~{positions_annotation_json_name}

    >>>

    runtime {
        docker: variants_nirvana_docker
        memory: "128 GB"
        cpu: "4"
        preemptible: 1
        maxRetries: 1
        disks: "local-disk 2000 HDD"
    }

    output {
        File genes_annotation_json = "~{gene_annotation_json_name}"
        File positions_annotation_json = "~{positions_annotation_json_name}"
        File monitoring_log = "monitoring.log"
    }
}

task PrepVtAnnotationJson {
    input {
        File positions_annotation_json
        String output_file_suffix
        String output_path
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Kick off the monitoring script
        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_vt_bqloadjson_from_annotations.py \
            --annotated_json ~{positions_annotation_json} \
            --output_vt_json ~{output_vt_json}

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'

    >>>

    runtime {
        docker: variants_docker
        memory: "16 GB"
        preemptible: 2
        cpu: "1"
        disks: "local-disk 500 HDD"
    }

    output {
        File vat_vt_json="~{output_vt_json}"
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}

task PrepGenesAnnotationJson {
    input {
        File genes_annotation_json
        String output_file_suffix
        String output_path
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_genes_gcp_path = output_path + 'genes/'

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_genes_bqloadjson_from_annotations.py \
            --annotated_json ~{genes_annotation_json} \
            --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

    >>>

    runtime {
        docker: variants_docker
        memory: "7 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 500 HDD"
    }

    output {
        File vat_genes_json="~{output_genes_json}"
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}


task BigQueryLoadJson {
    meta {
        # since the WDL will not see the updated data (its getting put in a gcp bucket)
        volatile: true
    }

    input {
        String base_vat_table_name
        File? nirvana_schema
        File? vt_schema
        File? genes_schema
        String project_id
        String dataset_name
        String output_path
        Array[Boolean] prep_vt_json_done
        Array[Boolean] prep_genes_json_done
        String cloud_sdk_docker
    }

    # This is the name of the vat table. Due to sharding (VS-1191) there may be some duplicated entries.
    # So we create it here, and then deduplicate it in a later step
    String vat_table_name = base_vat_table_name + "_w_dups"

    String variant_transcript_table = base_vat_table_name + "_variants"
    String genes_table = base_vat_table_name + "_genes"

    String vt_path = output_path + 'vt/*'
    String genes_path = output_path + 'genes/*'

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{vt_schema}
        fi

        echo "Loading data into a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        echo ~{vt_path}
        echo ~{genes_path}
        bq --apilog=false load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON ~{dataset_name}.~{variant_transcript_table} ~{vt_path}

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{genes_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
        fi

        echo "Loading data into a pre-vat table ~{dataset_name}.~{genes_table}"
        bq --apilog=false load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_path}

        set +e
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{vat_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the vat table ~{dataset_name}.~{vat_table_name}"
        else
            echo "Dropping and recreating the vat table ~{dataset_name}.~{vat_table_name}"
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table_name}
        fi

        CLUSTERING_STRING="--clustering_fields=contig"
        bq --apilog=false mk ${CLUSTERING_STRING} --expiration=$DATE --project_id=~{project_id} ~{dataset_name}.~{vat_table_name} ~{nirvana_schema}
        echo "Loading data into it"


        # Now we run a giant query in BQ to get this all in the right table and join the genes properly
        # Note the genes table join includes the group by to avoid the duplicates that get created from genes that span shards
        # Commented out columns in the query are to be added in the next release
        # We want the vat creation query to overwrite the destination table because if new data has been put into the pre-vat tables
        # and this workflow has been run an additional time, we dont want duplicates being appended from the original run

        # bq query --max_rows check: ok selecting into a table
        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{vat_table_name} --replace --project_id=~{project_id} \
        'SELECT
            v.vid,
            v.transcript,
            v.contig,
            v.position,
            v.ref_allele,
            v.alt_allele,
            v.gvs_all_ac,
            v.gvs_all_an,
            v.gvs_all_af,
            v.gvs_all_sc,
            v.gvs_max_af,
            v.gvs_max_ac,
            v.gvs_max_an,
            v.gvs_max_sc,
            v.gvs_max_subpop,
            v.gvs_afr_ac,
            v.gvs_afr_an,
            v.gvs_afr_af,
            v.gvs_afr_sc,
            v.gvs_amr_ac,
            v.gvs_amr_an,
            v.gvs_amr_af,
            v.gvs_amr_sc,
            v.gvs_eas_ac,
            v.gvs_eas_an,
            v.gvs_eas_af,
            v.gvs_eas_sc,
            v.gvs_eur_ac,
            v.gvs_eur_an,
            v.gvs_eur_af,
            v.gvs_eur_sc,
            v.gvs_mid_ac,
            v.gvs_mid_an,
            v.gvs_mid_af,
            v.gvs_mid_sc,
            v.gvs_oth_ac,
            v.gvs_oth_an,
            v.gvs_oth_af,
            v.gvs_oth_sc,
            v.gvs_sas_ac,
            v.gvs_sas_an,
            v.gvs_sas_af,
            v.gvs_sas_sc,
            v.gene_symbol,
            v.transcript_source,
            v.aa_change,
            v.consequence,
            v.dna_change_in_transcript,
            v.variant_type,
            v.exon_number,
            v.intron_number,
            v.genomic_location,
            # v.hgvsc AS splice_distance
            v.dbsnp_rsid,
            v.gene_id,
            # v.entrez_gene_id,
            # g.hgnc_gene_id,
            g.gene_omim_id,
            CASE WHEN ( v.transcript is not null and v.is_canonical_transcript is not True)
            THEN False WHEN ( v.transcript is not null and v.is_canonical_transcript is True) THEN True END AS is_canonical_transcript,
            v.gnomad_all_af,
            v.gnomad_all_ac,
            v.gnomad_all_an,
            v.gnomad_failed_filter,
            v.gnomad_max_af,
            v.gnomad_max_ac,
            v.gnomad_max_an,
            v.gnomad_max_subpop,
            v.gnomad_afr_ac,
            v.gnomad_afr_an,
            v.gnomad_afr_af,
            v.gnomad_amr_ac,
            v.gnomad_amr_an,
            v.gnomad_amr_af,
            v.gnomad_asj_ac,
            v.gnomad_asj_an,
            v.gnomad_asj_af,
            v.gnomad_eas_ac,
            v.gnomad_eas_an,
            v.gnomad_eas_af,
            v.gnomad_fin_ac,
            v.gnomad_fin_an,
            v.gnomad_fin_af,
            v.gnomad_nfe_ac,
            v.gnomad_nfe_an,
            v.gnomad_nfe_af,
            v.gnomad_sas_ac,
            v.gnomad_sas_an,
            v.gnomad_sas_af,
            v.gnomad_oth_ac,
            v.gnomad_oth_an,
            v.gnomad_oth_af,
            v.revel,
            v.splice_ai_acceptor_gain_score,
            v.splice_ai_acceptor_gain_distance,
            v.splice_ai_acceptor_loss_score,
            v.splice_ai_acceptor_loss_distance,
            v.splice_ai_donor_gain_score,
            v.splice_ai_donor_gain_distance,
            v.splice_ai_donor_loss_score,
            v.splice_ai_donor_loss_distance,
            g.omim_phenotypes_id,
            g.omim_phenotypes_name,
            v.clinvar_classification,
            v.clinvar_last_updated,
            v.clinvar_phenotype,
        FROM `~{dataset_name}.~{variant_transcript_table}` as v
            left join
        (SELECT gene_symbol, ANY_VALUE(gene_omim_id) AS gene_omim_id, ANY_VALUE(omim_phenotypes_id) AS omim_phenotypes_id, ANY_VALUE(omim_phenotypes_name) AS omim_phenotypes_name FROM `~{dataset_name}.~{genes_table}` group by gene_symbol) as g
        on v.gene_symbol = g.gene_symbol'
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    output {
        String vat_table = vat_table_name
        Boolean done = true
    }
}

task DeduplicateVatInBigQuery {
    meta {
        # since the WDL will not see the updated data (it's getting put in a gcp bucket)
        volatile: true
    }

    input {
        String input_vat_table_name
        String output_vat_table_name
        File? nirvana_schema

        String project_id
        String dataset_name
        String cloud_sdk_docker
    }


    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +e
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the final vat table ~{dataset_name}.~{output_vat_table_name}"
        else
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name}
        fi
        bq --apilog=false mk --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name} ~{nirvana_schema}
        echo "Loading data into it"

        # Now we query the original VAT table and recreate it, but remove any rows that appear twice.

        # bq query --max_rows check: ok selecting into a table
        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{output_vat_table_name} --replace --project_id=~{project_id} \
        ' SELECT * EXCEPT(row_number) FROM (
            SELECT
                *,
                row_number()
                    over (partition by vid, transcript)
                    row_number
            FROM
                `~{dataset_name}.~{input_vat_table_name}`
            )
            where row_number = 1'
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    output {
        String vat_table = output_vat_table_name
        Boolean done = true
    }
}
