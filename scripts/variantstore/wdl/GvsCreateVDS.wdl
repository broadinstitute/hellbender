version 1.0

# This WDL will create a VDS in Hail running in a Dataproc cluster.
import "GvsUtils.wdl" as Utils


workflow GvsCreateVDS {
    input {
        String avro_path
        String vds_destination_path

        String cluster_prefix = "vds-cluster"
        String gcs_subnetwork_name = "subnetwork"
        String? hail_temp_path
        Int? intermediate_resume_point
        String region = "us-central1"

        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end = false
        Float? master_memory_fraction
        Boolean use_classic_VQSR = false

        String? git_branch_or_tag
        String? hail_version
        File? hail_wheel
        String? basic_docker
        String? variants_docker
        String? workspace_bucket
        String? workspace_project
        String? cloud_sdk_slim_docker
    }

    parameter_meta {
        avro_path : {
            help: "Input location for the avro files"
        }
        vds_destination_path: {
            help: "Location for the final created VDS"
        }
        cluster_prefix: {
            help: "Prefix of the Dataproc cluster name"
        }
        gcs_subnetwork_name: {
            help: "Set to 'subnetwork' if running in Terra Cromwell"
        }
        hail_temp_path: {
            help: "Hail temp path to use, specify if resuming from a run that failed midway through creating intermediate VDSes."
        }
        intermediate_resume_point: {
            help: "Index at which to resume creating intermediate VDSes."
        }
        region: {
            help: "us-central1"
        }
        hail_version: {
            help: "Optional Hail version, defaults to 0.2.126. Cannot define both this parameter and `hail_wheel`."
        }
        hail_wheel: {
            help: "Optional Hail wheel. Cannot define both this parameter and `hail_version`."
        }
    }


    if (!defined(variants_docker) || !defined(basic_docker) || !defined(cloud_sdk_slim_docker) || !defined(workspace_bucket) || !defined(workspace_project) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_google_project = select_first([workspace_project, GetToolVersions.google_project])
    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])

    if (defined(hail_version) && defined(hail_wheel)) {
        call Utils.TerminateWorkflow as BothHailVersionAndHailWheelDefined {
            input:
                message = "Cannot define both `hail_version` and `hail_wheel`, exiting.",
                basic_docker = effective_basic_docker,
        }
    }

    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    if (defined(intermediate_resume_point) && !defined(hail_temp_path)) {
        call Utils.TerminateWorkflow as NeedHailTempPath {
            input:
                message = "GvsCreateVDS called with an intermediate resume point but no specified hail temp path from which to resume",
                basic_docker = effective_basic_docker,
        }
    }

    if (!defined(intermediate_resume_point) && defined(hail_temp_path)) {
        call Utils.TerminateWorkflow as NeedIntermediateResumePoint {
            input:
                message = "GvsCreateVDS called with no intermediate resume point but a specified hail temp path, which isn't a known use case",
                basic_docker = effective_basic_docker,
        }
    }

    call GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call CreateVds {
        input:
            prefix = cluster_prefix,
            vds_path = vds_destination_path,
            avro_path = avro_path,
            use_classic_VQSR = use_classic_VQSR,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            hail_temp_path = hail_temp_path,
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            gvs_import_script = GetHailScripts.gvs_import_script,
            hail_gvs_import_script = GetHailScripts.hail_gvs_import_script,
            intermediate_resume_point = intermediate_resume_point,
            workspace_project = effective_google_project,
            region = region,
            workspace_bucket = effective_workspace_bucket,
            gcs_subnetwork_name = gcs_subnetwork_name,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            master_memory_fraction = master_memory_fraction,
    }

    call ValidateVds {
        input:
            go = CreateVds.done,
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            vds_validation_script = GetHailScripts.vds_validation_script,
            prefix = cluster_prefix,
            vds_path = vds_destination_path,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            workspace_project = effective_google_project,
            region = region,
            workspace_bucket = effective_workspace_bucket,
            gcs_subnetwork_name = gcs_subnetwork_name,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            master_memory_fraction = master_memory_fraction,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
    }

    output {
        String create_cluster_name = CreateVds.cluster_name
        String validate_cluster_name = ValidateVds.cluster_name
        Boolean done = true
    }
}

task CreateVds {
    input {
        String prefix
        String vds_path
        String avro_path
        Boolean use_classic_VQSR
        Boolean leave_cluster_running_at_end
        File hail_gvs_import_script
        File gvs_import_script
        File run_in_hail_cluster_script
        String? hail_version
        File? hail_wheel
        String? hail_temp_path
        Int? intermediate_resume_point
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction

        String workspace_project
        String workspace_bucket
        String region
        String gcs_subnetwork_name

        String cloud_sdk_slim_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        pip3 install --upgrade pip

        if [[ ! -z "~{hail_wheel}" ]]
        then
            pip3 install ~{hail_wheel}
        else
            pip3 install hail~{'==' + hail_version}
        fi

        pip3 install --upgrade google-cloud-dataproc ijson

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | od -h -An | tr -d '[:space:]')-$(head -c2 < /dev/urandom | od -h -An | tr -d '[:space:]')"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        if [[ -z "~{hail_temp_path}" ]]
        then
            hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"
        else
            hail_temp_path="~{hail_temp_path}"
        fi

        # Set up the autoscaling policy
        cat > auto-scale-policy.yaml <<FIN
        workerConfig:
            minInstances: 2
            maxInstances: 2
        secondaryWorkerConfig:
            maxInstances: 500
        basicAlgorithm:
            cooldownPeriod: 120s
            yarnConfig:
                scaleUpFactor: 1.0
                scaleDownFactor: 1.0
                gracefulDecommissionTimeout: 120s
        FIN
        gcloud dataproc autoscaling-policies import gvs-autoscaling-policy --project=~{workspace_project} --source=auto-scale-policy.yaml --region=~{region} --quiet

        # construct a JSON of arguments for python script to be run in the hail cluster
        cat > script-arguments.json <<FIN
        {
            "vds-path": "~{vds_path}",
            "temp-path": "${hail_temp_path}",
            "avro-path": "~{avro_path}"
            ~{', "intermediate-resume-point": ' + intermediate_resume_point}
            ~{true=', "use-classic-vqsr": ""' false='' use_classic_VQSR}
        }
        FIN

        # Run the hail python script to make a VDS
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{hail_gvs_import_script} \
            --secondary-script-path-list ~{gvs_import_script} \
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
        docker: cloud_sdk_slim_docker
        bootDiskSizeGb: 10
    }

    output {
        String cluster_name = read_string("cluster_name.txt")
        Boolean done = true
    }
}

task ValidateVds {
    input {
        Boolean go
        Boolean leave_cluster_running_at_end
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction
        File run_in_hail_cluster_script
        File vds_validation_script
        String prefix
        String vds_path
        String? hail_version
        File? hail_wheel
        String workspace_project
        String workspace_bucket
        String region
        String gcs_subnetwork_name
        String cloud_sdk_slim_docker
    }

    meta {
        # should always be run
        volatile: true
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        pip3 install --upgrade pip

        if [[ ! -z "~{hail_wheel}" ]]
        then
            pip3 install ~{hail_wheel}
        else
            pip3 install hail~{'==' + hail_version}
        fi

        pip3 install --upgrade google-cloud-dataproc ijson

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | od -h -An | tr -d '[:space:]')-$(head -c2 < /dev/urandom | od -h -An | tr -d '[:space:]')"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt
        hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"

        # construct a JSON of arguments for python script to be run in the hail cluster
        cat > script-arguments.json <<FIN
        {
            "vds-path": "~{vds_path}",
            "temp-path": "${hail_temp_path}"
        }
        FIN

        # Run the hail python script to validate a VDS
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{vds_validation_script} \
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
        docker: cloud_sdk_slim_docker
        bootDiskSizeGb: 10
    }
    output {
        String cluster_name = read_string("cluster_name.txt")
    }
}

task GetHailScripts {
    input {
        String variants_docker
    }
    meta {
        # OK to cache this as the scripts are drawn from the stringified Docker image and as long as that stays the same
        # the script content should also stay the same.
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Not sure why this is required but without it:
        # Absolute path /app/run_in_hail_cluster.py doesn't appear to be under any mount points: local-disk 10 SSD

        mkdir app
        cp /app/*.py app
    >>>
    output {
        File run_in_hail_cluster_script = "app/run_in_hail_cluster.py"
        File hail_gvs_import_script = "app/hail_gvs_import.py"
        File gvs_import_script = "app/import_gvs.py"
        File vds_validation_script = "app/vds_validation.py"
    }
    runtime {
        docker: variants_docker
    }
}
