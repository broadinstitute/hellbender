version 1.0

task build_and_push {

    input {
        String repo_url
        String software_name
        String commit_hash
        String registry_host
    }

    parameter_meta {
        repo_url: "The url of the repository containing the project to build from"
        software_name: "The name that will be used to name the docker image"
        commit_hash: "The hash for the commit to build from; will also be used to tag the image"
        registry_host: "The docker repository to push the image to"
    }

    command {
        mkdir repo-folder
        cd repo-folder
        git clone ${repo_url} .
        git checkout ${commit_hash}
        echo -n "" >> .gcloudignore
        gcloud builds submit --tag ${registry_host}/${software_name}:${commit_hash} --timeout=24h
    }

    runtime {
        docker: "google/cloud-sdk:307.0.0"
    }

}

workflow docker_build {

    meta {
        description: "Build a docker image for a public git repository at a specified commit using GCLoud Build and push it to a GCR repository"
    }

    input {
        String repo_url = "https://github.com/broadinstitute/gatk.git"
        String software_name = "gatk"
        String commit_hash
        String registry_host
    }

    parameter_meta {
        repo_url: "The url of the repository containing the project to build from"
        software_name: "The name that will be used to name the docker image"
        commit_hash: "The hash for the commit to build from; will also be used to tag the image"
        registry_host: "The docker repository to push the image to"
    }

    call build_and_push {
        input:
            repo_url = repo_url,
            software_name = software_name,
            commit_hash = commit_hash,
            registry_host = registry_host
    }

    output {
        String image_url = registry_host + '/' + software_name + ':' + commit_hash
    }

}
