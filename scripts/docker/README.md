Docker for GATK4 (Developers)
--------------------------------------

*This docker script is meant for the GATK4 Broad dev team to create and push docker images.*

*Users should get gatk docker images from dockerhub (https://hub.docker.com/r/broadinstitute/gatk/)*

*All docker images should be created with the appropriate build script (see below)*

This doc does not give usage info on Docker and assumes you are already versed in its use.

Please see the excellent docker documentation if you need info about docker:  https://docs.docker.com/

Notes:
- HDF5 jni library is in ``/usr/lib/jni``, since this is the default for Ubuntu 16.04.  This has been added to the ``JAVA_LIBRARY_PATH`` environment variable.


## Overview

This repo contains the scripts for creating and pushing two docker images:
- gatkbase -- a basic docker image that we do not expect to change very often.  The GATK4 docker image uses this one (``FROM``)
- gatk4 -- the official docker image for GATK4.  The instructions in this document pertain to this image, unless otherwise stated.

``scripts/docker/gatkbase/build_docker_base_cloud.sh`` and ``scripts/docker/gatkbase/build_docker_base_locally.sh`` are scripts to create the gatkbase docker image. 
``build_docker.sh`` is a script to create the full gatk4 docker image.
``build_docker_remote.sh`` is a script to create the full gatk4 docker image using google cloud build remotely. This is useful if you can't build the docker image locally (for example if you have an M1 Macbook) NOTE: this requires the user first specify their project with the command `gcloud config set project <PROJECT>` to a project that has access to google cloud build. 

## GATK4 Docker image

This is the image where the actual GATK4 resides.  Location:  ``/gatk/gatk.jar``

#### Create GATK4 docker image and push the image to the gatk dockerhub (Seriously, only for Broad GATK4 dev team)

This allows you to create the "official" GATK4 docker image and push it to docker hub (if you have access) 

*Please allow 1 hour for completion*  This includes running the unit tests inside the docker image. 

```bash

# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_TAG=1.0.0.0-alpha1.2.1

# REPLACE with the directory where you'd like to clone the repo
export STAGING_DIR=~/tmp/tmp_build_docker_image/

sudo bash build_docker.sh -e ${GITHUB_TAG} -p -d ${STAGING_DIR}
```

#### Create GATK4 docker image and upload into GCR

This will attempt to build the docker remotely using google cloud build and tag it with the following this format: `us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:<USERNAME>-<GITHUB_TAG>-<GIT_HASH>`

From this directory, run:

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_TAG=4.2.6.1

# REPLACE with the directory where you'd like to clone the repo
export STAGING_DIR=~/tmp/tmp_build_docker_image/

sudo bash build_docker_remote.sh -e ${GITHUB_TAG} -d ${STAGING_DIR}

```

#### Create GATK4 docker image and upload into GCR (with a custom gcr repo tag)

This will attempt to build the docker remotely using google cloud build. From this directory, run:

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_TAG=4.2.6.1

# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GCR_REPO_TAG=us.gcr.io/myrepository:myname-master-4.2.6.1

# REPLACE with the directory where you'd like to clone the repo
export STAGING_DIR=~/tmp/tmp_build_docker_image/

sudo bash build_docker_remote.sh -e ${GITHUB_TAG} -t ${GCR_REPO_TAG} -d ${STAGING_DIR}

```

#### Create GATK4 docker image (do not push to dockerhub)

From this directory, run:

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_TAG=1.0.0.0-alpha1.2.1

# REPLACE with the directory where you'd like to clone the repo
export STAGING_DIR=~/tmp/tmp_build_docker_image/

sudo bash build_docker.sh -e ${GITHUB_TAG} -d ${STAGING_DIR}

```

#### Create GATK4 docker image from a github hash (pushing to dockerhub prohibited)

From this directory, run:

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_HASH=e454ac88c7791c0f8b385b3e82138ec52c61ef48

# REPLACE with the directory where you'd like to clone the repo
export STAGING_DIR=~/tmp/tmp_build_docker_image/

sudo bash build_docker.sh -e ${GITHUB_HASH} -s -d ${STAGING_DIR}
```

Alternatively, build a docker image based on the latest *push* of your current branch:
```bash
# REPLACE with the directory where you'd like to clone the repo
export STAGING_DIR=~/tmp/tmp_build_docker_image/

BRANCH=`git rev-parse --symbolic-full-name --abbrev-ref HEAD` 
sudo bash build_docker.sh -e `git rev-parse ${BRANCH}` -s -r -d ${STAGING_DIR}
```

```bash
# Same as above except do not run the unit tests.
sudo bash build_docker.sh -e `git rev-parse ${BRANCH}` -s -u -d ${STAGING_DIR}
```

#### Run GATK4 gradle tests in existing latest tagged docker image

*Note that the unit tests are run during a build of the Dockerfile*

Currently, this is a bit manual.

```bash
# bash is the default command for the docker image.
sudo docker pull broadinstitute/gatk
sudo docker run -i -t broadinstitute/gatk:latest

# On the docker prompt
cd /root/gatk
./gradlew test

# To leave the docker prompt:
exit
```

## gatkbase Docker image

This is a base image that does not require any files in the gatk repo (except the build script and Dockerfile, obviously).  GATK docker images are dependent on this one.

**IMPORTANT** 
- The gatkbase build scripts should be run from the ``scripts/docker/gatkbase`` directory.

#### Build the gatkbase docker image on your local machine:

```bash
build_docker_base_locally.sh <docker_image_version>
```

#### Build the gatkbase docker image remotely using Google Cloud Build:

```bash
build_docker_base_cloud.sh <docker_image_version>
```

#### Release a pre-built gatkbase image to the official repositories, after testing it:
 
```bash
release_prebuilt_base_image.sh <prebuilt_image> <version_number_for_release>
```


## Other useful Info
#### Delete all containers
```bash
# Credit:  http://stackoverflow.com/questions/17236796/how-to-remove-old-docker-containers
sudo docker rm `sudo docker ps --no-trunc -aq`
```