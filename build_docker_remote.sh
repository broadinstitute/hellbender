#!/usr/bin/env bash
#
# Build (and optionally push) a GATK docker image to GCR using Google Cloud Build. Images are built
# in the cloud rather than locally. Pushing to dockerhub is not supported by this script.
#
# By default the images are pushed to the following GCR repository:
#
# us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds
#
# with a name like "${YOUR_USERNAME}-${GITHUB_TAG}-${GIT_HASH_FOR_TAG}"
#
# This script should not be used to push to our release repositories. After
# you've built and staged an image, run the scripts/docker/release_prebuilt_docker_image.sh
# script to push to the release repositories.
#
# Usage: build_docker_remote.sh -e <GITHUB_TAG> -d <STAGING_DIRECTORY> [-str]
#
# where <GITHUB_TAG> is the github tag (or hash when -s is used) to use in building the docker image (e.g. 4.2.6.1)
# and <STAGING_DIRECTORY> is a directory in which to clone the repo and stage the build (DO NOT SPECIFY YOUR WORKING DIRECTORY)
# Optional arguments:
#   -s  The GITHUB_TAG (-e parameter) is actually a github hash, not tag.
#   -r  Build this image with the release flag set to true, causing the version number to not end with SNAPSHOT
#   -t <IMAGE_TAG>  The tag to assign image once it is finished constructing.  NOTE: currently this MUST be on either GCR or the Google Artifact Registry.
#

# Have script stop if there is an error
set -e

REPO=broadinstitute
PROJECT=gatk
GCR_REPO="us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds"
STAGING_CLONE_DIR=${PROJECT}_staging_temp

#################################################
# Parsing arguments
#################################################
while getopts "e:sd:t:r" option; do
	case "$option" in
		e) GITHUB_TAG="$OPTARG" ;;
		s) IS_HASH=true ;;
		d) STAGING_DIR="$OPTARG" ;;
		t) DOCKER_IMAGE_TAG="$OPTARG" ;;
		r) RELEASE=true ;;
	esac
done

function usage() {
  MESSAGE=$1
	printf "%s\n \
Usage: build_docker_remote.sh -e <GITHUB_TAG> -d <STAGING_DIRECTORY> [-str] \n \
where <GITHUB_TAG> is the github tag (or hash when -s is used) to use in building the docker image (e.g. 4.2.6.1)\n \
and <STAGING_DIRECTORY> is a directory in which to clone the repo and stage the build (DO NOT SPECIFY YOUR WORKING DIRECTORY)\n \
Optional arguments:  \n \
-s \t The GITHUB_TAG (-e parameter) is actually a github hash, not tag.  \n \
-r \t Build this image with the release flag set to true, causing the version number to not end with SNAPSHOT \n \
-t <IMAGE_TAG>\t  The tag to assign image once it is finished constructing.  NOTE: currently this MUST be on either GCR or the Google Artifact Registry. \n" "$MESSAGE"
}

if [ -z "$GITHUB_TAG" ]; then
	usage "Option -e (github tag) requires an argument."
	exit 1
fi

if [ -z "$STAGING_DIR" ]; then
	usage "Option -d (staging directory) requires an argument."
	exit 1
fi

# Output the parameters
echo -e "\n"
echo -e "github tag/hash: ${GITHUB_TAG}"
echo -e "github project: ${REPO}/${PROJECT}:${GITHUB_TAG}\n\n"
echo "Other options (Blank is false)"
echo "---------------"
echo "This is a git hash: ${IS_HASH}"
echo "Staging directory: ${STAGING_DIR}"
echo "Release: ${RELEASE}"

ORIGINAL_WORKING_DIRECTORY=$(pwd)

if [ -n "$STAGING_DIR" ]; then
    GITHUB_DIR="tags/"

    if [ -n "${IS_HASH}" ]; then
        GITHUB_DIR=" "
    fi

    mkdir -p ${STAGING_DIR}
    cd ${STAGING_DIR}
    set +e
    rm -Rf ${STAGING_CLONE_DIR}
    set -e
    GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/${REPO}/${PROJECT}.git ${STAGING_CLONE_DIR}
    cd ${STAGING_CLONE_DIR}
    STAGING_ABSOLUTE_PATH=$(pwd)

    echo "Now in $(pwd)"
    GIT_CHECKOUT_COMMAND="git checkout ${GITHUB_DIR}${GITHUB_TAG}"
    echo "${GIT_CHECKOUT_COMMAND}"
    ${GIT_CHECKOUT_COMMAND}
fi

GIT_HASH_FOR_TAG=$(git describe --tags --long)

## generate the tag if it wasn't explicitly specified
if [ -z "$DOCKER_IMAGE_TAG" ]; then
  DOCKER_IMAGE_TAG=${GCR_REPO}:$(whoami)-${GITHUB_TAG}-${GIT_HASH_FOR_TAG}
fi

echo "steps:" >> cloudbuild.yaml
echo "- name: 'gcr.io/cloud-builders/docker'" >> cloudbuild.yaml
if [ -n "$RELEASE" ]; then
  echo "  args: [ 'build', '-t', '${DOCKER_IMAGE_TAG}', '--build-arg', 'RELEASE=true', '.' ]" >> cloudbuild.yaml
else
  echo "  args: [ 'build', '-t', '${DOCKER_IMAGE_TAG}', '.' ]" >> cloudbuild.yaml
fi
echo "- name: 'gcr.io/cloud-builders/docker'" >> cloudbuild.yaml
echo "  args: [ 'push', '${DOCKER_IMAGE_TAG}' ]" >> cloudbuild.yaml

echo "Building image with the tag ${DOCKER_IMAGE_TAG}..."

SUBMIT_COMMAND="gcloud builds submit --config cloudbuild.yaml --timeout=24h --machine-type n1_highcpu_32"
echo "running the following gcloud command: ${SUBMIT_COMMAND}"
## We need to override the default .gcloudignore to preserve the .git directory which we need in order to download LFS files in the remote build.
echo -n "" >> .gcloudignore
${SUBMIT_COMMAND}

echo "Image successfully built and pushed to ${DOCKER_IMAGE_TAG}"

cd ${ORIGINAL_WORKING_DIRECTORY}
if [ -n "$STAGING_DIR" ] ; then
    rm -Rf ${STAGING_DIR}/${STAGING_CLONE_DIR}
fi

exit 0