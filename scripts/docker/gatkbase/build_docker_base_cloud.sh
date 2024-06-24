#!/bin/bash
#
# A script that builds the GATK base image using Google Cloud Build, and pushes it to
# a staging location at us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area
#
# Usage: build_docker_base_cloud.sh <docker_image_version>
#
# After staging the image, you should test it with GATK before actually releasing it
# using the release_prebuilt_base_image.sh script. You can test it by modifying the
# main GATK Dockerfile FROM clause to point temporarily at the staged base image, and
# submitting that as a PR to trigger a test suite run.
#

if [ $# -ne 1 ]; then
    echo "Usage: $0 <docker_image_version>"
    exit 1
fi

IMAGE_VERSION=$1
IMAGE_NAME="us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area"
DOCKER_IMAGE_TAG="${IMAGE_NAME}:gatkbase-${IMAGE_VERSION}"

gcloud builds submit --tag ${DOCKER_IMAGE_TAG} --timeout=24h --machine-type n1_highcpu_32

if [ $? -ne 0 ]; then
    echo "gcloud builds submit failed"
    exit 1
fi

echo "Successfully published image to staging area at ${DOCKER_IMAGE_TAG}"
exit 0