#!/bin/bash
#
# A script that builds the GATK base image locally (but does not push it anywhere).
#
# Usage: build_docker_base_locally.sh <docker_image_version>
#
# After building the image, you should test it with GATK before actually releasing it
# using the release_prebuilt_base_image.sh script. You can test it by modifying the
# main GATK Dockerfile FROM clause to point temporarily at the candidate base image, and
# submitting that as a PR to trigger a test suite run.
#

if [ $# -ne 1 ]; then
    echo "Usage: $0 <docker_image_version>"
    exit 1
fi

IMAGE_VERSION=$1
IMAGE_REPO="us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area"
IMAGE_FULL_TAG="${IMAGE_REPO}:gatkbase-${IMAGE_VERSION}"

# Build
echo "Building image to tag ${IMAGE_FULL_TAG}..."
docker build --squash -t "${IMAGE_FULL_TAG}" .

if [ $? -ne 0 ]; then
  echo "docker build failed"
  exit 1
fi

echo "docker build succeeded"
exit 0

