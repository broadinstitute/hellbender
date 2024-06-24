# How to update the GATK base docker image:

1. In a branch, make whatever updates you need to the base image Dockerfile in `scripts/docker/gatkbase`
2. Choose a new version number for the new base image
3. Build the GATK base image by running either `build_docker_base_cloud.sh` (to build in the cloud using Google Cloud Build) or `build_docker_base_locally.sh` (to build on your local machine) from within the `scripts/docker/gatkbase` directory in your GATK clone.
   * Both scripts take the version of the new image as their only argument. Eg., `build_docker_base_cloud.sh 3.0.0rc1` will create the remote image `us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area:gatkbase-3.0.0rc1`
4. If you built the image locally, you'll need to manually push it to the staging area at `us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area`
5. Test the new image using the GATK branch with your changes to the base image Dockerfile, by modifying the FROM clause at the top of the main GATK Dockerfile to point to the base image you staged, and submit a PR for the branch to let the test suite run.
6. Once tests pass and everything looks good, push the new base image to the release repositories using the `release_prebuilt_base_image.sh` script
   * The release script takes as arguments the prebuilt image you've tested, and a final version number for release. For example: `release_prebuilt_base_image.sh us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area:gatkbase-3.0.0rc1 3.0.0`
   * The release script looks for the image locally first, and if not found locally it pulls it from the remote repository before releasing.
7. Update the FROM clause in the main GATK Dockerfile in your GATK PR to point to the officially-released image, and merge it once tests pass.



