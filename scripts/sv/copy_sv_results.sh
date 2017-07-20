#!/bin/bash

# This script copies SV analysis results from a Google Dataproc cluster
# to an appropriate bucket/directory on GCS. It also uploads contents of
# local output logs

# terminate script on error or if a command fails before piping to another command
set -eu
set -o pipefail

if [[ "$#" -lt 3 ]]; then
    echo -e
"Please provide:
  [1] GCS project name (required)
  [2] GCS cluster name (required)
  [3] cluster output directory (required)
  [4] GCS user name (defaults to local user name)
  [5] path to local log file (default to empty, i.e. no log)
  [*] additional arguments that were passed to
      StructuralVariationDiscoveryPipelineSpark"
  exit 1
fi

PROJECT_NAME=$1
CLUSTER_NAME=$2
OUTPUT_DIR=$3
GCS_USER=${4:-${USER}}
LOCAL_LOG_FILE=${5:-"/dev/null"}

shift 5
SV_ARGS=${*:-${SV_ARGS:-""}}

# get appropriate ZONE for cluster
echo "CLUSTER_INFO=\$(gcloud dataproc clusters list --project=${PROJECT_NAME} --filter='clusterName=${CLUSTER_NAME}')"
CLUSTER_INFO=$(gcloud dataproc clusters list --project=${PROJECT_NAME} --filter="clusterName=${CLUSTER_NAME}" | tail -n 1 | tr -s ' ')
ZONE=$(echo "${CLUSTER_INFO}" | cut -d' ' -f 4)
if [ -z "${ZONE}" ]; then
    # cluster is down.
    echo "Cluster \"${CLUSTER_NAME}\" is down. Only log and command args will be uploaded"
    RESULTS_DIR=""
else
    # get the latest time-stamped results directory from the cluster
    # (may not be current date stamp if multiple jobs run on same cluster)
    MASTER="${CLUSTER_NAME}-m"
    RESULTS_DIR="$(dirname ${OUTPUT_DIR})"
    RESULTS_DIR=$(gcloud compute ssh ${MASTER} --zone ${ZONE} --command="hadoop fs -ls ${RESULTS_DIR} | tail -n 1")
    RESULTS_DIR=$(echo "${RESULTS_DIR}" | awk '{print $NF}' | sed -e 's/^\///')
fi

if [ -z "${RESULTS_DIR}" ]; then
    # directory is empty (presumably the job crashed). Use OUTPUT_DIR (without leading slash)
    RESULTS_DIR=$(echo "${OUTPUT_DIR}" | sed -e 's/^\///')
    echo "RESULTS_DIR=${RESULTS_DIR}" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    GCS_RESULTS_DIR="gs://${PROJECT_NAME}/${GCS_USER}/${RESULTS_DIR}"
else
    # copy the latest results to google cloud
    echo "RESULTS_DIR=${RESULTS_DIR}" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    GCS_RESULTS_DIR="gs://${PROJECT_NAME}/${GCS_USER}/${RESULTS_DIR}"
    # chose semi-optimal parallel args for distcp
    # 1) count number of files to copy
    COUNT_FILES_CMD="hadoop fs -count /${RESULTS_DIR}/ | tr -s ' ' | cut -d ' ' -f 3"
    NUM_FILES=$(gcloud compute ssh ${MASTER} --zone ${ZONE} --command="${COUNT_FILES_CMD}")
    # 2) get the number of instances
    NUM_INSTANCES=$(echo "${CLUSTER_INFO}" | cut -d' ' -f 2)
    # 3) choose number of maps as min of NUM_FILES or NUM_INSTANCES * MAPS_PER_INSTANCE
    MAPS_PER_INSTANCE=10
    NUM_MAPS=$((${NUM_INSTANCES} * ${MAPS_PER_INSTANCE}))
    NUM_MAPS=$((${NUM_FILES} < ${NUM_MAPS} ? ${NUM_FILES} : ${NUM_MAPS}))
    # 4) set arguments to distcp to prevent it from erroring out due to too many maps
    SPLIT_RATIO=$((${NUM_FILES}/${NUM_MAPS}))
    SPLIT_RATIO=$((${SPLIT_RATIO} < 3 ? ${SPLIT_RATIO} : 3))
    MAX_CHUNKS_IDEAL=$((${NUM_MAPS}*${SPLIT_RATIO}))
    MAX_CHUNKS_TOL=$((${MAX_CHUNKS_IDEAL} + ${NUM_MAPS}))
    DIST_CP_ARGS="-D distcp.dynamic.max.chunks.tolerable=${MAX_CHUNKS_TOL} -D distcp.dynamic.max.chunks.ideal=${MAX_CHUNKS_IDEAL} -D distcp.dynamic..min.records_per_chunk=0 -D distcp.dynamic.split.ratio=${SPLIT_RATIO}"
    CPY_CMD="hadoop distcp ${DIST_CP_ARGS} -m ${NUM_MAPS} -strategy dynamic /${RESULTS_DIR}/* ${GCS_RESULTS_DIR}/"
    echo "gcloud compute ssh ${MASTER} --zone ${ZONE} --command=\"${CPY_CMD}\" 2>&1 | tee -a ${LOCAL_LOG_FILE}" | tee -a ${LOCAL_LOG_FILE}
    gcloud compute ssh ${MASTER} --zone ${ZONE} --command="${CPY_CMD}" 2>&1 | tee -a ${LOCAL_LOG_FILE}
fi

# create file with command-line args. Always create (even if empty) so
# that it's unambiguous if args were logged or not (file is present),
# and if args were passed or not (file is non-empty)
echo "gsutil cp >(echo \"${SV_ARGS}\") \"${GCS_RESULTS_DIR}/sv-command-line-args.txt\"" | tee -a ${LOCAL_LOG_FILE}
gsutil cp <(echo "${SV_ARGS}") "${GCS_RESULTS_DIR}/sv-command-line-args.txt"

if [ "${LOCAL_LOG_FILE}" != "/dev/null" ]; then
    # copy log to google cloud
    REMOTE_LOG_FILE="${GCS_RESULTS_DIR}/sv-discovery.log"
    if gsutil -q stat ${REMOTE_LOG_FILE}; then
        # remote log file exists, append local to remote
        gsutil mv ${LOCAL_LOG_FILE} "${REMOTE_LOG_FILE}-append"
        gsutil compose ${REMOTE_LOG_FILE} "${REMOTE_LOG_FILE}-append"
        gsutil rm "${REMOTE_LOG_FILE}-append"
    else
        # mv log file from local to bucket
        gsutil mv ${LOCAL_LOG_FILE} ${REMOTE_LOG_FILE}
    fi
fi
