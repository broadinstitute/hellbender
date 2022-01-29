WF_ID=$1
TOKEN=$(gcloud auth print-access-token)

PARAMS="includeKey=start&includeKey=end&includeKey=executionStatus&includeKey=preemptible&includeKey=shardIndex&includeKey=wdl-task-name&includeKey=jes&includeKey=jobId&includeKey=executionEvents&includeKey=outputs"

# for BIG workflows where we are going to get PAPI info for start/stop times anyway
PARAMS="includeKey=start&includeKey=end&includeKey=executionStatus&includeKey=preemptible&includeKey=shardIndex&includeKey=wdl-task-name&includeKey=jes&includeKey=jobId"

curl -sL -X 'GET' \
  "https://api.firecloud.org/api/workflows/v1/${WF_ID}/metadata?${PARAMS}" \
  -H 'accept: application/json' \
  -H "Authorization: Bearer ${TOKEN}"
