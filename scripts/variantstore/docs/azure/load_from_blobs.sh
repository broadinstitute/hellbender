# Prepend date, time and pwd to xtrace log entries.
PS4='\D{+%F %T} \w $ '
set -o errexit -o nounset -o pipefail -o xtrace

usage() {
  echo "Usage: $(basename "$0") --group <Resource group name> --server <SQL Server> --database <SQL Database> --account <Storage account name> --container <Storage container name> --sas <Storage container SAS token> --password <Master key password>" 1>&2
  exit 1
}

VALID_ARGS=$(getopt --options g:s:d:a:c:t:p: --long group:,server:,database:,account:,container:,sas-token:,password: -- "$@")
if [[ $? -ne 0 ]]; then
    usage
fi

eval set -- "$VALID_ARGS"
while true ; do
  case "$1" in
    -g|--group)
      RESOURCE_GROUP="${OPTARG}"
      shift 2
      ;;
    -s|--server)
      SQL_SERVER="${OPTARG}"
      shift 2
      ;;
    -d|--database)
      SQL_DATABASE="${OPTARG}"
      shift 2
      ;;
    -a|--account)
      STORAGE_ACCOUNT_NAME="${OPTARG}"
      shift 2
      ;;
    -c|--container)
      STORAGE_CONTAINER_NAME="${OPTARG}"
      shift 2
      ;;
    -t|--sas-token)
      CSV_CONTAINER_SAS_TOKEN="${OPTARG}"
      shift 2
      ;;
    -p|--password)
      MASTER_KEY_PASSWORD="${OPTARG}"
      shift 2
      ;;
    --) shift;
      break
      ;;
  esac
done


az group create --location eastus --resource-group "${RESOURCE_GROUP}"

# Create the server (yes serverless Azure SQL Database requires a server).
az sql server create \
    --resource-group "${RESOURCE_GROUP}" \
    --enable-ad-only-auth \
    --enable-public-network \
    --location eastus \
    --external-admin-principal-type User \
    --external-admin-sid "$(az ad signed-in-user show | jq -r .id)" \
    --external-admin-name "$(az ad signed-in-user show | jq -r .mail)" \
    --name "${SQL_SERVER}"

# Get the external IP of this machine and create a firewall rule that allows it to connect to the Azure SQL Database.
# https://dev.to/0xbf/get-your-public-ip-address-from-command-line-57l4
EXTERNAL_IP="$(curl --silent ifconfig.me)"
az sql server firewall-rule create \
    --resource-group "${RESOURCE_GROUP}" \
    --server "${SQL_SERVER}" \
    --name AllowYourIp \
    --start-ip-address "${EXTERNAL_IP}" \
    --end-ip-address "${EXTERNAL_IP}"

# Note the database is specified as 'Serverless', not the server.
az sql db create \
    --resource-group "${RESOURCE_GROUP}" \
    --server "${SQL_SERVER}" \
    --name "${SQL_DATABASE}" \
    -e GeneralPurpose \
    -f Gen5 \
    -c 2 \
    --auto-pause-delay 10 \
    --compute-model Serverless


LOADER_SCRIPT="loader_script.sql"

cat > "${LOADER_SCRIPT}" <<FIN


CREATE MASTER KEY ENCRYPTION BY PASSWORD = '${MASTER_KEY_PASSWORD}';


CREATE DATABASE SCOPED CREDENTIAL MyAzureBlobStorageCredential
WITH IDENTITY = 'SHARED ACCESS SIGNATURE',
SECRET = '${CSV_CONTAINER_SAS_TOKEN}';


CREATE EXTERNAL DATA SOURCE MyAzureBlobStorage
WITH (
    TYPE = BLOB_STORAGE,
    LOCATION = 'https://${STORAGE_ACCOUNT_NAME}.blob.core.windows.net/${STORAGE_CONTAINER_NAME}',
    CREDENTIAL= MyAzureBlobStorageCredential
);


CREATE TABLE ref_ranges (
    location BIGINT,
    sample_id INT,
    length SMALLINT,
    state TINYINT
);


CREATE TABLE vets (
    location BIGINT,
    sample_id INT,
    ref VARCHAR(1024),
    alt VARCHAR(4096),
    AS_RAW_MQ VARCHAR(255),
    AS_RAW_MQRankSum VARCHAR(255),
    QUALapprox VARCHAR(255),
    AS_QUALapprox VARCHAR(255),
    AS_RAW_ReadPosRankSum VARCHAR(255),
    AS_SB_TABLE VARCHAR(255),
    AS_VarDP VARCHAR(255),
    call_GT VARCHAR(255),
    call_AD VARCHAR(255),
    call_GQ VARCHAR(255),
    call_PGT VARCHAR(255),
    call_PID VARCHAR(1024),
    call_PL VARCHAR(255)
);


FIN

for kind in vet ref_ranges
do
    list_file="${kind}_list.txt"
    az storage blob list \
        --account-name "${STORAGE_ACCOUNT_NAME}" \
        --container-name "${STORAGE_CONTAINER_NAME}" \
        --prefix "${kind}" \
        --sas-token "$SAS_TOKEN" | \
    jq -r '.[] | .name' | \
    grep -E '.csv$' > "${list_file}"

    while IFS= read -r file
    do
        if [[ "${kind}" == "vet" ]]
        then
            cat >> "${LOADER_SCRIPT}" <<FIN
INSERT INTO vets with (TABLOCK)
  SELECT
      CONVERT(BIGINT, location) AS location,
      CONVERT(INT, sample_id) AS sample_id,
      ref,
      alt,
      AS_RAW_MQ,
      AS_RAW_MQRankSum,
      QUALapprox,
      AS_QUALapprox,
      AS_RAW_ReadPosRankSum,
      AS_SB_TABLE,
      AS_VarDP,
      call_GT,
      call_AD,
      call_GQ,
      call_PGT,
      call_PID,
      call_PL
  FROM OPENROWSET(
              BULK '${file}',
              DATA_SOURCE = 'MyAzureBlobStorage',
              FORMAT ='CSV',
              FORMATFILE='vets_format.txt',
              FORMATFILE_DATA_SOURCE = 'MyAzureBlobStorage'
  ) as DataFile;


FIN
        else
            cat >> "${LOADER_SCRIPT}" <<FIN
  INSERT INTO ref_ranges with (TABLOCK)
  SELECT
      CONVERT(BIGINT, location) AS location,
      CONVERT(INT, sample_id) AS sample_id,
      CONVERT(SMALLINT, length) AS length,
      CONVERT(TINYINT, state)
  FROM OPENROWSET(
              BULK '${file}',
              DATA_SOURCE = 'MyAzureBlobStorage',
              FORMAT ='CSV',
              FORMATFILE='ref_ranges_format.txt',
              FORMATFILE_DATA_SOURCE = 'MyAzureBlobStorage'
  ) as DataFile;


FIN
        fi
    done < "${list_file}"
done

chmod +x "${LOADER_SCRIPT}"

get_db_token() {
    # Fetches a database access token with ~1 hour TTL, ASCII / UTF-8 encoded and newline terminated
    # https://learn.microsoft.com/en-us/sql/connect/odbc/linux-mac/connecting-with-sqlcmd?view=azuresqldb-current
    az account get-access-token --resource https://database.windows.net --query accessToken --output tsv
}

get_db_token | tr -d '\n' | iconv -f ascii -t UTF-16LE > token.txt

sqlcmd -S "tcp:${SQL_SERVER}.database.windows.net,1433" -d "${SQL_DATABASE}" -G -P token.txt -Q "$(cat "${LOADER_SCRIPT}")"
