#!/usr/bin/env bash

# Usage:
# Creates a sqlite3 database file based on a COSMIC download file.
# This resulting database file will only contain records that have either
# genome positions or protein positions so they can be matched to variants
# by Funcotator.
# To change the input file, change the `COSMIC_FILE` variable or pass a path
# in on the command-line.
# To change the output file, change the `OUT_DB_FILE` variable or pass a path
# in on the command-line for the COSMIC_FILE and a SECOND path for the OUT_DB_FILE

################################################################################

set -e

################################################################################

COSMIC_FILE=CosmicCompleteTargetedScreensMutantExport.tsv
OUT_DB_FILE="Cosmic.db"
OUT_TMP_FOLDER="~/tmp"

################################################################################

if [[ $# -gt 0 ]] ; then
    COSMIC_FILE=$1
fi

if [[ $# -gt 1 ]] ; then
    OUT_DB_FILE=$2
fi

if [[ $# -gt 2 ]] ; then
    OUT_TMP_FOLDER=$3
fi

if [ ! -f ${COSMIC_FILE} ] ; then
    echo "ERROR: Given COSMIC file does not exist: ${COSMIC_FILE}" 1>&2
    exit 1
fi

echo "Creating COSMIC database from file: ${COSMIC_FILE} -> ${OUT_DB_FILE}"

[ -f ${OUT_DB_FILE} ] && rm ${OUT_DB_FILE}

sqlite3 ${OUT_DB_FILE} <<EOF
.echo on
.mode tabs
.import ${COSMIC_FILE} RawCosmic
pragma temp_store_directory = ${OUT_TMP_FOLDER};
CREATE TABLE Cosmic AS SELECT * FROM RawCosmic WHERE ("Mutation AA" != "" OR "Mutation genome position" != "");
DROP TABLE RawCosmic;
UPDATE Cosmic SET "Mutation genome position" = "chr"||"Mutation genome position" WHERE "Mutation genome position" != "";
CREATE INDEX GeneIndex ON Cosmic("Gene name");
VACUUM;
EOF

echo 'DONE!'
