#!/usr/bin/env bash

################################################################################

set -e

################################################################################

SCRIPTDIR="$( cd -P "$( dirname "$0" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )

################################################################################

version="v98"

EMAIL=""
PASSWORD=""

################################################################################

function usage() {
	echo "This script creates the cosmic data sources for the Funcotator GATK tool."
	echo "It will download the raw data files from the COSMIC SFTP repository and "
	echo "create 3 folders in the current working directory, one for each data"
	echo "source:"
	echo "    cosmic"
	echo "    cosmic_tissue"
	echo "    cosmic_fusion"
	echo "Note: This script requires an account with COSMIC and will prompt the "
	echo "      user to input these credentials in order to download the raw data."
}

################################################################################

function createConfigFile() {

    dataSourceName=$1
    version=$2
    srcFile=$3
    originLocation=$4
    preprocessors=$5
    type=$6
    key=$7
    col=$8
    delim=$9

    echo "name = ${dataSourceName}"
    echo "version = ${version}"
    echo "src_file = ${srcFile}"
    echo "origin_location = ${originLocation}"
    echo "preprocessing_script = ${preprocessors}"
    echo ""
    echo "# Supported types:"
    echo "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID"
    echo "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location"
    echo "# gencode      -- Custom datasource class for GENCODE"
    echo "# cosmic       -- Custom datasource class for COSMIC"
    echo "type = ${type}"
    echo ""
    echo "# Required field for GENCODE files."
    echo "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:"
    echo "gencode_fasta_path ="
    echo ""
    echo "# Required field for simpleXSV files."
    echo "# Valid values:"
    echo "#     GENE_NAME"
    echo "#     TRANSCRIPT_ID"
    echo "xsv_key = ${key}"
    echo ""
    echo "# Required field for simpleXSV files."
    echo "# The 0-based index of the column containing the key on which to match"
    echo "xsv_key_column = ${col}"
    echo ""
    echo "# Required field for simpleXSV AND locatableXSV files."
    echo "# The delimiter by which to split the XSV file into columns."
    echo "xsv_delimiter = ${delim}"
    echo ""
    echo "# Required field for simpleXSV files."
    echo "# Whether to permissively match the number of columns in the header and data rows"
    echo "# Valid values:"
    echo "#     true"
    echo "#     false"
    echo "xsv_permissive_cols = ${PERMISSIVE_COLS}"
    echo ""
    echo "# Required field for locatableXSV files."
    echo "# The 0-based index of the column containing the contig for each row"
    echo "contig_column ="
    echo ""
    echo "# Required field for locatableXSV files."
    echo "# The 0-based index of the column containing the start position for each row"
    echo "start_column ="
    echo ""
    echo "# Required field for locatableXSV files."
    echo "# The 0-based index of the column containing the end position for each row"
    echo "end_column ="
    echo ""

}

################################################################################

if [[ $# -gt 0 ]] ; then
	usage
	exit 0
fi

################################################################################

echo "This script creates the cosmic data sources for the Funcotator GATK tool."
echo ""
echo "For usage information run with the '-h' option"
echo ""
echo "To retrieve the COSMIC data sources you must have a COSMIC account."
echo "Please enter your COSMIC account credentials:"

while [[ ${#EMAIL} -eq 0 ]] ; do
    echo -n "Enter your email address: "
    read EMAIL
done

while [[ ${#PASSWORD} -eq 0 ]] ; do
    echo -n "Enter your password: "
    read -s PASSWORD
done
echo

################################################################################
# Setup the Directory:

echo "Creating folders: ..."

mkdir -vp cosmic/hg19 cosmic/hg38 cosmic_fusion/hg19 cosmic_fusion/hg38 cosmic_tissue/hg19 cosmic_tissue/hg38

################################################################################
# Get the data files:

echo "Getting files ... "
AUTH_TOKEN=$(echo -n "$EMAIL:$PASSWORD" | base64)

curl -sS "$(curl -H "Authorization: Basic $AUTH_TOKEN"  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz | jq -r '.url')" -o cosmic/hg19/CosmicCompleteTargetedScreensMutantExport.tsv.gz
curl -sS "$(curl -H "Authorization: Basic $AUTH_TOKEN"  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/${version}/CosmicFusionExport.tsv.gz | jq -r '.url')" -o cosmic_fusion/hg19/CosmicFusionExport.tsv.gz
curl -sS "$(curl -H "Authorization: Basic $AUTH_TOKEN"  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz | jq -r '.url')" -o cosmic/hg38/CosmicCompleteTargetedScreensMutantExport.tsv.gz
curl -sS "$(curl -H "Authorization: Basic $AUTH_TOKEN"  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/${version}/CosmicFusionExport.tsv.gz | jq -r '.url')" -o cosmic_fusion/hg38/CosmicFusionExport.tsv.gz

echo "Retrieved COSMIC version ${version} on $(date) from sftp-cancer.sanger.ac.uk by: ${SCRIPTNAME}:" >  cosmic/metadata.txt
echo "User: ${EMAIL}" >> cosmic/metadata.txt
echo "    cosmic/grch37/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz" >> cosmic/metadata.txt
echo "    cosmic/grch38/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz" >> cosmic/metadata.txt
echo "" >> cosmic/metadata.txt

echo "Retrieved COSMIC version ${version} on $(date) from sftp-cancer.sanger.ac.uk by: ${SCRIPTNAME}:" >  cosmic_fusion/metadata.txt
echo "User: ${EMAIL}" >> cosmic_fusion/metadata.txt
echo "    cosmic/grch37/cosmic/${version}/CosmicFusionExport.tsv.gz" >> cosmic_fusion/metadata.txt
echo "    cosmic/grch38/cosmic/${version}/CosmicFusionExport.tsv.gz" >> cosmic_fusion/metadata.txt
echo "" >> cosmic_fusion/metadata.txt

echo "Retrieved COSMIC version ${version} on $(date) from sftp-cancer.sanger.ac.uk by: ${SCRIPTNAME}:" >  cosmic_tissue/metadata.txt
echo "User: ${EMAIL}" >> cosmic_tissue/metadata.txt
echo "    cosmic/grch37/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz" >> cosmic_tissue/metadata.txt
echo "    cosmic/grch38/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz" >> cosmic_tissue/metadata.txt
echo "" >> cosmic_tissue/metadata.txt

################################################################################
# Unzip the files, process into usable files for funcotator, and create
# config files for these data sources:

for d in "hg19" "hg38" ; do

    echo "Unzipping data for ${d} ..."

    gunzip cosmic/${d}/CosmicCompleteTargetedScreensMutantExport.tsv.gz
    gunzip cosmic_fusion/${d}/CosmicFusionExport.tsv.gz

    echo "Processing data source for Cosmic ${d} ..."
    ${SCRIPTDIR}/createSqliteCosmicDb.sh cosmic/${d}/CosmicCompleteTargetedScreensMutantExport.tsv cosmic/${d}/Cosmic.db

    echo "Processing data source for Cosmic Fusion ${d} ..."
    ${SCRIPTDIR}/createCosmicFusionGeneTsv.py cosmic_fusion/${d}/CosmicFusionExport.tsv cosmic_fusion/${d}/cosmic_fusion.tsv

    echo "Processing data source for Cosmic Tissue ${d} ..."
    ${SCRIPTDIR}/createCosmicGeneTsv.py cosmic/${d}/CosmicCompleteTargetedScreensMutantExport.tsv cosmic_tissue/${d}/cosmic_tissue.tsv

    ref="grch37"
    if [[ "${d}" == "hg38" ]] ; then
        ref="grch38"
    fi

    echo "Creating config file for Cosmic Fusion ${d} ..."
    PERMISSIVE_COLS="true"
    createConfigFile \
        "CosmicFusion" \
        "${version}" \
        "cosmic_fusion.tsv" \
        "sftp://sftp-cancer.sanger.ac.uk/cosmic/${ref}/cosmic/${version}/CosmicFusionExport.tsv.gz" \
        "${SCRIPTNAME}, createCosmicFusionGeneTsv.py" \
        "simpleXSV" \
        "GENE_NAME" \
        "0" \
        '\t' > cosmic_fusion/${d}/cosmic_fusion.config

    echo "Creating config file for Cosmic Tissue ${d} ..."
    PERMISSIVE_COLS="true"
    createConfigFile \
        "CosmicTissue" \
        "${version}" \
        "cosmic_tissue.tsv" \
        "sftp://sftp-cancer.sanger.ac.uk/cosmic/${ref}/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz" \
        "${SCRIPTNAME}, createCosmicGeneTsv.py" \
        "simpleXSV" \
        "GENE_NAME" \
        "0" \
        '\t' > cosmic_tissue/${d}/cosmic_tissue.config

    echo "Creating config file for Cosmic ${d} ..."
    PERMISSIVE_COLS=""
    createConfigFile \
        "Cosmic" \
        "${version}" \
        "Cosmic.db" \
        "sftp://sftp-cancer.sanger.ac.uk/cosmic/${ref}/cosmic/${version}/CosmicCompleteTargetedScreensMutantExport.tsv.gz" \
        "${SCRIPTNAME}, createSqliteCosmicDb.sh" \
        "cosmic" \
        "" \
        "" \
        ''  > cosmic/${d}/cosmic.config

    echo "Cleaning up raw files ..."
    rm -v cosmic/${d}/CosmicCompleteTargetedScreensMutantExport.tsv
    rm -v cosmic_fusion/${d}/CosmicFusionExport.tsv
done

echo "DONE."
