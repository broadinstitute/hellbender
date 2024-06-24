#!/usr/bin/env bash

################################################################################

# Creates a set of files that allows for creating annotations based off dbSNP
# Pulled directly from the ncbi FTP site.

################################################################################

#Setup variables for the script:
SCRIPTDIR="$( cd -P "$( dirname "$0" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=0

# Latest release numbers for our references.
# Update these numbers when a new Gencode is released.
LATEST_RELEASE=43

DATA_SOURCE_NAME="Gencode"
OUT_DIR_NAME='gencode'

# Check if we have samtools installed so we can index our fasta files:
HAS_SAMTOOLS=false
which samtools &> /dev/null
r=$?
[[ $r -eq 0 ]] && HAS_SAMTOOLS=true

################################################################################

set -e

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME [OPTIONS] ..."
  echo -e "Creates the data sources folder for gencode for the GATK Funcotator tool."
}

#Define a usage function:
function usage()
{
  simpleUsage
  echo -e ""
  echo -e "Will download all data sources directly from the Gencode website:"
  echo -e "    https://www.gencodegenes.org"
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  UNKNOWN ARGUMENT"
  echo -e "  4  BAD CHECKSUM"
  echo -e "  5  OUTPUT DIRECTORY ALREADY EXISTS"
  echo -e ""
}

#Display a message to std error:
function error()
{
  echo "$1" 2>&1
}

TMPFILELIST=''
function makeTemp()
{
  local f
  f=$( mktemp )
  TMPFILELIST="${TMPFILELIST} $f"
  echo $f
}

function cleanTempVars()
{
  rm -f ${TMPFILELIST}
}

function at_exit()
{
  cleanTempVars
}

################################################################################


function createConfigFile() {

  local dataSourceName=$1    # Name of the data source to use in the config file (by which funcotator will refer to this data source).
  local version=$2           # Numerical version of the Gencode files to download (e.g. 34)
  local refVersion=$3        # Reference version for which to download gencode (e.g. hg19/hg38)
  local srcFile=$4           # Name of the main file from which the data for this datasource will be pulled.
  local originLocation=$5    # Path specifier / URI of where the srcFile was originally hosted / pulled from
  local fastaPath=$6         # Name of the FASTA file icontaining transcript sequences corresponding to this version of Gencode.

  echo "name = ${dataSourceName}"
  echo "version = ${version}"
  echo "src_file = ${srcFile}"
  echo "origin_location = ${originLocation}"
  echo "preprocessing_script = ${SCRIPTNAME} , fixGencodeOrdering.py"
  echo ""
  echo "# Supported types:"
  echo "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID"
  echo "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location"
  echo "# gencode      -- Custom datasource class for GENCODE"
  echo "# cosmic       -- Custom datasource class for COSMIC"
  echo "# vcf          -- Custom datasource class for Variant Call Format (VCF) files"
  echo "type = gencode"
  echo ""
  echo "# Required field for GENCODE files."
  echo "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:"
  echo "gencode_fasta_path = ${fastaPath}"
  echo ""
  echo "# Required field for GENCODE files."
  echo "# NCBI build version (either hg19 or hg38):"
  echo "ncbi_build_version = ${refVersion}"
  echo ""
  echo "# Required field for simpleXSV files."
  echo "# Valid values:"
  echo "#     GENE_NAME"
  echo "#     TRANSCRIPT_ID"
  echo "xsv_key ="
  echo ""
  echo "# Required field for simpleXSV files."
  echo "# The 0-based index of the column containing the key on which to match"
  echo "xsv_key_column ="
  echo ""
  echo "# Required field for simpleXSV AND locatableXSV files."
  echo "# The delimiter by which to split the XSV file into columns."
  echo "xsv_delimiter ="
  echo ""
  echo "# Required field for simpleXSV files."
  echo "# Whether to permissively match the number of columns in the header and data rows"
  echo "# Valid values:"
  echo "#     true"
  echo "#     false"
  echo "xsv_permissive_cols ="
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

function getGencodeFiles()
{
  local version=$1          # Numerical version of the Gencode files to download (e.g. 34)
  local refVersion=$2        # Reference version for which to download gencode (e.g. hg19/hg38)

  echo "##################################################################################"
  echo "Processing ${refVersion} - Gencode ${version} ..."
  echo "===================================="

  mkdir -p ${OUT_DIR_NAME}/${refVersion}
  pushd ${OUT_DIR_NAME}/${refVersion} &> /dev/null

  fileRefVersion=${version}
  if [[ "${refVersion}" == "hg38" ]] ; then
    sourceUrl=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/gencode.v${version}.annotation.gtf.gz
    wget ${sourceUrl} 
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/gencode.v${version}.pc_transcripts.fa.gz
  else
    fileRefVersion="${version}lift37"
    sourceUrl=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/GRCh37_mapping/gencode.v${fileRefVersion}.annotation.gtf.gz
    wget ${sourceUrl}
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/GRCh37_mapping/gencode.v${fileRefVersion}.pc_transcripts.fa.gz
  fi

  gunzip gencode.v${fileRefVersion}.annotation.gtf.gz
  gunzip gencode.v${fileRefVersion}.pc_transcripts.fa.gz

  # We must fix the information in the gencode gtf file:
  echo "Reordering Gencode GTF data ..."
  ${SCRIPTDIR}/fixGencodeOrdering.py gencode.v${fileRefVersion}.annotation.gtf | sed -e '1,5s$^#\([a-z]\)$##\1$' > gencode.v${fileRefVersion}.annotation.REORDERED.gtf

  # Clean up original file:
  rm gencode.v${fileRefVersion}.annotation.gtf

  echo "Creating config file ..."
  createConfigFile "${DATA_SOURCE_NAME}" "${version}" "${refVersion}" "gencode.v${fileRefVersion}.annotation.REORDERED.gtf" "${sourceUrl}" "gencode.v${fileRefVersion}.pc_transcripts.fa" > gencode.config

  if $HAS_SAMTOOLS ; then
    echo "Indexing Fasta File: gencode.v${fileRefVersion}.pc_transcripts.fa"
    samtools faidx gencode.v${fileRefVersion}.pc_transcripts.fa
  fi

  echo
  popd > /dev/null
}

################################################################################

trap at_exit EXIT

################################################################################

#Check given arguments:
if [[ $# -gt $MAXARGS ]] ; then
  usage
  exit 1
elif [[ $# -lt $MINARGS ]] ; then
  usage
  exit 2
fi

################################################################################

#Read args:
while [ $# -gt 0 ] ; do

  case "$1" in
    -h|--h|--help|-help|help)
      usage;
      exit 0;
      ;;
    *)
      error "${SCRIPTNAME}: invalid option: $1"
      simpleUsage
      echo "Try \`$SCRIPTNAME --help' for more information."
      exit 3;
      ;;
  esac

  #Get next argument in $1:
  shift
done

################################################################################
# Do the work here:

# Make sure we don't have anything in our out folder already:
if [[ -d ${OUT_DIR_NAME} ]] ; then
  error "Output directory already exists: ${OUT_DIR_NAME} - aborting!"
  exit 5
fi

# Get the link for HG19:
getGencodeFiles $LATEST_RELEASE hg19

# Get the link for HG38:
getGencodeFiles $LATEST_RELEASE hg38

if ! $HAS_SAMTOOLS ; then
  echo -e "\033[1;33;40m##################################################################################\033[0;0m"
  echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
  echo -e "\033[1;33;40m# \033[1;5;37;41mWARNING\033[0;0m: You \033[4;37;40mMUST\033[0;0m index both Gencode Fasta files before using this data source \033[1;33;40m#\033[0;0m"
  echo -e "\033[1;33;40m#             Use samtools faidx <FASTA_FILE>                                    #\033[0;0m"
  echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
  echo -e "\033[1;33;40m##################################################################################\033[0;0m"
fi


echo -e "\033[1;33;40m##################################################################################\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m# \033[1;5;37;41mWARNING\033[0;0m: You \033[4;37;40mMUST\033[0;0m create sequence dictionaries for both Gencode Fasta files #\033[0;0m"
echo -e "\033[1;33;40m#             before using this data source.                                     #\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m#             Use gatk CreateSequenceDictionary <FASTA_FILE>                     #\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m# \033[1;5;37;41mWARNING\033[0;0m: You \033[4;37;40mMUST\033[0;0m ALSO index the Gencode GTF files for both HG19 and HG38 #\033[0;0m"
echo -e "\033[1;33;40m#             before using this data source.                                     #\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m#             Use gatk IndexFeatureFile <GTF_FILE>                               #\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m##################################################################################\033[0;0m"

exit 0

