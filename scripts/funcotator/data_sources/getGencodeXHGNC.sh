#!/usr/bin/env bash

# Creates a set of files that map records between GENCODE and HGNC.
# Pulled directly from the ensemble database.
# Currently not used by Funcotator.

outFileBaseName="gencode_xhgnc"
outExt=".tsv"

hg19db="homo_sapiens_core_75_37"
hg38db="homo_sapiens_core_110_38"

hg19FileName=${outFileBaseName}_v75_37.hg19${outExt}
hg38FileName=${outFileBaseName}_v110_38.hg38${outExt}

################################################################################

query="SELECT CONCAT(transcript.stable_id, '.', transcript.version) as transcript_id, xref.display_label as hgnc_id FROM translation, transcript, object_xref, xref, external_db WHERE transcript.transcript_id = translation.transcript_id AND translation.translation_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Translation' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = 'EMBL';"

echo "Getting HG19 gencode <=> hgnc..."
echo -e "transcript_id\thgnc_id" > ${hg19FileName} 
time mysql -u anonymous -h ensembldb.ensembl.org -e "use ${hg19db};${query}" | tail -n +2 | sort -n -k1 >> ${hg19FileName}

echo "Getting HG38 gencode <=> hgnc..."
echo -e "transcript_id\thgnc_id" > ${hg38FileName}
time mysql -u anonymous -h ensembldb.ensembl.org -e "use ${hg38db};${query}" | tail -n +2 | sort -n -k1 >> ${hg38FileName}

echo 'Done!'


