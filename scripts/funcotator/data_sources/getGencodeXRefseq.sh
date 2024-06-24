#!/usr/bin/env bash

# Creates a set of files that map records between GENCODE and RefSeq.
# Pulled directly from the ensemble database.
# Currently not used by Funcotator.

outFileBaseName="gencode_xrefseq"
outExt=".tsv"

hg19db="homo_sapiens_core_75_37"
hg38db="homo_sapiens_core_110_38"

hg19FileName=${outFileBaseName}_v75_37.hg19${outExt}
hg38FileName=${outFileBaseName}_v110_38.hg38${outExt}

################################################################################


# Create our query to the DB:
read -r -d '' query <<- ENDOFQUERYINPUT 
SELECT mrna.transcript_id as transcript_id, mRNA_id, prot_acc FROM 
	(
		SELECT CONCAT(transcript.stable_id, '.', transcript.version) AS transcript_id, xref.display_label AS mRNA_id  
			FROM transcript, object_xref, xref, external_db
			WHERE 
					transcript.transcript_id = object_xref.ensembl_id AND 
					object_xref.ensembl_object_type = 'Transcript'    AND 
					object_xref.xref_id = xref.xref_id                AND 
					xref.external_db_id = external_db.external_db_id  AND 
					external_db.db_name = 'RefSeq_mRNA'
	) AS mrna
	JOIN
	(
		SELECT CONCAT(transcript.stable_id, '.', transcript.version) AS transcript_id, xref.display_label AS prot_acc 
			FROM translation, transcript, object_xref, xref,external_db
			WHERE
				(
					transcript.transcript_id = translation.transcript_id  AND 
					translation.translation_id = object_xref.ensembl_id  AND 
					object_xref.ensembl_object_type = 'Translation'      AND 
					object_xref.xref_id = xref.xref_id                   AND 
					xref.external_db_id = external_db.external_db_id     AND 
					external_db.db_name = 'RefSeq_peptide'
				)
	) AS prot
	ON mrna.transcript_id = prot.transcript_id 
;
ENDOFQUERYINPUT

echo "Getting HG19 gencode <=> refseq..."
echo -e "transcript_id\tmRNA_id\tprot_acc" > ${hg19FileName} 
time mysql -u anonymous -h ensembldb.ensembl.org -e "use ${hg19db};${query}" | tail -n +2 | sort -n -k1 >> ${hg19FileName}

echo "Getting HG38 gencode <=> refseq..."
echo -e "transcript_id\tmRNA_id\tprot_acc" > ${hg38FileName} 
time mysql -u anonymous -h ensembldb.ensembl.org -e "use ${hg38db};${query}" | tail -n +2 | sort -n -k1 >> ${hg38FileName}

echo 'Done!'


