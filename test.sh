#!/bin/bash

#// TODO checkout table output format MAF

ref=/home/sshenker/GRCh38.fa
data="/u01/compbio/data/sshenker/funcotator/data-source/"
#gatk=/home/sshenker/git/gatk/build/libs/gatk-4.1.6.0-3-g43c4bc2-SNAPSHOT.jar
#java -Xmx8g -jar $gatk org.broadinstitute.hellbender.tools.funcotator.Funcotator \

./gatk Funcotator \
	   -R $ref \
	   -V broken.vcf \
	   -O out.vcf.gz \
	   --data-sources-path $data \
	   --output-file-format VCF \
	   --ref-version hg38
