#!/bin/bash

GATK_S3A_TEST_DIR=${GATK_S3A_TEST_DIR:-"/tmp/gatk-s3atest"}
CUR_DIR=$(realpath `dirname $0`)
HADOOP_VER=`grep "^final hadoopVersion" $CUR_DIR/../../build.gradle | sed -e "s/.*, '\([0-9.]\+\)')/\1/g"`

HCONF_DIR=${GATK_S3A_TEST_DIR}/etc/hadoop/
cat << EOF > ${HCONF_DIR}/core-site.xml
<configuration>
  <property>
    <name>fs.s3a.aws.credentials.provider</name>
    <value>org.apache.hadoop.fs.s3a.AnonymousAWSCredentialsProvider</value>
  </property>
</configuration>
EOF

# hadoop uses older apache-commons than gatk, so exclude them
HCORE_JAR=`find ${GATK_S3A_TEST_DIR}/share/hadoop/common/ -name '*.jar' | grep -v commons | paste -s -d :`
JAWS_JARS=`find ${GATK_S3A_TEST_DIR}/share/hadoop/tools/lib/ -name '*.jar' | grep -v commons | paste -s -d :`
JSR_JAR=`find ${GATK_S3A_TEST_DIR} -maxdepth 1 -name 'jsr203-s3a-*.jar'`
GATK_JAR=$(find ${CUR_DIR}/../../build/libs/ -name "gatk-package-*-local.jar")

CLASSPATH=${HCONF_DIR}:${HCORE_JAR}:${HAWS_JAR}:${JAWS_JARS}:${JSR_JAR}:${GATK_JAR}

# we do not use gatk command because we need to ensure the jsr203-s3a-*.jar is loaded before gatk-package-*.jar
CMD="java -cp ${CLASSPATH} \
-Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xms8g \
org.broadinstitute.hellbender.Main HaplotypeCaller \
-R s3a://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \
-I s3a://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam \
-O $CUR_DIR/out.vcf"

echo $CMD
$CMD
