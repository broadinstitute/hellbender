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
  <property>
    <name>fs.s3a.connection.maximum</name>
    <value>100</value>
  </property>
</configuration>
EOF

# hadoop uses older apache-commons, etc. than gatk and spark, so exclude them
HCORE_JAR=`find ${GATK_S3A_TEST_DIR}/share/hadoop/common/ -name '*.jar' | grep -v commons | grep -v jackson | grep -v avro | paste -s -d :`
JAWS_JARS=`find ${GATK_S3A_TEST_DIR}/share/hadoop/tools/lib/ -name '*.jar' | grep -v commons | grep -v jackson | grep -v avro | paste -s -d :`
JSR_JAR=`find ${GATK_S3A_TEST_DIR} -maxdepth 1 -name 'jsr203-s3a-*.jar'`
GATK_JAR=$(find ${CUR_DIR}/../../build/libs/ -name "gatk-package-*-local.jar")

CLASSPATH=${HCONF_DIR}:${HCORE_JAR}:${HAWS_JAR}:${JAWS_JARS}:${JSR_JAR}:${GATK_JAR}

# we do not use gatk command because we need to ensure the jsr203-s3a-*.jar is loaded before gatk-package-*.jar
CMD="java -cp ${CLASSPATH} \
org.broadinstitute.hellbender.Main CountReadsSpark \
--spark-master local[*] \
-I s3a://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam"

echo $CMD
$CMD
