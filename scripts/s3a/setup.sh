#!/bin/bash

CUR_DIR=$(dirname `realpath $0`)
GATK_S3A_TEST_DIR=${GATK_S3A_TEST_DIR:-"/tmp/gatk-s3atest"}
HADOOP_VER=`grep "^final hadoopVersion" $CUR_DIR/../../build.gradle | sed -e "s/.*, '\([0-9.]\+\)')/\1/g"`

mkdir -p $GATK_S3A_TEST_DIR
echo "Downloading hadoop-${HADOOP_VER}..."
curl -L https://archive.apache.org/dist/hadoop/core/hadoop-${HADOOP_VER}/hadoop-${HADOOP_VER}.tar.gz  | tar -xz -C $GATK_S3A_TEST_DIR --strip-components 1
echo "Downloading jsr203-s3a-0.0.2.jar..."
curl -L https://repo1.maven.org/maven2/net/fnothaft/jsr203-s3a/0.0.2/jsr203-s3a-0.0.2.jar -o ${GATK_S3A_TEST_DIR}/jsr203-s3a-0.0.2.jar