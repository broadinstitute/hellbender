#!/usr/bin/env bash

#NOTE: This script has been checked in to aid in the release process for future Funcotator datasource bundles.

echo "Making Tarballs of each Datasource Directory..."

tar -zcvf funcotator_dataSources.v1.8.hg38.20230908s.tar.gz funcotator_dataSources.v1.8.hg38.20230908s
tar -zcvf funcotator_dataSources.v1.8.hg38.20230908g.tar.gz funcotator_dataSources.v1.8.hg38.20230908g
tar -zcvf funcotator_dataSources.v1.8.hg19.20230908s.tar.gz funcotator_dataSources.v1.8.hg19.20230908s
tar -zcvf funcotator_dataSources.v1.8.hg19.20230908g.tar.gz funcotator_dataSources.v1.8.hg19.20230908g

echo "Making the various hashfiles for release"

find funcotator_dataSources.v1.8.hg38.20230908s -type f | xargs md5sum > funcotator_dataSources.v1.8.hg38.20230908s.dir.long.md5sum
md5sum funcotator_dataSources.v1.8.hg38.20230908s.tar.gz | awk '{print $1}' > funcotator_dataSources.v1.8.hg38.20230908s.dir.md5sum
sha256sum funcotator_dataSources.v1.8.hg38.20230908s.tar.gz > funcotator_dataSources.v1.8.hg38.20230908s.sha256

find funcotator_dataSources.v1.8.hg38.20230908g -type f | xargs md5sum > funcotator_dataSources.v1.8.hg38.20230908g.dir.long.md5sum
md5sum funcotator_dataSources.v1.8.hg38.20230908g.tar.gz | awk '{print $1}' > funcotator_dataSources.v1.8.hg38.20230908g.dir.md5sum
sha256sum funcotator_dataSources.v1.8.hg38.20230908g.tar.gz > funcotator_dataSources.v1.8.hg38.20230908g.sha256

find funcotator_dataSources.v1.8.hg19.20230908s -type f | xargs md5sum > funcotator_dataSources.v1.8.hg19.20230908s.dir.long.md5sum
md5sum funcotator_dataSources.v1.8.hg19.20230908s.tar.gz | awk '{print $1}' > funcotator_dataSources.v1.8.hg19.20230908s.dir.md5sum
sha256sum funcotator_dataSources.v1.8.hg19.20230908s.tar.gz > funcotator_dataSources.v1.8.hg19.20230908s.sha256

find funcotator_dataSources.v1.8.hg19.20230908g -type f | xargs md5sum > funcotator_dataSources.v1.8.hg19.20230908g.dir.long.md5sum
md5sum funcotator_dataSources.v1.8.hg19.20230908g.tar.gz | awk '{print $1}' > funcotator_dataSources.v1.8.hg19.20230908g.dir.md5sum
sha256sum funcotator_dataSources.v1.8.hg19.20230908g.tar.gz > funcotator_dataSources.v1.8.hg19.20230908g.sha256