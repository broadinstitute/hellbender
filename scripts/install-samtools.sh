#!/bin/sh
set -ex
#ubuntu specific
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install -y libncurses-dev libbz2-dev liblzma-dev

#install from the github tar
export SAMTOOLS_VERSION=1.21
wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar -xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION} && ./configure --prefix=/usr && make && sudo make install
