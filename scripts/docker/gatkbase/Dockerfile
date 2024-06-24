# Using OpenJDK 17
# This Dockerfile does not require any files that are in the GATK4 repo.
# NOTE: If you update the ubuntu version make sure to update the samtools/bcftools/bedtools versions in the README
FROM ubuntu:22.04

# Set environment variables.
# Avoid interactive prompts during apt installs/upgrades
ENV DEBIAN_FRONTEND="noninteractive" HOME="/root" JAVA_LIBRARY_PATH="/usr/lib/jni" DOWNLOAD_DIR="/downloads" CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-py310_23.10.0-1-Linux-x86_64.sh" CONDA_SHA256="c7a34df472feb69805b64df6e8db58363c5ccab41cd3b40b07e3e6dfb924359a" CONDA_PATH="/opt/miniconda" PATH="/opt/miniconda/bin:$PATH"

# Define working directory.
WORKDIR /root

#### Basic image utilities, google cloud support, and miniconda
RUN apt update && \
    apt full-upgrade -y && \
    apt install -y --no-install-recommends \
    python3 \
    python3-pip \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    samtools \
    bcftools \
    tabix \
    jq \
    git \
    gpg-agent \
    build-essential \
    openjdk-17-jdk \
    vim \
    software-properties-common && \
    apt -y clean  && \
    apt -y autoclean  && \
    apt -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    java -version && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt update -y && \
    apt install -y --no-install-recommends google-cloud-cli && \
    apt -y clean  && \
    apt -y autoclean  && \
    apt -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
#  Remove the anthos cli tool and related files since they are very large and we / anyone using the docker are unlikely to use them
#  Remove the bundled python because we have python installed separately
    rm -rf /usr/lib/google-cloud-sdk/bin/anthoscli /usr/lib/google-cloud-sdk/platform/anthoscli_licenses /usr/lib/google-cloud-sdk/platform/bundledpythonunix && \
    find / -wholename "*__pycache__/*.pyc" -exec rm {} + && \
    mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "$(sha256sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS=' ' -v ORS='' '{print $1}')" = "$CONDA_SHA256" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh && \
    conda clean -afy && \
    conda config --set auto_update_conda false && \
    conda config --set solver libmamba && \
    rm -rf /root/.cache/pip

CMD ["bash"]