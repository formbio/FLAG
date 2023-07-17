#Contains general highly used programs for bioinformatics
FROM ubuntu:20.04
LABEL maintainer="William Troy <William@formbio.com>"

#Included libraries are samtools, bcftools, htslib, picard, and bedtools
ENV SAMTOOLS_VER="1.17"
ENV HTSVER="1.17"
ENV PICARD_VER="3.0.0"
ENV BEDTOOLS_VER="2.31.0"


# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; apt-get upgrade; \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --allow-unauthenticated libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk r-base cpanminus pigz parallel

# Make destination and work directory
RUN mkdir -p /opt;
RUN mkdir -p /opt/bin;

# Install SamTools + BCFTools + htslib
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2; \
    tar -vxjf samtools-${SAMTOOLS_VER}.tar.bz2; rm samtools-${SAMTOOLS_VER}.tar.bz2; \
    cd samtools-${SAMTOOLS_VER}; \
    ./configure; make; make install; \
    mv samtools /opt/bin;
RUN cd /opt; wget https://github.com/samtools/htslib/releases/download/${HTSVER}/htslib-${HTSVER}.tar.bz2; \
    tar -vxjf htslib-${HTSVER}.tar.bz2; rm htslib-${HTSVER}.tar.bz2; \
    cd htslib-${HTSVER}; \
    ./configure; make; make install; \
    mv htslib /opt/bin;
RUN cd /opt; wget https://github.com/samtools/bcftools/releases/download/${HTSVER}/bcftools-${HTSVER}.tar.bz2; \
    tar -vxjf bcftools-${HTSVER}.tar.bz2; rm bcftools-${HTSVER}.tar.bz2; \
    cd bcftools-${HTSVER}; \
    ./configure; make; make install; \
    mv bcftools /opt/bin;

# Install Bedtools
RUN cd /opt; wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VER}/bedtools-${BEDTOOLS_VER}.tar.gz; \
    tar xvfz bedtools-${BEDTOOLS_VER}.tar.gz; rm bedtools-${BEDTOOLS_VER}.tar.gz; \
    cd bedtools2; \
    make; \
    cp -r bin/* /opt/bin;

RUN apt-get update; apt-get install -y build-essential checkinstall; apt-get upgrade; \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --allow-unauthenticated libc6-i386 libc6-x32

# Install Java
RUN cd /opt; wget https://download.oracle.com/java/17/latest/jdk-17_linux-x64_bin.deb; dpkg -i jdk-17_linux-x64_bin.deb

ENV PATH "/usr/local/bin:/usr/lib/jvm/jdk-17/bin/:$PATH"
ENV JAVA_HOME "/usr/lib/jvm/jdk-17"
# Install Picard
RUN cd /opt/bin; wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VER}/picard.jar;

RUN cp -R /opt/bin/* /usr/local/bin/; \
    chmod +x /usr/local/bin/picard.jar

ENV PICARD "/usr/local/bin"
ENV SNPEFF_HOME "/usr/local/bin/snpEff"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; DEBIAN_FRONTEND="noninteractive" apt-get install -y wget unzip git bzip2 curl gcc libbz2-dev libc6-dev libcurl4-gnutls-dev liblzma-dev make pigz zlib1g zlib1g-dev libncurses5-dev libncursesw5-dev  libvcflib-tools libvcflib-dev

RUN cd /usr/local; wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz; tar xvfz gffread-0.12.7.Linux_x86_64.tar.gz; mv gffread-0.12.7.Linux_x86_64/gffread bin

RUN curl https://rclone.org/install.sh | bash

RUN  cd /usr/local; wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2; tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
ENV FB_VER="1.3.6"
RUN cd /opt; git clone --recursive https://github.com/freebayes/freebayes.git; cp freebayes/scripts/* /usr/local/bin
RUN cd /opt; wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz; tar -xf seqkit_linux_amd64.tar.gz
ENV PATH "$PATH:/opt/"



# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
