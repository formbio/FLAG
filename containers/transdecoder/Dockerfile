FROM ubuntu:18.04

ARG version=2.13.0

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get -y install --no-install-recommends \
    automake \
    build-essential \
    ca-certificates \
    bzip2 \
    cmake \
    curl \
    g++ \
    gcc \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    git \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev \
    python \
    python-dev \
    python-pip \
    python-setuptools \
    hmmer
  
#perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm URI::Escape module

#install blast
RUN cd /opt; wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${version}/ncbi-blast-${version}+-x64-linux.tar.gz; tar xvfz ncbi-blast-2.13.0+-x64-linux.tar.gz; cp ncbi-blast-2.13.0+/bin/* /usr/local/bin

ENV PATH "$PATH:/opt/ncbi-blast-2.13.0+/bin"

#install DBs
RUN mkdir /seq; cd /seq; mkdir RNASEQ; cd RNASEQ; mkdir DBs; cd DBs; mkdir TRINOTATE_RESOURCES; cd TRINOTATE_RESOURCES; mkdir TRINOTATE_V3; cd TRINOTATE_V3

ADD uniprot_sprot.pep /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/uniprot_sprot.pep

ADD Pfam-A.hmm /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/Pfam-A.hmm

RUN cd /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/; makeblastdb -in uniprot_sprot.pep -dbtype prot; hmmpress Pfam-A.hmm

#install transcoder
RUN cd /opt; wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.5.0.tar.gz; tar -xf TransDecoder-v5.5.0.tar.gz; rm TransDecoder-v5.5.0.tar.gz; cd TransDecoder-TransDecoder-v5.5.0; make test
ENV PATH "$PATH:/opt/TransDecoder-TransDecoder-v5.5.0"


# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts


ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
