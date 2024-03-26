FROM ubuntu:20.04

RUN apt-get -qq update && DEBIAN_FRONTEND="noninteractive" apt-get -qq -y install --no-install-recommends \
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
    default-jre \
    python-setuptools \
    libdbi-perl \
    parallel \
    gffread 
  
#perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm URI::Escape module

#install EVM
RUN cd /opt; git clone https://github.com/EVidenceModeler/EVidenceModeler.git

#install conda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && chmod -R o+rX /opt/conda
# change conda for singularity
ENV PATH="/opt/conda/bin:${PATH}"
RUN conda init
RUN conda update -n base -c defaults conda
ADD env.yml /opt/env.yml
RUN conda env create -f /opt/env.yml

#install agat
RUN cpanm install Clone 
RUN cpanm install Graph::Directed 
RUN cpanm install LWP::UserAgent 
RUN cpanm install Carp 
RUN cpanm install Sort::Naturally File::Share File::ShareDir::Install Moose YAML LWP::Protocol::https
RUN cpanm install Term::ProgressBar
RUN cpanm install IPC::Run --force
RUN cpanm install Bio::DB::EUtilities
RUN cd /opt; git clone https://github.com/NBISweden/AGAT.git; cd AGAT; perl Makefile.PL; make; make install 

#begin gFACs install
RUN cpanm install Bio::Index::Fasta

RUN cd /opt; wget https://gitlab.com/PlantGenomicsLab/gFACs/-/archive/master/gFACs-master.tar.gz; tar -xf gFACs-master.tar.gz; cd gFACs-master

ENV PATH "$PATH:/opt/gFACs-master/:/opt/gFACs-master/format_scripts/:/opt/gFACs-master/support_scripts/:/opt/gFACs-master/task_scripts/"

RUN apt-get update; apt-get -y install samtools genometools 
# install stuff for bioconda genbank conversion
ADD biopython_environment.yml /opt/biopython_environment.yml
RUN conda env create -f /opt/biopython_environment.yml
RUN pip3 install pandas 
# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts


ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
