FROM gcr.io/bioinfo-devel/basebioinfo:latest

RUN apt-get -qq update && DEBIAN_FRONTEND="noninteractive" apt-get -qq -y install --no-install-recommends \
    automake \
    build-essential \
    ca-certificates \
    bzip2 \
    cmake \
    curl \
    g++ \
    gcc \
    git \
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
    perl \
    pkg-config \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev \
    default-jre \
    csh \
    ant \
    python3 \
    python3-pip \
    python3-setuptools \
    emboss \
    gffread
  
RUN pip3 install setuptools
RUN pip3 install networkx numpy biopython==1.76

RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install DB_File
RUN cpanm install URI::Escape

#install agat
RUN cpanm URI::Escape module
RUN cpanm install Clone 
RUN cpanm install Graph::Directed 
RUN cpanm install LWP::UserAgent 
RUN cpanm install Carp 
RUN cpanm install Sort::Naturally File::Share File::ShareDir::Install Moose YAML LWP::Protocol::https
RUN cpanm install Term::ProgressBar
RUN cpanm install IPC::Run --force
RUN cpanm install Bio::DB::EUtilities
RUN cd /opt; git clone https://github.com/NBISweden/AGAT.git; cd AGAT; perl Makefile.PL; make; make install 

RUN cd /opt; git clone https://github.com/harta55/EnTAP.git; cd EnTAP; cd libs/diamond-0.8.31; mkdir bin; cd bin; cmake ..; make install; make

#RUN cd /opt/EnTAP/libs/RSEM-1.3.3; make; make ebseq; make install
RUN cd /opt/EnTAP/libs/diamond-0.9.9; mkdir bin; cd bin; cmake ..; make install; make

RUN cd /opt/EnTAP; cmake CMakeLists.txt; make; make install

RUN cd /opt/EnTAP; rm entap_config.ini
ADD entap_config.ini /opt/EnTAP/


# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/

