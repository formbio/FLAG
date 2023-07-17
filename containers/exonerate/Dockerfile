FROM ubuntu:20.04

RUN apt update && DEBIAN_FRONTEND="noninteractive" apt install -y --no-install-recommends \
  autoconf \
  automake \
  build-essential \
  gcc \
  libglib2.0-dev \
  make \
  git \
  ca-certificates \
  parallel \
  wget \
  && rm -rf /var/lib/apt/lists/*
  
WORKDIR /opt
RUN git clone https://github.com/nathanweeks/exonerate.git
RUN cd exonerate \
  && autoreconf -i \
  && ./configure \
  && make -j \
  && make test \
  && make install

RUN cd /opt; wget https://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz; tar -xf gth-1.7.3-Linux_x86_64-64bit.tar.gz

ENV PATH "$PATH:/opt/gth-1.7.3-Linux_x86_64-64bit/bin/"

RUN cd /opt; git clone https://github.com/lh3/miniprot; cd miniprot; make
ENV PATH "$PATH:/opt/miniprot/"

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
