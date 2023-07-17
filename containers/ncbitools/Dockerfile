FROM ubuntu:20.04

ARG version=2.13.0

RUN apt-get update; apt-get -y -m update && DEBIAN_FRONTEND="noninteractive" apt-get install -y build-essential curl libidn11 libnet-perl perl-doc liblmdb-dev wget unzip git jq libdw-dev cpio pigz emboss
RUN apt-get update; apt-get -y -m update && DEBIAN_FRONTEND="noninteractive" apt-get install -y libgomp1 libnet-perl libxml-simple-perl libjson-perl perl-doc liblmdb-dev parallel vmtouch cpanminus curl python3-dev python3-pip fasta3 g++ automake cmake zlib1g-dev subversion libzstd-dev && rm -rf /var/lib/apt/lists/*  && cpanm HTML::Entities

RUN python3 -m pip install requests

RUN cd /opt; wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${version}/ncbi-blast-${version}+-x64-linux.tar.gz; tar xvfz ncbi-blast-${version}+-x64-linux.tar.gz; cp ncbi-blast-${version}+/bin/* /usr/local/bin
RUN cd /opt; sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"; mv /root/edirect/* /usr/local/bin; cd /usr/local/bin

RUN cd /opt; wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz; tar -xzf RpsbProc-x64-linux.tar.gz; cp RpsbProc-x64-linux/rpsbproc /usr/local/bin

RUN cd /opt; wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz; tar xvfz sratoolkit.2.11.3-ubuntu64.tar.gz; cp -R sratoolkit.2.11.3-ubuntu64/bin/* /usr/local/bin

#install NCBI Datasets CLI tool (v13) 
RUN cd /opt; curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v1/linux-amd64/datasets'; 
RUN cd /opt; curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v1/linux-amd64/dataformat'; 
RUN cd /opt; chmod +x datasets dataformat; cp datasets /usr/local/bin; cp dataformat /usr/local/bin; 

#install diamond
RUN cd /opt; wget http://github.com/bbuchfink/diamond/releases/download/v2.0.14/diamond-linux64.tar.gz; tar xzf diamond-linux64.tar.gz
ENV PATH "$PATH:/opt/"

RUN python3 -m pip install ffq
RUN python3 -m pip install snapgene_reader

RUN cd opt; wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz; tar xvfz seqkit_linux_amd64.tar.gz; cp seqkit /usr/local/bin/

#RUN mkdir -p /root/.ncbi
#COPY sratools.config /root/.ncbi/user-settings.mkfg

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/

