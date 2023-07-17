FROM pasapipeline/pasapipeline

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN apt-get update; apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libmysql++-dev libboost-iostreams-dev zlib1g-dev libbamtools-dev samtools libhts-dev cufflinks gffread genometools cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl

#add pasa pipeline perl script to path
RUN cd /opt; wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

ENV PATH="/root/miniconda3/bin:${PATH}"

Run conda create -n BUSCO -c conda-forge -c bioconda busco=5.3.2
RUN conda init
RUN apt-get install -y --no-install-recommends python3-biopython
ENV PATH "$PATH:/usr/local/src/PASApipeline"

# Set timezone in tzdata
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Berlin"

# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential rsync perl libipc-run-perl wget git autoconf cpanminus parallel
RUN apt-get install -y expat libexpat1-dev libcurl4-openssl-dev
RUN cpanm URI::Escape module
RUN cpanm install Clone
RUN cpanm install Graph::Directed
RUN cpanm install LWP::UserAgent
RUN cpanm install Carp
RUN cpanm install Sort::Naturally File::Share File::ShareDir::Install Moose YAML LWP::Protocol::https
RUN cpanm install Term::ProgressBar
RUN cpanm install IPC::Run --force
RUN cpanm install XML::LibXML --force
RUN cpanm install Bio::Root::IO
RUN cpanm install Bio::SeqIO
RUN cpanm install Bio::DB::EUtilities
#install blat
RUN cd /opt; mkdir blat; cd blat; rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

ENV PATH "$PATH:/opt/blat/blat:/opt/blat"


#install pblat
RUN cd /opt; git clone https://github.com/icebert/pblat.git; cd pblat; make;
ENV PATH "$PATH:/opt/pblat"


# Clone AUGUSTUS repository
Run git clone https://github.com/Gaius-Augustus/Augustus.git /root/augustus

# Build AUGUSTUS
WORKDIR "/root/augustus"
RUN make clean
RUN make
RUN make install
ENV PATH="/root/augustus/bin:/root/augustus/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH="/root/augustus/config/"
# Test AUGUSTUS
RUN make unit_test


RUN apt-get -y install cufflinks gffread samtools genometools

#clone in GALBA for the protein alignments 2 hints perl script
RUN cd /opt; git clone https://github.com/Gaius-Augustus/GALBA.git

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

# Install scipio
RUN mkdir -p /scipio;
ADD scipio /scipio

#Install modified autoAugTrain.pl script
RUN rm /root/augustus/scripts/autoAugTrain.pl
ADD autoAugTrain.pl /root/augustus/scripts/autoAugTrain.pl

#Install modified autoAug.pl script as autoaugpasa.pl
ADD autoAugpasa.pl /root/augustus/scripts/autoAugpasa.pl
RUN chmod 777 /root/augustus/scripts/autoAugpasa.pl

#Install modified E config file
RUN rm /root/augustus/config/extrinsic/extrinsic.E.cfg
ADD extrinsic.E.cfg /root/augustus/config/extrinsic/extrinsic.E.cfg
RUN chmod 777 /root/augustus/config/extrinsic/extrinsic.E.cfg

#change E config file to E.XNT config file
#RUN rm /root/augustus/config/extrinsic/extrinsic.E.cfg; cp /root/augustus/config/extrinsic/extrinsic.E.XNT.cfg /root/augustus/config/extrinsic/extrinsic.E.cfg
#RUN chmod 777 /root/augustus/config/extrinsic/extrinsic.E.cfg
RUN cd /scipio; chmod 777 *
ENV PATH "$PATH:/scipio"

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
