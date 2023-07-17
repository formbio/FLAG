#Contains liftoff
FROM ubuntu:18.04
LABEL maintainer="William Troy <William@formbio.com>"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; apt-get upgrade; \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --allow-unauthenticated libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk r-base cpanminus pigz parallel python3.7 python3-distutils curl

#install minimap2
RUN cd /opt; git clone https://github.com/lh3/minimap2; cd minimap2; make;
ENV PATH "$PATH:/opt/minimap2"

#install pip
RUN cd /opt; curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py; python3.7 get-pip.py; pip3.7 install six==1.12.0

#RUN cd /opt; git clone https://github.com/agshumate/Liftoff liftoff; cd liftoff; python3.7 setup.py install

RUN pip3.7 install Liftoff==1.6.3.2

#perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm URI::Escape module

#install agat
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
RUN cd /opt; git clone https://github.com/NBISweden/AGAT.git; cd AGAT; perl Makefile.PL; make; make install 

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
