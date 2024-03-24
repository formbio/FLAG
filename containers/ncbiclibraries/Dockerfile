FROM ubuntu:18.04

ARG version=2.13.0

RUN apt-get -y -m update; apt-get -y upgrade; apt-get install -y build-essential curl libidn11 libnet-perl perl-doc liblmdb-dev wget unzip git jq libdw-dev cpio
RUN apt-get -y -m update; apt-get -y upgrade; apt-get install -y libgomp1 autoconf libnet-perl libxml-simple-perl libncurses5-dev libncursesw5-dev libbz2-dev libjson-perl perl-doc liblmdb-dev parallel vmtouch cpanminus curl zlib1g make lmdb-utils mysql-server protobuf-compiler ccache wine64 libuv-dev nghttp2 libblas-dev liblapack-dev genometools && rm -rf /var/lib/apt/lists/*  && cpanm HTML::Entities

#RUN python3 -m pip install requests

RUN cd /opt; wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${version}/ncbi-blast-${version}+-x64-linux.tar.gz; tar xvfz ncbi-blast-${version}+-x64-linux.tar.gz; cp ncbi-blast-${version}+/bin/* /usr/local/bin
RUN cd /opt; sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"; mv /root/edirect/* /usr/local/bin; cd /usr/local/bin

RUN cd /opt; wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz; tar -xzf RpsbProc-x64-linux.tar.gz; cp RpsbProc-x64-linux/rpsbproc /usr/local/bin

RUN cd /opt; wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz; tar xvfz sratoolkit.2.11.3-ubuntu64.tar.gz; cp -R sratoolkit.2.11.3-ubuntu64/bin/* /usr/local/bin

#install odbc
RUN cd /opt; wget ftp://ftp.unixodbc.org/pub/unixODBC/unixODBC-2.3.9.tar.gz; tar xvzf unixODBC-2.3.9.tar.gz; cd unixODBC-2.3.9/; ./configure --prefix=/usr/local/unixODBC; make; make install

ENV PATH "$PATH:/usr/local/unixODBC/bin"

#install htslib
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2; \
    tar -vxjf samtools-1.15.1.tar.bz2; rm samtools-1.15.1.tar.bz2; \
    cd samtools-1.15.1; \
    ./configure; make; make install
RUN cd /opt; wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2; \
    tar -vxjf htslib-1.15.1.tar.bz2; rm htslib-1.15.1.tar.bz2; \
    cd htslib-1.15.1; \
    ./configure; make; make install

ENV PATH "$PATH:/opt/samtools-1.12"

#build Berkley db
RUN cd /opt; wget https://github.com/berkeleydb/libdb/releases/download/v5.3.28/db-5.3.28.tar.gz; tar -xf db-5.3.28.tar.gz; cd db-5.3.28/; cd build_unix; ../dist/configure; make; make install

#build boost libraries
RUN cd /opt; wget https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.gz; tar -xf boost_1_71_0.tar.gz; cd boost_1_71_0/; ./bootstrap.sh; ./b2

#install fastcgi
RUN apt-get install -y apache2 libapache2-mod-fcgid; apt-get update; apt-get -y upgrade; apt-get install -y apache2-api-20120211; cd /opt; wget https://mirrors.edge.kernel.org/ubuntu/pool/multiverse/liba/libapache-mod-fastcgi/libapache2-mod-fastcgi_2.4.7~0910052141-1.2_amd64.deb; dpkg -i libapache2-mod-fastcgi_2.4.7~0910052141-1.2_amd64.deb

#install libzip2
RUN apt-get install -y libbz2-dev cmake

ENV PATH "$PATH:/opt/boost_1_71_0:/opt/boost_1_71_0/stage/lib"

#libjpeg-turbo
RUN cd /opt; git clone https://github.com/libjpeg-turbo/libjpeg-turbo.git; cd libjpeg-turbo/; cmake -G"Unix Makefiles"; make;
#install libpgn, libtiff-dev, libungif-bin, libxml2, libxslt-dev
RUN apt-get install -y zlib1g-dev libpng-dev libtiff-dev libungif-bin libxml2 libxslt-dev

#install lzo
RUN cd /opt; wget https://www.oberhumer.com/opensource/lzo/download/lzo-2.10.tar.gz; tar -xf lzo-2.10.tar.gz; cd lzo-2.10; ./configure; make; make check; make test; make install

#install pcre
RUN cd /opt; wget https://sourceforge.net/projects/pcre/files/pcre/8.43/pcre-8.43.tar.gz; tar -xf pcre-8.43.tar.gz; cd pcre-8.43; ./configure; make; make check; make test; make install

#install pcre again
RUN /usr/local/src; wget https://sourceforge.net/projects/pcre/files/pcre/8.45/pcre-8.45.tar.gz; tar -xf pcre-8.45.tar.gz; cd pcre-8.45; ./configure --prefix=/usr/local/src/pcre-8.45; make; make install; ln -s /usr/local/src/pcre-8.45 /usr/local/pcre

ENV PATH "$PATH:/usr/local/pcre/bin"

#install sqlite
RUN cd /opt; DEBIAN_FRONTEND="noninteractive" apt-get install -y tcl tclsh tcl8.6-dev; git clone https://github.com/sqlite/sqlite.git; mkdir bld; cd bld;  ../sqlite/configure; make; make sqlite3.c

#skipping sybase because its not free

RUN apt-get -y update; apt-get -y upgrade; apt-get install -y libboost-all-dev sqlite3 sqlitebrowser libsqlite3-dev libzstd-dev libcurl4-openssl-dev

#reinstall Berkeley db. This one is actually necessary, the first install probably isn't but I just don't want to mess with it at this point
RUN cd /opt; wget http://download.oracle.com/berkeley-db/db-6.2.32.tar.gz; tar -xf db-6.2.32.tar.gz; cd db-6.2.32; cd build_unix/; ../dist/configure --prefix=/usr --enable-compat185 --enable-dbm --disable-static --enable-cxx; make; make docdir=/usr/share/doc/db-db-6.2.32 install; chown -v -R root:root /usr/bin/db_* /usr/include/db_185.h /usr/include/db_cxx.h /usr/lib/libdb*.so /usr/lib/libdb*.la

#install ncbic++ tools
RUN cd /opt; wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/ARCHIVE/2021/Sep_30_2021/ncbi_cxx--25_2_0.tar.gz; tar -xf ncbi_cxx--25_2_0.tar.gz; cd ncbi_cxx--25_2_0; ./configure --without-debug --with-optimization; cd GCC750-ReleaseMT64/build; make all_r

#make prosplign
RUN cd /opt/ncbi_cxx--25_2_0/GCC750-ReleaseMT64/build/algo/align/prosplign/demo; make

ENV PATH "$PATH:/opt/ncbi_cxx--25_2_0/GCC750-ReleaseMT64/bin"

# Install compart algos
RUN mkdir -p /compart_algos;
ADD compart_algos /compart_algos
RUN chmod 777 /compart_algos/*

ENV PATH "$PATH:/compart_algos"

#install diamond
RUN cd /opt; wget http://github.com/bbuchfink/diamond/releases/download/v2.0.14/diamond-linux64.tar.gz; tar xzf diamond-linux64.tar.gz
ENV PATH "$PATH:/opt/"

#install infernal
RUN cd /opt; wget http://eddylab.org/infernal/infernal-1.1.4-linux-intel-gcc.tar.gz; tar -xf infernal-1.1.4-linux-intel-gcc.tar.gz; cd infernal-1.1.4-linux-intel-gcc; ./configure; make; make check; make install

#tRNAscan-se
RUN cd /opt; wget http://trna.ucsc.edu/software/trnascan-se-2.0.10.tar.gz; tar -xf trnascan-se-2.0.10.tar.gz; cd tRNAscan-SE-2.0/; ./configure; make; make install

#upgrade cmake for skesa as it needs at least version 3.16 and the default is 3.10
RUN apt-get remove -y cmake; cd /opt; wget https://github.com/Kitware/CMake/releases/download/v3.25.3/cmake-3.25.3-linux-x86_64.tar.gz; tar -xf cmake-3.25.3-linux-x86_64.tar.gz
ENV PATH "$PATH:/opt/cmake-3.25.3-linux-x86_64/bin"

#skesa
RUN cd /opt; git clone https://github.com/ncbi/SKESA.git; cd SKESA; make
ENV PATH "$PATH:/opt/SKESA/"

RUN cd opt; wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz; tar xvfz seqkit_linux_amd64.tar.gz; cp seqkit /usr/local/bin/

RUN apt-get remove -y libpulse0

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts

ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
