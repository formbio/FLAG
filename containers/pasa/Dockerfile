FROM pasapipeline/pasapipeline

RUN apt-get update

#add pasa pipeline perl script to path
ENV PATH "$PATH:/usr/local/src/PASApipeline"

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts


ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
