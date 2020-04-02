FROM ubuntu:16.04
LABEL maintainer="Suzanne Jin <suzanne.jin@crg.es>"

ARG CACHEBUST=1

## Basic tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ed \
      less \
      vim-tiny \
      wget \
      git \
      python3 \ 
      build-essential \
      cmake \
      curl \
      libargtable2-0 \
      ca-certificates \
      libpng12-0 \
      python3-biopython \
      python3-numpy \ 
      python3-setuptools \
      default-jdk \
      libpng-dev \
      tcsh \ 
      unzip


# Change working directory, where the softwares will be installed
WORKDIR /usr/src/app

# Install PSIPRED
RUN git clone https://github.com/psipred/psipred.git && \
    checkout 5a5a6157c9011d30f43d33ffa41018c94025805c && \
    cd psipred && \
    cd src && \
    make && \
    make install 

# Install DeepMSA
RUN wget https://zhanglab.ccmb.med.umich.edu/DeepMSA/package.zip && \
    unzip package.zip && \
    rm package.zip && \
    cd hhsuite2 && \
    cat lib/hh/scripts/HHPaths.pm | sed "s/\/cluster\/toolkit\/production\/bioprogs\/psipred/\/usr\/src\/app\/psipred/g" | sed 's/\/cluster\/toolkit\/production\/bioprogs\/blast\/bin/\/usr\/src\/app\/hhsuite2\/bin/g'> a && \
    cat a > lib/hh/scripts/HHPaths.pm && \
    rm a

# Install QuantTest2 (original version)
# RUN wget http://bioinf.ucd.ie/quantest2.tar && \
#     tar -xvf quantest2.tar && \
#     rm quantest2.tar 

# Add modified scripts
RUN git clone https://github.com/suzannejin/mod_QuanTest2.git && \ 
    mv /usr/src/app/mod_QuanTest2/bin/runpsipred /usr/src/app/psipred/. && \
    cp /usr/src/app/mod_QuanTest2/bin/quantest2Local.py /usr/src/app/quantest2.py && \
    cp /usr/src/app/mod_QuanTest2/bin/quantest2Local.py /bin/quantest2 && \
    rm -rf mod_QuanTest2 


ENV HHLIB=/usr/src/app/hhsuite2

ENTRYPOINT ["quantest2"]
