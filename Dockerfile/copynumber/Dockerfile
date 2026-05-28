FROM ubuntu:20.04

LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="Copy number analysis container for PRCGAP pipeline"


RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

# Install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 && \
    tar jxvf minimap2-2.28_x64-linux.tar.bz2 

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/archive/v2.31.0.tar.gz && \
    tar zxvf v2.31.0.tar.gz && \
    rm v2.31.0.tar.gz && \
    cd bedtools2-2.31.0 && \
    make && \
    cp bin/* /usr/local/bin && \
    cd / && \
    rm -rf /opt/bedtools2-2.31.0

# Install Python dependencies
RUN pip3 install --no-cache-dir pysam

ENV PATH="/usr/local/bin:${PATH}"
ENV PATH="/opt/minimap2-2.28_x64-linux:${PATH}"

CMD ["/bin/bash"]
