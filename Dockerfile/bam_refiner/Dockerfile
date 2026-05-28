FROM ubuntu:20.04
LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="BAM Refiner container for PRCGAP pipeline"

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

WORKDIR /tools

RUN apt-get update && apt-get install -y \
    openssh-client \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    curl \
    python \
    python3 \
    python3-pip

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

RUN wget https://github.com/yos-sk/bam_refiner/archive/refs/tags/v0.3.6.tar.gz && \
    tar zxf v0.3.6.tar.gz && \
    cd  bam_refiner-0.3.6 && \
    ~/.cargo/bin/cargo build --release
RUN rm v0.3.6.tar.gz

RUN git clone https://github.com/yos-sk/split_bam.git && \
    cd split_bam && \
    ~/.cargo/bin/cargo build --release

RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 && \
    tar jxvf minimap2-2.28_x64-linux.tar.bz2 

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

RUN wget https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xJf meryl-1.4.1.Linux-amd64.tar.xz

ENV PATH $PATH:/tools/meryl-1.4.1/bin

ENV PATH $PATH:/tools/bam_refiner-0.3.6/target/release
ENV PATH $PATH:/tools/split_bam/target/release
ENV PATH $PATH:/tools/minimap2-2.28_x64-linux

CMD ["/bin/bash"]