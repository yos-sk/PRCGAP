FROM ubuntu:20.04

LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="Annotation container for PRCGAP pipeline (SV/SNV/INDEL annotation: pysam + samtools + coordconv + transanno)"

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Build deps + python + rust toolchain prerequisites. unzip is needed
# for the transanno release archive; the rest mirror the build deps used
# by Dockerfile/point_mutation_postprocess and Dockerfile/coordconv.
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    unzip \
    build-essential \
    curl \
    ca-certificates \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# samtools (+ htslib bgzip / tabix come along for free).
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    cd .. && rm -rf samtools-1.17*

# pysam (used by add_annotation_*.py, bed2vcf.py, filter_PRCGAP.py,
# reclassify_sv_type.py).
RUN pip3 install --no-cache-dir pysam

# Rust toolchain — required to build coordconv from source.
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH=/root/.cargo/bin:$PATH

WORKDIR /opt

# coordconv is in a private GitHub repo (github.com/yos-sk/coordconv).
# Provide a personal access token via Docker BuildKit secret:
#   docker build --secret id=token,src=token.txt ...
# Mirrors the existing Dockerfile/coordconv/Dockerfile install.
RUN wget https://github.com/yos-sk/coordconv/archive/refs/tags/v0.1.2.tar.gz && \
    tar zxf v0.1.2.tar.gz && \
    cd coordconv-0.1.2 && \
    cargo build --release
RUN rm v0.1.2.tar.gz

ENV PATH="/opt/coordconv-0.1.2/target/release:$PATH"

# transanno v0.4.1 (https://github.com/informationsea/transanno). The
# release page provides prebuilt linux binaries; download the musl
# static-linked tarball so we don't need to recompile. The exact asset
# name on the release page is the only thing likely to drift over time;
# adjust if v0.4.1 layout differs.
RUN wget https://github.com/informationsea/transanno/releases/download/v0.4.1/transanno-x86_64-unknown-linux-musl-v0.4.1.zip && \
    unzip transanno-x86_64-unknown-linux-musl-v0.4.1.zip && \
    mv transanno-x86_64-unknown-linux-musl-v0.4.1/transanno /usr/local/bin/transanno && \
    chmod +x /usr/local/bin/transanno && \
    rm -rf transanno-x86_64-unknown-linux-musl*

# ULTRA 1.2.1 (https://github.com/TravisWheelerLab/ULTRA) — tandem repeat
# annotator. Built from source with CMake (>=3.12 required; focal ships
# 3.16). `make install` puts the `ultra` binary in /usr/local/bin.
RUN wget https://github.com/TravisWheelerLab/ULTRA/archive/refs/tags/1.2.1.tar.gz -O ULTRA-1.2.1.tar.gz && \
    tar zxf ULTRA-1.2.1.tar.gz && \
    cd ULTRA-1.2.1 && cmake . && make && make install && \
    cd .. && rm -rf ULTRA-1.2.1 ULTRA-1.2.1.tar.gz

CMD ["/bin/bash"]
