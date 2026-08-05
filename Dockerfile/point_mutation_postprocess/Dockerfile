FROM ubuntu:20.04

LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="Point mutation postprocess container for PRCGAP pipeline"

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    curl \
    git \
    ca-certificates \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

# Install BWA
RUN wget -O bwa-0.7.18.tar.gz https://github.com/lh3/bwa/archive/refs/tags/v0.7.18.tar.gz && \
    tar xzf bwa-0.7.18.tar.gz && \
    cd bwa-0.7.18 && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa-0.7.18*

# Install Python dependencies
RUN pip3 install --no-cache-dir pysam

# Install rust (rustup); place cargo first on PATH so it shadows nothing else
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH=/root/.cargo/bin:$PATH

# Install mutation_postprocess (build under /opt so it matches the PATH below)
WORKDIR /opt
RUN wget https://github.com/yos-sk/mutation_postprocess/archive/refs/tags/v0.1.2.tar.gz && \
    tar zxf v0.1.2.tar.gz && \
    cd mutation_postprocess-0.1.2 && \
    cargo generate-lockfile && \
    cargo update -p hts-sys --precise 2.2.0 && \
    cargo build --release && \
    cd .. && rm v0.1.2.tar.gz


ENV PATH="/opt/mutation_postprocess-0.1.2/target/release:$PATH"

CMD ["/bin/bash"]ls
