FROM ubuntu:20.04

LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="NanoMonSV postprocess container for PRCGAP pipeline (postprocess, connect, merge)"

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
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

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH $PATH:/root/.cargo/bin

WORKDIR /opt

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

# Install nanomonsv_postprocess
RUN wget https://github.com/yos-sk/nanomonsv_postprocess/archive/refs/tags/v0.2.5.tar.gz && \
    tar zxf v0.2.5.tar.gz && \
    cd nanomonsv_postprocess-0.2.5 && \
    cargo build --release
RUN rm v0.2.5.tar.gz

ENV PATH="/opt/nanomonsv_postprocess-0.2.5/target/release:$PATH"

CMD ["/bin/bash"]
