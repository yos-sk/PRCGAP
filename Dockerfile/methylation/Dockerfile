FROM ubuntu:20.04

LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="Methylation calling container for PRCGAP pipeline"

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    wget \
    ca-certificates \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

# Install pb-CpG-tools (aligned_bam_to_cpg_scores)
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz && \
    tar zxvf pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz && \
    rm pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
ENV PATH="/opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin:${PATH}"

# Install modkit (precompiled binary for Linux x86_64)
# Extracted directory has a commit-hash suffix (dist_modkit_v0.4.1_<hash>/),
# so resolve the binary path with a glob and copy it to /usr/local/bin.
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.4.1/modkit_v0.4.1_u16_x86_64.tar.gz && \
    tar zxvf modkit_v0.4.1_u16_x86_64.tar.gz && \
    rm modkit_v0.4.1_u16_x86_64.tar.gz && \
    install -m 755 dist_modkit_v0.4.1_*/modkit /usr/local/bin/modkit && \
    rm -rf dist_modkit_v0.4.1_*

CMD ["/bin/bash"]
