# rocker/tidyverse ships R with tidyverse already built in, the Posit Public
# Package Manager preconfigured (binary CRAN packages, correct HTTPUserAgent),
# a compiler toolchain incl. gfortran, and the graphics-stack system libraries.
# This avoids both the slow tidyverse source build and installing R via the CRAN
# apt repository (which broke once Ubuntu focal reached EOL). The bioinformatics
# tools are added on top.
FROM rocker/tidyverse:4.5.2

LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="Copy number analysis container for PRCGAP pipeline"

# Avoid interactive prompts (e.g. tzdata asking for a geographic area) during apt installs.
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# System build tools and libraries for samtools/minimap2/bedtools, plus the
# graphics-stack dev libraries the R plotting packages link against.
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses-dev \
        libcurl4-openssl-dev \
        python3 \
        python3-pip \
        libxml2-dev libfontconfig1-dev libfreetype6-dev libpng-dev \
        libtiff-dev libjpeg-dev libwebp-dev libharfbuzz-dev libfribidi-dev libssl-dev \
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

# Install bedtools. Use the official prebuilt static binary: bedtools 2.31.0
# fails to compile from source under newer GCC (missing #include <cstdint>),
# and the static binary avoids the build entirely.
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static \
        -O /usr/local/bin/bedtools && \
    chmod +x /usr/local/bin/bedtools

# Install Python dependencies (Debian's Python is externally managed -> --break-system-packages).
RUN pip3 install --no-cache-dir --break-system-packages pysam

# Helvetica-compatible font for the plots. Helvetica itself is proprietary (no
# free package); fonts-urw-base35 provides Nimbus Sans, and fontconfig ships a
# built-in "Helvetica -> Nimbus Sans" alias, so plot_copy_number.R's
# family="Helvetica" resolves deterministically (consistent across rebuilds/runs)
# without changing the script.
RUN apt-get update && apt-get install -y --no-install-recommends \
        fonts-urw-base35 fontconfig && \
    fc-cache -f && \
    rm -rf /var/lib/apt/lists/*

# Install the remaining R packages. tidyverse (and its deps dplyr, scales) are
# already present in the base image; these come as precompiled binaries from the
# Posit Public Package Manager. DNAcopy is built from source via Bioconductor.
RUN install2.r --error --skipinstalled \
        optparse patchwork slider viridis ggnewscale BiocManager && \
    R -e "BiocManager::install('DNAcopy', update = FALSE, ask = FALSE)"

# Verify every R library used by the copy number scripts loads; fail the build
# if any is missing (keep this list in sync with the scripts' library() calls).
RUN R -e "pkgs <- c('DNAcopy','tidyverse','optparse','scales','dplyr','patchwork','slider','viridis','ggnewscale'); \
          ok <- sapply(pkgs, requireNamespace, quietly = TRUE); \
          if (any(!ok)) { stop('Missing R packages: ', paste(pkgs[!ok], collapse = ', ')); }; \
          cat('All R packages installed:', paste(pkgs, collapse = ', '), '\n')"

# Install mosdepth.
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.11/mosdepth \
        -O /usr/local/bin/mosdepth && \
    chmod +x /usr/local/bin/mosdepth

# Install seqtk. No release binary is published, so build from source; the
# Makefile's install target puts seqtk in /usr/local/bin.
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v1.4.tar.gz -O seqtk-1.4.tar.gz && \
    tar zxf seqtk-1.4.tar.gz && \
    cd seqtk-1.4 && make && make install && \
    cd .. && rm -rf seqtk-1.4 seqtk-1.4.tar.gz

ENV PATH="/usr/local/bin:${PATH}"
ENV PATH="/opt/minimap2-2.28_x64-linux:${PATH}"

CMD ["/bin/bash"]
