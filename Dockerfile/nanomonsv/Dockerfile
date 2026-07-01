# syntax=docker/dockerfile:1
#
# nanomonsv image for PRCGAP (paper setup: RepeatMasker 4.1.8 + Dfam 3.9).
#
# Dependencies mirror the bioconda `nanomonsv=0.8.0` meta-package, but RepeatMasker
# is pinned to 4.1.8 and nanomonsv itself is installed via pip (so the conda
# metapackage does not pull/pin its own RepeatMasker).
#
# The Dfam download and the RepeatMasker configure are kept as SEPARATE RUN
# layers: editing the configure step never invalidates the (large, cached)
# Dfam download layer above it.
FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Yoshitaka Sakamoto <yosakam2@gmail.com>"

USER root

# --- deps (bioconda nanomonsv=0.8.0 recipe) + RepeatMasker 4.1.8 --------------
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.10 pip \
        htslib \
        "mafft>=7.407" \
        "racon>=1.4.3" \
        "pysam>=0.19.1" \
        "numpy>=1.23.0" \
        "parasail-python>=1.2.4" \
        "minimap2>=2.24" \
        "bwa>=0.7.17" \
        "bedtools>=2.30.0" \
        "repeatmasker=4.1.8" \
        perl wget && \
    micromamba clean -afy
ENV PATH=/opt/conda/bin:$PATH

# --- nanomonsv (pip; no conda RepeatMasker pulled) ---------------------------
RUN pip install --no-cache-dir nanomonsv==0.8.0

# --- Dfam 3.9 famdb download: 0 = root (REQUIRED), 7 = Mammalia (human) -------
# Separate layer from configure. cache mount keeps the large downloads across
# rebuilds even if this layer is invalidated.
RUN --mount=type=cache,target=/opt/dfam_cache \
    RMLIB=/opt/conda/share/RepeatMasker/Libraries/famdb && \
    mkdir -p "$RMLIB" && \
    for p in 0 7; do \
        f=dfam39_full.${p}.h5; \
        [ -f /opt/dfam_cache/$f ] || { \
            wget -q "https://www.dfam.org/releases/Dfam_3.9/families/FamDB/${f}.gz" -O /opt/dfam_cache/${f}.gz && \
            gunzip -f /opt/dfam_cache/${f}.gz; }; \
        cp /opt/dfam_cache/$f "$RMLIB/"; \
    done

# --- RepeatMasker configure (SEPARATE layer; safe to re-edit) -----------------
# Run from inside the install dir and non-interactively (no TTY in Docker).
RUN cd /opt/conda/share/RepeatMasker && \
    perl ./configure \
        -default_search_engine rmblast \
        -trf_prgm    "$(which trf)" \
        -rmblast_dir "$(dirname "$(which rmblastn)")" \
        -libdir /opt/conda/share/RepeatMasker/Libraries \
        < /dev/null

CMD ["/bin/bash"]
