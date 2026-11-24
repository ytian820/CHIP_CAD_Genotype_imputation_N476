# CHIP-CAD genotype QC and imputation tools; independent linux/amd64 build.
# Build: docker build --platform linux/amd64 -t ytian820/chip_cad_genotype:1.0.0 docker

# Part 1. Build BCFtools and HTSlib from official release sources.
FROM ubuntu:24.04 AS builder
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates curl bzip2 make gcc libc6-dev \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /tmp/build
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/1.25/bcftools-1.25.tar.bz2 \
        | tar -xj \
    && cd bcftools-1.25 \
    && ./configure --prefix=/usr/local \
    && make -j2 \
    && make install \
    && cd htslib-1.25 \
    && make -j2 bgzip tabix \
    && install bgzip tabix /usr/local/bin/

# Part 2. Install runtime tools and R packages on a clean Ubuntu base.
FROM ubuntu:24.04
LABEL org.opencontainers.image.title="CHIP_CAD_genotype" \
      org.opencontainers.image.version="1.0.0"
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates curl unzip python3 \
        zlib1g libbz2-1.0 liblzma5 libcurl4t64 libssl3t64 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/ /usr/local/
RUN curl -fsSL https://s3.amazonaws.com/plink2-assets/beta1/plink2_linux_x86_64_20261107.zip \
        -o /tmp/plink2.zip \
    && unzip -j /tmp/plink2.zip plink2 -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        -o /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && echo 'deb [arch=amd64 signed-by=/etc/apt/trusted.gpg.d/cran_ubuntu_key.asc] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/' \
        > /etc/apt/sources.list.d/cran.list \
    && curl -fsSL https://r2u.stat.illinois.edu/ubuntu/dirk_eddelbuettel_pubkey.asc \
        -o /etc/apt/trusted.gpg.d/r2u_key.asc \
    && echo 'deb [arch=amd64 signed-by=/etc/apt/trusted.gpg.d/r2u_key.asc] https://r2u.stat.illinois.edu/ubuntu noble main' \
        > /etc/apt/sources.list.d/r2u.list \
    && printf 'Package: *\nPin: release o=CRAN-Apt Project\nPin-Priority: 700\n' \
        > /etc/apt/preferences.d/99r2u \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        r-base-core r-cran-tidyverse r-cran-data.table \
    && rm -rf /var/lib/apt/lists/*

# Part 3. Verify tools, plugins and R packages without research data.
RUN bcftools --version \
    && bcftools plugin -l | grep -Fx fixref \
    && bgzip --version \
    && tabix --version \
    && plink2 --version \
    && python3 --version \
    && Rscript -e 'library(tidyverse); library(data.table); print(R.version.string); print(packageVersion("tidyverse")); print(packageVersion("data.table"))'
CMD ["bash"]
