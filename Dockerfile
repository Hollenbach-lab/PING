FROM ubuntu:22.04

ENV TZ=UTC
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    gzip unzip gdebi-core wget \
    libtbb-dev libssl-dev libcurl4-openssl-dev \
    samtools bcftools \
    libxml2-dev pandoc git less curl bzip2 \
    build-essential gfortran tzdata \
    zlib1g zlib1g-dev libbz2-dev liblzma-dev \
    libpcre2-dev libudunits2-dev librsvg2-dev \
    liblapack3 libfreetype6-dev libfribidi-dev \
    libharfbuzz-dev libfontconfig1-dev libicu-dev \
    libjpeg-dev libpng-dev libtiff-dev make \
    libgit2-dev \
    libv8-dev libnode-dev libnlopt-dev \
    libx11-dev libxt-dev libcairo2-dev \
    libxrender-dev \
    pkg-config \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install R 4.3.1 from source, matching the Singularity recipe
WORKDIR /tmp
RUN wget -q https://cran.r-project.org/src/base/R-4/R-4.3.1.tar.gz \
    && tar -xf R-4.3.1.tar.gz \
    && cd R-4.3.1 \
    && ./configure --prefix=/opt/local/R --with-readline=no --with-x=no --enable-R-shlib --enable-BLAS-shlib --enable-LAPACK-shlib \
    && make -j"$(nproc)" \
    && make install \
    && make clean \
    && cd /tmp \
    && rm -rf /tmp/R-4.3.1 /tmp/R-4.3.1.tar.gz

RUN /opt/local/R/bin/R --slave -e '\
packages <- c( \
  "data.table","stringr","plotly","R.utils","gtools","zip", \
  "cluster","glue","argparser","pandoc","DescTools" \
); \
install.packages(packages, repos="https://cloud.r-project.org"); \
missing <- packages[!packages %in% installed.packages()[,"Package"]]; \
if(length(missing) > 0) stop(paste("Missing packages:", paste(missing, collapse=", ")))'

# Install bowtie2 2.4.1 from source, matching the Singularity recipe
WORKDIR /tmp
RUN wget -q -O bowtie2-source.zip \
      https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-source.zip/download \
    && unzip -q bowtie2-source.zip \
    && cd /tmp/bowtie2-2.4.1 \
    && sed -i 's/“/\"/g' processor_support.h \
    && sed -i 's/”/\"/g' processor_support.h \
    && make -j"$(nproc)" \
    && rm -f /tmp/bowtie2-source.zip

# Install PHASE, matching the Singularity recipe
RUN mkdir -p /opt/phase \
    && cd /opt/phase \
    && wget -q http://stephenslab.uchicago.edu/assets/software/phase/phasecode/phase.2.1.1.linux.tar.gz \
    && tar -xf phase.2.1.1.linux.tar.gz \
    && chown -R root:root /opt/phase/phase.2.1.1.linux \
    && chmod 755 /opt/phase/phase.2.1.1.linux \
    && chmod 755 /opt/phase/phase.2.1.1.linux/PHASE \
    && rm -f /opt/phase/phase.2.1.1.linux.tar.gz

ENV PATH="/opt/local/R/bin:/tmp/bowtie2-2.4.1:/opt/phase/phase.2.1.1.linux:${PATH}"