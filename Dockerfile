
FROM rocker/r-base:3.6.3

WORKDIR /usr/home
ENV CWD=/usr/home



RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    gzip \
    gdebi-core \
    wget \
    libssl-dev \
    libcurl4-openssl-dev \
    bowtie2 \
    samtools \
    bcftools


## copy files
COPY Resources/ Resources/
COPY test_sequence/ test_sequence/
COPY PING_run.R PING_run.R
COPY install_packages.R install_packages.R

RUN Rscript install_packages.R

CMD Rscript PING_run.R