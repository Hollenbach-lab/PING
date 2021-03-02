
FROM rocker/r-base:3.6.3

WORKDIR /usr/home

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    gzip \
    gdebi-core \
    wget \
    bedtools \
    libtbb-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    samtools \
    bcftools \
    libxml2-dev \
    pandoc

RUN wget --no-check-certificate https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip/download
RUN unzip download
RUN cd bowtie2-2.3.4.1 && make install
RUN cd && bowtie2 --help && rm -rf download
RUN mkdir -p /opt/PING/

## copy files
COPY Resources/ /opt/PING/Resources/
COPY test_sequence/ /opt/PING/test_sequence/
COPY PING_run.R /opt/PING/PING_run.R
COPY install_packages.R /opt/PING/install_packages.R

RUN Rscript /opt/PING/install_packages.R

#CMD Rscript PING_run.R
