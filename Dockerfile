
FROM rocker/r-base:3.6.3

WORKDIR /usr/home

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    gzip \
    gdebi-core \
    wget \
    libtbb-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    samtools \
    bcftools \
    libxml2-dev \
    pandoc

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip/download
RUN unzip download
RUN cd bowtie2-2.3.4.1 && make
ENV PATH=$PATH:/usr/home/bowtie2-2.3.4.1
RUN bowtie2 --help

ENV CWD=/usr/home
ENV RAW_FASTQ_DIR=test_sequence/
RUN mkdir results/
ENV RESULTS_DIR=results

## copy files
COPY Resources/ Resources/
COPY test_sequence/ test_sequence/
COPY PING_run.R PING_run.R
COPY install_packages.R install_packages.R

RUN Rscript install_packages.R

CMD Rscript PING_run.R
