# PING docker
An R-based bioinformatic pipeline to determine killer-cell immunoglobulin-like receptor (KIR) copy number and genotypes from short-read sequencing data.

This branch of PING is setup to be run using snakemake [https://snakemake.readthedocs.io/]


## Language
R


## System compatibility
* Linux (tested on Ubuntu, CentOS)
* OS X
* Windows (untested)


## System dependencies download and install
* conda 
* snakemake 

You can download and install conda [https://conda.io/en/latest/miniconda.html](here)

### 1. Install conda

```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3.sh  && bash ~/miniconda3.sh -b  && . ~/miniconda3/etc/profile.d/conda.sh && conda activate base && conda install mamba -n base -c conda-forge
```

### 2. Install snakemake 

```shell
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```


### 3. Download PING and run the test example


```shell
git clone https://github.com/wesleymarin/PING.git --single-branch --branch reborn_docker_copyOnly
```

### 4. Run the test example

```shell
cd PING/workflow && \
conda activate snakemake && \
snakemake -j 64 --use-conda --verbose --conda-frontend mamba # run snakemake with (up to) 64 threads, using conda for dependency management, and mamba as the conda frontend
```

### 5. Run on any of the 1000 genome samples

```shell
cd PING/workflow && \
conda activate snakemake && \
snakemake -j 64 --use-conda --verbose --conda-frontend mamba ../output/tgp/NA12004.tar.gz ../output/tgp/HG04200.tar.gz 
```

### 6. Run on your own data (starting with bam files)

Add all your bam files to the directory (within the `PING` directory) `input/data/bam`. Each file in the directory `input/data/bam` should correspond to a single sample, and have a file name with the pattern `{sample}.bam`, (do not use underscores in the file name).  To run the two samples `foo.bam` and `bar.bam`

```shell
cp /some_other_dir/foo.bar /some_other_dir/bar.bam input/data/bam/ &&\
cd PING/workflow && \
conda activate snakemake && \
snakemake -j 64 --use-conda --verbose --conda-frontend mamba ../output/data/foo.tar.gz ../output/data/bar.tar.gz 
```

snakemake will also automatically run PING for every bam file in the `input/data/bam` directory:

```shell
cp /some_other_dir/foo.bar /some_other_dir/bar.bam input/data/bam/ &&\
cd PING/workflow && \
conda activate snakemake && \
snakemake -j 64 --use-conda --verbose --conda-frontend mamba 
```
