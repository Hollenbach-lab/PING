FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="596f895364093d94dc95ab535d64a2d840d7fd9764851004ee19511b4d48f939"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: ../envs/ping.yaml
#   prefix: /conda-envs/9bde66ce6ea632074aec7d16e81308f2
#   name: ping
#   channels:
#     - r
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bcftools=1.11
#     - bedtools
#     - samtools
#     - bowtie2=2.3.4.1
#     - tbb=2020.2
#     - r=3.6.0
#     - r-stringr
#     - r-plotly
#     - r-r.utils
#     - r-data.table
#     - r-gtools
#     - r-zip
#     - r-pryr
#     - openjdk
#     - bazam
#   prefix: /Users/knoblaun/miniconda3/envs/ping
RUN mkdir -p /conda-envs/9bde66ce6ea632074aec7d16e81308f2
COPY ../envs/ping.yaml /conda-envs/9bde66ce6ea632074aec7d16e81308f2/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/9bde66ce6ea632074aec7d16e81308f2 --file /conda-envs/9bde66ce6ea632074aec7d16e81308f2/environment.yaml && \
    mamba clean --all -y
