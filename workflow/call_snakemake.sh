#!/bin/bash
export SINGULARITY_DOCKER_USERNAME='gitlab+deploy-token-1563'
export SINGULARITY_DOCKER_PASSWORD=MHooeydaC3NWYyktBTXf
#export SINGULARITY_DOCKER_PASSWORD=1E7ps5aCG-4ZThzErZMU
#1E7ps5aCG-4ZThzErZMU
snakemake --profile ../config/slurm "$@"
