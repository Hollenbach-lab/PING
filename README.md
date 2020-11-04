# PING
Use with the PING_run.R script.
The run script is meant to be run from top to bottom.

Change line 2 to wherever you download PING -> setwd('/home/LAB_PROJECTS/PING2_PAPER/PING')

Change lines 18-23 to fit your data/environment

Seeing the following message in the copy graphing module is normal:

'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'


PING submission notes

WRITE DOCUMENTATION

# PING
An R-based bioinformatic pipeline to determine killer-cell immunoglobulin-like receptor (KIR) copy number and genotypes from short-read sequencing data.

## Language, System compatibility
  R
  Linux (tested on Ubuntu, CentOs)
  Possible on OS X, but untested
  
## Dependencies download and install
  `sudo apt install bowtie2 gzip samtools bcftools`
  bowtie2
  samtools
  bcftools

## Required R packages:
    
    `install.packages(c("plotly","stringr","pryr","plotly","gtools","R.utils"),dependencies = T)`
    
    data.table
    stringr
    pryr
    plotly
    gtools
    R.utils

