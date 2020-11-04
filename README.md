# PING
Use with the PING_run.R script.
The run script is meant to be run from top to bottom.

Change line 2 to wherever you download PING -> setwd('/home/LAB_PROJECTS/PING2_PAPER/PING')

Change lines 18-23 to fit your data/environment

Seeing the following message in the copy graphing module is normal:

'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'



# PING
An R-based bioinformatic pipeline to determine killer-cell immunoglobulin-like receptor (KIR) copy number and genotypes from short-read sequencing data.

## Language
R

## System compatibility
Linux (tested on Ubuntu, CentOs)
Possible on OS X, but untested
  
## System dependencies download and install
* bowtie2 (tested with version 2.3.4.1) https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/
* samtools (tested with version 1.7) https://sourceforge.net/projects/samtools/files/samtools/1.7/
* bcftools (tested with version 1.7) https://sourceforge.net/projects/samtools/files/samtools/1.7/

### Command-line install of dependencies
`sudo apt install bowtie2 gzip samtools bcftools`

## Required R packages:    
* data.table 
* stringr
* pryr 
* plotly 
* gtools 
* R.utils
 
### R console dependency install command
`install.packages(c("plotly","stringr","pryr","plotly","gtools","R.utils"),dependencies = T)`

