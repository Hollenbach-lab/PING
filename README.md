# PING
An R-based bioinformatic pipeline to determine killer-cell immunoglobulin-like receptor (KIR) copy number and genotypes from short-read sequencing data.

## Language
R

## System compatibility
* Linux (tested on Ubuntu, CentOs)
* Possible on OS X, but untested
  
## System dependencies download and install
* bowtie2 (tested with version 2.3.4.1) https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/
* samtools (tested with version 1.7) https://sourceforge.net/projects/samtools/files/samtools/1.7/
* bcftools (tested with version 1.7) https://sourceforge.net/projects/samtools/files/samtools/1.7/
* R (tested wtih version 3.6.3) https://cran.r-project.org/src/base/R-3/R-3.6.3.tar.gz
* RStudio (used for running the script) https://rstudio.com/products/rstudio/download/
  - Either RStudio desktop free version, or RStudio server free version
  - It is possible to run the script without RStudio, but copy number thresholding will be more difficult

### Ubuntu 20.04 command-line install of dependencies
`sudo apt install bowtie2 gzip samtools bcftools r-base gdebi-core wget`
#### Download RStudio
`wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.3.1093-amd64.deb`
#### Install RStudio
`sudo gdebi rstudio-1.3.1093-amd64.deb`

## Required R packages:    
* data.table 
* stringr
* pryr 
* plotly 
* gtools 
* R.utils
 
### R dependency console install command
`install.packages(c("plotly","stringr","pryr","plotly","gtools","R.utils"),dependencies = T)`

## Setting up pipeline
### Downloading pipeline code
Option 1: Download zip file `wget https://github.com/wesleymarin/PING/archive/master.zip`

Option 2: Clone repository `git clone https://github.com/wesleymarin/PING.git`

### Environment setup
1. Change line 2 of PING_run.R to the local directory where you download PING
  - `2 setwd('/home/LAB_PROJECTS/PING_PAPER/PING')` Set this to your own PING working directory
2. Change lines 18-23 to fit your data/environment
  - `18 rawFastqDirectory <- '/home/LAB_PROJECTS/PING_PAPER/PING/test_sequence/'` Set to raw sequence directory or extracted fastq directory if extraction has already been performed
  - `19 fastqPattern <- '_KIR_'` Use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
  - `20 threads <- 4` Number of threads to use during bowtie2 alignments
  - `21 resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_test_sequence_results/'` Set the results directory, one will be created if it does not already exist (all pipeline output will be recorded here)
  - `22 shortNameDelim <- '_'` Set a delimiter to shorten sample ID's (ID will be characters before delim)
  - `23 minDP <- 10` Minimum depth for calling variants (set lower if using low-depth data, the default of 10 should work for most data)

## Running PING
Open PING_run.R in Rstudio
Set environment variables outlined in **Environment setup**
Run the script from top to bottom


Seeing the following message in the copy graphing module is normal:

'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'

### Running included test data
We have included 5 test sequences to run through the pipeline, they are located in the test_sequence/ directory. These samples were picked to cover a range of KIR haplotypes.


The default input settings (lines 18-23) are setup to run on this data, but must still be modified to fit your local environment. For example, line 18 `rawFastqDirectory` must be set to where test_sequence/ is located on your local machine. Additionally, line 21 `resultsDirectory` must be set to where you want the results output to be written.


The copy thresholding functionality will not work well for such a small cohort, but the pipeline can be run with only gene presence/absence information.

## Troubleshooting
Please save a copy of your R Console output and contact me through github or email at wesley.marin@ucsf.edu
