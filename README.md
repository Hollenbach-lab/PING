# PING (Pushing Immunogenetics to the Next Generation)
An R-based bioinformatic pipeline to determine killer-cell immunoglobulin-like receptor (KIR) copy number and high-resolution genotypes from short-read sequencing data.

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

## Ubuntu 20.04 command-line install of dependencies
#### System dependencies
`sudo apt install bowtie2 gzip samtools bcftools r-base gdebi-core wget libssl-dev libcurl4-openssl-dev`
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

#### R dependency console install command
`install.packages(c("data.table","plotly","stringr","pryr","gtools","R.utils"),dependencies = T)`

## Setting up pipeline
#### Downloading pipeline code
Option 1: Download zip file `wget https://github.com/wesleymarin/PING/archive/master.zip ; unzip master.zip`

Option 2: Clone repository `git clone https://github.com/wesleymarin/PING.git`

#### Environment setup
Change lines 33-39 to fit your data/environment, this does not need to be changed to run the included example dataset.
  - `33 rawFastqDirectory <- 'test_sequence/'` Set to raw sequence directory or extracted fastq directory if extraction has already been performed
  - `34 fastqPattern <- 'fastq'` Use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
  - `35 threads <- 4` Number of threads to use during bowtie2 alignments
  - `36 resultsDirectory <- '3_test_sequence_results/'` Set the results directory, one will be created if it does not already exist (all pipeline output will be recorded here)
  - `37 shortNameDelim <- ''` Set a delimiter to shorten sample ID's (ID will be characters before delim)
  - `38 setup.minDP <- 8` Minimum depth for calling variants used for genotype-matched references (set lower if using low-depth data, the default of 8 should work for most data)
  - `39 final.minDP <- 20` Minimum depth for calling variants used for final genotype determination (set lower if using low-depth data, the default of 20 should work for most data)

## Running PING
Open PING_run.R in Rstudio

Set environment variables outlined in **Environment setup**. Nothing needs to be changed to run the included example dataset.

Run the script from top to bottom

Seeing the following message in the copy graphing module is normal:

`'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'`

## Running included test data
We have included 5 test sequences to run through the pipeline, they are located in the test_sequence/ directory. These samples were picked to cover a range of KIR haplotypes.

The default input settings (lines 33-39) are setup to run on this data without modification.

The copy thresholding functionality will not work well for such a small cohort, but we included preset thresholds for the example dataset.

## Running your own data
Update line 33 `rawFastqDirectory <- 'test_sequence/'` to point to your own data directory.

Update line 34 `fastqPattern <- 'fastq'` to match your data naming if necessary. For example, if your sequencing data is named [SAMPLE_ID]\_R1_fq.gz, you would change line 34 to `fastqPattern <- 'fq'`.

Change line 84 option `use.threshFile=T` to `use.threshFile=F`, this will cause PING to prompt for copy number thresholding. If using Rstudio, copy number plots will be displayed in the 'Plots' panel, these plots can also be found in \[resultsDirectory\]/copyPlots/\[GENE\]\_copy\_number\_plot.html as html files to be opened with a web browser. 

It is recommended to only run the script to this line first if you are using your own data, then set the copy number thresholds, then run the rest of the script. 

## PING output
Copy number output can be found at `[resultsDirectory]/manualCopyNumberFrame.csv`

Genotype output can be found at `[resultsDirectory]/finalAlleleCalls.csv`

Aligned SNP tables can be found in `[resultsDirectory]/alignmentFiles/[sampleID]/iterAlign/`

Copy number graphs can be found in `[resultsDirectory]/copyPlots/`

#### Unresolved genotypes
If PING is unable to perfectly match aligned SNPs to known KIR allele sequences an unresolved call will be produced.

Unresolved genotype information can be found in `[resultsDirectory]/iterAlleleCalls.csv`, where the closest allele match is recorded along with the mismatched SNP information in the following format:
`[closest_matched_allele]$[exon]_[position].[nucleotide]`

Where closest matched allele is the allele genotyping that best matches the aligned SNPs, nucleotide denotes the mismatched nucleotide located at the indicated exon and position within the exon. Multiple mismatched SNPs are connected with the `^` symbol.

## Troubleshooting
Please save a copy of your R Console output and contact me through github or email at wesley.marin@ucsf.edu
