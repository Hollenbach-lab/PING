
setwd('/home/LAB_PROJECTS/PING2_PAPER/PING2') #Set this to your own PING2 working directory

# ---- DEPENDENCIES ----
' if any dependencies are missing, install with
install.packages("plotly",dependencies = T)
'
library(data.table)
library(stringr)
library(methods)
library(pryr)
library(plotly)
library(R.utils)
library(gtools)


# Initialization variables ------------------------------------------------
rawFastqDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/testSequence//' # can be set to raw sequence or extractedFastq directory
fastqPattern <- '_KIR_' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 30
resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_testSequence_results/' # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)
minDP <- 10


source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('Resources/ping_copy.R') # do not change
source('Resources/ping_allele.R') # do not change
source('Resources/ping_gc_align.R') # do not change
source('Resources/alleleCombine_functions.R') # do not change

# Preparation -------------------------------------------------------------
# Build up a list of sample objects
sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, resultsDirectory, shortNameDelim) # no need to change



# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING2 KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq') # no need to change
# Run PING2 extractor
sampleList <- extractor.run(sampleList,threads,extractedFastqDirectory,forceRun=F) # set forceRun=T if you want to force alignments


# PING2 gene content and copy number --------------------------------------
cat('\n\n----- Moving to PING2 gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=resultsDirectory,forceRun=F,onlyKFF=F) # set forceRun=T if you want to force alignments
sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=resultsDirectory) # this function sets copy thresholds

## Seeing the following message when generating copy plots is normal
'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'

## use this section for notes / recording suspicious samples
' 
KIR3DP1 [example]
  IND00001 [example]
'

# PING2 alignments and allele calling ----------------------------------------------
# Iter align workflow
source('Resources/genotype_alignment_functions.R') # do not change

# Alignment and allele calling workflow
for(currentSample in sampleList){
  currentSample <- ping_iter.run_alignments(currentSample, threads)
  currentSample <- ping_iter.allele(currentSample)
  
  currentSample <- ping_filter.run_alignments(currentSample, threads)
  currentSample <- ping_filter.allele(currentSample)
  currentSample <- post.combineGenos(currentSample, resultsDirectory)
}

