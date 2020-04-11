

setwd('/home/LAB_PROJECTS/PING2_PAPER/PING2') #Set this to your own PING2 working directory

source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('ping_copy.R') # do not change

# Initialization variables ------------------------------------------------
rawFastqDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/extractedFastq/' # can be set to raw sequence or extractedFastq directory
fastqPattern <- '_KIR_' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 12
resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/'
shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)


# Preparation -------------------------------------------------------------
# Build up a list of sample objects
sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, resultsDirectory, shortNameDelim)


# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING2 KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq')
# Run PING2 extractor
sampleList <- extractor.run(sampleList,threads,extractedFastqDirectory,forceRun=F)


# PING2 gene content and copy number --------------------------------------
cat('\n\n----- Moving to PING2 gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=resultsDirectory,forceRun=F)
sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=resultsDirectory) #use this function to set copy thresholds

## Seeing the following message when generating copy plots is normal
'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'

## use this section for notes / recording suspicious samples
' 
KIR3DP1 [example]
  IND00001 [example]
'



# PING2 haplotype alignments ----------------------------------------------




