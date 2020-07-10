
setwd('/home/LAB_PROJECTS/PING2_PAPER/PING2') #Set this to your own PING2 working directory

source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('Resources/ping_copy.R') # do not change

# Initialization variables ------------------------------------------------
rawFastqDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/extractedFastq/' # can be set to raw sequence or extractedFastq directory
fastqPattern <- '_KIR_' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 12
resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/' # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)


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
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=resultsDirectory,forceRun=F, onlyKFF=T) # set forceRun=T if you want to force alignments
sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=resultsDirectory) # this function sets copy thresholds

## Seeing the following message when generating copy plots is normal
'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'

## use this section for notes / recording suspicious samples
' 
KIR3DP1 [example]
  IND00001 [example]
'

# PING2 iter alignments ----------------------------------------------
# Iter align workflow

for(currentSample in sampleList){
  currentSample <- ping_iter.run_alignments(currentSample)
}

# Example wrapper function to run iter alignments for each sample ( should be moved to genotype_alignment_functions.R )
ping_iter.run_alignments <- function( currentSample ){
  
  cat('\nLoading ref DF')
  currentSample <- sampleObj.loadRefDF(currentSample, referenceAlleleDF) # Subset reference allele dataframe by present loci, save to sample object
  cat('\nWriting reference files')
  currentSample <- sampleObj.writeRefFastaBed(currentSample, locusRefList, alignmentFileDirectory) # Write fasta reference file for sample object based on refDF
  currentSample <- sampleObj.iterBowtie2Index(currentSample, bowtie2Build, threads) # Converts fasta file from previous line into a bowtie2 index
  currentSample <- sampleObj.iterBowtie2Align(currentSample, bowtie2, threads, deleteSam=F) # Align sample to bowtie2 index
  currentSample <- sampleObj.iterVCFGen(currentSample, samtools, bcftools, threads) # Convert SAM file into VCF
  
  return( currentSample )
}

# PING2 filter alignments ----------------------------------------------


