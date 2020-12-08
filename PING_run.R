
setwd('/home/wmarin/ping_reborn/PING/') #Set this to your own PING2 working directory

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
rawFastqDirectory <- '/home/wmarin/african_samples/1_raw_fastq/' # can be set to raw sequence or extractedFastq directory
fastqPattern <- 'fastq' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 36
resultsDirectory <- '/home/wmarin/african_samples/3_validation_gc_align_method/' # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)
setup.hetRatio <- 0.5
final.hetRatio <- 0.25
setup.minDP <- 6
final.minDP <- 20
setup.readBoost <- F
final.readBoost <- F
readBoost.thresh <- 6
allele.fullAlign <- F


source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('Resources/ping_copy.R') # do not change
source('Resources/ping_allele.R') # do not change
source('Resources/ping_gc_align.R') # do not change
source('Resources/alleleCombine_functions.R') # do not change
setDTthreads(threads)

# Preparation -------------------------------------------------------------
# Build up a list of sample objects
sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, resultsDirectory, shortNameDelim) # no need to change


# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq') # no need to change
# Run PING2 extractor
#ext.startTime <- Sys.time()
sampleList <- extractor.run(sampleList,threads,extractedFastqDirectory,forceRun=F) # set forceRun=T if you want to force alignments
#ext.endTime <- Sys.time()
#cat('\nExtractor code timing:',ext.endTime - ext.startTime)

#copy.startTime <- Sys.time()
# PING2 gene content and copy number --------------------------------------
cat('\n\n----- Moving to PING gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=resultsDirectory,forceRun=T,onlyKFF=F) # set forceRun=T if you want to force alignments
#copy.endTime <- Sys.time()
#cat('\nCopy code timing:',copy.endTime - copy.startTime)

' 6.12 hours for 50 samples at 40 threads, 150bp, 50dp
'

sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=resultsDirectory,use.threshFile = T) # this function sets copy thresholds
sampleList <- ping_copy.load_copy_results( sampleList, resultsDirectory )

# Fix for poor 2DL2 
#sapply(sampleList, function(x) x$copyNumber[['KIR2DL2']] <- as.character(2-as.integer(x$copyNumber[['KIR2DL3']])))

## Seeing the following message when generating copy plots is normal
'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'

## use this section for notes / recording suspicious samples (has no impact on pipeline operation or results)
' 
KIR3DP1 [example]
  IND00001 [example]
'

# PING2 alignments and allele calling ----------------------------------------------
# Iter align workflow
source('Resources/genotype_alignment_functions.R') # do not change
source('Resources/alleleSetup_functions.R')

# Alignment and allele calling workflow
for(currentSample in sampleList){
  
  currentSample <- alleleSetup.gc_matched_ref_alignment( currentSample, alleleSetupDirectory, alleleSeq.list, threads)
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$ASSamPath, delIndex.list, processSharedReads = setup.readBoost, readBoost.thresh )
  file.remove(currentSample$ASSamPath)
  
  currentSample <- alleleSetup.prep_results_directory( currentSample, alignmentFileDirectory )
  
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP )
  currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory )
  
  synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
  
  currentSample <- ping_iter.run_alignments(currentSample, threads)
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = final.readBoost, readBoost.thresh )
  currentSample <- pingAllele.generate_snp_df( currentSample,uniqueSamDT,setup.knownSnpDFList,'final', final.hetRatio, final.minDP )
  
  cat('\n\n\n----- Final allele calling -----')
  for( currentLocus in names( currentSample[['snpDFPathList']][['final']] )){
    cat('\n\t',currentLocus)
    currentSample <- pingAllele.call_final_alleles(currentSample, currentLocus, knownSnpDFList[[currentLocus]]$snpDF)
  }
  
  currentSample <- pingAllele.save_call( currentSample, alleleDFPathList$iter$alleleCallPath )
}

source('Resources/alleleFinalize_functions.R')
# ----- Formatting Results Genotypes (carry directly over to PING_run) -----

cat('\n\n ----- FINALIZING GENOTYPES ----- ')
finalCallPath <- pingFinalize.format_calls( resultsDirectory )
cat('\nFinal calls written to:',finalCallPath)


