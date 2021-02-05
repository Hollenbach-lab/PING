
setwd('/home/wmarin/PING/') #Set this to your own PING2 working directory
#cwd <- Sys.getenv("CWD", unset='~/PING')
#setwd(cwd) #Set this to your own PING working directory

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
library(zip)

# Initialization variables ------------------------------------------------
# rawFastqDirectory <- Sys.getenv("RAW_FASTQ_DIR", unset='~/PING/test_sequence/') # can be set to raw sequence or extractedFastq directory
# fastqPattern <- Sys.getenv("FASTQ_PATTERN", unset='fastq') # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
# threads <- Sys.getenv("THREADS", unset=4)
# resultsDirectory <- Sys.getenv("RESULTS_DIR", unset='~/3_test_sequence_results/') # Set the master results directory (all pipeline output will be recorded here)
# shortNameDelim <- Sys.getenv("SHORTNAME_DELIM", unset='_') # can set a delimiter to shorten sample ID's (ID will be characters before delim)
# minDP <- Sys.getenv("MIN_DP", unset=10)
# hetRatio <- Sys.getenv("HET_RATIO", unset=0.25)
# setup.readBoost <- Sys.getenv("SETUP_READBOOST", unset=T)
# final.readBoost <- Sys.getenv("FINAL_READBOOST", unset=F)
# readBoost.thresh <- Sys.getenv("READBOOST_THRESH", unset=6)
# allele.fullAlign <- Sys.getenv("ALLELE_FULLALIGN", unset=F)
# copy.fullAlign <- Sys.getenv("COPY_FULLALIGN", unset=T)

# Initialization variables ------------------------------------------------
rawFastqDirectory <- '~/synSeq_data_run4/reads/' # can be set to raw sequence or extractedFastq directory
fastqPattern <- 'fq' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 36
resultsDirectory <- '~/4_condensed_gc_align/' # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)
setup.hetRatio <- 0.25
final.hetRatio <- 0.25
setup.minDP <- 8
final.minDP <- 20
copy.readBoost <- T
setup.readBoost <- T
final.readBoost <- F
readBoost.thresh <- 2
allele.fullAlign <- F
copy.fullAlign <- F


source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('Resources/ping_copy.R') # do not change
source('Resources/ping_allele.R') # do not change
source('Resources/ping_gc_align.R') # do not change
source('Resources/alleleCombine_functions.R') # do not change
setDTthreads(threads)

# Preparation -------------------------------------------------------------
outDir <- pathObj(name='output_directory',path=resultsDirectory)
outDir$dirGen()

# Build up a list of sample objects
sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, outDir$path, shortNameDelim) # no need to change


# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq') # no need to change
outDir.extFqDir <- pathObj(name='extractedFqDir',path=extractedFastqDirectory)
outDir.extFqDir$dirGen()

# Run PING2 extractor
sampleList <- extractor.run(sampleList,threads,outDir.extFqDir$path,forceRun=F) # set forceRun=T if you want to force alignments

# PING2 gene content and copy number --------------------------------------
source('Resources/genotype_alignment_functions.R') # do not change
source('Resources/alleleSetup_functions.R')
cat('\n\n----- Moving to PING gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=outDir$path,forceRun=T,onlyKFF=T,fullAlign = copy.fullAlign) # set forceRun=T if you want to force alignments

#cat('\n\nZipping relevant results.')
#zip(file.path(resultsDirectory,'copy_output.zip'),c(file.path(resultsDirectory,'snp_output'), file.path(resultsDirectory,'locusCountFrame.csv'),file.path(resultsDirectory,'kffCountFrame.csv'), file.path(resultsDirectory,'copyPlots')))
#cat('\n\n----- Finished run, please find results at',paste0(resultsDirectory,'copy_output.zip'))
# 
# ' 6.12 hours for 50 samples at 40 threads, 150bp, 50dp
# '
sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=outDir$path,use.threshFile = T) # this function sets copy thresholds
sampleList <- ping_copy.load_copy_results( sampleList, outDir$path )
# 
# # Fix for poor 2DL2 
# #sapply(sampleList, function(x) x$copyNumber[['KIR2DL2']] <- as.character(2-as.integer(x$copyNumber[['KIR2DL3']])))
# 
# ## Seeing the following message when generating copy plots is normal
# 'A line object has been specified, but lines is not in the mode
# Adding lines to the mode...'
# 
# ## use this section for notes / recording suspicious samples (has no impact on pipeline operation or results)
# ' 
# KIR3DP1 [example]
#   IND00001 [example]
# '
# 
# PING2 alignments and allele calling ----------------------------------------------
# Iter align workflow


if(allele.fullAlign){
  as.list <- alleleSeq.list
}else{
  as.list <- compact.alleleSeq.list
}


probelistFile='probelist_2021_01_24.csv'
gcResourceDirectory <- normalizePath('Resources/gc_resources', mustWork = T)
cat('\n\nReading in the KFF probelist file: ', file.path(gcResourceDirectory, probelistFile))
probeDF <- read.csv(file.path(gcResourceDirectory, probelistFile), stringsAsFactors = F, check.names = F)
row.names(probeDF) <- probeDF$Name

# Alignment and allele calling workflow
for( currentSample in sampleList[1:length(sampleList)] ){
  
  currentSample <- alleleSetup.gc_matched_ref_alignment( currentSample, alleleSetupDirectory, as.list, threads)
  
  if( currentSample[['ASSamPath']] == 'failed' ){
    next
  }
  
  currentSample[['setCallList']] <- list()
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$ASSamPath, delIndex.list, processSharedReads = setup.readBoost, readBoost.thresh )
  file.remove(currentSample$ASSamPath)
  currentSample <- alleleSetup.prep_results_directory( currentSample, alignmentFileDirectory )
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP, includeAmb = T )
  currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory, setup.knownSnpDFList, alleleSetupRef.df, addFullyDefined = T)
  
  synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
  currentSample <- ping_iter.run_alignments(currentSample, threads, all.align=T)
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = T, 2 )
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP, includeAmb=F, ambScore=F, allPosScore = F, homScoreBuffer = 1, addDiversity = F)
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, final.hetRatio, final.minDP, includeAmb=T, ambScore=T, allPosScore = F, onlyExonScore = T, homScoreBuffer = 1, addDiversity=T,combineTypings=T, skipSnpGen=T)
  currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory, setup.knownSnpDFList, alleleSetupRef.df, addFullyDefined = F)
  
  synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
  currentSample <- ping_iter.run_alignments(currentSample, threads, all.align=T)
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = T, 2 )
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP, includeAmb=F, ambScore=F, allPosScore = F, homScoreBuffer = 1, addDiversity = F, skipSet = F)
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, final.hetRatio, final.minDP, includeAmb=T, ambScore=T, allPosScore = F, onlyExonScore = T, homScoreBuffer = 1, addDiversity=T,combineTypings=T,skipSnpGen=T, skipSet=F)
  currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory, setup.knownSnpDFList, alleleSetupRef.df, addFullyDefined = T)
  
  synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
  currentSample <- ping_iter.run_alignments(currentSample, threads, all.align=F)
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = F, 2 )
  currentSample <- pingAllele.generate_snp_df( currentSample,uniqueSamDT,currentSample[['iterRefDirectory']],setup.knownSnpDFList,'final', final.hetRatio, final.minDP )

  cat('\n\n\n----- Final allele calling -----')
  for( currentLocus in names( currentSample[['snpDFPathList']][['final']][['SNP']] )){
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
