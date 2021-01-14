
#setwd('/home/wmarin/ping_reborn/PING/') #Set this to your own PING2 working directory
cwd <- Sys.getenv("CWD", unset='~/PING')
setwd(cwd) #Set this to your own PING working directory

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
rawFastqDirectory <- Sys.getenv("RAW_FASTQ_DIR", unset='~/PING/test_sequence/') # can be set to raw sequence or extractedFastq directory
fastqPattern <- Sys.getenv("FASTQ_PATTERN", unset='fastq') # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- Sys.getenv("THREADS", unset=4)
resultsDirectory <- Sys.getenv("OUTDIR", unset='~/3_test_sequence_results/') # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- Sys.getenv("SHORTNAMEDELIM", unset='_') # can set a delimiter to shorten sample ID's (ID will be characters before delim)
setup.hetRatio <- Sys.getenv("SETUP_HETRATIO", unset=0.25)
final.hetRatio <- Sys.getenv("FINAL_HETRATIO", unset=0.25)
setup.minDP <- Sys.getenv("SETUP_MINDP", unset=8)
final.minDP <- Sys.getenv("FINAL_MINDP", unset=20)
copy.readBoost <- Sys.getenv("COPY_READBOOST", unset=T)
setup.readBoost <- Sys.getenv("SETUP_READBOOST", unset=T)
final.readBoost <- Sys.getenv("FINAL_READBOOST", unset=F)
readBoost.thresh <- Sys.getenv("READBOOST_THRESH", unset=2)
allele.fullAlign <- Sys.getenv("ALLELE_FULLALIGN", unset=F)
copy.fullAlign <- Sys.getenv("COPY_FULLALIGN", unset=F)

# Initialization variables ------------------------------------------------
# rawFastqDirectory <- 'test_sequence/' # can be set to raw sequence or extractedFastq directory
# fastqPattern <- 'fastq' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
# threads <- 36
# resultsDirectory <- '3_test_sequence_results/' # Set the master results directory (all pipeline output will be recorded here)
# shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)
# setup.hetRatio <- 0.5
# final.hetRatio <- 0.25
# setup.minDP <- 6
# final.minDP <- 20
# copy.readBoost <- T
# setup.readBoost <- T
# final.readBoost <- F
# readBoost.thresh <- 2
# allele.fullAlign <- F
# copy.fullAlign <- F


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
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=outDir$path,forceRun=F,onlyKFF=F,fullAlign = copy.fullAlign) # set forceRun=T if you want to force alignments

cat('\n\nZipping relevant results.')
zip(file.path(resultsDirectory,'copy_output.zip'),c(file.path(resultsDirectory,'snp_output'), file.path(resultsDirectory,'locusCountFrame.csv'),file.path(resultsDirectory,'kffCountFrame.csv'), file.path(resultsDirectory,'copyPlots')))
cat('\n\n----- Finished run, please find results at',paste0(resultsDirectory,'copy_output.zip'))
# 
# ' 6.12 hours for 50 samples at 40 threads, 150bp, 50dp
# '
# sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=outDir$path,use.threshFile = F) # this function sets copy thresholds
# sampleList <- ping_copy.load_copy_results( sampleList, outDir$path )
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

# 
# if(allele.fullAlign){
#   as.list <- alleleSeq.list
# }else{
#   as.list <- compact.alleleSeq.list
# }
# 
# # Alignment and allele calling workflow
# for(currentSample in sampleList){
#   
#   currentSample <- alleleSetup.gc_matched_ref_alignment( currentSample, alleleSetupDirectory, as.list, threads)
#   
#   if( currentSample[['ASSamPath']] == 'failed' ){
#     next
#   }
#   
#   uniqueSamDT <- alleleSetup.process_samDT( currentSample$ASSamPath, delIndex.list, processSharedReads = setup.readBoost, readBoost.thresh )
#   file.remove(currentSample$ASSamPath)
# 
#   currentSample <- alleleSetup.prep_results_directory( currentSample, alignmentFileDirectory )
# 
#   currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP )
#   currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory )
# 
#   synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
# 
#   currentSample <- ping_iter.run_alignments(currentSample, threads)
#   uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = final.readBoost, readBoost.thresh )
#   currentSample <- pingAllele.generate_snp_df( currentSample,uniqueSamDT,currentSample[['iterRefDirectory']],setup.knownSnpDFList,'final', final.hetRatio, final.minDP )
# 
#   cat('\n\n\n----- Final allele calling -----')
#   for( currentLocus in names( currentSample[['snpDFPathList']][['final']][['SNP']] )){
#     cat('\n\t',currentLocus)
#     currentSample <- pingAllele.call_final_alleles(currentSample, currentLocus, knownSnpDFList[[currentLocus]]$snpDF)
#   }
# 
#   currentSample <- pingAllele.save_call( currentSample, alleleDFPathList$iter$alleleCallPath )
# }
# 
# source('Resources/alleleFinalize_functions.R')
# # ----- Formatting Results Genotypes (carry directly over to PING_run) -----
# 
# cat('\n\n ----- FINALIZING GENOTYPES ----- ')
# finalCallPath <- pingFinalize.format_calls( resultsDirectory )
# cat('\nFinal calls written to:',finalCallPath)
# 
# 
