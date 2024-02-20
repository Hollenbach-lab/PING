#library("rstudioapi")
#setwd( dirname( getActiveDocumentContext()$path ) ) # This command should set the working directory to the PING_run.R directory
#setwd() # Set working directory manually if the above command does not work

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
library(parallel)
library(tidyverse)
library(cluster)
#library(factoextra)
library(glue)
library(argparser)
library(pandoc)


# DOCKER Initialization variables ------------------------------------------------
# rawFastqDirectory <- Sys.getenv("RAW_FASTQ_DIR", unset='test_sequence/') # can be set to raw sequence or extractedFastq directory
# fastqPattern <- Sys.getenv("FASTQ_PATTERN", unset='fastq') # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
# threads <- Sys.getenv("THREADS", unset=4)
# resultsDirectory <- Sys.getenv("RESULTS_DIR", unset='test_sequence_output/') # Set the master results directory (all pipeline output will be recorded here)
# shortNameDelim <- Sys.getenv("SHORTNAME_DELIM", unset='') # can set a delimiter to shorten sample ID's (ID will be characters before delim)
# minDP <- Sys.getenv("MIN_DP", unset=10)
# hetRatio <- Sys.getenv("HET_RATIO", unset=0.25)
# setup.readBoost <- Sys.getenv("SETUP_READBOOST", unset=T)
# final.readBoost <- Sys.getenv("FINAL_READBOOST", unset=F)
# readBoost.thresh <- Sys.getenv("READBOOST_THRESH", unset=6)
# allele.fullAlign <- Sys.getenv("ALLELE_FULLALIGN", unset=F)
# copy.fullAlign <- Sys.getenv("COPY_FULLALIGN", unset=T)

# ARGPARSER setup
p <- arg_parser("Run PING")
p <- add_argument(p, "--fqDirectory", help='The path to the directory holding your fastq dataThe path to the directory holding your fastq data')
p <- add_argument(p, "--fastqPattern", help='A string that is shared across all of your fastq file names (used to find fq files and match pairs), this is usually fq or fastq', default = 'fq')
p <- add_argument(p, "--resultsDirectory", help='The path to your desired output directory')
argv <- parse_args(p)

# RSTUDIO / RSCRIPT Initialization variables ------------------------------------------------
rawFastqDirectory <- argv$fqDirectory # can be set to raw sequence or extractedFastq directory
fastqPattern <- argv$fastqPattern # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 40
resultsDirectory <- argv$resultsDirectory # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- '' # can set a delimiter to shorten sample ID's (ID will be characters before delim, ID's must be unique or else there will be an error)
setup.minDP <- 8
final.minDP <- 20

# Run mode variables ------------------------------------------------
setup.hetRatio <- 0.25
final.hetRatio <- 0.25
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

failureLog <- file.path(outDir$path,'failure.log')
# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq') # no need to change
outDir.extFqDir <- pathObj(name='extractedFqDir',path=extractedFastqDirectory)
outDir.extFqDir$dirGen()

# Run PING2 extractor
sampleList <- extractor.run(sampleList,threads,outDir.extFqDir$path,forceRun=F) # set forceRun=T if you want to force alignments

# PING2 gene content and copy number --------------------------------------
source('Resources/genotype_alignment_functions.R')
source('Resources/alleleSetup_functions.R')
cat('\n\n----- Moving to PING gene content and copy determination -----')

sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=outDir$path,forceRun=F,onlyKFF=F,fullAlign = F) # set forceRun=T if you want to force alignments
# sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=outDir$path,use.threshFile = F) # this function sets copy thresholds
sampleList <- ping_copy.load_copy_results( sampleList, outDir$path )

## Fix for poor 2DL2  
#sapply(sampleList, function(x) x$copyNumber[['KIR2DL2']] <- as.character(2-as.integer(x$copyNumber[['KIR2DL3']])))

## Seeing the following message when generating copy plots is normal
'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'


# PING2 alignments and allele calling ----------------------------------------------

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
sampleList <- ping_allele(sampleList)

# ----- Formatting Results Genotypes -----
source('Resources/alleleFinalize_functions.R')
cat('\n\n ----- FINALIZING GENOTYPES ----- ')
finalCallPath <- pingFinalize.format_calls( resultsDirectory )
cat('\nFinal calls written to:',finalCallPath)
