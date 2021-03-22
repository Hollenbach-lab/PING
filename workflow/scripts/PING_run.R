# ---- DEPENDENCIES ----
library(data.table)
library(stringr)
library(methods)
library(pryr)
library(plotly)
library(R.utils)
library(gtools)
library(zip)


##' Checks that the input is well-formed and returns it or errors out
assert_not_null <- function(x){
  stopifnot(!is.null(x),length(x)>0)
  return(x)
}

workingDirectory <- assert_not_null(snakemake@params[["workingDirectory"]])
threads <- as.integer(snakemake@threads)
resultsDirectory <- assert_not_null(snakemake@params[["resultsDirectory"]])
shortNameDelim <- assert_not_null(snakemake@params[["shortNameDelim"]])
setup.hetRatio <- as.numeric(assert_not_null(snakemake@params[["setup_hetRatio"]]))
final.hetRatio <- as.numeric(assert_not_null(snakemake@params[["final_hetRatio"]]))
setup.minDP <- as.numeric(assert_not_null(snakemake@params[["setup_minDP"]]))
final.minDP <- as.numeric(assert_not_null(snakemake@params[["final_minDP"]]))
copy.readBoost <- (assert_not_null(snakemake@params[["copy_readBoost"]]))=="T"
setup.readBoost <- (assert_not_null(snakemake@params[["setup_readBoost"]]))=="T"
final.readBoost <- (assert_not_null(snakemake@params[["final_readBoost"]]))=="T"
readBoost.thresh <- as.numeric(assert_not_null(snakemake@params[["readBoost_thresh"]]))
allele.fullAlign <- assert_not_null(snakemake@params[["allele_fullAlign"]])=="T"
copy.fullAlign <- assert_not_null(snakemake@params[["copy_fullAlign"]])=="T"

output_f <- assert_not_null(unlist(snakemake@output))
inputf1 <- normalizePath(assert_not_null(snakemake@input[["fq1"]]),mustWork=TRUE)
inputf2 <- normalizePath(assert_not_null(snakemake@input[["fq2"]]),mustWork=TRUE)

owd <- getwd()
setwd(workingDirectory)


source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('Resources/ping_copy.R') # do not change
source('Resources/ping_allele.R') # do not change
source('Resources/ping_gc_align.R') # do not change
source('Resources/alleleCombine_functions.R') # do not change



sample <- setRefClass("sample",
                      fields=list(name='character',
                                  rawfastq1path='character',
                                  rawfastq2path='character',
                                  kirfastq1path='character',
                                  kirfastq2path='character',
                                  geneContent='list',
                                  kffHits='list',
                                  copyNumber='list',
                                  failed='logical',
                                  haploType='list',
                                  filterType='list',
                                  gzip='logical',
                                  samPath='character',
                                  bamPath='character'))


sessionInfo()
setDTthreads(threads)
                                        # Preparation -------------------------------------------------------------
outDir <- pathObj(name='output_directory',path=resultsDirectory)
outDir$dirGen()

                                        # Build up a list of sample objects


sampleList <- list(sample(name=assert_not_null(snakemake@params[["samplename"]]),
                          rawfastq1path=inputf1,
                          rawfastq2path=inputf2,
                          gzip=TRUE,
                          failed=FALSE))
names(sampleList) <- snakemake@params[["samplename"]]

                                        # PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING KIR extraction -----')
                                        # Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq') # no need to change
outDir.extFqDir <- pathObj(name='extractedFqDir',path=extractedFastqDirectory)
outDir.extFqDir$dirGen()

                                        # Run PING2 extractor

sampleList <- extractor.run(sampleList,threads,outDir.extFqDir$path,forceRun=TRUE) # set forceRun=T if you want to force alignments

                                        # PING2 gene content and copy number --------------------------------------
source('Resources/genotype_alignment_functions.R') # do not change
source('Resources/alleleSetup_functions.R')
cat('\n\n----- Moving to PING gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=outDir$path,forceRun=TRUE,onlyKFF=F,fullAlign = copy.fullAlign, hetRatio=setup.hetRatio, minDP = setup.minDP) # set forceRun=T if you want to force alignments
out_d <- normalizePath(outDir$path)
setwd(owd)
tar(output_f,out_d , compression = 'gzip', tar="tar")
unlink(out_d,recursive=TRUE)
