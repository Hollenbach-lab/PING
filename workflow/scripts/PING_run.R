# ---- DEPENDENCIES ----
  #' if any dependencies are missing, install with
  #install.packages("plotly",dependencies = T)
  #'
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
  rawFastqDirectory <- assert_not_null(snakemake@params[["rawFastqDirectory"]])
  fastqPattern <- assert_not_null(snakemake@params[["fastqPattern"]])
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
owd <- getwd()
setwd(workingDirectory)


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
  sampleList <- extractor.run(sampleList,threads,outDir.extFqDir$path,forceRun=T) # set forceRun=T if you want to force alignments

  # PING2 gene content and copy number --------------------------------------
  source('Resources/genotype_alignment_functions.R') # do not change
  source('Resources/alleleSetup_functions.R')
  cat('\n\n----- Moving to PING gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=outDir$path,forceRun=F,onlyKFF=F,fullAlign = copy.fullAlign, hetRatio=setup.hetRatio, minDP = setup.minDP) # set forceRun=T if you want to force alignments
out_d <- normalizePath(outDir$path)
setwd(owd)
tar(output_f,out_d , compression = 'gzip', tar="tar")

  #cat('\n\nZipping relevant results.')
  #zip(file.path(resultsDirectory,'copy_output.zip'),c(file.path(resultsDirectory,'snp_output'), file.path(resultsDirectory,'locusCountFrame.csv'),file.path(resultsDirectory,'kffCountFrame.csv'), file.path(resultsDirectory,'copyPlots')))
  #cat('\n\n----- Finished run, please find results at',paste0(resultsDirectory,'copy_output.zip'))
  # 
  # ' 6.12 hours for 50 samples at 40 threads, 150bp, 50dp
  # '
  # sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=outDir$path,use.threshFile = T) # this function sets copy thresholds
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
