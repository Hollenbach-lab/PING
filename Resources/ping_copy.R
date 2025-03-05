#.libPaths("/home/wmarin/R/x86_64-redhat-linux-gnu-library/3.4")
library(data.table)
library(stringr)
library(methods)
library(pryr)
library(plotly)
#library(randomForest)

#### 3/11/20 gc to-do notes
# 1. Create results subdirectory for bam files [x]
# 2. Place copy lines on graphs for manual copy setting [x]
# 3. Create system for marking samples for follow-up [x]

########## Development INPUT variables
#setwd('/home/LAB_PROJECTS/PING2_PAPER/PING2')
#sequenceDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/sequence_data/'
#fastqPattern <- 'fastq'
#threads <- 12
#resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/'
#KIR3DL3MinReadThreshold <- 100
#maxReadThreshold <- 50000
#probelistFile <- 'probelist_2018_08_02.csv'
#predictCopy <- T
###########
source('Resources/gc_functions.R')

ping_copy.version <- '2.1'
cat('\nPING_copy version:',ping_copy.version)

ping_copy.graph <- function(sampleList=list(),
                    threads=4,
                    resultsDirectory,
                    KIR3DL3MinReadThreshold=500,
                    maxReadThreshold=50000,
                    probelistFile='probelist_2021_01_24.csv',
                    onlyKFF=F,
                    forceRun=F,
                    fullAlign=F){
  predictCopy=T
  kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                    'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                    'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')
  
  cat('Current working directory:', getwd(),'\n')
  cat('Max read count set at:',maxReadThreshold,'\n')
  #cat('Sequence directory set to:',sequenceDirectory,'\n')
  #cat('FASTQ pattern set to:',fastqPattern,'\n')
  #cat('Results directory set at:',resultsDirectory,'\n')
  
  
  setDTthreads(threads)
  
  ### Set up directory paths, make sure they exist
  resultsDirectory <- normalizePath(resultsDirectory)
  ## Create the results directory if it does not exist
  if(!file.exists(resultsDirectory)){
    dir.create(resultsDirectory, recursive = T)
  }
  
  ## Create directory for holding copy plots
  plotDirectory <- file.path(resultsDirectory,'copyPlots')
  if(!file.exists(plotDirectory)){
    dir.create(plotDirectory)
  }
  
  ## Create directory for holding bam files
  bamDirectory <- file.path(resultsDirectory,'gc_bam_files')
  if(!file.exists(bamDirectory)){
    dir.create(bamDirectory)
  }
  
  ## Create directory for holding the SNP output files
  snpDir <- pathObj(name='snp_output',path=file.path(resultsDirectory,'snp_output'))
  snpDir$dirGen()
  
  #sequenceDirectory <- normalizePath(sequenceDirectory, mustWork=T)
  gcResourceDirectory <- normalizePath('Resources/gc_resources', mustWork = T)
  ### /Set up
  
  if( fullAlign ){
    ### Read in reference files
    kirReferenceFasta <- normalizePath(file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled.fasta'), mustWork=T)
    kirReferenceIndex <- file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled')
  }else{
    kirReferenceFasta <- normalizePath(file.path(gcResourceDirectory,'filled_kir_reference','KIR_compact_filled.fasta'), mustWork=T)
    kirReferenceIndex <- file.path(gcResourceDirectory,'filled_kir_reference','KIR_compact_filled')
  }
  ### /Read in
  
  ### Initialize lists of kir alleles at different resolutions
  kirAlleleList <- read.kir_allele_list_from_reference_fasta(kirReferenceFasta)
  kirAlleleListRes3 <- unique(unlist(lapply(kirAlleleList, kir.allele_resolution, 3)))
  kirAlleleListRes5 <- unique(unlist(lapply(kirAlleleList, kir.allele_resolution, 5)))
  ### /Initialize  
  
  ### Read in the probelist CSV file as a dataframe
  cat('\n\nReading in the KFF probelist file: ', file.path(gcResourceDirectory, probelistFile))
  probeDF <- read.csv(file.path(gcResourceDirectory, probelistFile), stringsAsFactors = F, check.names = F)
  row.names(probeDF) <- probeDF$Name
  ### /Read in  
  
  ### Pull out all the probe names with '>'
  kffPresenceProbeNameList <- grep('>', probeDF$Name, fixed=T, value=T)
  
  ### Split the pulled out probe names by '>', then grab all the unique locus names
  kffLociList <- unique(tstrsplit(kffPresenceProbeNameList,'>',fixed=T)[[2]])
  
  ### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
  #bowtie2Build <- system2('which', c('bowtie2-build'), stdout=T, stderr=T)
  #check.system2_output(bowtie2Build, 'bowtie2-build not found')
  
  ## Creqte a bowtie2 index for the kir_reference.fasta file <- only needed when building a new reference index
  #createIndex <- system2(bowtie2Build, c(fullKirReferenceFasta, fullKirReferenceIndex))
  #check.system2_output(createIndex, 'bowtie2 index building failed')
  
  ### Building a list of sample objects from files in sampleDirectory that match fastqPattern
  #sampleList <- build.paired_sample_objects(sequenceDirectory,fastqPattern,resultsDirectory)
  
  ### Check to make sure bowtie2is accessible
  #bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
  #check.system2_output(bowtie2, 'bowtie2 not found')
  
  ### Define paths for output files
  kffCountDFFile <- file.path(resultsDirectory, 'kffCountFrame.csv')
  kffNormDFFile <- file.path(resultsDirectory, 'kffNormFrame.csv')
  kffPresenceDFFile <- file.path(resultsDirectory, 'kffPresenceFrame.csv')
  locusCountDFFile <- file.path(resultsDirectory, 'locusCountFrame.csv')
  #alleleCountDFFile <- file.path(resultsDirectory, 'alleleCountFrame.csv')
  ### /Define
  
  ## Set up dataframe to store copy thresholds
  threshPath <- file.path(resultsDirectory, 'manualCopyThresholds.csv')
  
  if(!file.exists(threshPath)){
    ## Set up dataframe to store copy thresholds
    thresholdCols <- c('0-1','1-2','2-3','3-4','4-5','5-6')
    thresholdDF <- data.frame(matrix(NA,nrow=length(kirLocusList),ncol=length(thresholdCols)),stringsAsFactors = F)
    rownames(thresholdDF) <- kirLocusList
    colnames(thresholdDF) <- thresholdCols
  }else{
    ## Load in copy threshold file
    threshPath <- normalizePath(threshPath, mustWork=T)
    cat(paste0('\nFound manualCopyThresholds.csv in ',threshPath,'. Loading these results.'))
    thresholdDF <- read.csv(threshPath, stringsAsFactors=F, check.names=F,row.names=1)
  }
  
  if( grepl('test_sequence_output',resultsDirectory,fixed=T) & grepl('IND00001', sampleList[[1]]$name, fixed=T) ){
    cat('\nUsing example threshold data from Resources/gc_resources/manualCopyThresholds_example.csv')
    thresholdDF <- read.csv('Resources/gc_resources/manualCopyThresholds_example.csv',stringsAsFactors = F, check.names=F,row.names=1)
  }
  
  # write.csv(thresholdDF, file = threshPath)
  ## /Threshold df setup
  
  previousIDVect <- c()
  ### Check to see if the output files exist in the results directory
  if(!forceRun & all(file.exists(c(kffCountDFFile, kffNormDFFile, kffPresenceDFFile, locusCountDFFile)))){
    
    ## Load in found count file
    cat(paste0('\n\nFound kffCountFrame.csv in ', resultsDirectory, '. Loading these results.'))
    kffCountDF <- read.csv(kffCountDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
    
    ## Load in found count file
    cat(paste0('\nFound kfNormFrame.csv in ', resultsDirectory, '. Loading these results.'))
    kffNormDF <- read.csv(kffNormDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
    
    ## Load in found count file
    cat(paste0('\nFound kffPresenceFrame.csv in ', resultsDirectory, '. Loading these results.'))
    kffPresenceDF <- read.csv(kffPresenceDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
    
    ## Load in found count file
    cat(paste0('\nFound locusCountFrame.csv in ', resultsDirectory, '. Loading these results.'))
    locusCountDF <- read.csv(locusCountDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
    
    ## Load in found count file
    #cat(paste0('\nFound alleleCountFrame.csv in ', resultsDirectory, '. Loading these results.'))  
    #alleleCountDF <- read.csv(alleleCountDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
    
    sampleStart <- sum(apply(locusCountDF, 1, sum) > 0)
    previousIDVect <- intersect(rownames(locusCountDF)[apply(locusCountDF,1,sum) > 0],names(sampleList))
  }else{
    
    ## Initialize a dataframe for counting KFF probe matches
    kffCountDF <- data.frame(matrix(0, length(sampleList), length(probeDF$Name)),row.names=names(sampleList),check.names=F,stringsAsFactors=F)
    colnames(kffCountDF) <- probeDF$Name
    
    ## Initialize a dataframe for storing the normalized probe matches
    kffNormDF <- data.frame(matrix(0, length(sampleList), length(kffLociList)),row.names=names(sampleList),check.names=F,stringsAsFactors=F)
    colnames(kffNormDF) <- kffLociList
    
    ## Initialize a dataframe for storing kff determined locus presence/absence values
    kffPresenceDF <- data.frame(matrix(0, length(sampleList), length(kffLociList)),row.names=names(sampleList),check.names=F,stringsAsFactors=F)
    colnames(kffPresenceDF) <- kffLociList
    
    ## Initialize a dataframe for counting how many reads uniquely align to a locus
    locusCountDF <- data.frame(matrix(0, length(sampleList), length(kirLocusList)),row.names=names(sampleList),check.names=F,stringsAsFactors=F)
    colnames(locusCountDF) <- kirLocusList
    
    ## Initialize a dataframe for counting how many reads uniquely align to an allele (at protein coding resolution)
    #alleleCountDF <- data.frame(matrix(0, length(sampleList), length(kirAlleleListRes3)), row.names=names(sampleList),check.names=F,stringsAsFactors=F)
    #colnames(alleleCountDF) <- kirAlleleListRes3
    
    sampleStart <- 1
  }
  ### /Check
  
  ## Load results into sample objects for samples that have already been run
  if(length(previousIDVect) > 0){
    for(sampleID in previousIDVect){
      ## Transfer the KFF hit data
      currentSampleProbeList <- as.list(kffCountDF[sampleID,])
      
      currentSampleGeneContent <- as.list(kffPresenceDF[sampleID,])
      
      ## Consolidate probe ID's
      probeIDVect <- unique(tstrsplit(names(currentSampleProbeList),'_rc')[[1]])
      
      ## Fill out the probe ID list
      currentSampleProbeIDList <- list()
      for(probeID in probeIDVect){
        currentProbeVect <- grep(probeID,names(currentSampleProbeList),value=T,fixed=T)
        currentSampleProbeIDList[[probeID]] <- sum(unlist(currentSampleProbeList[currentProbeVect]))
      }
    
      sampleList[[sampleID]]$kffHits <- currentSampleProbeIDList
      sampleList[[sampleID]]$geneContent <- currentSampleGeneContent
    }
  }
  
  ## Run all samples through bowtie2 gc alignment
  for(currentSample in sampleList[!(names(sampleList) %in% previousIDVect)]){
    
    if(currentSample$failed){next}
    
    cat('\n\nProcessing', currentSample$name)
    cat('\n------------------------------------')
    
    if( file.size(currentSample$kirfastq1path) < 10000 ){
      message('\n-- Sample failure --')
      message(paste0(currentSample$name,': extracted fastq size smaller than 10KB'))
      currentSample$failed <- T
      currentSample[['failureMessage']] <- paste0(currentSample$name,': Extracted fastq smaller than 10KB')
      cat(currentSample[['failureMessage']], file = failureLog, sep = "\n", append = TRUE)
      next
    }
    
    cat('\n\nCounting KFF primer matches.')
    
    ## Count KFF probe matches for the currentSample
    kffCountList <- run.count_kff_probes(currentSample, probeDF, maxReadThreshold)
    kffCountDF[currentSample$name,names(kffCountList)] <- kffCountList 
    
    ## Transfer the KFF hit data
    currentSampleProbeList <- as.list(kffCountDF[currentSample$name,])
    
    ## Consolidate probe ID's
    probeIDVect <- unique(tstrsplit(names(currentSampleProbeList),'_rc')[[1]])
    
    ## Fill out the probe ID list
    currentSampleProbeIDList <- list()
    for(probeID in probeIDVect){
      currentProbeVect <- grep(probeID,names(currentSampleProbeList),value=T,fixed=T)
      currentSampleProbeIDList[[probeID]] <- sum(unlist(currentSampleProbeList[currentProbeVect]))
    }
    
    currentSample$kffHits <- currentSampleProbeIDList
    
    ## Write the results to a csv file
    write.csv(kffCountDF, file = kffCountDFFile)
    
    cat('\nNormalizing KFF primer matches.')
    
    ## Normalize the KFF probe matches by KIR3DL3
    kffNormList <- run.reduce_and_normalize_kff_probes(kffCountList, kffLociList)
    kffNormDF[currentSample$name,names(kffNormList)] <- kffNormList
    
    ## Write the results to a csv file
    write.csv(kffNormDF, file = kffNormDFFile)
    
    cat('\nDetermining KIR locus presence/absence')
    
    ## Determine locus presence/absence
    kffPresenceList <- run.kff_determine_presence_from_norm_values(kffNormList, kffThreshold=0.2)
    
    ## Mutate kff presence results to the same format as the copy results
    if(any(sapply(kffPresenceList, is.na))){
      kffPresenceList[names(kffPresenceList)] <- 'failed'
    }else{
      kffPresenceList[names(kffPresenceList)] <- as.character(kffPresenceList)
    }
    
    kffPresenceDF[currentSample$name,names(kffPresenceList)] <- kffPresenceList
    currentSample$geneContent <- kffPresenceList
    
    ## Write the results to a csv file
    write.csv(kffPresenceDF, file = kffPresenceDFFile)
    
    cat('\n\nFinished with presence/absence determination, moving to copy number determination.')
    
    if(onlyKFF){
      cat('\n\nonlyKFF set to TRUE, skipping copy number determination.')
      next
    }
    
    ## Fill in the path to the alignment file (it may or may not be present)
    currentSample$samPath <- file.path(bamDirectory,paste0(currentSample$name,'.sam'))
    
    ## Fill in the path to the alignment file (it may or may not be present)
    currentSample$bamPath <- file.path(bamDirectory,paste0(currentSample$name,'.bam'))
    
    ## If the alignment file does not exist, then run bowtie2 alignment, otherwise continue
    if(!file.exists(currentSample$bamPath)){
      cat('\n\nCurrent used memory: ', mem_used())
      cat('\n\nPerforming bowtie2 alignment for this sample.')
      currentSample <- run.bowtie2_gc_alignment(bowtie2, kirReferenceIndex, threads, currentSample, bamDirectory)
      currentSample <- samtools.sam_to_bam(samtools, currentSample, bamDirectory, threads)
    }else{
      cat('\n\nFound a previous alignment file for this sample, skipping bowtie2 alignment.')
      currentSample <- samtools.bam_to_sam(samtools, currentSample, bamDirectory, threads)
    }
    
    samDT <- alleleSetup.process_samDT( currentSample$samPath, delIndex.list,processSharedReads = copy.readBoost, readBoost.thresh = 2 )
    
    ## Add the counts to the appropriate count dataframe
    locusCount.tab <- table( samDT$locus )
    locusCountDF[currentSample$name,names(locusCount.tab)] = as.integer( locusCount.tab )
    locusCountDF[is.na(locusCountDF)] <- 0
    
    currentSample <- pingAllele.generate_snp_df( currentSample, samDT, snpDir$path, setup.knownSnpDFList,'copy', setup.hetRatio, setup.minDP)
    
    # 
    # ## Count how many header lines there in the SAM file so they can be skipped during read in
    # headerLineCountInt <- samfile.count_header_lines(currentSample)
    # 
    # ## Read in the SAM file to analyze where the reads are aligning
    # samTable <- read.bowtie2_sam_nohd(currentSample, rows_to_skip=headerLineCountInt)
    # 
    file.remove(currentSample$samPath) ## Remove the SAM file to save space
    
    # cat('\nCounting reads that align uniquely to a locus or allele ')
    # 
    # ## Count how many reads align uniquely to a locus or allele
    # countList <- run.count_kir_read_matches(currentSample, samTable, maxReadThreshold, kirLocusList, kirAlleleListRes3)
    # 
    # ## Add the counts to the appropriate count dataframe
    # locusCountDF[currentSample$name,names(countList$locusMatches)] = countList$locusMatches
    # #alleleCountDF[currentSample$name,names(countList$alleleMatches)] = countList$alleleMatches
    
    ## Write the results to a csv file
    write.csv(locusCountDF, file = locusCountDFFile)
    #write.csv(alleleCountDF, file = file.path(resultsDirectory, 'alleleCountFrame.csv'))
    
    rm(samDT)
  }
  
  if(onlyKFF){
    return(sampleList)
  }
  
  cat('\n\n----- Finished with alignment! -----')
  cat('\n\nMoving on to copy number graphing.')
  
  ## Read in the csv results
  locusCountDF <- read.csv(locusCountDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
  kffPresenceDF <- read.csv(kffPresenceDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
  
  ## Initialize a copy number frame
  copyNumberDF <- data.frame(locusCountDF)
  copyNumberDF[,] <- 0
  
  ## Initialize a list of sample names
  sampleNameList <- names(sampleList)
  
  ## Only analyze samples that have at least 'KIR3DL3MinReadThreshold' number of unique KIR3DL3 reads
  goodRows <- rownames(locusCountDF[sampleNameList,])[apply(locusCountDF[sampleNameList,], 1, function(x) x['KIR3DL3']>=KIR3DL3MinReadThreshold)]
  
  ## Keep track of what samples are being discarded
  badRows <- rownames(locusCountDF[sampleNameList,])[apply(locusCountDF[sampleNameList,], 1, function(x) x['KIR3DL3']<KIR3DL3MinReadThreshold)]
  cat('\nSkipping', length(badRows), 'samples that had fewer than',KIR3DL3MinReadThreshold,'KIR3DL3 reads.')
  
  ## Mark bad samples as failed
  if(length(badRows) > 0){
    for(sampleID in badRows){
      if(!sampleList[[sampleID]]$failed){
        sampleList[[sampleID]]$failed <- T
        sampleList[[sampleID]][['failureMessage']] <- paste0(sampleID,': Fewer than ',KIR3DL3MinReadThreshold,' KIR3DL3 reads.')
        cat(sampleList[[sampleID]][['failureMessage']], file = failureLog, sep = "\n", append = TRUE)
        message('\n-- Sample failure --')
        message(paste0(sampleID,': Fewer than ',KIR3DL3MinReadThreshold,' KIR3DL3 reads.'))
      }
    }
  }
  
  if( length(goodRows) == 0){
    stop('No samples passed copy determination quality control.')
  }else if( length(goodRows) == 1 ){
    locusRatioDF <- as.data.frame( locusCountDF[goodRows,] / locusCountDF[goodRows,'KIR3DL3'] )
  }else{
    ## Subset the count dataframe by the samples that were determined to be good, then normalize each locus unique read count by KIR3DL3
    locusRatioDF <- apply(locusCountDF[goodRows,], 2, function(x) x / locusCountDF[goodRows,'KIR3DL3'])
    locusRatioDF <- as.data.frame(locusRatioDF)
  }
  
  ## Write the locus ratio results to a csv file
  write.csv(locusRatioDF, file = file.path(resultsDirectory, 'locusRatioFrame.csv'))
  
  if(predictCopy){
    ## Use the random forest models to predict copy number
    cat('\nPredicting copy number... ')
    # copyNumberDF <- run.predict_copy(locusRatioDF, locusCountDF, copyNumberDF, goodRows, resultsDirectory, rfAllPathList)
    # thresholdDF <- run.get_threshold(locusRatioDF, locusCountDF, copyNumberDF, goodRows, resultsDirectory, rfAllPathList)  
    cnThresDF <- run.predict_copy_get_threshold(locusRatioDF, locusCountDF, copyNumberDF, goodRows, resultsDirectory, rfAllPathList)
    copyNumberDF <- cnThresDF$dataframe1
    thresholdDF <- cnThresDF$dataframe2
    
    ## Write the results to a csv file
    cat('\nFinished with copy predictions.')
    write.csv(copyNumberDF, file = file.path(resultsDirectory, 'predictedCopyNumberFrame.csv'))

    ## Generate ratio graphs and color according to predicted copy number
    cat('\nGenerating predicted copy number graphs...')
    # run.generate_predicted_copy_number_graphs(locusRatioDF, kffPresenceDF, kirLocusList, plotDirectory, locusCountDF, copyNumberDF)
    run.generate_copy_number_graphs(locusRatioDF, kffPresenceDF, kirLocusList, plotDirectory, locusCountDF, thresholdDF)
  }else{
    ## Generate ratio graphs and color according to kff presence/absence
    cat('\nGenerating copy number graphs... ')
    run.generate_copy_number_graphs(locusRatioDF, kffPresenceDF, kirLocusList, plotDirectory, locusCountDF, thresholdDF)
  }

cat('\n\nALL FINISHED!!')
return(sampleList)
}

ping_copy.manual_threshold <- function(sampleList=list(),resultsDirectory='',use.threshFile=F){
  kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                    'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                    'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')
  
  cat('Current working directory: ', getwd(),'\n')
  
  ### Set up directory paths, make sure they exist
  resultsDirectory <- normalizePath(resultsDirectory, mustWork=T)
  cat('Running manual thresholding on:',resultsDirectory)
  
  plotDirectory <- normalizePath(file.path(resultsDirectory,'copyPlots'))
  
  ## Load in KFF presence file
  kffPresenceDFFile <- normalizePath(file.path(resultsDirectory, 'kffPresenceFrame.csv'), mustWork=T)
  cat(paste0('\n\nFound kffPresenceFrame.csv in ', resultsDirectory, '. Loading these results.'))
  kffPresenceDF <- read.csv(kffPresenceDFFile, stringsAsFactors = F, check.names = F, row.names = 1)
  
  ## Load in locus count file
  locusCountDFFile <- normalizePath(file.path(resultsDirectory, 'locusCountFrame.csv'), mustWork=T)
  cat(paste0('\nFound locusCountFrame.csv in ',resultsDirectory,'. Loading these results.'))
  locusCountDF <- read.csv(locusCountDFFile, stringsAsFactors=F, check.names=F,row.names=1)
  
  ## Load in locus ratio file
  locusRatioDFFile <- normalizePath(file.path(resultsDirectory, 'locusRatioFrame.csv'), mustWork=T)
  cat(paste0('\nFound locusRatioFrame.csv in ',resultsDirectory,'. Loading these results.'))
  locusRatioDF <- read.csv(locusRatioDFFile, stringsAsFactors=F, check.names=F,row.names=1)
  
  threshPath <- file.path(resultsDirectory, 'manualCopyThresholds.csv')
  
  threshFile.bool <- file.exists(threshPath)
  
  if(!threshFile.bool){
    ## Set up dataframe to store copy thresholds
    thresholdCols <- c('0-1','1-2','2-3','3-4','4-5','5-6')
    thresholdDF <- data.frame(matrix(NA,nrow=length(kirLocusList),ncol=length(thresholdCols)),stringsAsFactors = F)
    rownames(thresholdDF) <- kirLocusList
    colnames(thresholdDF) <- thresholdCols
  }else{
    ## Load in copy threshold file
    threshPath <- normalizePath(threshPath, mustWork=T)
    cat(paste0('\nFound manualCopyThresholds.csv in ',threshPath,'. Loading these results.'))
    thresholdDF <- read.csv(threshPath, stringsAsFactors=F, check.names=F,row.names=1)
    thresholdNA.bool <- all( apply(thresholdDF, 1, function(x) all( is.na(x)) ) )
  }
  
  
  ## Generate ratio graphs and color according to kff presence/absence
  cat('\nGenerating copy number graphs... ')
  run.generate_copy_number_graphs(locusRatioDF, kffPresenceDF, kirLocusList, plotDirectory, locusCountDF, thresholdDF)
  
  ## Initialize a copy number frame
  copyNumberDF <- data.frame(locusCountDF)
  copyNumberDF[,] <- 0
  
  
  if( (use.threshFile & !threshFile.bool) | (use.threshFile & thresholdNA.bool) ){
    
    cat("\nuse.threshFile set to TRUE but manualCopyThresholds.csv is either not found or not filled out. Terminating run.")
    stop()
    
  }else if( use.threshFile==T & threshFile.bool ){
    
    cat("\nuse.threshFile set to TRUE, loading results from manualCopyThresholds.csv")
    copyResultsList <- copy.set_copy_from_threshFile(copyNumberDF, locusRatioDF, locusCountDF, thresholdDF)
    
  }else{
  
    copyResultsList <- run.set_copy(kirLocusList=kirLocusList,
                                 copyNumberDF=copyNumberDF,
                                 locusRatioDF=locusRatioDF,
                                 locusCountDF=locusCountDF,
                                 thresholdDF=thresholdDF)
    
  }
  
  copyNumberDF <- copyResultsList$copyDF
  thresholdDF <- copyResultsList$threshDF
  
  ## Save copy results to sample objects
  for(sampleID in names(sampleList)){
    sampleList[[sampleID]]$copyNumber <- as.list(copyNumberDF[sampleID,])
  }
  
  ## Generate ratio graphs and color according to kff presence/absence
  cat('\nAdding thresholds to copy number graphs... ')
  run.generate_copy_number_graphs(locusRatioDF, kffPresenceDF, kirLocusList, plotDirectory, locusCountDF, thresholdDF)
  
  ## Write the copy number results to a csv file
  manualCopyPath <- file.path(resultsDirectory, 'manualCopyNumberFrame.csv')
  write.csv(copyNumberDF, file = manualCopyPath)
  
  write.csv(thresholdDF, file = threshPath)
  
  cat('\nFinished with copy number setting. Results are written to',manualCopyPath,'. Thresholds recorded at',threshPath)
  return(sampleList)
}

ping_copy.load_copy_results <- function( sampleList, resultsDirectory ){
  
  cat('\n\n ----- Loading up results from ping_copy -----')
  kffCountDFFile <- file.path(resultsDirectory, 'kffCountFrame.csv')
  # copyDFFile <- file.path(resultsDirectory, 'manualCopyNumberFrame.csv')
  copyDFFile <- file.path(resultsDirectory, 'predictedCopyNumberFrame.csv')
  
  if( all(file.exists(c(kffCountDFFile, copyDFFile))) ){
    
    cat(paste0('\n\nFound kffCountFrame.csv in ', resultsDirectory, '. Loading these results.'))
    kffCountDF <- read.csv(kffCountDFFile, stringsAsFactors = F, check.names = F, row.names = 1 )
    
    cat(paste0('\n\nFound manualCopyNumberFrame.csv in ', resultsDirectory, '. Loading these results.'))
    copyDF <- read.csv(copyDFFile , stringsAsFactors = F, check.names = F, row.names = 1)
    
  }else{
    cat('\n',kffCountDFFile,'\n',copyDFFile,'\n')
    stop('Copy results files absent in:',resultsDirectory)
  }
  
  for( sampleID in names(sampleList) ){
    currentSampleProbeList <- as.list( kffCountDF[sampleID,] )
    currentSampleCopyList <- as.list( copyDF[sampleID,] )
    
    ## Consolidate probe ID's
    probeIDVect <- unique(tstrsplit(names(currentSampleProbeList),'_rc')[[1]])
    
    ## Fill out the probe ID list
    currentSampleProbeIDList <- list()
    for(probeID in probeIDVect){
      currentProbeVect <- grep(probeID,names(currentSampleProbeList),value=T,fixed=T)
      currentSampleProbeIDList[[probeID]] <- sum(unlist(currentSampleProbeList[currentProbeVect]))
    }
    
    sampleList[[sampleID]]$kffHits <- currentSampleProbeIDList
    sampleList[[sampleID]]$copyNumber <- currentSampleCopyList
  }
  cat('\n\nFinished.')
  return(sampleList)
}
