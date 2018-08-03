library(methods)
library(reshape2)

## This function checks to make sure the output of system2 is valid
check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}

## This function finds samples in sampleDirectory and turns them into objects for downstream use
build.paired_sample_objects <- function(sample_directory, fastq_pattern='fastq'){
  #####
  ## This function takes in a directory and a file name pattern and attempts to pair fastq files
  ## Returns a list of sample objects that contain the paired fastq file names
  #####
  
  cat("\nAttempting automatic fastq pairing in", sample_directory, "using", fastq_pattern)
  
  ## Find all the files in sampleDirectory that match fastqPattern
  unpairedFastqList <- list.files(path=sample_directory, pattern=fastq_pattern)
  
  ## To pair reads, we will split the file names by fastqPattern, then continuously chop a
  ## character off the end of each name until the number of unique names is exactly half of the total names
  
  ## Setting up an initial fastq list that splits the files names by fastqPattern
  strList <- sapply(unpairedFastqList, function(x) str_split(x, fastq_pattern)[[1]][1])
  
  ## Setting the maximum number of times to chop off the last character to the length of the shortest name
  maxChop <- min(sapply(strList, nchar))
  
  ## Iterate from 0 to maxChop, cutting off i characters from the file names each time
  for(i in 0:maxChop){
    
    ## In each iteration, the file names are reset. There is no particular reason I implemented it this way
    subStrList <- strList
    
    ## Cut off i characters from the end of each fastq file name
    subStrList <- sapply(subStrList, function(x) substr(x, 1, nchar(x)-i))
    
    ## After cutting, determine the unique names
    uniqueFastqList <- unique(subStrList)
    
    ## If the number of unique names is exactly half of the total names, then cutting should be finished
    ## and it is time to move on to matching!
    if(length(uniqueFastqList) == (length(subStrList)/2)){
      break
    }
    
    ## Pairing failed if i reaches maxChop. This will raise a fatal error.
    if(i == maxChop){
      stop("Was not able to pair fastq file names, please check that fastqPattern and sampleDirectory are set correctly.")
    }
  }
  
  ## Initialize a list for storing the paired fastq file names
  pairedFastqList <- list()
  
  ## Iterate through the unique fastq names to make pairings
  for(fastqName in uniqueFastqList){
    
    ## Pull out the matches for fastqName in the subStrList
    fastqMatches <- subStrList[fastqName == subStrList]
    
    ## Determine how many matches there are for fastqName in the subStrList
    matchCount <- length(fastqMatches)
    
    ## Stop the program if more or less than 2 matches are found for fastqName
    if(matchCount != 2){
      cat('\n',names(fastqMatches))
      stop('Auto fastq matching failed due to an improper number of matches for ',fastqName)
    }
    
    ## Save the file names in the pairedFastqList under the unique name
    pairedFastqList[[fastqName]] <- names(fastqMatches)
  }
  
  cat("\nFound", length(uniqueFastqList), "samples in", sample_directory)
  
  ## Creating the sample object class
  sample <- setRefClass("sample",
                        fields=list(name='character',
                                    fastq1path='character',
                                    fastq2path='character',
                                    gzip='logical',
                                    gcSamPath='character'))
  
  ## Initializing a sample object list. This will be returned
  output.sampleList <- list()
  
  for(i in 1:length(pairedFastqList)){
    
    ## Pulling the current working element out of the list
    pairedFastq <- pairedFastqList[i]
    
    ## Creating an absolute path to the first fastq file
    fastq1path <- normalizePath(file.path(sample_directory, pairedFastq[[1]][1]), mustWork=T)
    
    ## Creating a absolute path to the second fastq file
    fastq2path <- normalizePath(file.path(sample_directory, pairedFastq[[1]][2]), mustWork=T)
    
    ## Checking if the first fastq file is gzipped
    gzip <- substr(fastq1path, nchar(fastq1path)-2, nchar(fastq1path)) == '.gz'
    
    ## Building a sample object and adding it to sampleList
    output.sampleList[[names(pairedFastq)]] <- sample(name=names(pairedFastq),fastq1path=fastq1path,fastq2path=fastq2path,gzip=gzip)
  }
  
  cat("\nAll samples were successfully paired")
  return(output.sampleList)
}

## This function runs a bowtie2 alignment looking for exact matches to all KIR references in subsetKirReference
run.bowtie2_gc_alignment <- function(bowtie2_command, reference_index, threads, current_sample){
  ## Intitialize an output path for the SAM file
  current_sample$gcSamPath <- file.path(resultsDirectory,paste0(current_sample$name,'.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-5 0', '-3 6', '-N 0', '--end-to-end', paste0('-p ',threads), '--score-min "C,-2,0"',
                      paste0('-1 ',current_sample$fastq1path),
                      paste0('-2 ',current_sample$fastq2path),
                      '--no-unal','--no-hd','-a','--np 0', '--mp 2,2', '--rdg 1,1', '--rfg 1,1',
                      paste0('-S ',current_sample$gcSamPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  current_sample$gcSamPath <- normalizePath(current_sample$gcSamPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',current_sample$name,'to',reference_index)
  
  return(output.sampleAlign)
}

## This function reads in a SAM file with no header to a data.table
read.bowtie2_sam_nohd <- function(sam_path){
  
  ## Make sure the SAM file can be read in
  sam_path <- normalizePath(sam_path, mustWork=T)
  
  ## SAM files can have a variable number of column names, the col.names=1:25 is a possible
  ## point of failure if the SAM file has more than 25 columns
  output.samTable <- read.table(sam_path, sep='\t', col.names=1:25, stringsAsFactors=F, check.names=F, fill=T)
  
  ## Convert the dataframe to a datatable for faster access
  output.samTable <- as.data.table(output.samTable)
  
  ## Name the columns that are used for downstream analysis
  colnames(output.samTable)[1] <- 'read_name'
  colnames(output.samTable)[3] <- 'reference_name'
  
  ## Convert the alignment scores into integers
  alignmentScoreList <- as.integer(tstrsplit(output.samTable$`12`, ':', fixed=TRUE)[[3]])
  
  ## Save alignment scores to their own columns
  output.samTable$alignment_score <- alignmentScoreList
  
  cat('\n\nRemoving read alignments that do not include alignment scores')
  
  ## Discard any rows of the samTable that do not have alignment scores (can happen for various reasons)
  output.samTable <- output.samTable[!is.na(alignmentScoreList),]
  
  ## Check if there were any formatting errors with the alignment score conversions
  if(any(is.na(output.samTable$alignment_score))){
    stop("NA found in SAM file alignment scores for ", currentSample$name)
  }
  
  return(output.samTable)
}

## This function returns the list of alleles found in the reference fasta
read.kir_allele_list_from_reference_fasta <- function(fasta_path){
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  output.alleleList <- c()
  
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      output.alleleList <- c(output.alleleList, strsplit(currentLine, '>',fixed=TRUE)[[1]][2])
    }
  }
  
  return(output.alleleList)
}

## This function counts how many reads map to a unique locus or allele
run.count_kir_read_matches <- function(currentSample, samTable, maxReadThreshold){
  
  ## Pull out the unique read names
  uniqueReadNames <- unique(samTable$read_name)
  
  ## Randomize the read name order
  set.seed(001) # just to make it reproducible
  randomUniqueReadNames <- sample(uniqueReadNames)
  
  ## Check if there are more reads than the threshold and take some out if so
  if(length(randomUniqueReadNames) > maxReadThreshold){
    randomUniqueReadNames <- randomUniqueReadNames[1:maxReadThreshold]
  }
  
  ## Initialize the list for storing reference matches
  uniqueLocusMatchList <- as.list(kirLocusList)
  names(uniqueLocusMatchList) <- kirLocusList
  uniqueLocusMatchList[kirLocusList] <- 0
  
  ## Initialize the list for storing reference matches
  uniqueAlleleMatchList <- as.list(kirAlleleListRes3)
  names(uniqueAlleleMatchList) <- kirAlleleListRes3
  uniqueAlleleMatchList[kirAlleleListRes3] <- 0
  
  ## Initialize variables to check on the progress of the next for loop
  i = 1
  max_i = length(randomUniqueReadNames)
  checkAmount = ceiling(max_i/10)
  j=0
  
  ## Find the reference matches for each unique read name
  for(currentReadName in randomUniqueReadNames){
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    ## Find the best alignment score for this read
    maxAlignmentScore <- max(samSubsetTable$alignment_score)
    
    ## Pull out the current read name alignments that have the best alignment score
    samSubsetTable <- samSubsetTable[alignment_score == maxAlignmentScore]
    
    ## Read in the matched references to a list for further processing
    matchedAlleleList <- samSubsetTable$reference_name
    
    ## Pull out the res 3 allele names
    matchedAlleleList <- unique(unlist(lapply(matchedAlleleList, kir.allele_resolution, 3)))
    
    ## Pull out the unique locus names from the matched allele list
    matchedLocusList <- unique(unlist(lapply(matchedAlleleList, kir.allele_resolution, 0)))
    
    ## If there is only 1 unique locus, then add 1 to the unique match count for that locus
    if(length(matchedLocusList) <= 2){
      
      if(length(matchedLocusList) == 1){
        uniqueLocusMatchList[matchedLocusList] = uniqueLocusMatchList[matchedLocusList][[1]] + 1
        
        ## Adding in some special processing for 2DL5 because of the A/B naming scheme
      }else if('KIR2DL5A' %in% matchedLocusList & 'KIR2DL5B' %in% matchedLocusList){
        uniqueLocusMatchList['KIR2DL5'] = uniqueLocusMatchList['KIR2DL5'][[1]] + 1
      }
      
      ## If there is only a single matching reference allele, then iterate the count of that allele
      if(length(matchedAlleleList) == 1){
        uniqueAlleleMatchList[matchedAlleleList] = uniqueAlleleMatchList[matchedAlleleList][[1]] + 1
      }
    }
    
    ## Will display the percent completion every 10%
    i = i+1
    if(i%%checkAmount == 0){
      j = j+10
      cat(paste0(j,'% ', collapse = ''))
    }
  }
  
  cat("\n\nFinished counting!")
  
  ## Cutting the function short for now to test what a good threshold would be
  return(list(locusMatches = uniqueLocusMatchList, alleleMatches = uniqueAlleleMatchList))
  #locusRatio <- sapply(uniqueLocusMatchList, function(x) x/uniqueLocusMatchList$KIR3DL3)
}

## This function shortens KIR allele names
kir.allele_resolution <- function(allele_name, res){
  
  ## Split the locus name from the allele number
  alleleLocusNumber <- str_split(allele_name, fixed('*'))[[1]]
  
  alleleLocus <- alleleLocusNumber[1]
  alleleNumber <- alleleLocusNumber[2]
  
  ## Only return the alleleLocus if the resolution is 0
  if(res == 0){
    return(alleleLocus)
  }
  
  ## Subset the allele_number based on res
  shortAlleleNumber <- substr(alleleNumber, 1, res)
  
  ## Return the shortened allele name
  return(paste0(alleleLocus,'*',shortAlleleNumber))
}

## This function generates copy number graphs
run.generate_copy_number_graphs <- function(countRatioDF, kffDF){
  
  ## Iterate over all KIR loci, create a plot for each one
  for(currentLocus in kirLocusList){
    
    ## Determine the rank of the ratios (x-axis order from lowest to highest)
    countRatioDF$ratioRank <- rank(countRatioDF[,currentLocus], ties.method = 'first')
    
    ## Special dual locus graphing for KIR loci in tight LD
    if(currentLocus == 'KIR2DL1'){
      
      ## Graph 2DL1 and 2DP1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DL1', 'KIR2DP1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR2DL1' = '#000000', 'KIR2DP1' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR2DL1'=20,'KIR2DP1'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR2DP1'){
      
      ## Graph 2DP1 and 2DL1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DP1', 'KIR2DL1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR2DP1' = '#000000', 'KIR2DL1' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR2DP1'=20,'KIR2DL1'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR2DL2'){
      
      ## Graph 2DL2 and 2DL3 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DL2', 'KIR2DL3', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR2DL2' = '#000000', 'KIR2DL3' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR2DL2'=20,'KIR2DL3'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR2DL3'){
      
      ## Graph 2DL3 and 2DL2 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DL3', 'KIR2DL2', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR2DL3' = '#000000', 'KIR2DL2' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR2DL3'=20,'KIR2DL2'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR3DL1'){
      
      ## Graph 3DL1 and 3DS1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR3DL1', 'KIR3DS1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR3DL1' = '#000000', 'KIR3DS1' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR3DL1'=20,'KIR3DS1'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR3DS1'){
      
      ## Graph 3DS1 and 3DL1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR3DS1', 'KIR3DL1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR3DS1' = '#000000', 'KIR3DL1' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR3DS1'=20,'KIR3DL1'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else{
      
      ## Graph all other loci by themselves
      print(ggplot(countRatioDF, aes(x=ratioRank, y=countRatioDF[[currentLocus]])) +
              geom_point() +
              scale_color_manual(values=c(currentLocus = '#000000')) + 
              scale_shape_manual(values=c(currentLocus = 20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF[[currentLocus]]))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }
    
    ## Save each plot
    ggsave(file.path(resultsDirectory, paste0(currentLocus, '_copy_number_plot.pdf')), width=10, height=10)
  }
}

## This function performs a nucleotide string search and count to determine
run.count_kff_probes <- function(currentSample, probelistDF, maxReadThreshold){
  
  ## Initialize the list for kff probe match counts
  kffProbeMatchList <- as.list(probelistDF$Name)
  names(kffProbeMatchList) <- probelistDF$Name
  kffProbeMatchList[probelistDF$Name] <- 0
  
  ## Read in the first paired-end fastq file
  if(currentSample$gzip){
    fileContents <- fread(paste("zcat", currentSample$fastq1path), sep='\n', header=F, nrows=maxReadThreshold*4)
  }else{
    fileContents <- fread(currentSample$fastq1path, sep='\n', header=F, nrows=maxReadThreshold*4)
  }
  
  ## Pull out the read rows from the file object
  fileContents <- fileContents[seq(2, length(fileContents[[1]]), 4)]
  
  ## Count the number of grep hits for the probe sequence in the file object
  for(probeName in probelistDF$Name){
    kffProbeMatchList[probeName] <- length(grep(probelistDF[probeName,'Sequence'], fileContents[[1]], fixed=T)) + kffProbeMatchList[probeName][[1]]
  }
  
  ## Clean up
  remove(fileContents)
  
  ## Read in the second paired-end fastq file
  if(currentSample$gzip){
    fileContents <- fread(paste("zcat", currentSample$fastq2path), sep='\n', header=F, nrows=maxReadThreshold*4)
  }else{
    fileContents <- fread(currentSample$fastq2path, sep='\n', header=F, nrows=maxReadThreshold*4)
  }
  
  ## Pull out the read rows from the file object
  fileContents <- fileContents[seq(2, length(fileContents[[1]]), 4)]
  
  ## Count the number of grep hits for the probe sequence in the file object
  for(probeName in probelistDF$Name){
    kffProbeMatchList[probeName] <- length(grep(probelistDF[probeName,'Sequence'], fileContents[[1]], fixed=T)) + kffProbeMatchList[probeName][[1]]
  }
  
  ## Clean up
  remove(fileContents)
  
  return(kffProbeMatchList)
}

## This function normalizes the KFF probe counts by the number of probes and by the KIR3DL3 probe count
run.reduce_and_normalize_kff_probes <- function(probeCountList, probeLociList){

  ## Initialize the list for locus match counts
  locusKffRatioList <- as.list(probeLociList)
  names(locusKffRatioList) <- probeLociList
  locusKffRatioList[probeLociList] <- 0
  
  ## Pull out the KIR3DL3 specific probes
  allKIR3DL3Probes <- grep('>KIR3DL3>', names(kffCountList), fixed=T, value=T)
  
  ## Sum the KIR3DL3 probe counts, and divide by the number of probes to establish a KIR3DL3 normalization value
  KIR3DL3NormValue <- sum(as.integer(kffCountList[allKIR3DL3Probes]))/length(allKIR3DL3Probes)
  
  for(currentLocus in probeLociList){
    ## pull out the currentLocus specific probes
    allCurrentLocusProbes <- grep(paste0('>',currentLocus,'>'), names(kffCountList), fixed=T, value=T)
    
    ## Sum the currentLocus probe counts and divide by the number of probes to establish a currentLocus normalization value
    currentLocusNormValue <- sum(as.integer(kffCountList[allCurrentLocusProbes]))/length(allCurrentLocusProbes)
    
    ## Normalize the normalized value of the current locus by the KIR3DL3 normalized value
    locusKffRatioList[currentLocus] <- currentLocusNormValue/KIR3DL3NormValue
  }
  
  return(locusKffRatioList)
}

## This function checks the kff normalized ratios against kffThreshold to determine locus presence/absence
run.kff_determine_presence_from_norm_values <- function(kffNormList, kffThreshold){

  ## Initialize the list for locus match counts
  locusPresenceList <- kffNormList
  locusPresenceList[names(locusPresenceList)] <- 0
  
  ## Check each locus to see if the ratio is greater than or equal to kffThreshold
  for(currentLocus in names(kffNormList)){
    locusPresenceList[currentLocus] <- (kffNormList[[currentLocus]]>=kffThreshold)*1
  }
  
  return(locusPresenceList)
}
