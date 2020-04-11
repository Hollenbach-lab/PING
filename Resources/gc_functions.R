#library(methods)
#library(data.table)
#library(ggplot2)
#library(plotly)
kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

## This function checks to make sure the output of system2 is valid
#check.system2_output <- function(system2_output, system2_error){
#  
#  ## Checking if the attributes of system2_command are NULL, if not the command was not found
#  if(!is.null(attributes(system2_output))){
##    cat('\n',system2_output)
#    stop(system2_error, '. Stopping program.')
#  }
#}

## Check to make sure samtools is accessible
#samtools <- system2('which', c('samtools'), stdout=T, stderr=T)
#check.system2_output(samtools, 'samtools not found')

## This function finds samples in sampleDirectory and turns them into objects for downstream use
# build.paired_sample_objects <- function(sample_directory, fastq_pattern='fastq', resultsDirectory){
#   #####
#   ## This function takes in a directory and a file name pattern and attempts to pair fastq files
#   ## Returns a list of sample objects that contain the paired fastq file names
#   #####
#   
#   cat("\nAttempting automatic fastq pairing in", sample_directory, "using", fastq_pattern)
#   
#   ## Find all the files in sampleDirectory that match fastqPattern
#   unpairedFastqList <- list.files(path=sample_directory, pattern=fastq_pattern)
#   
#   ## To pair reads, we will split the file names by fastqPattern, then continuously chop a
#   ## character off the end of each name until the number of unique names is exactly half of the total names
#   
#   ## Setting up an initial fastq list that splits the files names by fastqPattern
#   strList <- sapply(unpairedFastqList, function(x) str_split(x, fastq_pattern)[[1]][1])
#   
#   ## Setting the maximum number of times to chop off the last character to the length of the shortest name
#   maxChop <- min(sapply(strList, nchar))
#   
#   ## Iterate from 0 to maxChop, cutting off i characters from the file names each time
#   for(i in 0:maxChop){
#     
#     ## In each iteration, the file names are reset. There is no particular reason I implemented it this way
#     subStrList <- strList
#     
#     ## Cut off i characters from the end of each fastq file name
#     subStrList <- sapply(subStrList, function(x) substr(x, 1, nchar(x)-i))
#     
#     ## After cutting, determine the unique names
#     uniqueFastqList <- unique(subStrList)
#     
#     ## If the number of unique names is exactly half of the total names, then cutting should be finished
#     ## and it is time to move on to matching!
#     if(length(uniqueFastqList) == (length(subStrList)/2)){
#       break
#     }
#     
#     ## Pairing failed if i reaches maxChop. This will raise a fatal error.
#     if(i == maxChop){
#       stop("Was not able to pair fastq file names, please check that fastqPattern and sampleDirectory are set correctly.")
#     }
#   }
#   
#   ## Initialize a list for storing the paired fastq file names
#   pairedFastqList <- list()
#   
#   ## Iterate through the unique fastq names to make pairings
#   for(fastqName in uniqueFastqList){
#     
#     ## Pull out the matches for fastqName in the subStrList
#     fastqMatches <- subStrList[fastqName == subStrList]
#     
#     ## Determine how many matches there are for fastqName in the subStrList
#     matchCount <- length(fastqMatches)
#     
#     ## Stop the program if more or less than 2 matches are found for fastqName
#     if(matchCount != 2){
#       cat('\n',names(fastqMatches))
#       stop('Auto fastq matching failed due to an improper number of matches for ',fastqName)
#     }
#     
#     ## Save the file names in the pairedFastqList under the unique name
#     pairedFastqList[[fastqName]] <- names(fastqMatches)
#   }
#   
#   cat("\nFound", length(uniqueFastqList), "samples in", sample_directory)
#   
#   ## Creating the sample object class
#   sample <- setRefClass("sample",
#                         fields=list(name='character',
#                                     fastq1path='character',
#                                     fastq2path='character',
#                                     gzip='logical',
#                                     samPath='character',
#                                     bamPath='character',
#                                     assembledRefPath='character',
#                                     assembledRefIndex='character',
#                                     assembledSamPath='character',
#                                     assembledBedPath='character',
#                                     assembledDelIndex='character'))
#   
#   ## Initializing a sample object list. This will be returned
#   output.sampleList <- list()
#   
#   for(i in 1:length(pairedFastqList)){
#     
#     ## Pulling the current working element out of the list
#     pairedFastq <- pairedFastqList[i]
#     
#     ## Creating an absolute path to the first fastq file
#     fastq1path <- normalizePath(file.path(sample_directory, pairedFastq[[1]][1]), mustWork=T)
#     
#     ## Creating a absolute path to the second fastq file
#     fastq2path <- normalizePath(file.path(sample_directory, pairedFastq[[1]][2]), mustWork=T)
#     
#     ## Checking if the first fastq file is gzipped
#     gzip <- substr(fastq1path, nchar(fastq1path)-2, nchar(fastq1path)) == '.gz'
#     
#     ## Fill in the path to the alignment file (it may or may not be present)
#     samPath <- file.path(resultsDirectory,paste0(names(pairedFastq),'.sam'))
#     
#     ## Fill in the path to the alignment file (it may or may not be present)
#     bamPath <- file.path(resultsDirectory,paste0(names(pairedFastq),'.bam'))
#     
#     ## Building a sample object and adding it to sampleList
#     output.sampleList[[names(pairedFastq)]] <- sample(name=names(pairedFastq),fastq1path=fastq1path,fastq2path=fastq2path,gzip=gzip,samPath=samPath,bamPath=bamPath)
#   }
#   
#   cat("\nAll samples were successfully paired")
#   return(output.sampleList)
# }

## This function runs a bowtie2 alignment looking for exact matches to all KIR references in subsetKirReference
run.bowtie2_gc_alignment <- function(bowtie2_command, reference_index, threads, current_sample, bamDirectory){
  ## Intitialize an output path for the SAM file
  current_sample$samPath <- file.path(bamDirectory,paste0(current_sample$name,'.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-5 0', '-3 6', '-N 0', '--end-to-end', paste0('-p ',threads), '--score-min "C,-2,0"',
                      '-I 75', '-X 1000',
                      paste0('-1 ',current_sample$kirfastq1path),
                      paste0('-2 ',current_sample$kirfastq2path),
                      '--no-unal','-a','--np 0', '--mp 2,2', '--rdg 1,1', '--rfg 1,1',
                      paste0('-S ',current_sample$samPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 failed, retrying alignment...')
    output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  }
  
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  current_sample$samPath <- normalizePath(current_sample$samPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',current_sample$name,'to',reference_index)
  
  return(current_sample)
}

## This function counts the number of header lines in a SAM file (using '@SQ')
samfile.count_header_lines <- function(currentSample){
  if(!file.exists(currentSample$samPath)){
    stop('This sam file does not exist')
  }
  
  headerLines <- 0
  con = file(currentSample$samPath, "r")
  
  while(TRUE){
    
    line = readLines(con, n = 1)
    #if(strsplit(line,'\t')[[1]][1] != "@SQ"){
    #  break
    #}
    
    if(substr(line,1,1) != '@'){
      break
    }
    
    headerLines <- headerLines + 1
  }
  close(con)
  
  return(headerLines)
}

## This function runs a samtools sam to bam conversion
samtools.sam_to_bam <- function(samtools_command, currentSample, bamDirectory, threads){
  
  ## Initialize an output path for the BAM file
  currentSample[['bamPath']] <- file.path(bamDirectory,paste0(currentSample$name,'.bam'))
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),
                      currentSample$samPath, '-o', currentSample$bamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamConv, 'samtools sam to bam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['bamPath']] <- normalizePath(currentSample$bamPath, mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample$samPath,'to',currentSample$bamPath)
  
  return(currentSample)
}


## This function runs a samtools bam to sam conversion
samtools.bam_to_sam <- function(samtools_command, currentSample, bamDirectory, threads){
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),'-h',
                      currentSample$bamPath, '-o', currentSample$samPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.samConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.samConv, 'samtools bam to sam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.samConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['samPath']] <- normalizePath(currentSample$samPath, mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample$bamPath,'to',currentSample$samPath)
  
  return(currentSample)
}

## This function reads in a SAM file with no header to a data.table
read.bowtie2_sam_nohd <- function(currentSample, allele_alignment=F, rows_to_skip=0){
  cat("\n\nReading in",currentSample$samPath)
  
  ## Make sure the SAM file can be read in
  sam_path <- normalizePath(currentSample$samPath, mustWork=T)
  
  ## SAM files can have a variable number of column names, the col.names=1:25 is a possible
  ## point of failure if the SAM file has more than 25 columns
  #output.samTable <- read.table(sam_path, sep='\t', col.names=1:25, stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip, nrows = 30000000)
  output.samTable <- read.table(sam_path, sep='\t', stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip,col.names=1:25)
  
  ## Convert the dataframe to a datatable for faster access
  output.samTable <- as.data.table(output.samTable)
  
  ## Name the columns that are used for downstream analysis
  colnames(output.samTable)[1] <- 'read_name'
  colnames(output.samTable)[3] <- 'reference_name'
  colnames(output.samTable)[4] <- 'ref_pos'
  colnames(output.samTable)[10] <- 'read_seq'
  
  ## Convert the alignment scores into integers
  alignmentScoreList <- as.integer(tstrsplit(output.samTable$`12`, ':', fixed=TRUE)[[3]])
  
  ## Save alignment scores to their own columns
  cat('\nSetting alignment scores.')
  output.samTable[,alignment_score := alignmentScoreList]
  
  ## Convert the allele names into locus names
  if(allele_alignment){
    locusList <- tstrsplit(output.samTable$reference_name, '_', fixed=T)[[1]]
  }else{
    locusList <- tstrsplit(output.samTable$reference_name, '*', fixed=TRUE)[[1]]
  }
  locusList <- tstrsplit(output.samTable$reference_name, '*', fixed=TRUE)[[1]]
  locusList[locusList %in% 'KIR2DL5A'] = 'KIR2DL5'
  locusList[locusList %in% 'KIR2DL5B'] = 'KIR2DL5'

  ## Save locus names to their own columns
  cat('\nSetting locus names.')
  output.samTable[,locus:=locusList]
  
  ## Create new column names to store universal coordinates
  output.samTable$startPos <- 0
  output.samTable$endPos <- 0
  
  ## Create a new column to store read lengths
  cat('\nSetting read lengths.')
  output.samTable[,readLen := nchar(read_seq)]
  
  cat('\nRemoving read alignments that do not include alignment scores.')
  
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

## This function returns a list of dataframes of allele sequences found in the reference fasta
read.kir_allele_dataframe_from_reference_fasta <- function(fasta_path, kirLocusList){
  
  ## Make sure the fasta can be found
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  ## Initialize a list to store the sequence strings from the file
  alleleSeqList <- list()
  
  ## Read in the fasta file and store the allele names and sequences
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      alleleName <- strsplit(currentLine, '>',fixed=TRUE)[[1]][2]
    }else{
      alleleSeq <- currentLine
      alleleSeqList[[alleleName]] <- alleleSeq
    }
  }
  
  ## Initialize the output list
  output.alleleSeqDFList <- list()
  
  cat('\n\tProcessing:')
  
  ## Input the allele sequence strings to a dataframe for each locus
  for(currentLocus in kirLocusList){
    
    cat('',currentLocus)
    
    if(currentLocus == 'KIR2DL5'){
      tempCurrentLocus <- c('KIR2DL5A', 'KIR2DL5B')
    }else{
      tempCurrentLocus <- currentLocus
    }
    
    ## Pull out all of the allele names for the current locus
    currentLocusAlleleNames <- names(alleleSeqList)[tstrsplit(names(alleleSeqList), '*', fixed=T)[[1]] %in% tempCurrentLocus]
    
    ## Check the size of the largest and smalles allele
    largestAlleleSize <- max(sapply(alleleSeqList[currentLocusAlleleNames], nchar))
    minAlleleSize <- min(sapply(alleleSeqList[currentLocusAlleleNames], nchar))
    
    ## Make sure the largest and smallest allele are the same sizes
    if(largestAlleleSize != minAlleleSize){
      stop(currentLocus,' has alleles of different sizes. Cannot transform into dataframe.')
    }
    
    ## Initialize a dataframe for storing the allele sequence strings
    alleleSeqDF <- data.frame(matrix(0, length(currentLocusAlleleNames), largestAlleleSize),row.names=currentLocusAlleleNames,check.names=F,stringsAsFactors=F)
    
    ## For each allele for the current locus, input the sequence into the dataframe
    for(currentAlleleName in currentLocusAlleleNames){
      alleleSeqDF[currentAlleleName,] <- strsplit(alleleSeqList[[currentAlleleName]],'')[[1]]
    }
    
    ## Save the dataframe to the output list
    output.alleleSeqDFList[[currentLocus]] <- alleleSeqDF
  }
  
  return(output.alleleSeqDFList)
}

## This function counts how many reads map to a unique locus or allele
run.count_kir_read_matches <- function(currentSample, samTable, maxReadThreshold, kirLocusList, kirAlleleList){
  
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
  #uniqueAlleleMatchList <- as.list(kirAlleleList)
  #names(uniqueAlleleMatchList) <- kirAlleleList
  #uniqueAlleleMatchList[kirAlleleList] <- 0
  
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
    #matchedAlleleList <- samSubsetTable$reference_name
    
    ## Pull out the res 3 allele names
    #matchedAlleleList <- unique(unlist(lapply(matchedAlleleList, kir.allele_resolution, 3)))
    
    ## Pull out the unique locus names from the matched allele list
    matchedLocusList <- unique(samSubsetTable$locus)
    
    ###### This section will count all locus matches (as opposed to only unique locus matches)
    #for(matchedLocus in matchedLocusList){
    #  uniqueLocusMatchList[matchedLocus] = uniqueLocusMatchList[matchedLocus][[1]] + 1
    #}
    #if('KIR2DL5A' %in% matchedLocusList & 'KIR2DL5B' %in% matchedLocusList){
    #  uniqueLocusMatchList['KIR2DL5'] = uniqueLocusMatchList['KIR2DL5'][[1]] + 1
    #}
    ###### /s
    
    ###### This section will count only unique locus matches (as opposed to all locus matches)
    ## If there is only 1 unique locus, then add 1 to the unique match count for that locus
    
    if(length(matchedLocusList) == 1){
      uniqueLocusMatchList[matchedLocusList] = uniqueLocusMatchList[matchedLocusList][[1]] + 1
        
      ## This will count all allele matches
      #for(matchedAllele in matchedAlleleList){
      #  uniqueAlleleMatchList[matchedAllele] = uniqueAlleleMatchList[matchedAllele][[1]] + 1
      #}
      
      ## If there is only a single matching reference allele, then iterate the count of that allele
      #if(length(matchedAlleleList) == 1){
      #  uniqueAlleleMatchList[matchedAlleleList] = uniqueAlleleMatchList[matchedAlleleList][[1]] + 1
      #}
    }
    ###### /s
    
    ## Will display the percent completion every 10%
    i = i+1
    if(i%%checkAmount == 0){
      j = j+10
      cat(paste0(j,'% ', collapse = ''))
    }
  }
  cat("\n\nFinished counting!")
  
  #return(list(locusMatches = uniqueLocusMatchList, alleleMatches = uniqueAlleleMatchList))
  return(list(locusMatches=uniqueLocusMatchList))
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

## Generates copy number graphs with threshold lines
run.generate_copy_number_graphs <- function(countRatioDF, kffDF, kirLocusList, plotDirectory, countDF, thresholdDF){
  
  ## Check to see what samples are in both data frames
  samplesInBoth <- intersect(row.names(countRatioDF), row.names(kffDF))
  
  ## If there are no samples in both data frames, stop the script
  if(length(samplesInBoth) < 1){
    stop('\n\nThe kff presence table does not match up with the locus count table.')
  }
  
  ## Subset both data frames by the samples found in both
  countRatioDF <- countRatioDF[samplesInBoth,]
  kffDF <- kffDF[samplesInBoth,]
  countDF <- countDF[samplesInBoth,]
  
  ## Give sample id their own column
  countRatioDF$id <- rownames(countRatioDF)
  
  ## Iterate over all KIR loci, create a plot for each one
  for(currentLocus in kirLocusList){
    
    if(currentLocus =='KIR3DL3'){
      next
    }
    
    ## Format threshold stuff
    currentLocusThresholds <- thresholdDF[currentLocus,]
    setThresholdCols <- colnames(thresholdDF)[!is.na(currentLocusThresholds)]
    setThresholdValList <- as.list(currentLocusThresholds[setThresholdCols])
    

    ## Determine the rank of the ratios (x-axis order from lowest to highest)
    countRatioDF$ratioRank <- rank(countRatioDF[,currentLocus], ties.method = 'first')
    
    ## Initialize columns for storing kff presence
    countRatioDF$kffPresence <- 0
    countRatioDF$altKffPresence <- 0
    
    ## Initialize neg locus name for currentLocus
    currentNeg <- paste0(currentLocus,'_kff_neg')
    
    ## Initialize 3DL3 ratio to max 3DL3 count
    countRatioDF$overall3DL3Ratio <- 0
    
    ## Set the currentLocus kff presence values
    countRatioDF[samplesInBoth,'kffPresence'] <- kffDF[samplesInBoth,currentLocus]
    
    ## Set the 3DL3 ratio to max 3DL3 count
    countRatioDF[samplesInBoth,'overall3DL3Ratio'] <- countDF[samplesInBoth,'KIR3DL3']/max(countDF[samplesInBoth,'KIR3DL3'])
    
    ## Find the sample id's that match each condition for currentLocus
    currentLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'kffPresence'] > 0]
    currentLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'kffPresence'] == 0]
    
    ## Special dual locus graphing for KIR loci in tight LD
    if(currentLocus == 'KIR2DL1'){
      
      ## Initialze the alternate locus to compare to
      altLocus<-'KIR2DP1'
      
      ## Initialize neg locus names for alt locus
      altNeg <- paste0(altLocus,'_kff_neg')
      
      ## Set the alternate locus kff presence values
      countRatioDF[samplesInBoth,'altKffPresence'] <- kffDF[samplesInBoth,altLocus]
      
      ## Find the sample id's that match each condition for altLocus
      altLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] > 0]
      altLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] == 0]
      
    }else if(currentLocus == 'KIR2DP1'){
      
      ## Initialze the alternate locus to compare to
      altLocus<-'KIR2DL1'
      
      ## Initialize neg locus names for alt locus
      altNeg <- paste0(altLocus,'_kff_neg')
      
      ## Set the alternate locus kff presence values
      countRatioDF[samplesInBoth,'altKffPresence'] <- kffDF[samplesInBoth,altLocus]
      
      ## Find the sample id's that match each condition for altLocus
      altLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] > 0]
      altLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] == 0]
      
    }else if(currentLocus == 'KIR2DL2'){
      
      ## Initialze the alternate locus to compare to
      altLocus<-'KIR2DL3'
      
      ## Initialize neg locus names for alt locus
      altNeg <- paste0(altLocus,'_kff_neg')
      
      ## Set the alternate locus kff presence values
      countRatioDF[samplesInBoth,'altKffPresence'] <- kffDF[samplesInBoth,altLocus]
      
      ## Find the sample id's that match each condition for altLocus
      altLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] > 0]
      altLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] == 0]
      
    }else if(currentLocus == 'KIR2DL3'){
      
      ## Initialze the alternate locus to compare to
      altLocus<-'KIR2DL2'
      
      ## Initialize neg locus names for alt locus
      altNeg <- paste0(altLocus,'_kff_neg')
      
      ## Set the alternate locus kff presence values
      countRatioDF[samplesInBoth,'altKffPresence'] <- kffDF[samplesInBoth,altLocus]
      
      ## Find the sample id's that match each condition for altLocus
      altLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] > 0]
      altLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] == 0]
      
    }else if(currentLocus == 'KIR3DL1'){
      
      ## Initialze the alternate locus to compare to
      altLocus<-'KIR3DS1'
      
      ## Initialize neg locus names for alt locus
      altNeg <- paste0(altLocus,'_kff_neg')
      
      ## Set the alternate locus kff presence values
      countRatioDF[samplesInBoth,'altKffPresence'] <- kffDF[samplesInBoth,altLocus]
      
      ## Find the sample id's that match each condition for altLocus
      altLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] > 0]
      altLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] == 0]
      
    }else if(currentLocus == 'KIR3DS1'){
      
      ## Initialze the alternate locus to compare to
      altLocus<-'KIR3DL1'
      
      ## Initialize neg locus names for alt locus
      altNeg <- paste0(altLocus,'_kff_neg')
      
      ## Set the alternate locus kff presence values
      countRatioDF[samplesInBoth,'altKffPresence'] <- kffDF[samplesInBoth,altLocus]
      
      ## Find the sample id's that match each condition for altLocus
      altLocusPresent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] > 0]
      altLocusAbsent <- samplesInBoth[countRatioDF[samplesInBoth,'altKffPresence'] == 0]
      
    }else{
      ## Set no alternative loci
      altLocus=''
      altNeg=''
      
      ## Set no samples names for alternative loci
      altLocusPresent <- currentLocus[0]
      altLocusAbsent <- currentLocus[0]
      kir3DL3 <- samplesInBoth[countRatioDF[samplesInBoth,'']]
    }
    kir3DL3RatioToMax <- 'KIR3DL3_ratio_to_max3DL3'
    
    ## Initialize color scheme for plots
    pal1 <- c("black","red")
    pal1 <- setNames(pal1, c(currentLocus, currentNeg))
    
    pal2 <- c("#d3dcff", "gray", "pink", "red", "black")
    pal2 <- setNames(pal2, c(kir3DL3RatioToMax, altLocus, altNeg, currentNeg, currentLocus))
    
    palThresholds <- c('#bb00ff','#2eee00','#007aff','#ffb200','#ff1aa7','#dbdd07')
    palThresholds <- setNames(palThresholds, colnames(thresholdDF))
    
    ## Set labels for the positive and negative points
    posPointText <- paste('Sample ID:',countRatioDF[currentLocusPresent,'id'],'$<br>Ratio:',countRatioDF[currentLocusPresent,currentLocus],'$<br>3DL3Ratio:',countRatioDF[currentLocusPresent,'overall3DL3Ratio'])
    negPointText <- paste('Sample ID:',countRatioDF[currentLocusAbsent,'id'],'$<br>Ratio:',countRatioDF[currentLocusAbsent,currentLocus],'$<br>3DL3Ratio:',countRatioDF[currentLocusAbsent,'overall3DL3Ratio'])
    altPosPointText <- paste('Sample ID:',countRatioDF[altLocusPresent,'id'],'$<br>Ratio:',countRatioDF[altLocusPresent,altLocus])
    altNegPointText <- paste('Sample ID:',countRatioDF[altLocusAbsent,'id'],'$<br>Ratio:',countRatioDF[altLocusAbsent,altLocus])
    kir3DL3RatioText <- paste('Sample ID:',countRatioDF[samplesInBoth,'id'],'$<br>Ratio:',countRatioDF[samplesInBoth,'overall3DL3Ratio'])
    
    ## Set a title for the graph
    currentLocusTitle <- currentLocus
    
    ## Fix for graph breaking when there are no positive points for current locus
    if(length(currentLocusPresent)==0){
      currentLocusTitle=currentLocusTitle[0]
      posPointText=posPointText[0]
    }
    
    ## Fix for graph breaking when there are no negative points for current locus
    if(length(currentLocusAbsent)==0){  
      currentNeg=currentNeg[0]
      negPointText=negPointText[0]
    }
    
    ## Fix for graph breaking when there are no positive points for alt locus
    if(length(altLocusPresent)==0){
      altLocus=altLocus[0]
      altPosPointText=altPosPointText[0]
    }
    
    ## Fix for graph breaking when there are no negative points for alt locus
    if(length(altLocusAbsent)==0){
      altNeg=altNeg[0]
      altNegPointText=altNegPointText[0]
    }
    
    ## Snippet of code to remove list elements with empty names from the color pallete
    i=1
    removeList = c()
    for(keyName in names(pal1)){
      if(is.na(keyName)){
        removeList <- c(removeList, i)
      }else if(nchar(keyName)==0){
        removeList <- c(removeList, i)
      }
      i = i+1
    }
    if(length(removeList>0)){
      pal1 <- pal1[-removeList]
    }
    
    i=1
    removeList = c()
    for(keyName in names(pal2)){
      if(is.na(keyName)){
        removeList <- c(removeList, i)
      }else if(nchar(keyName)==0){
        removeList <- c(removeList, i)
      }
      i = i+1
    }
    if(length(removeList>0)){
      pal2 <- pal2[-removeList]
    }
    
    maxY <- max(c(1,unlist(countRatioDF[samplesInBoth,c(currentLocus,altLocus)])))+0.2
    
    ## Create plot
    
    ## Add histogram data
    p1 <- plot_ly(colors=pal1) %>%
      add_histogram(y=~countRatioDF[currentLocusPresent,currentLocusTitle],name=currentLocusTitle,
                    color=currentLocusTitle) %>%
      add_histogram(y=~countRatioDF[currentLocusAbsent,currentLocus],name=currentNeg,
                    color=currentNeg) %>%
      layout(title='KIR3DL3',
             xaxis=list(title='density',showgrid=F),
             yaxis=list(title=paste(currentLocus,'/ KIR3DL3 Ratio'),range=c(0,maxY)))
    
    ## Add ratio data
    p2 <- plot_ly(colors=pal2) %>%
      add_trace(x=countRatioDF[samplesInBoth,'ratioRank'],
                y=countRatioDF[samplesInBoth,'overall3DL3Ratio'],
                mode='markers',type='scatter',name=kir3DL3RatioToMax,color=kir3DL3RatioToMax,
                text=kir3DL3RatioText) %>%
      add_trace(x=countRatioDF[altLocusPresent, 'ratioRank'],
                y=countRatioDF[altLocusPresent,altLocus], 
                mode='markers',type='scatter',name=altLocus,color=altLocus,
                text=altPosPointText) %>%
      add_trace(x=countRatioDF[altLocusAbsent,'ratioRank'],
                y=countRatioDF[altLocusAbsent,altLocus],
                mode='markers',type='scatter',name=altNeg,color=altNeg,
                text=altNegPointText) %>%
      add_trace(x=countRatioDF[currentLocusAbsent, 'ratioRank'],
                y=countRatioDF[currentLocusAbsent,currentLocus], 
                mode='markers',type='scatter',name=currentNeg,color=currentNeg,
                text=negPointText) %>%
      add_trace(x=countRatioDF[currentLocusPresent,'ratioRank'],
                y=countRatioDF[currentLocusPresent,currentLocusTitle],
                mode='markers',type='scatter',name=currentLocusTitle,color=currentLocusTitle,
                text=posPointText)
      
      ## Add theshold data
      for(currentThreshold in names(setThresholdValList)){
        ## Pull out the current threshold value
        thresholdVal <- as.numeric(setThresholdValList[[currentThreshold]])
        
        ## Initialize overlay text
        currentThreshText <- paste('Threshold:',currentThreshold,
                                   '$<br>Value:',thresholdVal)
        
        ## Add trace for current threhsold
        p2 <- p2 %>% add_trace(p2,
                        x=1:length(samplesInBoth),
                        y=thresholdVal,
                        type='scatter',
                        mode='line',
                        line=list(color=palThresholds[currentThreshold]),
                        name=currentThreshold)
      }
      
      ## Format axis layout
      p2 <- p2 %>% layout(title=currentLocus,
                          xaxis = list(title='Sample rank',showgrid=F),
                          yaxis = list(title='',range=c(0,maxY),showgrid=T))
    
    p <- subplot(p1, p2, shareY=T, widths = c(0.15,0.85),titleX=T)
    print(p)
    
    ## Save each plot
    htmlwidgets::saveWidget(p, file=file.path(plotDirectory,paste0(currentLocus, '_copy_number_plot.html')))
    
  }
}

## Sets copy number manually
run.set_copy <- function(kirLocusList, 
                         copyNumberDF, 
                         locusRatioDF,
                         locusCountDF,
                         thresholdDF){
  
  ## Set aside the samples that we do not know the copy
  unsetCopySampleList <- row.names(locusRatioDF)
  failedSampleList <- setdiff(row.names(locusCountDF),unsetCopySampleList)
  
  copyNumberDF[failedSampleList,] <- 'failed'
  
  for(kirLocus in kirLocusList){
    
    cat('\n\nProcessing',kirLocus)
    
    ## Set aside the samples that we do not know the copy
    unsetCopySampleList <- row.names(locusRatioDF)
    
    ## Skip KIR3DL3 (since that is the normalization locus)
    if(kirLocus == 'KIR3DL3'){
      cat('\nAuto-setting all KIR3DL3 copy to 2.')
      copyNumberDF[unsetCopySampleList,kirLocus] <- 2
      next
    }
    
    maxCopyInt <- readMaxCopyInt(kirLocus)
    
    if(maxCopyInt == 0){
      copyNumberDF[,kirLocus] <- 0
    }else{
      for(topCopy in 1:maxCopyInt){
        copyThresholdDouble <- readThresholdFloat(topCopy)
        
        ## save threshold, indexing could also be done using paste0(topCopy,'-',topCopy-1)
        thresholdDF[kirLocus,colnames(thresholdDF)[topCopy]] <- copyThresholdDouble
        
        ## Subset the sample list by the names of samples that fall in the lower copy group
        lowerSampleList <- unsetCopySampleList[locusRatioDF[unsetCopySampleList,kirLocus] < copyThresholdDouble]
        
        ## Subset the sample list by the names of samples that fall in the upper copy group
        upperSampleList <- unsetCopySampleList[locusRatioDF[unsetCopySampleList,kirLocus] >= copyThresholdDouble]
        
        ## Set copy number for samples that fall under the threshold
        copyNumberDF[lowerSampleList,kirLocus] <- topCopy-1
        
        ## Set copy number for samples that fall above the threshold
        copyNumberDF[upperSampleList,kirLocus] <- topCopy
        
        ## Remove sample names that had copy set
        unsetCopySampleList <- setdiff(unsetCopySampleList,lowerSampleList)
        
        ## Break out of the loop if we run out of samples to set
        if(length(unsetCopySampleList) == 0){
          break
        }
      }
    }
    
  }
  return(list('copyDF'=copyNumberDF,'threshDF'=thresholdDF))
}

## Helper function for run.set_copy()
readMaxCopyInt <- function(kirLocus){ 
  n <- readline(prompt=paste("Please enter the maximum copy number for",kirLocus,"as an integer: "))
  if(!grepl("^[0-9]+$",n)){
    return(readMaxCopyInt(kirLocus))
  }
  return(as.integer(n))
}

## Helper function for run.set_copy()
readThresholdFloat <- function(topCopy){
  n <- readline(prompt=paste("Please enter the threshold value that separates copy groups",topCopy-1,"and",topCopy,": "))
  n <- as.double(n)
  if(is.na(n)){
    n <- readThresholdFloat(topCopy)
  }
  return(n)
}

## This function generates copy number graphs
run.generate_predicted_copy_number_graphs <- function(countRatioDF, kirLocusList, plotDirectory, countDF, copyDF){
  
  ## Check to see what samples are in both data frames
  samplesInBoth <- intersect(row.names(countRatioDF), row.names(copyDF))
  
  ## If there are no samples in both data frames, stop the script
  if(length(samplesInBoth) < 1){
    stop('\n\nThe kff presence table does not match up with the locus count table.')
  }
  
  ## Subset both data frames by the samples found in both
  countRatioDF <- countRatioDF[samplesInBoth,]
  countDF <- countDF[samplesInBoth,]
  copyDF <- copyDF[samplesInBoth,]
  
  ## Give sample id their own column
  countRatioDF$id <- rownames(countRatioDF)
  
  ## Initialize 3DL3 ratio to max 3DL3 count
  countRatioDF$overall3DL3Ratio <- 0
  
  ## Set the 3DL3 ratio to max 3DL3 count
  countRatioDF[samplesInBoth,'overall3DL3Ratio'] <- countDF[samplesInBoth,'KIR3DL3']/max(countDF[samplesInBoth,'KIR3DL3'])
  
  ## Iterate over all KIR loci, create a plot for each one
  for(currentLocus in kirLocusList){
    
    if(currentLocus =='KIR3DL3'){
      next
    }
    
    ## Determine the rank of the ratios (x-axis order from lowest to highest)
    countRatioDF$ratioRank <- sample(rank(countRatioDF[,currentLocus], ties.method = 'first'))
    
    ## Set the predicted copy number for the current sample
    countRatioDF$copyNumber <- copyDF[,currentLocus]
    
    ## Set the color scheme
    pal <- c("#ff0000", "#000cff", "#ff00ff", "#00d6dd", "#cda425")
    pal <- setNames(pal, c('0', '1', '2', '3', '4'))
    pal <- pal[as.character(unique(countRatioDF$copyNumber))]
    pal[['KIR3DL3']] <- '#cdcdcd' 
    
    ## Set the maximum Y value
    maxY <- max(c(1,unlist(countRatioDF[samplesInBoth,currentLocus])))+0.2
    
    ## Initialize the plot
    p1 <- plot_ly(colors=pal)

    ## Plot the copy number predictions
    for(copyNumber in unique(countRatioDF$copyNumber)){
      

      ## Pull out the current copy sample names
      currentCopySamples <- row.names(countRatioDF)[countRatioDF[,'copyNumber'] == copyNumber]
      
      ## Initialize the overlay text
      currentCopyText <- paste('Sample ID:',countRatioDF[currentCopySamples,'id'],
                               '$<br>Ratio:',countRatioDF[currentCopySamples,currentLocus],
                               '$<br>KIR3DL3_ratio:',countRatioDF[currentCopySamples,'overall3DL3Ratio'])
      
      ## Add the trace for the current copy information
      #p1 <- add_trace(p1, 
      #                x=countRatioDF[currentCopySamples,'ratioRank'], 
      #                y=countRatioDF[currentCopySamples,'overall3DL3Ratio'],
      #                type="scatter",mode="markers",
      #                showlegend=F,
      #                color='KIR3DL3',
      #                text=currentCopyText)
      
      ## Add the trace for the current copy information
      p1 <- add_trace(p1, 
                      x=countRatioDF[currentCopySamples,'ratioRank'], 
                      y=countRatioDF[currentCopySamples,currentLocus],
                      type="scatter",mode="markers",
                      name=as.character(copyNumber),
                      color=as.character(copyNumber),
                      text=currentCopyText)
    }
    
    ## Format the layout of the graph
    p1 <- layout(p1,
                 title=currentLocus,
                 xaxis = list(title='Sample rank',showgrid=F),
                 yaxis = list(title='',range=c(0,maxY),showgrid=T))
    
    ## Print the graph
    print(p1)
    
    ## Save each plot
    htmlwidgets::saveWidget(p1, file=file.path(plotDirectory,paste0(currentLocus, '_predicted_copy_number_plot.html')))
  }
}

## This function performs a nucleotide string search and count to determine
run.count_kff_probes <- function(currentSample, probelistDF, maxReadThreshold){
  
  ## Initialize the list for kff probe match counts
  kffProbeMatchList <- as.list(probelistDF$Name)
  names(kffProbeMatchList) <- probelistDF$Name
  kffProbeMatchList[probelistDF$Name] <- 0
  
  ## Read in the first paired-end fastq file
  fileContents <- fread(cmd=(paste("zcat", currentSample$kirfastq1path)), sep='\n', header=F, nrows=maxReadThreshold*4)
  
  ## Pull out the read rows from the file object
  fileContents <- fileContents[seq(2, length(fileContents[[1]]), 4)]
  
  ## Count the number of grep hits for the probe sequence in the file object
  for(probeName in probelistDF$Name){
    kffProbeMatchList[probeName] <- length(grep(probelistDF[probeName,'Sequence'], fileContents[[1]], fixed=T)) + kffProbeMatchList[probeName][[1]]
  }
  
  ## Clean up
  remove(fileContents)
  
  ## Read in the second paired-end fastq file
  fileContents <- fread(cmd=(paste("zcat", currentSample$kirfastq2path)), sep='\n', header=F, nrows=maxReadThreshold*4)
  
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
  allKIR3DL3Probes <- grep('>KIR3DL3>', names(probeCountList), fixed=T, value=T)
  
  ## Sum the KIR3DL3 probe counts, and divide by the number of probes to establish a KIR3DL3 normalization value
  KIR3DL3NormValue <- sum(as.integer(probeCountList[allKIR3DL3Probes]))/length(allKIR3DL3Probes)
  
  for(currentLocus in probeLociList){
    ## pull out the currentLocus specific probes
    allCurrentLocusProbes <- grep(paste0('>',currentLocus,'>'), names(probeCountList), fixed=T, value=T)
    
    ## Sum the currentLocus probe counts and divide by the number of probes to establish a currentLocus normalization value
    currentLocusNormValue <- sum(as.integer(probeCountList[allCurrentLocusProbes]))/length(allCurrentLocusProbes)
    
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

## This function reads through the SAM table and assigns reads to different lists for downstream analysis
run.sam_read_assignment <- function(samTable, sortedReadNames, currentPresentLoci, multiReferenceSorting=F){
  cat('\nAssigning reads to paired/unpaired/multiRef/multiLocus... 0% ')
  
  ## Initialize list to store reads that match multiple loci
  multiLocusReads <- c()
  
  ## Initialze list to store reads that match absent loci
  #absentLocusReads <- c()
  
  ## Initialize list to store reads that are paired and match a single locus
  singleLocusPairedReads <- c()
  
  ## Initialize list to store reads that are unpaired and match a single locus but multiple references
  singleLocusMultiReferenceReads <- c()
  
  ## Initialize list tto store reads that are unpaired and match a single locu
  singleLocusUnpairedReads <- c()
  
  ## Initialize variables to check on the progress of the next for loop
  i = 0
  max_i = length(sortedReadNames)
  checkAmount = ceiling(max_i/10)
  j=0
  
  ## Iterate through the sorted read names
  for(currentReadName in sortedReadNames[(i+1):length(sortedReadNames)]){
    i = i + 1
    
    ## Will display the percent completion every 10%
    if(i%%checkAmount == 0){
      j = j+10
      if(j != 100){
        cat(paste0(j,'% ', collapse = ''))
      }
    }
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    ## Pull out the loci matched for this read
    currentReadLocusList <- unique(samSubsetTable$locus)
    
    ## Special KIR2Dl5 finagling
    #if('KIR2DL5A' %in% currentReadLocusList | 'KIR2DL5B' %in% currentReadLocusList){
    #  currentReadLocusList <- c(currentReadLocusList, 'KIR2DL5')
    #}
    
    ## Get rid of the extra 2DL5A/B loci names
    #currentReadLocusList <- currentReadLocusList[currentReadLocusList %in% kirLocusList]
    
    ## Only continue if there is only one matched locus
    if(length(currentReadLocusList) > 1){
      multiLocusReads <- c(multiLocusReads, currentReadName)
      next
    }
    
    ## Pull out the current read name alignments that match present loci
    #samSubsetTable <- samSubsetTable[locus %in% currentPresentLoci]
    
    ## Skip this iteration if the samSubsetTable is empty
    #if(nrow(samSubsetTable)==0){
    #  absentLocusReads <- c(absentLocusReads, currentReadName)
    #  next
    #}
    
    ## Interpret the SAM flags
    samFlagTable <- sapply(samSubsetTable$'2', sam_flag_convert)
    
    ## Pull out reads that have proper pair mapping
    #samSubsetTable <- samSubsetTable[unlist(samFlagTable['properPairMapping',]),]
    unpairedSubsetTable <- samSubsetTable[!unlist(samFlagTable['mateUnmapped',]),]
    
    ## If all of the rows of the subset table are removed, then move on to the next read
    if(nrow(unpairedSubsetTable)==0){
      
      if(multiReferenceSorting){
        currentReadRef <- unique(samSubsetTable$reference_name)
        if(length(currentReadRef) > 1){
          singleLocusMultiReferenceReads <- c(singleLocusMultiReferenceReads, currentReadName)
        }
      }
      
      singleLocusUnpairedReads <- c(singleLocusUnpairedReads, currentReadName)
      next
    }
    
    singleLocusPairedReads <- c(singleLocusPairedReads, currentReadName)
    next
  }
  cat(' 100%')
  cat('\n\tSingle locus paired-end reads:',length(singleLocusPairedReads))
  cat('\n\tSingle locus unpaired reads:',length(singleLocusUnpairedReads))
  cat('\n\tSingle locus multi reference reads:',length(singleLocusMultiReferenceReads))
  cat('\n\tMulti locus reads:',length(multiLocusReads))
  
  return(list(singleLocusMultiReferenceReads=singleLocusMultiReferenceReads,
              multiLocusReads=multiLocusReads,
              singleLocusPairedReads=singleLocusPairedReads,
              singleLocusUnpairedReads=singleLocusUnpairedReads))
}

## Assign single locus matched paired-end reads to their mapped postions in the assembly table
old.run.assemble_paired_reads <- function(samTable, singleLocusPairedReads, assembledSeqList){
  
  cat('\nProcessing single locus paired-end reads... 0% ')
  
  ## Initialize variables to check on the progress of the next for loop
  i = 0
  max_i = length(singleLocusPairedReads)
  checkAmount = ceiling(max_i/10)
  j=0
  
  for(currentReadName in singleLocusPairedReads){
    i = i+1
   
     ## Will display the percent completion every 10%
    if(i%%checkAmount == 0){
      j = j+10
      cat(paste0(j,'% ', collapse = ''))
    }
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    ## Pull out the loci matched for this read
    currentReadLocusList <- unique(samSubsetTable$locus)
    
    ## Special KIR2Dl5 finagling
    if('KIR2DL5A' %in% currentReadLocusList | 'KIR2DL5B' %in% currentReadLocusList){
      currentReadLocusList <- c(currentReadLocusList, 'KIR2DL5')
    }
    
    ## Only consider loci that are present in the sample
    currentReadLocusList <- currentReadLocusList[currentReadLocusList %in% currentPresentLoci]
    
    ## Pull out the current read name alignments that match present loci
    samSubsetTable <- samSubsetTable[locus %in% currentReadLocusList]
    
    ## Pull out the first reference allele name 
    refAlleleName <- samSubsetTable$reference_name[[1]]
    
    ## Pull out the current locus
    currentLocus <- kir.allele_resolution(refAlleleName, 0)
    
    if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
      currentLocus <- 'KIR2DL5'
    }
    
    ## Subset the table by the current reference allele
    currentLocusSubsetTable <- samSubsetTable[reference_name == refAlleleName]
    
    ## Check to make sure the expected number of rows are in the table
    if(nrow(currentLocusSubsetTable) != 2){
      stop('currentLocusSubsetTable has a different number of rows than expected.')
    }
    
    ## Initialize a list to store the full index of both reads of the pair
    fullIndexList <- c()
    
    ## Initialize a list to store the full sequence of both reads of the pair
    seqListList <- c()
    
    ## Iterate over both paired reads and figure out their positions and sequence
    for(currentRowIndex in 1:2){
      
      ## Pull the reference allele position from the SAM table
      refAllelePos <- as.integer(currentLocusSubsetTable[currentRowIndex,4][[1]])
      
      ## Pull the current read sequence from the SAM table
      currentReadSeq <- currentLocusSubsetTable[currentRowIndex,10][[1]]
      
      ## Determine the current read length
      currentReadLen <- nchar(currentReadSeq)
      
      ## If there is more than 1 reference allele, then we have a problem
      if(length(refAllelePos) > 1){
        stop('More than 1 reference allele was found in the currentLocusSubsetTable')
      }
      
      ## Establish the starting position, accounting for deletions
      delStartOffset <- which(kirAlleleDFList[[currentLocus]][refAlleleName,] != '.')[refAllelePos]
      
      ## Establish the ending position, accounting for deletions
      delEndOffset <- which(kirAlleleDFList[[currentLocus]][refAlleleName,] != '.')[refAllelePos+currentReadLen-1]
      
      ## Testing if there deletions found in the reference for the current read
      if((delEndOffset-delStartOffset+1) > currentReadLen){
        
        ## Picking out the indices of deletions
        delIndexList <- which(kirAlleleDFList[[currentLocus]][refAlleleName,(delStartOffset):(delEndOffset)] == '.')
        
        ## If the index list is empty, something went wrong
        if(length(delIndexList) == 0){
          stop('Something weird happened with the deletion index.')
        }
        
        ## Split apart the read sequence and add in deletion symbols
        subStringList <- c()
        for(delIndex in delIndexList){
          preString <- substr(currentReadSeq, 1, delIndex-1)
          postString <- substr(currentReadSeq, delIndex, nchar(currentReadSeq))
          
          currentReadSeq <- paste0(preString, '.', postString)
        }
      }
      
      ## Create a full index that cooresponds to the sequence nucleotides and where they go
      fullIndex <- delStartOffset:delEndOffset
      
      ## Add the current index to the list of indices
      fullIndexList <- c(fullIndexList, fullIndex)
      
      ## Turn the sequence string into a list
      seqList <- strsplit(currentReadSeq,'')[[1]]
      
      ## Add the current read nucleotide list to the paired-read nucleotide list
      seqListList <- c(seqListList,seqList)
      
      ## Name the nucleotides by their position
      names(seqListList) <- fullIndexList
    }
    
    ## Cut down the paired-read nucleotide list by unique positions
    sequenceToInsert <- seqListList[unique(names(seqListList))]
    
    ## Create a depth table for each position (positions overlapping in the paired-end reads will have depth = 2)
    depthTable <- table(fullIndexList)
    
    ## If the length of the depth table and the sequnce list are not the same, then we have a problem
    if(length(sequenceToInsert) != length(depthTable)){
      stop('The length of the depth table and the sequence table do not match up')
    }
    
    ## If the name of both are not the same, then something is not running correctly
    if(!all(names(sequenceToInsert) %in% names(depthTable))){
      stop('The names of the depth table and sequence table do not match up')
    }
    
    ## This section compares the current read seqence to the assembly frame and determines which row it should be assigned to
    ## If the full sequence matches the sequence already found in a row, discounting unfilled positions (0's),
    ## then place the sequence in that row, otherwise go to the next row
    rowIndex <- -1
    consensusBool <- FALSE
    while(!consensusBool){
      ## Iterate by 2 row indices each iteration since 1 row is for depth
      rowIndex <- rowIndex + 2
      
      ## Pull out the existing assembled sequence for this region
      assembledSequencePortion <- assembledSeqList[[currentLocus]][rowIndex,names(sequenceToInsert)]
      
      ## Figure out which postions are currently unfilled (0's)
      zeroPositions <- names(assembledSequencePortion)[assembledSequencePortion == 0]
      
      ## Figure out which positions are filled
      filledPositions <- setdiff(names(sequenceToInsert), zeroPositions)
      
      ## If one or more positions are filled, then compare the new sequence to the assembled sequence and see if they match
      if(length(filledPositions) > 0){ 
        
        ## Determine if we break out of the loop based on if the sequences match or not
        consensusBool <- all(assembledSeqList[[currentLocus]][rowIndex,filledPositions] == sequenceToInsert[filledPositions])
      }else{
        
        ## If all current positions are unfilled, then break out of the loop
        consensusBool <- TRUE
      }
      
      ## If we have not broken out of the loop once we have exceeded the size of the assembly table, then throw an error
      if(consensusBool == FALSE & rowIndex > nrow(assembledSeqList[[currentLocus]])){
        consensusBool <- TRUE
        stop('Exceeded all alternative sequence slots for this region')
      }
    }
    
    ## Place the new sequence into the assembly table
    assembledSeqList[[currentLocus]][rowIndex,names(sequenceToInsert)] <- sequenceToInsert
    
    ## Increase the depth of the current positions based on the depth table
    assembledSeqList[[currentLocus]][(rowIndex+1),names(depthTable)] <- as.integer(assembledSeqList[[currentLocus]][(rowIndex+1),names(depthTable)]) + depthTable
    
    ## If we are close to maxing out the assembly table, then stop the program.
    if(rowIndex > 98){
      stop()
    }
  }
  
  return(assembledSeqList)
}

## This function adds single locus paired-end reads to the assembly
run.assemble_paired_reads <- function(samTable, singleLocusPairedReads, assembledNucList, deletionIndexList, kirAlleleDFList){
  cat('\nProcessing single locus paired-end reads... 0% ')
  
  ## Randomize the read name order
  set.seed(001) ## Just to make it reproducible
  randomPairedReads <- sample(singleLocusPairedReads)
  
  ## Check if there are more reads than the threshold and take some out if so
  if(length(randomPairedReads ) > 10000){
    randomPairedReads <- randomPairedReads[1:10000]
  }
  
  ## Initialize variables to check on the progress of the next for loop
  i = 0
  max_i = length(singleLocusPairedReads)
  checkAmount = ceiling(max_i/10)
  j=0
  
  for(currentReadName in randomPairedReads){
    i = i+1
    
    ## Will display the percent completion every 10%
    if(i%%checkAmount == 0){
      j = j+10
      if(j != 100){
        cat(paste0(j,'% ', collapse = ''))
      }
    }
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    ## Pull out the loci matched for this read
    currentReadLocusList <- unique(samSubsetTable$locus)
    
    ## Special KIR2Dl5 finagling
    if('KIR2DL5A' %in% currentReadLocusList | 'KIR2DL5B' %in% currentReadLocusList){
      currentReadLocusList <- c(currentReadLocusList, 'KIR2DL5')
    }
    
    ## Only consider loci that are present in the sample
    currentReadLocusList <- currentReadLocusList[currentReadLocusList %in% currentPresentLoci]
    
    ## Pull out the current read name alignments that match present loci
    samSubsetTable <- samSubsetTable[locus %in% currentReadLocusList]
    
    ## Pull out the first reference allele name 
    refAlleleName <- samSubsetTable$reference_name[[1]]
    
    ## Pull out the current locus
    currentLocus <- kir.allele_resolution(refAlleleName, 0)
    
    if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
      currentLocus <- 'KIR2DL5'
    }
    
    ## Subset the table by the current reference allele
    currentLocusSubsetTable <- samSubsetTable[reference_name == refAlleleName]
    
    ## Check to make sure the expected number of rows are in the table
    if(nrow(currentLocusSubsetTable) != 2){
      next
      #stop('currentLocusSubsetTable has a different number of rows than expected.')
    }
    
    ## Initialize a list to store the full index of both reads of the pair
    fullIndexList <- c()
    
    ## Initialize a list to store the full sequence of both reads of the pair
    seqListList <- c()
    
    ## Iterate over both paired reads and figure out their positions and sequence
    for(currentRowIndex in 1:2){
      
      ## Pull the reference allele position from the SAM table
      refAllelePos <- as.integer(currentLocusSubsetTable[currentRowIndex,4][[1]])
      
      ## Pull the current read sequence from the SAM table
      currentReadSeq <- currentLocusSubsetTable[currentRowIndex,10][[1]]
      
      ## Determine the current read length
      currentReadLen <- nchar(currentReadSeq)
      
      ## If there is more than 1 reference allele, then we have a problem
      if(length(refAllelePos) > 1){
        stop('More than 1 reference allele was found in the currentLocusSubsetTable')
      }
      
      ## Establish the starting position, accounting for deletions
      delStartOffset <- deletionIndexList[[refAlleleName]][refAllelePos]
      
      ## Establish the ending position, accounting for deletions
      delEndOffset <- deletionIndexList[[refAlleleName]][refAllelePos+currentReadLen-1]
      
      ## Testing if there deletions found in the reference for the current read
      if((delEndOffset-delStartOffset+1) > currentReadLen){
        
        ## Picking out the indices of deletions
        delIndexList <- which(kirAlleleDFList[[currentLocus]][refAlleleName,(delStartOffset):(delEndOffset)] == '.')
        
        ## If the index list is empty, something went wrong
        if(length(delIndexList) == 0){
          stop('Something weird happened with the deletion index.')
        }
        
        ## Split apart the read sequence and add in deletion symbols
        subStringList <- c()
        for(delIndex in delIndexList){
          preString <- substr(currentReadSeq, 1, delIndex-1)
          postString <- substr(currentReadSeq, delIndex, nchar(currentReadSeq))
          
          currentReadSeq <- paste0(preString, '.', postString)
        }
      }
      
      ## Create a full index that cooresponds to the sequence nucleotides and where they go
      fullIndex <- delStartOffset:delEndOffset
      
      ## Add the current index to the list of indices
      fullIndexList <- c(fullIndexList, fullIndex)
      
      ## Turn the sequence string into a list
      seqList <- strsplit(currentReadSeq,'')[[1]]
      
      seqListList <- c(seqListList,seqList)
      
      names(seqListList) <- fullIndexList
    }
    
    ## Cut down the paired-read nucleotide list by unique positions
    sequenceToInsert <- seqListList[unique(names(seqListList))]
    
    ## Create a depth table for each position (positions overlapping in the paired-end reads will have depth = 2)
    depthTable <- table(fullIndexList)
    
    ## If the length of the depth table and the sequnce list are not the same, then we have a problem
    if(length(sequenceToInsert) != length(depthTable)){
      stop('The length of the depth table and the sequence table do not match up')
    }
    
    ## If the name of both are not the same, then something is not running correctly
    if(!all(names(sequenceToInsert) %in% names(depthTable))){
      stop('The names of the depth table and sequence table do not match up')
    }
    
    ## Add the depth to the corresponding row in the assembledNucList frame
    for(pos in names(sequenceToInsert)){
      rowNum <- nucListConv[[sequenceToInsert[pos]]]
      assembledNucList[[currentLocus]][rowNum,pos] <- depthTable[pos] + as.integer(assembledNucList[[currentLocus]][rowNum,pos])
    }
  }
  cat('100%')
  return(assembledNucList)
}

## This function adds single locus unpaired reads to the assembly
run.assemble_unpaired_reads <- function(samTable, singleLocusUnpairedReads, assembledNucList, deletionIndexDFList, inverseDeletionIndexDFList, kirAlleleDFList, nucListConv, kirLocusList){
  cat('\nProcessing single locus unpaired reads... 0% ')
  
  ## Randomize the read name order
  set.seed(001) ## Just to make it reproducible
  randomUnpairedReads <- sample(singleLocusUnpairedReads)
  
  ## Check if there are more reads than the threshold and take some out if so
  if(length(randomUnpairedReads ) > 10000){
    randomUnpairedReads <- randomUnpairedReads[1:10000]
  }
  
  ## Initialize variables to check on the progress of the next for loop
  i = 0
  max_i = length(randomUnpairedReads)
  checkAmount = ceiling(max_i/10)
  j=0
  
  for(currentReadName in randomUnpairedReads){
    i = i+1
    
    ## Will display the percent completion every 10%
    if(i%%checkAmount == 0){
      j = j+10
      if(j != 100){
        cat(paste0(j,'% ', collapse = ''))
      }
    }
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    ## Pull out the loci matched for this read
    currentReadLocusList <- unique(samSubsetTable$locus)
    
    ## Special KIR2Dl5 finagling
    if('KIR2DL5A' %in% currentReadLocusList | 'KIR2DL5B' %in% currentReadLocusList){
      currentReadLocusList <- c(currentReadLocusList, 'KIR2DL5')
      currentReadLocusList <- currentReadLocusList[currentReadLocusList %in% kirLocusList]
    }
    
    ## Only consider loci that are present in the sample
    #currentReadLocusList <- currentReadLocusList[currentReadLocusList %in% currentPresentLoci]
    
    ## Pull out the current read name alignments that match present loci
    #samSubsetTable <- samSubsetTable[locus %in% currentReadLocusList]
    
    ## Pull out the first reference allele name 
    #refAlleleName <- samSubsetTable$reference_name[[1]]
    
    ## Pull out the current locus
    #currentLocus <- kir.allele_resolution(refAlleleName, 0)
    
    #if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
    #  currentLocus <- 'KIR2DL5'
    #}
    
    ## Subset the table by the current reference allele
    #currentLocusSubsetTable <- samSubsetTable[reference_name == refAlleleName][1,,drop=F]
    
    ## Add the read to the assembly
    assembledNucList <- build.add_read_to_assembly(assembledNucList,samSubsetTable,delIndexList,inverseDeletionIndexDFList,kirAlleleDFList,nucListConv,currentReadLocusList)
  }
  cat('100%')
  return(assembledNucList)
}

## This function adds multi locus reads to the assembly (comparing them to the determined good allele list to sort out which locus they belong to)
run.assemble_multi_reads <- function(samTable, multiLocusReads, assembledNucList, deletionIndexDFList, inverseDeletionIndexDFList, kirAlleleDFList, nucListConv, allGoodAlleles, kirLocusList){
  
  cat('\nProcessing multi locus reads... 0% ')
  
  ## Randomize the read name order
  set.seed(001) ## Just to make it reproducible
  randomMultiReads <- sample(multiLocusReads)
  
  ## Check if there are more reads than the threshold and take some out if so
  if(length(randomMultiReads) > 10000){
    randomMultiReads <- randomMultiReads[1:10000]
  }
  
  ## Initialize variables to check on the progress of the next for loop
  i = 0
  max_i = length(randomMultiReads)
  checkAmount = ceiling(max_i/10)
  j=0
  
  for(currentReadName in randomMultiReads){
    i = i+1
    
    ## Will display the percent completion every 10%
    if(i%%checkAmount == 0){
      j = j+10
      if(j != 100){
        cat(paste0(j,'% ', collapse = ''))
      }
    }
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    ## Pull out the lines of the subset sam table that have references that match the determined good alleles
    goodAlleleSubsetTable <- samSubsetTable[reference_name %in% allGoodAlleles]
    
    ## Pull out the loci matched for this read
    currentReadLocusList <- unique(goodAlleleSubsetTable$locus)
    
    ## Special KIR2Dl5 finagling
    if('KIR2DL5A' %in% currentReadLocusList | 'KIR2DL5B' %in% currentReadLocusList){
      currentReadLocusList <- c(currentReadLocusList, 'KIR2DL5')
      currentReadLocusList <- currentReadLocusList[currentReadLocusList %in% kirLocusList]
    }
    
    if(length(currentReadLocusList) != 1){
      next
    }else{
      ## Pull out the first reference allele name 
      refAlleleName <- goodAlleleSubsetTable$reference_name[[1]]
      
      ## Pull out the current locus
      currentLocus <- kir.allele_resolution(refAlleleName, 0)
      
      if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
        currentLocus <- 'KIR2DL5'
      }
      
      ## Subset the table by the current reference allele
      currentLocusSubsetTable <- goodAlleleSubsetTable[reference_name == refAlleleName][1,,drop=F]
      
      assembledNucList <- build.add_read_to_assembly(assembledNucList,currentLocusSubsetTable,deletionIndexDFList,inverseDeletionIndexDFList,kirAlleleDFList,nucListConv,currentLocus)
    }
  }
  
  cat('100%')
  return(assembledNucList)
}

## This function runs the count data through a trained random forest model for predicting copy number
run.predict_copy <- function(locusRatioDF, locusCountDF, copyNumberDF, goodRows, resultsDirectory, rfAllPathList){
  
  ### Prepare data for merging
  colnames(locusCountDF) <- paste0('locusCount_', colnames(locusCountDF))
  colnames(locusRatioDF) <- paste0('locusNorm_', colnames(locusRatioDF))
  ###
  
  locusCountDF <- locusCountDF[goodRows,]
  locusRatioDF <- locusRatioDF[goodRows,]
  
  ### Merge the data used for prediction
  mergedData <- merge(locusCountDF[,'locusCount_KIR3DL3',drop=F], locusRatioDF, by=0, all=T)
  rownames(mergedData) <- mergedData$Row.names
  mergedData <- mergedData[,2:ncol(mergedData)]
  ###
  
  ### Predict
  for(currentLocus in names(rfAllPathList)){
    load(rfAllPathList[[currentLocus]])
    predData <- predict(model1, mergedData, type='class')
    copyNumberDF[names(predData),currentLocus] <- as.numeric(levels(predData))[predData]
  }
  copyNumberDF[,'KIR3DL3'] <- 2
  ###
  
  return(copyNumberDF)
}

## This function prompts for user input for determining copy number thresholds
readRatioPrompt <- function(copyNumber){ 
  n <- readline(prompt=paste0('Enter the lowest ratio for copy number ', copyNumber, ': '))
  n <- suppressWarnings(as.numeric(n))
  if(is.na(n)){
    cat('Only numeric input is accepted. Please try again.\n')
    return(readRatio(copyNumber))
  }
  return(n)
}

## This function prompts for user input for determining the maximum copy number for a locus
copyNumberPrompt <- function(currentLocus,minCopy,maxCopy){ 
  n <- readline(prompt=paste0('Please input the max copy number for ', currentLocus, ': '))
  n <- suppressWarnings(as.integer(n))
  if(is.na(n) | n<minCopy | n>maxCopy){
    cat(paste0('Only numeric input from ',minCopy,'-',maxCopy,' is allowed. Please try again.\n'))
    return(copyNumberPrompt(currentLocus,minCopy,maxCopy))
  }
  return(n)
}

## This function builds a list of deletion indices for each kir allele
build.deletion_index_list <- function(kirAlleleList, kirAlleleDFList){
  cat('\nBuilding up a deletion index for faster lookup during read assignment.')
  
  ## Initialize a list for storing deletion position conversions
  deletionIndexDFList <- list()
  for(currentAllele in kirAlleleList){
    ## Pull out the current locus
    currentLocus <- kir.allele_resolution(currentAllele, 0)
    
    if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
      currentLocus <- 'KIR2DL5'
    }
    
    deletionIndexDFList[[currentAllele]] <- which(kirAlleleDFList[[currentLocus]][currentAllele,] != '.')
  }
  return(deletionIndexDFList)
}

## This function builds the inverse of the deletion indices for each kir allele
build.inverse_deletion_index_list <- function(kirAlleleList, kirAlleleDFList){
  cat('\nBuilding up an inverse deletion index for faster lookup during read assignment.')
  
  ## Initialize a list for storing deletion position conversions
  deletionIndexDFList <- list()
  for(currentAllele in kirAlleleList){
    ## Pull out the current locus
    currentLocus <- kir.allele_resolution(currentAllele, 0)
    
    if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
      currentLocus <- 'KIR2DL5'
    }
    
    deletionIndexDFList[[currentAllele]] <- which(kirAlleleDFList[[currentLocus]][currentAllele,] == '.')
  }
  return(deletionIndexDFList)
}

## This function initializes a list of dataframes for storing the assembled nucleotides for each kir locus
build.assembled_nuc_list <- function(kirLocusList, kirAlleleDFList){
  assembledNucList <- list()
  for(currentLocus in kirLocusList){
    refSeqDF <- data.frame(matrix(0, 5, ncol(kirAlleleDFList[[currentLocus]])), check.names=F, stringsAsFactors=F)
    assembledNucList[[currentLocus]] <- as.data.table(refSeqDF)
  }
  return(assembledNucList)
}

## This function initializes a list of dataframes for storing the assembled nucleotides for each kir locus
build.promoter_nuc_list <- function(kirLocusList, kirExonCoords){
  assembledNucList <- list()
  for(currentLocus in kirLocusList){
    
    ## Pull out up through the first intron, and then the last intron to the end
    posVect <- unlist(kirExonCoords[[currentLocus]][c('5UTR','E1','I1', tail(names(kirExonCoords[[currentLocus]]), 3))], use.names=F)
    
    refSeqDF <- data.frame(matrix(0, 5, length(posVect)), check.names=F, stringsAsFactors=F)
    colnames(refSeqDF) <- as.character(posVect)
    assembledNucList[[currentLocus]] <- refSeqDF
  }
  return(assembledNucList)
}

## This function adds single read information to the assembled nuc list
build.add_read_to_assembly <- function(assembledNucList, samSubsetTable, delIndexDFList, inverseDeletionIndexDFList, kirAlleleDFList, nucListConv, currentLocus){

  ## Pull out the first reference allele name 
  refAlleleName <- samSubsetTable$reference_name[[1]]

  ## Pull out the the first reference allele position
  refAllelePos <- as.integer(samSubsetTable[1,4][[1]])

  ## Pull out the current read sequence
  currentReadSeq <- samSubsetTable[1,10][[1]]

  ## Store the length of the sequence
  currentReadLen <- nchar(currentReadSeq)

  ## Establish the starting position, accounting for deletions
  delStartOffset <- deletionIndexDFList[[refAlleleName]][refAllelePos]
  
  ## Establish the ending position, accounting for deletions
  delEndOffset <- deletionIndexDFList[[refAlleleName]][refAllelePos+currentReadLen-1]

  ## Testing if there deletions found in the reference for the current read
  if((delEndOffset-delStartOffset+1) > currentReadLen){

    delIndexList <- which(delStartOffset:delEndOffset %in% inverseDeletionIndexDFList[[refAlleleName]])
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    subStringList <- c()
    for(delIndex in delIndexList){
      preString <- substr(currentReadSeq, 1, delIndex-1)
      postString <- substr(currentReadSeq, delIndex, nchar(currentReadSeq))
      
      currentReadSeq <- paste0(preString, '.', postString)
    }
  }

  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- delStartOffset:delEndOffset
  
  ## Turn the sequence string into a list
  seqList <- strsplit(currentReadSeq,'')[[1]]
  
  names(seqList) <- fullIndex
  
  sequenceToInsert <- seqList[unique(names(seqList))]
  depthTable <- table(fullIndex)
  
  if(length(sequenceToInsert) != length(depthTable)){
    stop('The length of the depth table and the sequence table do not match up')
  }
  
  if(!all(names(sequenceToInsert) %in% names(depthTable))){
    stop('The names of the depth table and sequence table do not match up')
  }
  
  for(pos in names(sequenceToInsert)){
    rowNum <- nucListConv[[sequenceToInsert[pos]]]
    assembledNucList[[currentLocus]][rowNum,pos] <- as.integer(assembledNucList[[currentLocus]][rowNum,pos]) + 1
  }
  return(assembledNucList)
}

## This function takes in the assembly and determines which positions pass the depthThreshold
build.good_coord_bad_coord_list <- function(kirLocusList,assembledNucList,depthThreshold){
  
  ## Figure out which positions are above or below the set depthThreshold for each locus
  cat('\nDetermining what proportion of each locus is above the depth threshold.')
  allLocusGoodCoords <- list()
  allLocusBadCoords <- list()
  for(currentLocus in kirLocusList){
    
    ## Pull out the assembly for the current locus
    currentLocusAssembly <- assembledNucList[[currentLocus]]
    
    ## Mark above depth coords as positions that have depth above or equal to the depthThreshold
    #aboveDepthCoords <- colnames(currentLocusAssembly)[apply(currentLocusAssembly, 2, sum) >= depthThreshold]
    aboveDepthCoords <- colnames(currentLocusAssembly)[apply(currentLocusAssembly, 2, function(x){
      any(x>=depthThreshold)
    })]
    #aboveDepthCoords <- intersect(aboveDepthCoords, colnames(kirAlleleDFList[[currentLocus]]))
    
    ## Do the inverse for below depth coords
    #belowDepthCoords <- colnames(currentLocusAssembly)[apply(currentLocusAssembly, 2, sum) < depthThreshold]
    belowDepthCoords <- colnames(currentLocusAssembly)[apply(currentLocusAssembly, 2, function(x){
      all(x<depthThreshold)
    })]
    #belowDepthCoords <- intersect(belowDepthCoords, colnames(kirAlleleDFList[[currentLocus]]))
    
    ## Put the above depth coords into a list
    allLocusGoodCoords[[currentLocus]] <- aboveDepthCoords
    
    ## Put the below depth coords into a list
    allLocusBadCoords[[currentLocus]] <- belowDepthCoords
    filledRatio <- round(length(aboveDepthCoords)/ncol(currentLocusAssembly)*100,1)
    
    ## Print statements!
    cat(paste0('\n\n\t',currentLocus,': ',filledRatio, '% filled'))
    
    ## Figuring out the fill percentage by region
    for(currentExon in names(kirExonCoords[[currentLocus]])){
      currentExonCoords <- kirExonCoords[[currentLocus]][[currentExon]]
      totalExonLength <- length(currentExonCoords)
      fillExonLength <- sum(aboveDepthCoords %in% currentExonCoords)
      fillExonRatio <- round((fillExonLength/totalExonLength)*100,1)
      
      ## Double tab for exons, triple for introns and UTR's
      if(grepl(pattern='E', x = currentExon, fixed = T)){
        cat(paste0('\n\t\t[',currentExon,']',fillExonRatio,'% '))
      }else{
        cat(paste0('\n\t\t\t[',currentExon,']',fillExonRatio,'% '))
      }
      
    }
    
    #cat(paste0('\n',currentLocus,': ',filledRatio, '% filled'))
  }
  return(list('badCoords'=allLocusBadCoords,'goodCoords'=allLocusGoodCoords))
}

## This function takes in the assembly and converts depth values to nucleotides
build.nuc_assembly_frame <- function(assembledNucList, kirLocusList, goodBadCoordList, hetRatio, depthThreshold){
  
  allLocusGoodCoords <- goodBadCoordList$goodCoords
  
  allLocusAssemblyList <- list()
  
  cat('\n\tProcessing:')
  
  for(currentLocus in names(assembledNucList)){
    #currentLocus <- 'KIR3DL2'
    cat('',currentLocus)
    
    ## Pull out the assembly for the current locus
    currentLocusAssembly <- as.data.frame(assembledNucList[[currentLocus]])
    
    #currentLocusColNames <- colnames(kirAlleleDFList[[currentLocus]])
    #currentLocusAssembly <- currentLocusAssembly[,currentLocusColNames]
    
    ## Initialize a dataframe for storing the assembled nucleotides
    currentLocusAssemblyNucDF <- data.frame(matrix('',5,ncol(currentLocusAssembly)),stringsAsFactors=F,check.names=F)
    colnames(currentLocusAssemblyNucDF) <- colnames(currentLocusAssembly)
    
    ## For each coordinate that passes the good call threshold, transfer the nucleotide value(s) to the assembled nucleotide dataframe
    for(assemblyPos in allLocusGoodCoords[[currentLocus]]){
      #cat('\n',assemblyPos)
      
      ## Pull out the depth for the current position
      currentPosReadCount <- currentLocusAssembly[,assemblyPos,drop=F]
      
      ## Detemine which nucs pass the hetRatio threshold
      goodSnps <- (currentPosReadCount/sum(currentPosReadCount))>hetRatio & currentPosReadCount >= depthThreshold
      
      ### Fill in the nuc value
      currentLocusAssemblyNucDF[1:sum(goodSnps),assemblyPos] <- rownames(currentPosReadCount)[goodSnps]
    }
    
    ## If there are any positions with 3 SNP calls, stop the program for investigation
    if(sum(currentLocusAssemblyNucDF[3,] != '')>0){
      #stop('Found a position with more than 3 SNP calls')
    }
    
    ## Cut the assembly dataframe down to 2 rows
    currentLocusAssemblyNucDF <- currentLocusAssemblyNucDF[1:2,]
    
    ## Add the current nuc assembly dataframe to the output list of dataframes 
    allLocusAssemblyList[[currentLocus]] <- currentLocusAssemblyNucDF
  }
  return(allLocusAssemblyList)
}

## Use the determined good alleles to finish off filling in the allLocusNucAssembly
build.fill_in_remaining_positions_with_good_alleles <- function(currentPresentLoci, allGoodAlleles, kirAlleleDFList,goodBadCoordList, allLocusNucAssembly,copyNumberList){
  
  cat('\nUsing the determined good alleles to fill in any remaining unfilled positions.')
  
  for(currentLocus in currentPresentLoci){
    
    if(currentLocus == 'KIR2DL5'){
      currentLocusGoodAlleles <- allGoodAlleles[tstrsplit(allGoodAlleles,'*',fixed=T)[[1]] %in% c('KIR2DL5A','KIR2DL5B')]
    }else{
      currentLocusGoodAlleles <- allGoodAlleles[tstrsplit(allGoodAlleles,'*',fixed=T)[[1]] == currentLocus]
    }
    
    if(length(currentLocusGoodAlleles)==0){
      stop('No good alleles found for this locus')
    }
    
    copyNumber <- copyNumberList[[currentLocus]]
    
    goodAlleleSequenceDF <- kirAlleleDFList[[currentLocus]][currentLocusGoodAlleles,goodBadCoordList$badCoords[[currentLocus]]]
    
    goodAlleleNAllDF <- goodAlleleSequenceDF[,apply(goodAlleleSequenceDF,2,function(x)all(x=='N')),drop=F]
    
    if(ncol(goodAlleleNAllDF) != 0){
      cat('\n',colnames(goodAlleleNAllDF))
      stop('good alleles all have N for at least one position')
    }
    
    for(colPos in colnames(goodAlleleSequenceDF)){
      uniqueNucs <- unique(goodAlleleSequenceDF[,colPos])
      uniqueNucs <- uniqueNucs[uniqueNucs != 'N']
      
      if(length(uniqueNucs) > copyNumber){
        nucTable <- table(goodAlleleSequenceDF[,colPos])
        uniqueNucs <- names(nucTable)[order(nucTable, decreasing=T)][1:copyNumber]
      }
      
      allLocusNucAssembly[[currentLocus]][1:length(uniqueNucs),colPos] <- uniqueNucs
    }
    
    currentLocusFilledSum <- sum(allLocusNucAssembly[[currentLocus]][1,] != '')
    currentLocusAllSum <- ncol(allLocusNucAssembly[[currentLocus]])
    
    if(currentLocusFilledSum != currentLocusAllSum){
      stop('Not all positions were filled in!')
    }
  }
  
  return(allLocusNucAssembly)
}


new.run.assemble_unpaired_reads <- function(currentReadName, samTable, deletionIndexDFList, inverseDeletionIndexDFList){
  samSubsetTable <- samTable[read_name == currentReadName]
  currentReadLocus <- samSubsetTable$locus[[1]]
  
  return(new.build.add_read_to_assembly(samSubsetTable, deletionIndexDFList, inverseDeletionIndexDFList, currentReadLocus))
}

new.build.add_read_to_assembly <- function(samSubsetTable, deletionIndexDFList, inverseDeletionIndexDFList, currentReadLocus){
  ## Pull out the first reference allele name 
  refAlleleName <- samSubsetTable$reference_name[[1]]
  
  ## Pull out the the first reference allele position
  #refAllelePos <- as.integer(samSubsetTable[1,4][[1]])
  #refAllelePos <- as.integer(samSubsetTable$ref_pos[[1]])
  
  ## Pull out the current read sequence
  #currentReadSeq <- samSubsetTable[1,10][[1]]
  currentReadSeq <- samSubsetTable$read_seq[[1]]
  
  ## Store the length of the sequence
  #currentReadLen <- nchar(currentReadSeq)
  currentReadLen <- as.integer(samSubsetTable$readLen[[1]])
  
  ## Establish the starting position, accounting for deletions
  #delStartOffset <- deletionIndexDFList[[refAlleleName]][refAllelePos]
  #samTable[read_name == samSubsetTable$read_name[1], startPos := delStartOffset]
  delStartOffset <- as.integer(samSubsetTable$startPos[[1]])
  
  ## Establish the ending position, accounting for deletions
  #delEndOffset <- deletionIndexDFList[[refAlleleName]][refAllelePos+currentReadLen-1]
  #samTable[read_name == samSubsetTable$read_name[1], endPos := delEndOffset]
  delEndOffset <- as.integer(samSubsetTable$endPos[[1]])
  
  ## Testing if there deletions found in the reference for the current read
  if((delEndOffset-delStartOffset+1) > currentReadLen){
    
    delIndexList <- which(delStartOffset:delEndOffset %in% inverseDeletionIndexDFList[[refAlleleName]])
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    subStringList <- c()
    for(delIndex in delIndexList){
      preString <- substr(currentReadSeq, 1, delIndex-1)
      postString <- substr(currentReadSeq, delIndex, nchar(currentReadSeq))
      
      currentReadSeq <- paste0(preString, '.', postString)
    }
  }
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- delStartOffset:delEndOffset
  
  ## Turn the sequence string into a list
  seqList <- strsplit(currentReadSeq,'')[[1]]
  
  names(seqList) <- fullIndex
  
  #sequenceToInsert <- seqList[unique(names(seqList))]
  #depthTable <- table(fullIndex)
  seqListList <- list(seqList)
  return(setNames(seqListList, currentReadLocus))
}

new.run.assemble_paired_reads <- function(currentReadName, samTable, deletionIndexDFList, inverseDeletionIndexDFList){
  samSubsetTable <- samTable[read_name == currentReadName]
  currentReadLocus <- samSubsetTable$locus[[1]]
  
  samFlagTable <- sapply(samSubsetTable$'2', sam_flag_convert)
  
  samSubsetTableFirst <- samSubsetTable[unlist(samFlagTable['firstInPair',]),]
  samSubsetTableSecond <- samSubsetTable[unlist(samFlagTable['secondInPair',]),]
  
  return(sapply(list(samSubsetTableFirst, samSubsetTableSecond), new.build.add_read_to_assembly, deletionIndexDFList, inverseDeletionIndexDFList, currentReadLocus))
}

run.condense_assembly_list <- function(assemblyReadList, assembledNucList, nucListConv, kirAlleleDFList){
  
  ## Initialize a deep list for counting the depth of each nuc at each position for each locus
  allLocusPosCountList <- list()
  
  ## Unlist the read list, each read should be named by locus
  allUnpairedByLocus <- unlist(assemblyReadList, recursive=F)
  
  cat('\n\tDetermining nucleotide depth by position.')
  
  cat('\n\t\tProcessing:')
  ## Iterate through each locus matched in the read list
  for(currentLocus in unique(names(allUnpairedByLocus))){
    cat('',currentLocus)
    
    ## Set the current locus element for the output list
    allLocusPosCountList[[currentLocus]] <- list()
    
    ## Pull out and unlist all of the reads that match the current locus (creates a single named vector)
    allUnpairedCurrentLocus <- unlist(allUnpairedByLocus[names(allUnpairedByLocus) == currentLocus])
    
    ## Split the names of the current locus vector to only include the position (take out the locus name)
    names(allUnpairedCurrentLocus) <- tstrsplit(names(allUnpairedCurrentLocus),'.',fixed=T)[[2]]
    
    ## Fix for reads that extend over the reference end coordinate for the current locus
    ## Define the last position for this locus
    lastPos <- as.integer(tail(colnames(kirAlleleDFList[[currentLocus]]),1))
    
    ## Figure out which read positions are <= the last position
    goodNames <- names(allUnpairedCurrentLocus)[as.integer(names(allUnpairedCurrentLocus)) <= lastPos]
    
    ## Subset the read positions to only include ones <= the last position
    allUnpairedCurrentLocus <- allUnpairedCurrentLocus[names(allUnpairedCurrentLocus) %in% goodNames]
    
    ## Iterate over each possible nucleotide (so max 5 iterations for A, C, T, G, .)
    for(currentNuc in unique(allUnpairedCurrentLocus)){
      
      ## Initialize a list for storing the read positions that match the current nucleotide
      allLocusPosCountList[[currentLocus]][[currentNuc]] <- list()
      
      ## Pull out the read positions that match the current nucleotide
      currentNucPosVect <- names(allUnpairedCurrentLocus)[allUnpairedCurrentLocus == currentNuc]
      
      ## Create a table that shows the number of matches for the current nuc at each unique position (get the depth)
      currentNucPosTable <- table(currentNucPosVect)
      
      ## For each unique nuc depth value, store the positions that match the nuc and depth
      for(i in unique(currentNucPosTable)){
        allLocusPosCountList[[currentLocus]][[currentNuc]][[as.character(i)]] <- names(currentNucPosTable)[currentNucPosTable == i]
      }
    }
  }
  cat('\n\tFinished.')
  cat('\n\tFilling out the nucleotide data.table.')
  ## Now we want to use the allLocusPosCountList to build a nucleotide dataframe for each locus
  cat('\n\t\tProcessing:')
  for(currentLocus in names(allLocusPosCountList)){
    cat('',currentLocus)
    
    for(currentNuc in names(allLocusPosCountList[[currentLocus]])){
      rowInd <- unlist(nucListConv[currentNuc],use.names=F)
      for(depthValue in names(allLocusPosCountList[[currentLocus]][[currentNuc]])){
        
        columnList <- allLocusPosCountList[[currentLocus]][[currentNuc]][[depthValue]]
        depthInt <- as.integer(depthValue) + as.integer(unlist(assembledNucList[[currentLocus]][rowInd, (columnList), with=F], use.names=F))
        
        assembledNucList[[currentLocus]][rowInd, columnList := as.list(depthInt), with=F]
      }
    }
  }
  cat('\n\tFinished.')
  return(assembledNucList)
}


new.run.find_allele_matches_for_assembly <- function(allLocusNucAssembly, currentPresentLoci, kirAlleleDFList, copyNumberList){
  allGoodAlleleList <- list()
  for(currentLocus in currentPresentLoci){
    cat('\n\n',currentLocus)
    copyNumber <- copyNumberList[[currentLocus]]
    currentLocusAssemblyNucDF <- allLocusNucAssembly[[currentLocus]]
    
    ## Pull out the heterozygous coordinates
    currentLocusHetCoords <- colnames(currentLocusAssemblyNucDF)[currentLocusAssemblyNucDF[2,] != '']
    
    ## Pull out the homozygous coordinates
    currentLocusHomCoords <- setdiff(colnames(currentLocusAssemblyNucDF)[currentLocusAssemblyNucDF[1,] != ''], currentLocusHetCoords)
    
    ## Create a dataframe of the het coords
    currentLocusHetSnpDF <- currentLocusAssemblyNucDF[,currentLocusHetCoords,drop=F]
    
    ## Create a dataframe of the hom coords
    currentLocusHomSnpDF <- currentLocusAssemblyNucDF[,currentLocusHomCoords,drop=F]
    
    cat('\nNumber of het positions:',length(currentLocusHetCoords))
    cat('\nNumber of hom positions:',length(currentLocusHomCoords))
    
    ## Pull out the sequence information for all known allele for this locus
    allCurrentLocusAlleles <- kirAlleleDFList[[currentLocus]]
    
    ## Generate a unique position frame to analyze the homozygous positions (or else it takes forever)
    allCurrentLocusAllelesUniquePos <- make_unique_pos_frame(allCurrentLocusAlleles)
    
    ## Initializing a list of possible alleles / allele combinations (for heterozygous loci)
    possibleAlleles <- rownames(kirAlleleDFList[[currentLocus]])
    
    possibleAlleleDF <- data.frame(matrix(0,length(possibleAlleles),2),stringsAsFactors=F,row.names=possibleAlleles,check.names=F)
    colnames(possibleAlleleDF) <- c('hetDistance','homDistance')
    
    currentLocusHomPosVect <- intersect(colnames(currentLocusHomSnpDF), colnames(allCurrentLocusAlleles))
    
    currentLocusHomSnpVect <- apply(currentLocusHomSnpDF[1,currentLocusHomPosVect,drop=F],2,function(x){names(nucListConv)[as.integer(x)]})
    
    ## Determining the homozygous position distance (calculated for all present loci)
    for(singleAllele in possibleAlleles){
      allCurrentLocusAllelesHomSnpVect <- unlist(allCurrentLocusAlleles[singleAllele,currentLocusHomPosVect])
      Nsum <- 0# sum(allCurrentLocusAlleles[singleAllele,currentLocusHomPosVect] == 'N')
      possibleAlleleDF[singleAllele,'homDistance'] <- sum(currentLocusHomSnpVect != allCurrentLocusAllelesHomSnpVect) - Nsum
    }
    
    ## Determining the het position distance (only calculated for copy 2 loci)
    cat('\nMin distance:',min(possibleAlleleDF$homDistance))
    possibleAlleleHits <- rownames(possibleAlleleDF)[possibleAlleleDF$homDistance == min(possibleAlleleDF$homDistance)]
    cat('\nNumber of alleles matching min distance:',length(possibleAlleleHits))
    allGoodAlleleList[[currentLocus]] <- possibleAlleleHits
  }
  return(allGoodAlleleList)
}

## This function compares the nuc assembly to known alleles to try and find the closest allele matches for each locus
run.find_allele_matches_for_assembly <- function(allLocusNucAssembly,currentPresentLoci,kirAlleleDFList,copyNumberList){
  
  cat('\nUsing the assembly to find the closest matching alleles for each locus')
  
  ## Initializing lists for storing the mismatch flexibility for each locus for both heterozygous and homozygous position allele matching
  locusMinHetDistanceAddition <- list('KIR3DL3' = 4,'KIR2DS2'=4, 'KIR2DL2'=4,'KIR2DL3' = 4,'KIR2DL1' = 2,'KIR2DL4'=4,'KIR3DL1'=4,'KIR3DL2'=4,'KIR2DP1'=3,'KIR3DS1'=2,'KIR2DL5'=4,'KIR2DS5'=4,
                                      'KIR2DS1' = 4, 'KIR2DS4'=4,'KIR3DP1'=4,'KIR2DS3'=4)
  
  locusMinHomDistanceAddition <- list('KIR3DL3' = 4,'KIR2DS2'=4, 'KIR2DL2'=4, 'KIR2DL3' = 4,'KIR2DL1' = 4,'KIR2DL4'=4,'KIR3DL1'=4,'KIR3DL2'=4,'KIR2DP1'=3,'KIR3DS1'=4,
                                      'KIR2DL5'=99,'KIR2DS5'=9,'KIR2DS1'=1,'KIR2DS4'=4,'KIR3DP1'=4, 'KIR2DS3'=4)
  
  allGoodAlleleVector <- c() 
  for(currentLocus in currentPresentLoci){
    cat('\n\n',currentLocus)
    copyNumber <- copyNumberList[[currentLocus]]
    currentLocusAssemblyNucDF <- allLocusNucAssembly[[currentLocus]]
    
    ## Pull out the heterozygous coordinates
    currentLocusHetCoords <- colnames(currentLocusAssemblyNucDF)[currentLocusAssemblyNucDF[2,] != '']
    
    ## Pull out the homozygous coordinates
    currentLocusHomCoords <- setdiff(colnames(currentLocusAssemblyNucDF)[currentLocusAssemblyNucDF[1,] != ''], currentLocusHetCoords)
    
    ## Create a dataframe of the het coords
    currentLocusHetSnpDF <- currentLocusAssemblyNucDF[,currentLocusHetCoords,drop=F]
    
    ## Create a dataframe of the hom coords
    currentLocusHomSnpDF <- currentLocusAssemblyNucDF[,currentLocusHomCoords,drop=F]
    
    cat('\nNumber of het positions:',length(currentLocusHetCoords))
    cat('\nNumber of hom positions:',length(currentLocusHomCoords))
    
    ## Initializing a list of possible alleles / allele combinations (for heterozygous loci)
    possibleAlleles <- rownames(kirAlleleDFList[[currentLocus]])
    #possibleAlleles <- rownames(kirAlleleDFList[[currentLocus]])[apply(kirAlleleDFList[[currentLocus]]=='N', 1, sum) < 4000]
    
    
    #for(currentLocus in kirLocusList){
    #  cat('\n',currentLocus)
    #  print(table(apply(kirAlleleDFList[[currentLocus]]=='N', 1, sum)))
    #}
    
    ## If there are more than 1 copies of this locus, then generate all possible 2-allele combinations
    if(copyNumber>1){
      
      ## Generate all 2-allele combinations
      possibleAlleleComboDF <- combinations(length(possibleAlleles),2,possibleAlleles)
      
      ## Initialize a list for storing the combinations once they have been merged
      possibleAlleleList <- c()
      
      ## For each allele combination, merged the two columns into a single string
      for(rowNum in 1:nrow(possibleAlleleComboDF)){
        possibleAlleleList <- c(possibleAlleleList, paste0(possibleAlleleComboDF[rowNum,],collapse='+'))
      }
    }else{
      
      ## If there is only 1 copy of this locus, then set the current locus allele list as all possible alleles
      possibleAlleleList <- possibleAlleles
    }
    
    ## Initialize a dataframe of the possible alleles, this will be used to judge each possible alleles distance to the assembly
    possibleAlleleDF <- data.frame(matrix(0,length(possibleAlleleList),2),stringsAsFactors=F,row.names=possibleAlleleList,check.names=F)
    colnames(possibleAlleleDF) <- c('hetDistance','homDistance')
    
    ## Pull out the sequence information for all known allele for this locus
    allCurrentLocusAlleles <- kirAlleleDFList[[currentLocus]]
    
    ## Generate a unique position frame to analyze the homozygous positions (or else it takes forever)
    allCurrentLocusAllelesUniquePos <- make_unique_pos_frame(allCurrentLocusAlleles)
    
    ## intersect the coordinates between the newly generated unique position frame and the homozygous coords from the assembly
    uniquePosHetCoords <- intersect(colnames(allCurrentLocusAllelesUniquePos),currentLocusHetCoords)
    
    cat('\nComparing',length(uniquePosHetCoords),'het positions from the assembly to',nrow(possibleAlleleDF),'known alleles/allele combinations')
    for(alleleCombo in rownames(possibleAlleleDF)){
      
      if(copyNumber>1){
        alleleVector <- strsplit(alleleCombo,'+',fixed=T)[[1]]
      }else{
        alleleVector <- alleleCombo
      }
      
      for(hetPos in uniquePosHetCoords){
        possibleSnpVect <- allCurrentLocusAlleles[alleleVector,hetPos]
        assembledSnpVect <- currentLocusHetSnpDF[,hetPos]
        
        possibleAlleleDF[alleleCombo,'hetDistance'] <- as.integer(possibleAlleleDF[alleleCombo,'hetDistance']) + sum(!(assembledSnpVect %in% possibleSnpVect))
      }
    }
    
    ## Keep all allele combos that score within 5 of the best score (lower is better)
    possibleAlleleList <- rownames(possibleAlleleDF)[possibleAlleleDF$hetDistance <= min(possibleAlleleDF$hetDistance)+locusMinHetDistanceAddition[[currentLocus]]]
    
    ## Subset the possible allele DF by these new allele combos
    possibleAlleleDF <- possibleAlleleDF[possibleAlleleList,]
    
    ## intersect the coordinates between the newly generated unique position frame and the homozygous coords from the assembly
    uniquePosHomCoords <- intersect(colnames(allCurrentLocusAllelesUniquePos),currentLocusHomCoords)
    
    cat('\nComparing',length(uniquePosHomCoords),'hom positions from the assembly to',nrow(possibleAlleleDF),'known alleles/allele combinations')
    for(alleleCombo in rownames(possibleAlleleDF)){
      
      if(copyNumber > 1){
        alleleVector <- strsplit(alleleCombo,'+',fixed=T)[[1]]
      }else{
        alleleVector <- alleleCombo
      }
      
      #alleleComboHetCoords <- colnames(allCurrentLocusAlleles)[apply(allCurrentLocusAlleles[alleleVector,],2,num_unique_nuc)>1]
      if(copyNumber>1){
        possibleAlleleDF[alleleCombo,2] <- sum(apply(allCurrentLocusAllelesUniquePos[alleleVector,uniquePosHomCoords] != currentLocusAssemblyNucDF[c(1,1),uniquePosHomCoords],2,sum))
      }else{
        possibleAlleleDF[alleleCombo,2] <- sum(apply(allCurrentLocusAllelesUniquePos[alleleVector,uniquePosHomCoords] != currentLocusAssemblyNucDF[1,uniquePosHomCoords],2,sum))
      }
    }
    
    ## Keep all allele combos that score within 5 of the best score (lower is better)
    possibleAlleleList <- rownames(possibleAlleleDF)[possibleAlleleDF$homDistance <= min(possibleAlleleDF$homDistance)+locusMinHomDistanceAddition[[currentLocus]]]
    
    ## Subset the possible allele DF by these new allele combos
    possibleAlleleDF <- possibleAlleleDF[possibleAlleleList,]
    
    ## Pull out the remaining alleles and separate them if needed
    if(copyNumber > 1){
      currentLocusGoodAlleles <- unique(unlist(sapply(rownames(possibleAlleleDF),tstrsplit,'+',fixed=T)))
    }else{
      currentLocusGoodAlleles <- rownames(possibleAlleleDF)
    }
    
    ## Add the good alleles to the output list
    allGoodAlleleVector <- c(allGoodAlleleVector, currentLocusGoodAlleles)
  }
  
  return(allGoodAlleleVector)
}

## This function writes the final assembly to a fasta file
run.convert_assembly_to_fasta <- function(currentPresentLoci, copyNumberList, allLocusNucAssembly, currentSample, assembledReferenceDirectory){
  cat('\nWriting the assembly to a fasta file.')
  
  ## Setting all of the important file paths
  currentSample$assembledRefPath <- file.path(assembledReferenceDirectory,paste0(currentSample$name,'assembled.fasta'))
  currentSample$assembledRefIndex <- file.path(assembledReferenceDirectory,paste0(currentSample$name,'assembled'))
  currentSample$assembledBedPath <- file.path(assembledReferenceDirectory,paste0(currentSample$name,'assembled.bed'))
  currentSample$assembledDelIndex <- file.path(assembledReferenceDirectory, paste0(currentSample$name,'assembled_deletion.index'))
  
  ## Initializing a list to store the deletion indices (which we will be removing when converting to fasta)
  assemblyDelIndexList <- list()
  
  ## Initializing a dataframe to store all the sequences that will be written to the fasta file
  allLocusSequenceDF <- data.frame(matrix('',sum(unlist(copyNumberList)),2),stringsAsFactors=F,check.names=F)
  
  ## Initializing a counter to know which row of the dataframe we are writing to
  currentRow <- 1
  for(currentLocus in currentPresentLoci){
    
    ## Pull out the copy number for the current locus
    copyNumber <- copyNumberList[[currentLocus]]
    
    ## Pull out the nucleotide assembly for the current locus
    currentLocusAssembly <- allLocusNucAssembly[[currentLocus]]
    
    ## If there is more than 1 copy for this locus, then we will write 2 sequences for this locus
    if(copyNumber > 1){
      
      ## Figure out which columns of row 2 are empty (this initially held het nucs, but now we want it all filled out)
      emptyCols <- colnames(currentLocusAssembly)[currentLocusAssembly[2,] == '']
      
      ## Fill in the empty row 2 columns with the nucs from row 1
      currentLocusAssembly[2,emptyCols] <- currentLocusAssembly[1,emptyCols]
      
      ## Figure out which columns in row 1 have deletions
      delIndices1 <- colnames(currentLocusAssembly[1,currentLocusAssembly[1,] == '.',drop=F])
      
      ## Figure out which columns in row 2 have deletions
      delIndices2 <- colnames(currentLocusAssembly[2,currentLocusAssembly[2,] == '.',drop=F])
      
      ## If there are deletions in row 1, then we want to store those positions in the delIndexList
      if(length(delIndices1) > 0){
        assemblyDelIndexList[[paste0(currentLocus,'_build1')]] <- c(assemblyDelIndexList[[currentLocus]], delIndices1)
        
        ## Write over the deletion positions with empty space
        currentLocusAssembly[1,currentLocusAssembly[1,] == '.'] <- ''
      }
      
      ## If there are deletions in row 2, then we want to store those positions in the delIndexList
      if(length(delIndices2) > 0){
        assemblyDelIndexList[[paste0(currentLocus,'_build2')]] <- c(assemblyDelIndexList[[currentLocus]], delIndices2)
        
        ## Write over the deletion positions with empty space
        currentLocusAssembly[2,currentLocusAssembly[2,] == '.'] <- ''
      }
      
      ## Merge the nucs from row 1 into a single string
      sequence1 <- paste0(currentLocusAssembly[1,], collapse='')
      
      ## Merge the nucs from row 2 into a single string
      sequence2 <- paste0(currentLocusAssembly[2,], collapse='')
      
      ## Give the nuc string a name according to what locus it is, and store it in the dataframe
      allLocusSequenceDF[currentRow,1] <- paste0('>',currentLocus,'_build1')
      allLocusSequenceDF[currentRow,2] <- sequence1
      
      ## Do the same for the second sequence
      allLocusSequenceDF[(currentRow+1),1] <- paste0('>',currentLocus,'_build2')
      allLocusSequenceDF[(currentRow+1),2] <- sequence2
      
      ## iterate the rows by 2 since we placed 2 sequence builds in the dataframe
      currentRow <- currentRow + 2
    }else{
      
      ## Figure out which columns in row 1 have deletions
      delIndices <- colnames(currentLocusAssembly[1,currentLocusAssembly[1,] == '.',drop=F])
      
      ## If there are deletions in row 1, then we want to store those positions in the delIndexList
      if(length(delIndices) > 0){
        assemblyDelIndexList[[paste0(currentLocus,'_build1')]] <- c(assemblyDelIndexList[[currentLocus]], delIndices)
        
        ## Overwrite the del positions with empty space
        currentLocusAssembly[1,delIndices] <- ''
      }
      
      ## Collapse the sequence into a single string
      sequence1 <- paste0(currentLocusAssembly[1,], collapse='')
      
      ## Give the sequence a name based on the current locus and write both to the dataframe
      allLocusSequenceDF[currentRow,1] <- paste0('>',currentLocus,'_build1')
      allLocusSequenceDF[currentRow,2] <- sequence1
      
      ## Iterate the row count by 1 since we placed 1 sequence build in the dataframe
      currentRow <- currentRow + 1
    }
  }
  
  ## Use the sequence string information to make a bed file
  allLocusBedDF <- data.frame(matrix('',nrow(allLocusSequenceDF),3),stringsAsFactors=F,check.names=F)
  allLocusBedDF[,1] <- allLocusSequenceDF[,1]
  allLocusBedDF[,2] <- 0
  allLocusBedDF[,3] <- nchar(allLocusSequenceDF[,2])
  
  allLocusBedDF[,1] <- tstrsplit(allLocusBedDF[,1],'>',fixed=T)[[2]]
  
  ## Write the deletion indices to a text file for easy read in later
  run.write_list(assemblyDelIndexList, currentSample$assembledDelIndex)
  
  ## Write the fasta file and bed file
  write.table(allLocusSequenceDF,file=currentSample$assembledRefPath,col.names=F,row.names=F,sep='\n',quote=F)
  write.table(allLocusBedDF,file=currentSample$assembledBedPath,col.names=F,row.names=F,sep=' ',quote=F)
  
  ## Return the sample object with the new file paths filled out
  return(currentSample)
}

## This function runs a bowtie2 alignment
run.bowtie2_assembly_alignment <- function(bowtie2_command, reference_index, threads, currentSample, assembledReferenceDirectory){
  ## Intitialize an output path for the SAM file
  currentSample$assembledSamPath <- file.path(assembledReferenceDirectory,paste0(currentSample$name,'.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-5 3', '-3 7', '--end-to-end', paste0('-p ',threads), '--score-min "L,-0.0,-0.2"',
                      '-I 75', '-X 1200','--very-sensitive',
                      paste0('-1 ',currentSample$fastq1path),
                      paste0('-2 ',currentSample$fastq2path),
                      '--no-unal','-a',
                      paste0('-S ',currentSample$assembledSamPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  currentSample$assembledSamPath <- normalizePath(currentSample$assembledSamPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',currentSample$name,'to',reference_index)
  
  return(output.sampleAlign)
}


## This function takes in a list and writes it to a text file
run.write_list <- function(listToWrite, filePath){
  firstLineBool = T
  
  for(listElemName in names(listToWrite)){
    listElem <- listToWrite[[listElemName]]
    
    write(listElemName, file=filePath, ncolumns=1, append=!firstLineBool, sep=' ')
    
    if(firstLineBool){
      firstLineBool = F
    }
    
    write(listElem, file=filePath, ncolumns=length(listElem), append=!firstLineBool, sep=' ')
  }
  
  return(filePath)
}

## This function takes in a text file reads it in as a list (its really only meant to work with run.write_list)
run.read_list <- function(filePath){
  #filePath <- file.path(assembledReferenceDirectory, paste0(currentSample$name,'assembled_deletion.index'))
  
  listToReturn <- list()
  
  conn <- file(filePath,open="r")
  linn <-readLines(conn)
  
  for (i in 1:length(linn)){
    if(i%%2 == 1){
      elemName <- linn[i]
    }else{
      elemValues <- linn[i]
      listToReturn[[elemName]] <- strsplit(elemValues, ' ', fixed=T)[[1]]
    }
  }
  close(conn)
  
  return(listToReturn)
}

### This function attempts to phase het SNPs using unpaired and paired-end reads
allele.snp_phaser <- function(samTable, currentLocus, currentLocusHetSnpDF, readAssignmentList, nucConv = T, hetRatio){
  cat('\nPhasing',ncol(currentLocusHetSnpDF),'het SNPs')
  ## Initialize a list for storing the phased snps (new elements are added when phasing cannot be completed)
  phasedList <- list()
  
  ## Initialze a count for keeping track of what list element we are on
  i <- 1
  
  ## Initialize the previous snp position
  prevPos <- 0
  
  ## Loop over all het snps
  for(snpPos in colnames(currentLocusHetSnpDF)){
    ## Pull out the nucleotides for the current position
    if(nucConv){
      snpVect <- names(nucListConv)[as.numeric(currentLocusHetSnpDF[,snpPos])]
    }else{
      snpVect <- as.vector(currentLocusHetSnpDF[,snpPos])
    }
    
    ## If this is the first iteration, then initialize a dataframe for storing the phased Snps
    if(prevPos == 0){
      ## Data frame with two rows for storing strand1 and strand2 SNPs (s1 and s2 from different list elements do not match up)
      phasedList[[i]] <- data.frame(matrix('',nrow=2,ncol=1),row.names=c('s1','s2'))
      
      ## Name the dataframe colums based on the current snp position
      colnames(phasedList[[i]]) <- snpPos
      
      ## Save the dataframe to the current list element
      phasedList[[i]][,snpPos] <- snpVect
    }else{
      
      ## Pull out all unpaired reads that overlap prevPos and snpPos for the current locus
      spanningUnpairedReads <- unique(samTable[read_name %in% readAssignmentList$singleLocusUnpairedReads][locus == currentLocus][startPos <= as.numeric(prevPos) & endPos >= as.numeric(snpPos)]$read_name)
      
      ## Pull out all paired reads that start before prevPos
      spanningPairedReads <- unique(samTable[read_name %in% readAssignmentList$singleLocusPairedReads][locus == currentLocus][startPos <= as.numeric(prevPos)]$read_name)# & endPos >= as.numeric(snpPos)]$read_name)
      
      ## Subset prevPos paired reads by the paired reads that also end after snpPos
      spanningPairedReads <- unique(samTable[read_name %in% spanningPairedReads][locus == currentLocus][endPos >= as.numeric(snpPos)]$read_name)
      
      ## Combine the unpaired and paired reads
      allSpanningReads <- c(spanningUnpairedReads, spanningPairedReads)
      
      ## If there are no reads that span both positions, then move on to the next list element and initalize a new phased dataframe
      ## else attempt to phase the prevPos snps with the snpPos snps using the overlapping reads
      if(length(allSpanningReads) == 0){
        i <- i + 1
        phasedList[[i]] <- data.frame(matrix('',nrow=2,ncol=1),row.names=c('s1','s2'))
        colnames(phasedList[[i]]) <- snpPos
        phasedList[[i]][,snpPos] <- snpVect
      }else{
        
        ## Pull out the prevPos snps
        s1Snp1 <- phasedList[[i]]['s1',prevPos]
        s2Snp1 <- phasedList[[i]]['s2',prevPos]
        
        ## Set the snpPos snps
        s1Snp2 <- snpVect[1]
        s2Snp2 <- snpVect[2]
        
        ## Initialize a new table based on the reads that overlap both positions from the SAM table
        miniTable <- samTable[read_name %in% allSpanningReads][locus == currentLocus]
        
        ## Subset the new table based on unique read sequences (get rid of duplicated entries)
        miniTable <- unique(miniTable, by=c('read_seq'))
        
        ## Apply the read_formatter function to the overlapping reads (creates a named vector of the read, adds in deletion positions)
        miniTable[,read_table:=mapply(allele.read_formatter, read_seq, startPos, endPos, reference_name)]
        
        ## Paired-end read combining and formatting
        miniTable[,read_table := allele.paired_read_combiner(read_name, startPos, endPos, read_table, spanningPairedReads), by=list(read_name)]
        
        ## Set the start position for paired end reads to be the minimum start position for both reads of the pair
        miniTable[,startPos := min(as.integer(startPos)), by=list(read_name)]
        
        ## Set the end position for paired end reads to be the maximum end position for both reads of the pair
        miniTable[,endPos := max(as.integer(endPos)), by=list(read_name)]
        
        ## Drop duplicate reads from the table (the information has been consolidated at this point)
        miniTable <- miniTable[startPos <= as.numeric(prevPos) & endPos >= as.numeric(snpPos)]
        miniTable <- unique(miniTable, by=c('startPos', 'endPos', 'locus'))
        
        ## Apply the snp_phaser_bool function to each named read vector to see if the s1-Snp1 and s1-Snp2 snps are in phase
        miniTable[,phasedS1:=mapply(allele.snp_phaser_bool, read_table, prevPos, s1Snp1, snpPos, s1Snp2)]
        miniTable[,phasedS2:=mapply(allele.snp_phaser_bool, read_table, prevPos, s2Snp1, snpPos, s2Snp2)]
        
        ## Save the results to a bool vector (TRUE/FALSE if the snps are in phase for each read)
        case1 <- miniTable$phasedS1 | miniTable$phasedS2
        
        ## Apply the snp_phaser_bool function to each named read vector to see if the s2-Snp1 and s1-Snp2 snps are in phase
        miniTable[,phasedS1:=mapply(allele.snp_phaser_bool, read_table, prevPos, s2Snp1, snpPos, s1Snp2)]
        miniTable[,phasedS2:=mapply(allele.snp_phaser_bool, read_table, prevPos, s1Snp1, snpPos, s2Snp2)]
        
        ## Save the results to a bool vector (TRUE/FALSE if the snps are in phase for each read)
        case2 <- miniTable$phasedS1 | miniTable$phasedS2
        
        ## If case1 has more TRUE values than case2, then s1-Snp1 and s1-Snp2 are in phase and saved to the phased list
        ## Else if case2 hase more TRUE values than case1, then s2-Snp1 and s1-Snp2 are in phased and saved to the phased list
        ## else no phased could be determined and we move on to the next list element and initialze a new phased dataframe
        if(((sum(case1)/nrow(miniTable)) > hetRatio) & (sum(case1) > sum(case2))){
          phasedList[[i]]['s1',snpPos] <- s1Snp2
          phasedList[[i]]['s2',snpPos] <- s2Snp2
        }else if(((sum(case2)/nrow(miniTable)) > hetRatio) & (sum(case2) > sum(case1))){
          phasedList[[i]]['s1',snpPos] <- s2Snp2
          phasedList[[i]]['s2',snpPos] <- s1Snp2
        }else{
          i <- i + 1
          phasedList[[i]] <- data.frame(matrix('',nrow=2,ncol=1),row.names=c('s1','s2'))
          colnames(phasedList[[i]]) <- snpPos
          phasedList[[i]][,snpPos] <- snpVect
        }
      }
    }
    
    ## Set the prevPos to be the current snpPos, then move on to the next snpPos
    prevPos <- snpPos
  }
  return(phasedList)
}

## This function combines paired-end read named vectors into a single named vector
allele.paired_read_combiner <- function(readName, startPos, endPos, pairedList, spanningPairedReads){
  
  ## If there are two startPos values for this entry and the read name is in the spanningPairedReads list, then continue with combining
  ## else return the current named vector with no editing
  if(length(startPos) == 2 & readName %in% spanningPairedReads){
    
    ## Unlist the named paired-end read vectors
    combinedList <- unlist(pairedList)
    
    ## Pull out the position names
    realPositionList <- names(combinedList)
    
    ## Format the starting position for the paired-end read
    minPosition <- min(as.numeric(startPos))
    
    ## Format the ending position for the paired-end read
    maxPosition <- max(as.numeric(endPos))
    
    ## Set a full spanning position list that goes from the min to the max position
    fullPositionList <- as.character(minPosition:maxPosition)
    
    ## Compare the fullPosition list to the realPositionList to see what positions are not actually found in the reads (non-overlapping paired-end reads)
    absentPositionList <- setdiff(fullPositionList, realPositionList)
    
    ## If there are absent positions, then we want to set those positions to 'N' in the named vector we will return
    if(length(absentPositionList) > 0){
      
      ## Initialize the named vector
      absentPositionVect <- vector(mode="character", length=length(absentPositionList))
      
      ## Set absent positions to 'N'
      absentPositionVect[1:length(absentPositionVect)] <- 'N'
      
      ## Set the names of the absent position vector 
      names(absentPositionVect) <- as.character(absentPositionList)
      
      ## Add the named absent position vector to the named paired-end read vector
      combinedList <- c(combinedList, absentPositionVect)
    }
    
    ## Return the named combined list vector
    return(list(list(combinedList[fullPositionList])))
  }else{
    return(list(pairedList))
  }
}

## Format reads for snp phasing (fill in deletion positions with '.')
allele.read_formatter <- function(readSeq, startPos, endPos, refAlleleName){
  
  ## Format the start position
  startPos <- as.numeric(startPos)
  
  ## Format the end position
  endPos <- as.numeric(endPos)
  
  ## Testing if there deletions found in the reference for the current read
  if((endPos-startPos+1) > nchar(readSeq)){
    
    delIndexList <- which(startPos:endPos %in% currentDelIndex[[refAlleleName]])
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    subStringList <- c()
    for(delIndex in delIndexList){
      preString <- substr(readSeq, 1, delIndex-1)
      postString <- substr(readSeq, delIndex, nchar(readSeq))
      
      readSeq <- paste0(preString, '.', postString)
    }
  }
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- startPos:endPos
  
  ## Turn the sequence string into a list
  seqList <- strsplit(readSeq,'')[[1]]
  
  if(length(seqList) != length(fullIndex)){
    cat('\n',startPos,endPos,refAlleleName)
  }
  
  names(seqList) <- fullIndex
  
  seqListList <- list(seqList)
  return(seqListList)
}

## Return if the two snps are present in the same read
allele.snp_phaser_bool <- function(readTab, snp1Pos, snp1Nuc, snp2Pos, snp2Nuc, read_name){
  
  ## If snp1 is found in the read and snp2 is found in the read, then the snps are phased and return TRUE
  return(readTab[[snp1Pos]] == snp1Nuc & readTab[[snp2Pos]] == snp2Nuc)
}

## Fill in positions from fillNucList into allLocusNucAssembly
run.fill_nucs <- function(allLocusNucAssembly, currentLocus, nucListConv, fillNucList){
  for(fillPos in names(fillNucList)){
    fillNucs <- fillNucList[[fillPos]]
    if(length(fillNucs) > 2){
      fillNucs <- fillNucs[1:2]
      cat('\nFound 3 nucs for this locus-------------------------')
    }
    allLocusNucAssembly[[currentLocus]][1:length(fillNucs),fillPos] <- unlist(nucListConv[fillNucs], use.names=F)
  }
  return(allLocusNucAssembly)
}

allele.format_vcf <- function(currentVcf, currentDelIndex, depthThreshold, validatedDelList){
  ## Initialize global positions
  currentVcf$pos <- currentVcf$V2
  
  ## Pull out the locus names
  currentVcf$locus <- tstrsplit(currentVcf$V1, '_', fixed=T)[[1]]
  
  ## Set up the depth statistic
  currentVcf$depth <- tstrsplit(currentVcf$V8, ';', fixed=T)[[1]]
  
  ## Discard Indels
  currentVcf <- currentVcf[depth != 'INDEL']
  
  ## Set the depth field
  currentVcf$depth <- as.integer(tstrsplit(currentVcf$depth, '=', fixed=T)[[2]])
  
  ## Drop snps that do not pass the depth threshold
  currentVcf <- currentVcf[depth >= depthThreshold]
  
  ## Update the global positions based on the deletion index
  currentVcf[,pos := mapply(allele.add_del_index, pos, V1)]
  
  ## Add the deletion positions to the vcf frame so they can be used for allele calling
  for(currentBuild in names(validatedDelList)){
    currentLocus <- strsplit(currentBuild,'_',fixed=T)[[1]][1]
    delIndex <- currentDelIndex[[currentBuild]]
    
    for(currentPos in delIndex){
      currentPos <- as.integer(currentPos)
      addList <- data.table('V1'=currentBuild,'V2'=0,'V3'='.','V4'='.','V5'='.',
                      'V6'=0,'V7'='.','V8'='.','V9'='.','V10'='0/0',
                      'pos'=currentPos,'locus'=currentLocus,'depth'=depthThreshold)
      currentVcf <- rbindlist(list(currentVcf, addList))
      
    }
  }
  
  ## Initialize snp1 calls
  currentVcf$snp1 <- ''
  
  ## Initialize snp2 calls
  currentVcf$snp2 <- ''
  
  ## Format the geno field call
  currentVcf$V10 <- tstrsplit(currentVcf$V10, ':', fixed=T)[[1]]
  
  ## Discard no calls
  currentVcf <- currentVcf[V10 != './.']
  
  ## Set up the geno calls for snp1 and snp2
  snp1Call <- as.integer(tstrsplit(currentVcf$V10, '/', fixed=T)[[1]])+1
  snp2Call <- as.integer(tstrsplit(currentVcf$V10, '/', fixed=T)[[2]])+1
  
  ## Assign nucs to the snp columns based on the geno calls
  callMatrix <- cbind(currentVcf$V4, currentVcf$V5)
  currentVcf$snp1 <- sapply(1:length(snp1Call), function(x){callMatrix[x,snp1Call[x]]})
  currentVcf$snp2 <- sapply(1:length(snp2Call), function(x){callMatrix[x,snp2Call[x]]})
  
  ## Set the regions based on the kirExonCoords list
  currentVcf[,region := mapply(allele.set_vcf_region, locus, pos)]
  
  return(currentVcf)
}

allele.set_vcf_region <- function(locus, pos){
  return(names(kirExonCoords[[locus]])[sapply(kirExonCoords[[locus]],function(x) pos %in% x)])
}

allele.add_del_index <- function(pos, currentBuild){
  sumAmount <- sum(pos >= currentDelIndex[[currentBuild]])
  
  newSumAmount <- sum(sumAmount+pos >= currentDelIndex[[currentBuild]])
  
  while(sumAmount != newSumAmount){
    sumAmount <- newSumAmount
    newSumAmount <- sum(sumAmount+pos >= currentDelIndex[[currentBuild]])
  }
  
  return(as.integer(sumAmount + pos))
}

allele.find_hom_region_alleles <- function(homRegionVect, kirExonCoords, currentLocusSnpFrame, kirAlleleDFList, onlyExons=F){
  cat('\nFinding homozygous region allele matches...')
  
  if(onlyExons){
    homRegionVect <- homRegionVect[sapply(homRegionVect, grepl, pattern='E', fixed=T)]
  }
  
  homRegionAlleleList <- list()
  cat('\n\tExamining:')
  for(currentRegion in homRegionVect){
    cat('', currentRegion)
    currentRegionPosVect <- kirExonCoords[[currentLocus]][[currentRegion]]
    currentRegionPosVect <- intersect(as.character(currentRegionPosVect), colnames(currentLocusSnpFrame))
    
    if(length(currentRegionPosVect) == 0){
      next
    }
    
    currentRegionAlleleFrame <- kirAlleleDFList[[currentLocus]][,currentRegionPosVect,drop=F]
    possibleAlleleList <- rownames(currentRegionAlleleFrame)
    
    for(currentPos in currentRegionPosVect){
      possibleAlleleList <- possibleAlleleList[sapply(currentRegionAlleleFrame[possibleAlleleList,currentPos,drop=F], function(x){x %in% c(currentLocusSnpFrame[,currentPos], 'N')})]
      if(length(possibleAlleleList) == 0){
        break
      }
    }
    
    if(length(possibleAlleleList) > 0){
      homRegionAlleleList[[currentRegion]] <- possibleAlleleList
    }
  }
  
  cat('\nCondensing all homozygous region allele matches.')
  
  if(onlyExons){
    onlyExonNames <- names(homRegionAlleleList)[sapply(names(homRegionAlleleList), grepl, pattern='E', fixed=T)]
    
    if(length(onlyExonNames) == 0){
      cat('\nUh oh! The onlyExons paramater was set to TRUE, but no homozygous exons were found.')
      stop()
    }
    homRegionAlleleVect <- Reduce(intersect, homRegionAlleleList[onlyExonNames])
  }else{
    homRegionAlleleVect <- Reduce(intersect, homRegionAlleleList)
  }
  
  cat('\n\tMatched alleles: ', homRegionAlleleVect)
  
  return(homRegionAlleleVect)
}

allele.find_het_region_alleles <- function(hetRegionVect, kirExonCoords, s1NucVect, s2NucVect, currentLocusSnpFrame, homRegionAlleleVect, kirAlleleDFList, onlyExons=F){
  cat('\nFinding heterozygous region allele matches...')
  
  if(onlyExons){
    hetRegionVect <- hetRegionVect[sapply(hetRegionVect, grepl, pattern='E', fixed=T)]
  }
  
  #names(s1RegionAlleleList)[sapply(names(s1RegionAlleleList), grepl, pattern='E', fixed=T)]
  s1RegionAlleleList <- list()
  s2RegionAlleleList <- list()
  
  cat('\n\tExamining:')
  for(currentRegion in hetRegionVect){
    cat('',currentRegion)
    currentRegionPosVect <- kirExonCoords[[currentLocus]][[currentRegion]]
    
    currentRegionHetPosVect <- intersect(as.character(currentRegionPosVect), names(s1NucVect))
    currentRegionHomPosVect <- intersect(as.character(currentRegionPosVect), colnames(currentLocusSnpFrame))
    currentRegionHomPosVect <- setdiff(currentRegionHomPosVect, currentRegionHetPosVect)
    
    currentRegionPosVect <- union(currentRegionHetPosVect, currentRegionHomPosVect)
    
    if(length(currentRegionPosVect) == 0){
      next
    }
    
    currentRegionAlleleFrame <- kirAlleleDFList[[currentLocus]][,currentRegionPosVect,drop=F]
    
    possibleAlleleList <- intersect(rownames(currentRegionAlleleFrame), homRegionAlleleVect)
    #possibleAlleleList <- rownames(currentRegionAlleleFrame)
    
    noMatchingAlleles = F
    for(currentPos in currentRegionHomPosVect){
      possibleAlleleList <- possibleAlleleList[sapply(currentRegionAlleleFrame[possibleAlleleList,currentPos,drop=F], function(x){x %in% c(currentLocusSnpFrame[1,currentPos], 'N')})]
      if(length(possibleAlleleList) == 0){
        noMatchingAlleles = T
        break
      }
    }
    
    if(noMatchingAlleles){
      cat(' No perfect allele matches for this region.')
      next
    }
    
    s1PossibleAlleleList <- possibleAlleleList
    s2PossibleAlleleList <- possibleAlleleList
    
    for(currentPos in currentRegionHetPosVect){
      s1PossibleAlleleList <- s1PossibleAlleleList[sapply(currentRegionAlleleFrame[s1PossibleAlleleList,currentPos,drop=F], function(x){x %in% c(s1NucVect[currentPos], 'N')})]
      if(length(s1PossibleAlleleList) == 0){
        break
      }
    }
    
    for(currentPos in currentRegionHetPosVect){
      s2PossibleAlleleList <- s2PossibleAlleleList[sapply(currentRegionAlleleFrame[s2PossibleAlleleList,currentPos,drop=F], function(x){x %in% c(s2NucVect[currentPos], 'N')})]
      if(length(s2PossibleAlleleList) == 0){
        break
      }
    }
    
    if(length(s1PossibleAlleleList) > 0){
      s1RegionAlleleList[[currentRegion]] <- s1PossibleAlleleList
    }
    if(length(s2PossibleAlleleList) > 0){
      s2RegionAlleleList[[currentRegion]] <- s2PossibleAlleleList
    }
  }
  
  cat('\nCondensing all heterozygous region allele matches.')
  
  if(onlyExons){
    s1OnlyExonNames <- names(s1RegionAlleleList)[sapply(names(s1RegionAlleleList), grepl, pattern='E', fixed=T)]
    s2OnlyExonNames <- names(s2RegionAlleleList)[sapply(names(s2RegionAlleleList), grepl, pattern='E', fixed=T)]
    
    if(length(s1OnlyExonNames) == 0 | length(s2OnlyExonNames) == 0){
      cat('\nUh oh! The onlyExons paramater was set to TRUE, but no heterozygous exons were found.')
      stop()
    }
    s1RegionAlleleVect <- Reduce(intersect, s1RegionAlleleList[s1OnlyExonNames])
    s2RegionAlleleVect <- Reduce(intersect, s2RegionAlleleList[s2OnlyExonNames])
  }else{
    s1RegionAlleleVect <- Reduce(intersect, s1RegionAlleleList)
    s2RegionAlleleVect <- Reduce(intersect, s2RegionAlleleList)
  }
  
  cat('\n\tStrand1 matched alleles: ', s1RegionAlleleVect)
  cat('\n\tStrand2 matched alleles: ', s2RegionAlleleVect)
  
  return(list('s1'=s1RegionAlleleVect, 's2'=s2RegionAlleleVect))
}

allele.KIR2DL4_fix <- function(samTable, alleleBuildList, readAssignmentList, nucListConv){
  
  ##### Fix 1
  # pos 9929 has lots of reads that do not match up together
  # Filter the pos 9929 reads based on pos 9930 'G' and 9931 'T', which are invariant in 2DL4
  
  cat('\n\tPos 9929')
  
  ## First we pull out the reads that align to KIR2DL4 from position 9929-9931
  samSubsetTable <- samTable[locus == 'KIR2DL4'][read_name %in% readAssignmentList$singleLocusPairedReads][startPos <= 9929][endPos >= 9931]
  samSubsetTable[,read_table := mapply(allele.read_formatter, read_seq, startPos, endPos, reference_name)]
  
  ## If there is more than 0 matching reads, then continue with the fix
  if(nrow(samSubsetTable) > 0){
    
    ## Pull out the pos 9929 nuc from reads that have 'G' at pos 9930 and 'T' at pos 9931
    posNucList <- sapply(samSubsetTable$read_table, function(x){
      if(x['9930'] == 'G' & x['9931'] == 'T'){
        return(x['9929'])
      }else{
        return(NULL)
      }
    })
    
    ## Convert the list into a vector without NULL values
    posNucVect <- unlist(posNucList)
    
    ## If there are more than 0 nucs pulled from pos 9929, then continue
    if(length(posNucVect) > 0){
      
      ## Make a table to easily count up the number of nuc calls
      posNucTable <- table(posNucVect)
      
      ## Turn the table names into a row index vector
      rowNameVect <- unlist(nucListConv[names(posNucTable)], use.names=F)
      
      ## Use the row index vector to insert the nuc counts into the KIR2DL4 allele build
      alleleBuildList[['KIR2DL4']][rowNameVect,'9929'] <- as.vector(posNucTable)
    }else{
      
      ## IF ther are 0 nucs pulled form pos 9929, then write 0's to the nuc counts of the KIR2DL4 allele build
      alleleBuildList[['KIR2DL4']][,'9929'] <- c(0,0,0,0,0)
    }
  }
  return(alleleBuildList)
}

allele.KIR3DS1_fix <- function(allLocusNucAssembly){
  
  ##### Fix 1
  # pos 484 has >100 bp of surrounding homology, making this a really hard het call to filter out
  # The fix is to remove the SNP calls from this position if it is a het call
  # There should still be enough allele differentiating snps to call specific alleles
  
  cat('\n\tPos 484')
  if(allLocusNucAssembly[['KIR3DS1']]['2','484'] != ''){
    allLocusNucAssembly[['KIR3DS1']][,'484'] <- c('','')
  } 
  
  return(allLocusNucAssembly)
}

sam_flag_convert <- function(samFlag){
  flagList <- list(readPaired=F,properPairMapping=F,readUnmapped=F,mateUnmapped=F,
                   readReverseStrand=F,mateReverseStrand=F,firstInPair=F,secondInPair=F,
                   notPrimaryAlignment=F,readFailsQualityChecks=F,readIsPcrOrOpticalDuplicate=F,
                   supplementaryAlignment=F)
  
  #read paired (0x1)
  #read mapped in proper pair (0x2)
  #read unmapped (0x4)
  #mate unmapped (0x8)
  #read reverse strand (0x10)
  #mate reverse strand (0x20)
  #first in pair (0x40)
  #second in pair (0x80)
  #not primary alignment (0x100)
  #read fails platform/vendor quality checks (0x200)
  #read is PCR or optical duplicate (0x400)
  #supplementary alignment (0x800)
  
  
  hexInt <- as.integer(samFlag)
  
  if(hexInt >= 2048){
    hexInt <- hexInt - 2048
    flagList$supplementaryAlignment <- T
  }
  
  if(hexInt >= 1024){
    hexInt <- hexInt - 1024
    flagList$readIsPcrOrOpticalDuplicate <- T
  }
  
  if(hexInt >= 512){
    hexInt <- hexInt - 512
    flagList$readFailsQualityChecks <- T
  }
  
  if(hexInt >= 256){
    hexInt <- hexInt - 256
    flagList$notPrimaryAlignment <- T
  }
  
  if(hexInt >= 128){
    hexInt <- hexInt - 128
    flagList$secondInPair <- T
  }
  
  if(hexInt >= 64){
    hexInt <- hexInt - 64
    flagList$firstInPair <- T
  }
  
  if(hexInt >= 32){
    hexInt <- hexInt - 32
    flagList$mateReverseStrand <- T
  }
  
  if(hexInt >= 16){
    hexInt <- hexInt - 16
    flagList$readReverseStrand <- T
  }
  
  if(hexInt >= 8){
    hexInt <- hexInt - 8
    flagList$mateUnmapped <- T
  }
  
  if(hexInt >= 4){
    hexInt <- hexInt - 4
    flagList$readUnmapped <- T
  }
  
  if(hexInt >= 2){
    hexInt <- hexInt - 2
    flagList$properPairMapping <- T
  }
  
  if(hexInt >= 1){
    hexInt <- hexInt - 1
    flagList$readPaired <- T
  }
  
  return(flagList)
}

## Condenses allele fram into only positions that have multiple nucleotides
make_unique_pos_frame <- function(pos_frame){
  num_unique <- lapply(pos_frame, num_unique_nuc)
  unique_pos_frame <- pos_frame[,names(num_unique[num_unique > 1]), drop=F]
  return(unique_pos_frame)
}

## Returns boolean
is_nuc <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G', '.')
  return(as.character(chr) %in% nuc_list)
}

## Counts how many unique nucleotides are in a vector
num_unique_nuc <- function(vector_to_check){
  return(sum(is_nuc(names(table(vector_to_check)))))
}

nucMatching <- function(x,assemblyRownames){
  return(assemblyRownames[x])
}

nucListConv <- list('A'=1,
                    'T'=2,
                    'C'=3,
                    'G'=4,
                    '.'=5)

hetNucConv <- list('R'=c('A','G'),
                   'Y'=c('C','T'),
                   'K'=c('G','T'),
                   'M'=c('A','C'),
                   'S'=c('G','C'),
                   'W'=c('A','T'),
                   'E'=c('A','.'),
                   'F'=c('T','.'),
                   'I'=c('C','.'),
                   'J'=c('G','.'))

strandConvList <- list('s1'='s2',
                       's2'='s1')

kirExonCoords <- list('KIR3DP1'=list('5UTR'=1:267,'E1'=268:301,'I1'=302:1323,'E2'=1324:1359,
                                     'I2'=1360:2117,'E3'=2118:2402,'I3'=2403:3852,'E4'=3853:4152,
                                     'I4'=4153:5679,'E5'=5680:5973,'3UTR'=5974:5982),
                      'KIR2DS5'=list('5UTR'=1:268,'E1'=269:304,'I1'=305:1799,'E2'=1800:1835,
                                     'I2'=1836:2562,'E3'=2563:2844,'I3'=2845:4279,'E4'=4280:4579,
                                     'I4'=4580:6107,'E5'=6108:6401,'I5'=6402:9554,'E6'=9555:9605,
                                     'I6'=9606:13870,'E7'=13871:13975,'I7'=13976:14437,'E8'=14438:14490,
                                     'I8'=14491:14588,'E9'=14589:14630,'3UTR'=14631:15275),
                      'KIR2DL3'=list('5UTR'=1:268,'E1'=269:302,'I1'=303:1229,'E2'=1230:1265,
                                     'I2'=1266:1992,'E3'=1993:2274,'I3'=2274:3702,'E4'=3703:4002,
                                     'I4'=4003:5520,'E5'=5521:5814,'I5'=5815:9095,'E6'=9096:9146,
                                     'I6'=9147:13411,'E7'=13412:13516,'I7'=13517:13978,'E8'=13979:14031,
                                     'I8'=14032:14129,'E9'=14130:14282,'3UTR'=14283:14816),
                      'KIR2DP1'=list('5UTR'=1:267,'E1'=268:301,'I1'=302:1358,'E2'=1359:1394,
                                     'I2'=1395:2119,'E3'=2120:2401,'I3'=2402:3854,'E4'=3855:4153,
                                     'I4'=4154:5679,'E5'=5680:5973,'I5'=5974:9141,'E6'=9142:9192,
                                     'I6'=9193:11728,'E7'=11729:11833,'I7'=11834:12295,'E8'=12296:12348,
                                     'I8'=12349:12446,'E9'=12447:12623,'3UTR'=12624:13133),
                      'KIR2DS3'=list('5UTR'=1:300,'E1'=301:334,'I1'=335:1639,'E2'=1640:1675,
                                     'I2'=1676:2402,'E3'=2403:2684,'I3'=2685:4118,'E4'=4119:4418,
                                     'I4'=4419:5937,'E5'=5938:6231,'I5'=6232:9384,'E6'=9385:9435,
                                     'I6'=9436:13700,'E7'=13701:13805,'I7'=13806:14267,'E8'=14268:14320,
                                     'I8'=14321:14418,'E9'=14419:14460,'3UTR'=14461:15105),
                      'KIR2DS2'=list('5UTR'=1:300,'E1'=301:334,'I1'=335:1126,'E2'=1127:1162,
                                     'I2'=1163:1889,'E3'=1890:2170,'I3'=2171:3598,'E4'=3599:3898,
                                     'I4'=3899:5416,'E5'=5417:5710,'I5'=5711:8879,'E6'=8880:8930,
                                     'I6'=8931:13196,'E7'=13197:13301,'I7'=13302:13763,'E8'=13764:13816,
                                     'I8'=13817:13914,'E9'=13915:13956,'3UTR'=13957:14580),
                      'KIR2DL4'=list('5UTR'=1:267,'E1'=268:307,'I1'=308:506,'E2'=507:542,
                                     'I2'=543:1476,'E3'=1477:1761,'I3'=1762:2497,'I4'=2498:2640,
                                     'E5'=2641:2934,'I5'=2935:5529,'E6'=5530:5580,'I6'=5581:9824,
                                     'E7'=9825:9929,'I7'=9930:10390,'E8'=10391:10443,'I8'=10444:10542,
                                     'E9'=10543:10812,'3UTR'=10813:11219),
                      'KIR3DL3'=list('5UTR'=1:321,'E1'=322:357,'I1'=358:1071,'E2'=1072:1107,
                                     'I2'=1108:1877,'E3'=1878:2162,'I3'=2163:3435,'E4'=3436:3735,
                                     'I4'=3736:5317,'E5'=5318:5611,'I5'=5612:7212,'I6'=7213:11086,
                                     'E7'=11087:11191,'I7'=11192:11653,'E8'=11654:11706,'I8'=11707:11804,
                                     'E9'=11805:11930,'3UTR'=11931:12491),
                      'KIR3DL1'=list('5UTR'=1:267,'E1'=268:301,'I1'=302:1323,'E2'=1324:1359,
                                     'I2'=1360:2104,'E3'=2105:2389,'I3'=2390:3502,'E4'=3503:3802,
                                     'I4'=3803:5354,'E5'=5355:5648,'I5'=5649:8820,'E6'=8821:8871,
                                     'I6'=8871:15550,'E7'=15551:15655,'I7'=15656:16118,'E8'=16119:16171,
                                     'I8'=16172:16290,'E9'=16291:16467,'3UTR'=16468:16977),
                      'KIR3DS1'=list('5UTR'=1:479,'E1'=480:513,'I1'=514:1513,'E2'=1514:1549,
                                     'I2'=1550:2293,'E3'=2294:2578,'I3'=2579:4066,'E4'=4067:4366,
                                     'I4'=4367:5946,'E5'=5947:6240,'I5'=6241:9410,'E6'=9411:9461,
                                     'I6'=9462:13741,'E7'=13742:13847,'I7'=13848:14309,'E8'=14310:14360,
                                     'I8'=14361:14458,'E9'=14459:14466,'3UTR'=14467:15145),
                      'KIR2DL2'=list('5UTR'=1:300,'E1'=301:334,'I1'=335:1242,'E2'=1243:1278,
                                     'I2'=1279:2015,'E3'=2016:2296,'I3'=2297:3724,'E4'=3725:4024,
                                     'I4'=4025:5538,'E5'=5539:5832,'I5'=5833:9111,'E6'=9112:9162,
                                     'I6'=9163:13427,'E7'=13428:13529,'I7'=13530:13991,'E8'=13992:14044,
                                     'I8'=14045:14142,'E9'=14143:14319,'3UTR'=14320:14829),
                      'KIR3DL2'=list('5UTR'=1:268,'E1'=269:302,'I1'=303:1031,'E2'=1032:1067,
                                     'I2'=1068:1809,'E3'=1810:2094,'I3'=2095:3558,'E4'=3559:3858,
                                     'I4'=3859:5437,'E5'=5438:5731,'I5'=5732:8898,'E6'=8899:8949,
                                     'I6'=8950:15632,'E7'=15633:15737,'I7'=15738:16197,'E8'=16198:16250,
                                     'I8'=16251:16349,'E9'=16350:16559,'3UTR'=16560:17044),
                      'KIR2DS4'=list('5UTR'=1:267,'E1'=268:301,'I1'=302:2683,'E2'=2684:2719,
                                     'I2'=2720:3446,'E3'=3447:3728,'I3'=3729:5180,'E4'=5181:5480,
                                     'I4'=5481:7032,'E5'=7033:7326,'I5'=7327:10495,'E6'=10496:10546,
                                     'I6'=10547:14813,'E7'=14814:14918,'I7'=14919:15380,'E8'=15381:15433,
                                     'I8'=15434:15531,'E9'=15532:15573,'3UTR'=15574:16197),
                      'KIR2DL1'=list('5UTR'=1:268,'E1'=269:302,'I1'=303:1266,'E2'=1267:1302,
                                     'I2'=1303:2030,'E3'=2031:2312,'I3'=2313:3753,'E4'=3754:4053,
                                     'I4'=4054:5587,'E5'=5588:5881,'I5'=5882:9037,'E6'=9038:9088,
                                     'I6'=9089:13357,'E7'=13358:13459,'I7'=13460:13921,'E8'=13922:13974,
                                     'I8'=13975:14072,'E9'=14073:14249,'3UTR'=14250:14759),
                      'KIR2DS1'=list('5UTR'=1:267,'E1'=268:301,'I1'=302:1265,'E2'=1266:1301,
                                     'I2'=1302:2028,'E3'=2029:2310,'I3'=2311:3749,'E4'=3750:4049,
                                     'I4'=4050:5574,'E5'=5575:5868,'I5'=5869:9022,'E6'=9023:9073,
                                     'I6'=9074:13336,'E7'=13337:13441,'I7'=13442:13903,'E8'=13904:13956,
                                     'I8'=13957:14054,'E9'=14055:14096,'3UTR'=14097:14720),
                      'KIR2DL5'=list('5UTR'=1:547,'E1'=548:581,'I1'=582:1451,'E2'=1452:1487,
                                     'I2'=1488:2249,'E3'=2250:2534,'I3'=2535:3277,'I4'=3278:3409,
                                     'E5'=3410:3703,'I5'=3704:5875,'E6'=5876:5926,'I6'=5927:8692,
                                     'E7'=8693:8797,'I7'=8797:9259,'E8'=9260:9312,'I8'=9313:9412,
                                     'E9'=9413:9682,'3UTR'=9683:10097))


### Setting paths to random forest models
rfSaveDirectory <- 'Resources/gc_resources/rf_models/'
rfKIR2DL1path <- file.path(rfSaveDirectory,'rfKIR2DL1.Rdata')
rfKIR2DL2path <- file.path(rfSaveDirectory,'rfKIR2DL2.Rdata')
rfKIR2DL3path <- file.path(rfSaveDirectory,'rfKIR2DL3.Rdata')
rfKIR2DL4path <- file.path(rfSaveDirectory,'rfKIR2DL4.Rdata')
rfKIR2DL5path <- file.path(rfSaveDirectory,'rfKIR2DL5.Rdata')
rfKIR2DS1path <- file.path(rfSaveDirectory,'rfKIR2DS1.Rdata')
rfKIR2DS2path <- file.path(rfSaveDirectory,'rfKIR2DS2.Rdata')
rfKIR2DS3path <- file.path(rfSaveDirectory,'rfKIR2DS3.Rdata')
rfKIR2DS4path <- file.path(rfSaveDirectory,'rfKIR2DS4.Rdata')
rfKIR2DS5path <- file.path(rfSaveDirectory,'rfKIR2DS5.Rdata')
rfKIR2DP1path <- file.path(rfSaveDirectory,'rfKIR2DP1.Rdata')
rfKIR3DL1path <- file.path(rfSaveDirectory,'rfKIR3DL1.Rdata')
rfKIR3DS1path <- file.path(rfSaveDirectory,'rfKIR3DS1.Rdata')
rfKIR3DP1path <- file.path(rfSaveDirectory,'rfKIR3DP1.Rdata')
rfKIR3DL2path <- file.path(rfSaveDirectory,'rfKIR3DL2.Rdata')
### /Paths

### Initialize a list for storing the rf model paths
rfAllPathList <- list()
rfAllPathList[['KIR2DL1']] <- rfKIR2DL1path
rfAllPathList[['KIR2DL2']] <- rfKIR2DL2path
rfAllPathList[['KIR2DL3']] <- rfKIR2DL3path
rfAllPathList[['KIR2DL4']] <- rfKIR2DL4path
rfAllPathList[['KIR2DL5']] <- rfKIR2DL5path
rfAllPathList[['KIR2DS1']] <- rfKIR2DS1path
rfAllPathList[['KIR2DS2']] <- rfKIR2DS2path
rfAllPathList[['KIR2DS3']] <- rfKIR2DS3path
rfAllPathList[['KIR2DS4']] <- rfKIR2DS4path
rfAllPathList[['KIR2DS5']] <- rfKIR2DS5path
rfAllPathList[['KIR2DP1']] <- rfKIR2DP1path
rfAllPathList[['KIR3DL1']] <- rfKIR3DL1path
rfAllPathList[['KIR3DS1']] <- rfKIR3DS1path
rfAllPathList[['KIR3DP1']] <- rfKIR3DP1path
rfAllPathList[['KIR3DL2']] <- rfKIR3DL2path
### /Initialize
































  