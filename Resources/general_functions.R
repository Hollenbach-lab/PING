
## This function checks to make sure the output of system2 is valid
check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}

## Check to make sure samtools is accessible
samtools <- system2('which', c('samtools'), stdout=T, stderr=T)
check.system2_output(samtools, 'samtools not found')

## Check to make sure samtools is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')

## This function finds samples in sampleDirectory and turns them into objects for downstream use
general.paired_sample_objects <- function(rawFastqDirectory,fastqPattern,resultsDirectory,shortNameDelim){
  #####
  ## This function takes in a directory and a file name pattern and attempts to pair fastq files
  ## Returns a list of sample objects that contain the paired fastq file names
  #####
  
  cat("\nAttempting automatic fastq pairing in", rawFastqDirectory, "using", fastqPattern)
  
  ## Find all the files in sampleDirectory that match fastqPattern
  unpairedFastqList <- list.files(path=rawFastqDirectory, pattern=fastqPattern)
  
  ## To pair reads, we will split the file names by fastqPattern, then continuously chop a
  ## character off the end of each name until the number of unique names is exactly half of the total names
  
  ## Setting up an initial fastq list that splits the files names by fastqPattern
  strList <- sapply(unpairedFastqList, function(x) str_split(x, fastqPattern)[[1]][1])
  
  ## Shorten names based on delim
  if(shortNameDelim != ''){
    strList <- sapply(strList, function(x) str_split(x, fixed(shortNameDelim))[[1]][1])
  }
  
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
  
  cat("\nFound", length(uniqueFastqList), "samples in", rawFastqDirectory)
  
  ## Creating the sample object class
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
  
  ## Initializing a sample object list. This will be returned
  output.sampleList <- list()
  
  for(i in 1:length(pairedFastqList)){
    
    ## Pulling the current working element out of the list
    pairedFastq <- pairedFastqList[i]
    
    ## Creating an absolute path to the first fastq file
    rawfastq1path <- normalizePath(file.path(rawFastqDirectory, pairedFastq[[1]][1]), mustWork=T)
    
    ## Creating a absolute path to the second fastq file
    rawfastq2path <- normalizePath(file.path(rawFastqDirectory, pairedFastq[[1]][2]), mustWork=T)
    
    ## Checking if the first fastq file is gzipped
    gzip <- substr(rawfastq1path, nchar(rawfastq1path)-2, nchar(rawfastq1path)) == '.gz'
    
    ## Building a sample object and adding it to sampleList
    output.sampleList[[names(pairedFastq)]] <- sample(name=names(pairedFastq),
                                                      rawfastq1path=rawfastq1path,
                                                      rawfastq2path=rawfastq2path,
                                                      gzip=gzip,
                                                      failed=FALSE)
  }
  
  cat("\nAll samples were successfully paired")
  return(output.sampleList)
}
