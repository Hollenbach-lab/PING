# ----- Generate results directory -----
## Create the results directory if it does not exist
if(!file.exists(resultsDirectory)){
  dir.create(resultsDirectory, recursive = T)
}

## This function checks to make sure the output of system2 is valid
check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}

# ----- Getting alignment tools ready -----
### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
bowtie2Build <- system2('which', c('bowtie2-build'), stdout=T, stderr=T)
check.system2_output(bowtie2Build, 'bowtie2-build not found')

### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
samtools <- system2('which', c('samtools'), stdout=T, stderr=T)
check.system2_output(samtools, 'samtools')

### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
bcftools <- system2('which', c('bcftools'), stdout=T, stderr=T)
check.system2_output(bcftools, 'bcftools')

## Check to make sure samtools is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')


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

## Returns boolean [copied from gc_functions.R]
#is_nuc <- function(chr){
#  nuc_list <- c('A', 'T', 'C', 'G', '.')
#  return(as.character(chr) %in% nuc_list)
#}

is_nodel_nuc <- function(chr){
  nuc_list <- c('A','T','C','G')
  return( as.character(chr) %in% nuc_list )
}

## Counts how many unique nucleotides are in a vector
num_unique_nuc <- function(vector_to_check){
  return(sum(is_nuc(names(table(vector_to_check)))))
}

## Condenses allele fram into only positions that have multiple nucleotides
make_unique_pos_frame <- function(pos_frame){
  num_unique <- lapply(pos_frame, num_unique_nuc)
  unique_pos_frame <- pos_frame[,names(num_unique[num_unique > 1]), drop=F]
  return(unique_pos_frame)
}

kir.locus.vect <- c("KIR3DP1","KIR2DS5","KIR2DL3","KIR2DP1","KIR2DS3","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL1","KIR3DS1","KIR2DL2","KIR3DL2","KIR2DS4","KIR2DL1","KIR2DS1","KIR2DL5")

kirLocusFeatureNameList <- list()
kirLocusFeatureNameList[['KIR2DL1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL2']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL3']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL4']] <- c("5UTR","E1","I1","E2","I2","E3","I3/4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL5']] <- c("5UTR","E1","I1","E2","I2","E3","I3/4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DP1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS2']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS3']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS4']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS5']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL2']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL3']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5/6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DP1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","3UTR")
kirLocusFeatureNameList[['KIR3DS1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")

## This function checks to make sure the output of system2 is valid
check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}

# Write allele to FASTA file
general.write_fasta <- function(fastaCon, alleleName, alleleSeq){
  writeLines(paste0('>',alleleName,'\n',alleleSeq), fastaCon)
}


## This function returns the list of alleles found in the reference fasta
general.read_fasta <- function(fasta_path){
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  output.alleleList <- list()
  
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      currentAllele <- strsplit(currentLine, '>',fixed=TRUE)[[1]][2]
      output.alleleList[[ currentAllele ]] <- ''
    }else{
      output.alleleList[[ currentAllele ]] <- paste0( output.alleleList[[currentAllele]], currentLine )
    }
  }
  
  return(output.alleleList)
}



# Write feature coordinates to BED file
general.write_bed <- function(bedCon, alleleName, startPos, endPos, featureName,utr5ExtraLen=0,utr3ExtraLen=0){
  
  if(featureName == '5UTR'){
    writeLines(paste(alleleName,startPos,(endPos+utr5ExtraLen),featureName,sep='\t'), bedCon)
  }else if(featureName == '3UTR'){
    writeLines(paste(alleleName,(startPos+utr5ExtraLen),(endPos+utr5ExtraLen+utr3ExtraLen),featureName,sep='\t'),bedCon)
  }else{
    writeLines(paste(alleleName,(startPos+utr5ExtraLen),(endPos+utr5ExtraLen),featureName,sep='\t'),bedCon)
  }
}

## Finds samples in sampleDirectory and turns them into objects for downstream use
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

## Initializes a locus reference object, which stores reference allele information
general.initialize_locus_ref_object <- function(kirLocusVect = kir.locus.vect){
  
  ## Creating the sample object class
  locusRef <- setRefClass("locusRef",
                          fields=list(name='character',
                                      rawPath='character',
                                      formattedPath='character',
                                      filledPath='character',
                                      alleleBedList='list',
                                      alleleStrList='list'))
  
  ## Initialize locus ref list
  locusRefList <- list()
  
  ## Add locus names to ref list
  for(kirLocus in kirLocusVect){
    locusRefList[[kirLocus]] <- locusRef(name=kirLocus)
  }
  
  return(locusRefList)
}

## Read MSF files from copied_msf directory into a list of locus reference objects
initLocusRef.read_raw_msf <- function(locusRefList, copiedMsfDirectory, kirLocusVect = kir.locus.vect){
  
  ## Make sure the MSF directory is accessible
  copiedMsfDirectory <- normalizePath(copiedMsfDirectory,mustWork = T)
  
  cat('Checking',copiedMsfDirectory,'for IPD-KIR MSF files\n')
  
  ## List all of the correctly named raw.msf files
  rawMSFfileVect <- list.files(path=copiedMsfDirectory,pattern='raw.msf')
  
  ## Match the MSF files to the kir.locus.vect
  locusMatchedRawFileList <- sapply(names(locusRefList), function(x) return(c(rawMSFfileVect[grepl(x,rawMSFfileVect,fixed=T)])))
  
  ## Storing ref info
  for(locusName in names(locusRefList)){
    locusRefList[[locusName]]$rawPath <- file.path(copiedMsfDirectory, locusMatchedRawFileList[[locusName]])
  }
  
  cat('\nReference files found: ')
  cat(kirLocusVect[kirLocusVect %in% names(locusRefList)])
  
  cat('\nReference files not found: ')
  cat(kirLocusVect[!(kirLocusVect %in% names(locusRefList))])
  
  cat('\n\nPulling allele sequence from IPD-KIR reference files\n')
  cat('Processing:')
  
  for(locusRef in locusRefList){
    cat('',locusRef$name)
    
    ## Formatting locus name to match IPD-KIR
    altLocusNameStr <- strsplit(locusRef$name,'KIR',fixed=T)[[1]][2]
    
    ## Initialzing list to store allele strings
    locusAlleleStrList <- list()
    
    con = file(locusRef$rawPath, "r")
    while(TRUE){
      line = readLines(con, n = 1)
      
      ## Break out of the loop when end of file is reached
      if (length(line) == 0){
        break
      }
      
      ## Split the line by spaces
      lineVect <- unlist(strsplit(line,' ',fixed=T))
      
      ## Remove all empty elements
      lineVect <- lineVect[lineVect != '']
      
      ## If this is a sequence line, then process
      if(any(grepl(altLocusNameStr, lineVect)) & length(lineVect) > 1){
        
        ## Pull out the allele name, format for PING
        alleleNameStr <- paste0('KIR',lineVect[grepl(altLocusNameStr, lineVect)])
        
        ## Pull out the allele sequence
        alleleSeqStr <- paste0(lineVect[!grepl(altLocusNameStr, lineVect)],collapse='')
        
        ## Initialize a list element for this allele if it does not already exist
        if( !(alleleNameStr %in% names(locusAlleleStrList)) ){
          locusAlleleStrList[[alleleNameStr]] <- ''
        }
        
        ## Add the allele sequence to the allele string list
        locusAlleleStrList[[alleleNameStr]] <- paste0(locusAlleleStrList[[alleleNameStr]],alleleSeqStr)
      }
      
    }
    
    close(con)
    
    locusRef$alleleStrList <- locusAlleleStrList
  }
  
  return(locusRefList)
}

## Determine reference feature coordinates, write to BED file and store in locus reference objects
initLocusRef.create_bed <- function(locusRefList, referenceResourceDirectory, kirLocusFeatureNameList, writeBed=T){
  
  cat('\n\nWriting allele sequence as PING reference files\n')
  
  if(writeBed){
    kirBedPath <- file.path(referenceResourceDirectory, 'allKIR.bed')
    con = file(kirBedPath, "w")
  }
  
  cat('Processing:')
  
  for(locusRef in locusRefList){
    cat('',locusRef$name)
    
    locusBedList <- list()
    
    for( alleleNameStr in names(locusRef$alleleStrList) ){
      alleleStr <- locusRef$alleleStrList[[alleleNameStr]]
      
      ## Initialize list for storing BED information
      locusBedList[[alleleNameStr]] <- list()
      
      ## Pull out intron/exon sequences
      alleleFeatSeqVect <- strsplit(alleleStr,'|',fixed=T)[[1]]
      
      ## Save allele seq without boundary symbols back to alleleStrList
      locusRef$alleleStrList[[alleleNameStr]] <- paste0(unlist(alleleFeatSeqVect), collapse='')
      
      ## Pull out the locus feature names
      featureNameVect <- kirLocusFeatureNameList[[locusRef$name]]
      
      if( length(alleleFeatSeqVect) != length(featureNameVect) ){
        stop('Mismatched features for',locusRef$name)
      }
      
      ## Process each gene feature for this allele
      preBoundaryInt <- 0
      for( iterInt in 1:length(alleleFeatSeqVect) ){
        
        ## Pull out information for this iteration
        featSeqStr <- alleleFeatSeqVect[iterInt]
        boundaryInt <- nchar(featSeqStr)+preBoundaryInt
        delCount <- str_count(featSeqStr,pattern = fixed('.'))
        boundaryInt <- boundaryInt - delCount
        featNameStr <- featureNameVect[iterInt]
        
        delIndex <- grep('.',unlist(strsplit(featSeqStr,'')),fixed=T)
        
        noDelFeatSeqStr <- gsub('.','',featSeqStr,fixed=T)
        
        snpVect <- strsplit(featSeqStr,'',fixed=T)[[1]]
        names(snpVect) <- paste0(featNameStr,'_',1:length(snpVect))
        
        #### BED coordinates should be preBoundaryInt : boundaryInt
        ## Save gene feature coordinates to list
        locusBedList[[alleleNameStr]][[featNameStr]] <- list(alleleName=alleleNameStr,
                                                             startPos=preBoundaryInt,
                                                             endPos=boundaryInt,
                                                             featName=featNameStr,
                                                             featSeq=noDelFeatSeqStr,
                                                             featDelIndex=delIndex,
                                                             snpVect=snpVect)
        
        if(writeBed){
          ## Write feature coordinates to BED file
          general.write_bed(con, alleleNameStr, preBoundaryInt, boundaryInt, featNameStr)
          #writeLines(paste(alleleNameStr,preBoundaryInt,boundaryInt,featNameStr,sep='\t'), con)
        }
        
        ## Reset previous coordinate
        preBoundaryInt <- boundaryInt
      }
    }
    
    locusRef$alleleBedList <- locusBedList
  }
  
  if(writeBed){
    close(con)
  }
  cat('\n\nProcessing complete.')
  return(locusRefList)
}









