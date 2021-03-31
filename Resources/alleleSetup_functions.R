cigarOperationVect <- c('M','I','D')

## Load in the allele SNP dataframes
#annotatedAlleleDirectory <- file.path('/home/LAB_PROJECTS/PING2_PAPER/4_synSeq_dp50rl150_results/','alleleFiles')
annotatedAlleleDirectory <- file.path('Resources/genotype_resources/extended_SNP_files')
setup.knownSnpDFList <- list()
alleleSeq.list <- list()
for(locus in kirLocusList){
  alleleSnpDF <- read.csv( normalizePath(
    file.path(annotatedAlleleDirectory,
              paste0(locus, '_alleleSNPs.csv')
    )), check.names=F, row.names=1, stringsAsFactors = F,colClasses = 'character')
  
  alleleSeq.list <- c(alleleSeq.list, apply(alleleSnpDF,1, function(x) {
    seqX <- paste0(x,collapse='')
    seqX <- gsub('.','',seqX,fixed=T)
    seqX <- gsub('*','N',seqX,fixed=T)
    return(seqX)
  }) )
  
  setup.knownSnpDFList[[locus]] <- alleleSnpDF#[,grep('UTR',colnames(alleleSnpDF),invert = T,value = T)]
}

alleleSetupRef.df <- read.table('Resources/genotype_resources/gc_allele_reference.csv',stringsAsFactors = F, sep=',', check.names=F,row.names=1)

compact.alleleSeq.list <- alleleSeq.list[ names(alleleSeq.list) %in% unlist(alleleSetupRef.df) ]

delIndex.list <- list()
for( locus in names(setup.knownSnpDFList) ){
  snpDF <- setup.knownSnpDFList[[locus]]
  
  if(locus == 'KIR2DS1'){
    delIndex.list[[locus]] <- sapply(rownames(snpDF), function(x) {
      which(snpDF[x,] == '.')
    })
  }else{
    delIndex.list[[locus]] <- apply(snpDF, 1, function(x) {
      which(x == '.')
    })
  }
}

alleleSetupDirectory <- file.path(resultsDirectory,'allele_setup_files')
if(!file.exists(alleleSetupDirectory)){
  dir.create(alleleSetupDirectory)
}

KIR2DL4.10A.vect <- rownames(setup.knownSnpDFList$KIR2DL4)[ setup.knownSnpDFList$KIR2DL4[,'E7_105',drop=F] == 'A' ]
KIR2DL4.9A.vect <- rownames(setup.knownSnpDFList$KIR2DL4)[ setup.knownSnpDFList$KIR2DL4[,'E7_105',drop=F] != 'A' ]


# ----- PING allele setup -----
alleleSetup.write_fasta <- function( alleleSeq.list, fastaPath){
  i <- 1
  for( alleleID in names(alleleSeq.list) ){
    if( i == 1){
      cat(paste0('>',alleleID,'\n',alleleSeq.list[[alleleID]],'\n'),file=fastaPath,append=F)
    }else{
      cat(paste0('>',alleleID,'\n',alleleSeq.list[[alleleID]],'\n'),file=fastaPath,append=T)
    }
    i <- i+1
  }
  #cat('\n',file=fastaPath,append=T)
}
alleleSetup.write_sample_gc_fasta <- function( currentSample, alleleSetupDirectory, alleleSeq.list ){
  
  if( any( currentSample$copyNumber == 'failed' ) | any( is.na(currentSample$copyNumber) ) ){
    currentSample[['ASFastaPath']] <- 'failed'
    return(currentSample)
  }
  currentSample[['ASFastaPath']] <- file.path(alleleSetupDirectory,paste0(currentSample$name,'.fasta'))
  
  presentLociVect <- names(currentSample$copyNumber)[as.integer(currentSample$copyNumber) > 0]
  
  presentBoolVect <- sapply( names(alleleSeq.list), function(x){
    locus <- tstrsplit(x,'*',fixed=T)[[1]]
    if( 'KIR2DL5A' == locus | 'KIR2DL5B' == locus ){
      locus <- 'KIR2DL5'
    }
    
    return(locus %in% presentLociVect)
  })
  
  alleleSetup.write_fasta(alleleSeq.list[presentBoolVect], currentSample[['ASFastaPath']])
  
  return(currentSample)
}
alleleSetup.write_bt2_index <- function( currentSample, alleleSetupDirectory, bowtie2Build, threads ){
  
  if( currentSample[['ASFastaPath']] == 'failed' ){
    currentSample[['ASIndexPath']] <- 'failed'
    return(currentSample)
  }
  
  currentSample[['ASIndexPath']] <- file.path(alleleSetupDirectory,paste0(currentSample$name,'.alleleSetupRef'))
  
  ## Creqte a bowtie2 index for the kir_reference.fasta file <- only needed when building a new reference index
  createIndex <- system2(bowtie2Build, c(currentSample[['ASFastaPath']], 
                                         currentSample[['ASIndexPath']],'--quiet',paste('--threads', threads)))
  check.system2_output(createIndex, 'bowtie2 index building failed')
  
  return(currentSample)
}
alleleSetup.bt2_align <- function( currentSample, alleleSetupDirectory, bowtie2, threads ){
  
  if( currentSample[['ASIndexPath']] == 'failed' ){
    currentSample[['ASSamPath']] <- 'failed'
    return(currentSample)
  }
  currentSample[['ASSamPath']] <- file.path(alleleSetupDirectory,paste0(currentSample$name,'.sam'))
  
  ### 1. Align KIR extracted reads to haplo-reference
  #bt2_p <- paste0("-p", threads)
  #bt2_5 <- "-5 3"
  #bt2_3 <- "-3 7"
  #bt2_i <- "-i S,1,0.5"
  #bt2_min_score <- "--score-min L,0,-0.187"
  #bt2_I <- "-I 75"
  #bt2_X <- "-X 1500"
  #bt2_noUnal <- '--no-unal'
  
  #bt2_x <- paste0("-x ", currentSample[['ASIndexPath']] )
  #bt2_1 <- paste0('-1 ',currentSample$kirfastq1path)
  #bt2_2 <- paste0('-2 ',currentSample$kirfastq2path)
  
  #bt2_stream <- paste0("-S ", currentSample[['ASSamPath']] )    
  #bt2_al_conc <- paste0("--al-conc-gz ", fastqBase, "_%.fastq.gz")
  #bt2_un <- "--un dump.me"
  
  #optionsCommand <- c(bt2_p, bt2_5, bt2_3, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_noUnal, bt2_stream, bt2_al_conc, bt2_un)
  
  
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',currentSample[['ASIndexPath']]),
                      '-5 0', '-3 6', '-N 0', '--end-to-end', paste0('-p ',threads), '--score-min L,-2,-0.08',
                      '-I 75', '-X 1000',
                      paste0('-1 ',currentSample$kirfastq1path),
                      paste0('-2 ',currentSample$kirfastq2path),
                      '--no-unal',
                      '-a','--np 1', '--mp 2,2', '--rdg 1,1', '--rfg 1,1',
                      paste0('-S ',currentSample[['ASSamPath']]))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 failed, retrying alignment...')
    output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  }
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 failed, marking sequence as failure and moving on')
    currentSample[['ASSamPath']] <- 'failed'
    return(currentSample)
  }
  
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  currentSample[['ASSamPath']] <- normalizePath(currentSample[['ASSamPath']], mustWork=T)
  
  cat('\n\nSuccessfully aligned',currentSample$name)
  
  return(currentSample)
}
## This function runs a samtools sam to bam conversion
alleleSetup.sam_to_bam <- function(samtools_command, currentSample, bamDirectory, threads){
  
  if( currentSample[['ASSamPath']] == 'failed' ){
    currentSample[['ASBamPath']] <- 'failed'
    return(currentSample)
  }
  
  ## Initialize an output path for the BAM file
  currentSample[['ASBamPath']] <- file.path(bamDirectory,paste0(currentSample$name,'.bam'))
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),
                      currentSample$ASSamPath, '-o', currentSample[['ASBamPath']])
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamConv, 'samtools sam to bam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['bamPath']] <- normalizePath(currentSample[['ASBamPath']], mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample[['ASSamPath']],'to',currentSample[['ASBamPath']])
  
  return(currentSample)
}
## This function runs a samtools bam to sam conversion
alleleSetup.bam_to_sam <- function(samtools_command, currentSample, bamDirectory, threads){
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),'-h',
                      currentSample$ASBamPath, '-o', currentSample$ASSamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.samConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.samConv, 'samtools bam to sam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.samConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['ASSamPath']] <- normalizePath(currentSample$ASSamPath, mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample$bamPath,'to',currentSample$ASSamPath)
  
  return(currentSample)
}
## Cleanup function
alleleSetup.remove_index <- function( currentSample, alleleSetupDirectory ){
  removeVect <- list.files( alleleSetupDirectory, paste0(currentSample$name, '.alleleSetupRef'), full.names = T )
  file.remove( removeVect )
  return(NULL)
}

alleleSetup.gc_matched_ref_alignment <- function( currentSample, alleleSetupDirectory, alleleSeq.list, threads ){
  
  currentSample <- alleleSetup.write_sample_gc_fasta( currentSample, alleleSetupDirectory, alleleSeq.list )
  
  currentSample <- alleleSetup.write_bt2_index( currentSample, alleleSetupDirectory, bowtie2Build, threads )
  
  currentSample <- alleleSetup.bt2_align( currentSample, alleleSetupDirectory, bowtie2, threads )
  
  alleleSetup.remove_index( currentSample, alleleSetupDirectory )
  
  currentSample <- alleleSetup.sam_to_bam(samtools, currentSample, alleleSetupDirectory, threads)
  
  return( currentSample )
}

alleleSetup.del_start_offset <- function( refPos, currentDelIndex, delRegionIndexVect ){
  #cat('\n',refPos,'',length(currentDelIndex))
  
  if(length( delRegionIndexVect ) == 0 ){
    delRegionIndexVect <- length(currentDelIndex)
  }
  
  if(length(currentDelIndex) == 0){
    return(refPos)
  }else{
    
    #oldRefPos <- refPos
    delOffsetBool <- any( refPos >= currentDelIndex[1:delRegionIndexVect[1]] )
    #delOffset <- length(currentDelIndex)
    
    if( !delOffsetBool ){
      return(refPos)
    }
    
    delOffset <- length(currentDelIndex[1:delRegionIndexVect[1]])
    refPos <- refPos + delOffset
    
    if(length( currentDelIndex ) > delOffset ){
      currentDelIndex <- currentDelIndex[ (delOffset+1):length(currentDelIndex) ]
      delRegionIndexVect <- delRegionIndexVect - delOffset
      return( alleleSetup.del_start_offset( refPos, currentDelIndex, delRegionIndexVect[-1]) )
    }else{
      return(refPos)
    }
    
  }
  
}

alleleSetup.del_end_offset <- function( startPos, seqLength, endDelIndex, delRegionIndexVect ){
  #cat('\n',startPos,'',seqLength,'',length(endDelIndex))
  endPos <- startPos + seqLength
  
  if(length( delRegionIndexVect ) == 0 ){
    delRegionIndexVect <- length(endDelIndex)
  }
  
  if(length(endDelIndex) == 0){
    return( seqLength )
  }else{
    
    delOffsetBool <- any( endPos >= endDelIndex )
    
    if( !delOffsetBool ){
      return( seqLength )
    }
    
    delOffset <- length( endDelIndex[1:delRegionIndexVect[1]] )
    
    seqLength <- seqLength + delOffset
    
    if( length(endDelIndex) > delOffset ){
      endDelIndex <- endDelIndex[ (delOffset+1):length(endDelIndex) ]
      delRegionIndexVect <- delRegionIndexVect - delOffset
      return( alleleSetup.del_end_offset( startPos, seqLength, endDelIndex, delRegionIndexVect[-1] ) )
    }else{
      return( seqLength )
    }
  }
}

## Format reads for snp phasing (fill in deletion positions with '.')
alleleSetup.read_formatter <- function( readName, currentSeq, refPos, currentLocus, currentRef, cigarStr ){
  
  #cat('\n',readName,currentLocus,currentRef)
  cigarMod.del.bool <- grepl('D', cigarStr)
  cigarMod.ins.bool <- grepl('I', cigarStr)
  
  if( cigarMod.del.bool | cigarMod.ins.bool ){
    resultList <- samFormat.processCigarStr( cigarStr, currentSeq )
    currentSeq <- resultList$readSeq
    insIndex <- resultList$insIndex
  }
  
  currentDelIndex <- delIndex.list[[currentLocus]][[currentRef]]
  
  #cat('',length(currentDelIndex))
  
  seqLength <- nchar(currentSeq)
  delBool <- length(currentDelIndex) > 0
  
  delRegionIndexVect <- which( diff( currentDelIndex ) > 1 )
  
  startPos <- alleleSetup.del_start_offset( refPos, currentDelIndex, delRegionIndexVect )
  #cat('','passed startPos')
  ## Format the start position
  #startPos <- as.numeric(refPos)
  
  ## Format the end position
  endDelIndex <- currentDelIndex[ currentDelIndex > startPos ]
  delRegionIndexVect <- which( diff( endDelIndex ) > 1 )
  endPos <- startPos + alleleSetup.del_end_offset( startPos, seqLength, endDelIndex, delRegionIndexVect ) - 1
  
  ## Testing if there deletions found in the reference for the current read
  if( (endPos-startPos + 1) > seqLength ){
    
    delIndexList <- which(startPos:endPos %in% currentDelIndex)
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    for(delIndex in delIndexList){
      seqLength <- nchar(currentSeq)
      preString <- substr(currentSeq, 1, delIndex-1)
      postString <- substr(currentSeq, delIndex, seqLength)
      
      currentSeq <- paste0(preString, '.', postString)
    }
  }
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- startPos:endPos
  
  ## Turn the sequence string into a list
  seqList <- strsplit(currentSeq,'')[[1]]
  
  if(length(seqList) != length(fullIndex)){
    cat('\n',startPos,endPos,currentRef,length(seqList),length(fullIndex))
  }
  
  names(seqList) <- fullIndex
  
  if( cigarMod.ins.bool ){
    for( insPos in names( insIndex ) ){
      insPosInt <- as.integer(insPos)
      #cat('\n',readName)
      # Mod to assign insertion sequence over subsequet del characters
      if( seqList[ (as.integer(insPos)+1) ] == '.' ){
        
        insSeqVect <- strsplit( insIndex[[insPos]], '')[[1]]
        insSeqVect <- insSeqVect[2:length(insSeqVect)]
        insSeqPosVect <- (insPosInt+1):(length(insSeqVect)+insPosInt)
        
        readDelIndex <- which( seqList == '.' )
        
        assign.index <- insSeqPosVect %in% readDelIndex
        
        if( any(assign.index) ){
          seqList[ insSeqPosVect[ assign.index ] ] <- insSeqVect[ assign.index ]
        }
        
        # If there is leftover insertion sequence, append it to the last assigned ins character
        if( any(!assign.index) ){
          lastInsPos <- insSeqPosVect[ assign.index ][ sum( assign.index ) ]
          seqList[ lastInsPos ] <- paste0( c(seqList[ lastInsPos ], insSeqVect[ !assign.index ]), collapse='')
        }
      }else{
        seqList[ insPosInt ] <- insIndex[[insPos]]
      }
    }
  }
  
  seqListList <- list( list(seqList) )
  return(seqListList)
}

alleleSetup.process_samDT <- function( samPath, delIndex.list, processSharedReads, readBoost.thresh){
  
  # Read in SAM file, then delete SAM to save space
  nHead <- general.count_sam_header( samPath )
  samTable <- general.read_sam( samPath , rows_to_skip = nHead )
  
  # Pull out reads that align to s single locus
  isUniqueDT <- samTable[, length( unique(locus) ) == 1, by=read_name ]
  
  uniqueReadVect <- isUniqueDT$read_name[isUniqueDT$V1]
  sharedReadVect <- isUniqueDT$read_name[!isUniqueDT$V1]
  
  if(processSharedReads){
    cat("\n\tRunning cross-mapped read separation with readBoost.thresh (alignment score buffer) =",readBoost.thresh)
    
    sharedSamTable <- samTable[read_name %in% sharedReadVect]
    sharedSamTable$asSum <- 0
    sharedSamTable[,asSum:= sum(alignment_score), by=c('read_name','reference_name')]
    
    passedThresh.index <- sharedSamTable[,asSum >= (max(asSum) - readBoost.thresh),by=c('read_name')]$V1
    passedThresh.sharedSamTable <- sharedSamTable[passedThresh.index,]
    
    unique.passedThresh.subDT <- passedThresh.sharedSamTable[,length(unique(locus)) == 1,by='read_name']
    unique.sharedReadVect <- unique.passedThresh.subDT$read_name[unique.passedThresh.subDT$V1]
    
    addSamTable <- passedThresh.sharedSamTable[ read_name %in% unique.sharedReadVect ]
    addSamTable <- unique(addSamTable, by=c('read_seq'))
    addSamTable[,asSum:=NULL]
    cat("\n\tRescued",nrow(addSamTable),'reads.')
  }
  
  uniqueSamDT <- copy( samTable[read_name %in% uniqueReadVect] )
  uniqueSamDT <- unique(uniqueSamDT, by=c('read_seq'))
  
  if(processSharedReads){
    #uniqueSamDT <- merge( uniqueSamDT, addSamTable, all=T, by=colnames(uniqueSamDT) )
    uniqueSamDT <- rbindlist( list(uniqueSamDT, addSamTable) )
    rm(addSamTable)
  }
  
  rm(samTable)
  
  cat('\nFormatting reads.')
  uniqueSamDT[, 'readTable' := alleleSetup.read_formatter(read_name, read_seq, ref_pos, locus, reference_name, cigarStr), by=seq_len(nrow(uniqueSamDT))]
  
  return(uniqueSamDT)
}

alleleSetup.generate_snpDF <- function( rawSnpTbl ){
  
  snpPosVect <- unique( tstrsplit( names(rawSnpTbl), '^', fixed=T )[[1]] )
  
  snpDF <- as.data.frame( matrix('',nrow=3,ncol=length(snpPosVect) ), stringsAsFactors = F )
  rownames(snpDF) <- paste0('SNP_',1:3)
  colnames(snpDF) <- snpPosVect
  
  for( snpPos in snpPosVect ){
    snpVect <- grep( paste0(snpPos, '^'), names(rawSnpTbl), fixed=T, value=T )
    snpVect <- tstrsplit(snpVect,'^',fixed=T)[[2]]
    
    if(length(snpVect) == 1){
      snpDF[1,snpPos] <- snpVect
      snpDF[2:3,snpPos] <- 'N'
    }else if(length(snpVect) == 2){
      snpDF[1:2,snpPos] <- snpVect
      snpDF[3,snpPos] <- 'N'
    }else if(length(snpVect) == 3){
      snpDF[,snpPos] <- snpVect
    }
  }
  
  return(snpDF)
  
}

alleleSetup.generate_rawSnpTbl <- function( currentLocus, currentLocusSnpDF, currentADSnpDF, uniqueSamDT, setup.minDP ){
  
  adCols <- colnames(currentADSnpDF)
  
  colIndexVect <- as.character( which( colnames(currentLocusSnpDF) %in% adCols ) )
  colIndexConv.list <- as.list(colnames(currentLocusSnpDF)[as.integer(colIndexVect)])
  names(colIndexConv.list) <- colIndexVect
  
  cat('\n\t\tGenerating list of aligned AD SNPs.')
  adSnpList <- sapply( uniqueSamDT[ locus == currentLocus ]$readTable, function(x){
    boolVect <- colIndexVect %in% names(x)
    if( any(boolVect) ){
      return( x[colIndexVect[boolVect]] )
    }else{
      return(NA)
    }
  })
  
  adSnpList <- adSnpList[!is.na(adSnpList)]
  
  cat('\n\t\tFormatting',length(adSnpList),'AD SNPs.')
  
  formattedADSnp.list <- sapply( adSnpList, function(x){
    sapply( names(x), function(y){
      snpNuc <- x[[y]]
      return( paste0(colIndexConv.list[[y]],'^',snpNuc) )
    })
  })
  formattedADSnp.list <- unlist( formattedADSnp.list, use.names = F )
  
  rawSnpTbl <- table(unlist(formattedADSnp.list))
  rawSnpTbl <- rawSnpTbl[rawSnpTbl>=setup.minDP]
  
  if(length( rawSnpTbl ) <= 1 ){
    cat('\n\t\tNo SNP positions passing setup.minDP =',paste0(setup.minDP,'.'),' Retrying with setup.minDP = 1.')
    rawSnpTbl <- table(unlist(formattedADSnp.list))
    rawSnpTbl <- rawSnpTbl[rawSnpTbl>=1]
    cat('\n\t\tNow',length(rawSnpTbl),'SNP positions are passing.')
  }
  
  return(rawSnpTbl)
}

alleleSetup.call_allele <- function( currentSample, currentLocus, currentSnpDF, currentADSnpDF, homScoreBuffer=1, kirRes=7, ambScore=F ){
  
  if( currentLocus == 'KIR3DL1S1' ){
    locusCopy <- as.integer(currentSample$copyNumber[['KIR3DS1']]) + as.integer(currentSample$copyNumber[['KIR3DL1']])
  }else if( currentLocus == 'KIR2DL23' ){
    locusCopy <- as.integer(currentSample$copyNumber[['KIR2DL2']]) + as.integer(currentSample$copyNumber[['KIR2DL3']])
  }else if( currentLocus == 'KIR2DS35' ){
    locusCopy <- as.integer(currentSample$copyNumber[['KIR2DS3']]) + as.integer(currentSample$copyNumber[['KIR2DS5']])
  }else{
    locusCopy <- as.integer( currentSample$copyNumber[[currentLocus]] )
  }
  
  locusCopy <- min( locusCopy, 3 )
  
  if( ncol(currentSnpDF) == 0 ){
    return(NULL)
  }
  
  # Attempt to cut down the rows of the currentSnpDF to save processing time
  maxVariants <- max( apply(currentSnpDF, 2, function(x){ num_unique_nuc(x) } ) )
  if( maxVariants < 3 ){
    currentSnpDF <- currentSnpDF[1:2,]
  }
  
  # Format sample DF and transform into data table
  currentSnpDT <- as.data.table(currentSnpDF)
  currentSnpDT$name <- rownames(currentSnpDF)
  setkey(currentSnpDT, name)
  
  # Subset by aligned SNP positions, remove duplicate alleles
  adSnpDF <- unique(currentADSnpDF[,colnames(currentSnpDF)])
  
  # Reset the adSnpDF back to input if no AD cols
  if( ncol(adSnpDF) == 0 ){
    adSnpDF <- currentADSnpDF[,colnames(currentSnpDF)]
  }
  
  '
  1.21677 mins (nameVect)
  1.19399 mins (boolVect)
  0.5976 secs  (matrix/vector comparison)
  '
  # Generate excluded allele list (this cuts down processing identical alleles, but must be added back in)
  tempADSnpDF <- currentADSnpDF[,colnames(currentSnpDF)]
  excludedAlleleVect <- setdiff( rownames(currentADSnpDF), rownames(adSnpDF) )
  excludedAlleleList <- sapply(excludedAlleleVect, function(x){
    nonAmbColVect <- unlist( tempADSnpDF[x,,drop=T] != '*' )
    
    adMat <- as.matrix( tempADSnpDF[,nonAmbColVect,drop=F] )
    names( which( apply( t(t(adMat) == adMat[x,]), 1, all ) ) )
  })
  
  
  adSnpDT <- as.data.table(adSnpDF)
  adSnpDT$alleleName <- rownames(adSnpDF)
  setkey(adSnpDT, alleleName)
  
  # Pull out AD positions
  adCols <- colnames(adSnpDF)
  
  # Identify het/hom positions (should be length 0 for homozygous)
  namedHetVect <- which( apply(currentSnpDF, 2, function(x){ num_unique_nuc(x) }) >= 2 )
  hetPosVect <- names(namedHetVect)
  homPosVect <- setdiff(colnames(currentSnpDF), hetPosVect)
  
  # Set if het allele calling should be enabled
  hetBool <- length(hetPosVect) > 0
  homBool <- length(homPosVect) > 0
  
  
  # Bool check is fix for no AD homozygous positions
  if( homBool ){
    cat('\n\t\tScoring hom positions.')
    
    scoreDF <- sapply( homPosVect, function(x){
      #cat('\n',x)
      outVect1 <- (adSnpDT[,..x] != currentSnpDT[1,..x][[1]])*1
      
      if( !ambScore ){
        ## 11/30/2020 Change to now penalize ambiguous characters (since we are using filled SNP DF matching)
        outVect2 <- (adSnpDT[,..x] != '*')*1
        finalOutVect <- pmin(outVect1, outVect2)
        names(finalOutVect) <- adSnpDT$alleleName
      }else{
        finalOutVect <- outVect1
        names(finalOutVect) <- adSnpDT$alleleName
      }
      
      return(finalOutVect)
    })
    
    homScoreList <- apply(scoreDF, 1, sum)
  }
  
  #homScoreBuffer <- 2
  if( !hetBool ){
    bestScoreInt <- min(homScoreList)
    #bestMatchIndex <- which( homScoreList == bestScoreInt )
    bestMatchAlleleVect <- names( which( homScoreList == bestScoreInt ) )
  }else{
    
    # Fix for no AD homozygous positions
    if( !homBool ){
      homAlleleVect <- adSnpDT$alleleName
    }else{
      #homAlleleVect <- names(homScoreList)[homScoreList <= ( min(homScoreList)+homScoreBuffer ) ] # Cut any alleles with more than min mismatches
      homAlleleVect <- names( which( homScoreList <= ( min(homScoreList)+homScoreBuffer ) ) )
    }
    
    possAllelePairMat <- combinations(length(homAlleleVect),locusCopy,homAlleleVect,repeats.allowed = T)
    
    # Set up possible allele data table
    possAlleleDT <- as.data.table(possAllelePairMat)
    possAlleleDT$alleleVect <- apply(possAllelePairMat, 1, list)
    
    # initialize column for storing distance scores
    possAlleleDT$distance <- sapply( possAlleleDT$alleleVect, function(x) sum( homScoreList[x[[1]]] ) )
    
    # Score het positions for allele pairs (this can take awhile if there are many pairs)
    cat('\n\t\tScoring',nrow(possAlleleDT),'allele pairings...')
    
    '
    alleleSNP (TCG) == alignedSNP (TCC)
    alignedSNP
    1 2 3
    allele1 1 T F F
    allele2 2 F T T
    allele3 3 F F F
    
    score <- max( sum(!rowBoolVect), sum(!colBoolVect) )
    '
    
    sapply( hetPosVect, function(currentPos){
      alignedSnpVect <- currentSnpDT[,..currentPos][[1]]
      possAlleleDT[,'distance' := alleleSetup.geno_score_calc( alleleVect, distance, adSnpDT, alignedSnpVect, currentPos, ambScore ), by=seq_len(nrow(possAlleleDT)) ]
      return(NULL)
    })
    
    bestScoreInt <- min(possAlleleDT$distance)
    bestMatchIndex <- which(possAlleleDT$distance == bestScoreInt)
    
    bestMatchAlleleVect <- sapply( unlist( possAlleleDT$alleleVect[bestMatchIndex], recursive=F ), paste0, collapse='+' )
  }
  
  ## Add back in exactly matching excluded alleles + cut allele calls to desired res
  formattedAlleleVect <- alleleSetup.format_call( bestMatchAlleleVect, excludedAlleleList, hetBool, kirRes )
  
  # Print out the best score and cooresponding allele matches
  allMatchingAlleleStr <- paste0(formattedAlleleVect, collapse=' ')
  cat('\n\t\t\tBest score:',bestScoreInt)
  cat('\n\t\t\tAllele match(es):',allMatchingAlleleStr)
  
  return(list('alleleVect'=formattedAlleleVect,'scoreInt'=bestScoreInt))
}

alleleSetup.geno_score_calc <- function( alleleVect, distance, adSnpDT, alignedSnpVect, currentPos, ambScore ){
  alleleVect <- alleleVect[[1]]
  
  genoSnpVect <- adSnpDT[alleleVect,..currentPos][[1]]
  
  snpMat <- t( sapply(genoSnpVect, function(x){x == alignedSnpVect}) )
  
  if( !ambScore ){
    ## 11/30/2020 change to penalize/accept ambiguous character matching
    ambMatch <- genoSnpVect == '*'
    if( any(ambMatch) ){
      snpMat[ambMatch,] <- TRUE
    }
  }
  
  ## 12/04/2020 change to accept 'N' matching
  nMatch <- alignedSnpVect == 'N'
  if( any( nMatch ) ){
    snpMat[,nMatch] <- TRUE
  }
  
  rowSum <- sum( !apply(snpMat, 1, any) )
  colSum <- sum( !apply(snpMat, 2, any) )
  score <- max( rowSum, colSum )
  
  return(distance + score)
}

alleleSetup.pair_score_calc <- function( allele1, allele2, distance, adSnpDT, sampleSnpDT, adCols ){
  
  scoreList <- sapply(adCols, function(curPos){
    
    pos1 <- adSnpDT[ c(allele1, allele2), ..curPos][1,]
    pos2 <- adSnpDT[ c(allele1, allele2), ..curPos][2,]
    
    con1 <- !all( sampleSnpDT[, ..curPos][1,] == pos1 )
    con2 <- !all( sampleSnpDT[, ..curPos][2,] == pos2 )
    
    con3 <- !all( sampleSnpDT[, ..curPos][1,] == pos2 )
    con4 <- !all( sampleSnpDT[, ..curPos][2,] == pos1 )
    
    if( pos1 == '*' ){
      con1=F
      con4=F
    }
    
    if( pos2 == '*' ){
      con2=F
      con3=F
    }
    
    return( min((con1*1+con2*1),(con3*1+con4*1)) )
  })
  
  distanceInt <- sum(scoreList)
  
  return( distanceInt + distance )
}

alleleSetup.format_call <- function( bestMatchAlleleVect, excludedAlleleList, hetBool, kirRes){
  output.alleleVect <- c()
  
  if( length(excludedAlleleList) == 0 ){
    excludedAlleleList[[1]] <- 'dummy'
  }
  
  if(hetBool == F){
    
    splitAlleleVect <- tstrsplit( bestMatchAlleleVect, ' ', fixed=T )[[1]]
    homAlleleVect <- c()
    
    for( alleleMatch in splitAlleleVect ){
      homAlleleVect <- c( alleleMatch, homAlleleVect )
      homAlleleVect <- c( homAlleleVect, names( which( sapply( excludedAlleleList, function(x) alleleMatch %in% x ) ) ) )
    }
    
    output.alleleVect <- unique( unlist( lapply( homAlleleVect, kir.allele_resolution, kirRes ) ) )
    #output.alleleVect <- unique( homAlleleVect )
  }else{
    
    splitAllelePairVect <- tstrsplit( bestMatchAlleleVect, ' ', fixed=T )[[1]]
    hetAlleleVect <- c()
    
    for(allelePair in splitAllelePairVect){
      
      splitAlleleVect <- strsplit( allelePair, '+', fixed=T )[[1]]
      
      ambList <- lapply(1:length(splitAlleleVect), function(x){
        alleleName <- splitAlleleVect[x]
        c(alleleName, names( which( sapply( excludedAlleleList, function(y) alleleName %in% y ) ) ) )
      })
      
      ambList <- lapply(ambList, function(x){
        unique( unlist(lapply(x, kir.allele_resolution, kirRes)) )
      })
      
      ambMat <- expand.grid(ambList)
      
      hetAlleleVect <- c( hetAlleleVect, apply( ambMat, 1, paste0, collapse='+' ) )
    }
    
    output.alleleVect <- unique( hetAlleleVect )
  }
  
  return( output.alleleVect )
}


alleleSetup.call_setup_alleles <- function( currentSample, uniqueSamDT, setup.knownSnpDFList, hetRatio, setup.minDP, ambScore=T, allPosScore=F, onlyExonScore=F, homScoreBuffer=4, includeAmb=F, addDiversity=F, combineTypings=F, skipSnpGen=F, skipSet=F){
  
  '
  pingAllele.generate_snp_df could be sped up by only processing AD SNP positions, would save time and storage without impacting results
  '
  if( !skipSnpGen ){
    currentSample <- pingAllele.generate_snp_df( currentSample, uniqueSamDT, currentSample[['iterRefDirectory']], setup.knownSnpDFList, 'setup', hetRatio, setup.minDP )
  }

  
  # ----- CALL ALLELES -----
  presentLociVect <- names(currentSample$copyNumber)[currentSample$copyNumber > 0]
  currentAlleleCall.list <- list()
  currentScore.list <- list()
  setCall.list <- currentSample[['setCallList']]
  
  cat('\nProcessing:')
  for(currentLocus in presentLociVect){
    cat('\n\n\t',currentLocus)
    
    currentLocusSnpDF <- setup.knownSnpDFList[[currentLocus]]
    
    if( currentLocus %in% names(currentSample[['setCallList']]) & !skipSet ){
      #currentLocusSnpDF <- currentLocusSnpDF[ currentSample[['setCallList' ]], ]
      
      currentAlleleCall.list[[currentLocus]] <- currentSample[['setCallList']][[currentLocus]]
      currentScore.list[[currentLocus]] <- 0
      next
    }
    
    if( allPosScore ){
      currentADSnpDF <- currentLocusSnpDF
      callRes <- 7
    }else if(onlyExonScore){
      exonNameVect <- grep('E',kirLocusFeatureNameList[[currentLocus]],value=T)
      exonNameVect <- grep('PE',exonNameVect,value=T,invert = T)
      sampleCols <- tstrsplit(colnames(currentLocusSnpDF),'_',fixed=T)[[1]] 
      sampleCols <- colnames(currentLocusSnpDF)[ sampleCols %in% exonNameVect ]
      currentADSnpDF <- currentLocusSnpDF[,sampleCols]
      callRes <- 5
    }else{
      currentADSnpDF <- make_unique_pos_frame( currentLocusSnpDF )
      callRes <- 7
    }
    
    # rawSnpTbl <- alleleSetup.generate_rawSnpTbl( currentLocus , currentLocusSnpDF, currentADSnpDF, uniqueSamDT, setup.minDP )
    # currentSnpDF <- alleleSetup.generate_snpDF( rawSnpTbl )
    
    currentSnpDF <- read.csv(file.path(currentSample$snpDFPathList[['setup']][['SNP']][[currentLocus]]),
                             check.names=F,stringsAsFactors = F,row.names=1,header = T,colClasses = c("character"))
    
    
    #nonADCol.vect <- colnames(currentSnpDF)[!(colnames(currentSnpDF) %in% colnames(currentADSnpDF))]
    
    #dropCol.vect <- names( which( apply( currentSnpDF[,nonADCol.vect], 2, function(x) num_unique_nuc(x) > 1) ) )
    keepCol.vect <- colnames(currentSnpDF)[(colnames(currentSnpDF) %in% colnames(currentADSnpDF))]
    
    if( length(keepCol.vect) < 2 ){
      cat('\n\t\tFewer than 2 positions for calling, setting call to all alleles for this locus')
      currentAlleleCall.list[[currentLocus]] <- rownames( setup.knownSnpDFList[[currentLocus]] )
      currentScore.list[[currentLocus]] <- ncol(currentADSnpDF)
      next
    }
    
    currentSnpDF <- currentSnpDF[,keepCol.vect,drop=F]
    
    totalCols <- ncol(currentADSnpDF)
    presentCols <- sum( colnames(currentADSnpDF) %in% colnames(currentSnpDF) )
    posPerc <- paste0( round( presentCols / totalCols, 3)*100,'%')
    
    if( allPosScore ){
      cat('\n\t\tCalling on',posPerc,'of all positions')
    }else if( onlyExonScore ){
      cat('\n\t\tCalling on',posPerc,'of exon positions')
    }else{
      cat('\n\t\tCalling on',posPerc,'of AD positions')
    }
    
    callList <- alleleSetup.call_allele( currentSample, currentLocus, currentSnpDF, currentADSnpDF, homScoreBuffer=homScoreBuffer, kirRes=callRes, ambScore=ambScore )
    
    if( onlyExonScore ){
      out.upres.vect <- c()
      for( allelePair in callList$alleleVect ){
        compAllele.vect <- strsplit( allelePair, '+' ,fixed=T)[[1]]
        compAllele.vect <- as.vector( sapply(compAllele.vect, function(x) general.kir_upres(x,setup.knownSnpDFList)[[1]]  ) )
        out.upres.vect <- c(out.upres.vect, paste0(compAllele.vect, collapse='+'))
      }
      callList$alleleVect <- out.upres.vect
    }
    
    if( combineTypings ){
      currentAlleleCall.list[[currentLocus]] <- unique( c(currentSample$setupAlleleList[[currentLocus]], callList$alleleVect) )
      currentScore.list[[currentLocus]] <- callList$scoreInt
    }else{
      currentAlleleCall.list[[currentLocus]] <- callList$alleleVect
      currentScore.list[[currentLocus]] <- callList$scoreInt
    }
  }
  
  # ----- REF ALLELE SELECTION -----
  '
  Reference alleles are chosen based on which ambiguity has the most defined positions
  '
  
  lockedLocus.vect <- intersect( names( which( currentScore.list == 0 ) ), names( which( sapply(currentAlleleCall.list, length ) == 1) ) )
  if( length(lockedLocus.vect) > 0 ){
    addLocus.vect <- setdiff( lockedLocus.vect, names(setCall.list) )

    for( locus in addLocus.vect ){
      setCall.list[[locus]] <- currentAlleleCall.list[[locus]]
    }
  }
  
  setCall.locusVect <- intersect( names(setCall.list), names(currentAlleleCall.list) )
  for( locus in setCall.locusVect ){

    currentCall.alleleVect <- currentAlleleCall.list[[locus]]
    setCall.alleleVect <- setCall.list[[locus]]

    if( ! (all(currentCall.alleleVect %in% setCall.alleleVect) & all(setCall.alleleVect %in% currentCall.alleleVect)) ){
      
      if( length(currentCall.alleleVect) == 1){
        cat('\n\n-- Adding',paste0(setCall.alleleVect,collapse=' '),'to',paste0(currentCall.alleleVect,collapse=' '),'---')
        
        currentCall.alleleVect <- unique( unlist( strsplit( currentCall.alleleVect, '+', fixed=T) ) )
        setCall.alleleVect <- unique( unlist( strsplit( setCall.alleleVect, '+', fixed=T) ) )
        
        currentAlleleCall.list[[locus]] <- paste0( unique(c(setCall.alleleVect,currentCall.alleleVect)), collapse='+' )
        setCall.list[[locus]] <- paste0( unique(c(setCall.alleleVect,currentCall.alleleVect)), collapse='+' )
      }
      # else{
      #   cat('\n\n-- Replacing',paste0(currentCall.alleleVect,collapse=' '),'with',paste0(setCall.alleleVect,collapse=' '),'---')
      #   currentAlleleCall.list[[locus]] <- setCall.list[[locus]]
      # }
    }
  }
  
  # lockedLocus.vect <- names( which( currentScore.list == 0 ) )
  # if( length(lockedLocus.vect) > 0 ){
  #   
  #   for( locus in lockedLocus.vect){
  #     
  #     if( !(locus %in% names(setCall.list)) ){
  #       setCall.list[[locus]] <- currentAlleleCall.list[[locus]]
  #     }else{
  #       
  #       setCall.vect <- setCall.list[[locus]] 
  #       currentCall.vect <- currentAlleleCall.list[[locus]]
  #       
  #       if( length(setCall.vect) > length(currentCall.vect) ){
  #         setCall.list[[locus]] <- currentCall.vect
  #       }
  #     }
  #   }
  # }
  
  if( includeAmb ){
    cat('\nIncluding all genotype ambiguity in reference.')
    nonAmbAlleleCall.list <- currentAlleleCall.list
  }else{
    nonAmbAlleleCall.list <- alleleSetup.most_defined_typing( currentAlleleCall.list )
  }
  
  if( addDiversity ){
    
    addDiversityLoci.vect <- names( which( currentScore.list > 0 ) )
    cat('\nAdding diverse allele set to:')
    for( addLocus in addDiversityLoci.vect ){
      cat('',addLocus)
      nonAmbAlleleCall.list[[addLocus]] <- paste0( alleleSetupRef.df[addLocus,], collapse='+' )
    }
  }
  # cat('\n\nReducing genotype ambiguity by selecting most defined sequences')
  # nonAmbAlleleCall.list <- sapply( names(currentAlleleCall.list), function(x){
  #   alleleCallVect <- currentAlleleCall.list[[x]]
  #   
  #   if(length(alleleCallVect) == 1){
  #     return(alleleCallVect)
  #   }
  #   
  #   alleleCallList <- strsplit(alleleCallVect, '+', fixed=T)
  #   
  #   nonAmbSumVect <- sapply( alleleCallList, function(y){
  #     sum( sapply( y, function(z){
  #       sum( setup.knownSnpDFList[[x]][z,] != '*' )
  #     }) )
  #   })
  #   
  #   selectedCall <- which( nonAmbSumVect == max( nonAmbSumVect ) )[1]
  #   return( alleleCallVect[selectedCall] )
  # })
  

  if( 'KIR2DL1' %in% names(nonAmbAlleleCall.list) ){
    pos4710 <- as.integer( currentSample$kffHits[['*KIR2DL1*4710b']] ) > 10
    pos47 <- as.integer( currentSample$kffHits[['*KIR2DL1*4and7']] ) > 10
    pos7 <- as.integer( currentSample$kffHits[['*KIR2DL1*7only']] ) > 10
    pos003 <- as.integer( currentSample$kffHits[['*KIR2DL1*003']] ) > 10
    
    calledAlleleVect <- unique( unlist( strsplit( nonAmbAlleleCall.list[['KIR2DL1']], '+', fixed=T) ) )
    
    if( pos4710 & !pos47 & !pos7 ){
      addAllele <- 'KIR2DL1*010'
      cat('\n\nAdjusting KIR2DL1 call to include',addAllele,'based on KFF')
    }else if( pos4710 & pos47 & !pos7 ){
      addAllele <- 'KIR2DL1*0040101'
      cat('\n\nAdjusting KIR2DL1 call to include',addAllele,'based on KFF')
    }else if( pos4710 & pos47 & pos7 ){
      addAllele <- 'KIR2DL1*007'
      cat('\n\nAdjusting KIR2DL1 call to include',addAllele,'based on KFF')
    }else if( pos4710 | pos47 | pos7 ){
      addAllele <- 'KIR2DL1*0040101'
      cat('\n\nAdjusting KIR2DL1 call to include',addAllele,'based on KFF')
    }else{
      addAllele <- integer(0)
    }
    
    calledAlleleVect <- c(calledAlleleVect, addAllele)
    addAllele <- integer(0)
    
    if( pos003 ){
      probeSeq <- probeDF['*KIR2DL1*003','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- 'KIR2DL1*0030201'
        cat('\n\nAdjusting KIR2DL1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    nonAmbAlleleCall.list[['KIR2DL1']] <- paste0(calledAlleleVect,collapse='+')
    
  }
  
  if( 'KIR2DL2' %in% names(nonAmbAlleleCall.list) ){
    posE5292T <- as.integer( currentSample$kffHits[['*KIR2DL2*E5292T']] ) > 10
    pos00103 <- as.integer( currentSample$kffHits[['*KIR2DL2*00103']] ) > 10
    pos004 <- as.integer( currentSample$kffHits[['*KIR2DL2*004']] ) > 10
    
    calledAlleleVect <- unique( unlist( strsplit( nonAmbAlleleCall.list[['KIR2DL2']], '+', fixed=T) ) )
    
    if( posE5292T ){
      probeSeq <- probeDF['*KIR2DL2*E5292T','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DL2 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( pos00103 ){
      probeSeq <- probeDF['*KIR2DL2*00103','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DL2 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }else{
      probeSeq <- probeDF['*KIR2DL2*00103','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      
      if( any( calledAlleleVect %in% hitAllele.vect ) ){
        removeAllele <- "KIR2DL2*00103"
        cat('\n\nAdjusting KIR2DL2 call to exclude',removeAllele,'based on KFF')
        calledAlleleVect <- setdiff( calledAlleleVect, removeAllele )
        
        if( length(calledAlleleVect) == 0 ){
          calledAlleleVect <- 'KIR2DL2*0010101'
        }
      }
    }
    
    if( pos004 ){
      probeSeq <- probeDF['*KIR2DL2*004','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DL2 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    nonAmbAlleleCall.list[['KIR2DL2']] <- paste0(calledAlleleVect,collapse='+')
    
  }
  
  if( 'KIR2DL4' %in% names(nonAmbAlleleCall.list) ){
    pos10A <- as.integer( currentSample$kffHits[['*KIR2DL4*10A64']] ) > 10
    pos9A <- as.integer( currentSample$kffHits[['*KIR2DL4*9A64']] ) > 10
    
    calledAlleleVect <- unique( unlist( strsplit( nonAmbAlleleCall.list[['KIR2DL4']], '+', fixed=T) ) )
    
    if( pos10A & pos9A ){
      addAllele.10A <- 'KIR2DL4*0010201'
      addAllele.9A <- 'KIR2DL4*0080101'
      removeAllele <- ''
    }else if( pos10A & !pos9A ){
      addAllele.10A <- 'KIR2DL4*0010201'
      addAllele.9A <- ''
      removeAllele <- ''#calledAlleleVect[ calledAlleleVect %in% KIR2DL4.9A.vect ]
      if(length(removeAllele) == 0) removeAllele <- ''
    }else if( !pos10A & pos9A ){
      addAllele.10A <- ''
      addAllele.9A <- 'KIR2DL4*0080101'
      removeAllele <- ''#calledAlleleVect[ calledAlleleVect %in% KIR2DL4.10A.vect ]
      if(length(removeAllele) == 0) removeAllele <- ''
    }else{
      addAllele.10A <- ''
      addAllele.9A <- ''
      removeAllele <- ''
    }
    
    probeSeq <- probeDF['*KIR2DL4*10A64','Sequence']
    KIR2DL4.10A.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
    
    probeSeq <- probeDF['*KIR2DL4*9A64','Sequence']
    KIR2DL4.9A.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
    
    if( addAllele.10A != '' & !any(calledAlleleVect %in% KIR2DL4.10A.vect) ){
      cat('\n\nAdjusting KIR2DL4 call to include',addAllele.10A,'based on KFF')
      calledAlleleVect <- c(calledAlleleVect, addAllele.10A)
    }
    
    if( addAllele.9A != '' & !any(calledAlleleVect %in% KIR2DL4.9A.vect) ){
      cat('\n\nAdjusting KIR2DL4 call to include',addAllele.9A,'based on KFF')
      calledAlleleVect <- c(calledAlleleVect, addAllele.9A)
    }
    
    if( nchar(removeAllele[1]) != 0 ){
      cat('\n\nAdjusting KIR2DL4 call to remove',paste(removeAllele,collapse=' '),'based on KFF')
      calledAlleleVect <- calledAlleleVect[ !calledAlleleVect %in% removeAllele ]
    }
    
    if( length(calledAlleleVect) == 0 ){
      stop('HELP, removed all KIR2DL4 reference alleles when performing KFF modifications.')
    }
    
    nonAmbAlleleCall.list[['KIR2DL4']] <- paste0(calledAlleleVect,collapse='+')
  }
  
  if( 'KIR2DS1' %in% names(nonAmbAlleleCall.list) ){
    pos001 <- as.integer( currentSample$kffHits[['*KIR2DS1*001']] ) > 10
    pos00202.00302 <- as.integer( currentSample$kffHits[['*KIR2DS1*00202^00302']] ) > 10
    pos005 <- as.integer(currentSample$kffHits[['*KIR2DS1*005']]) > 10
    pos003 <- as.integer(currentSample$kffHits[['*KIR2DS1*00301']]) > 10
    pos006 <- as.integer(currentSample$kffHits[['*KIR2DS1*006']]) > 10
    pos004.008 <- as.integer(currentSample$kffHits[['*KIR2DS1*004^008']]) > 10
    pos009 <- as.integer(currentSample$kffHits[['*KIR2DS1*009']]) > 10
    #pos010 <- as.integer(currentSample$kffHits[['*KIR2DS1*010']]) > 10
    #pos011 <- as.integer(currentSample$kffHits[['*KIR2DS1*011']]) > 10
    
    calledAlleleVect <- unique( unlist( strsplit( nonAmbAlleleCall.list[['KIR2DS1']], '+', fixed=T) ) )
    
    if( pos001 ){
      probeSeq <- probeDF['*KIR2DS1*001','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    
    if( pos00202.00302 ){
      probeSeq <- probeDF['*KIR2DS1*00202^00302','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( pos005 ){
      probeSeq <- probeDF['*KIR2DS1*005','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( pos003 ){
      probeSeq <- probeDF['*KIR2DS1*00301','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( pos006 ){
      probeSeq <- probeDF['*KIR2DS1*006','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( pos004.008 ){
      probeSeq <- probeDF['*KIR2DS1*004^008','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( pos009 ){
      probeSeq <- probeDF['*KIR2DS1*009','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR2DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    # if( pos010 ){
    #   probeSeq <- probeDF['*KIR2DS1*010','Sequence']
    #   hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
    #   if( !any( calledAlleleVect %in% hitAllele.vect ) ){
    #     addAllele.010 <- hitAllele.vect[1]
    #   }else{
    #     addAllele.010 <- integer(0)
    #   }
    # }else{
    #   addAllele.010 <- integer(0)
    # }
    
    nonAmbAlleleCall.list[['KIR2DS1']] <- paste0(calledAlleleVect,collapse='+')
    
    
  }
  
  if( 'KIR3DP1' %in% names(nonAmbAlleleCall.list) ){
    
    pos3DP1del <- as.integer( currentSample$kffHits[['*KIR3DP1*E2del']] ) > 10
    neg3DP1del <- as.integer( currentSample$kffHits[['*KIR3DP1*E2nodel']] ) > 10
    
    calledAlleleVect <- unique( unlist( strsplit( nonAmbAlleleCall.list[['KIR3DP1']], '+', fixed=T) ) )
    
    if( pos3DP1del ){
      probeSeq <- probeDF['*KIR3DP1*E2del','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- 'KIR3DP1*0030101'
        cat('\n\nAdjusting KIR3DP1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    if( neg3DP1del ){
      probeSeq <- probeDF['*KIR3DP1*E2nodel','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- 'KIR3DP1*0090101'
        cat('\n\nAdjusting KIR3DP1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    nonAmbAlleleCall.list[['KIR3DP1']] <- paste0(calledAlleleVect,collapse='+')
  }
  
  if( 'KIR3DS1' %in% names(nonAmbAlleleCall.list) ){
    pos049 <- as.integer( currentSample$kffHits[['*KIR3DS1*049']] ) > 10
    
    calledAlleleVect <- unique( unlist( strsplit( nonAmbAlleleCall.list[['KIR3DS1']], '+', fixed=T) ) )
    
    if( !any(grepl( 'KIR3DS1*01301', calledAlleleVect, fixed=T)) ){
      cat('\n\nAdjusting KIR3DS1 call to include KIR3DS1*0130101')
      calledAlleleVect <- c(calledAlleleVect, 'KIR3DS1*0130101')
    }
    
    if( pos049 ){
      probeSeq <- probeDF['*KIR3DS1*049','Sequence']
      hitAllele.vect <- names(alleleSeq.list)[ grepl(probeSeq,alleleSeq.list,fixed = T) ]
      if( !any( calledAlleleVect %in% hitAllele.vect ) ){
        addAllele <- hitAllele.vect[1]
        cat('\n\nAdjusting KIR3DS1 call to include',addAllele,'based on KFF')
        calledAlleleVect <- c(calledAlleleVect,addAllele)
        addAllele <- integer(0)
      }
    }
    
    nonAmbAlleleCall.list[['KIR3DS1']] <- paste0(calledAlleleVect,collapse='+')
  }
  
  currentSample[['setCallList']] <- setCall.list
  currentSample[['setupAlleleList']] <- nonAmbAlleleCall.list
  return(currentSample)
}


# ----- REF INFO WRITING FUNCTIONS -----
alleleSetup.write_sample_ref_info <- function( currentSample, alleleSetupDirectory, setup.knownSnpDFList, referenceAlleleDF, addFullyDefined=F){
  
  refInfoPath <- file.path(alleleSetupDirectory,paste0(currentSample$name,'.refInfo.txt'))
  currentSample[['refInfoPath']] <- alleleSetup.initialize_info_file( refInfoPath )
  
  formattedAlleleCall.vect <- as.vector( unlist( sapply( currentSample$setupAlleleList, function(x) unlist( tstrsplit(x, '+', fixed=T )) ) ) )
  
  formattedAlleleCall.vect <- formattedAlleleCall.vect[ !is.na( formattedAlleleCall.vect ) ]
  formattedAlleleCall.vect <- unique( formattedAlleleCall.vect )
  
  lociVect <- names(currentSample$setupAlleleList)
  
  if( addFullyDefined ){
    for(locus in lociVect){
      locus.alleleVect <- grep( locus, formattedAlleleCall.vect, value=T, fixed=T)
      locus.undefinedRegion.boolVect <- apply( setup.knownSnpDFList[[locus]][locus.alleleVect,], 1, function(x) any( !is_nuc(x) ) )
      
      if( all(locus.undefinedRegion.boolVect) ){
        formattedAlleleCall.vect <- c( formattedAlleleCall.vect, referenceAlleleDF[locus,1] )
      }
    }
  }
  
  alleleSetup.write_info( currentSample$name, formattedAlleleCall.vect, currentSample$refInfoPath)
  return(currentSample)
}

alleleSetup.write_info <- function(synSeq.ID, allele.vect, infoPath){
  alleleStr <- paste0(allele.vect,collapse='_')
  cat(paste0(synSeq.ID,'\t',alleleStr,'\n'),file=infoPath,append=T)
}

alleleSetup.initialize_info_file <- function(refInfoPath){
  textStr <- ''
  cat(textStr, file=refInfoPath)
  return(normalizePath(refInfoPath,mustWork = T))
}

# ----- Final geno determination functions -----
samFormat.processCigarStr <- function( cigarStr, readSeq ){
  # INDEL handling
  # deletions - add . character to read for deletion positions marked in CIGAR string
  # insertions - remove insertion characters from read, send read through the indexed vector formatting step, 
  #              then replace the insertion position with the insertion string ( will need to add the current
  #              read del index to correctly place )
  # reads with both - should work following the above steps
  
  ## Split the cigar string up into a vector
  cigarVect <- strsplit(cigarStr,'')[[1]]
  
  ## Determine which position of the vector denote cigar operators
  operationPosVect <- which(cigarVect %in% cigarOperationVect)
  
  ## Initialize a list for storing the different operations and their positions
  operationList <- list()
  prevPos <- 1
  countInt <- 1
  
  ## Convert the cigar vect into the operation list
  for(operationPos in operationPosVect){
    operationList[[countInt]] <- cigarVect[prevPos:operationPos]
    prevPos <- operationPos+1
    countInt <- countInt+1
  }
  
  ## Initialize a start position for CIGAR operations
  prevPos <- 1
  outStr <- ''
  out.insIndex <- list()
  
  #operationVect <- operationList[[6]]
  for(operationVect in operationList){
    
    ## Pull out the end position of the current operation
    operationEndPos <- as.integer(paste0(operationVect[1:(length(operationVect)-1)], collapse='')) + prevPos - 1
    
    ## Pull out the operation type
    operationTypeChar <- operationVect[length(operationVect)]
    
    if(operationTypeChar == 'M'){
      prevPos <- operationEndPos + 1
    }else if(operationTypeChar == 'I'){
      
      preString <- substr(readSeq, 1, prevPos-1)
      insStr <- substr(readSeq, prevPos-1, operationEndPos)
      postString <- substr(readSeq, operationEndPos+1, nchar(readSeq) )
      
      out.insIndex[[as.character( prevPos-1 )]] <- insStr
      readSeq <- paste0( preString, postString )
      
      ## We make no changes to the position variables since the INS seq was removed
      #operationOffset <- operationOffset + nchar(insStr) - 1
      
      #prevPos <- operationEndPos + 1
    }else if(operationTypeChar == 'D'){
      
      seqLength <- nchar(readSeq)
      preString <- substr(readSeq, 1, prevPos-1)
      postString <- substr(readSeq, prevPos, seqLength)
      readSeq <- paste0(c( preString, rep('.', length( prevPos:operationEndPos )), postString), collapse='')
      
      prevPos <- operationEndPos + 1
      
    }else{
      stop('samFormat.adjustReadSeq something weird happened for ',readSeq)
    }
    
  }
  
  return(list('readSeq'=readSeq,'insIndex'=out.insIndex))
}

pingAllele.generate_snp_df <- function( currentSample, uniqueSamDT, output.dir, setup.knownSnpDFList, workflow, hetRatio, minDP){
  
  cat('\n\nSetting read start positions.')
  uniqueSamDT$startPos <- as.integer( apply( uniqueSamDT, 1, function(x) names( x$readTable )[1] ) )
  
  cat('\n\nSetting read end positions.')
  uniqueSamDT$endPos <- as.integer( apply( uniqueSamDT, 1, function(x) names(x$readTable)[length(names(x$readTable))] ) )
  
  cat('\nSorting read table by locus and alignment coordinates.')
  uniqueSamDT <- uniqueSamDT[order(locus, startPos)]
  
  currentSample[['snpDFPathList']][[workflow]] <- list('DP'=list(),'SNP'=list())
  
  cat('\n\nGenerating SNP tables for',paste0(currentSample$name,':'))
  for(currentLocus in unique(uniqueSamDT$locus)){
    cat('',currentLocus)
    
    locusSnpDF <- setup.knownSnpDFList[[currentLocus]]
    locusSamDT <- uniqueSamDT[locus == currentLocus]
    currentDepthDF <- as.data.frame( matrix(data=0,nrow=6,ncol=ncol(locusSnpDF)) )
    
    locusSamDT <- locusSamDT[locusSamDT$endPos < ncol( locusSnpDF ),]
    
    namedReadVect <- unlist( locusSamDT$readTable )
    
    # cat('\n\t\tProcessing:')
    
    # First process insertions
    insIndex <- which( nchar( namedReadVect ) != 1 ,useNames = F)
    if( length(insIndex) > 0){
      #  cat('','INS')
      insTab <- table( as.integer(names(insIndex)) )
      currentDepthDF[ 6, as.integer(names(insTab)) ] <- as.vector( insTab )
      namedReadVect.noIns <- namedReadVect[ -as.vector(insIndex) ]
    }else{
      namedReadVect.noIns <- namedReadVect
    }
    remove(namedReadVect)
    
    # Then process all other nucleotides and deletion positions
    uniqueNucVect <- intersect( unique(namedReadVect.noIns), names(nucListConv) )
    for( nuc in uniqueNucVect ){
      #cat('',nuc)
      nucIndex <- which( namedReadVect.noIns == nuc, useNames=F )
      nucTab <- table( as.integer(names(nucIndex)) )
      currentDepthDF[ nucListConv[[nuc]], as.integer(names(nucTab)) ] <- as.vector( nucTab )
    }
    remove(nucIndex)
    remove(nucTab)
    remove(namedReadVect.noIns)
    
    rownames(currentDepthDF) <- c(names(nucListConv),'INS')
    colnames(currentDepthDF) <- colnames(locusSnpDF)
    
    write.csv(currentDepthDF, file.path(output.dir,paste0(workflow,'_',currentLocus,'_',currentSample$name,'_DP.csv')))
    currentSample[['snpDFPathList']][[workflow]][['DP']][[currentLocus]] <- file.path(output.dir,paste0(workflow,'_',currentLocus,'_',currentSample$name,'_DP.csv'))
    #cat('\n\t\tCompleted depth file generation')
    
    currentSnpDF <- as.data.frame( matrix('',nrow=3,ncol=ncol(currentDepthDF)),stringsAsFactors=F)
    passedMinDP.index <- apply(currentDepthDF,2,sum)>=minDP
    
    currentSnpDF <- currentSnpDF[,passedMinDP.index]
    currentDepthDF <- currentDepthDF[,passedMinDP.index]
    
    currentSnpDF[,] <- apply( currentDepthDF, 2, function(x){
      snpIndex <- which( ( x / max(x) ) > hetRatio )
      snpVect <- names(nucListConvWIns)[snpIndex]
      addVect <- rep('N',max((3-length(snpIndex)),0))
      returnVect <- c(snpVect,addVect)
      return(returnVect[1:3])
    })
    
    colnames(currentSnpDF) <- colnames(currentDepthDF)
    rownames(currentSnpDF) <- c('SNP_1','SNP_2','SNP_3')
    
    write.csv(currentSnpDF, file.path(output.dir,paste0(workflow,'_',currentLocus,'_',currentSample$name,'_SNP.csv')))
    #cat('\n\t\tCompleted SNP file generation')
    currentSample[['snpDFPathList']][[workflow]][['SNP']][[currentLocus]] <- file.path(output.dir,paste0(workflow,'_',currentLocus,'_',currentSample$name,'_SNP.csv'))
  }
  return(currentSample)
}
pingAllele.call_final_alleles <- function( currentSample, currentLocus, refSnpDF ){
  #refSnpDF <- knownSnpDFList[[currentLocus]]$snpDF
  currentADSnpDF <- refSnpDF
  
  exonNameVect <- grep('E',kirLocusFeatureNameList[[currentLocus]],value=T)
  exonNameVect <- grep('PE',exonNameVect,value=T,invert = T)
  
  currentSnpDF <- read.csv(file.path(currentSample[['snpDFPathList']][['final']][['SNP']][[currentLocus]]),
                           check.names=F,stringsAsFactors = F,row.names=1,header = T,colClasses = c("character"))
  
  if( currentLocus != 'KIR2DS35' ){
    sampleCols <- tstrsplit(colnames(currentSnpDF),'_',fixed=T)[[1]] 
    sampleCols <- colnames(currentSnpDF)[ sampleCols %in% exonNameVect ]
    currentSnpDF <- currentSnpDF[,sampleCols]
    
    totalCols <- tstrsplit(colnames(currentADSnpDF),'_',fixed=T)[[1]]
    totalCols <- colnames(currentADSnpDF)[ totalCols %in% exonNameVect ]
    
    exonPosPerc <- paste0( round( length(sampleCols) / length(totalCols), 3)*100,'%')
    
    cat('\n\t\tCalling on',exonPosPerc,'of exon positions')
  }
  
  #currentSample[['callList']][[currentLocus]] <- alleleSetup.call_allele( currentSample, currentLocus, currentSnpDF, currentADSnpDF, homScoreBuffer=1, kirRes=5, ambScore=F )
  
  currentSample[['callList']][[currentLocus]]  <- tryCatch({
    alleleSetup.call_allele( currentSample, currentLocus, currentSnpDF, currentADSnpDF, homScoreBuffer=1, kirRes=5, ambScore=F )
  },
  error=function(cond) {
    message(paste("Allele calling error for:", currentLocus))
    message("Error message:")
    message(cond)
    # Choose a return value in case of error
    return(list('alleleVect'='failed','scoreInt'=0))
  })
  
  #unresLocus.vect <- which( sapply( currentSample[['callList']], function(x) x$scoreInt > 0 ) )
  
  if( currentSample[['callList']][[currentLocus]]$scoreInt > 0 ){
    currentSample <- pingAllele.addUnresSNPs( currentSample, currentADSnpDF, currentSnpDF, currentLocus )
  }
  
  return(currentSample)
}

pingAllele.save_call <- function( currentSample, acDFPath ){
  
  # Failure condition
  #if( currentSample[[ paste0(workflow,'RefDirectory') ]] == 'failed' ){
  #  currentSample[[ paste0(workflow,'AlleleCallDF') ]] <- 'failed'
  #  return(currentSample)
  #}
  
  cat('\n\nSaving allele calls to',acDFPath)
  acDF <- read.table(acDFPath, check.names = F, stringsAsFactors = F, sep=',',colClasses = 'character',comment.char='')
  acDF[currentSample$name,] <- NA
  
  for( currentLocus in names( currentSample$callList ) ){
    alleleCall <- paste(currentSample$callList[[currentLocus]]$alleleVect, collapse=' ')
    callScore <- currentSample$callList[[currentLocus]]$scoreInt
    
    #if( callScore > 0 ){
    #  
    #  alleleCall <- 'unresolved_mismatch'
    #}
    
    acDF[currentSample$name,currentLocus] <- alleleCall
  }
  
  write.table( acDF, file=acDFPath, sep=',', quote=F, row.names=T, col.names=T )
  
  return(currentSample)
}

alleleSetup.prep_results_directory <- function(currentSample, alignmentFileDirectory){
  currentSample[['snpDFPathList']] <- list()
  # Create a directory for the current iteration reference files
  currentSampleResultsDirectory <- file.path(alignmentFileDirectory,currentSample$name,'iterAlign')
  if(!file.exists(currentSampleResultsDirectory)){
    dir.create(currentSampleResultsDirectory, recursive=T)
  }
  currentSample[['iterRefDirectory']] <- currentSampleResultsDirectory
  
  return(currentSample)
}

general.kir_upres <- function( alleleType, setup.knownSnpDFList ){
  locus <- strsplit(alleleType,'*',fixed=T)[[1]][1]
  
  if(locus == 'KIR2DL5A' | locus == 'KIR2DL5B'){
    locus <- 'KIR2DL5'
  }
  
  return( grep( alleleType, rownames(  setup.knownSnpDFList[[locus]]), fixed=T, value=T) )
}

alleleSetup.most_defined_typing <- function( currentAlleleCall.list ){
  
  cat('\n\nReducing genotype ambiguity by selecting most defined sequences')
  nonAmbAlleleCall.list <- sapply( names(currentAlleleCall.list), function(x){
    alleleCallVect <- currentAlleleCall.list[[x]]
    
    if(length(alleleCallVect) == 1){
      return(alleleCallVect)
    }
    
    alleleCallList <- strsplit(alleleCallVect, '+', fixed=T)
    
    nonAmbSumVect <- sapply( alleleCallList, function(y){
      sum( sapply( y, function(z){
        sum( setup.knownSnpDFList[[x]][z,] != '*' )
      }) )
    })
    
    selectedCall <- which( nonAmbSumVect == max( nonAmbSumVect ) )[1]
    return( alleleCallVect[selectedCall] )
  })
  
  return(nonAmbAlleleCall.list)
}

alleleSetup.finalize_set_calls <- function( currentSample ){
  setCall.list <- currentSample$setCallList
  setupCall.list <- currentSample$setupAlleleList
  
  outCall.list <- list()
  
  if(length(setCall.list) > 0){
    
    setCall.list <- alleleSetup.most_defined_typing(setCall.list)
    setLocus.vect <- intersect( names(setCall.list), names(setupCall.list) )
    
    for( currentLocus in setLocus.vect ){
      outCall.list[[currentLocus]] <- setCall.list[[currentLocus]]
    }
  }
  
  unsetLocus.vect <- setdiff( names(setupCall.list), names(setCall.list) )
  
  for( currentLocus in unsetLocus.vect ){
   outCall.list[[currentLocus]] <- setupCall.list[[currentLocus]]
  }
  
  currentSample[['setupAlleleList']] <- outCall.list
  return(currentSample)
}

pingAllele.addUnresSNPs <- function( currentSample, currentADSnpDF, currentSnpDF, currentLocus ){
  rawAllele.vect <- currentSample[['callList']][[currentLocus]]$alleleVect
  
  out.modAllele.vect <- c()
  for( rawGenoCall in rawAllele.vect ){
    
    alleleVect <- strsplit(rawGenoCall, '+', fixed=T)[[1]]
    alleleVect <- unique( unlist( sapply(alleleVect, function(x) grep( x, rownames(currentADSnpDF), fixed=T, value=T)[1] ) ) )
    
    colID.vect <- colnames(currentSnpDF)
    unresolvedAllele.list <- list()
    for( allele in alleleVect ){
      unresolvedAllele.list[[allele]] <- integer(0)
    }
    
    for( colID in colID.vect ){
      snpVect <- currentSnpDF[,colID]
      snpVect <- snpVect[ is_nuc( snpVect ) ]
      
      refSnpVect <- unique( currentADSnpDF[alleleVect, colID] )
      
      if( !all(snpVect %in% refSnpVect) ){
        #cat('',colID)
        
        diffSnp <- setdiff( snpVect, refSnpVect )
        addSnp <- paste0(colID,'.',diffSnp)
        unresolvedAllele.list <- lapply( unresolvedAllele.list, function(x) c(x,addSnp))
      }
      
      if( !all(refSnpVect %in% snpVect) ){
        #cat('\n',colID)
        mismatchAllele.vect <- alleleVect[ !(currentADSnpDF[alleleVect, colID, drop=T] %in% snpVect) ]
        for( mismatchAllele in mismatchAllele.vect ){
          diffSnp <- currentADSnpDF[mismatchAllele,colID,drop=T]
          addSnp <- paste0(colID,'.',diffSnp)
          unresolvedAllele.list[[mismatchAllele]] <- c(unresolvedAllele.list[[mismatchAllele]], addSnp)
        }
      }
    }
    
    unresolvedAllele.list <- lapply( unresolvedAllele.list, paste0, collapse='^' )
    
    unresAllele.vect <- c()
    for( alleleID in names(unresolvedAllele.list) ){
      if( nchar( unresolvedAllele.list[[alleleID]] ) == 0 ){
        unresAllele.vect <- c(unresAllele.vect, alleleID)
      }else{
        unresAllele.vect <- c(unresAllele.vect, paste0(alleleID,'$',unresolvedAllele.list[[alleleID]]) )
      }
    }
    
    unresGenoCall <- paste0( unresAllele.vect, collapse='+' )
    out.modAllele.vect <- c(out.modAllele.vect, unresGenoCall)
  }
  
  currentSample[['callList']][[currentLocus]]$alleleVect <- out.modAllele.vect
  return(currentSample)
}

ping_allele <- function(sampleList){
  
  for( currentSample in sampleList ){
    
    currentSample <- alleleSetup.gc_matched_ref_alignment( currentSample, alleleSetupDirectory, as.list, threads)
    
    currentSample[['setCallList']] <- list()
    uniqueSamDT <- alleleSetup.process_samDT( currentSample$ASSamPath, delIndex.list, processSharedReads = setup.readBoost, readBoost.thresh )
    file.remove(currentSample$ASSamPath)
    currentSample <- alleleSetup.prep_results_directory( currentSample, alignmentFileDirectory )
    currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP, includeAmb = T )
    currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory, setup.knownSnpDFList, alleleSetupRef.df, addFullyDefined = F)
    
    synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
    currentSample <- ping_iter.run_alignments(currentSample, threads, all.align=T, synSeq.key)
    uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = T, 2 )
    currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP, includeAmb=T, ambScore=F, allPosScore = F, homScoreBuffer = 1, addDiversity = F)
    currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, final.hetRatio, final.minDP, includeAmb=T, ambScore=T, allPosScore = F, onlyExonScore = T, homScoreBuffer = 1, addDiversity=T,combineTypings=T, skipSnpGen=T)
    currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory, setup.knownSnpDFList, alleleSetupRef.df, addFullyDefined = T)
    
    synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
    currentSample <- ping_iter.run_alignments(currentSample, threads, all.align=F, synSeq.key)
    uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = F, 2 )
    currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.hetRatio, setup.minDP, includeAmb=T, ambScore=F, allPosScore = F, homScoreBuffer = 1, addDiversity = F, skipSet = F)
    currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, final.hetRatio, final.minDP, includeAmb=T, ambScore=T, allPosScore = F, onlyExonScore = T, homScoreBuffer = 1, addDiversity=F,combineTypings=T,skipSnpGen=T, skipSet=F)
    currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory, setup.knownSnpDFList, alleleSetupRef.df, addFullyDefined = T)
    
    synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
    currentSample <- ping_iter.run_alignments(currentSample, threads, all.align=F, synSeq.key)
    uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list, processSharedReads = F, 2 )
    currentSample <- pingAllele.generate_snp_df( currentSample,uniqueSamDT,currentSample[['iterRefDirectory']],setup.knownSnpDFList,'final', final.hetRatio, final.minDP )
    
    cat('\n\n\n----- Final allele calling -----')
    for( currentLocus in names( currentSample[['snpDFPathList']][['final']][['SNP']] )){
      cat('\n\t',currentLocus)
      currentSample <- pingAllele.call_final_alleles(currentSample, currentLocus, knownSnpDFList[[currentLocus]]$snpDF)
    }
    
    currentSample <- pingAllele.save_call( currentSample, alleleDFPathList$iter$alleleCallPath )
  }
  
  return(sampleList)
}




