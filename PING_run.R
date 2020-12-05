
setwd('/home/LAB_PROJECTS/PING2_PAPER/PING2') #Set this to your own PING2 working directory

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


# Initialization variables ------------------------------------------------
rawFastqDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/sequence_data/' # can be set to raw sequence or extractedFastq directory
fastqPattern <- '_KIR_' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
threads <- 90
resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_validation_gc_align_method/' # Set the master results directory (all pipeline output will be recorded here)
shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)
minDP <- 20
setup.minDP <- 6


source('Resources/general_functions.R') # do not change
source('Resources/extractor_functions.R') # do not change
source('Resources/ping_copy.R') # do not change
source('Resources/ping_allele.R') # do not change
source('Resources/ping_gc_align.R') # do not change
source('Resources/alleleCombine_functions.R') # do not change
setDTthreads(threads)

# Preparation -------------------------------------------------------------
# Build up a list of sample objects
sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, resultsDirectory, shortNameDelim) # no need to change


# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq') # no need to change
# Run PING2 extractor
#ext.startTime <- Sys.time()
sampleList <- extractor.run(sampleList,threads,extractedFastqDirectory,forceRun=F) # set forceRun=T if you want to force alignments
#ext.endTime <- Sys.time()
#cat('\nExtractor code timing:',ext.endTime - ext.startTime)

#copy.startTime <- Sys.time()
# PING2 gene content and copy number --------------------------------------
cat('\n\n----- Moving to PING gene content and copy determination -----')
sampleList <- ping_copy.graph(sampleList=sampleList,threads=threads,resultsDirectory=resultsDirectory,forceRun=T,onlyKFF=F) # set forceRun=T if you want to force alignments
#copy.endTime <- Sys.time()
#cat('\nCopy code timing:',copy.endTime - copy.startTime)

' 6.12 hours for 50 samples at 40 threads, 150bp, 50dp
'

sampleList <- ping_copy.manual_threshold(sampleList=sampleList,resultsDirectory=resultsDirectory,use.threshFile = T) # this function sets copy thresholds

## Seeing the following message when generating copy plots is normal
'A line object has been specified, but lines is not in the mode
Adding lines to the mode...'

## use this section for notes / recording suspicious samples
' 
KIR3DP1 [example]
  IND00001 [example]
'

# ----- PING allele setup -----

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
alleleSeq.list <- alleleSeq.list[ names(alleleSeq.list) %in% unlist(alleleSetupRef.df) ]

delIndex.list <- list()
for( locus in names(setup.knownSnpDFList) ){
  snpDF <- setup.knownSnpDFList[[locus]]
  
  if(locus == 'KIR2DS1'){
  delIndex.list[[locus]] <- sapply(rownames(snpDF), function(x) {
    which(snpDF[x,] == '.')
  })
  }else{
    odelIndex.list[[locus]] <- apply(snpDF, 1, function(x) {
      which(x == '.')
    })
  }
}

alleleSetupDirectory <- file.path(resultsDirectory,'allele_setup_files')
if(!file.exists(alleleSetupDirectory)){
  dir.create(alleleSetupDirectory)
}


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
  
  currentSample[['ASIndexPath']] <- file.path(alleleSetupDirectory,paste0(currentSample$name,'.alleleSetupRef'))
  
  ## Creqte a bowtie2 index for the kir_reference.fasta file <- only needed when building a new reference index
  createIndex <- system2(bowtie2Build, c(currentSample[['ASFastaPath']], 
                                         currentSample[['ASIndexPath']],'--quiet',paste('--threads', threads)))
  check.system2_output(createIndex, 'bowtie2 index building failed')
  
  return(currentSample)
}
alleleSetup.bt2_align <- function( currentSample, alleleSetupDirectory, bowtie2, threads ){
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

alleleSetup.process_samDT <- function( samPath, delIndex.list ){
  
  #currentSample <- samtools.bam_to_sam(samtools, currentSample, '/home/LAB_PROJECTS/PING2_PAPER/3_validation_gc_align_method/gc_bam_files/', threads)
  
  # Read in SAM file, then delete SAM to save space
  nHead <- general.count_sam_header( samPath )
  samTable <- general.read_sam( samPath , rows_to_skip = nHead )
  
  #nHead <- general.count_sam_header( currentSample$samPath )
  #samTable <- general.read_sam( currentSample$samPath, rows_to_skip = nHead )
  
  
  # Separate KIR2DL1 & KIR2DS1 reads
  # if( 'KIR2DL1' %in% samTable$locus & 'KIR2DS1' %in% samTable$locus & sepL1S1 ){
  #   L1.readVect <- unique( samTable[locus == 'KIR2DL1']$read_name )
  #   S1.readVect <- unique( samTable[locus == 'KIR2DS1']$read_name )
  #   L1S1.readVect <- intersect( L1.readVect, S1.readVect )
  #   
  #   subSamTable <- samTable[read_name %in% L1S1.readVect][,sum(alignment_score),by=c('read_name','reference_name','locus')]
  #   
  #   maxScore.subSamTable <- subSamTable[subSamTable[,V1 >= (max(V1)),by=c('read_name')]$V1,]
  #   maxScore.isUniqueDT <- maxScore.subSamTable[,length(unique(locus)) == 1,by='read_name']
  #   maxScore.uniqueReadVect <- maxScore.isUniqueDT$read_name[maxScore.isUniqueDT$V1]
  #   
  #   L1.selectedReadVect <- unique( maxScore.subSamTable[ read_name %in% maxScore.uniqueReadVect & locus == 'KIR2DL1']$read_name )
  #   S1.selectedReadVect <- unique( maxScore.subSamTable[ read_name %in% maxScore.uniqueReadVect & locus == 'KIR2DS1']$read_name )
  #   
  #   L1.samTable <- samTable[ read_name %in% L1.selectedReadVect & locus == 'KIR2DL1' ]
  #   S1.samTable <- samTable[ read_name %in% S1.selectedReadVect & locus == 'KIR2DS1' ]
  #   
  #   L1S1.samTable <- merge(L1.samTable, S1.samTable, all=T)
  #   
  #   L1S1.selectedReadVect <- unique( c(L1.selectedReadVect, S1.selectedReadVect) )
  #   
  #   samTable <- merge( samTable[ !(read_name %in% L1S1.selectedReadVect) ], L1S1.samTable, all=T, by=colnames(L1S1.samTable))
  # }
  
  samTable$combLocus <- samTable$locus
  
  # if(combS35){
  #   if( all(c('KIR2DS3','KIR2DS5') %in% unique(samTable$locus) )){
  #     samTable[ locus %in% c('KIR2DS3','KIR2DS5')]$combLocus <- 'KIR2DS35'
  #   }
  # }

  # Pull out reads that align to s single locus
  isUniqueDT <- samTable[, length( unique(combLocus) ) == 1, by=read_name ]
  
  uniqueReadVect <- isUniqueDT$read_name[isUniqueDT$V1]
  sharedReadVect <- isUniqueDT$read_name[!isUniqueDT$V1]
  
  uniqueSamDT <- copy( samTable[read_name %in% uniqueReadVect] )
  uniqueSamDT <- unique(uniqueSamDT, by=c('read_seq'))
  
  rm(samTable)
  
  #uniqueSamDT$readTable <- list()
  
  #cigarMod.indel.index <- grepl('D', uniqueSamDT$cigarStr) & grepl('I', uniqueSamDT$cigarStr)
  
  cat('\nFormatting reads.')
  cigarOperationVect <- c('M','I','D')
  #uniqueSamDT[cigarMod.del.index,'read_seq' := samFormat.adjustReadSeq( cigarStr, read_seq), by=seq_len(nrow(uniqueSamDT))]
  #uniqueSamDT$readLen <- nchar(uniqueSamDT$read_seq)
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

alleleSetup.generate_rawSnpTbl <- function( currentLocus, currentLocusSnpDF, currentADSnpDF, uniqueSamDT, minDP ){
  
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
  rawSnpTbl <- rawSnpTbl[rawSnpTbl>=minDP]
  
  if(length( rawSnpTbl ) <= 1 ){
    rawSnpTbl <- table(unlist(formattedADSnp.list))
    rawSnpTbl <- rawSnpTbl[rawSnpTbl>=1]
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
    
    adMat <- as.matrix( tempADSnpDF[,nonAmbColVect] )
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

alleleSetup.call_setup_alleles <- function( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.minDP ){
  
  # ----- CALL ALLELES -----
  presentLociVect <- names(currentSample$copyNumber)[currentSample$copyNumber > 0]
  currentAlleleCall.list <- list()
  cat('\nProcessing:')
  for(currentLocus in presentLociVect){
    cat('\n\n\t',currentLocus)
    
    currentLocusSnpDF <- setup.knownSnpDFList[[currentLocus]]
    currentADSnpDF <- make_unique_pos_frame( currentLocusSnpDF )
    
    rawSnpTbl <- alleleSetup.generate_rawSnpTbl( currentLocus , currentLocusSnpDF, currentADSnpDF, uniqueSamDT, setup.minDP )
    currentSnpDF <- alleleSetup.generate_snpDF( rawSnpTbl )
    
    callList <- alleleSetup.call_allele( currentSample, currentLocus, currentSnpDF, currentADSnpDF, homScoreBuffer=4, kirRes=7, ambScore=T )
    
    # if(callList$scoreInt > 100){
    #   cat('\n\t\tRetrying without ambiguity scoring.')
    #   secondCallList <- alleleSetup.call_allele( currentSample, currentLocus, currentSnpDF, currentADSnpDF, homScoreBuffer=2, kirRes=7, ambScore=F )
    # 
    #   firstCallScore <- callList$scoreInt
    #   secondCallScore <- secondCallList$scoreInt
    # 
    #   if( secondCallScore < firstCallScore ){
    #     cat('\n\t\tUsing second call.')
    #     callList <- secondCallList
    #   }else{
    #     cat('\n\t\tUsing first call.')
    #   }
    # }
    
    currentAlleleCall.list[[currentLocus]] <- callList$alleleVect
  }
  
  # ----- REF ALLELE SELECTION -----
  '
  Reference alleles are chosen based on which ambiguity has the most defined positions
  '
  
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
  
  if( 'KIR2DL1' %in% names(nonAmbAlleleCall.list) ){
    pos4710 <- as.integer( currentSample$kffHits[['*KIR2DL1*4710b']] ) > 10
    pos47 <- as.integer( currentSample$kffHits[['*KIR2DL1*4and7']] ) > 10
    pos7 <- as.integer( currentSample$kffHits[['*KIR2DL1*7only']] ) > 10
    
    if( pos4710 & !pos47 & !pos7 ){
      addAllele <- 'KIR2DL1*010'
    }else if( pos4710 & pos47 & !pos7 ){
      addAllele <- 'KIR2DL1*0040101'
    }else if( pos4710 & pos47 & pos7 ){
      addAllele <- 'KIR2DL1*007'
    }else if( pos4710 | pos47 | pos7 ){
      addAllele <- 'KIR2DL1*0040101'
    }else{
      addAllele <- ''
    }
    
    if(nchar(addAllele) > 0){
      if( !grepl( kir.allele_resolution(addAllele,5), nonAmbAlleleCall.list[['KIR2DL1']], fixed=T) ){
        cat('\n\nAdjusting KIR2DL1 call to include',addAllele,'based on KFF')
        nonAmbAlleleCall.list[['KIR2DL1']] <- paste0(nonAmbAlleleCall.list[['KIR2DL1']],'+',addAllele)
      }
    }
    
  }
  
  #if( 'KIR2DL4' %in% names(nonAmbAlleleCall.list) ){
  #  pos10A <- as.integer( currentSample$kffHits[['*KIR2DL4*10A64']] ) > 10
  #  pos9A <- as.integer( currentSample$kffHits[['*KIR2DL4*9A64']] ) > 10
  #}
  
  currentSample[['setupAlleleList']] <- nonAmbAlleleCall.list
  return(currentSample)
}


# ----- REF INFO WRITING FUNCTIONS -----
alleleSetup.write_sample_ref_info <- function( currentSample, alleleSetupDirectory ){
  
  refInfoPath <- file.path(alleleSetupDirectory,paste0(currentSample$name,'.refInfo.txt'))
  currentSample[['refInfoPath']] <- alleleSetup.initialize_info_file( refInfoPath )
  
  formattedAlleleCall.vect <- unlist( sapply( currentSample$setupAlleleList, function(x) unlist( tstrsplit(x, '+', fixed=T )) ) )
  
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


# ----- Final geno determination functions

samFormat.processCigarStr <- function( cigarStr, readSeq ){
  
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

# INDEL handling
# deletions - add . character to read for deletion positions marked in CIGAR string
# insertions - remove insertion characters from read, send read through the indexed vector formatting step, 
#              then replace the insertion position with the insertion string ( will need to add the current
#              read del index to correctly place )
# reads with both - should work following the above steps

pingAllele.generate_snp_df <- function( currentSample, uniqueSamDT, setup.knownSnpDFList, hetRatio ){
  
  cat('\n\nSetting read start positions.')
  uniqueSamDT$startPos <- as.integer( apply( uniqueSamDT, 1, function(x) names( x$readTable )[1] ) )
  
  cat('\nSorting read table by locus and alignment coordinates.')
  uniqueSamDT <- uniqueSamDT[order(locus, startPos)]
  
  currentSample[['snpDFPathList']] <- list()
  
  cat('\n\nGenerating SNP tables for',currentSample$name)
  for(currentLocus in unique(uniqueSamDT$locus)){
    cat('\n\n\t',currentLocus)
    
    locusSnpDF <- setup.knownSnpDFList[[currentLocus]]
    locusSamDT <- uniqueSamDT[locus == currentLocus]
    currentDepthDF <- as.data.frame( matrix(data=0,nrow=6,ncol=ncol(locusSnpDF)) )
    
    namedReadVect <- unlist( locusSamDT$readTable )
    
    # First process insertions
    insIndex <- which( nchar( namedReadVect ) != 1 ,useNames = F)
    if( length(insIndex) > 0){
      insTab <- table( as.integer(names(insIndex)) )
      currentDepthDF[ 6, as.integer(names(insTab)) ] <- as.vector( insTab )
      namedReadVect.noIns <- namedReadVect[ -as.vector(insIndex) ]
    }else{
      namedReadVect.noIns <- namedReadVect
    }
    remove(namedReadVect)
    
    cat('\n\t\tProcessing:')
    # Then process all other nucleotides and deletion positions
    uniqueNucVect <- intersect( unique(namedReadVect.noIns), names(nucListConv) )
    for( nuc in uniqueNucVect ){
      cat('',nuc)
      nucIndex <- which( namedReadVect.noIns == nuc, useNames=F )
      nucTab <- table( as.integer(names(nucIndex)) )
      currentDepthDF[ nucListConv[[nuc]], as.integer(names(nucTab)) ] <- as.vector( nucTab )
    }
    remove(nucIndex)
    remove(nucTab)
    remove(namedReadVect.noIns)
    
    rownames(currentDepthDF) <- c(names(nucListConv),'INS')
    colnames(currentDepthDF) <- colnames(locusSnpDF)
    
    write.csv(currentDepthDF, file.path(currentSample$iterRefDirectory,paste0(currentLocus,'_DP.csv')))
    cat('\n\t\tCompleted depth file generation')
    
    currentSnpDF <- as.data.frame( matrix('',nrow=3,ncol=ncol(currentDepthDF)),stringsAsFactors=F)
    passedMinDP.index <- apply(currentDepthDF,2,sum)>=minDP
    
    currentSnpDF <- currentSnpDF[,passedMinDP.index]
    currentDepthDF <- currentDepthDF[,passedMinDP.index]
    
    currentSnpDF[,] <- apply( currentDepthDF, 2, function(x){
      snpIndex <- which( ( x / max(x) ) > 0.25 )
      snpVect <- names(nucListConvWIns)[snpIndex]
      addVect <- rep(snpVect[1],max((3-length(snpIndex)),0))
      returnVect <- c(snpVect,addVect)
      return(returnVect[1:3])
    })
    
    colnames(currentSnpDF) <- colnames(currentDepthDF)
    rownames(currentSnpDF) <- c('SNP_1','SNP_2','SNP_3')
    
    write.csv(currentSnpDF, file.path(currentSample$iterRefDirectory,paste0(currentLocus,'_SNP.csv')))
    cat('\n\t\tCompleted SNP file generation')
    currentSample[['snpDFPathList']][[currentLocus]] <- file.path(currentSample$iterRefDirectory,paste0(currentLocus,'_SNP.csv'))
  }
  return(currentSample)
}

pingAllele.call_final_alleles <- function( currentSample, currentLocus, refSnpDF ){
  currentADSnpDF <- refSnpDF
  
  
  exonNameVect <- grep('E',kirLocusFeatureNameList[[currentLocus]],value=T)
  exonNameVect <- grep('PE',exonNameVect,value=T,invert = T)
  
  currentSnpDF <- read.csv(file.path(currentSample$iterRefDirectory,paste0(currentLocus,'_SNP.csv')),
                           check.names=F,stringsAsFactors = F,row.names=1,header = T,colClasses = c("character"))
  
  sampleCols <- tstrsplit(colnames(currentSnpDF),'_',fixed=T)[[1]] 
  sampleCols <- colnames(currentSnpDF)[ sampleCols %in% exonNameVect ]
  currentSnpDF <- currentSnpDF[,sampleCols]
  
  totalCols <- tstrsplit(colnames(currentADSnpDF),'_',fixed=T)[[1]]
  totalCols <- colnames(currentADSnpDF)[ totalCols %in% exonNameVect ]
  
  exonPosPerc <- paste0( round( length(sampleCols) / length(totalCols), 3)*100,'%')
  
  cat('\n\t\tCalling on',exonPosPerc,'of exon positions')
  
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
  
  return(currentSample)
}

pingAllele.save_call <- function( currentSample, acDFPath ){
  
  # Failure condition
  #if( currentSample[[ paste0(workflow,'RefDirectory') ]] == 'failed' ){
  #  currentSample[[ paste0(workflow,'AlleleCallDF') ]] <- 'failed'
  #  return(currentSample)
  #}
  
  cat('\n\nSaving allele calls to',acDFPath)
  acDF <- read.table(acDFPath, check.names = F, stringsAsFactors = F, sep=',')
  acDF[currentSample$name,] <- NA
  
  for( currentLocus in names( currentSample$callList ) ){
    alleleCall <- paste(currentSample$callList[[currentLocus]]$alleleVect, collapse=' ')
    callScore <- currentSample$callList[[currentLocus]]$scoreInt
    
    if( callScore > 0 ){
      alleleCall <- 'unresolved_mismatch'
    }
    
    acDF[currentSample$name,currentLocus] <- alleleCall
  }
  
  write.table( acDF, file=acDFPath, sep=',', quote=F, row.names=T, col.names=T )
  
  return(currentSample)
}


# PING2 alignments and allele calling ----------------------------------------------
# Iter align workflow
source('Resources/genotype_alignment_functions.R') # do not change

# Alignment and allele calling workflow
for(currentSample in sampleList[70:length(sampleList)]){
  
  if(currentSample$name %in% c('IND00003','IND00007','IND00074','IND00075')){
    next
  }
  
  currentSample <- alleleSetup.gc_matched_ref_alignment( currentSample, alleleSetupDirectory, alleleSeq.list, threads)
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$ASSamPath, delIndex.list )
  #uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list )
  #file.remove(currentSample$ASSamPath)
  currentSample <- alleleSetup.call_setup_alleles( currentSample, uniqueSamDT, setup.knownSnpDFList, setup.minDP )
  currentSample <- alleleSetup.write_sample_ref_info( currentSample, alleleSetupDirectory )
  
  synSeq.key <- alleleSetup.readAnswerKey( currentSample$refInfoPath )
  
  currentSample <- ping_iter.run_alignments(currentSample, threads)
  
  uniqueSamDT <- alleleSetup.process_samDT( currentSample$iterSamPathList[[1]], delIndex.list )
  
  
  currentSample <- pingAllele.generate_snp_df( currentSample, uniqueSamDT, setup.knownSnpDFList, 0.25 )
  
  
  
  cat('\n\n\n----- Final allele calling -----')
  for( currentLocus in names( currentSample[['snpDFPathList']] )){
    cat('\n\t',currentLocus)
    currentSample <- pingAllele.call_final_alleles(currentSample, currentLocus, knownSnpDFList[[currentLocus]]$snpDF)
  }
  
  currentSample <- pingAllele.save_call( currentSample, alleleDFPathList$iter$alleleCallPath )

}

