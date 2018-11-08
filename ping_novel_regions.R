library(data.table)
library(ggplot2)
library(stringr)
library(methods)
library(plotly)

source('Resources/gc_functions.R')


########## INPUT variables
setwd('/home/wmarin/PING_projects/PING2/')
sampleDirectory <- '/home/wmarin/PING_projects/2_copy_1_sequences/'
fastqPattern <- 'fastq'
threads <- 12
resultsDirectory <- 'copy_number_1_samples'
KIR3DL3MinReadThreshold <- 100
maxReadThreshold <- 30000
probelistFile <- 'probelist_2018_08_02.csv'
depthThreshold <- 8
hetRatio <- 0.25
###########

#allCopyNumberList <- list('IND00109_S26_L001__KIR_'=list('KIR3DP1'=2,'KIR2DL3'=1,'KIR2DP1'=1,'KIR2DS2'=1,'KIR2DL4'=2,'KIR3DL3'=2,'KIR3DL1'=2,'KIR2DL2'=1,'KIR3DL2'=2,'KIR2DS4'=2,'KIR2DL1'=1),
#                          'IND00110_S34_L001__KIR_'=list('KIR3DP1'=2,'KIR2DL3'=1,'KIR2DP1'=1,'KIR2DS2'=1,'KIR2DL4'=2,'KIR3DL3'=2,'KIR3DL1'=2,'KIR2DL2'=1,'KIR3DL2'=2,'KIR2DS4'=2,'KIR2DL1'=1),
#                          'IND00112_S48_L001__KIR_'=list('KIR3DP1'=2,'KIR2DS5'=1,'KIR2DL3'=2,'KIR2DP1'=2,'KIR2DL4'=2,'KIR3DL3'=2,'KIR3DL1'=1,'KIR3DS1'=1,'KIR3DL2'=2,'KIR2DS4'=1,'KIR2DL1'=2, 'KIR2DS1'=1, 'KIR2DL5'=1),
#                          'IND00114_S63_L001__KIR_'=list('KIR3DP1'=2,'KIR2DS5'=1,'KIR2DL3'=1,'KIR2DP1'=1,'KIR2DS2'=1,'KIR2DL4'=2,'KIR3DL3'=2,'KIR3DL1'=1,'KIR3DS1'=1,'KIR2DL2'=1,'KIR3DL2'=2,'KIR2DS4'=1,'KIR2DL1'=1, 'KIR2DS1'=1, 'KIR2DL5'=2),
#                          'IND00115_S71_L001__KIR_'=list('KIR3DP1'=2,'KIR2DL3'=2,'KIR2DP1'=2,'KIR2DL4'=2,'KIR3DL3'=2,'KIR3DL1'=2,'KIR3DL2'=2,'KIR2DS4'=2,'KIR2DL1'=2))


kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

cat('Current working directory: ', getwd(),'\n')

resultsDirectory <- normalizePath(resultsDirectory, mustWork=T)
sampleDirectory <- normalizePath(sampleDirectory, mustWork=T)
assembledReferenceDirectory <- normalizePath(file.path(resultsDirectory,'assembled_references'), mustWork=T)
gcResourceDirectory <- normalizePath('Resources/gc_resources', mustWork = T)

kirReferenceFasta <- normalizePath(file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled.fasta'), mustWork=T)
kirReferenceAlignedFasta <- normalizePath(file.path(gcResourceDirectory, 'filled_kir_reference', 'KIR_gen_onelines_aligned_filled.fasta'), mustWork=T)
kirReferenceIndex <- file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled')

kirAlleleList <- read.kir_allele_list_from_reference_fasta(kirReferenceFasta)
kirAlleleListRes3 <- unique(unlist(lapply(kirAlleleList, kir.allele_resolution, 3)))

## Set file names
kffPresenceDFFile <- file.path(resultsDirectory, 'kffPresenceFrame.csv')
alleleCountDFFile <- file.path(resultsDirectory, 'alleleCountFrame.csv')
copyNumberDFFile <- file.path(resultsDirectory, 'copyNumber.csv')

## Load in found count file
cat(paste0('\nFound alleleCountFrame.csv in ', resultsDirectory, '. Loading these results.'))  
alleleCountDF <- read.csv(alleleCountDFFile, stringsAsFactors = F, check.names = F, row.names = 1)

## Load in found count file
cat(paste0('\nFound kffPresenceFrame.csv in ', resultsDirectory, '. Loading these results.'))
kffPresenceDF <- read.csv(kffPresenceDFFile, stringsAsFactors = F, check.names = F, row.names = 1)

## Load in copy number file
cat(paste0('\nFound copyNumber.csv in ',resultsDirectory,'. Loading these results.'))
copyNumberDF <- read.csv(copyNumberDFFile, stringsAsFactors=F,check.names=F,row.names=1)

## Building a list of sample objects from files in sampleDirectory that match fastqPattern
sampleList <- build.paired_sample_objects(sampleDirectory,fastqPattern,resultsDirectory)

## Read in the aligned fasta and convert it into a list of dataframes (a dataframe for each locus)
kirAlleleDFList <- read.kir_allele_dataframe_from_reference_fasta(kirReferenceAlignedFasta, kirLocusList)

## Check to make sure bowtie2is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')

#currentSample <- sampleList[[5]]

for(currentSample in sampleList){
  
if(!grepl('IND',currentSample$name)){
  next
}
  
### Fill in the path to the alignment file (it may or may not be present)
currentSample$gcSamPath <- file.path(resultsDirectory,paste0(currentSample$name,'.sam'))

## Run the alignment
sampleAlign <- run.bowtie2_gc_alignment(bowtie2, kirReferenceIndex, threads, currentSample, resultsDirectory)

## Read in the SAM file to analyze where the reads are aligning
samTable <- read.bowtie2_sam_nohd(currentSample$gcSamPath)

## Determine which loci are present for this sample
currentPresentLoci <- colnames(kffPresenceDF)[kffPresenceDF[currentSample$name,]>0]
currentPresentLoci <- currentPresentLoci[currentPresentLoci %in% kirLocusList]

#copyNumberList <- allCopyNumberList[[currentSample$name]]
copyNumberList <- as.list(copyNumberDF[currentSample$name,])

## Initialize a list of dataframes for storing the read information
assembledNucList <- build.assembled_nuc_list(kirLocusList, kirAlleleDFList)

## Build up the deletion index list to speed up read assignment
#deletionIndexDFList <- build.deletion_index_list(kirAlleleList, kirAlleleDFList)

## Build up the inverse deletion index list to speed up read assignment
#inverseDeletionIndexDFList <- build.inverse_deletion_index_list(kirAlleleList, kirAlleleDFList)

## Pull out the unique read names
uniqueReadNames <- unique(samTable$read_name)

## Figure out which reads map to multiple loci / map to absent loci / have unpaired mappings / have paired mappings
readAssignmentList <- run.sam_read_assignment(samTable, uniqueReadNames, currentPresentLoci)

## Assign paired-end reads to the assembledNucList
#assembledNucList <- run.assemble_paired_reads(samTable, readAssignmentList$singleLocusPairedReads, assembledNucList, deletionIndexDFList, kirAlleleDFList)

## Assign unpaired reads to the assembledNucList
#assembledNucList <- run.assemble_unpaired_reads(samTable, readAssignmentList$singleLocusUnpairedReads, assembledNucList, deletionIndexDFList, inverseDeletionIndexDFList, kirAlleleDFList, nucListConv, kirLocusList)


new.run.assemble_unpaired_reads <- function(currentReadName, samTable, deletionIndexDFList, inverseDeletionIndexDFList){
  samSubsetTable <- samTable[read_name == currentReadName]
  currentReadLocus <- samSubsetTable$locus[[1]]
  
  return(new.build.add_read_to_assembly(samSubsetTable, deletionIndexDFList, inverseDeletionIndexDFList, currentReadLocus))
}

new.build.add_read_to_assembly <- function(samSubsetTable, deletionIndexDFList, inverseDeletionIndexDFList, currentReadLocus){
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


cat('\nProcessing',length(readAssignmentList$singleLocusPairedReads),'single locus paired reads...')
pairedAssemblyList <- lapply(readAssignmentList$singleLocusPairedReads, new.run.assemble_paired_reads, samTable, deletionIndexDFList, inverseDeletionIndexDFList)
pairedAssemblyUnlistedList <- list()
i = 1
for(subList in pairedAssemblyList){
  pairedAssemblyUnlistedList[[i]] <- subList[1]
  i = i + 1
  pairedAssemblyUnlistedList[[i]] <- subList[2]
  i = i + 1
}
cat('\nFinished.')

cat('\nProcessing',length(readAssignmentList$singleLocusUnpairedReads),'single locus unpaired reads...')
unpairedAssemblyList <- lapply(readAssignmentList$singleLocusUnpairedReads, new.run.assemble_unpaired_reads, samTable, deletionIndexDFList, inverseDeletionIndexDFList)
cat('\nFinished.')

emptyAssembledNucList <- assembledNucList

cat('\nAdding paired end reads to assembly...')
assembledNucList <- run.condense_assembly_list(pairedAssemblyUnlistedList, assembledNucList, kirLocusList, emptyAssembledNucList)
cat('\nFinished.')

cat('\nAdding unpaired end reads to assembly...')
assembledNucList <- run.condense_assembly_list(unpairedAssemblyList, assembledNucList, kirLocusList, emptyAssembledNucList)
cat('\nFinished.')


run.condense_assembly_list <- function(assemblyReadList, assembledNucList, kirLocusList, emptyAssembledNucList, nucListConv){
  
  for(currentLocus in kirLocusList){
    cat('\nProcessing',currentLocus)
    currentLocusAssembledNucDT <- assembledNucList[[currentLocus]]
    
    currentLocusIndex <- lapply(assemblyReadList,function(x, currentLocus){return(names(x)==currentLocus)}, currentLocus)
    currentLocusAssemblyList <- assemblyReadList[unlist(currentLocusIndex)]
    lapply(currentLocusAssemblyList[1:10], run.convert_assembly_list_to_frames, currentLocusAssembledNucDT, nucListConv)
    
    for(subFrame in convertedCurrentLocusAssemblyFrames){
      currentLocusAssembledNucList[,colnames(subFrame)] <- currentLocusAssembledNucList[,colnames(subFrame)] + subFrame
    }
    
    assembledNucList[[currentLocus]] <- currentLocusAssembledNucList
  }
  return(assembledNucList)
}

run.convert_assembly_list_to_frames <- function(singleListElement, currentLocusAssembledNucDT, nucListConv){
  singleListElement <- singleListElement[[1]]
  
  for(nuc in names(nucListConv)){
    currentNucElementList <- singleListElement[singleListElement == nuc]
    
    colIndex <- names(currentNucElementList)
    
    if(length(colIndex) == 0 ){
      next
    }
    
    currentLocusAssembledNucDT[nucListConv[nuc][[1]],colIndex] <- currentLocusAssembledNucDT[nucListConv[nuc][[1]],colIndex,with=FALSE] + 1
  }
  
  return()
}





## Determine which coordinates of the assembly pass the depthThreshold
goodBadCoordList <- build.good_coord_bad_coord_list(kirLocusList,assembledNucList,depthThreshold)

## Convert above threshold positions to their nucleotides
allLocusNucAssembly <- build.nuc_assembly_frame(assembledNucList,kirLocusList,goodBadCoordList,hetRatio)

#copyNumberList <- list('KIR3DL3'=2,'KIR2DL3'=2,'KIR2DP1'=2,'KIR2DL1'=2,'KIR3DP1'=2,'KIR2DL4'=2,
#                       'KIR3DL1'=1,'KIR3DS1'=1,'KIR2DL5'=1,'KIR2DS5'=1,'KIR2DS1'=1,'KIR2DS4'=1,'KIR3DL2'=2)

## Figure out the closest matching alleles for each locus
allGoodAlleles <- run.find_allele_matches_for_assembly(allLocusNucAssembly,currentPresentLoci,kirAlleleDFList,copyNumberList)

## Use the good allele information to assign multi locus reads to a single locus
assembledNucList <- run.assemble_multi_reads(samTable,readAssignmentList$multiLocusReads,assembledNucList,deletionIndexDFList,inverseDeletionIndexDFList,kirAlleleDFList,nucListConv,allGoodAlleles,kirLocusList)

## Determine which coordinates of the assembly pass the depthThreshold
goodBadCoordList <- build.good_coord_bad_coord_list(kirLocusList,assembledNucList,depthThreshold)

## Convert above threshold positions to their nucleotides
allLocusNucAssembly <- build.nuc_assembly_frame(assembledNucList,kirLocusList,goodBadCoordList,hetRatio)

## Again figure out the closest matching alleles for each locus
allGoodAlleles <- run.find_allele_matches_for_assembly(allLocusNucAssembly,currentPresentLoci,kirAlleleDFList,copyNumberList)

## Use the good alleles to fill in the rest of the assembly
allLocusNucAssembly <- build.fill_in_remaining_positions_with_good_alleles(currentPresentLoci,allGoodAlleles,kirAlleleDFList,goodBadCoordList,allLocusNucAssembly,copyNumberList)

currentSample <- run.convert_assembly_to_fasta(currentPresentLoci, copyNumberList, allLocusNucAssembly, currentSample, assembledReferenceDirectory)

currentDelIndex <- run.read_list(currentSample$assembledDelIndex)



## Check to make sure bowtie2-build is accessible
bowtie2Build <- system2('which', c('bowtie2-build'), stdout=T, stderr=T)
check.system2_output(bowtie2Build, 'bowtie2-build not found')

## Check to make sure bowtie2is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')

## Creqte a bowtie2 index for the kir_reference.fasta file
createIndex <- system2(bowtie2Build, c(currentSample$assembledRefPath, currentSample$assembledRefIndex))
check.system2_output(createIndex, 'bowtie2 index building failed')

sampleAlign <- run.bowtie2_assembly_alignment(bowtie2,currentSample$assembledRefIndex,8,currentSample,assembledReferenceDirectory)


### 2. Convert SAM to BAM file
samtParam <- "-b "
bamFile <- file.path(assembledReferenceDirectory, paste0(currentSample$name, ".bam"))  # BAM output file
samtOut <- paste0("-o ", bamFile)
cat("samtools view", samtParam, currentSample$assembledSamPath, samtOut)
system2("samtools", c("view",samtParam, currentSample$assembledSamPath, samtOut))
cat("\n\n")

### 3. Sort BAM file
samtT <- "-T temp "
sortedBamFile <- file.path(assembledReferenceDirectory, paste0(currentSample$name, ".sorted.bam"))   # Sorted BAM output file
samtOutSorted <- paste0("-o ", sortedBamFile)
cat("samtools sort", samtOutSorted, samtT, bamFile)
system2("samtools", c("sort", samtOutSorted, samtT, bamFile))
cat("\n\n")

### 4. Generate VCF file
mpM <- "-m 0"
mpF <- "-F 0.0002"
mpU <- "-u"
mpf <- paste0("-f ", currentSample$assembledRefPath)
mpL <- paste0("-l ", currentSample$assembledBedPath)

vcfFilepath <- file.path(assembledReferenceDirectory,paste0(currentSample$name,'assembled.vcf'))
mpo <- paste0("-o ", vcfFilepath)
mpO <- "-O v"

cat("samtools", "mpileup", mpM, mpF, mpU, mpf, sortedBamFile, mpL, "|", "bcftools", "call", "--multiallelic-caller", mpO, mpo)
system2("samtools", c("mpileup", mpM, mpF, mpU, mpf, sortedBamFile, mpL, "|", "bcftools", "call", "--multiallelic-caller", mpO, mpo))
cat("\n\n")
}

currentVcf <- as.data.table(read.table(vcfFilepath,stringsAsFactors=F,check.names=F))



currentLocus <- 'KIR3DL3'
#belowDepthCoords <- colnames(assembledNucListPaired[[currentLocus]])[apply(assembledNucListPaired[[currentLocus]], 2, sum) < 1]


positionDepth <- apply(assembledNucList[[currentLocus]], 2, sum)
positionDominantDepth <- apply(assembledNucList[[currentLocus]], 2, max)

positionCoord <- colnames(assembledNucList[[currentLocus]])

plot(positionCoord, positionDominantDepth, main=paste(currentLocus, 'depth per position'),
     xlab='Position',ylab='Depth',pch=20,cex=0.1, col=2)
points(positionCoord, positionDepth, pch=20, cex=0.1, col=1)



#currentLocusAssembly <- assembledSeqList[[currentLocus]]
#currentLocusAssembly <- currentLocusAssembly[apply(currentLocusAssembly != 0, 1, any),]

#frameNRows <- nrow(currentLocusAssembly)

#filledPositions <- colnames(currentLocusAssembly)[apply(currentLocusAssembly[seq(1,frameNRows,2),] != 0, 2, any)]
#emptyPositions <- colnames(currentLocusAssembly)[apply(currentLocusAssembly[seq(1,frameNRows,2),] == 0, 2, all)]

#varReadsPerPos <- apply(currentLocusAssembly[seq(1,frameNRows,2),] != 0, 2, sum)
#totalDepthPerPos <- apply(apply(currentLocusAssembly[seq(2,frameNRows,2),],2,as.integer),2,sum)

#positionCoord <- colnames(currentLocusAssembly)

#varPos <- colnames(currentLocusAssembly)[sapply(sapply(apply(currentLocusAssembly[seq(1,frameNRows,2),],2,unique),setdiff,c(0)),length)>1]


#plot(positionCoord, varReadsPerPos, main=paste(currentLocus, "variable reads per position"),
#     xlab="Position", ylab="Depth", pch=20, cex=0.1)

#plot(positionCoord, totalDepthPerPos, main=paste(currentLocus, "depth per position"),
#     xlab="Position", ylab="Depth", pch=20, cex=0.1)

variantPoints <- as.integer(variantPositions)
variantPoints[1:length(variantPoints)] <- as.integer(max(positionDepth)/3*2)
points(as.integer(variantPositions), variantPoints, col=2, pch=20, cex=0.3)

divergentPoints <- as.integer(divergentPositions)
divergentPoints[1:length(divergentPoints)] <- as.integer(max(positionDepth)/3*1)
points(as.integer(divergentPositions), divergentPoints, col=2, pch=20, cex=0.3)





