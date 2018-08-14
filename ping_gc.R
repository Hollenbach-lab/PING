library(data.table)
library(ggplot2)
library(stringr)
library(methods)
library(plotly)


source('Resources/gc_functions.R')


########## INPUT variables
#setwd('/home/wmarin/PING_projects/PING2/')
#sampleDirectory <- '/home/wmarin/african_samples/2_extracted_kir/'
#fastqPattern <- 'fastq'
#threads <- 10
#resultsDirectory <- 'indigo_filled_one_mismatch_kir_results'
#KIR3DL3MinReadThreshold <- 100
#maxReadThreshold <- 30000
#probelistFile <- 'probelist_2018_08_02.csv'
###########

ping_gc <- function(sampleDirectory='/home/common_arse/INDIGO/2_KIR_runs_extracted/indigo_plates1-4_extracted_8-5-17/',
                    fastqPattern='fastq',
                    threads=24,
                    resultsDirectory='indigo_filled_one_mismatch_kir_results',
                    KIR3DL3MinReadThreshold=100,
                    maxReadThreshold=30000,
                    probelistFile='probelist_2018_08_02.csv'){

kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

cat('Current working directory: ', getwd(),'\n')

resultsDirectory <- normalizePath(resultsDirectory, mustWork=T)
sampleDirectory <- normalizePath(sampleDirectory, mustWork=T)
gcResourceDirectory <- normalizePath('Resources/gc_resources', mustWork = T)

kirReferenceFasta <- normalizePath(file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled.fasta'), mustWork=T)
kirReferenceIndex <- file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled')

kirAlleleList <- read.kir_allele_list_from_reference_fasta(kirReferenceFasta)
kirAlleleListRes3 <- unique(unlist(lapply(kirAlleleList, kir.allele_resolution, 3)))
  
## Read in the probelist CSV file as a dataframe
cat('\n\nReading in the KFF probelist file: ', file.path(gcResourceDirectory, probelistFile))
probeDF <- read.csv(file.path(gcResourceDirectory, probelistFile), stringsAsFactors = F, check.names = F)
row.names(probeDF) <- probeDF$Name
  
## Pull out all the probe names with '>'
kffPresenceProbeNameList <- grep('>', probeDF$Name, fixed=T, value=T)
  
## Split the pulled out probe names by '>', then grab all the unique locus names
kffLociList <- unique(tstrsplit(kffPresenceProbeNameList,'>',fixed=T)[[2]])
  
## Check to make sure bowtie2-build is accessible
#bowtie2Build <- system2('which', c('bowtie2-build'), stdout=T, stderr=T)
#check.system2_output(bowtie2Build, 'bowtie2-build not found')
  
## Creqte a bowtie2 index for the kir_reference.fasta file
#createIndex <- system2(bowtie2Build, c(fullKirReferenceFasta, fullKirReferenceIndex))
#check.system2_output(createIndex, 'bowtie2 index building failed')
  
## Building a list of sample objects from files in sampleDirectory that match fastqPattern
sampleList <- build.paired_sample_objects(sampleDirectory,fastqPattern)
  
## Check to make sure bowtie2is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')
  
  
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
alleleCountDF <- data.frame(matrix(0, length(sampleList), length(kirAlleleListRes3)), row.names=names(sampleList),check.names=F,stringsAsFactors=F)
colnames(alleleCountDF) <- kirAlleleListRes3
  
  
## Run all samples through bowtie2 gc alignment
for(currentSample in sampleList[1:length(sampleList)]){
  cat('\n\nProcessing', currentSample$name)
  cat('\n------------------------------------')
  
  cat('\n\nCounting KFF primer matches.')
    
  ## Count KFF probe matches for the currentSample
  kffCountList <- run.count_kff_probes(currentSample, probeDF, maxReadThreshold)
  kffCountDF[currentSample$name,names(kffCountList)] <- kffCountList 
    
  ## Write the results to a csv file
  write.csv(kffCountDF, file = file.path(resultsDirectory, 'kffCountFrame.csv'))
    
  cat('\nNormalizing KFF primer matches.')
    
  ## Normalize the KFF probe matches by KIR3DL3
  kffNormList <- run.reduce_and_normalize_kff_probes(kffCountList, kffLociList)
  kffNormDF[currentSample$name,names(kffNormList)] <- kffNormList
    
  ## Write the results to a csv file
  write.csv(kffNormDF, file = file.path(resultsDirectory, 'kffNormFrame.csv'))
  
  cat('\nDetermining KIR locus presence/absence')
  
  ## Determine locus presence/absence
  kffPresenceList <- run.kff_determine_presence_from_norm_values(kffNormList, kffThreshold=0.2)
  kffPresenceDF[currentSample$name,names(kffPresenceList)] <- kffPresenceList
  
  ## Write the results to a csv file
  write.csv(kffPresenceDF, file = file.path(resultsDirectory, 'kffPresenceFrame.csv'))
    
  cat('\n\nFinished with presence/absence determination, moving to copy number determination.')
    
  ## Fill in the path to the alignment file (it may or may not be present)
  currentSample$gcSamPath <- file.path(resultsDirectory,paste0(currentSample$name,'.sam'))
    
  ## If the alignment file does not exist, then run bowtie2 alignment, otherwise continue
  if(!file.exists(currentSample$gcSamPath)){
    cat('\n\nPerforming bowtie2 alignment for this sample.')
    sampleAlign <- run.bowtie2_gc_alignment(bowtie2, kirReferenceIndex, threads, currentSample)
  }else{
    cat('\n\nFound a previous alignment file for this sample, skipping bowtie2 alignment.')
  }
    
  cat("\nReading in",currentSample$gcSamPath)
    
  ## Read in the SAM file to analyze where the reads are aligning
  samTable <- read.bowtie2_sam_nohd(currentSample$gcSamPath)
    
  cat('\nCounting reads that align uniquely to a locus or allele ')
  
  ## Count how many reads align uniquely to a locus or allele
  countList <- run.count_kir_read_matches(currentSample, samTable, maxReadThreshold, kirLocusList)
  
  ## Add the counts to the appropriate count dataframe
  locusCountDF[currentSample$name,names(countList$locusMatches)] = countList$locusMatches
  alleleCountDF[currentSample$name,names(countList$alleleMatches)] = countList$alleleMatches
    
  ## Write the results to a csv file
  write.csv(locusCountDF, file = file.path(resultsDirectory, 'locusCountFrame.csv'))
  write.csv(alleleCountDF, file = file.path(resultsDirectory, 'alleleCountFrame.csv'))
}

cat('\n\n----- Finished with alignment! -----')
cat('\n\nMoving on to copy number graphing.')

## Read in the csv results
locusCountDF <- read.csv(file.path(resultsDirectory, 'locusCountFrame.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
kffPresenceDF <- read.csv(file.path(resultsDirectory, 'kffPresenceFrame.csv'), stringsAsFactors = F, check.names = F, row.names = 1)

## Only analyze samples that have at least 'KIR3DL3MinReadThreshold' number of unique KIR3DL3 reads
goodRows <- rownames(locusCountDF)[apply(locusCountDF, 1, function(x) x['KIR3DL3']>=KIR3DL3MinReadThreshold)]
  
## Keep track of what samples are being discarded
badRows <- rownames(locusCountDF)[apply(locusCountDF, 1, function(x) x['KIR3DL3']<KIR3DL3MinReadThreshold)]
cat('\nSkipping', length(badRows), 'samples that had fewer than',KIR3DL3MinReadThreshold,'KIR3DL3 reads.')
  
## Subset the count dataframe by the samples that were determined to be good, then normalize each locus unique read count by KIR3DL3
locusRatioDF <- apply(locusCountDF[goodRows,], 2, function(x) x / locusCountDF[goodRows,'KIR3DL3'])
locusRatioDF <- as.data.frame(locusRatioDF)

cat('\nGenerating copy number graphs... ')
run.generate_copy_number_graphs(locusRatioDF, kffPresenceDF, kirLocusList)
}
