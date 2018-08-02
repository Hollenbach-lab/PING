library(data.table)
library(ggplot2)
library(stringr)
library(methods)

source('Resources/gc_functions.R')
setwd('/home/wmarin/PING_projects/PING2/')

########## INPUT variables
sampleDirectory <- '/home/wmarin/african_samples/2_extracted_kir/'
fastqPattern <- 'fastq'
threads <- 10
resultsDirectory <- 'african_filled_one_mismatch_3DL3_limit_1000_kir_results'
###########
kirLocusList <- c('KIR3DP1','KIR2DL5A','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3','KIR2DL5B',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

cat('Current working directory: ', getwd(),'\n')

sampleDirectory <- normalizePath(sampleDirectory, mustWork=T)
gcResourceDirectory <- normalizePath('Resources/gc_resources', mustWork = T)

#kirReferenceFasta <- normalizePath(file.path(gcResourceDirectory,'subset_kir_reference','subset_kir_reference.fasta'), mustWork=T)
#kirReferenceIndex <- file.path(gcResourceDirectory,'subset_kir_reference','subset_kir_reference')

kirReferenceFasta <- normalizePath(file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled.fasta'), mustWork=T)
kirReferenceIndex <- file.path(gcResourceDirectory,'filled_kir_reference','KIR_gen_onelines_filled')

kirAlleleList <- read.kir_allele_list_from_reference_fasta(kirReferenceFasta)
kirAlleleListRes3 <- unique(unlist(lapply(kirAlleleList, kir.allele_resolution, 3)))

## Check to make sure bowtie2-build is accessible
bowtie2Build <- system2('which', c('bowtie2-build'), stdout=T, stderr=T)
check.system2_output(bowtie2Build, 'bowtie2-build not found')

## Creqte a bowtie2 index for the kir_reference.fasta file
#createIndex <- system2(bowtie2Build, c(fullKirReferenceFasta, fullKirReferenceIndex))
#check.system2_output(createIndex, 'bowtie2 index building failed')

## Building a list of sample objects from files in sampleDirectory that match fastqPattern
sampleList <- build.paired_sample_objects(sampleDirectory,fastqPattern)

## Check to make sure bowtie2is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')

locusCountDF <- data.frame(matrix(0, length(sampleList), length(kirLocusList)),row.names=names(sampleList),check.names=F,stringsAsFactors=F)
colnames(locusCountDF) <- kirLocusList

alleleCountDF <- data.frame(matrix(0, length(sampleList), length(kirAlleleListRes3)), row.names=names(sampleList),check.names=F,stringsAsFactors=F)
colnames(alleleCountDF) <- kirAlleleListRes3

## Run all samples through bowtie2 gc alignment (for now just testing with a couple)
for(currentSample in sampleList[1:length(sampleList)]){
  cat('\n\nProcessing', currentSample$name)
  #sampleAlign <- run.bowtie2_gc_alignment(bowtie2, kirReferenceIndex, threads, currentSample)
  currentSample$gcSamPath <- file.path(resultsDirectory,paste0(currentSample$name,'.sam'))
  
  cat("\n\nReading in",currentSample$gcSamPath)
  samTable <- read.bowtie2_sam_nohd(currentSample$gcSamPath)
  countList <- run.determine_kir_presence(currentSample, samTable, max_3DL3=1000)
  
  locusCountDF[currentSample$name,names(countList$locusMatches)] = countList$locusMatches
  alleleCountDF[currentSample$name,names(countList$alleleMatches)] = countList$alleleMatches
}

write.csv(locusCountDF, file = file.path(resultsDirectory, 'locusCountFrame.csv'))
write.csv(alleleCountDF, file = file.path(resultsDirectory, 'alleleCountFrame.csv'))

### Accuracy testing code
goodRows <- rownames(locusCountDF)[apply(locusCountDF, 1, function(x) x['KIR3DL3']>100)]

locusRatioDF <- apply(locusCountDF[goodRows,], 2, function(x) x / locusCountDF[goodRows,'KIR3DL3'])
locusRatioDF <- as.data.frame(locusRatioDF)

currentLocus <- 'KIR2DL3'
for(currentLocus in kirLocusList){

  threshold <- 0.1
  
  ggplot(locusRatioDF, aes(x=1:nrow(locusRatioDF), y=locusRatioDF[order(locusRatioDF[,currentLocus]),currentLocus])) + 
    geom_point() +
    scale_y_continuous(name=currentLocus, limits=c(0, max(locusRatioDF[,currentLocus])))
  
  ggsave(file.path(resultsDirectory, paste0(currentLocus, '_plot.pdf')), width=6.72, height=6.72)
}

absentSamples <- rownames(locusRatioDF)[locusRatioDF[,currentLocus] < threshold]
pingAbsentSamples <- rownames(pingGCDF)[pingGCDF[,unlist(str_split(currentLocus, 'KIR'))[2]] == 0]

write.csv(locusCountDF, file = file.path(resultsDirectory, 'locusCountFrame.csv'))
write.csv(alleleCountDF, file = file.path(resultsDirectory, 'alleleCountFrame.csv'))

#pingGCDF <- read.csv(file = file.path(resultsDirectory, 'GC_table_all_batches.csv'), check.names=F)
#rownames(pingGCDF) <- unlist(lapply(rownames(pingGCDF), function(x) paste0(x, '_')))
#inBoth <- intersect(rownames(locusRatioDF), rownames(pingGCDF))
#pingGCDF <- pingGCDF[inBoth,]
#locusRatioDF <- locusRatioDF[inBoth,]
