library(data.table)
library(ggplot2)
library(stringr)
library(methods)
library(plotly)

source('Resources/gc_functions.R')


########## INPUT variables
setwd('/home/wmarin/PING_projects/PING2/')
africanResultsDirectory <- '/home/wmarin/PING_projects/PING2/african_filled_one_mismatch_kir_results/'
indigoResultsDirectory <- '/home/wmarin/PING_projects/PING2/indigo_filled_one_mismatch_kir_results/'
fastqPattern <- 'fastq'
resultsDirectory <- 'kff_threshold_graphs'
KIR3DL3MinReadThreshold <- 100
maxReadThreshold <- 30000
###########

kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

## Make sure the needed directories can be accessed
africanResultsDirectory <- normalizePath(africanResultsDirectory, mustWork=T)
indigoResultsDirectory <- normalizePath(indigoResultsDirectory, mustWork=T)
resultsDirectory <- normalizePath(resultsDirectory, mustWork=T)

## Define the file paths to the normalized KFF frames
africanKffFile <- file.path(africanResultsDirectory, 'kffNormFrame.csv')
indigoKffFile <- file.path(indigoResultsDirectory, 'kffNormFrame.csv')

## Feed in the normalized KFF frames
africanKffRatioDF <- read.csv(africanKffFile, stringsAsFactors = F, check.names = F, row.names = 1)
indigoKffRatioDF <- read.csv(indigoKffFile, stringsAsFactors = F, check.names = F, row.names = 1)

## Define the file paths to the read mapping count frames
africanLocusCountFile <- file.path(africanResultsDirectory, 'locusCountFrame.csv')
indigoLocusCountFile <- file.path(indigoResultsDirectory, 'locusCountFrame.csv')

## Read in the read count frames
africanLocusCountDF <- read.csv(africanLocusCountFile, stringsAsFactors = F, check.names = F, row.names = 1)
indigoLocusCountDF <- read.csv(indigoLocusCountFile, stringsAsFactors = F, check.names = F, row.names = 1)

## Define good sample ID's as those having over or equal to 'KIR3DL3MinReadThreshold' read counts for KIR3DL3
africanGoodSamples <- rownames(africanLocusCountDF)[apply(africanLocusCountDF, 1, function(x) x['KIR3DL3']>=KIR3DL3MinReadThreshold)]
indigoGoodSamples <- rownames(indigoLocusCountDF)[apply(indigoLocusCountDF, 1, function(x) x['KIR3DL3']>=KIR3DL3MinReadThreshold)]

## Define bad sample ID's as those having under 'KIR3DL3MinReadThreshold' read counts for KIR3DL3
africanBadSamples <- rownames(africanLocusCountDF)[apply(africanLocusCountDF, 1, function(x) x['KIR3DL3']<KIR3DL3MinReadThreshold)]
indigoBadSamples <- rownames(indigoLocusCountDF)[apply(indigoLocusCountDF, 1, function(x) x['KIR3DL3']<KIR3DL3MinReadThreshold)]

## Set the good sample ID's as the ones that will be used for further analysis
africanSampleIds <- africanGoodSamples
indigoSampleIds <- indigoGoodSamples

## Keep track of what samples are being discarded
cat('\nSkipping', length(africanBadSamples), 'KhoeSan samples that had fewer than',KIR3DL3MinReadThreshold,'KIR3DL3 reads.')
cat('\nSkipping', length(indigoBadSamples), 'Indigo samples that had fewer than',KIR3DL3MinReadThreshold,'KIR3DL3 reads.')

## Subset the KFF ratio dataframe by the good sample ID's
africanKffRatioDF <- africanKffRatioDF[africanSampleIds,]
indigoKffRatioDF <- indigoKffRatioDF[indigoSampleIds,]

## For each KIR locus, print a graph of the KFF ratio's for INDIGO and KhoeSan samples
for(currentLocus in kirLocusList){
  if(!(currentLocus %in% colnames(africanKffRatioDF))){
    next
  }
  
  ## Determine the rank of the ratios (x-axis order from lowest to highest)
  africanKffRatioDF$ratioRank <- rank(africanKffRatioDF[,currentLocus], ties.method = 'first')
  indigoKffRatioDF$ratioRank <- rank(indigoKffRatioDF[,currentLocus], ties.method = 'first')
  
  p <- plot_ly() %>%
    add_trace(x=africanKffRatioDF[africanSampleIds,'ratioRank'],
              y=africanKffRatioDF[africanSampleIds,currentLocus],
              mode='markers',type='scatter',name='KhoeSan') %>%
    add_trace(x=indigoKffRatioDF[indigoSampleIds,'ratioRank'],
              y=indigoKffRatioDF[indigoSampleIds,currentLocus],
              mode='markers',type='scatter',name='indigo') %>%
    layout(title=currentLocus,
           xaxis = list(title='Sample rank'),
           yaxis = list(title=paste(currentLocus,'/ KIR3DL3 Ratio'),rangemode='tozero'))
    print(p)
    
    htmlwidgets::saveWidget(p, file=file.path(resultsDirectory,paste0(currentLocus, '_kff_plot.html')))
}



