library(data.table)
library(stringr)
library(methods)
library(pryr)
library(plotly)

source('Resources/gc_functions.R')

kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

gcResourceDirectory <- normalizePath('Resources/gc_resources', mustWork = T)

kirReferenceAlignedFasta <- normalizePath(file.path(gcResourceDirectory, 'filled_kir_reference', 'KIR_gen_onelines_aligned_filled.fasta'), mustWork=T)

## Read in the aligned fasta and convert it into a list of dataframes (a dataframe for each locus)
kirAlleleDFList <- read.kir_allele_dataframe_from_reference_fasta(kirReferenceAlignedFasta, kirLocusList)

