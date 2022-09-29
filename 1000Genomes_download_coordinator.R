'
This script downloads 30X WGS data from the 1000Genomes project and converts it into FQ data for input into PING
PING must be run as a seperate step after downloading and converting the data.

It is advised to run PING on only 1 population at a time, because the allelic makeup differences between
populations can skew copy number thresholding.

Example use of this script is is set data_outputDir to whereever you would like the fq data to be sent,
setting the number of threads to match your system, then running the full script.
'
# RSTUDIO / RSCRIPT Initialization variables ------------------------------------------------
data_outputDir <- 'tgpData/' # <-- set this to whatever path you would like the data to be output to
threads <- 40                # <-- Set this to the maximum number of threads you want the script to use

## --- ##
'After setting the above two variables, run the full script to download the FQ data'
## --- ##

# ---- DEPENDENCIES ----
' if any dependencies are missing, install with
install.packages("plotly",dependencies = T)
'
library(parallel)
source('Resources/general_functions.R')
source('Resources/tgp_functions.R')


# Preparation ------------------------------------------------------------------------------
dir.create(data_outputDir, showWarnings = FALSE)
threads <- min(c(detectCores(),threads))
resourceDir <- 'Resources/1000Genomes_resources/'
awsUrl <- 'http://s3.amazonaws.com/1000genomes'


# Read in the data manifests ---------------------------------------------------------------
tgpManifest <- read.table(file.path(resourceDir,'manifest.tsv'), sep='\t',stringsAsFactors = F, header=T, check.names=F)
urlManifest <- read.table(file.path(resourceDir,'tgp_full_30x.tsv'), sep='\t',stringsAsFactors = F,col.names=c('Sample name','url'), check.names=F)
dataManifest <- merge(tgpManifest, urlManifest, by="Sample name")
rownames(dataManifest) <- dataManifest$`Sample name`


# Download the reference if needed ---------------------------------------------------------
refPath <- retrieve_ref(resourceDir)
bedPath <- file.path(resourceDir,'kir_regions.bed')


# Download 1000Genomes CRAMs and convert to FQ  --------------------------------------------
pop_vect <- unique(dataManifest$`Population code`)
for( pop in pop_vect ){
  popData_dir <- file.path(data_outputDir,pop)
  dir.create(popData_dir,showWarnings = F)
  
  popManifest <- dataManifest[dataManifest$`Population code` == pop,]
  
  popData_cramDir <- file.path(popData_dir,'cram')
  dir.create(popData_cramDir,showWarnings = F)
  popData_fqDir <- file.path(popData_dir,'fq')
  dir.create(popData_fqDir,showWarnings = F)
  
  cat(paste('\n\n--> Starting download of',pop,'KIR aligned reads to',popData_cramDir,'<--'))
  cat(paste('\nThere are',length(popManifest$`Sample name`),'samples being downloaded.'))
  cat(paste('\nAfter download they will be converted to fq format.\n'))
  fq_path_list <- mclapply(popManifest$`Sample name`, function(sampleID){
    cram_path <- retrieve_cram(sampleID, popData_cramDir, refPath, bedPath, resourceDir)
    fq_paths <- cram_to_fq(sampleID, cram_path, popData_fqDir)
    file.remove(paste0(sampleID,'.final.cram.crai'))
    return(fq_paths)
  }, mc.cores=threads, mc.silent=F)
  cat('\n\n-- DONE! --')
}
