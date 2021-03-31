library(data.table)
library(stringr)
library(methods)

#synSeqAnswerKey <- '/home/LAB_PROJECTS/PING2_PAPER/figure_scripts/synSeq_data_run3/synSeq.info.txt'
# Reads in answer key (synSeq.info.txt)
synSeq.readAnswerKey <- function(keyFile){
  
  out.list <- list()
  con  <- file(keyFile, open = "r")
  
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    lineVect <- unlist((strsplit(oneLine, "\t",fixed=T)))
    sampleID <- lineVect[1]
    genoID <- lineVect[2]
    genoStr <- lineVect[3]
    
    genoVect <- unlist(strsplit(genoStr,'_',fixed=T))
    out.list[[sampleID]] <- list('genoID'=genoID,'genoVect'=genoVect)
  }
  
  close(con)
  
  return(out.list)
}

alleleSetup.readAnswerKey <- function(keyFile){
  out.list <- list()
  con  <- file(keyFile, open = "r")
  
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    lineVect <- unlist((strsplit(oneLine, "\t",fixed=T)))
    sampleID <- lineVect[1]
    genoStr <- lineVect[2]
    
    genoVect <- unlist(strsplit(genoStr,'_',fixed=T))
    out.list[[sampleID]] <- list('genoVect'=genoVect)
  }
  
  close(con)
  
  return(out.list)
}

#synSeq.key <- synSeq.readAnswerKey(synSeqAnswerKey)

# rawFastqDirectory
# fastqPattern
# threads <- 60
# resultsDirectory
# shortNameDelim

# Alignment directory structure
# alignmentFiles
# - sampleName
#   - iterAlign
#     - iter_1...iter_5
#       - bed
#       - fasta
#       - vcf
#   - filterAlign
#     - [gene].vcf

# alignmentFileDirectory <- file.path(resultsDirectory,'alignmentFiles')
# if(!file.exists(alignmentFileDirectory)){
#   dir.create(alignmentFileDirectory)
# }

#referenceAlleleDF <- read.csv('Resources/genotype_resources/master_haplo_iteration_testing_v10.csv',row.names=1,stringsAsFactors = F)


### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
#bowtie2Build <- system2('which', c('bowtie2-build'), stdout=T, stderr=T)
#check.system2_output(bowtie2Build, 'bowtie2-build not found')

### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
#samtools <- system2('which', c('samtools'), stdout=T, stderr=T)
#check.system2_output(samtools, 'samtools')

### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
#bcftools <- system2('which', c('~/tools/samtools_update/bcftools/bcftools'), stdout=T, stderr=T)
#check.system2_output(bcftools, 'bcftools')

# load reference alleles for present loci into a sample object 'currentSample$refAlleleDF'
sampleObj.loadRefDF <- function(currentSample, referenceAlleleDF){
  '
  UPDATED FOR (ONLY GC / ONLY COPY) RUNNING
  '
  # Skip this sample if either gene content or copy number detemrination failed
  if( 'failed' %in% c(currentSample$geneContent, currentSample$copyNumber) ){
    
    currentSample[['refAlleleDF']] <- 'failed'
    
  }else{
    
    presentLociVect <- c()
    
    # pull out copy number and gene content info for current sample
    currentSampleCopy <- as.list(sapply(currentSample$copyNumber, as.numeric))
    currentSampleGC <- as.list(sapply(currentSample$geneContent, as.numeric))
    
    # Intersect all loci found across both copy and gene content
    #allCopyLoci <- intersect(names(currentSampleCopy),names(currentSampleGC)) # 2020/09/02 only GC fix
    allCopyLoci <- unique( c(names(currentSampleCopy), names(currentSampleGC)) )
    
    if( length(currentSampleCopy) > 0 ){
      # pull out all present loci determined by the copy module
      #presentCopyLoci <- allCopyLoci[currentSampleCopy[allCopyLoci] > 0]
      
      temp.allCopyLoci <- intersect(names(currentSampleCopy), allCopyLoci)
      
      presentLociVect <- c(presentLociVect, temp.allCopyLoci[currentSampleCopy[temp.allCopyLoci] > 0])
    }
    
    if( length(currentSampleGC) > 0 ){
      # pull out all present loci determined by the gc module
      #presentGCLoci <- allCopyLoci[currentSampleGC[allCopyLoci] > 0]
      
      temp.allCopyLoci <- intersect(names(currentSampleGC), allCopyLoci)
      
      presentLociVect <- c(presentLociVect, temp.allCopyLoci[currentSampleGC[temp.allCopyLoci] > 0])
    }
    # combine the two to cover all possible present loci
    #presentLociVect <- unique(c(presentCopyLoci,presentGCLoci))
    presentLociVect <- intersect( rownames(referenceAlleleDF), unique(presentLociVect) )
    
    if(length(presentLociVect) == 0){
      currentSample[['refAlleleDF']] <- 'failed'
    }else{
      # Load reference allele vect into the sample object
      currentSample[['refAlleleDF']] <- referenceAlleleDF[presentLociVect,]
    }
    
  }
  
  return(currentSample)
}

#UTRextList <- general.read_fasta('/home/LAB_PROJECTS/PING2_PAPER/PING2/Resources/genotype_resources/KIR_UTR_ext.fasta')

# Write fasta and bed files for dynamic reference building
sampleObj.writeRefFastaBed <- function(currentSample, locusRefList, alignmentFileDirectory){
  currentRefAlleleDF <- currentSample$refAlleleDF
  
  if('failed' %in% currentRefAlleleDF){
    currentSample[['iterRefDirectory']] <- 'failed'
    currentSample[['refIterVect']] <- 'failed'
    return(currentSample)
  }
  
  # pull out the number of iterations from the ref allele DF
  iterVect <- colnames(currentRefAlleleDF)
  
  # Create a directory for the current iteration reference files
  currentSampleResultsDirectory <- file.path(alignmentFileDirectory,currentSample$name,'iterAlign')
  if(!file.exists(currentSampleResultsDirectory)){
    dir.create(currentSampleResultsDirectory, recursive=T)
  }
  
  # For each iteration, write a fasta and bed file
  for(currentIter in iterVect){
    
    # Pull out ref alleles for current iteration
    currentIterDF <- currentRefAlleleDF[,currentIter,drop=F]
    
    # Create a directory for the current iteration reference files
    currentIterResultsDirectory <- file.path(currentSampleResultsDirectory,currentIter)
    if(!file.exists(currentIterResultsDirectory)){
      dir.create(currentIterResultsDirectory)
    }
    
    # Open fasta file connection
    fastaPath <- file.path(currentIterResultsDirectory, 'alleleReference.fasta')
    fastaCon <- file(fastaPath, 'w')
    
    # Open bed file connection
    bedPath <- file.path(currentIterResultsDirectory, 'alleleReference.bed')
    bedCon <- file(bedPath, 'w')
    
    # Pull out loci to create references for
    currentIterLociVect <- rownames(currentIterDF)
    
    # For each locus, write a reference allele and bed coordinates
    #for(currentLocus in currentIterLociVect){
    #  currentRefAllele <- currentRefAlleleDF[currentLocus,currentIter]
    for(currentRefAllele in refAlleleVect){  
      cat('\n',currentRefAllele)
      
      currentLocus <- strsplit(currentRefAllele,'*',fixed=T)[[1]][1]
      
      if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
        currentLocus <- 'KIR2DL5'
      }
      
      # Remove deletion positions from the allele string
      noDelAlleleStr <- gsub('.','',locusRefList[[currentLocus]]$alleleStrList[[currentRefAllele]],fixed = T)
      
      # Add in 5UTR extension (by replacing current 5UTR with 1000bp 5UTR)
      UTRextName <- paste0(currentLocus,'_5UTR')
      UTRextStr5 <- UTRextList[[UTRextName]]
      noDelAlleleStr <- gsub(paste0(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`5UTR`$featSeq,
                                    locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`E1`$featSeq),
                             paste0(UTRextStr5,
                                    locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`E1`$featSeq),
                             noDelAlleleStr,
                             fixed=T)
      
      # Add in 3UTR extension (by replacing current 5UTR with 1000bp 5UTR)
      lastExonLab <- names(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]])[ (length(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]])-1) ]
      UTRextName <- paste0(currentLocus,'_3UTR')
      UTRextStr3 <- UTRextList[[UTRextName]]
      noDelAlleleStr <- gsub(paste0(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]][[lastExonLab]]$featSeq,
                                    locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`3UTR`$featSeq),
                             paste0(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]][[lastExonLab]]$featSeq,
                                    UTRextStr3),
                             noDelAlleleStr,
                             fixed=T)
      
      # Replace all '*' with 'N'
      noDelAlleleStr <- gsub('*','N',noDelAlleleStr,fixed=T)
      
      # Pull out bed coordinate list
      alleleBed <- locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]
      
      # Modify allele bed to account for 5' and 3' UTR extensions
      utr5ExtraLen <- nchar(UTRextStr5) - nchar(alleleBed$`5UTR`$featSeq)
      utr3ExtraLen <- nchar(UTRextStr3) - nchar(alleleBed$`3UTR`$featSeq)
      
      # Count the number of ambiguous characters in this reference allele
      ambChrCount <- str_count(noDelAlleleStr, pattern=fixed('N'))
      if(ambChrCount > 0){
        message('\nAmbiguous characters found in reference allele\nAllele ID:\t',currentRefAllele,
                '\nNumber:\t',ambChrCount)
      }
      
      # Catch if the last bed coord and the no del reference allele length do not match
      if( (as.numeric(alleleBed$`3UTR`$endPos) + utr5ExtraLen + utr3ExtraLen) != nchar(noDelAlleleStr) ){
        message('\nMismatch between BED coordinates to reference allele length\nAllele ID:\t',currentRefAllele,
                '\nlast BED coord:\t',alleleBed$`3UTR`$endPos,
                '\nallele length:\t',nchar(noDelAlleleStr),
                '\niter:\t\t',currentIter)
        close(bedCon)
        close(fastaCon)
        stop()
      }
      
      # Write allele seq to fasta file
      general.write_fasta(fastaCon, currentRefAllele, noDelAlleleStr)
      
      # Write feature coordinates to bed file
      sapply(alleleBed, function(x){
        general.write_bed(bedCon, x$alleleName, x$startPos, x$endPos, x$featName, utr5ExtraLen, utr3ExtraLen)
      })
      
    }
    
    # close connections
    close(bedCon)
    close(fastaCon)
    
    currentSample[['refIterVect']] <- iterVect
  }
  
  currentSample[['iterRefDirectory']] <- currentSampleResultsDirectory
  
  return(currentSample)
}

# Write fasta and bed files for dynamic reference building
new.sampleObj.writeRefFastaBed <- function(currentSample, locusRefList, alignmentFileDirectory, synSeq.key){
  currentRefAlleleDF <- currentSample$refAlleleDF
  
  if('failed' %in% currentRefAlleleDF){
    currentSample[['iterRefDirectory']] <- 'failed'
    currentSample[['refIterVect']] <- 'failed'
    return(currentSample)
  }
  
  # pull out the number of iterations from the ref allele DF
  #iterVect <- colnames(currentRefAlleleDF)
  
  # Create a directory for the current iteration reference files
  currentSampleResultsDirectory <- file.path(alignmentFileDirectory,currentSample$name,'iterAlign')
  if(!file.exists(currentSampleResultsDirectory)){
    dir.create(currentSampleResultsDirectory, recursive=T)
  }
  
  refAlleleVect <- synSeq.key[[currentSample$name]]$genoVect
  iterVect <- 'iter_1'
  
  # For each iteration, write a fasta and bed file
  for(currentIter in iterVect){
    
    # Pull out ref alleles for current iteration
    #currentIterDF <- currentRefAlleleDF[,currentIter,drop=F]
    
    # Create a directory for the current iteration reference files
    currentIterResultsDirectory <- file.path(currentSampleResultsDirectory,currentIter)
    if(!file.exists(currentIterResultsDirectory)){
      dir.create(currentIterResultsDirectory)
    }
    
    # Open fasta file connection
    fastaPath <- file.path(currentIterResultsDirectory, 'alleleReference.fasta')
    fastaCon <- file(fastaPath, 'w')
    
    # Open bed file connection
    bedPath <- file.path(currentIterResultsDirectory, 'alleleReference.bed')
    bedCon <- file(bedPath, 'w')
    
    # Pull out loci to create references for
    #currentIterLociVect <- rownames(currentIterDF)
    
    # For each locus, write a reference allele and bed coordinates
    #for(currentLocus in currentIterLociVect){
    #  currentRefAllele <- currentRefAlleleDF[currentLocus,currentIter]
    for(currentRefAllele in refAlleleVect){  
      cat('\n',currentRefAllele)
      
      currentLocus <- strsplit(currentRefAllele,'*',fixed=T)[[1]][1]
      
      if(currentLocus == 'KIR2DL5A' | currentLocus == 'KIR2DL5B'){
        currentLocus <- 'KIR2DL5'
      }
      
      # Remove deletion positions from the allele string
      noDelAlleleStr <- gsub('.','',locusRefList[[currentLocus]]$alleleStrList[[currentRefAllele]],fixed = T)
      
      # Add in 5UTR extension (by replacing current 5UTR with 1000bp 5UTR)
      #UTRextName <- paste0(currentLocus,'_5UTR')
      #UTRextStr5 <- UTRextList[[UTRextName]]
      #noDelAlleleStr <- gsub(paste0(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`5UTR`$featSeq,
      #                              locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`E1`$featSeq),
      #                       paste0(UTRextStr5,
      #                              locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`E1`$featSeq),
      #                       noDelAlleleStr,
      #                       fixed=T)
      
      # Add in 3UTR extension (by replacing current 5UTR with 1000bp 5UTR)
      #lastExonLab <- names(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]])[ (length(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]])-1) ]
      #UTRextName <- paste0(currentLocus,'_3UTR')
      #UTRextStr3 <- UTRextList[[UTRextName]]
      #noDelAlleleStr <- gsub(paste0(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]][[lastExonLab]]$featSeq,
      #                              locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]$`3UTR`$featSeq),
      #                       paste0(locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]][[lastExonLab]]$featSeq,
      #                              UTRextStr3),
      #                       noDelAlleleStr,
      #                       fixed=T)
      
      # Replace all '*' with 'N'
      noDelAlleleStr <- gsub('*','N',noDelAlleleStr,fixed=T)
      
      # Pull out bed coordinate list
      alleleBed <- locusRefList[[currentLocus]]$alleleBedList[[currentRefAllele]]
      
      # Modify allele bed to account for 5' and 3' UTR extensions
      #utr5ExtraLen <- nchar(UTRextStr5) - nchar(alleleBed$`5UTR`$featSeq)
      #utr3ExtraLen <- nchar(UTRextStr3) - nchar(alleleBed$`3UTR`$featSeq)
      
      # Count the number of ambiguous characters in this reference allele
      ambChrCount <- str_count(noDelAlleleStr, pattern=fixed('N'))
      if(ambChrCount > 0){
        message('\nAmbiguous characters found in reference allele\nAllele ID:\t',currentRefAllele,
                '\nNumber:\t',ambChrCount)
      }
      
      # Catch if the last bed coord and the no del reference allele length do not match
      #if( (as.numeric(alleleBed$`3UTR`$endPos) + utr5ExtraLen + utr3ExtraLen) != nchar(noDelAlleleStr) ){
      #  message('\nMismatch between BED coordinates to reference allele length\nAllele ID:\t',currentRefAllele,
      #          '\nlast BED coord:\t',alleleBed$`3UTR`$endPos,
      #          '\nallele length:\t',nchar(noDelAlleleStr),
      #          '\niter:\t\t',currentIter)
      #  close(bedCon)
      #  close(fastaCon)
      #  stop()
      #}
      
      # Write allele seq to fasta file
      general.write_fasta(fastaCon, currentRefAllele, noDelAlleleStr)
      
      # Write feature coordinates to bed file
      sapply(alleleBed, function(x){
        general.write_bed(bedCon, x$alleleName, x$startPos, x$endPos, x$featName)
      })
      
    }
    
    # close connections
    close(bedCon)
    close(fastaCon)
    
    currentSample[['refIterVect']] <- iterVect
  }
  
  currentSample[['iterRefDirectory']] <- currentSampleResultsDirectory
  
  return(currentSample)
}


# Write bowtie2 index
sampleObj.iterBowtie2Index <- function(currentSample, bowtie2Build, threads){
  
  if(currentSample$iterRefDirectory == 'failed'){
    return(currentSample)
  }
  
  cat('\n\tBuilding Bowtie2 indexes for:',currentSample$name)
  
  currentSample[['iterIndexPathList']] <- list()
  
  for(currentIter in currentSample$refIterVect){
    iterDir <- normalizePath(file.path(currentSample$iterRefDirectory,currentIter), mustWork=T)
    
    fastaPath <- normalizePath(file.path(iterDir,'alleleReference.fasta'), mustWork=T)
    indexPath <- file.path(iterDir,'alleleReference')
    
    currentSample[['iterIndexPathList']][[currentIter]] <- indexPath
    
    # vcfPath <- file.path(iterDir,paste0(currentSample$name, '.vcf'))
    # if( file.exists(vcfPath) ){
    #   cat('\nFound VCF, skipping bowtie2-build')
    #   next
    # }
    
    ## Creqte a bowtie2 index for the kir_reference.fasta file <- only needed when building a new reference index
    createIndex <- system2(bowtie2Build, c(fastaPath, indexPath,'--quiet',paste('--threads', threads)))
    check.system2_output(createIndex, 'bowtie2 index building failed')
  }
  
  return(currentSample)
}

# Generate VCF from BAM
sampleObj.iterVCFGen <- function(currentSample, samtools, bcftools, threads, deleteBam=T, deleteIndex=T){
  
  if(currentSample$iterRefDirectory == 'failed'){
    return(currentSample)
  }
  
  # bcftools parameters
  mp_m <- '-m 3'
  mp_F <- "-F 0.0002"
  mp_u <- "-u "
  mp_O <- "-O v"
  mp_threads <- paste0('--threads ',threads)
  
  currentSample[['iterVCFPathList']] <- list()
  
  for(currentIter in currentSample$refIterVect){
    
    iterDir <- normalizePath(file.path(currentSample$iterRefDirectory,currentIter), mustWork=T)
    
    vcfPath <- file.path(iterDir,paste0(currentSample$name, '.vcf'))
    
    # Intitialize an output path for the sorted BAM file and VCF file
    currentSample[['iterSortedBamPathList']][[currentIter]] <- file.path(iterDir,paste0(currentSample$name,'.sorted.bam'))
    currentSample[['iterVCFPathList']][[currentIter]] <- vcfPath
    
    if( file.exists(vcfPath) ){
      cat("\nFound VCF, skipping mpileup")
      next
    }
    
    # set samtools output
    sam_sort_output <- paste0('-o ', currentSample[['iterSortedBamPathList']][[currentIter]])
    
    sam_input <- currentSample$iterBamPathList[[currentIter]]
    
    optionsCommand <- c('sort', paste0('-@', threads), sam_sort_output, sam_input)
    
    ## Run the samtools BAM sort command
    cat('\n\n',samtools,optionsCommand)
    output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
    
    check.system2_output(output.sampleAlign, 'samtools bam sort failed')
    
    fastaPath <- normalizePath(file.path(iterDir,'alleleReference.fasta'), mustWork=T)
    bedPath <- normalizePath(file.path(iterDir,'alleleReference.bed'), mustWork=T)
    
    mp_f <- paste0('-f ', fastaPath)
    mp_l <- paste0('-T ', bedPath)
    mp_o <- paste0('-o ', vcfPath)
    #mp_A <- '-A'
    #mp_B <- '-B'
    #mp_m,
    #mp_F,
    optionsCommand <- c('mpileup',mp_m,mp_F,mp_threads,'-d 2000', mp_f, currentSample[['iterSortedBamPathList']][[currentIter]], mp_l,
                        '|', bcftools ,'call','-m','-M',mp_O, mp_o)
    
    ## Run the samtools | bcftools command
    cat('\n\n',bcftools,optionsCommand)
    output.sampleAlign <- system2(bcftools, optionsCommand, stdout=T, stderr=T)
    
    check.system2_output(output.sampleAlign, 'bcftools mpileup | call failed')
    
    if(deleteBam){
      file.remove(currentSample[['iterSortedBamPathList']][[currentIter]])
      file.remove(currentSample$iterBamPathList[[currentIter]])
    }
    
    if(deleteIndex){
      lapply(list.files(iterDir,'.bt2',full.names = T), file.remove)
      lapply(list.files(iterDir,'.fai',full.names = T), file.remove)
    }
  }
  
  return(currentSample)
}

## Samtools sam to bam conversion
iterAlign.sam_to_bam <- function(samtools, samPath, bamPath, threads){
  
  sam_b <- '-b'
  #sam_q <- '-q 10'
  sam_q <- '-q 0'
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads), sam_b, sam_q,
                      samPath, '-o', bamPath)
  
  cat('\n\n',samtools, optionsCommand)
  output.bamConv <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamConv, 'samtools sam to bam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamConv), collapse='\n')
  
  cat('\n\nSuccessfully converted',samPath,'to',bamPath)
  
  return()
}

# bowtie2-align
sampleObj.iterBowtie2Align <- function(currentSample, bowtie2, threads, deleteSam=F, all.align=F, forceRun=F){
  
  if(currentSample$iterRefDirectory == 'failed'){
    return(currentSample)
  }
  
  cat('\n\t-----',currentSample$name,'Bowtie2 alignment -----')
  
  currentSample[['iterSamPathList']] <- list()
  currentSample[['iterBamPathList']] <- list()
  
  for(currentIter in currentSample$refIterVect){
    
    iterDir <- normalizePath(file.path(currentSample$iterRefDirectory,currentIter), mustWork=T)
    
    #-- catch for skipping alingment if VCF file exists
    # vcfPath <- file.path(iterDir,paste0(currentSample$name, '.vcf'))
    # 
    # vcfExists.bool <- file.exists(vcfPath)
    # 
    # if( vcfExists.bool ){
    #   cat('\nFound VCF, skipping bowtie2 alignment')
    #   ## Intitialize an output path for the SAM file
    #   currentSample[['iterSamPathList']][[currentIter]] <- file.path(iterDir,paste0(currentSample$name,'.sam'))
    #   currentSample[['iterBamPathList']][[currentIter]] <- file.path(iterDir,paste0(currentSample$name,'.bam'))
    #   next
    # }
    
    ## Intitialize an output path for the SAM file
    currentSample[['iterSamPathList']][[currentIter]] <- file.path(iterDir,paste0(currentSample$name,'.sam'))
    currentSample[['iterBamPathList']][[currentIter]] <- file.path(iterDir,paste0(currentSample$name,'.bam'))
    
    ### 1. Align KIR extracted reads to haplo-reference
    bt2_p <- paste0("-p", threads)
    bt2_5 <- "-5 3"
    bt2_3 <- "-3 7"
    bt2_i <- "-i S,1,0.5"
    bt2_min_score <- "--score-min L,0,-0.187"
    bt2_I <- "-I 75"
    bt2_X <- "-X 1500"
    bt2_noUnal <- '--no-unal'
    if( all.align ){
      bt2_a <- '-a'
    }
    
    bt2_x <- paste0("-x ", currentSample$iterIndexPathList[[currentIter]])
    bt2_1 <- paste0('-1 ',currentSample$kirfastq1path)
    bt2_2 <- paste0('-2 ',currentSample$kirfastq2path)
    
    bt2_stream <- paste0("-S ", currentSample[['iterSamPathList']][[currentIter]])    
    fastqBase <- file.path(iterDir, currentSample$name)
    bt2_al_conc <- paste0("--al-conc-gz ", fastqBase, "_%.fastq.gz")
    bt2_un <- "--un dump.me"
    
    if( all.align ){
      optionsCommand <- c(bt2_p, bt2_5, bt2_3, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_noUnal, bt2_stream, bt2_al_conc, bt2_un, bt2_a)
    }else{
      optionsCommand <- c(bt2_p, bt2_5, bt2_3, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_noUnal, bt2_stream, bt2_al_conc, bt2_un)
    }
    
    ## Run the bowtie2 alignment command
    cat('\n\n',bowtie2,optionsCommand)
    output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
    
    check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
    
    ## Print the bowtie2 output
    cat('\n',paste0(output.sampleAlign, collapse='\n'))
    
    iterAlign.sam_to_bam(samtools, 
                         currentSample[['iterSamPathList']][[currentIter]], 
                         currentSample[['iterBamPathList']][[currentIter]],
                         threads)
    
    if(deleteSam){
      file.remove(currentSample[['iterSamPathList']][[currentIter]])
    }
    file.remove('dump.me')
    file.remove(paste0(fastqBase,'_1.fastq.gz'))
    file.remove(paste0(fastqBase,'_2.fastq.gz'))
  }
  
  return(currentSample)
}

#skipSampleVect <- c('IND00003','IND00007','IND00026','IND00074','IND00075')

# # Iter align workflow
# sampleList[16:length(sampleList)] <- sapply(sampleList[16:length(sampleList)], function(x){
#   
#   if(x$name %in% skipSampleVect){
#     return(x)
#   }
#   
#   cat('\nLoading ref DF')
#   x <- sampleObj.loadRefDF(x, referenceAlleleDF) # Subset reference allele dataframe by present loci, save to sample object
#   cat('\nWriting reference files')
#   x <- sampleObj.writeRefFastaBed(x, locusRefList, alignmentFileDirectory) # Write fasta reference file for sample object based on refDF
#   x <- sampleObj.iterBowtie2Index(x, bowtie2Build, threads) # Converts fasta file from previous line into a bowtie2 index
#   x <- sampleObj.iterBowtie2Align(x, bowtie2, threads, deleteSam=F) # Align sample to bowtie2 index
#   x <- sampleObj.iterVCFGen(x, samtools, bcftools, threads) # Convert SAM file into VCF
#   
#   x <- allele.iter_alignments_to_snp_dfs(x, locusRefList, referenceAlleleDF, minDP, kirLocusFeatureNameList)
#   x <- allele.combine_iter_snps(x, snpDFList)
#   
#   if( 'KIR2DL2' %in% rownames(x$refAlleleDF) & 'KIR2DL3' %in% rownames(x$refAlleleDF) ){
#     x <- allele.iter_combine_KIR2DL23( x, knownSnpDFList, alleleFileDirectory )
#   }
#   
#   x <- allele.setup_iter_allele_call_df(x)
#   
#   cat('\n\nFinding allele matches for',x$name)
#   for(currentLocus in colnames(x[['iterAlleleCallDF']])){
#     x <- allele.call_allele(x, currentLocus, alleleFileDirectory, knownSnpDFList, alleleDFPathList$iter$newAllelePath, filterLocusConv, 'iter')
#   }
#   
#   cat('\nWriting allele matches to', alleleDFPathList$iter$alleleCallPath )
#   x <- allele.save_call( x, alleleDFPathList$iter$alleleCallPath, 'iter' )
# })


# Filter align workflow
# Create filter align directory structure
' ADDS THESE ATTRIBUTES TO SAMPLE OBJECTS
sample
  $filterRefDirectory
    [RESULTSDIR]/alignmentFiles/[NAME]/filterAlign/
  $filterVCFList
    $[LOCUS]
      [NAME]_[LOCUS].vcf
  $filterBEDList
    $[LOCUS]
      [LOCUS]***.bed
  
'
sampleObj.filterAlign.setup <- function(currentSample, alignmentFileDirectory){
  
  if('failed' %in% c(currentSample$geneContent, currentSample$copyNumber)){
    currentSample[['filterRefDirectory']] <- 'failed'
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  # Create a directory for the current iteration reference files
  currentSampleFilterResultsDirectory <- file.path(alignmentFileDirectory,currentSample$name,'filterAlign')
  if(!file.exists(currentSampleFilterResultsDirectory)){
    dir.create(currentSampleFilterResultsDirectory, recursive=T)
  }
  
  currentSample[['filterRefDirectory']] <- currentSampleFilterResultsDirectory
  currentSample[['filterVCFList']] <- list()
  currentSample[['filterBEDList']] <- list()
  
  return(currentSample)
}

# only used for converting PING1 filter references to PING2
filterAlign.refConversion <- function(locusRefList){
  
  filterRefList <- list('KIR2DL1'='2DL1AClong.fas',
                        'KIR2DL2'='2DL2long.fas',
                        'KIR2DL3'='2DL3long.fas',
                        'KIR3DL3'='3DL3long.fas',
                        'KIR3DL2'='3DL2long.fas',
                        'KIR2DL4'='2DL4FH5.fas',
                        'KIR2DL5'='2DL5B.fas',
                        'KIR2DP1'='2DP1.fas',
                        'KIR2DS4'='2DS4long.fas',
                        'KIR3DL1'='3DL1longCAT.fas',
                        'KIR3DS1'='3DS1longCAT.fas',
                        'KIR2DS5'='2DS5long.fas')
  
  currentLocus <- 'KIR2DL3'
  currentRefFile <- list.files('old_stuff/PING/Resources/caller_resources/Filters/2DL23/',full.names = T, pattern = filterRefList[[currentLocus]])
  
  currentRefAllele <- general.read_fasta(currentRefFile)
  
  cat(names(currentRefAllele),'\n')
  
  featNameVect <- grep('UTR',kirLocusFeatureNameList[[currentLocus]],value=T,invert=T)
  #featNameVect <- c(featNameVect, grep('I',kirLocusFeatureNameList[[currentLocus]],value=T))
  #exonNameVect <- grep('PE',exonNameVect,value=T,invert = T)
  
  sapply(locusRefList[[currentLocus]]$alleleBedList, function(x){
    sapply(x[featNameVect], function(y){
      featSeq <- y[['featSeq']]
      grepl(featSeq, currentRefAllele[1],fixed=T)
    })
  })
  
  alleleMatch <- names(which(sapply(locusRefList[[currentLocus]]$alleleBedList, function(x){
    all(sapply(x[featNameVect], function(y){
      featSeq <- y[['featSeq']]
      grepl(featSeq, currentRefAllele[1],fixed=T)
    }))
  })))
  
  cat(alleleMatch,'\n')
  
  bedMat <- t(sapply(locusRefList[[currentLocus]]$alleleBedList[[alleleMatch]][featNameVect], function(x){
    str_locate(currentRefAllele[[1]], x[['featSeq']])
  }))
  
  bedMat[,1] <- bedMat[,1]-1
  
  bedMat
  
  filterDirectory <- file.path('Resources/genotype_resources/filters/', 'KIR2DL23')
  kirBedPath <- file.path(filterDirectory, '2DL3longtail.bed.new')
  bedCon = file(kirBedPath, "w")
  
  for(i in 1:nrow(bedMat)){
    general.write_bed(bedCon, names(currentRefAllele), bedMat[i,1], bedMat[i,2], rownames(bedMat)[i])
  }
  
  close(bedCon)
  
  
}

general.read_fasta.2 <- function(fastaFile){
  output.fastaList <- list()
  
  con  <- file(fastaFile, open = "r")
  
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    alleleNameBool <- grepl(pattern = fixed('>'), x = oneLine)
    
    if(alleleNameBool){
      alleleName <- gsub('>','',oneLine,fixed=T)
      output.fastaList[[alleleName]] <- ''
    }else{
      output.fastaList[[alleleName]] <- paste0(output.fastaList[[alleleName]], oneLine)
    }
  } 
  close(con)
  
  return(output.fastaList)
}

sampleObj.filterAlign.KIR3DL3 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # Step 1: positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL3/All3DL3'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # Step 2: negative filter
  bt2_score_min <- '--score-min L,0,-0.2'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL3/not3DL3'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # Step 3: final alignment
  bt2_score_min <- '--score-min L,0,-0.35'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL3/3DL3long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_3DL3.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_3DL3.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_3DL3.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_3DL3.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_3DL3.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR3DL3/3DL3long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR3DL3/3DL3long.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR3DL3.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR3DL3']] <- vcfPath
  currentSample[['filterBEDList']][['KIR3DL3']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR3DL2 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL2/All3DL2'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.2'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL2/not3DL2'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL2/3DL2long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_3DL2.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_3DL2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_3DL2.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_3DL2.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_3DL2.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR3DL2/3DL2long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR3DL2/3DL2long.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR3DL2.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR3DL2']] <- vcfPath
  currentSample[['filterBEDList']][['KIR3DL2']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR3DS1 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/All3DL1andS1'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.17'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/not3DL1S1'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_local <- '--local'
  bt2_N <- '-N 1'
  bt2_score_min <- '--score-min L,1,0.6'
  bt2_X <- '-X 750'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_3DS1.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_3DS1.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_local, bt2_N, bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_3DS1.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_3DS1.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_3DS1.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR3DS1.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR3DS1']] <- vcfPath
  currentSample[['filterBEDList']][['KIR3DS1']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR3DL1S1 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/All3DL1andS1'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq'))) # 3DL12in.fastq
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2 3DL1: Negative filter
  bt2_score_min <- '--score-min L,0,-0.2'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/not3DL1S1'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.L1.inter2.fastq'))) # _KIR3DL1S1.fastq
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3 3DL1: Final alignment
  bt2_local <- '--local'
  bt2_N <- '-N 1'
  bt2_score_min <- '--score-min L,1,0.6'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/3DL1longCAT'
  bt2_X <- '-X 750'
  bt2_no_unal <- '--no-unal'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.L1.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.L1.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_3DL1.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_3DL1.fastq')))
  #bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_local, bt2_N, bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam, bt2_no_unal)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2 3DS1: Negative filter
  bt2_score_min <- '--score-min L,0,-0.17'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/not3DL1S1'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.S1.inter2.fastq'))) # _KIR3DL1S1.fastq
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3 3DS1: Final alignment
  bt2_local <- '--local'
  bt2_N <- '-N 1'
  bt2_score_min <- '--score-min L,1,0.6'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT'
  bt2_no_unal <- '--no-unal'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.S1.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.S1.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_3DS1.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_3DS1.fastq')))
  
  optionsCommand <- c(bt2_local, bt2_N, bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_no_unal, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # ------ 3DL1 VCF generation
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_3DL1.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_3DL1.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_3DL1.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR3DL1/3DL1longCAT.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR3DL1/3DL1longCAT.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR3DL1.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  currentSample[['filterVCFList']][['KIR3DL1het']] <- vcfPath
  currentSample[['filterBEDList']][['KIR3DL1het']] <- bedPath
  
  # ------ 3DS1 VCF generation
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_3DS1.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_3DS1.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_3DS1.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR3DS1.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  currentSample[['filterVCFList']][['KIR3DS1het']] <- vcfPath
  currentSample[['filterBEDList']][['KIR3DS1het']] <- bedPath
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  return(currentSample)
}

sampleObj.filterAlign.KIR3DL1 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/All3DL1andS1'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.2'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/not3DL1S1'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_local <- '--local'
  bt2_N <- '-N 1'
  bt2_score_min <- '--score-min L,1,0.6'
  bt2_X <- '-X 750'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR3DL1/3DL1longii'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_3DL1.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_3DL1.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_local, bt2_N, bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_3DL1.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_3DL1.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_3DL1.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR3DL1/3DL1longii.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR3DL1/3DL1longii.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR3DL1.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR3DL1']] <- vcfPath
  currentSample[['filterBEDList']][['KIR3DL1']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DS35 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS35/All2DS35_Gen'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.18'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS35/not2DS35'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.5'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS35/2DS5long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DS35.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DS35.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DS35.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DS35.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DS35.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DS35/2DS5long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DS35/2DS5long.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DS35.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DS35']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DS35']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DS4 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS4/All2DS4'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.1'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS4/not2DS4'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.17'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS4/2DS4long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DS4.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DS4.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DS4.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DS4.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DS4.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DS4/2DS4long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DS4/2DS4noE5.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DS4.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DS4']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DS4']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DS3 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS35/All2DS35_Gen'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.153'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS35/not2DS35'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DS35/2DS3'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DS3.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DS3.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DS3.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DS3.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DS3.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DS35/2DS3.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DS35/2DS3.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DS3.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DS3']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DS3']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DP1 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DP1/All2DP1gen'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.2'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DP1/not2DP1'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DP1/2DP1'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DP1.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DP1.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DP1.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DP1.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DP1.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DP1/2DP1.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DP1/2DP1.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DP1.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DP1']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DP1']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DL5 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL5/All2DL5'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  #bt2_score_min <- '--score-min L,0,-0.1'
  #bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL5/not2DL5'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  #bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  #optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
  #                    bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  #cat('\n\n',bowtie2,optionsCommand)
  #output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  #check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  #message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL5/2DL5B'
  #bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  #bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL5.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL5.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DL5.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DL5.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DL5.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DL5/2DL5B.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DL5/2DL5B.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DL5.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DL5']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DL5']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DL4 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL4/All2DL4gen'
  bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_score_min <- '--score-min L,0,-0.2'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL4/not2DL4'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: Final alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL4/2DL4FH5'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL4.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL4.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DL4.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DL4.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DL4.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DL4/2DL4FH5.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DL4/2DL4FH5.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DL4.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DL4']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DL4']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DL23 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  # step 1: Positive filter
  bt2_threads <- paste0('-p ',threads)
  bt2_5 <- '-5 3'
  bt2_3 <- '-3 7'
  bt2_L <- '-L 20'
  bt2_i <- '-i S,1,0.5'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_I <- '-I 75'
  bt2_X <- '-X 1000'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/All2DL23'
  bt2_sequence_1 <- paste0('-1 ', currentSample$kirfastq1path)
  bt2_sequence_2 <- paste0('-2 ', currentSample$kirfastq2path)
  
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 2: Negative filter
  bt2_3 <- '-3 10'
  bt2_score_min <- '--score-min L,0,-0.135'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/Not2DL23b'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 3: 2DS2 negative filter
  bt2_3 <- '-3 35'
  bt2_score_min <- '--score-min L,0,-0.09'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/All2DS2'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 4: remove L3 tail
  bt2_3 <- '-3 7'
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/2DL3long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.noL3tail.inter4.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter4 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 5: final L2 alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/2DL2long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.noL3tail.inter4.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.noL3tail.inter4.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL2.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL2.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter5 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # Step 6: final L23 alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/2DL3long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL23.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL23.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter6 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # Step 7: remove L2 tail
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/2DL2long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.noS2.inter3.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.noL2tail.inter4.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter7 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # step 8: final L3 alignment
  bt2_score_min <- '--score-min L,0,-0.187'
  bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL23/2DL3long'
  bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.noL2tail.inter4.1.fastq')))
  bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.noL2tail.inter4.2.fastq')))
  bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL3.sam')))
  bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL3.fastq')))
  bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
  
  optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                      bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2,optionsCommand)
  output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'bowtie2 filter iter8 alignment failed')
  message(paste0(output.sampleAlign, collapse='\n'))
  
  # 2DL2 VCF gen ----- 
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DL2.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DL2.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DL2.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DL23/2DL2long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DL23/2DL2longtail.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DL2.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  currentSample[['filterVCFList']][['KIR2DL2']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DL2']] <- bedPath
  
  # 2DL23 VCF gen ----- 
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DL23.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DL23.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DL23.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DL23/2DL3long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DL23/2DL3longshort.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DL23.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  currentSample[['filterVCFList']][['KIR2DL23']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DL23']] <- bedPath
  
  # 2DL3 VCF gen ----- 
  # SAM to BAM conversion
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DL3.sam'))
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DL3.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DL3.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DL23/2DL3long.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DL23/2DL3longtail.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DL3.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DL3']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DL3']] <- bedPath
  return(currentSample)
}

sampleObj.filterAlign.KIR2DL1 <- function(currentSample, bowtie2, samtools, bcftools, threads, removeTemp=T){
  
  if(currentSample$filterRefDirectory == 'failed'){
    currentSample[['filterVCFList']] <- 'failed'
    currentSample[['filterBEDList']] <- 'failed'
    return(currentSample)
  }
  
  tempDir <- file.path('tempFiles')
  if(!file.exists(tempDir)){
    dir.create(tempDir)
  }
  
  KIR2DS1Bool <- F
  if( length(currentSample$copyNumber) > 0 ){
    KIR2DS1Bool <- as.numeric(currentSample$copyNumber$KIR2DS1) == 0
  }
  if( length(currentSample$geneContent) > 0){
    KIR2DS1Bool <- as.numeric(currentSample$geneContent$KIR2DS1) == 0
  }
  
  if( KIR2DS1Bool ){
    
    # step 1: Positive filter
    bt2_threads <- paste0('-p ',threads)
    bt2_5 <- '-5 3'
    bt2_3 <- '-3 7'
    bt2_L <- '-L 20'
    bt2_i <- '-i S,1,0.5'
    bt2_score_min <- '--score-min L,0,-0.187'
    bt2_I <- '-I 75'
    bt2_X <- '-X 1000'
    bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/All2DL1S1gen'
    bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
    bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
    
    bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
    bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
    
    optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                        bt2_sequence_2, bt2_al_conc, bt2_sam)
    
    ## Run the bowtie2 alignment command
    cat('\n\n',bowtie2,optionsCommand)
    output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
    
    check.system2_output(output.sampleAlign, 'bowtie2 filter iter1 alignment failed')
    message(paste0(output.sampleAlign, collapse='\n'))
    
    # step 2: Negative filter
    bt2_score_min <- '--score-min L,0,-0.155'
    bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/not2DL1S1'
    bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
    bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
    bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
    
    optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                        bt2_sequence_2, bt2_un_conc, bt2_sam)
    
    ## Run the bowtie2 alignment command
    cat('\n\n',bowtie2,optionsCommand)
    output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
    
    check.system2_output(output.sampleAlign, 'bowtie2 filter iter2 alignment failed')
    message(paste0(output.sampleAlign, collapse='\n'))
    
    # step 3: Final alignment
    bt2_score_min <- '--score-min L,0,-0.55'
    bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/2DL1AClong'
    bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
    bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
    bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL1.sam')))
    bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL1.fastq')))
    bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
    
    optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                        bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
    
    ## Run the bowtie2 alignment command
    cat('\n\n',bowtie2,optionsCommand)
    output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
    
    check.system2_output(output.sampleAlign, 'bowtie2 filter iter3 alignment failed')
    message(paste0(output.sampleAlign, collapse='\n'))
    
  }else{
    
    if( as.numeric(currentSample$kffHits$`*KIR2DL1*4and7`) >= 10 | as.numeric(currentSample$kffHits$`*KIR2DL1*4710b`) >= 10 ){
      
      # step 1: Positive filter
      bt2_threads <- paste0('-p ',threads)
      bt2_5 <- '-5 3'
      bt2_3 <- '-3 7'
      bt2_L <- '-L 20'
      bt2_i <- '-i S,1,0.5'
      bt2_score_min <- '--score-min L,0,-0.187'
      bt2_I <- '-I 75'
      bt2_X <- '-X 1000'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/All2DL1S1gen'
      bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
      bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
      
      bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
      bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_al_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter4 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
      # step 2: Negative filter
      bt2_score_min <- '--score-min L,0,-0.155'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/not2DL1S1'
      bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
      bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
      bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_un_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter5 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
      # step 3: Negative filter
      bt2_score_min <- '--score-min L,0,-0.153'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/2DS1gen'
      bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
      bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
      bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter3.fastq')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_un_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter6 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
      # step 4: Final alignment
      bt2_score_min <- '--score-min L,0,-0.55'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/2DL1AClong'
      bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter3.1.fastq')))
      bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter3.2.fastq')))
      bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL1.sam')))
      bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL1.fastq')))
      bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter7 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
    }else{
      
      # step 1: Positive filter
      bt2_threads <- paste0('-p ',threads)
      bt2_5 <- '-5 3'
      bt2_3 <- '-3 7'
      bt2_L <- '-L 20'
      bt2_i <- '-i S,1,0.5'
      bt2_score_min <- '--score-min L,0,-0.187'
      bt2_I <- '-I 75'
      bt2_X <- '-X 1000'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/All2DL1S1gen'
      bt2_sequence_1 <- paste0('-1', currentSample$kirfastq1path)
      bt2_sequence_2 <- paste0('-2', currentSample$kirfastq2path)
      
      bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name,'.inter1.fastq')))
      bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'.temp')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_al_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter7 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
      # step 2: Negative filter
      bt2_score_min <- '--score-min L,0,-0.155'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/not2DL1S1'
      bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter1.1.fastq')))
      bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter1.2.fastq')))
      bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter2.fastq')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_un_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter8 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
      # step 3: Negative filter
      bt2_score_min <- '--score-min L,0,-0.18'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/2DS1gen'
      bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter2.1.fastq')))
      bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter2.2.fastq')))
      bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name,'.inter3.fastq')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_un_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter9 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
      # step 4: Final alignment
      bt2_score_min <- '--score-min L,0,-0.55'
      bt2_index <- '-x Resources/genotype_resources/filters/KIR2DL1/2DL1AClong'
      bt2_sequence_1 <- paste0('-1 ', file.path(tempDir, paste0(currentSample$name,'.inter3.1.fastq')))
      bt2_sequence_2 <- paste0('-2 ', file.path(tempDir, paste0(currentSample$name,'.inter3.2.fastq')))
      bt2_sam <- paste0('-S ', file.path(tempDir, paste0(currentSample$name,'_2DL1.sam')))
      bt2_al_conc <- paste0('--al-conc ', file.path(tempDir, paste0(currentSample$name, '_2DL1.fastq')))
      bt2_un_conc <- paste0('--un-conc ', file.path(tempDir, paste0(currentSample$name, '_notmapped.fastq')))
      
      optionsCommand <- c(bt2_threads, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1,
                          bt2_sequence_2, bt2_un_conc, bt2_al_conc, bt2_sam)
      
      ## Run the bowtie2 alignment command
      cat('\n\n',bowtie2,optionsCommand)
      output.sampleAlign <- system2(bowtie2, optionsCommand, stdout=T, stderr=T)
      
      check.system2_output(output.sampleAlign, 'bowtie2 filter iter10 alignment failed')
      message(paste0(output.sampleAlign, collapse='\n'))
      
    }
    
  }
  
  samPath <- file.path(tempDir, paste0(currentSample$name,'_2DL1.sam'))
  
  KIR2DL1.filter_contam_reads_from_sam_file(samPath, currentSample$name)
  
  # SAM to BAM conversion
  bamPath <- file.path(tempDir, paste0(currentSample$name,'_2DL1.bam'))
  
  iterAlign.sam_to_bam(samtools, samPath, bamPath, threads)
  
  sortedBamPath <- file.path(tempDir, paste0(currentSample$name,"_2DL1.sorted.bam"))
  st_out <- paste0("-o ", sortedBamPath)
  st_in <- bamPath
  
  optionsCommand <- c('sort', paste0('-@', threads), st_out, st_in)
  
  ## Run the samtools BAM sort command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools sort failed')
  
  # VCF generation
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/genotype_resources/filters/KIR2DL1/2DL1AClong.fas"
  st_in <- sortedBamPath
  bedPath <- normalizePath('Resources/genotype_resources/filters/KIR2DL1/2DL1AClong.bed',mustWork=T)
  st_l <- paste0('-l ', bedPath)
  st_break <- "|"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  vcfPath <- file.path(currentSample$filterRefDirectory,paste0(currentSample$name,'_KIR2DL1.vcf'))
  bcf_out <- paste0("-o ", vcfPath)
  
  
  optionsCommand <- c('mpileup', st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcftools, 'call', bcf_multi_al, bcf_O, bcf_out)
  
  ## Run the samtools | bcftools command
  cat('\n\n',samtools,optionsCommand)
  output.sampleAlign <- system2(samtools, optionsCommand, stdout=T, stderr=T)
  
  check.system2_output(output.sampleAlign, 'samtools mpileup | bcftools call failed')
  
  if(removeTemp){
    unlink(tempDir, recursive = T)
  }
  
  currentSample[['filterVCFList']][['KIR2DL1']] <- vcfPath
  currentSample[['filterBEDList']][['KIR2DL1']] <- bedPath
  return(currentSample)
}

# # Filter align workflow
# sampleList <- sapply(sampleList, function(x){
#   x <- sampleObj.filterAlign.setup(x, alignmentFileDirectory)
#   
#   if(x[['filterRefDirectory']] == 'failed'){
#     return(x)
#   }
#   
#   locusPresenceList <- sapply(intersect(names(x$geneContent), names(x$copyNumber)), function(locusName){
#     as.numeric(x$copyNumber[[locusName]]) > 0 | as.numeric(x$geneContent[[locusName]]) > 0
#   })
#   
#   x[[ 'samplePresentLocusVect' ]] <- c()
#   
#   x <- sampleObj.filterAlign.KIR3DL3(x, bowtie2, samtools, bcftools, threads)
#   x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR3DL3' )
#   
#   x <- sampleObj.filterAlign.KIR3DL2(x, bowtie2, samtools, bcftools, threads)
#   x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR3DL2' )
#   
#   if( locusPresenceList[['KIR3DL1']] & locusPresenceList[['KIR3DS1']] ){
#     x <- sampleObj.filterAlign.KIR3DL1S1(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR3DL1het', 'KIR3DS1het' )
#   }
#   
#   if( locusPresenceList[['KIR3DL1']] & !locusPresenceList[['KIR3DS1']] ){
#     x <- sampleObj.filterAlign.KIR3DL1(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR3DL1' )
#   }
#   
#   if( locusPresenceList[['KIR3DS1']] & !locusPresenceList[['KIR3DL1']] ){
#     x <- sampleObj.filterAlign.KIR3DS1(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR3DS1' )
#   }
#   
#   # If 2DS5 present regardless of 2DS3
#   if( locusPresenceList[['KIR2DS5']] ){
#     x <- sampleObj.filterAlign.KIR2DS35(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DS35' )
#   }
#   
#   # If 2DS3 present and not 2DS5
#   if( locusPresenceList[['KIR2DS3']] & !locusPresenceList[['KIR2DS5']] ){
#     x <- sampleObj.filterAlign.KIR2DS3(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DS3' )
#   }
#   
#   if( locusPresenceList[['KIR2DS4']] ){
#     x <- sampleObj.filterAlign.KIR2DS4(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DS4' )
#   }
#   
#   if( locusPresenceList[['KIR2DP1']] ){
#     x <- sampleObj.filterAlign.KIR2DP1(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DP1' )
#   }
#   
#   x <- sampleObj.filterAlign.KIR2DL23(x, bowtie2, samtools, bcftools, threads)
#   x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DL23', 'KIR2DL2', 'KIR2DL3' )
#   
#   if( locusPresenceList[['KIR2DL1']] ){
#     x <- sampleObj.filterAlign.KIR2DL1(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DL1' )
#   }
#   
#   if( locusPresenceList[['KIR2DL5']] ){
#     x <- sampleObj.filterAlign.KIR2DL5(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DL5' )
#   }
#   
#   if( locusPresenceList[['KIR2DL4']] ){
#     x <- sampleObj.filterAlign.KIR2DL4(x, bowtie2, samtools, bcftools, threads)
#     x[[ 'samplePresentLocusVect' ]] <- c( x[[ 'samplePresentLocusVect' ]], 'KIR2DL4' )
#   }
#   
#   return(x)
# })