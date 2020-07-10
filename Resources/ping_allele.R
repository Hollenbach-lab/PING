library(data.table)
library(stringr)
library(methods)
library(plotly)
library(R.utils)

rawFastqDirectory
fastqPattern
threads
resultsDirectory
shortNameDelim
minDP <- 10

alignmentFileDirectory

currentSample <- sampleList[[1]]

currentSample$iterVCFPathList
currentSample$iterRefDirectory
currentSample$refIterVect

# ---------- ALLELE CALLING -----------

' TO DO
filter SNP processing
'

referenceAlleleDF <- read.csv('Resources/genotype_resources/master_haplo_iteration_testing_v4.csv',row.names=1,stringsAsFactors = F)

# Create directory for storing files relating to allele calling
alleleFileDirectory <- file.path(resultsDirectory,'alleleFiles')
if(!file.exists(alleleFileDirectory)){
  dir.create(alleleFileDirectory)
}

filterAlleleFileDirectory <- file.path(alleleFileDirectory,'filterFiles')
if(!file.exists(filterAlleleFileDirectory)){
  dir.create(filterAlleleFileDirectory)
}

# ---- Function base -----
# Initialize SNP dataframes for storing aligned SNPs
allele.initialize_SNP_tables <- function(alleleFileDirectory, locusRefList, referenceAlleleDF){
  
  output.snpDFList <- list()
  
  cat('\nInitializing SNP tables...')
  
  # Use the reference allele dataframe to pull out locus and allele names
  locusVect <- rownames(referenceAlleleDF)
  alleleVect <- referenceAlleleDF[locusVect,'iter_1']
  
  alleleList <- list()
  alleleList[locusVect] <- alleleVect
  
  # Create a SNP df for each locus
  for( currentLocus in names(alleleList) ){
    cat('',currentLocus)
    currentAllele <- alleleList[[currentLocus]]
    
    # Use the iter_1 reference allele to pull out all feature lengths
    featLenList <- sapply(locusRefList[[currentLocus]]$alleleBedList[[currentAllele]], function(x){
      featSeq <- x$featSeq
      delIndex <- x$featDelIndex
      seqLength <- nchar(general.del_insert(featSeq, delIndex))
      return(list( seqLength ))
    })
    
    # Pull out exon names
    exonFeatNameVect <- grep('E',names(featLenList),fixed=T,value=T)
    exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
    
    # Pull out other feature names
    otherFeatNameVect <- setdiff(names(featLenList), exonFeatNameVect)
    
    # Separate exon features from other features
    exonLenList <- featLenList[exonFeatNameVect]
    otherLenList <- featLenList[otherFeatNameVect]
    
    # Append the feature label to each position within the feature
    exonLabList <- sapply(names(exonLenList), function(x){
      featEnd <- exonLenList[[x]]
      return(paste0(x,'_',1:as.numeric(featEnd)))
    })
    
    otherLabList <- sapply(names(otherLenList), function(x){
      featEnd <- otherLenList[[x]]
      return(paste0(x,'_',1:as.numeric(featEnd)))
    })
    
    # Unlist the appended labels
    exonLabVect <- unlist(exonLabList, use.names = F)
    otherLabVect <- unlist(otherLabList, use.names = F)
    
    # Initialize CSV paths for saving the dataframe
    exonDFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_exonSNPs.csv'))
    intronDFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_intronSNPs.csv'))
    
    # Create the exon dataframe and save it
    exonSnpDF <- data.frame(matrix(NA,nrow=1,ncol=length(exonLabVect)), check.names = F, stringsAsFactors = F)
    rownames(exonSnpDF) <- 'coord'
    exonSnpDF['coord',] <- exonLabVect
    colnames(exonSnpDF) <- exonLabVect
    write.csv(exonSnpDF, exonDFPath)
    
    # Create the intron (+UTR & PE) dataframe and save it
    otherSnpDF <- data.frame(matrix(NA,nrow=1,ncol=length(otherLabVect)), check.names = F, stringsAsFactors = F)
    rownames(otherSnpDF) <- 'coord'
    otherSnpDF['coord',] <- otherLabVect
    colnames(otherSnpDF) <- otherLabVect
    write.csv(otherSnpDF, intronDFPath)
    
    # Save the CSV paths and dataframes to a return list
    output.snpDFList[[currentLocus]] <- list()
    output.snpDFList[[currentLocus]][['exonSNPs']] <- list('csvPath'=exonDFPath)
    output.snpDFList[[currentLocus]][['intronSNPs']] <- list('csvPath'=intronDFPath)
    output.snpDFList[[currentLocus]][['failure']] <- FALSE
  }
  
  cat('\nFinished. Exon and intron SNP tables located in',alleleFileDirectory)
  return(output.snpDFList)
}

# Returns list of BED features
allele.read_bed_w_del_index <- function(bedFile, locusRefList){
  output.bedList <- list()
  
  con  <- file(bedFile, open = "r")
  
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    
    bedLineVect <- strsplit(oneLine,'\t',fixed=T)[[1]]
    alleleName <- bedLineVect[1]
    featStart <- as.numeric(bedLineVect[2])+1
    featEnd <- as.numeric(bedLineVect[3])
    
    if(featEnd < featStart){
      next
    }
    
    featName <- bedLineVect[4]
    
    locusName <- tstrsplit(alleleName, '*', fixed=T)[[1]]
    
    if( locusName %in% c('KIR2DL5A','KIR2DL5B') ){
      locusName <- 'KIR2DL5'
    }
    
    featObject <- locusRefList[[locusName]]$alleleBedList[[alleleName]][[featName]]
    
    if(!alleleName %in% names(output.bedList)){
      output.bedList[[alleleName]] <- list()
    }
    
    delFeatSeq <- general.del_insert(featObject$featSeq, featObject$featDelIndex)
    
    featLength <- nchar(delFeatSeq)
    
    posVect <- as.character( featStart:featEnd )
    
    featPosVect <- as.character( 1:featLength )
    
    if( length(featObject$featDelIndex) > 0 ){
      featPosVect <- featPosVect[ -featObject$featDelIndex ]
    }
    
    #cat('\n',alleleName,'-',featName,'-',featLength,'-',length(posVect),'-',length(featPosVect))
    
    output.bedList[[alleleName]][ posVect ] <- sapply(featPosVect, function(x){
      paste0(featName,'_',x)
    })
    
    #names(output.bedList[[alleleName]][[featName]]$geneCoords) <- output.bedList[[alleleName]][[featName]]$bedCoords
    
  } 
  close(con)
  
  return(output.bedList)
}

'
Reads VCFs and outputs intron/exon SNP tables for each locus, as well as any INDELs
Adds currentSample$iterSnpDFList
'
allele.iter_alignments_to_snp_dfs <- function(currentSample, locusRefList, referenceAlleleDF, minDP, kirLocusFeatureNameList){
  
  if(currentSample$iterRefDirectory == 'failed'){
    currentSample[['iterSnpDFList']] <- 'failed'
    return(currentSample)
  }
  
  cat('\n\nProcessing alignment files for',currentSample$name,'.\n')
  
  # Initialize SNP dataframes
  sampleSnpDFList <- allele.initialize_SNP_tables(currentSample$iterRefDirectory, locusRefList, currentSample$refAlleleDF)
  cat('\n')
  for(refIter in currentSample$refIterVect){
    cat('\n\t',refIter)
    bedPath <- file.path(currentSample$iterRefDirectory, refIter,'alleleReference.bed')
    bedList <- allele.read_bed_w_del_index(bedPath, locusRefList)
    
    # Pull out the deletion index for the reference allele set
    bedDelIndex <- sapply(names(bedList), function(x){
      
      locusName <- strsplit(x,'*',fixed=T)[[1]][1]
      
      sapply(locusRefList[[locusName]]$alleleBedList[[x]], function(y){
        y$featDelIndex
      })
      
    })
    
    vcfPath <- normalizePath(currentSample$iterVCFPathList[[refIter]], mustWork=T)
    
    cat('\n\t\tReading VCF')
    vcfDT <- general.read_VCF(vcfPath)
    
    cat('\n\t\tProcessing VCF')
    vcfDT <- allele.iter_VCF_process(vcfDT, currentSample, refIter, minDP, bedList)
    
    # Pull out all loci 
    locusVect <- unique(vcfDT$Locus)
    
    cat('\n\t\tWriting SNP table(s)')
    
    for(currentLocus in locusVect){
      cat('',currentLocus)
      
      exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[currentLocus]],fixed=T,value=T)
      exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
      
      otherFeatNameVect <- setdiff(kirLocusFeatureNameList[[currentLocus]], exonFeatNameVect)
      
      exonDT <- vcfDT[ Locus == currentLocus ][ featLab %in% exonFeatNameVect ]
      
      # If there are no exonic positions passing qc checks, then mark this locus as a failure
      if( nrow(exonDT) == 0 ){
        cat('\n----- insufficient depth for this locus, marking as failure -----\n')
        sampleSnpDFList[[currentLocus]]$failure <- TRUE
        next
      }
      
      # Split up SNPs from INDELs
      exonDTList <- general.VCF_sep_INDEL(exonDT)
      exonDT <- exonDTList$nodelDT
      exonIndelDT <- exonDTList$indelDT
      
      exonIndelBoolVect <- !sapply(exonIndelDT$genoVect, function(x) all(x == 1))
      
      # Write exon INDEL table
      if( any(exonIndelBoolVect) ){
        exonIndelDT <- exonIndelDT[exonIndelBoolVect,]
        exonIndelDT$refIter <- refIter
        
        exonIndelRepList <- lapply(1:nrow(exonIndelDT), function(i){
          x <- exonIndelDT[i,]
          allele.formatIndelSnps(x$REF, x$SNP1, x$SNP2, x$featLab, x$featCoord)
        })
        
        indelPath <- file.path(currentSample$iterRefDirectory, paste0(refIter,'_exonINDELs.tsv'))
        
        if(refIter == 'iter_1' | !file.exists(indelPath)){
          appendBool <- FALSE
        }else{
          appendBool <- TRUE
        }
        
        write.table(exonIndelDT[,c('CHROM','POS','ID',
                                   'REF','ALT','QUAL',
                                   'FILTER','INFO','FORMAT',
                                   'GENO','DP','SNP1',
                                   'SNP2','feat','featLab',
                                   'featCoord','Locus','refIter')], 
                    file=indelPath, sep='\t', quote = F, 
                    row.names = F, append=appendBool, col.names = !appendBool)
      }
      
      exonSnpDF <- read.csv(sampleSnpDFList[[currentLocus]]$exonSNPs$csvPath, 
                            stringsAsFactors = F, 
                            check.names = F, 
                            row.names = 1)
      
      currentRow <- paste0(refIter,'_SNP1')
      if(!currentRow %in% rownames(exonSnpDF)){
        exonSnpDF[nrow(exonSnpDF)+1,] <- NA
        rownames(exonSnpDF)[nrow(exonSnpDF)] <- paste0(refIter,'_SNP1')
        exonSnpDF[nrow(exonSnpDF)+1,] <- NA
        rownames(exonSnpDF)[nrow(exonSnpDF)] <- paste0(refIter,'_SNP2')
      }
      
      if( any(exonIndelBoolVect) ){
        lapply(exonIndelRepList, function(x){
        snp1List <- x$SNP1
        snp2List <- x$SNP2
        
        combNames <- unique(names(snp1List), names(snp2List))
        
        # Remove INDEL positions from intron datatable
        exonDT <<- exonDT[ !feat %in% combNames ]
        
        # Write INDEL positions to SNP DF
        exonSnpDF[paste0(refIter,'_SNP1'),names(snp1List)] <<- unlist(snp1List)
        exonSnpDF[paste0(refIter,'_SNP2'),names(snp2List)] <<- unlist(snp2List)
        
        return(NULL)
      })
      }
      
      # Write the exon SNPs to the dataframe
      exonSnpDF[paste0(refIter,'_SNP1'), unlist(exonDT$feat)] <- unlist(exonDT$SNP1)
      exonSnpDF[paste0(refIter,'_SNP2'), unlist(exonDT$feat)] <- unlist(exonDT$SNP2)
      
      snp1Row <- paste0(refIter,'_SNP1')
      snp2Row <- paste0(refIter,'_SNP2')
      
      exonSnpDF <- allele.process_del_index(exonSnpDF, bedDelIndex, unique(exonDT$CHROM), snp1Row, snp2Row)
      
      write.csv(exonSnpDF, sampleSnpDFList[[currentLocus]]$exonSNPs$csvPath)
      
      intronDT <- vcfDT[ Locus == currentLocus ][ featLab %in% otherFeatNameVect ]
      
      # Split up SNPs from INDELs
      intronDTList <- general.VCF_sep_INDEL(intronDT)
      intronDT <- intronDTList$nodelDT
      intronIndelDT <- intronDTList$indelDT
      
      intronIndelBoolVect <- !sapply(intronIndelDT$genoVect, function(x) all(x == 1))
      
      # Write intron INDEL table
      if( any(intronIndelBoolVect) ){
        intronIndelDT <- intronIndelDT[intronIndelBoolVect,]
        intronIndelDT$refIter <- refIter
        
        intronIndelRepList <- lapply(1:nrow(intronIndelDT), function(i){
          x <- intronIndelDT[i,]
          allele.formatIndelSnps(x$REF, x$SNP1, x$SNP2, x$featLab, x$featCoord)
        })
        
        indelPath <- file.path(currentSample$iterRefDirectory, paste0(currentLocus,'_intronINDELs.tsv'))
        
        if(refIter == 'iter_1' | !file.exists(indelPath)){
          appendBool <- FALSE
        }else{
          appendBool <- TRUE
        }
        
        write.table(intronIndelDT[,c('CHROM','POS','ID',
                                     'REF','ALT','QUAL',
                                     'FILTER','INFO','FORMAT',
                                     'GENO','DP','SNP1',
                                     'SNP2','feat','featLab',
                                     'featCoord','Locus','refIter')], 
                    file=indelPath, sep='\t', quote = F, 
                    row.names = F, append=appendBool, col.names = !appendBool)
      }
      
      
      intronSnpDF <- read.csv(sampleSnpDFList[[currentLocus]]$intronSNPs$csvPath, 
                              stringsAsFactors = F, 
                              check.names = F, 
                              row.names = 1)
      
      currentRow <- paste0(refIter,'_SNP1')
      if(!currentRow %in% rownames(intronSnpDF)){
        intronSnpDF[nrow(intronSnpDF)+1,] <- NA
        rownames(intronSnpDF)[nrow(intronSnpDF)] <- paste0(refIter,'_SNP1')
        intronSnpDF[nrow(intronSnpDF)+1,] <- NA
        rownames(intronSnpDF)[nrow(intronSnpDF)] <- paste0(refIter,'_SNP2')
      }
      
      if( any(intronIndelBoolVect) ){
        lapply(intronIndelRepList, function(x){
        snp1List <- x$SNP1
        snp2List <- x$SNP2
        
        combNames <- unique(names(snp1List), names(snp2List))
        
        # Remove INDEL positions from intron datatable
        intronDT <<- intronDT[ !feat %in% combNames ]
        
        # Write INDEL positions to SNP DF
        intronSnpDF[paste0(refIter,'_SNP1'),names(snp1List)] <<- unlist(snp1List)
        intronSnpDF[paste0(refIter,'_SNP2'),names(snp2List)] <<- unlist(snp2List)
        
        return(NULL)
      })
      }
      
      # Write the intron SNPs to the dataframe
      intronSnpDF[paste0(refIter,'_SNP1'), unlist(intronDT$feat)] <- unlist(intronDT$SNP1)
      intronSnpDF[paste0(refIter,'_SNP2'), unlist(intronDT$feat)] <- unlist(intronDT$SNP2)
      
      intronSnpDF <- allele.process_del_index(intronSnpDF, bedDelIndex, unique(exonDT$CHROM), snp1Row, snp2Row)
      
      write.csv(intronSnpDF, sampleSnpDFList[[currentLocus]]$intronSNPs$csvPath)
    }
    cat('\n')
  }
  
  currentSample[['iterSnpDFList']] <- sampleSnpDFList
  return(currentSample)
}

allele.filter.process_del_index <- function(snpDF, bedDelIndex, snp1Row, snp2Row){
  
  for(curPos in bedDelIndex){
    
    if( !curPos %in% colnames(snpDF) ){
      next
    }
    
    spltVect <- strsplit(curPos, '_', fixed=T)[[1]]
    featName <- spltVect[1]
    curDelPos <- as.numeric(spltVect[2])
    
    prevPos <- paste0(featName,'_',curDelPos-1)
    
    # If the previous position passed QC checks and the current position is undefined, then replace NA's with '.'
    if( !any(is.na(snpDF[c(snp1Row, snp2Row),prevPos])) ){
      
      if( all(is.na(snpDF[c(snp1Row, snp2Row), curPos])) ){
        
        snpDF[c(snp1Row, snp2Row), curPos] <- '.'
      }
    }
  }
  
  return(snpDF)
}

allele.process_del_index <- function(snpDF, bedDelIndex, alleleName, snp1Row, snp2Row){
  
  for( featName in names(bedDelIndex[[alleleName]]) ){
    
    curDelIndex <- bedDelIndex[[alleleName]][[featName]]
    
    if( length(curDelIndex) == 0 ){
      next
    }
    
    for(curDelPos  in curDelIndex){
      curPos <- paste0(featName,'_',curDelPos)
      
      if( !curPos %in% colnames(snpDF) ){
        next
      }
      
      prevPos <- paste0(featName,'_',curDelPos-1)
      
      # If the previous position passed QC checks and the current position is undefined, then replace NA's with '.'
      if( !any(is.na(snpDF[c(snp1Row, snp2Row),prevPos])) ){
        
        if( all(is.na(snpDF[c(snp1Row, snp2Row), curPos])) ){
          
          snpDF[c(snp1Row, snp2Row), curPos] <- '.'
        }
      }
    }
    
  }
  
  return(snpDF)
}

allele.formatIndelSnps <- function(REF, SNP1, SNP2, featLab, featCoord){
  output.indelList <- list('SNP1'=list(),'SNP2'=list())
  
  featCoord <- as.numeric(featCoord)
  
  refLen <- nchar(REF)
  snp1Len <- nchar(SNP1)
  snp2Len <- nchar(SNP2)
  
  if(refLen > snp1Len){
    # del scenario
    # SNP1+[.]*(refLen-snp1Len)
    
    indelVect <- c(strsplit(SNP1,'')[[1]], rep('.',(refLen-snp1Len)))
    names(indelVect) <- paste0(featLab,'_',featCoord:(featCoord+refLen-1))
    indelList <- as.list(indelVect)
    
  }else if(refLen < snp1Len){
    # ins scenario
    # SNP1 replaces REF+(snp1Len-refLen)
    
    indelVect <- strsplit(SNP1,'')[[1]]
    names(indelVect) <- paste0(featLab,'_',featCoord:(featCoord+snp1Len-1))
    indelList <- as.list(indelVect)
  }else{
    # snp1 = ref scenario
    # Same logic as ins scenario, should not add extra positions
    
    indelVect <- strsplit(SNP1,'')[[1]]
    names(indelVect) <- paste0(featLab,'_',featCoord:(featCoord+snp1Len-1))
    indelList <- as.list(indelVect)
  }
  
  output.indelList$SNP1 <- indelList
  
  if(refLen > snp2Len){
    # del scenario
    # SNP1+[.]*(refLen-snp2Len)
    
    indelVect <- c(strsplit(SNP2,'')[[1]], rep('.',(refLen-snp2Len)))
    names(indelVect) <- paste0(featLab,'_',featCoord:(featCoord+refLen-1))
    indelList <- as.list(indelVect)
    
  }else if(refLen < snp2Len){
    # ins scenario
    # SNP1 replaces REF+(snp2Len-refLen)
    
    indelVect <- strsplit(SNP2,'')[[1]]
    names(indelVect) <- paste0(featLab,'_',featCoord:(featCoord+snp2Len-1))
    indelList <- as.list(indelVect)
  }else{
    # snp1 = ref scenario
    # Same logic as ins scenario, should not add extra positions
    indelVect <- strsplit(SNP2,'')[[1]]
    names(indelVect) <- paste0(featLab,'_',featCoord:(featCoord+snp2Len-1))
    indelList <- as.list(indelVect)
  }
  
  output.indelList$SNP2 <- indelList
  
  return(output.indelList)
}

# Filters VCF data table by minimum depth, adds SNP columns and feature/position labels
allele.iter_VCF_process <- function(vcfDT, currentSample, refIter, minDP, bedList){
  # Remove low quality positions
  vcfDT <- vcfDT[GENO != './.']
  
  # Pull out Depth values
  dpScore <- unlist(lapply(strsplit(vcfDT$INFO,fixed(';')), function(x){
    grep('DP=',x,value=T)
  }))
  
  # Check that only depth values were pulled
  dpCheckBool <- all(grepl('DP=',dpScore))
  
  # Stop program if a non depth value was pulled
  if(!dpCheckBool){
    message('\nNon DP element found.\t',currentSample$name,'\t',refIter)
    stop()
  }
  
  # Set a depth column
  vcfDT$DP <- as.numeric(tstrsplit(dpScore, fixed('DP='))[[2]])
  
  # Split up the geno call into a vector
  vcfDT$genoVect <- lapply(strsplit(tstrsplit(vcfDT$GENO, fixed(':'))[[1]], fixed('/')), function(x){
    return(as.numeric(x)+1)
  })
  
  # Filter out low depth positions
  vcfDT <- vcfDT[DP >= minDP]
  
  # Set important features as new data table columns
  vcfDT[, c('SNP1','SNP2') := general.vcfDT_set_SNPs( REF, ALT, genoVect ), by = 1:nrow(vcfDT) ]
  vcfDT[, feat:= bedList[[CHROM]][[POS]] , by = 1:nrow(vcfDT) ]
  vcfDT$featLab <- tstrsplit(vcfDT$feat, '_', fixed=T)[[1]]
  vcfDT$featCoord <- tstrsplit(vcfDT$feat, '_', fixed=T)[[2]]
  vcfDT$Locus <- tstrsplit(vcfDT$CHROM, '*', fixed=T)[[1]]
  
  vcfDT$Locus[vcfDT$Locus %in% c('KIR2DL5A','KIR2DL5B')] <- 'KIR2DL5'
  
  return(vcfDT)
}

allele.set_DF_snps <- function(snpDF, snpVect, coord, iterLab){
  coordCol <- which(snpDF['coord',] == coord)
  
  if(length(coordCol) != 1){
    stop('coordinate not found')
  }
  
  snp1Row <- paste0(iterLab,'_SNP1')
  snp2Row <- paste0(iterLab,'_SNP2')
  
  snpDF[snp1Row,coordCol] <- snpVect[1]
  snpDF[snp2Row,coordCol] <- snpVect[2]
}

# 
general.vcfDT_set_SNPs <- function(snp1, snp2, genoVect){
  
  if( grepl(fixed(','),snp2) ){
    rawSnpVect <- unlist(strsplit(paste(snp1, snp2, sep=','), fixed(',')))
  }else{
    rawSnpVect <- c(snp1, snp2)
  }
  
  snpVect <- rawSnpVect[genoVect[[1]]]

  return(list(snpVect[1],snpVect[2]))
}

# Inserts deletion characters into sequence string using a deletion index vector
general.del_insert <- function(noDelStr,delIndexVect){
  
  for(currentDelIndex in delIndexVect){
    
    noDelChrVect <- strsplit(noDelStr, '')[[1]]
    
    if( currentDelIndex > 1 ){
      firstHalf <- noDelChrVect[1:(currentDelIndex-1)]
    }else{
      firstHalf <- ''
    }
    
    if( currentDelIndex < (length(noDelChrVect)+1) ){
      secondHalf <- noDelChrVect[(currentDelIndex):length(noDelChrVect)]
    }else{
      secondHalf <- ''
    }
    
    noDelStr <- paste0(paste0(firstHalf, collapse=''), '.', paste0(secondHalf, collapse=''))
    
  }
  
  return(noDelStr)
}

# Returns data table
general.read_VCF <- function(vcfFile){
  vcfDT <- fread(vcfFile)
  
  colnames(vcfDT) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GENO')
  
  return(vcfDT)
}

# Returns list with 'nodelDT' and 'indelDT'
general.VCF_sep_INDEL <- function(vcfDT){
  indelIndex <- which(grepl('INDEL', tstrsplit(vcfDT$INFO, fixed(';'))[[1]]))
  
  if(length(indelIndex) > 0){
    nodelDT <- vcfDT[-indelIndex,]
  }else{
    nodelDT <- vcfDT
  }
  
  indelDT <- vcfDT[indelIndex,]
  return(list('indelDT'=indelDT,'nodelDT'=nodelDT))
}

allele.combine_iter_snps <- function(currentSample, snpDFList){
  
  if(currentSample$iterRefDirectory == 'failed'){
    return(currentSample)
  }
  
  cat('\nCombining iter snps for',currentSample$name)
  
  locusVect <- names(currentSample$iterSnpDFList)
  
  for(currentLocus in locusVect){
    
    if( currentSample$iterSnpDFList[[currentLocus]]$failure == TRUE ){
      cat('\n----- skipping locus -----\n')
      next
    }
      
    cat('',currentLocus)
    masterExonSnpDFPath <- snpDFList[[currentLocus]]$exonSNPs$csvPath
    masterIntronSnpDFPath <- snpDFList[[currentLocus]]$intronSNPs$csvPath
    
    sampleExonSnpDFPath <- currentSample$iterSnpDFList[[currentLocus]]$exonSNPs$csvPath
    sampleIntronSnpDFPath <- currentSample$iterSnpDFList[[currentLocus]]$intronSNPs$csvPath
    
    # Intron Processing
    sampleIntronSnpDF <- read.csv(sampleIntronSnpDFPath,
                                stringsAsFactors = F,
                                check.names = F,
                                row.names = 1)
    
    masterIntronSnpDF <- read.csv(masterIntronSnpDFPath,
                                stringsAsFactors = F,
                                check.names = F,
                                row.names = 1)
    
    snp1Row <- paste0(currentSample$name,'_SNP1')
    snp2Row <- paste0(currentSample$name,'_SNP2')
    if(!snp1Row %in% rownames(masterIntronSnpDF)){
      masterIntronSnpDF[nrow(masterIntronSnpDF)+1,] <- NA
      rownames(masterIntronSnpDF)[nrow(masterIntronSnpDF)] <- paste0(currentSample$name,'_SNP1')
      masterIntronSnpDF[nrow(masterIntronSnpDF)+1,] <- NA
      rownames(masterIntronSnpDF)[nrow(masterIntronSnpDF)] <- paste0(currentSample$name,'_SNP2')
    }
    
    for(currentCol in colnames(masterIntronSnpDF)){
      snpVect <- unique(sampleIntronSnpDF[2:nrow(sampleIntronSnpDF),currentCol])
      
      snpVect <- snpVect[!is.na(snpVect)]
      
      if(length(snpVect) > 2){
        cat('\n--- Triple SNP found ---',currentCol)
        masterIntronSnpDF[c(snp1Row, snp2Row),currentCol] <- '3SNP'
      }else if( length(snpVect) == 0 ){
        masterIntronSnpDF[c(snp1Row, snp2Row),currentCol] <- NA
      }else{
        masterIntronSnpDF[c(snp1Row, snp2Row),currentCol] <- snpVect
      }
    }
    
    write.csv(masterIntronSnpDF, masterIntronSnpDFPath)
    
    # Exon processing
    sampleExonSnpDF <- read.csv(sampleExonSnpDFPath,
                                stringsAsFactors = F,
                                check.names = F,
                                row.names = 1)
    
    masterExonSnpDF <- read.csv(masterExonSnpDFPath,
                                stringsAsFactors = F,
                                check.names = F,
                                row.names = 1)
    
    snp1Row <- paste0(currentSample$name,'_SNP1')
    snp2Row <- paste0(currentSample$name,'_SNP2')
    if(!snp1Row %in% rownames(masterExonSnpDF)){
      masterExonSnpDF[nrow(masterExonSnpDF)+1,] <- NA
      rownames(masterExonSnpDF)[nrow(masterExonSnpDF)] <- paste0(currentSample$name,'_SNP1')
      masterExonSnpDF[nrow(masterExonSnpDF)+1,] <- NA
      rownames(masterExonSnpDF)[nrow(masterExonSnpDF)] <- paste0(currentSample$name,'_SNP2')
    }
    
    for(currentCol in colnames(masterExonSnpDF)){
      snpVect <- unique(sampleExonSnpDF[2:nrow(sampleExonSnpDF),currentCol])
      
      snpVect <- snpVect[!is.na(snpVect)]
      
      if(length(snpVect) > 2){
        cat('\n--- Triple SNP found ---',currentCol)
        masterExonSnpDF[c(snp1Row, snp2Row),currentCol] <- '3SNP'
      }else if( length(snpVect) == 0 ){
        masterExonSnpDF[c(snp1Row, snp2Row),currentCol] <- NA
      }else{
        masterExonSnpDF[c(snp1Row, snp2Row),currentCol] <- snpVect
      }
    }
    
    write.csv(masterExonSnpDF, masterExonSnpDFPath)
  }
  return(currentSample)
}

# Isolate allele differentiating SNP positions
## Use the locusRefList to generate known allele SNP files
allele.create_allele_resources <- function(locusRefList, alleleFileDirectory){
  
  cat('\nWriting allele SNP resource files.')
  output.alleleSnpDFList <- list()
  
  for( currentLocus in names(locusRefList) ){
    cat('',currentLocus)
    
    alleleVect <- names(locusRefList[[currentLocus]]$alleleBedList)
    
    # initialize DF for storing SNPs
    locusSnpDF <- data.frame(matrix(NA,nrow=length(alleleVect),ncol=0),stringsAsFactors = F)
    row.names(locusSnpDF) <- alleleVect
    
    featVect <- names(locusRefList[[currentLocus]]$alleleBedList[[1]])
    
    #currentFeat <- featVect[1]
    for( currentFeat in featVect ){
      
      featSnpMat <- sapply(locusRefList[[currentLocus]]$alleleBedList, function(x){
        x[[currentFeat]]$snpVect
      })
      
      # Transpose SNP matrix
      featSnpMat <- t(featSnpMat)
      
      # Append SNP matrix to the locus SNP df
      locusSnpDF <- cbind(locusSnpDF, featSnpMat)
    }
    
    output.alleleSnpDFList[[currentLocus]] <- list('csvPath' = file.path(alleleFileDirectory, paste0(currentLocus,'_alleleSNPs.csv')),
                                                   'snpDF' = locusSnpDF)
      
      
    write.csv(x = locusSnpDF, file = output.alleleSnpDFList[[currentLocus]]$csvPath )
  }
  
  cat('\nFinished. Allele SNPs written to',alleleFileDirectory,'as [LOCUS]_alleleSNPs.csv')
  return(output.alleleSnpDFList)
}

'
Below function should be modified to return sample object marking a new snp/allele
'
# Writes any newly found snps to alleleFileDirectory/newSNPs/[locus]_newSNPs.csv
allele.save_new_snps <- function(newSnpLocus, newSnpMat, knownSnpDF, alleleFileDirectory, sampleSnpDF, nonAdCols){
  
  # Create directory for storing files relating to allele calling
  newSnpDirectory <- file.path(alleleFileDirectory,'newSNPs')
  
  if(!file.exists(newSnpDirectory)){
    dir.create(newSnpDirectory)
  }
  
  newSnpPath <- file.path(newSnpDirectory,paste0(newSnpLocus,'_newSNPs.tsv'))
  
  # Pull out new Snp indices
  newSnpIndMat <- which( newSnpMat, arr.ind=T )
  
  # Save each new SNP
  for( newSnp in rownames(newSnpIndMat) ){
  
    snpRow <- newSnpIndMat[newSnp,1]
    snpPos <- nonAdCols[newSnpIndMat[newSnp,2]]
    snpNuc <- sampleSnpDF[newSnp,snpPos]
    
    cat('\t', snpPos)
    
    # Pull out known nucleotides for this position
    knownNuc <- paste0(unique(knownSnpDF[,snpPos]), collapse='')
    
    if( !file.exists(newSnpPath) ){
      appendBool <- F
    }else{
      appendBool <- T
    }
    
    newSnpDF <- data.frame(matrix(NA,nrow=1,ncol=3), row.names=newSnp)
    colnames(newSnpDF) <- c('Pos','Nuc','Known')
    
    newSnpDF[newSnp,c('Pos','Nuc','Known')] <- c(snpPos,snpNuc,knownNuc)
    
    write.table( newSnpDF, file=newSnpPath, 
                 sep='\t', quote=F, row.names=T, 
                 append=appendBool, col.names=!appendBool)
  }
  
  cat('\nFinished saving',newSnpLocus,'new SNP(s)')
  return(length(rownames(newSnpIndMat)))
}


# Allele calling
allele.setup_allele_call_df <- function(currentSample){
  currentSampleLociVect <- rownames(currentSample$refAlleleDF)
  currentSample[['iterAlleleCallDF']] <- data.frame(matrix(NA,nrow=3,ncol=length(currentSampleLociVect)))
  colnames(currentSample[['iterAlleleCallDF']]) <- currentSampleLociVect
  rownames(currentSample[['iterAlleleCallDF']]) <- c('allele_call','mismatch_score','new_snps')
  
  return(currentSample)
}

allele.call_allele <- function(currentSample, currentLocus, alleleFileDirectory, knownSnpDFList, newAlleleDFPath, workflow){
  
  if(currentSample$iterRefDirectory == 'failed'){
    currentSample[['iterAlleleCallDF']] <- 'failed'
    return(currentSample)
  }
  
  
  if( currentSample$iterSnpDFList[[currentLocus]]$failure == TRUE ){
    cat('\n----- skipping locus -----\n')
    currentSample[['iterAlleleCallDF']][,currentLocus] <- 'failed'
    return(currentSample)
  }
  
  
  cat('\n\nMatching alleles for',currentLocus)
  
  # load up the aligned SNP df for the current locus
  locusExonSnpDF <- read.csv(file.path(alleleFileDirectory,paste0(currentLocus,'_exonSNPs.csv')), check.names=F, stringsAsFactors = F, row.names = 1)
  
  # Define rownames for the current sample
  sampleRows <- paste0(currentSample$name,c('_SNP1','_SNP2'))
  
  # Pull out the aligned SNPs for the current sample
  sampleSnpDF <- locusExonSnpDF[sampleRows,]
  
  # Format sample DF and transform into data table
  sampleSnpDT <- as.data.table(sampleSnpDF)
  sampleSnpDT$name <- rownames(sampleSnpDF)
  setkey(sampleSnpDT, name)
  
  # Pull out the aligned positions that passed QC checks
  definedCols <- colnames(sampleSnpDF)[apply(sampleSnpDF, 2, function(x){ all(is.nuc(x)) })]
  
  # Subset the known SNP df by the aligned SNP positions
  knownSnpDF <- knownSnpDFList[[currentLocus]]$snpDF[,definedCols]
  
  # Generate allele differentiating (AD) SNP df
  adSnpDF <- make_unique_pos_frame(knownSnpDF)
  
  # Remove duplicate alleles
  adSnpDF <- unique(adSnpDF)
  
  adSnpDT <- as.data.table(adSnpDF)
  adSnpDT$alleleName <- rownames(adSnpDF)
  setkey(adSnpDT, alleleName)
  
  # Pull out AD positions
  adCols <- colnames(adSnpDF)
  
  # Identify non variable positions (based on known alleles)
  nonAdCols <- setdiff(colnames(knownSnpDF),adCols)
  
  # Identify new snps
  newSnpMat <- sapply(nonAdCols, function(x){ 
    apply(sampleSnpDF[,x,drop=F],1,function(y){
      !any(y == knownSnpDF[,x,drop=F])
    })
  })
  
  # Process new snps
  if( any(newSnpMat) ){
    cat('\n----- Found new SNPs, saving to file -----')
    newSnpInt <- allele.save_new_snps(currentLocus, newSnpMat, knownSnpDF, alleleFileDirectory, sampleSnpDF, nonAdCols) # Should add return to mark this sample as having a new allele/SNPs
  }else{
    newSnpInt <- 0
  }
  
  # Identify het positions (should be length 0 for homozygous)
  namedHetVect <- which( apply(sampleSnpDF, 2, function(x){ num_unique_nuc(x) }) == 2 )
  
  # Process het SNPs
  hetPosVect <- names(namedHetVect)
  hetPosVect <- intersect(hetPosVect, adCols) # This removes new snp hets which mess up allele calling
  homPosVect <- setdiff(adCols, hetPosVect)
  
  # Set if het allele calling should be enabled
  hetBool <- length(hetPosVect) > 0
  
  # Narrow down adSnpDT
  homScoreList <- sapply(rownames(adSnpDF), function(x){
    sum(sapply(homPosVect, function(y){
      (adSnpDT[x, ..y] != sampleSnpDT[1, ..y])*1
    }))
  })
  
  if( !hetBool ){
    bestScoreInt <- min(homScoreList)
    bestMatchIndex <- which( homScoreList == bestScoreInt )
    bestMatchAlleleVect <- rownames(adSnpDF)[bestMatchIndex]
    bestMatchAlleleMat <- combinations(length(bestMatchAlleleVect), 2, bestMatchAlleleVect, repeats.allowed = T)
    allMatchingAlleleVect <- apply(bestMatchAlleleMat,1,paste0,collapse='+')
  }else{
    homAlleleVect <- rownames(adSnpDF)[homScoreList < 2] # Cut any alleles with 2 or more mismatches
    
    # Generate all possible allele pairings
    possAllelePairMat <- combinations(length(homAlleleVect),2,homAlleleVect,repeats.allowed = T)
    
    # Set up possible allele data table
    possAlleleDT <- as.data.table(possAllelePairMat)
    possAlleleDT$allelePair <- apply(possAllelePairMat,1,paste0,collapse='+')
    setkey(possAlleleDT, allelePair)
    colnames(possAlleleDT)[1:2] <- c('allele1','allele2')
    
    # initialize column for storing distance scores
    possAlleleDT$distance <- 0
    
    # First score homozygous positions for allele pairs (this lowers the computational load for the het position scoring)
    possAlleleDT[, 'distance' := allele.add_hom_score( allele1, allele2, distance, homScoreList ), by=allelePair]
    possAlleleDT <- possAlleleDT[possAlleleDT$distance < 2,] # Remove all allele pairings that have more than 1 mismatch
    
    # Score het positions for allele pairs (this can take awhile if there are many pairs)
    cat('\nScoring',nrow(possAlleleDT),'allele pairings...')
    possAlleleDT[, 'distance' := allele.pair_score_calc( allele1, allele2, distance, adSnpDT, sampleSnpDT, hetPosVect ), by=allelePair ]
    
    bestScoreInt <- min(possAlleleDT$distance)
    bestMatchIndex <- which(possAlleleDT$distance == bestScoreInt)
    
    allMatchingAlleleVect <- possAlleleDT$allelePair[bestMatchIndex]
  }
  
  allMatchingAlleleStr <- paste0(allMatchingAlleleVect, collapse=' ')
  
  bestScoreInt <- bestScoreInt+newSnpInt
  
  # New allele catch
  if( bestScoreInt > 0 ){
    allele.save_new_allele( currentSample, currentLocus, knownSnpDF, sampleSnpDF, allMatchingAlleleStr, newAlleleDFPath )
  }
  
  # Print out the best score and cooresponding allele matches
  cat('\nBest score:',bestScoreInt)
  cat('\nAllele match(es):',allMatchingAlleleStr)
  
  if( workflow == 'iter' ){
    currentSample[['iterAlleleCallDF']]['allele_call',currentLocus] <- allMatchingAlleleStr
    currentSample[['iterAlleleCallDF']]['mismatch_score',currentLocus] <- bestScoreInt
    currentSample[['iterAlleleCallDF']]['new_snps',currentLocus] <- newSnpInt
  }else if( workflow == 'filter' ){
    currentSample[['filterAlleleCallDF']]['allele_call',currentLocus] <- allMatchingAlleleStr
    currentSample[['filterAlleleCallDF']]['mismatch_score',currentLocus] <- bestScoreInt
    currentSample[['filterAlleleCallDF']]['new_snps',currentLocus] <- newSnpInt
  }
  
  return(currentSample)
}

allele.save_new_allele <- function( currentSample, currentLocus, knownSnpDF, sampleSnpDF, allMatchingAlleleStr, newAlleleDFPath ){
  
  cat('\nSaving new allele call.')
  
  newDF <- read.table(newAlleleDFPath, check.names=F, stringsAsFactors = F, sep=',', header = T)
  
  mismatchPosVect <- sapply(tstrsplit(allMatchingAlleleStr, ' '), function(x){
    checkCols <- colnames(sampleSnpDF)
    
    alleleVect <- strsplit(x,'+',fixed=T)[[1]]
    
    scoreList <- sapply(checkCols, function(curPos){
      
      con1 <- !all(sampleSnpDF[, curPos][1] == knownSnpDF[ alleleVect, curPos][1])
      con2 <- !all(sampleSnpDF[, curPos][2] == knownSnpDF[ alleleVect, curPos][2])
      
      con3 <- !all(sampleSnpDF[, curPos][1] == knownSnpDF[ alleleVect, curPos][2])
      con4 <- !all(sampleSnpDF[, curPos][2] == knownSnpDF[ alleleVect, curPos][1])
      
      return( min((con1*1+con2*1),(con3*1+con4*1)) )
    })
    
    mismatchIndex <- which( scoreList > 0 )
    
    if(length(mismatchIndex) == 0){
      cat('\n\tIn new allele calling but no mismatches found.')
      stop()
    }
    
    return(names(mismatchIndex))
  })
  
  mismatchPosVect <- unique(unlist(mismatchPosVect))
  mismatchPosStr <- paste0(mismatchPosVect, collapse=' ')
  
  fun.output <- c(currentSample$name, currentLocus, mismatchPosStr, allMatchingAlleleStr, snpScore>0)
  
  newDF[ (nrow(newDF)+1) , c('sampleName','locusName','mimatchPos','bestAlleleMatch','newSnpFound') ] <- fun.output
  
  write.table( newDF, file=newAlleleDFPath, sep=',', quote=F, row.names=T, col.names=T )
  
}

allele.add_hom_score <- function( allele1, allele2, distance, homScoreList ){
  return( homScoreList[[allele1]] + homScoreList[[allele2]] + distance )
}

allele.pair_score_calc <- function( allele1, allele2, distance, adSnpDT, sampleSnpDT, adCols ){

  scoreList <- sapply(adCols, function(curPos){
    
    con1 <- !all(sampleSnpDT[, ..curPos][1,] == adSnpDT[ c(allele1, allele2), ..curPos][1,])
    con2 <- !all(sampleSnpDT[, ..curPos][2,] == adSnpDT[ c(allele1, allele2), ..curPos][2,])
    
    con3 <- !all(sampleSnpDT[, ..curPos][1,] == adSnpDT[ c(allele1, allele2), ..curPos][2,])
    con4 <- !all(sampleSnpDT[, ..curPos][2,] == adSnpDT[ c(allele1, allele2), ..curPos][1,])
    
    return( min((con1*1+con2*1),(con3*1+con4*1)) )
  })
  
  distanceInt <- sum(scoreList)

  return( distanceInt + distance )
}

allele.save_call <- function( currentSample, acDFPath, newDFPath ){

  if(currentSample$iterRefDirectory == 'failed'){
    currentSample[['iterAlleleCallDF']] <- 'failed'
    return(currentSample)
  }
  
  acDF <- read.table(acDFPath, check.names = F, stringsAsFactors = F, sep=',')
  
  for(currentLocus in colnames(currentSample[['iterAlleleCallDF']])){
    alleleCall <- currentSample[['iterAlleleCallDF']]['allele_call',currentLocus]
    callScore <- currentSample[['iterAlleleCallDF']]['mismatch_score',currentLocus]
    snpScore <- currentSample[['iterAlleleCallDF']]['new_snps',currentLocus]
    
    if(callScore > 0){
      alleleCall <- 'unresolved_mismatch'
    }
    
    if( snpScore > 0 ){
      alleleCall <- 'unresolved_newSnp(s)'
    }
    
    acDF[currentSample$name,currentLocus] <- alleleCall
    
  }
  
  write.table( acDF, file=acDFPath, sep=',', quote=F, row.names=T, col.names=T )
  
  return(currentSample)
}

allele.setup_results_df <- function( locusRefList, resultsDirectory, sampleList, workflow ){
  locusVect <- names(locusRefList)
  
  acDF <- data.frame(matrix(NA,nrow=length(sampleList),ncol=length(locusVect)), stringsAsFactors = F)
  colnames(acDF) <- locusVect
  rownames(acDF) <- unlist(lapply(sampleList, function(x) x$name))
  
  alleleCallPath <- file.path(resultsDirectory, paste0(workflow,'AlleleCalls.csv'))
  write.table( acDF, file=alleleCallPath, sep=',', quote=F, row.names=T )
  
  newDF <- data.frame(matrix(NA,nrow=0,ncol=5), stringsAsFactors = F)
  colnames(newDF) <- c('sampleName','locusName','mimatchPos','bestAlleleMatch','newSnpFound')
  
  newCallPath <- file.path(resultsDirectory, paste0(workflow,'NewAlleles.csv'))
  write.table( newDF, file= newCallPath, sep=',', quote=F, row.names=T)
  
  return( list('alleleCallPath'=alleleCallPath, 'newAllelePath'=newCallPath) )
}


# ----- Run workflow -----
# Initialze data frames for consolidating SNPs across all samples for each gene
snpDFList <- allele.initialize_SNP_tables(alleleFileDirectory, locusRefList, referenceAlleleDF)

# Generate known SNP df's for allele calling
knownSnpDFList <- allele.create_allele_resources(locusRefList, alleleFileDirectory)

alleleDFPathList <- allele.setup_results_df( locusRefList, resultsDirectory, sampleList, 'iter')

filterAlleleDFPathList <- allele.setup_results_df( locusRefList, resultsDirectory, sampleList, 'filter' )

# Iter SNP consolidation workflow + allele calling
sampleList[1:10] <- sapply(sampleList[1:10], function(x){
  x <- allele.iter_alignments_to_snp_dfs(x, locusRefList, referenceAlleleDF, minDP, kirLocusFeatureNameList)
  x <- allele.combine_iter_snps(x, snpDFList)
  x <- allele.setup_allele_call_df(x)
  
  cat('\n\nFinding allele matches for',x$name)
  for(currentLocus in rownames(x$refAlleleDF)){
    x <- allele.call_allele(x, currentLocus, alleleFileDirectory, knownSnpDFList, alleleDFPathList$newAllelePath, 'iter')
  }
  cat('\nWriting allele matches to', alleleDFPathList$alleleCallPath )
  x <- allele.save_call( x, alleleDFPathList$alleleCallPath )
})


filterSnpDFList <- allele.initialize_filter_SNP_tables(filterAlleleFileDirectory, locusRefList, names(filterLocusConv), filterLocusConv)

# Filter SNP processing workflow + allele calling
sampleList[4] <- sapply(sampleList[4], function(x){
  cat('\n\nProcessing filtration alignments for',x$name)
  x <- allele.filter_alignments_to_snp_dfs(x, locusRefList, minDP, kirLocusFeatureNameList, filterRefFastaList, filterLocusConv, filterSnpDFList, knownSnpDFList)
})

# Processes VCF data and outputs SNP dataframes to pass to allele calling
allele.filter_alignments_to_snp_dfs <- function(currentSample, locusRefList, minDP, kirLocusFeatureNameList, filterRefFastaList, filterLocusConv, filterSnpDFList, knownSnpDFList){
  
  if('failed' %in% c(currentSample$geneContent, currentSample$copyNumber)){
    return(currentSample)
  }
  
  cat('\nConverting VCF files to SNP dataframes...')
  
  for( currentLocus in names(currentSample$filterVCFList) ){
    
    cat('\n',currentLocus)
    
    vcfPath <- currentSample$filterVCFList[[currentLocus]]
    bedPath <- currentSample$filterBEDList[[currentLocus]]
    
    cat('\n\tReading BED')
    bedList <- allele.convert_filter_bed(bedPath, currentLocus, filterRefFastaList, locusRefList, filterLocusConv)
    
    realLocus <- filterLocusConv[[currentLocus]]
    
    # Pull out the deletion index for the reference locus
    bedDelIndex <- names( which( apply(knownSnpDFList[[realLocus]]$snpDF, 2, function(x){ '.' %in% x }) ) )
    
    exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[realLocus]],fixed=T,value=T)
    exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
    
    otherFeatNameVect <- setdiff(kirLocusFeatureNameList[[realLocus]], exonFeatNameVect)
    
    cat('\n\tReading VCF')
    vcfDT <- general.read_VCF(vcfPath)
    
    vcfDT <- allele.filter_VCF_process(vcfDT, currentSample, minDP, bedList, realLocus)
    
    # ----- EXON processing -----
    cat('\n\tProcessing exons')
    exonDT <- vcfDT[ Locus == realLocus ][ featLab %in% exonFeatNameVect ]
    
    # Split up SNPs from INDELs
    exonDTList <- general.VCF_sep_INDEL(exonDT)
    exonDT <- exonDTList$nodelDT
    exonIndelDT <- exonDTList$indelDT
    
    exonIndelBoolVect <- !sapply(exonIndelDT$genoVect, function(x) all(x == 1))
    
    # Write exon INDEL table
    if( any(exonIndelBoolVect) ){
      exonIndelDT <- exonIndelDT[exonIndelBoolVect,]
      
      exonIndelList <- lapply(1:nrow(exonIndelDT), function(i){
        x <- exonIndelDT[i,]
        allele.formatIndelSnps(x$REF, x$SNP1, x$SNP2, x$featLab, x$featCoord)
      })
      
      indelPath <- file.path(currentSample$filterRefDirectory, paste0(currentLocus,'_exonINDELs.tsv'))
      
      appendBool <- FALSE
      
      write.table(exonIndelDT[,c('CHROM','POS','ID',
                                   'REF','ALT','QUAL',
                                   'FILTER','INFO','FORMAT',
                                   'GENO','DP','SNP1',
                                   'SNP2','feat','featLab',
                                   'featCoord','Locus')], 
                  file=indelPath, sep='\t', quote = F, 
                  row.names = F, append=appendBool, col.names = !appendBool)
    }
    
    
    exonSnpDF <- read.csv(filterSnpDFList[[currentLocus]]$exonSNPs$csvPath, 
                          stringsAsFactors = F, 
                          check.names = F, 
                          row.names = 1)
    
    snp1Row <- paste0(currentSample$name,'_SNP1')
    snp2Row <- paste0(currentSample$name,'_SNP2')
    
    if(!snp1Row %in% rownames(exonSnpDF)){
      exonSnpDF[nrow(exonSnpDF)+1,] <- NA
      rownames(exonSnpDF)[nrow(exonSnpDF)] <- snp1Row
      exonSnpDF[nrow(exonSnpDF)+1,] <- NA
      rownames(exonSnpDF)[nrow(exonSnpDF)] <- snp2Row
    }
    
    if( any(exonIndelBoolVect) ){
      lapply(exonIndelList, function(x){
        snp1List <- x$SNP1
        snp2List <- x$SNP2
        
        combNames <- unique(names(snp1List), names(snp2List))
        
        # Remove INDEL positions from exon datatable
        exonDT <<- exonDT[ !feat %in% combNames ]
        
        # Write INDEL positions to SNP DF
        exonSnpDF[snp1Row,names(snp1List)] <<- unlist(snp1List)
        exonSnpDF[snp2Row,names(snp2List)] <<- unlist(snp2List)
        
        return(NULL)
      })
    }
    
    # Write the exon SNPs to the dataframe
    exonSnpDF[snp1Row, unlist(exonDT$feat)] <- unlist(exonDT$SNP1)
    exonSnpDF[snp2Row, unlist(exonDT$feat)] <- unlist(exonDT$SNP2)
    
    exonSnpDF <- allele.filter.process_del_index(exonSnpDF, bedDelIndex, snp1Row, snp2Row)
    
    write.csv(exonSnpDF, filterSnpDFList[[currentLocus]]$exonSNPs$csvPath)
    
    # ----- INTRON processing -----
    cat('\n\tProcessing introns')
    intronDT <- vcfDT[ Locus == realLocus ][ featLab %in% otherFeatNameVect ]
    
    if(nrow(intronDT) == 0){
      next
    }
    
    # Split up SNPs from INDELs
    intronDTList <- general.VCF_sep_INDEL(intronDT)
    intronDT <- intronDTList$nodelDT
    intronIndelDT <- intronDTList$indelDT
    
    intronIndelBoolVect <- !sapply(intronIndelDT$genoVect, function(x) all(x == 1))
    
    # Write intron INDEL table
    if( any(intronIndelBoolVect) ){
      intronIndelDT <- intronIndelDT[intronIndelBoolVect,]
      
      intronIndelList <- lapply(1:nrow(intronIndelDT), function(i){
        x <- intronIndelDT[i,]
        allele.formatIndelSnps(x$REF, x$SNP1, x$SNP2, x$featLab, x$featCoord)
      })
      
      indelPath <- file.path(currentSample$filterRefDirectory, paste0(currentLocus,'_intronINDELs.tsv'))
      
      appendBool <- FALSE
      
      write.table(intronIndelDT[,c('CHROM','POS','ID',
                                   'REF','ALT','QUAL',
                                   'FILTER','INFO','FORMAT',
                                   'GENO','DP','SNP1',
                                   'SNP2','feat','featLab',
                                   'featCoord','Locus')], 
                  file=indelPath, sep='\t', quote = F, 
                  row.names = F, append=appendBool, col.names = !appendBool)
    }
    
    intronSnpDF <- read.csv(filterSnpDFList[[currentLocus]]$intronSNPs$csvPath, 
                          stringsAsFactors = F, 
                          check.names = F, 
                          row.names = 1)
    
    snp1Row <- paste0(currentSample$name,'_SNP1')
    snp2Row <- paste0(currentSample$name,'_SNP2')
    
    if(!snp1Row %in% rownames(intronSnpDF)){
      intronSnpDF[nrow(intronSnpDF)+1,] <- NA
      rownames(intronSnpDF)[nrow(intronSnpDF)] <- snp1Row
      intronSnpDF[nrow(intronSnpDF)+1,] <- NA
      rownames(intronSnpDF)[nrow(intronSnpDF)] <- snp2Row
    }
    
    if( any(intronIndelBoolVect) ){
      lapply(intronIndelList, function(x){
        snp1List <- x$SNP1
        snp2List <- x$SNP2
        
        combNames <- unique(names(snp1List), names(snp2List))
        
        # Remove INDEL positions from intron datatable
        intronDT <<- intronDT[ !feat %in% combNames ]
        
        # Write INDEL positions to SNP DF
        intronSnpDF[snp1Row,names(snp1List)] <<- unlist(snp1List)
        intronSnpDF[snp2Row,names(snp2List)] <<- unlist(snp2List)
        
        return(NULL)
      })
    }
    
    # Write the intron SNPs to the dataframe
    intronSnpDF[snp1Row, unlist(intronDT$feat)] <- unlist(intronDT$SNP1)
    intronSnpDF[snp2Row, unlist(intronDT$feat)] <- unlist(intronDT$SNP2)
    
    intronSnpDF <- allele.filter.process_del_index(intronSnpDF, bedDelIndex, snp1Row, snp2Row)
    
    write.csv(intronSnpDF, filterSnpDFList[[currentLocus]]$intronSNPs$csvPath)
  }
  
  return(currentSample)
}

allele.convert_filter_bed <- function(bedPath, currentLocus, filterRefFastaList, locusRefList, filterLocusConv){
  
  ## Process reference fasta file
  fastaPath <- filterRefFastaList[[currentLocus]]
  
  output.bedList <- list()
  
  ## Read the reference fasta file
  fastaList <- general.read_fasta(fastaPath)
    
  con  <- file(bedPath, open = "r")
    
  # Process each line of the bed file
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
      
    if(grepl('\t',oneLine)){
      lineVect <- strsplit(oneLine, '\t')[[1]]
    }else{
      lineVect <- strsplit(oneLine, ' ')[[1]]
    }
      
    # Name the attributes
    lineVect <- lineVect[ nchar(lineVect) > 0 ]
    names(lineVect) <- c('alleleName','startPos','endPos','featName')[1:length(lineVect)]
    
    realLocus <- filterLocusConv[[currentLocus]]
    
    # pull out the reference sequence for the current bed line
    featSeq <- substr( fastaList[[ lineVect[['alleleName']] ]], 
                       ( as.numeric( lineVect[['startPos']] ) + 1 ) , 
                       as.numeric( lineVect[['endPos']] ) )
      
    # Match the reference sequence against known allele sequence, report all allele matches
    matchedAlleleList <- sapply(locusRefList[[realLocus]]$alleleBedList, function(x){
      x[[ lineVect[[ 'featName' ]] ]]$featSeq == featSeq
    })
      
    matchedAlleleVect <- names( which( matchedAlleleList ) )
      
    # If there are matches alleles, then move on to further processiong
    if( length(matchedAlleleVect) > 0 ){
      
      # pull out the first matched allele (they should all be the same for the current feature so the exact allele doesnt matter)
      matchedAllele <- matchedAlleleVect[1]
        
      featObject <- locusRefList[[realLocus]]$alleleBedList[[matchedAllele]][[ lineVect[[ 'featName' ]] ]]
      
      snpVect <- featObject$snpVect
        
      if(length(featObject$featDelIndex) > 0){
        snpVect <- snpVect[ -featObject$featDelIndex ]
      }
        
      posVect <- as.character( (as.numeric(lineVect[[ 'startPos' ]])+1): as.numeric(lineVect[[ 'endPos' ]]) )
      
      output.bedList[[ lineVect[[ 'alleleName' ]] ]][ posVect ] <- names(snpVect)
        
    }
      
  }
    
  close(con)
  
  return(output.bedList)
}

# Filters VCF data table by minimum depth, adds SNP columns and feature/position labels
allele.filter_VCF_process <- function(vcfDT, currentSample, minDP, bedList, realLocus){
  # Remove low quality positions
  vcfDT <- vcfDT[GENO != './.']
  
  # Pull out Depth values
  dpScore <- unlist(lapply(strsplit(vcfDT$INFO,fixed(';')), function(x){
    grep('DP=',x,value=T)
  }))
  
  # Check that only depth values were pulled
  dpCheckBool <- all(grepl('DP=',dpScore))
  
  # Stop program if a non depth value was pulled
  if(!dpCheckBool){
    message('\nNon DP element found.\t',currentSample$name,'\t','filter')
    stop()
  }
  
  # Set a depth column
  vcfDT$DP <- as.numeric(tstrsplit(dpScore, fixed('DP='))[[2]])
  
  # Split up the geno call into a vector
  vcfDT$genoVect <- lapply(strsplit(tstrsplit(vcfDT$GENO, fixed(':'))[[1]], fixed('/')), function(x){
    return(as.numeric(x)+1)
  })
  
  # Filter out low depth positions
  vcfDT <- vcfDT[DP >= minDP]
  
  # Set important features as new data table columns
  vcfDT[, c('SNP1','SNP2') := general.vcfDT_set_SNPs( REF, ALT, genoVect ), by = 1:nrow(vcfDT) ]
  
  vcfDT <- vcfDT[ vcfDT$POS %in% names(bedList[[1]]), ]
  
  vcfDT[, feat:= bedList[[CHROM]][[as.character(POS)]] , by = 1:nrow(vcfDT) ]
  
  vcfDT$featLab <- tstrsplit(vcfDT$feat, '_', fixed=T)[[1]]
  vcfDT$featCoord <- tstrsplit(vcfDT$feat, '_', fixed=T)[[2]]
  vcfDT$Locus <- realLocus
  
  return(vcfDT)
}

# Initialize SNP dataframes for storing aligned SNPs
allele.initialize_filter_SNP_tables <- function(alleleFileDirectory, locusRefList, tableNameVect, filterLocusConv){
  
  output.snpDFList <- list()
  
  cat('\nInitializing SNP tables...')
  
  # Create a SNP df for each locus
  for( currentLocus in tableNameVect ){
    cat('',currentLocus)
    
    realLocus <- filterLocusConv[[currentLocus]]
    
    # Use the iter_1 reference allele to pull out all feature lengths
    featLenList <- sapply(locusRefList[[realLocus]]$alleleBedList[[1]], function(x){
      featSeq <- x$featSeq
      delIndex <- x$featDelIndex
      seqLength <- nchar(general.del_insert(featSeq, delIndex))
      return(list( seqLength ))
    })
    
    # Pull out exon names
    exonFeatNameVect <- grep('E',names(featLenList),fixed=T,value=T)
    exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
    
    # Pull out other feature names
    otherFeatNameVect <- setdiff(names(featLenList), exonFeatNameVect)
    
    # Separate exon features from other features
    exonLenList <- featLenList[exonFeatNameVect]
    otherLenList <- featLenList[otherFeatNameVect]
    
    # Append the feature label to each position within the feature
    exonLabList <- sapply(names(exonLenList), function(x){
      featEnd <- exonLenList[[x]]
      return(paste0(x,'_',1:as.numeric(featEnd)))
    })
    
    otherLabList <- sapply(names(otherLenList), function(x){
      featEnd <- otherLenList[[x]]
      return(paste0(x,'_',1:as.numeric(featEnd)))
    })
    
    # Unlist the appended labels
    exonLabVect <- unlist(exonLabList, use.names = F)
    otherLabVect <- unlist(otherLabList, use.names = F)
    
    # Initialize CSV paths for saving the dataframe
    exonDFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_exonSNPs.csv'))
    intronDFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_intronSNPs.csv'))
    
    # Create the exon dataframe and save it
    exonSnpDF <- data.frame(matrix(NA,nrow=1,ncol=length(exonLabVect)), check.names = F, stringsAsFactors = F)
    rownames(exonSnpDF) <- 'coord'
    exonSnpDF['coord',] <- exonLabVect
    colnames(exonSnpDF) <- exonLabVect
    write.csv(exonSnpDF, exonDFPath)
    
    # Create the intron (+UTR & PE) dataframe and save it
    otherSnpDF <- data.frame(matrix(NA,nrow=1,ncol=length(otherLabVect)), check.names = F, stringsAsFactors = F)
    rownames(otherSnpDF) <- 'coord'
    otherSnpDF['coord',] <- otherLabVect
    colnames(otherSnpDF) <- otherLabVect
    write.csv(otherSnpDF, intronDFPath)
    
    # Save the CSV paths and dataframes to a return list
    output.snpDFList[[currentLocus]] <- list()
    output.snpDFList[[currentLocus]][['exonSNPs']] <- list('csvPath'=exonDFPath)
    output.snpDFList[[currentLocus]][['intronSNPs']] <- list('csvPath'=intronDFPath)
    output.snpDFList[[currentLocus]][['failure']] <- FALSE
  }
  
  cat('\nFinished. Exon and intron SNP tables located in',alleleFileDirectory)
  return(output.snpDFList)
}


filterLocusConv <- list(
  'KIR3DL3'='KIR3DL3',
  'KIR3DL2'='KIR3DL2',
  'KIR3DS1'='KIR3DS1',
  'KIR3DS1het'='KIR3DS1',
  'KIR3DL1het'='KIR3DL1',
  'KIR3DL1'='KIR3DL1',
  'KIR2DS35'='KIR2DS5',
  'KIR2DS4'='KIR2DS4',
  'KIR2DS3'='KIR2DS3',
  'KIR2DP1'='KIR2DP1',
  'KIR2DL5'='KIR2DL5',
  'KIR2DL4'='KIR2DL4',
  'KIR2DL2'='KIR2DL2',
  'KIR2DL23'='KIR2DL3',
  'KIR2DL3'='KIR2DL3',
  'KIR2DL1'='KIR2DL1'
)

filterRefFastaList <- list(
  'KIR3DL3'='Resources/genotype_resources/filters/KIR3DL3/3DL3long.fas',
  'KIR3DL2'='Resources/genotype_resources/filters/KIR3DL2/3DL2long.fas',
  'KIR3DS1'='Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT.fas',
  'KIR3DS1het'='Resources/genotype_resources/filters/KIR3DL1/3DS1longCAT.fas',
  'KIR3DL1het'='Resources/genotype_resources/filters/KIR3DL1/3DL1longCAT.fas',
  'KIR3DL1'='Resources/genotype_resources/filters/KIR3DL1/3DL1longii.fas',
  'KIR2DS35'='Resources/genotype_resources/filters/KIR2DS35/2DS5long.fas',
  'KIR2DS4'='Resources/genotype_resources/filters/KIR2DS4/2DS4long.fas',
  'KIR2DS3'='Resources/genotype_resources/filters/KIR2DS35/2DS3.fas',
  'KIR2DP1'='Resources/genotype_resources/filters/KIR2DP1/2DP1.fas',
  'KIR2DL5'='Resources/genotype_resources/filters/KIR2DL5/2DL5B.fas',
  'KIR2DL4'='Resources/genotype_resources/filters/KIR2DL4/2DL4FH5.fas',
  'KIR2DL2'='Resources/genotype_resources/filters/KIR2DL23/2DL2long.fas',
  'KIR2DL23'='Resources/genotype_resources/filters/KIR2DL23/2DL3long.fas',
  'KIR2DL3'='Resources/genotype_resources/filters/KIR2DL23/2DL3long.fas',
  'KIR2DL1'='Resources/genotype_resources/filters/KIR2DL1/2DL1AClong.fas'
)

'
no current logic for overlapping indels

finished with filter SNP processing, next move to sample looping
  - KIR2DL4 filter SNP df writing failed
  - KIR2DL1 allele calling not optimal ( IND00001 allele mismatch by 1 SNP )

  - KIR3DL3 filter allele calling appears working
  - KIR3DL2 filter allele calling appears working
  - KIR3DL1hom filter allele calling appears working
  - KIR2DP1 allele calling appears working
  - KIR2DS4 allele calling appears working
  - KIR2DL2hom allele calling appears working
  - KIR2DL3hom allele calling appears working
  - KIR2DL5 allele calling appears working (not exact match tho)
  - KIR2DS35 allele calling appears working
  - KIR2DS3 working but added new SNP that did not appear in original PING results
  - KIR2DL23 allele calling appears working (need to modify allele call workflow for compatibility with combined genes)

  - Sill left to test
    - KIR3DL1S1het

'

allele.generate_L23_SNP_df <- function( kirLocusFeatureNameList, knownSnpDFList ){
  
  ## Pull out L2 exon position names
  currentLocus <- 'KIR2DL2'
  exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[ currentLocus ]],fixed=T,value=T)
  exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
  
  otherFeatNameVect <- setdiff(kirLocusFeatureNameList[[ currentLocus ]], exonFeatNameVect)
  
  colLabels <- tstrsplit(colnames(knownSnpDFList[[ currentLocus ]]$snpDF), '_', fixed=T )[[1]]
  exonLabelIndex <- which( colLabels %in% exonFeatNameVect )
  
  exonLabelVect2DL2 <- colnames( knownSnpDFList[[ currentLocus ]]$snpDF )[exonLabelIndex]
  
  ## Pull out L3 exon position names
  currentLocus <- 'KIR2DL3'
  exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[ currentLocus ]],fixed=T,value=T)
  exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
  
  otherFeatNameVect <- setdiff(kirLocusFeatureNameList[[ currentLocus ]], exonFeatNameVect)
  
  colLabels <- tstrsplit(colnames(knownSnpDFList[[ currentLocus ]]$snpDF), '_', fixed=T )[[1]]
  exonLabelIndex <- which( colLabels %in% exonFeatNameVect )
  
  exonLabelVect2DL3 <- colnames( knownSnpDFList[[ currentLocus ]]$snpDF )[exonLabelIndex]
  
  # Isolate L2 and L3 exonic SNPs
  L2DF <- knownSnpDFList[[ 'KIR2DL2' ]]$snpDF[,exonLabelVect2DL2]
  L3DF <- knownSnpDFList[[ 'KIR2DL3' ]]$snpDF[,exonLabelVect2DL3]
  
  # KIR2DL2 positions to add (post positions get shifted by 3)
  L2RemCols <- c('E7_58','E7_59','E7_60')
  L2AddCols <- c('E7_103','E7_104','E7_105')
  
  # KIR2DL3 positions to add
  L3AddCols <- c('E9_154', 'E9_155', 'E9_156', 'E9_157', 'E9_158', 'E9_159', 'E9_160',
    'E9_161', 'E9_162', 'E9_163', 'E9_164', 'E9_165', 'E9_166', 'E9_167',
    'E9_168', 'E9_169', 'E9_170', 'E9_171', 'E9_172', 'E9_173', 'E9_174',
    'E9_175', 'E9_176', 'E9_177')
  
  # Add in L3 alleles to L23 dataframe
  L23DF <- cbind(L3DF,L2DF[1,L3AddCols])
  
  # Process L2 cols to add
  L2DF[,L2AddCols] <- NA
  oldL2Cols <- paste0('E7_',58:102)
  newL2Cols <- paste0('E7_',61:105)
  
  L2DF[,newL2Cols] <- L2DF[,oldL2Cols]
  L2DF[,L2RemCols] <- L3DF[1,L2RemCols]
  
  # Add L2 alleles to L23 dataframe
  L23DF <- rbind(L23DF, L2DF)
  
  knownSnpDFList[['KIR2DL23']] <- list(snpDF=L23DF, csvPath='')
  return( knownSnpDFList )
}

allele.generate_S35_SNP_df <- function( kirLocusFeatureNameList, knownSnpDFList ){
  ## Pull out S3 exon position names
  currentLocus <- 'KIR2DS3'
  exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[ currentLocus ]],fixed=T,value=T)
  exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
  
  otherFeatNameVect <- setdiff(kirLocusFeatureNameList[[ currentLocus ]], exonFeatNameVect)
  
  colLabels <- tstrsplit(colnames(knownSnpDFList[[ currentLocus ]]$snpDF), '_', fixed=T )[[1]]
  exonLabelIndex <- which( colLabels %in% exonFeatNameVect )
  
  exonLabelVect2DS3 <- colnames( knownSnpDFList[[ currentLocus ]]$snpDF )[exonLabelIndex]
  
  ## Pull out S5 exon position names
  currentLocus <- 'KIR2DS5'
  exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[ currentLocus ]],fixed=T,value=T)
  exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
  
  otherFeatNameVect <- setdiff(kirLocusFeatureNameList[[ currentLocus ]], exonFeatNameVect)
  
  colLabels <- tstrsplit(colnames(knownSnpDFList[[ currentLocus ]]$snpDF), '_', fixed=T )[[1]]
  exonLabelIndex <- which( colLabels %in% exonFeatNameVect )
  
  exonLabelVect2DS5 <- colnames( knownSnpDFList[[ currentLocus ]]$snpDF )[exonLabelIndex]
  
  # Add L2 alleles to L23 dataframe
  S35DF <- rbind( knownSnpDFList[[ 'KIR2DS3' ]]$snpDF[,exonLabelVect2DS3],
                  knownSnpDFList[[ 'KIR2DS5' ]]$snpDF[,exonLabelVect2DS5] )
  
  knownSnpDFList[['KIR2DS35']] <- list(snpDF=S35DF, csvPath='')
  return( knownSnpDFList )
}

# Add KIR2DL23 combined SNP df to the knownSnpDFList
knownSnpDFList <- allele.generate_L23_SNP_df( kirLocusFeatureNameList, knownSnpDFList )

knownSnpDFList <- allele.generate_S35_SNP_df( kirLocusFeatureNameList, knownSnpDFList )
