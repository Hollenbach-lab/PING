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
minDP <- 2

alignmentFileDirectory

currentSample <- sampleList[[1]]

currentSample$iterVCFPathList
currentSample$iterRefDirectory
currentSample$refIterVect

# ---------- ALLELE CALLING -----------

' TO DO
Deletion position processing
Intron SNP consolidation
'

referenceAlleleDF <- read.csv('Resources/genotype_resources/master_haplo_iteration_testing_v4.csv',row.names=1,stringsAsFactors = F)

# Create directory for storing files relating to allele calling
alleleFileDirectory <- file.path(resultsDirectory,'alleleFiles')
if(!file.exists(alleleFileDirectory)){
  dir.create(alleleFileDirectory)
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

allele.call_allele <- function(currentSample, currentLocus, alleleFileDirectory, knownSnpDFList, newAlleleDFPath ){
  
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
  
  cat('\nBest score:',bestScoreInt)
  cat('\nAllele match(es):',allMatchingAlleleStr)
  
  currentSample[['iterAlleleCallDF']]['allele_call',currentLocus] <- allMatchingAlleleStr
  currentSample[['iterAlleleCallDF']]['mismatch_score',currentLocus] <- bestScoreInt
  currentSample[['iterAlleleCallDF']]['new_snps',currentLocus] <- newSnpInt
  
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

# Iter SNP consolidation workflow
sampleList[1:10] <- sapply(sampleList[1:10], function(x){
  x <- allele.iter_alignments_to_snp_dfs(x, locusRefList, referenceAlleleDF, minDP, kirLocusFeatureNameList)
  x <- allele.combine_iter_snps(x, snpDFList)
  x <- allele.setup_allele_call_df(x)
  
  cat('\n\nFinding allele matches for',x$name)
  for(currentLocus in rownames(x$refAlleleDF)){
    x <- allele.call_allele(x, currentLocus, alleleFileDirectory, knownSnpDFList, alleleDFPathList$newAllelePath )
  }
  cat('\nWriting allele matches to', alleleDFPathList$alleleCallPath )
  x <- allele.save_call( x, alleleDFPathList$alleleCallPath )
})

