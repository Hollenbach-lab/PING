library(data.table)
library(stringr)
library(methods)
library(R.utils)
library(gtools)

# rawFastqDirectory
# fastqPattern
# threads
# resultsDirectory
# shortNameDelim
# minDP <- 10
# 
# alignmentFileDirectory
# 
# currentSample <- sampleList[[6]]
# 
# currentSample$iterVCFPathList
# currentSample$iterRefDirectory
# currentSample$refIterVect

# ---------- ALLELE CALLING -----------
# maskedPosList <- list('KIR2DS4'=paste0('E5_',84:105),'KIR2DL4'='E7_105')

# referenceAlleleDF <- read.csv('Resources/genotype_resources/master_haplo_iteration_testing_v10.csv',row.names=1,stringsAsFactors = F)

# Create directory for storing files relating to allele calling
# alleleFileDirectory <- file.path(resultsDirectory,'alleleFiles')
# if(!file.exists(alleleFileDirectory)){
#   dir.create(alleleFileDirectory)
# }
# 
# filterAlleleFileDirectory <- file.path(alleleFileDirectory,'filterFiles')
# if(!file.exists(filterAlleleFileDirectory)){
#   dir.create(filterAlleleFileDirectory)
# }

# ---- Function base -----
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

# Initialize SNP dataframes for storing aligned SNPs
allele.initialize_SNP_tables <- function(alleleFileDirectory, locusRefList, referenceAlleleDF, workflow){
  
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
    
    # Initialize CSV paths for saving the dataframe
    exonDFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_exonSNPs.csv'))
    intronDFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_intronSNPs.csv'))
    
    ## 7/26/2020 addition - depth dataframes (meant for iter alignments)
    # Initialize CSV paths for saving the dataframe
    exon.DP4.DFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_exonDP4.csv'))
    intron.DP4.DFPath <- file.path(alleleFileDirectory,paste0(currentLocus,'_intronDP4.csv'))
    
    # Save the CSV paths and dataframes to a return list
    output.snpDFList[[currentLocus]] <- list()
    output.snpDFList[[currentLocus]][['exonSNPs']] <- list('csvPath'=exonDFPath)
    output.snpDFList[[currentLocus]][['exonDP4']] <- list('csvPath'=exon.DP4.DFPath)
    output.snpDFList[[currentLocus]][['intronSNPs']] <- list('csvPath'=intronDFPath)
    output.snpDFList[[currentLocus]][['intronDP4']] <- list('csvPath'=intron.DP4.DFPath)
    output.snpDFList[[currentLocus]][['failure']] <- FALSE
    
    resultExists.bool <- file.exists(exonDFPath) & file.exists(intronDFPath) & file.exists(exon.DP4.DFPath) & file.exists(intron.DP4.DFPath)
    
    if( resultExists.bool ){
      cat('\nFound SNP dataframes in',alleleFileDirectory,'~ skipping dataframe generation.\n')
      next
    }
    
    currentAllele <- alleleList[[currentLocus]]
    
    # Use the iter_1 reference allele to pull out all feature lengths
    featLenList <- sapply(locusRefList[[currentLocus]]$alleleBedList[[currentAllele]], function(x){
      featSeq <- x$featSeq
      delIndex <- x$featDelIndex
      seqLength <- nchar(general.del_insert(featSeq, delIndex))
      
      if( workflow == 'iter' ){
        if( x$featName == '5UTR' | x$featName == '3UTR' ){
          seqLength <- 1000
        }
      }
      
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
    
    # Create the exon dataframe and save it
    exonSnpDF <- data.frame(matrix(NA,nrow=1,ncol=length(exonLabVect)), check.names = F, stringsAsFactors = F)
    rownames(exonSnpDF) <- 'coord'
    exonSnpDF['coord',] <- exonLabVect
    colnames(exonSnpDF) <- exonLabVect
    write.csv(exonSnpDF, exonDFPath)
    
    write.csv(exonSnpDF, exon.DP4.DFPath)
    
    # Create the intron (+UTR & PE) dataframe and save it
    otherSnpDF <- data.frame(matrix(NA,nrow=1,ncol=length(otherLabVect)), check.names = F, stringsAsFactors = F)
    rownames(otherSnpDF) <- 'coord'
    otherSnpDF['coord',] <- otherLabVect
    colnames(otherSnpDF) <- otherLabVect
    write.csv(otherSnpDF, intronDFPath)
    
    write.csv(otherSnpDF, intron.DP4.DFPath)
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
    
    if(featName == '5UTR' | featName == '3UTR'){ # hack for UTR extensions
      posVect <- as.character( featStart:featEnd )
      featLength <- length(posVect)
      featPosVect <- as.character( 1:featLength )
    }else{
      delFeatSeq <- general.del_insert(featObject$featSeq, featObject$featDelIndex)
      
      featLength <- nchar(delFeatSeq)
      
      posVect <- as.character( featStart:featEnd )
      
      featPosVect <- as.character( 1:featLength )
      
      if( length(featObject$featDelIndex) > 0 ){
        featPosVect <- featPosVect[ -featObject$featDelIndex ]
      }
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
  sampleSnpDFList <- allele.initialize_SNP_tables(currentSample$iterRefDirectory, locusRefList, currentSample$refAlleleDF, 'iter')
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
      cat('\n',currentLocus)
      
      cat('\tProcessing EXONs')
      
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
        cat('\tINDELs')
        exonIndelDT <- exonIndelDT[exonIndelBoolVect,]
        exonIndelDT$refIter <- refIter
        
        exonIndelRepList <- lapply(1:nrow(exonIndelDT), function(i){
          x <- exonIndelDT[i,]
          allele.formatIndelSnps(x$REF, x$SNP1, x$SNP2, x$featLab, x$featCoord)
        })
        
        indelPath <- file.path(currentSample$iterRefDirectory, paste0(currentLocus,'_exonINDELs.tsv'))
        
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
      
      exon.DP4.DF <- read.csv(sampleSnpDFList[[currentLocus]]$exonDP4$csvPath, 
                             stringsAsFactors = F, 
                             check.names = F, 
                             row.names = 1)
      
      currentRow <- paste0(refIter,'_SNP1')
      if(!currentRow %in% rownames(exonSnpDF)){
        exonSnpDF[nrow(exonSnpDF)+1,] <- NA
        rownames(exonSnpDF)[nrow(exonSnpDF)] <- paste0(refIter,'_SNP1')
        exonSnpDF[nrow(exonSnpDF)+1,] <- NA
        rownames(exonSnpDF)[nrow(exonSnpDF)] <- paste0(refIter,'_SNP2')
        
        exon.DP4.DF[nrow(exon.DP4.DF)+1,] <- NA
        rownames(exon.DP4.DF)[nrow(exon.DP4.DF)] <- paste0(refIter,'_SNP1')
        exon.DP4.DF[nrow(exon.DP4.DF)+1,] <- NA
        rownames(exon.DP4.DF)[nrow(exon.DP4.DF)] <- paste0(refIter,'_SNP2')
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
      
      # Write the exon DP4 to the dataframe
      exon.DP4.DF[paste0(refIter,'_SNP1'), unlist(exonDT$feat)] <- unlist(exonDT$SNP1.DP4)
      exon.DP4.DF[paste0(refIter,'_SNP2'), unlist(exonDT$feat)] <- unlist(exonDT$SNP2.DP4)
      
      snp1Row <- paste0(refIter,'_SNP1')
      snp2Row <- paste0(refIter,'_SNP2')
      
      exonSnpDF <- allele.process_del_index(exonSnpDF, bedDelIndex, unique(exonDT$CHROM), snp1Row, snp2Row)
      
      write.csv(exonSnpDF, sampleSnpDFList[[currentLocus]]$exonSNPs$csvPath)
      write.csv(exon.DP4.DF, sampleSnpDFList[[currentLocus]]$exonDP4$csvPath)
      
      cat('\tINTRONs')
      intronDT <- vcfDT[ Locus == currentLocus ][ featLab %in% otherFeatNameVect ]
      
      # Split up SNPs from INDELs
      intronDTList <- general.VCF_sep_INDEL(intronDT)
      intronDT <- intronDTList$nodelDT
      intronIndelDT <- intronDTList$indelDT
      
      intronIndelBoolVect <- !sapply(intronIndelDT$genoVect, function(x) all(x == 1))
      
      # Write intron INDEL table
      if( any(intronIndelBoolVect) ){
        cat('\tINDELs')
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
      
      intron.DP4.DF <- read.csv(sampleSnpDFList[[currentLocus]]$intronDP4$csvPath, 
                              stringsAsFactors = F, 
                              check.names = F, 
                              row.names = 1)
      
      currentRow <- paste0(refIter,'_SNP1')
      if(!currentRow %in% rownames(intronSnpDF)){
        intronSnpDF[nrow(intronSnpDF)+1,] <- NA
        rownames(intronSnpDF)[nrow(intronSnpDF)] <- paste0(refIter,'_SNP1')
        intronSnpDF[nrow(intronSnpDF)+1,] <- NA
        rownames(intronSnpDF)[nrow(intronSnpDF)] <- paste0(refIter,'_SNP2')
        
        intron.DP4.DF[nrow(intron.DP4.DF)+1,] <- NA
        rownames(intron.DP4.DF)[nrow(intron.DP4.DF)] <- paste0(refIter,'_SNP1')
        intron.DP4.DF[nrow(intron.DP4.DF)+1,] <- NA
        rownames(intron.DP4.DF)[nrow(intron.DP4.DF)] <- paste0(refIter,'_SNP2')
      }
      
      if( any(intronIndelBoolVect) ){
        lapply(intronIndelRepList, function(x){
        snp1List <- x$SNP1
        snp2List <- x$SNP2
        
        x.feat <- unique( tstrsplit(names(snp2List),'_',fixed=T)[[1]] )
        if(x.feat == '3UTR'){
          x.1.pos <- tstrsplit(names(snp1List),'_',fixed=T)[[2]]
          x.2.pos <- tstrsplit(names(snp2List),'_',fixed=T)[[2]]
          
          # Skip INDEL processing that happens at the end of 3'UTR
          if( any( as.numeric( unique(c(x.1.pos, x.2.pos)) ) > 950 ) ) {
            return(NULL)
          }
        }
        
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
      
      # Write the intron DP4 to the dataframe
      intron.DP4.DF[paste0(refIter,'_SNP1'), unlist(intronDT$feat)] <- unlist(intronDT$SNP1.DP4)
      intron.DP4.DF[paste0(refIter,'_SNP2'), unlist(intronDT$feat)] <- unlist(intronDT$SNP2.DP4)
      
      intronSnpDF <- allele.process_del_index(intronSnpDF, bedDelIndex, unique(exonDT$CHROM), snp1Row, snp2Row)
      
      write.csv(intronSnpDF, sampleSnpDFList[[currentLocus]]$intronSNPs$csvPath)
      write.csv(intron.DP4.DF, sampleSnpDFList[[currentLocus]]$intronDP4$csvPath)
    }
    cat('\n')
  }
  
  currentSample[['iterSnpDFList']] <- sampleSnpDFList
  return(currentSample)
}

allele.iter_combine_KIR2DL23 <- function( currentSample, knownSnpDFList, alleleFileDirectory, snpDFList ){

  snpFilePath <- file.path(alleleFileDirectory,'KIR2DL23_exonSNPs.csv')
  
  # Create the KIR2DL23 snp file if it does not exist
  if( !file.exists(snpFilePath) ){
    
    KIR2DL23SnpDF <-  knownSnpDFList$KIR2DL23$snpDF[1,,drop=F]
    KIR2DL23SnpDF[1,] <- as.character( colnames(KIR2DL23SnpDF) )
    rownames(KIR2DL23SnpDF)[1] <- 'coord'
    
    write.csv(x = KIR2DL23SnpDF, file = snpFilePath )
  
  }
  
  KIR2DL23SnpDF <- read.csv(snpFilePath, 
                            stringsAsFactors = F, 
                            check.names=F,
                            row.names=1 )
  
  row1 <- paste0(currentSample$name,'_SNP1')
  row2 <- paste0(currentSample$name,'_SNP2')
  rowVect <- c(row1, row2)

  if( !row1 %in% rownames(KIR2DL23SnpDF) ){
    KIR2DL23SnpDF[nrow(KIR2DL23SnpDF)+1,] <- NA
    rownames(KIR2DL23SnpDF)[nrow(KIR2DL23SnpDF)] <- row1
    KIR2DL23SnpDF[nrow(KIR2DL23SnpDF)+1,] <- NA
    rownames(KIR2DL23SnpDF)[nrow(KIR2DL23SnpDF)] <- row2
  }
  
  # Define the KIR2DL2 exon SNPs DF path
  KIR2DL2SnpDFPath <- snpDFList[['KIR2DL2']]$exonSNPs$csvPath
  
  # Define the KIR2DL3 exon SNPs DF path
  KIR2DL3SnpDFPath <- snpDFList[['KIR2DL3']]$exonSNPs$csvPath
  
  # Read in the KIR2DL2 exon SNPs
  KIR2DL2SnpDF <- read.csv(KIR2DL2SnpDFPath,
                                stringsAsFactors = F,
                                check.names = F,
                                row.names = 1)
    
  # Read in the KIR2DL2 exon SNPs
  KIR2DL3SnpDF <- read.csv(KIR2DL3SnpDFPath,
                           stringsAsFactors = F,
                           check.names = F,
                           row.names = 1)
  
  if( !all( rowVect %in% rownames(KIR2DL2SnpDF) ) ){
    stop('Did not find [',paste(row1,row2),'] in ',KIR2DL2SnpDFPath)
  }
  
  if( !all( rowVect %in% rownames(KIR2DL3SnpDF) ) ){
    stop('Did not find [',paste(row1,row2),'] in ',KIR2DL3SnpDFPath)
  }
  
  oldL2Cols <- paste0('E7_',58:102)
  newL2Cols <- paste0('E7_',61:105)
  
  # KIR2DL3 positions to add
  L3AddCols <- c('E9_154', 'E9_155', 'E9_156', 'E9_157', 'E9_158', 'E9_159', 'E9_160',
                 'E9_161', 'E9_162', 'E9_163', 'E9_164', 'E9_165', 'E9_166', 'E9_167',
                 'E9_168', 'E9_169', 'E9_170', 'E9_171', 'E9_172', 'E9_173', 'E9_174',
                 'E9_175', 'E9_176', 'E9_177')
  
  # ----- KIR2DL3 will be row1
  KIR2DL23SnpDF[row1,L3AddCols] <- knownSnpDFList$KIR2DL23$snpDF[1,L3AddCols]
  
  cat('\n\tAdding KIR2DL3 SNPs')
  onlyL3Cols <- paste0('E7_',103:105)
  for(colID in colnames(KIR2DL3SnpDF)){
    
    nonNAIndex <- !is.na( KIR2DL3SnpDF[rowVect,colID] )
    
    if( !any( nonNAIndex ) ){
      KIR2DL23SnpDF[rowVect,colID] <- NA
    }else{
      L3Snps <- unique( KIR2DL3SnpDF[ rowVect[nonNAIndex], colID] )
      
      KIR2DL23SnpDF[rowVect[1:length(L3Snps)],colID] <- L3Snps
    }
    
  }
  
  # ----- Add in KIR2DL2 SNPs
  cat('\n\tAdding KIR2DL2 SNPs')
  oldL2Cols <- paste0('E7_',58:102)
  newL2Cols <- paste0('E7_',61:105)
  KIR2DL23SnpDF[row2,c('E7_58','E7_59','E7_60')] <- c('T','C','C')
  for(colID in colnames(KIR2DL2SnpDF)){
    
    oldColID <- colID
    
    if(colID %in% oldL2Cols){
      idIndex <- which( colID == oldL2Cols )
      oldColID <- colID
      colID <- newL2Cols[idIndex]
    }
    
    nonNAIndex <- !is.na(KIR2DL2SnpDF[rowVect,oldColID])
      
    if( any( nonNAIndex ) ){
      nonNAL23Index <- !is.na(KIR2DL23SnpDF[rowVect,colID])
      
      L23Snps <- unique( KIR2DL2SnpDF[ rowVect[nonNAIndex], oldColID] )
      L23Snps <- unique( c( L23Snps, KIR2DL23SnpDF[ rowVect[nonNAL23Index], colID] ) )
      
      KIR2DL23SnpDF[rowVect,colID] <- L23Snps
    }
  }
  
  write.csv(x = KIR2DL23SnpDF, file = snpFilePath )
  
  currentSample[['iterSnpDFList']][['KIR2DL23']] <- list()
  currentSample[['iterSnpDFList']][['KIR2DL23']][['exonSNPs']] <- list('csvPath'=snpFilePath)
  currentSample[['iterSnpDFList']][['KIR2DL23']][['failure']] <- FALSE
  
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
  
  # Set the DP4
  cat('\n\t\t\tSetting DP4 sums')
  vcfDT[, c('SNP1.DP4','SNP2.DP4') := general.vcfDT_set_DP4( INFO ), by = 1:nrow(vcfDT) ]
  
  # Split up the geno call into a vector
  vcfDT$genoVect <- lapply(strsplit(tstrsplit(vcfDT$GENO, fixed(':'))[[1]], fixed('/')), function(x){
    return(as.numeric(x)+1)
  })
  
  
  # Filter out low depth positions
  vcfDT <- vcfDT[DP >= minDP]
  
  # Filter out bad REF call positions
  vcfDT <-  vcfDT[ !(SNP2.DP4 > SNP1.DP4 & GENO == '0/0') ]
  
  cat('\n\t\t\tSetting SNP calls')
  
  # Set important features as new data table columns
  #vcfDT[, c('SNP1','SNP2') := general.vcfDT_set_SNPs( REF, ALT, genoVect ), by = 1:nrow(vcfDT) ]
  
  tempSnpList <- lapply(1:nrow(vcfDT),function(x){
    
    if( grepl( fixed(','), vcfDT$ALT[[x]] ) ){
      tempSnpVect <- unlist( strsplit(paste(vcfDT$REF[[x]], vcfDT$ALT[[x]], sep=','), fixed(',')) )
    }else{
      tempSnpVect <- c(vcfDT$REF[[x]],vcfDT$ALT[[x]])
    }
    
    tempSnpVect[ vcfDT$genoVect[[x]] ]
  })
  
  vcfDT$SNP1 <- unlist(lapply(1:nrow(vcfDT), function(x){
    tempSnpList[[x]][1]
  }))
  
  vcfDT$SNP2 <- unlist(lapply(1:nrow(vcfDT), function(x){
    tempSnpList[[x]][2]
  }))
  
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
  
  if( grepl( fixed(','), snp2 ) ){
    rawSnpVect <- unlist( strsplit(paste(snp1, snp2, sep=','), fixed(',')) )
  }else{
    rawSnpVect <- c(snp1, snp2)
  }
  
  snpVect <- rawSnpVect[ genoVect[[1]] ]

  return( list( snpVect[1], snpVect[2] ) )
}

# set dp4 in vcf data table
general.vcfDT_set_DP4 <- function( INFO ){
  dp4Str <- tstrsplit(grep('DP4', tstrsplit(INFO, fixed(';')),value=T), fixed('DP4='))[[2]]
  
  dp4Vect <- as.numeric(strsplit(dp4Str,fixed(','))[[1]])
  
  return( list( sum(dp4Vect[1:2]), sum(dp4Vect[3:4]) ) )
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
  
  #vcfDT <- fread(vcfFile)
  
  vcfDT <- tryCatch({
    fread(vcfFile)
  },
  error=function(cond){
    message('failure to read VCF')
    return('failure')
  },
  warning=function(cond){
    message('failure to read VCF')
    return('failure')
  }
  )
  
  if(length(vcfDT) == 1){
    if(vcfDT == 'failure'){
      return(vcfDT)
    }
  }
  
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

allele.combine_iter_snps <- function(currentSample, snpDFList, dp4SnpRatio=5){
  
  if(currentSample$iterRefDirectory == 'failed'){
    return(currentSample)
  }
  
  cat('\nCombining iter snps for',currentSample$name)
  
  locusVect <- names(currentSample$iterSnpDFList)
  
  if( 'KIR2DL23' %in% locusVect ){
    locusVect <- locusVect[ locusVect != 'KIR2DL23' ]
  }
  
  for(currentLocus in locusVect){
    
    if( currentSample$iterSnpDFList[[currentLocus]]$failure == TRUE ){
      cat('\n----- skipping locus -----\n')
      next
    }
      
    cat('',currentLocus)
    masterExonSnpDFPath <- snpDFList[[currentLocus]]$exonSNPs$csvPath
    masterIntronSnpDFPath <- snpDFList[[currentLocus]]$intronSNPs$csvPath
    
    sampleExonSnpDFPath <- currentSample$iterSnpDFList[[currentLocus]]$exonSNPs$csvPath
    sampleExonDP4DFPath <- currentSample$iterSnpDFList[[currentLocus]]$exonDP4$csvPath
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
    
    sampleExonDP4DF <- read.csv(sampleExonDP4DFPath,
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
      
      snpVect <- sampleExonSnpDF[2:nrow(sampleExonSnpDF),currentCol]
      snpVect <- snpVect[!is.na(snpVect)]
      
      dp4Vect <- sampleExonDP4DF[2:nrow(sampleExonSnpDF),currentCol]
      dp4Vect <- as.numeric( dp4Vect[!is.na(dp4Vect)] )
      
      uniqueSnpVect <- unique(snpVect)
      
      if( length(dp4Vect) > 0 && length( uniqueSnpVect ) > 1 ){
        #maxDP4 <- max( dp4Vect )
        #dp4Thresh <- as.integer( maxDP4*dp4SnpRatio )
        dp4Thresh <- dp4SnpRatio
        
        snpVect <- unique( snpVect[ dp4Vect >= dp4Thresh ] )
      }else{
        snpVect <- uniqueSnpVect
      }
      
      #snpVect <- uniqueSnpVect
      
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
      locusSnpDF <- cbind.data.frame(locusSnpDF, featSnpMat, stringsAsFactors=F)
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
allele.setup_iter_allele_call_df <- function(currentSample){
  
  currentSampleLociVect <- rownames(currentSample$refAlleleDF)
  
  if( 'KIR2DL2' %in% currentSampleLociVect & 'KIR2DL3' %in% currentSampleLociVect ){
    currentSampleLociVect <- c(currentSampleLociVect, 'KIR2DL23')
  }
  
  currentSample[['iterAlleleCallDF']] <- data.frame(matrix(NA,nrow=3,ncol=length(currentSampleLociVect)))
  
  colnames(currentSample[['iterAlleleCallDF']]) <- currentSampleLociVect
  
  rownames(currentSample[['iterAlleleCallDF']]) <- c('allele_call','mismatch_score','new_snps')
  
  return(currentSample)
}

allele.setup_filter_allele_call_df <- function(currentSample){
  
  currentSampleLociVect <- currentSample[[ 'samplePresentLocusVect' ]]
  currentSampleLociVect <- gsub('KIR3DL1het', 'KIR3DL1', currentSampleLociVect)
  currentSampleLociVect <- gsub('KIR3DS1het', 'KIR3DS1', currentSampleLociVect)
  
  currentSample[['filterAlleleCallDF']] <- data.frame(matrix(NA,nrow=3,ncol=length(currentSampleLociVect)))
  
  colnames(currentSample[['filterAlleleCallDF']]) <- currentSampleLociVect
  
  rownames(currentSample[['filterAlleleCallDF']]) <- c('allele_call','mismatch_score','new_snps')
  
  return(currentSample)
}

allele.call_allele <- function(currentSample, currentLocus, alleleFileDirectory, knownSnpDFList, newAlleleDFPath, filterLocusConv, workflow){
  
  # Fix for KIR3DL1het and KIR3DS1het VCF files, which are treated as KIR3DL1 and KIR3DS1 in allele calling
  if( currentLocus == 'KIR3DL1het' | currentLocus == 'KIR3DS1het' ){
    currentLocus <- filterLocusConv[[currentLocus]]
  }
  
  # Failure conditions -----  
  
  if( currentSample[[ paste0(workflow,'RefDirectory') ]] == 'failed' ){
    currentSample[[ paste0(workflow,'AlleleCallDF') ]] <- 'failed'
    return(currentSample)
  }
  
  if( currentSample[[ paste0(workflow,'SnpDFList') ]][[ currentLocus ]]$failure == TRUE ){
    cat('\n----- skipping',currentLocus,'-----\n')
    currentSample[[ paste0(workflow,'AlleleCallDF') ]][,currentLocus] <- 'failed'
    return(currentSample)
  }
  
  cat('\n\nMatching alleles for',currentLocus)
  
  # load up the aligned SNP df for the current locus
  locusExonSnpDF <- read.csv(file.path(alleleFileDirectory,paste0(currentLocus,'_exonSNPs.csv')), check.names=F, stringsAsFactors = F, row.names = 1)
  
  # Define rownames for the current sample
  sampleRows <- paste0(currentSample$name,c('_SNP1','_SNP2'))
  
  # Pull out the aligned SNPs for the current sample
  sampleSnpDF <- locusExonSnpDF[sampleRows,]
  
  if( workflow == 'filter' & currentLocus %in% c( 'KIR2DL1', 'KIR2DL5' ) ){
    sampleSnpDF <- filter_contam_snps( currentSample, currentLocus, sampleSnpDF )
  }
  
  # Pull out the aligned positions that passed QC checks
  definedCols <- colnames(sampleSnpDF)[apply(sampleSnpDF, 2, function(x){ all(is_nuc(x)) })]
  
  # Take out any masked positions ( usually deletion positions with poor alignment accuracy )
  if( currentLocus %in% names(maskedPosList) & workflow == 'iter' ){
    definedCols <- definedCols[ !definedCols %in% maskedPosList[[currentLocus]] ]
    
    sampleSnpDF <- sampleSnpDF[,definedCols]
  }
  
  # Format sample DF and transform into data table
  sampleSnpDT <- as.data.table(sampleSnpDF)
  sampleSnpDT$name <- rownames(sampleSnpDF)
  setkey(sampleSnpDT, name)
  
  # Pull out any triple SNPs
  tripleSnpCols <- colnames(sampleSnpDF)[ apply(sampleSnpDF, 2, function(x){ any(x=='3SNP', na.rm=T) }) ]
  
  # Subset the known SNP df by the aligned SNP positions
  knownSnpDF <- knownSnpDFList[[currentLocus]]$snpDF[,definedCols]
  
  # Generate allele differentiating (AD) SNP df
  adSnpDF <- make_unique_pos_frame(knownSnpDF)
  
  # failure condition for all alleles being taken out by make_unique_pos_frame()
  if( ncol(adSnpDF) == 0 ){
    adSnpDF <- knownSnpDF
  }
  
  # Remove duplicate alleles
  adSnpDF <- unique(adSnpDF)
  
  excludedAlleleVect <- setdiff( rownames(knownSnpDF), rownames(adSnpDF) )
  
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
  homBool <- length(homPosVect) > 0
  
  if( !hetBool & !homBool ){
    cat('\nNo SNPs found, returning NULL call (check copy results to verify).')
   
     if( workflow == 'iter' ){
      currentSample[['iterAlleleCallDF']]['allele_call',currentLocus] <- paste0(currentLocus,'*NULL')
      currentSample[['iterAlleleCallDF']]['mismatch_score',currentLocus] <- 1
      currentSample[['iterAlleleCallDF']]['new_snps',currentLocus] <- 0
    }else if( workflow == 'filter' ){
      currentSample[['filterAlleleCallDF']]['allele_call',currentLocus] <- paste0(currentLocus,'*NULL')
      currentSample[['filterAlleleCallDF']]['mismatch_score',currentLocus] <- 1
      currentSample[['filterAlleleCallDF']]['new_snps',currentLocus] <- 0
    }
    
    return(currentSample)
  }
  
  # Bool check is fix for no AD homozygous positions
  if( homBool ){
    cat('\nScoring hom positions.')
    # Narrow down adSnpDT
    homScoreList <- sapply(rownames(adSnpDF), function(x){
      sum(sapply(homPosVect, function(y){
        addScore <- (adSnpDT[x, ..y] != sampleSnpDT[1, ..y])*1
        subScore <- (adSnpDT[x, ..y] == '*')*1
        max(0,addScore-subScore)
      }))
    })
  }
  
  homScoreBuffer <- 1
  if( !hetBool ){
    bestScoreInt <- min(homScoreList)
    bestMatchIndex <- which( homScoreList == bestScoreInt )
    bestMatchAlleleVect <- rownames(adSnpDF)[bestMatchIndex]
    
    #bestMatchAlleleMat <- combinations(length(bestMatchAlleleVect), 2, bestMatchAlleleVect, repeats.allowed = T)
    #allMatchingAlleleVect <- apply(bestMatchAlleleMat,1,paste0,collapse='+')
  }else{
    
    # Fix for no AD homozygous positions
    if( !homBool ){
      homAlleleVect <- rownames(adSnpDF)
    }else{
      homAlleleVect <- rownames(adSnpDF)[homScoreList <= ( min(homScoreList)+homScoreBuffer ) ] # Cut any alleles with more than min mismatches
    }
    
    # Generate all possible allele pairings
    possAllelePairMat <- combinations(length(homAlleleVect),2,homAlleleVect,repeats.allowed = T)
    
    # Set up possible allele data table
    possAlleleDT <- as.data.table(possAllelePairMat)
    possAlleleDT$allelePair <- apply(possAllelePairMat,1,paste0,collapse='+')
    setkey(possAlleleDT, allelePair)
    colnames(possAlleleDT)[1:2] <- c('allele1','allele2')
    
    # initialize column for storing distance scores
    possAlleleDT$distance <- 0
    
    if( homBool ){
      # First score homozygous positions for allele pairs (this lowers the computational load for the het position scoring)
      possAlleleDT[, 'distance' := allele.add_hom_score( allele1, allele2, distance, homScoreList ), by=allelePair]
      possAlleleDT <- possAlleleDT[possAlleleDT$distance <= (min(possAlleleDT$distance)+homScoreBuffer),] # Remove all allele pairings that have more than the min mismatch
    }
    
    # Score het positions for allele pairs (this can take awhile if there are many pairs)
    cat('\nScoring',nrow(possAlleleDT),'allele pairings...')
    possAlleleDT[, 'distance' := allele.pair_score_calc( allele1, allele2, distance, adSnpDT, sampleSnpDT, hetPosVect ), by=allelePair ]
    
    bestScoreInt <- min(possAlleleDT$distance)
    bestMatchIndex <- which(possAlleleDT$distance == bestScoreInt)
    
    bestMatchAlleleVect <- possAlleleDT$allelePair[bestMatchIndex]
  }
  
  formattedAlleleVect <- allele.format_call(bestMatchAlleleVect, knownSnpDF, hetBool, excludedAlleleVect )
  
  if( currentLocus == 'KIR2DL1' & workflow == 'filter' ){
    formattedAlleleVect <- allele.custom_2DL1_allele_filter( currentSample, formattedAlleleVect )
  }
  
  allMatchingAlleleStr <- paste0(formattedAlleleVect, collapse=' ')
  
  bestScoreInt <- bestScoreInt+newSnpInt+length(tripleSnpCols)
  
  # New allele catch, the !unresolved_alleleFilter flag is to bypass unresolved calls added by the custom allele filter
  if( bestScoreInt > 0 & !grepl('unresolved_alleleFilter', allMatchingAlleleStr) ){
    bestMatchAlleleStr <- paste0( bestMatchAlleleVect, collapse= ' ' )
    allele.save_new_allele( currentSample, currentLocus, knownSnpDF, sampleSnpDF, bestMatchAlleleStr, newAlleleDFPath, hetBool, newSnpInt, allMatchingAlleleStr )
  }
  
  # Print out the best score and cooresponding allele matches
  cat('\n\t\tBest score:',bestScoreInt)
  cat('\n\t\tAllele match(es):',allMatchingAlleleStr)
  
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

allele.save_new_allele <- function( currentSample, currentLocus, knownSnpDF, sampleSnpDF, allMatchingAlleleStr, newAlleleDFPath, hetBool, snpScore, formAlleleStr ){
  
  cat('\nSaving new allele call.')
  
  newDF <- read.table(newAlleleDFPath, check.names=F, stringsAsFactors = F, sep=',', header = T)
  
  if( hetBool ){
    snpInt <- 2
  }else{
    snpInt <- 1
  }
  
  mismatchPosVect <- sapply(tstrsplit(allMatchingAlleleStr, ' ',fixed=T), function(x){
    checkCols <- colnames(sampleSnpDF)
    
    alleleVect <- strsplit(x,'+',fixed=T)[[1]]
    
    scoreList <- sapply(checkCols, function(curPos){
      
      if( all( is.na(sampleSnpDF[,curPos]) ) ){
        return(0)
      }
      
      if(any(sampleSnpDF[,curPos] == '3SNP', na.rm = T)){
        return(1)
      }
      
      pos1 <- knownSnpDF[ alleleVect, curPos][1]
      pos2 <- knownSnpDF[ alleleVect, curPos][snpInt]
      
      con1 <- !all( sampleSnpDF[, curPos][1] == pos1 )
      con2 <- !all( sampleSnpDF[, curPos][2] == pos2 )
      
      con3 <- !all( sampleSnpDF[, curPos][1] == pos2 )
      con4 <- !all( sampleSnpDF[, curPos][2] == pos1 )
      
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
    
    mismatchIndex <- which( scoreList > 0 )
    
    if(length(mismatchIndex) == 0){
      cat('\n\tIn new allele calling but no mismatches found.')
      return('unknown')
    }
    
    return(names(mismatchIndex))
  })
  
  mismatchPosVect <- unique(unlist(mismatchPosVect))
  mismatchPosStr <- paste0(mismatchPosVect, collapse=' ')
  
  fun.output <- c(currentSample$name, currentLocus, mismatchPosStr, formAlleleStr, snpScore>0)
  
  newDF[ (nrow(newDF)+1) , c('sampleName','locusName','mimatchPos','bestAlleleMatch','newSnpFound') ] <- fun.output
  
  write.table( newDF, file=newAlleleDFPath, sep=',', quote=F, row.names=T, col.names=T )
  
}

allele.add_hom_score <- function( allele1, allele2, distance, homScoreList ){
  return( homScoreList[[allele1]] + homScoreList[[allele2]] + distance )
}

allele.pair_score_calc <- function( allele1, allele2, distance, adSnpDT, sampleSnpDT, adCols ){

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

allele.save_call <- function( currentSample, acDFPath, workflow ){
  
  # Failure condition
  if( currentSample[[ paste0(workflow,'RefDirectory') ]] == 'failed' ){
    currentSample[[ paste0(workflow,'AlleleCallDF') ]] <- 'failed'
    return(currentSample)
  }
  
  acDF <- read.table(acDFPath, check.names = F, stringsAsFactors = F, sep=',')
  
  
  
  for( currentLocus in colnames( currentSample[[ paste0(workflow,'AlleleCallDF') ]] ) ){
    alleleCall <- currentSample[[ paste0(workflow,'AlleleCallDF') ]]['allele_call',currentLocus]
    callScore <- currentSample[[ paste0(workflow,'AlleleCallDF') ]]['mismatch_score',currentLocus]
    snpScore <- currentSample[[ paste0(workflow,'AlleleCallDF') ]]['new_snps',currentLocus]
    
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

allele.setup_results_df <- function( locusRefList, filterLocusConv, resultsDirectory, sampleList ){
  
  output.list <- list()
  
  # set up allele call dataframes for both workflows
  for( workflow in c('filter','iter') ){
    
    # pull out all locus names specific to the workflow
    if( workflow == 'filter' ){
      locusVect <- c(unique( unlist(filterLocusConv, use.names=F) ), 'KIR2DL23', 'KIR2DS35')
    }else if ( workflow == 'iter' ){
      locusVect <- c(names(locusRefList), 'KIR2DL23')
    }
    
    # Set up dataframe for storing allele calls across all samples and loci
    acDF <- data.frame(matrix(NA,nrow=length(sampleList),ncol=length(locusVect)), stringsAsFactors = F)
    colnames(acDF) <- locusVect
    rownames(acDF) <- unlist(lapply(sampleList, function(x) x$name))
    
    alleleCallPath <- file.path(resultsDirectory, paste0(workflow,'AlleleCalls.csv'))
    write.table( acDF, file=alleleCallPath, sep=',', quote=F, row.names=T )
    
    newDF <- data.frame(matrix(NA,nrow=0,ncol=5), stringsAsFactors = F)
    colnames(newDF) <- c('sampleName','locusName','mimatchPos','bestAlleleMatch','newSnpFound')
    
    newCallPath <- file.path(resultsDirectory, paste0(workflow,'NewAlleles.csv'))
    write.table( newDF, file= newCallPath, sep=',', quote=F, row.names=T)
    
    output.list[[ workflow ]] <- list('alleleCallPath'=alleleCallPath, 'newAllelePath'=newCallPath)
  }
  
  return( output.list )
}

allele.format_call <- function( bestMatchAlleleVect, knownSnpDF, hetBool, excludedAlleleVect ){
  
  output.alleleVect <- c()
  
  if(hetBool == F){
    
    splitAlleleVect <- tstrsplit( bestMatchAlleleVect, ' ', fixed=T )[[1]]
    homAlleleVect <- c()
    
    for( alleleMatch in splitAlleleVect ){
      homAlleleVect <- c( alleleMatch, homAlleleVect )
      homAlleleVect <- c( homAlleleVect, names( which( apply( knownSnpDF[excludedAlleleVect,], 1, function(x){ all( x == knownSnpDF[alleleMatch,] ) } ) ) ) )
    }
    
    output.alleleVect <- unique( unlist( sapply( homAlleleVect, kir.allele_resolution, 5 ), use.names=F ) )
    
  }else{
    
    splitAllelePairVect <- tstrsplit( bestMatchAlleleVect, ' ', fixed=T )[[1]]
    hetAlleleVect <- c()
    
    for(allelePair in splitAllelePairVect){
      splitAlleleVect <- strsplit( allelePair, '+', fixed=T )[[1]]
      
      firstMatchVect <- names( which( apply( knownSnpDF[excludedAlleleVect,], 1, function(x) { all( x == knownSnpDF[splitAlleleVect[1],] ) } ) ) )
      firstMatchVect <- unique( c( kir.allele_resolution(splitAlleleVect[1], 5), 
                           unique( unlist( sapply( firstMatchVect, kir.allele_resolution, 5 ), use.names=F ) ) ) )
      
      secondMatchVect <- names( which( apply( knownSnpDF[excludedAlleleVect,], 1, function(x) { all( x == knownSnpDF[splitAlleleVect[2],] ) } ) ) )
      secondMatchVect <- unique( c( kir.allele_resolution(splitAlleleVect[2], 5), 
                            unique( unlist( sapply( secondMatchVect, kir.allele_resolution, 5 ), use.names=F ) ) ) )
      
      matchMat <- expand.grid(firstMatchVect,secondMatchVect)
      hetAlleleVect <- c( hetAlleleVect, apply( matchMat, 1, paste0, collapse='+' ) )
    }
    
    output.alleleVect <- hetAlleleVect
    
  }

  return( output.alleleVect )  
}

# ----- Run workflow -----
# Initialze data frames for consolidating SNPs across all samples for each gene
#snpDFList <- allele.initialize_SNP_tables(alleleFileDirectory, locusRefList, referenceAlleleDF, workflow='iter')

# # Generate known SNP df's for allele calling
# knownSnpDFList <- allele.create_allele_resources(locusRefList, alleleFileDirectory)
# 
# # Add KIR2DL23 combined SNP df to the knownSnpDFList [filter specific]
# knownSnpDFList <- allele.generate_L23_SNP_df( kirLocusFeatureNameList, knownSnpDFList )
# 
# # Add KIR2DS35 combined SNP df to the knownSnpDFList [filter specific]
# knownSnpDFList <- allele.generate_S35_SNP_df( kirLocusFeatureNameList, knownSnpDFList )

# # Iter SNP consolidation workflow + allele calling
# sampleList <- sapply(sampleList, function(x){
#   
#   if(x$name %in% skipSampleVect){
#     return(x)
#   }
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

# Initialize dataframes for consolidating SNPs across all samples for each gene
#filterSnpDFList <- allele.initialize_filter_SNP_tables(filterAlleleFileDirectory, locusRefList, names(filterLocusConv), filterLocusConv)

# # Filter SNP processing workflow + allele calling
# sampleList <- sapply(sampleList, function(x){
#   cat('\n\nProcessing filtration alignments for',x$name)
#   x <- allele.filter_alignments_to_snp_dfs( x, locusRefList, minDP, kirLocusFeatureNameList, filterRefFastaList, filterLocusConv, filterSnpDFList, knownSnpDFList )
#   
#   x <- allele.setup_filter_allele_call_df( x )
#   
#   if( x[['filterRefDirectory']] == 'failed' ){
#     return(x)
#   }
#   
#   for( currentLocus in x$samplePresentLocusVect ){
#     x <- allele.call_allele( x, currentLocus, filterAlleleFileDirectory, knownSnpDFList, alleleDFPathList$filter$newAllelePath, filterLocusConv, 'filter' )
#   }
#   
#   x <- allele.save_call( x, alleleDFPathList$filter$alleleCallPath, 'filter' )
# })

# Processes VCF data and outputs SNP dataframes to pass to allele calling
allele.filter_alignments_to_snp_dfs <- function(currentSample, locusRefList, minDP, kirLocusFeatureNameList, 
                                                filterRefFastaList, filterLocusConv, filterSnpDFList, knownSnpDFList){
  
  if('failed' %in% c(currentSample$geneContent, currentSample$copyNumber)){
    return(currentSample)
  }
  
  cat('\nConverting VCF files to SNP dataframes...')
  
  currentSample[[ 'filterSnpDFList' ]] <- filterSnpDFList
  
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
    
    # VCF failure condition ( most likely due to insufficient depth )
    if( length(vcfDT) == 1 ){
      cat('\t-- VCF READ FAILURE --')
      currentSample[[ 'filterSnpDFList' ]][[currentLocus]]$failure <- TRUE
      next
    }
    
    vcfDT <- allele.filter_VCF_process(vcfDT, currentSample, minDP, bedList, realLocus)
    
    # VCF failure condition ( most likely due to insufficient depth )
    if( length(vcfDT) == 1 ){
      cat('\t-- VCF READ FAILURE --')
      currentSample[[ 'filterSnpDFList' ]][[currentLocus]]$failure <- TRUE
      next
    }
    
    # Fix for KIR3DL1het and KIR3DS1het VCF files, which are treated as KIR3DL1 and KIR3DS1 in allele calling
    if( currentLocus == 'KIR3DL1het' | currentLocus == 'KIR3DS1het' ){
      currentLocus <- filterLocusConv[[currentLocus]]
    }
    
    # ----- EXON processing -----
    cat('\n\tProcessing EXONs')
    exonDT <- vcfDT[ Locus == realLocus ][ featLab %in% exonFeatNameVect ]
    
    # VCF failure condition ( most likely due to insufficient depth )
    if( nrow(exonDT) < 1 ){
      cat('\t-- EXON READ FAILURE --')
      currentSample[[ 'filterSnpDFList' ]][[currentLocus]]$failure <- TRUE
      next
    }
    
    # Split up SNPs from INDELs
    exonDTList <- general.VCF_sep_INDEL(exonDT)
    exonDT <- exonDTList$nodelDT
    exonIndelDT <- exonDTList$indelDT
    
    # VCF failure condition ( most likely due to insufficient depth )
    if( nrow(exonDT) < 1 ){
      cat('\t-- EXON READ FAILURE --')
      currentSample[[ 'filterSnpDFList' ]][[currentLocus]]$failure <- TRUE
      next
    }
    
    exonIndelBoolVect <- !sapply(exonIndelDT$genoVect, function(x) all(x == 1))
    
    # Write exon INDEL table
    if( any(exonIndelBoolVect) ){
      cat('\tINDELSs')
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
    cat('\n\tProcessing INTRONs')
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
      cat('\tINDELs')
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
        
        x.feat <- unique( tstrsplit(names(snp2List),'_',fixed=T)[[1]] )
        if(x.feat == '3UTR'){
          x.1.pos <- tstrsplit(names(snp1List),'_',fixed=T)[[2]]
          x.2.pos <- tstrsplit(names(snp2List),'_',fixed=T)[[2]]
          
          x.featLen <- length( locusRefList[[currentLocus]]$alleleBedList[[1]]$`3UTR`$snpVect ) - 50
          
          # Skip INDEL processing that happens at the end of 3'UTR
          if( any( as.numeric( unique(c(x.1.pos, x.2.pos)) ) > x.featLen ) ) {
            return(NULL)
          }
        }
        
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
  
  '
  KIR2DL3 exon9 is missed by this method because of alterations to the fasta, so the fasta sequence does not match known allele sequence
  '
  
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
  
  if( nrow(vcfDT) < 10 ){
    return('failure')
  }
  
  # Set important features as new data table columns
  vcfDT[, c('SNP1','SNP2') := general.vcfDT_set_SNPs( REF, ALT, genoVect ), by = 1:nrow(vcfDT) ]
  
  vcfDT <- vcfDT[ vcfDT$POS %in% names(bedList[[1]]), ]
  
  if( nrow(vcfDT) < 10 ){
    return('failure')
  }
  
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
    
    # Fix for KIR3DL1het and KIR3DS1het files, which are treated as KIR3DL1 and KIR3DS1 in allele calling
    if( currentLocus == 'KIR3DL1het' | currentLocus == 'KIR3DS1het' ){
      next
    }
    
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
    
    # Save the CSV paths and dataframes to a return list
    output.snpDFList[[currentLocus]] <- list()
    output.snpDFList[[currentLocus]][['exonSNPs']] <- list('csvPath'=exonDFPath)
    output.snpDFList[[currentLocus]][['intronSNPs']] <- list('csvPath'=intronDFPath)
    output.snpDFList[[currentLocus]][['failure']] <- FALSE
    
    resultExists.bool <- file.exists(exonDFPath) & file.exists(intronDFPath)
    
    if( resultExists.bool ){
      cat('\nFound SNP dataframes in',alleleFileDirectory,'~ skipping dataframe generation.\n')
      next
    }
    
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
    #output.snpDFList[[currentLocus]] <- list()
    #output.snpDFList[[currentLocus]][['exonSNPs']] <- list('csvPath'=exonDFPath)
    #output.snpDFList[[currentLocus]][['intronSNPs']] <- list('csvPath'=intronDFPath)
    #output.snpDFList[[currentLocus]][['failure']] <- FALSE
  }
  
  cat('\nFinished. Exon and intron SNP tables located in',alleleFileDirectory)
  return(output.snpDFList)
}

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
  L23DF <- cbind.data.frame(L3DF,L2DF[1,L3AddCols],stringsAsFactors=F)
  
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


# filterLocusConv <- list(
#   'KIR3DL3'='KIR3DL3',
#   'KIR3DL2'='KIR3DL2',
#   'KIR3DS1'='KIR3DS1',
#   'KIR3DS1het'='KIR3DS1',
#   'KIR3DL1het'='KIR3DL1',
#   'KIR3DL1'='KIR3DL1',
#   'KIR2DS35'='KIR2DS5',
#   'KIR2DS4'='KIR2DS4',
#   'KIR2DS3'='KIR2DS3',
#   'KIR2DP1'='KIR2DP1',
#   'KIR2DL5'='KIR2DL5',
#   'KIR2DL4'='KIR2DL4',
#   'KIR2DL2'='KIR2DL2',
#   'KIR2DL23'='KIR2DL3',
#   'KIR2DL3'='KIR2DL3',
#   'KIR2DL1'='KIR2DL1'
# )

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





