new.initLocusRef.read_snp_df <- function(locusRefList, snpDFDirectory, kirLocusVect = kir.locus.vect){
  
  knownSnpDFList <- list()
  cat('\n\nReading in:')
  for(locus in kirLocusVect){
    cat('',locus)
    alleleSnpDF <- read.csv( normalizePath(
      file.path(snpDFDirectory,
                paste0(locus, '_alleleSNPs.csv')
      )), check.names=F, row.names=1, stringsAsFactors = F,colClasses = 'character')
    
    knownSnpDFList[[locus]] <- alleleSnpDF
  }
  
  cat('\n\nFilling invariant positions:')
  for( locus in kirLocusVect ){
    cat('',locus)
    snpDF <- knownSnpDFList[[locus]]
    numUniqueNuc.byPosition.list <- apply( snpDF, 2, num_unique_nuc )
    invariantColVect <- names( which( numUniqueNuc.byPosition.list == 1 ) )
    
    invariantNucVect <- unlist( apply( snpDF[,invariantColVect], 2, function(x) unique( x[is_nuc(x)] ) ) )
    
    for(i in 1:nrow(snpDF)){
      snpDF[i,invariantColVect] <- invariantNucVect
    }
    
    knownSnpDFList[[locus]] <- snpDF
  }
  
  return(knownSnpDFList)
}
new.initLocusRef.extend_5UTR <- function( filled.snpDFList, UTRextList ){
  
  cat('\n\nExtending 5UTR:')
  out.list <- list()
  for( currentLocus in names(filled.snpDFList) ){
    cat('',currentLocus)
    snpDF <- filled.snpDFList[[currentLocus]]
    utr5.ext.str <- UTRextList[[paste0(currentLocus,'_5UTR')]]
    utr5.ext.str.length <- nchar(utr5.ext.str)
    
    utr5.colVect <- grep('5UTR', colnames(snpDF), value=T)
    utr5.indexVect <- grep('5UTR', colnames(snpDF))
    utr5.length <- length(utr5.colVect)
    
    utr5.addition.length <- utr5.ext.str.length - utr5.length
    
    utr5.addition.str <- substr(utr5.ext.str, 1, utr5.addition.length)
    
    
    # Replace 5UTR naming by offset (to account for added positions)
    utr5.replacement.colVect <- unlist( sapply(utr5.colVect, function(x) {
      xVect <- strsplit(x,'_',fixed=T)[[1]]
      xPos <- as.integer(xVect[2]) + utr5.addition.length
     return( paste0(xVect[1],'_',xPos) )
      }) )
    
    colnames(snpDF)[utr5.indexVect] <- utr5.replacement.colVect
    
    utr5.add.df <- as.data.frame( matrix('',nrow=nrow(snpDF),ncol=utr5.addition.length), check.names=F, stringsAsFactors = F)
    rownames(utr5.add.df) <- rownames(snpDF)
    colnames(utr5.add.df) <- paste0('5UTR_',1:utr5.addition.length)
    
    for(i in 1:nrow(utr5.add.df)){
      utr5.add.df[i,] <- strsplit( utr5.addition.str, '')[[1]]
    }
    
    snpDF <- cbind(utr5.add.df, snpDF)
    
    out.list[[currentLocus]] <- snpDF
  }
  
  return(out.list)
}
new.initLocusRef.extend_3UTR <- function( filled.snpDFList, UTRextList ){
  
  cat('\n\nExtending 3UTR:')
  out.list <- list()
  for( currentLocus in names(filled.snpDFList) ){
    cat('',currentLocus)
    snpDF <- filled.snpDFList[[currentLocus]]
    utr3.ext.str <- UTRextList[[paste0(currentLocus,'_3UTR')]]
    
    utr3.ext.str.length <- nchar(utr3.ext.str)
    
    utr3.colVect <- grep('3UTR',colnames(snpDF),value=T)
    utr3.indexVect <- grep('3UTR',colnames(snpDF))
    utr3.length <- length(utr3.colVect)
    
    utr3.addition.length <- utr3.ext.str.length - utr3.length
    
    utr3.addition.str <- substr(utr3.ext.str, (utr3.length+1), utr3.ext.str.length)
    
    utr3.add.df <- as.data.frame( matrix('',nrow=nrow(snpDF),ncol=utr3.addition.length), check.names=F, stringsAsFactors = F)
    rownames(utr3.add.df) <- rownames(snpDF)
    colnames(utr3.add.df) <- paste0('3UTR_',(utr3.length+1):utr3.ext.str.length)
    
    for(i in 1:nrow(utr3.add.df)){
      utr3.add.df[i,] <- strsplit( utr3.addition.str, '')[[1]]
    }
    
    snpDF <- cbind(snpDF, utr3.add.df)
    
    out.list[[currentLocus]] <- snpDF
  }
  
  return(out.list)
}
new.initLocusRef.snpDFtoLocusRefAlleleSeq <- function( filled.snpDFList, locusRefList ){
  cat('\n\nAdding allele sequences to locus reference object:')
  for(locusRef in locusRefList){
    cat('',locusRef$name)
    locusRef$alleleStrList <- as.list( apply(filled.snpDFList[[locusRef$name]],1,function(x) paste0(x,collapse='') ) )
  }
  
  return(locusRefList)
}
new.initLocusRef.snpDFtoLocusRefBed <- function( filled.snpDFList, locusRefList, kirLocusFeatureNameList ){
  cat('\n\nAdding feature names and coordinates to locus reference object:')
  for(locusRef in locusRefList){
    cat('',locusRef$name)
    
    snpDF <- filled.snpDFList[[locusRef$name]]
    
    ## pull out the locus feature names from the SNP DF and predefined
    locusFeatVect <- unique( tstrsplit( colnames(snpDF), '_', fixed=T)[[1]] )
    featureNameVect <- kirLocusFeatureNameList[[locusRef$name]]
    if( length(locusFeatVect) != length(featureNameVect) ){
      stop('Mismatched features for',locusRef$name)
    }
    
    locusBedList <- list()
    for( alleleNameStr in names(locusRef$alleleStrList) ){
      alleleStr <- locusRef$alleleStrList[[alleleNameStr]]
      
      ## Initialize list for storing BED information
      locusBedList[[alleleNameStr]] <- list()
      
      totalDelOffset <- 0
      ## Process each gene feature for this allele
      for(featNameStr in featureNameVect){
        
        featIndexVect <- grep(paste0(featNameStr,'_'),colnames(snpDF),fixed=T)
        featColVect <- grep(paste0(featNameStr,'_'),colnames(snpDF),fixed=T,value=T)
        
        featSeqStr <- paste0( snpDF[alleleNameStr,featColVect], collapse='' )
        snpVect <- unlist( snpDF[alleleNameStr,featColVect] )
        preBoundaryInt <- as.integer( featIndexVect[1]-1 ) - totalDelOffset
        boundaryInt <- as.integer( featIndexVect[length(featIndexVect)] ) - totalDelOffset
        
        delCount <- str_count(featSeqStr,pattern = fixed('.'))
        totalDelOffset <- delCount + totalDelOffset
        
        boundaryInt <- boundaryInt - delCount
        
        delIndex <- grep('.',unlist(strsplit(featSeqStr,'')),fixed=T)
        noDelFeatSeqStr <- gsub('.','',featSeqStr,fixed=T)
        
        #### BED coordinates should be preBoundaryInt : boundaryInt
        ## Save gene feature coordinates to list
        locusBedList[[alleleNameStr]][[featNameStr]] <- list(alleleName=alleleNameStr,
                                                             startPos=preBoundaryInt,
                                                             endPos=boundaryInt,
                                                             featName=featNameStr,
                                                             featSeq=noDelFeatSeqStr,
                                                             featDelIndex=delIndex,
                                                             snpVect=snpVect)
        
        preBoundaryInt <- boundaryInt
      }
    }
    
    locusRef$alleleBedList <- locusBedList
  }
  
  return(locusRefList)
}

update_gc_reference_csv <- function(referenceCSVPath,
                                    filled.snpDFList,
                                    outputCSVPath = referenceCSVPath,
                                    logPath = NULL,
                                    strict = TRUE) {

  refDF <- read.csv(
    referenceCSVPath,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"iter_1" %in% colnames(refDF))
    stop("gc_allele_reference.csv must contain iter_1...iter_5 columns")

  logLines <- character()

  for (locus in rownames(refDF)) {

    if (!locus %in% names(filled.snpDFList)) {
      msg <- paste0("Locus not found in filled.snpDFList: ", locus)

      if (strict) stop(msg)
      warning(msg)
      next
    }

    validAlleles <- rownames(filled.snpDFList[[locus]])

    if (length(validAlleles) == 0) {
      msg <- paste0("No alleles found for locus: ", locus)

      if (strict) stop(msg)
      warning(msg)
      next
    }

    for (iterCol in colnames(refDF)) {

      refAllele <- refDF[locus, iterCol]

      if (is.na(refAllele) || !nzchar(refAllele))
        next

      ## Exact match
      if (refAllele %in% validAlleles)
        next

      ## Try to find the successor allele
      replacement <- find_successor_allele(
        refAllele = refAllele,
        validAlleles = validAlleles
      )

      if (length(replacement) == 0) {

        msg <- paste0(
          "No replacement found for ",
          locus,
          " / ",
          iterCol,
          ": ",
          refAllele
        )

        if (strict) {
          stop(msg)
        } else {
          warning(msg)
          next
        }
      }

      if (replacement != refAllele) {

        message(
          sprintf(
            "[gc_allele_reference] %s %s: %s -> %s",
            locus,
            iterCol,
            refAllele,
            replacement
          )
        )

        logLines <- c(
          logLines,
          sprintf(
            "%s,%s,%s,%s",
            locus,
            iterCol,
            refAllele,
            replacement
          )
        )

        refDF[locus, iterCol] <- replacement
      }
    }
  }

  write.csv(
    refDF,
    file = outputCSVPath,
    quote = FALSE
  )

  if (!is.null(logPath)) {
    writeLines(logLines, logPath)
  }

  invisible(refDF)
}

find_successor_allele <- function(refAllele, validAlleles) {

  escape_regex <- function(x) {
    gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", x)
  }
  ## 1. Exact match
  if (refAllele %in% validAlleles)
    return(refAllele)

  ## 2. Direct prefix match
  matches <- grep(
    paste0("^", escape_regex(refAllele)),
    validAlleles,
    value = TRUE
  )

  if (length(matches) > 0)
    return(sort(matches)[1])

  ## 3. Try inserting 01 before an expression suffix
  ##
  ##   *049N  -> *04901N
  ##   *027L  -> *02701L
  ##
  expanded <- sub(
    "([0-9]+)([NLSQCA])$",
    "\\101\\2",
    refAllele
  )

  if (expanded != refAllele) {

    matches <- grep(
      paste0("^", escape_regex(expanded)),
      validAlleles,
      value = TRUE
    )

    if (length(matches) > 0)
      return(sort(matches)[1])
  }

  ## 4. Try appending 0101
  ##
  ##   *004 -> *0040101
  ##
  expanded <- paste0(refAllele, "0101")

  matches <- grep(
    paste0("^", escape_regex(expanded)),
    validAlleles,
    value = TRUE
  )

  if (length(matches) > 0)
    return(sort(matches)[1])

  ## 5. Fall back to any prefix beginning with the allele number
  matches <- grep(
    paste0("^", escape_regex(refAllele)),
    validAlleles,
    value = TRUE
  )

  if (length(matches) > 0)
    return(sort(matches)[1])

  return(character(0))
}

kir.locus.vect <- c("KIR3DP1","KIR2DS5","KIR2DL3","KIR2DP1","KIR2DS3","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL1","KIR3DS1","KIR2DL2","KIR3DL2","KIR2DS4","KIR2DL1","KIR2DS1","KIR2DL5")

