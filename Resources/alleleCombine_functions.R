
post.combineGenos <- function( currentSample, resultsDirectory ){
  currentID <- currentSample$name
  
  # ---------- DATA READ IN
  DFList <- post.preGenoDFRead( resultsDirectory, currentID )
  post.iterAlleleDF <- DFList$iter
  post.filterAlleleDF <- DFList$filter
  
  # ---------- Initialize final genotype DF
  # Create the final allele call dataframe and save it
  colNameVect <- c("KIR3DP1","KIR2DS5","KIR2DP1","KIR2DS3","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL1","KIR3DS1",
                   "KIR3DL2","KIR2DS4","KIR2DL1","KIR2DS1","KIR2DL5","KIR2DL23")
  
  finalGenoPath <- file.path(resultsDirectory,'finalAlleleCalls.csv')
  post.finalGenoDF <- post.finalDFRead( finalGenoPath, colNameVect )
  
  # Add current sample ID to the final genotype dataframe
  if(!currentID %in% rownames(post.finalGenoDF)){
    post.finalGenoDF[nrow(post.finalGenoDF)+1,] <- NA
    rownames(post.finalGenoDF)[nrow(post.finalGenoDF)] <- currentID
  }else{
    post.finalGenoDF[currentID,] <- NA
  }
  
  # Setting up KIR2DL23 adjustments (turning 2DL2, 2DL3 calls into 2DL23)
  iter.L23BoolDF <- !is.na( post.iterAlleleDF[currentID,c('KIR2DL2','KIR2DL3','KIR2DL23'),drop=F] )
  post.L23Bool <- iter.L23BoolDF[currentID,'KIR2DL23']
  post.L2Bool <- iter.L23BoolDF[currentID,'KIR2DL2']
  post.L3Bool <- iter.L23BoolDF[currentID,'KIR2DL3']
  
  # ----- Create locus processing loop
  locusVect <- colnames( post.finalGenoDF )
  
  otherLocusVect <- c('KIR3DP1','KIR2DS1','KIR2DS2','KIR2DL2','KIR2DL3','KIR2DL23')
  sharedLocusVect <- setdiff(locusVect, otherLocusVect)
  
  # KIR3DP1 - DONE
  iter.3DP1Result <- post.iterAlleleDF[currentID,'KIR3DP1']
  if( !is.na(iter.3DP1Result) ){
    
    iterNewBool <- grepl('unresolved',post.iterAlleleDF[currentID,'KIR3DP1'],fixed=T)
    
    if(iterNewBool){
      iter.3DP1Result.str <- paste0( post.iterAlleleDF[currentID,'KIR3DP1'],'_ITER' )
    }else{
      # only using ITER workflow results (no FILTER results)
      iterTypeStr <- post.iterAlleleDF[currentID,'KIR3DP1']
      iter.3DP1Result <- post.expandTypeStringToDF( iterTypeStr)
      # Convert the intersection type DF to a genotype string and record to the final genotype DF
      iter.3DP1Result.str <- post.collapseDFToStr( iter.3DP1Result )
    }
    post.finalGenoDF[currentID,'KIR3DP1'] <- iter.3DP1Result.str
  }
  
  # KIR2DS1 - DONE
  iter.2DS1Result <- post.iterAlleleDF[currentID,'KIR2DS1']
  if( !is.na(iter.2DS1Result) ){
    
    iterNewBool <- grepl('unresolved',post.iterAlleleDF[currentID,'KIR2DS1'],fixed=T)
    
    if(iterNewBool){
      iter.2DS1Result.str <- paste0( post.iterAlleleDF[currentID,'KIR2DS1'],'_ITER' )
    }else{
      # only using ITER workflow results (no FILTER results)
      iterTypeStr <- post.iterAlleleDF[currentID,'KIR2DS1']
      iter.2DS1Result <- post.expandTypeStringToDF( iterTypeStr)
      # Convert the intersection type DF to a genotype string and record to the final genotype DF
      iter.2DS1Result.str <- post.collapseDFToStr( iter.2DS1Result )
    }
    post.finalGenoDF[currentID,'KIR2DS1'] <- iter.2DS1Result.str
  }
  
  # KIR2DS2 - DONE
  iter.2DS2Result <- post.iterAlleleDF[currentID,'KIR2DS2']
  if( !is.na(iter.2DS2Result) ){
    
    iterNewBool <- grepl('unresolved',post.iterAlleleDF[currentID,'KIR2DS2'],fixed=T)
    
    if(iterNewBool){
      iter.2DS2Result.str <- paste0( post.iterAlleleDF[currentID,'KIR2DS2'],'_ITER' )
    }else{
      # only using ITER workflow results (no FILTER results)
      iterTypeStr <- post.iterAlleleDF[currentID,'KIR2DS2']
      iter.2DS2Result <- post.expandTypeStringToDF( iterTypeStr)
      # Convert the intersection type DF to a genotype string and record to the final genotype DF
      iter.2DS2Result.str <- post.collapseDFToStr( iter.2DS2Result )
    }
    post.finalGenoDF[currentID,'KIR2DS2'] <- iter.2DS2Result.str
  }
  
  if(post.L23Bool){
    sharedLocusVect <- c(sharedLocusVect,'KIR2DL23')
  }else if(post.L2Bool){
    sharedLocusVect <- c(sharedLocusVect,'KIR2DL2')
  }else if(post.L3Bool){
    sharedLocusVect <- c(sharedLocusVect,'KIR2DL3')
  }
  
  for(currentLocus in sharedLocusVect){
    
    #iterLocusBool <- currentLocus %in% colnames(post.iterAlleleDF)
    #filterLocusBool <- currentLocus %in% colnames(post.filterAlleleDF)
    
    iterNewBool <- grepl('unresolved',post.iterAlleleDF[currentID,currentLocus],fixed=T)
    filterNewBool <- grepl('unresolved',post.filterAlleleDF[currentID,currentLocus],fixed=T)
    
    # ---------- Locus processing -----------
    
    if( currentLocus %in% c('KIR2DL2','KIR2DL3') ){
      
      # prioritize 2DL23 over 2DL2+2DL3
      if( post.L23Bool ){
        currentLocus <- 'KIR2DL23'
      }
      
    }
    
    iterNA.bool <- is.na( post.iterAlleleDF[currentID,currentLocus] ) # flag to mark if NA result
    if( !iterNA.bool ){
      iterTypeStr <- post.iterAlleleDF[currentID,currentLocus]
      post.iterTypeDF <- post.expandTypeStringToDF( iterTypeStr)
    }
    
    filterNA.bool <- is.na( post.filterAlleleDF[currentID,currentLocus] ) # flag to mark if NA result
    if( !filterNA.bool ){
      filterTypeStr <- post.filterAlleleDF[currentID,currentLocus]
      post.filterTypeDF <- post.expandTypeStringToDF( filterTypeStr )
    }
    
    # Condition to catch if one result is NA and the other is defined
    if( iterNA.bool != filterNA.bool ){
      if(iterNA.bool){
        iterNewBool <- T
        filterNewBool <- F
      }else if(filterNA.bool){
        filterNewBool <- T
        iterNewBool <- F
      }
    }
    
    if( iterNewBool & filterNewBool ){
      
      post.iterTypeDF[,'allele1'] <- paste0( post.iterTypeDF[,'allele1'], '_ITER' )
      post.filterTypeDF[,'allele2'] <- paste0( post.filterTypeDF[,'allele2'], '_FILTER' )
      
      post.interTypeDF <- cbind( post.iterTypeDF[,'allele1',drop=F], post.filterTypeDF[,'allele2',drop=F] )
      
    }else if( iterNewBool & !filterNewBool ){
      
      #### Condition if ITER is 'new' allele and FILTER is known
      ## Set the INTERSECTED result to be the FILTER result
      post.interTypeDF <- post.filterTypeDF
      
    }else if( filterNewBool & !iterNewBool ){
      
      #### Condition if FILTER is 'new' allele and ITER is known
      ## Set the INTERSECTED result to be the ITER result
      post.interTypeDF <- post.iterTypeDF
      
    }else{
      
      if( filterNA.bool & iterNA.bool ){
        next
      }
      
      if( currentLocus == 'KIR3DP1' ){
        
        # only using ITER workflow results (no FILTER results)
        post.interTypeDF <- post.iterTypeDF
        
      }
      
      'KIR2DS3 and KIR2DS5 needs to be separated in filter results'
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR2DS3' ){
        
        # intersect results, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
      }
      
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR2DS5' ){
        
        # intersect results, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
      }
      
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR2DL3' | currentLocus == 'KIR2DL2' | currentLocus == 'KIR2DL23' ){
        
        # intersect results first
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
        currentLocus <- 'KIR2DL23' # Make sure the results are saved under KIR2DL23
      }
      
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR2DP1' ){
        
        # Intersect, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
      }
      
      # only ITER
      if( currentLocus == 'KIR2DS2' ){
        
        # only using ITER workflow results (no FILTER results)
        post.interTypeDF <- post.iterTypeDF
        
      }
      
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR2DL4' ){
        
        # intersect results, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
      }
      
      # intersect results, bias towards FILTER
      if( currentLocus == 'KIR3DL3' ){
        
        # intersect results, bias towards FILTER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.filterTypeDF
        }
        
      }
      
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR3DL1' ){
        
        # intersect results, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
      }
      
      # only ITER
      if( currentLocus == 'KIR3DS1' ){
        
        # only using ITER workflow results (no FILTER results)
        post.interTypeDF <- post.iterTypeDF
        
      }
      
      # intersect results, bias towards ITER
      if( currentLocus == 'KIR3DL2' ){
        
        # intersect results, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.iterTypeDF
        }
        
      }
      
      # only FILTER
      if( currentLocus == 'KIR2DS4' ){
        
        # only using FILTER workflow results
        post.interTypeDF <- post.filterTypeDF
        
      }
      
      # intersect results, bias towards FILTER
      if( currentLocus == 'KIR2DL1' ){
        
        # intersect results, bias towards ITER
        post.interTypeDF <- post.intersectITERwithFILTER( post.iterTypeDF, post.filterTypeDF )
        
        # Bias for ITER results if no intersection
        if( nrow(post.interTypeDF) == 0 ){
          post.interTypeDF <- post.filterTypeDF
        }
        
      }
      
      # only ITER
      if( currentLocus == 'KIR2DS1' ){
        
        # only using ITER workflow results (no FILTER results)
        post.interTypeDF <- post.iterTypeDF
        
      }
      
      # only FILTER
      if( currentLocus == 'KIR2DL5' ){
        
        # only using FILTER workflow results
        post.interTypeDF <- post.filterTypeDF
        
      }
      
    }
    
    if( !(filterNA.bool & iterNA.bool) ){
      
      # Convert the intersection type DF to a genotype string and record to the final genotype DF
      post.finalGenoDF[currentID,currentLocus] <- post.collapseDFToStr( post.interTypeDF )
      
    }
  }
  
  write.csv(post.finalGenoDF, finalGenoPath)
  return( currentSample )
}

# ----- support functions -----
post.expandTypeStringToDF <- function( currentLocusType ){
  
  ## Split ambiguous typings
  currentLocusTypeVect <- strsplit(currentLocusType,' ',fixed=T)[[1]]
  
  ## Expanding typings that need to be expanded
  if(all(!grepl('+',currentLocusTypeVect,fixed=T)) & length(currentLocusTypeVect)>1){
    
    combMat <- combinations(length(currentLocusTypeVect),2,currentLocusTypeVect,repeats.allowed = T)
    currentLocusTypeVect <- apply(combMat, 1, paste0, collapse='+')
  }
  
  ## Initialize a dataframe for storing allele typings, 1 row for each ambiguity
  alleleTypeFrame <- data.frame(matrix('',nrow=length(currentLocusTypeVect),ncol=2),
                                stringsAsFactors = F)
  colnames(alleleTypeFrame) <- c('allele1','allele2')
  
  i <- 0
  ## For each ambiguity, record the typings in the dataframe
  for(singleType in currentLocusTypeVect){
    i <- i+1
    con2 <- grepl('+',singleType,fixed=T)
    
    if(con2){
      hetTypeVect <- unlist(tstrsplit(singleType,'+',fixed=T))
      alleleTypeFrame[as.character(i),] <- hetTypeVect
    }else{
      alleleTypeFrame[as.character(i),] <- singleType
    }
  }
  
  return(alleleTypeFrame)
}

post.intersectITERwithFILTER <- function( iterTypeDF, filterTypeDF ){
  
  # Intersecting ITER results with FILTER results
  interBoolVectFW <- apply(iterTypeDF, 1, function(x){
    return(any(apply(filterTypeDF, 1, function(y){
      return(x[2] %in% y[2] & x[1] %in% y[1])
    })))
  })
  
  interBoolVectRW <- apply(iterTypeDF, 1, function(x){
    return(any(apply(filterTypeDF, 1, function(y){
      return(x[1] %in% y[2] & x[2] %in% y[1])
    })))
  })
  
  interBoolVect <- interBoolVectFW | interBoolVectRW
  
  interResultDF <- iterTypeDF[interBoolVect,]
  
  return(interResultDF)
}

post.collapseDFToStr <- function( alleleDF ){
  return( paste0( apply( alleleDF, 1, paste0, collapse='+' ), collapse=' ' ) )
}

post.filter.addKIR2DS5 <- function( typeDF, currentID ){
  
    
  index2DS5Int <- which(grepl('KIR2DS5',typeDF[currentID,,drop=F],fixed=T))
  
  if( length(index2DS5Int) == 0 ){
    newBool <- grepl('unresolved',typeDF[currentID,'KIR2DS35',drop=F],fixed=T)
    
    if( newBool ){
      type2DS5Str <- 'unresolved'
    }else{
      type2DS5Str <- ''
    }
    
  }else{
    type2DS5Str <- typeDF[currentID,'KIR2DS35'][[1]]
  }
  
  vectType2DS5 <- tstrsplit(type2DS5Str,' ',fixed=T)
  listType2DS5 <- list()
  j <- 1
  for(sub2DS5Str in vectType2DS5){
    
    ## KIR2DS3 split
    if(any(grepl('KIR2DS3', sub2DS5Str, fixed=T))){
      splitVect <- strsplit(sub2DS5Str,'+',fixed=T)[[1]]
      sub2DS5Str <- splitVect[grepl('KIR2DS5',splitVect,fixed=T)]
      listType2DS5[[j]] <- sub2DS5Str
      j <- j+1
    }else{
      listType2DS5[[j]] <- sub2DS5Str
      j <- j+1
    }
    
  }
  
  typeStr <- paste0(listType2DS5, collapse=' ')
  
  if( typeStr == '' ){
    typeStr <- NA
  }
  
  typeDF[currentID,'KIR2DS5'] <- typeStr
  return(typeDF)
}

post.filter.addKIR2DS3 <- function( typeDF, currentID ){
  
  ## KIR2DS3 isolation in PING typings
    
  index2DS3IntVect <- which(grepl('KIR2DS3',typeDF[currentID,,drop=F],fixed=T))
  
  if( length(index2DS3IntVect) == 0 ){
    
    newBool <- any(grepl('unresolved',typeDF[currentID,c('KIR2DS3','KIR2DS35'),drop=F],fixed=T))
    
    if( newBool ){
      type2DS3Str <- 'unresolved'
    }else{
      type2DS3Str <- ''
    }
  }else{
    type2DS3Str <- typeDF[currentID,index2DS3IntVect][[1]]
  }
  
  vectType2DS3 <- tstrsplit(type2DS3Str,' ',fixed=T)
  listType2DS3 <- list()
  j <- 1
  for(sub2DS3Str in vectType2DS3){
    
    ## KIR2DS3 split
    if(any(grepl('KIR2DS5', sub2DS3Str, fixed=T))){
      splitVect <- strsplit(sub2DS3Str,'+',fixed=T)[[1]]
      sub2DS3Str <- splitVect[grepl('KIR2DS3',splitVect,fixed=T)]
      listType2DS3[[j]] <- sub2DS3Str
      j <- j+1
    }else{
      listType2DS3[[j]] <- sub2DS3Str
      j <- j+1
    }
    
  }
  
  type2DS3Str <- paste0(listType2DS3, collapse=' ')
  
    #if(length(index2DS3IntVect) > 1){
    #  cat('\n',type2DS3Str)
    #}
    
  if( type2DS3Str == '' ){
    type2DS3Str <- NA
  }
    
  typeDF[currentID,'KIR2DS3'] <- type2DS3Str
    
  return(typeDF)
}

post.preGenoDFRead <- function( resultsDirectory, currentID ){
  
  post.iterAlleleDF <- read.csv(file.path(resultsDirectory,'iterAlleleCalls.csv'),
                                stringsAsFactors = F,check.names=F,row.names=1)
  
  post.filterAlleleDF <- read.csv(file.path(resultsDirectory,'filterAlleleCalls.csv'),
                                  stringsAsFactors=F, check.names=F, row.names=1)
  
  post.iterAlleleDF <- post.iterAlleleDF[currentID,,drop=F]
  post.filterAlleleDF <- post.filterAlleleDF[currentID,,drop=F]
  
  post.filterAlleleDF$KIR2DS5 <- NA # initialize KIR2DS5 allele call slot for FILTER ganotypes
  post.filterAlleleDF <- post.filter.addKIR2DS5( post.filterAlleleDF, currentID ) # separate KIR2DS5 from KIR2DS35
  post.filterAlleleDF <- post.filter.addKIR2DS3( post.filterAlleleDF, currentID ) # separate KIR2DS3 from KIR2DS35
  
  return( list('iter'=post.iterAlleleDF, 'filter'=post.filterAlleleDF) )
}

post.finalDFRead <- function( finalGenoPath, colNameVect ){
  
  if( !file.exists(finalGenoPath) ){
    post.finalGenoDF <- data.frame(matrix(NA,nrow=0,ncol=length(colNameVect)), check.names = F, stringsAsFactors = F)
    colnames(post.finalGenoDF) <- colNameVect
    write.csv(post.finalGenoDF, finalGenoPath)
  }
  
  # Read final genotype dataframe
  post.finalGenoDF <- read.csv(finalGenoPath, 
                               stringsAsFactors = F, 
                               check.names = F, 
                               row.names = 1)
  
  return(post.finalGenoDF)
  
}
