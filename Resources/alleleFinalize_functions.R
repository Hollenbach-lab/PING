
pingFinalize.combineL23 <- function( iterCall.df ){
  
  L23.call.list <- sapply( rownames( iterCall.df ), function(x){
    L2.call <- iterCall.df[x,'KIR2DL2']
    L3.call <- iterCall.df[x,'KIR2DL3']
    
    L2.call.bool <- !is.na(L2.call) & L2.call != 'failed'
    L2.call.unresolved <- grepl('$',L2.call,fixed=T)
    
    L3.call.bool <- !is.na(L3.call) & L3.call != 'failed'
    L3.call.unresolved <- grepl('$',L3.call,fixed=T)
    
    if( L2.call.unresolved ){
      L2.call <- paste0('KIR2DL2*unresolved')
    }
    if( L3.call.unresolved ){
      L3.call <- paste0('KIR2DL3*unresolved')
    }
    
    if( L2.call.bool & L3.call.bool ){
      L2.call.vect <- strsplit(L2.call,' ',fixed=T)[[1]]
      L3.call.vect <- strsplit(L3.call,' ',fixed=T)[[1]]
      L23.call.mat <- expand.grid(L2.call.vect,L3.call.vect)
      L23.call.vect <- apply( L23.call.mat, 1, paste0, collapse='+' )
    }else if( L2.call.bool & !L3.call.bool ){
      
      L2.hetBool <- grepl('+',L2.call,fixed=T)
      
      if(!L2.hetBool){
        L2.callVect <- strsplit(L2.call,' ',fixed=T)[[1]]
        
        '
        Combinations could cause recursion limit problems for highly ambiguous calls
        '
        L23.call.mat <- combinations(length(L2.callVect),2,L2.callVect,repeats.allowed = T)
        L23.call.vect <- apply( L23.call.mat, 1, paste0, collapse='+')
      }else{
        L23.call.vect <- unlist( tstrsplit(L2.call,' ',fixed=T)[[1]] )
      }
      
    }else if( !L2.call.bool & L3.call.bool ){
      L3.hetBool <- grepl('+',L3.call,fixed=T)
      
      if(!L3.hetBool){
        L3.callVect <- strsplit(L3.call,' ',fixed=T)[[1]]
        L23.call.mat <- combinations(length(L3.callVect),2,L3.callVect,repeats.allowed = T)
        L23.call.vect <- apply( L23.call.mat, 1, paste0, collapse='+')
      }else{
        L23.call.vect <-  unlist( tstrsplit(L3.call,' ',fixed=T)[[1]] )
      }
    }else{
      L23.call.vect <- 'failed'
    }
    
    if( length(L23.call.vect) > 64 ){
      L23.call.str <- 'failed'
    }else{
      L23.call.str <- paste0(L23.call.vect, collapse=' ')
    }
    return(L23.call.str)
  })
  
  iterCall.df$KIR2DL23 <- NA
  iterCall.df[names(L23.call.list),'KIR2DL23'] <- unlist( L23.call.list )
  
  return(iterCall.df)
}

pingFinalize.combineS35 <- function( iterCall.df, copyCall.df ){
  S35.call.list <- sapply( rownames( iterCall.df ), function(x){
    #cat('\n',x)
    S3.call <- iterCall.df[x,'KIR2DS3']
    S5.call <- iterCall.df[x,'KIR2DS5']
    
    S3.copy <- as.integer( copyCall.df[x,'KIR2DS3'] )
    S5.copy <- as.integer( copyCall.df[x,'KIR2DS5'] )
    
    S3.call.bool <- !is.na(S3.call) & S3.call != 'failed'
    S3.call.unresolved <- grepl('$',S3.call,fixed=T)
    
    S5.call.bool <- !is.na(S5.call) & S5.call != 'failed'
    S5.call.unresolved <- grepl('$',S5.call,fixed=T)
    
    if( S3.call.unresolved ){
      S3.call <- paste0('KIR2DS3*unresolved')
    }
    if( S5.call.unresolved ){
      S5.call <- paste0('KIR2DS5*unresolved')
    }
    
    if( S3.call.bool ){
      S3.het.bool <- grepl('+',S3.call,fixed=T)
    }
    
    if( S5.call.bool ){
      S5.het.bool <- grepl('+',S5.call,fixed=T)
    }
    
    if( S3.call.bool & S5.call.bool ){
      
      if( !S3.het.bool ){
        S3.allele1.vect <- strsplit(S3.call,' ',fixed=T)[[1]]
        
        if( S3.copy > 1 ){
          S3.call.mat <- combinations(length(S3.allele1.vect), S3.copy, S3.allele1.vect,repeats.allowed = T)
          S3.call.vect <- apply( S3.call.mat, 1, paste0, collapse='+')
        }else{
          S3.call.vect <- S3.allele1.vect
        }
      }else{
        S3.hetCall.vect <- tstrsplit(S3.call,' ',fixed=T)
        #S3.hetCall.list <- tstrsplit(S3.hetCall.vect,'+',fixed=T)
        #S3.allele1.vect <- S3.hetCall.list[[1]]
        #S3.allele2.vect <- S3.hetCall.list[[2]]
        
        #S3.call.vect <- apply( unique( expand.grid(S3.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        S3.call.vect <- unlist( S3.hetCall.vect )
      }
      
      if( !S5.het.bool ){
        S5.allele1.vect <- strsplit(S5.call,' ',fixed=T)[[1]]
        if( S5.copy > 1 ){
          S5.call.mat <- combinations(length(S5.allele1.vect), S5.copy, S5.allele1.vect,repeats.allowed = T)
          S5.call.vect <- apply( S5.call.mat, 1, paste0, collapse='+')
        }else{
          S5.call.vect <- S5.allele1.vect
        }
      }else{
        S5.hetCall.vect <- tstrsplit(S5.call,' ',fixed=T)
        #S5.hetCall.list <- tstrsplit(S5.hetCall.vect,'+',fixed=T)
        #S5.allele1.vect <- S5.hetCall.list[[1]]
        #S5.allele2.vect <- S5.hetCall.list[[2]]
        
        #S5.call.vect <- apply( unique( expand.grid(S5.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        S5.call.vect <- unlist( S5.hetCall.vect )
      }
      
      S35.call.mat <- expand.grid(S3.call.vect,S5.call.vect)
      S35.call.vect <- apply( S35.call.mat, 1, paste0, collapse='+')
    }else if( S3.call.bool & !S5.call.bool ){
      if( !S3.het.bool ){
        S3.allele1.vect <- strsplit(S3.call,' ',fixed=T)[[1]]
        
        if( S3.copy > 1 ){
          S3.call.mat <- combinations(length(S3.allele1.vect), S3.copy, S3.allele1.vect,repeats.allowed = T)
          S35.call.vect <- apply( S3.call.mat, 1, paste0, collapse='+')
        }else{
          S35.call.vect <- apply( expand.grid( S3.allele1.vect, 'KIR2DS35*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S3.hetCall.vect <- tstrsplit(S3.call,' ',fixed=T)
        #S3.hetCall.list <- tstrsplit(S3.hetCall.vect,'+',fixed=T)
        
        #S35.call.vect <- apply( unique( expand.grid(S3.hetCall.list) ), 1, function(x) paste0(x,collapse='+') )
        S35.call.vect <- unlist( S3.hetCall.vect )
      }
    }else if(!S3.call.bool & S5.call.bool){
      if( !S5.het.bool ){
        S5.allele1.vect <- strsplit(S5.call,' ',fixed=T)[[1]]
        
        if( S5.copy > 1 ){
          S5.call.mat <- combinations(length(S5.allele1.vect), S5.copy, S5.allele1.vect,repeats.allowed = T)
          S35.call.vect <- apply( S5.call.mat, 1, paste0, collapse='+')
        }else{
          S35.call.vect <- apply( expand.grid( S5.allele1.vect, 'KIR2DS35*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S5.hetCall.vect <- tstrsplit(S5.call,' ',fixed=T)
        #S5.hetCall.list <- tstrsplit(S5.hetCall.vect,'+',fixed=T)
        
        #S35.call.vect <- apply( unique( expand.grid(S5.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        #S35.call.vect <- paste(S35.call.vect)
        S35.call.vect <- unlist( S5.hetCall.vect )
      }
    }else if(S3.copy == 0 & S5.copy == 0) {
      S35.call.vect <- 'KIR2DS35*null+KIR2DS35*null'
    }else{
      S35.call.vect <- 'failed'
    }
    
    if( length(S35.call.vect) > 64 ){
      S35.call.str <- 'failed'
    }else{
      S35.call.str <- paste0(S35.call.vect, collapse=' ')
    }
  })
  
  iterCall.df$KIR2DS35 <- NA
  iterCall.df[names(S35.call.list),'KIR2DS35'] <- unlist( S35.call.list )
  
  return(iterCall.df)
}

pingFinalize.combineL1S1 <- function( iterCall.df, copyCall.df ){
  L1S1.call.list <- sapply( rownames( iterCall.df ), function(x){
    S1.call <- iterCall.df[x,'KIR3DS1']
    L1.call <- iterCall.df[x,'KIR3DL1']
    
    S1.copy <- as.integer( copyCall.df[x,'KIR3DS1'] )
    L1.copy <- as.integer( copyCall.df[x,'KIR3DL1'] )
    
    S1.call.bool <- !is.na(S1.call) & S1.call != 'failed'
    S1.call.unresolved <- grepl('$',S1.call,fixed=T)
    
    L1.call.bool <- !is.na(L1.call) & L1.call != 'failed'
    L1.call.unresolved <- grepl('$',L1.call,fixed=T)
    
    if( S1.call.unresolved ){
      S1.call <- paste0('KIR3DS1*unresolved')
    }
    if( L1.call.unresolved ){
      L1.call <- paste0('KIR3DL1*unresolved')
    }
    
    if( S1.call.bool ){
      S1.het.bool <- grepl('+',S1.call,fixed=T)
    }
    
    if( L1.call.bool ){
      L1.het.bool <- grepl('+',L1.call,fixed=T)
    }
    
    if( S1.call.bool & L1.call.bool ){
      
      if( !S1.het.bool ){
        S1.allele1.vect <- strsplit(S1.call,' ',fixed=T)[[1]]
        
        if( S1.copy > 1 ){
          S1.call.mat <- combinations(length(S1.allele1.vect), S1.copy, S1.allele1.vect,repeats.allowed = T)
          S1.call.vect <- apply( S1.call.mat, 1, paste0, collapse='+')
        }else{
          S1.call.vect <- S1.allele1.vect
        }
      }else{
        S1.hetCall.vect <- tstrsplit(S1.call,' ',fixed=T)
        S1.call.vect <- unlist( S1.hetCall.vect )
      }
      
      if( !L1.het.bool ){
        L1.allele1.vect <- strsplit(L1.call,' ',fixed=T)[[1]]
        if( L1.copy > 1 ){
          L1.call.mat <- combinations(length(L1.allele1.vect), L1.copy, L1.allele1.vect,repeats.allowed = T)
          L1.call.vect <- apply( L1.call.mat, 1, paste0, collapse='+')
        }else{
          L1.call.vect <- L1.allele1.vect
        }
      }else{
        L1.hetCall.vect <- tstrsplit(L1.call,' ',fixed=T)
        L1.call.vect <- unlist( L1.hetCall.vect )
      }
      
      L1S1.call.mat <- expand.grid(S1.call.vect,L1.call.vect)
      L1S1.call.vect <- apply( L1S1.call.mat, 1, paste0, collapse='+')
    }else if( S1.call.bool & !L1.call.bool ){
      if( !S1.het.bool ){
        S1.allele1.vect <- strsplit(S1.call,' ',fixed=T)[[1]]
        
        if( S1.copy > 1 ){
          S1.call.mat <- combinations(length(S1.allele1.vect), S1.copy, S1.allele1.vect,repeats.allowed = T)
          L1S1.call.vect <- apply( S1.call.mat, 1, paste0, collapse='+')
        }else{
          L1S1.call.vect <- apply( expand.grid( S1.allele1.vect, 'KIR3DL1S1*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S1.hetCall.vect <- tstrsplit(S1.call,' ',fixed=T)
        L1S1.call.vect <- unlist( S1.hetCall.vect )
      }
    }else if(!S1.call.bool & L1.call.bool){
      if( !L1.het.bool ){
        L1.allele1.vect <- strsplit(L1.call,' ',fixed=T)[[1]]
        
        if( L1.copy > 1 ){
          L1.call.mat <- combinations(length(L1.allele1.vect), L1.copy, L1.allele1.vect,repeats.allowed = T)
          L1S1.call.vect <- apply( L1.call.mat, 1, paste0, collapse='+')
        }else{
          L1S1.call.vect <- apply( expand.grid( L1.allele1.vect, 'KIR3DL1S1*null' ), 1, paste0, collapse='+' )
        }
      }else{
        L1.hetCall.vect <- tstrsplit(L1.call,' ',fixed=T)
        L1S1.call.vect <- unlist( L1.hetCall.vect )
      }
    }else if(S1.copy == 0 & L1.copy == 0) {
      L1S1.call.vect <- 'KIR3DL1S1*null+KIR3DL1S1*null'
    }else{
      L1S1.call.vect <- 'failed'
    }
    
    if( length(L1S1.call.vect) > 64 ){
      L1S1.call.str <- 'failed'
    }else{
      L1S1.call.str <- paste0(L1S1.call.vect, collapse=' ')
    }
    
    return(L1S1.call.str)
    
  })
  
  iterCall.df$KIR3DL1S1 <- NA
  iterCall.df[names(L1S1.call.list),'KIR3DL1S1'] <- unlist( L1S1.call.list )
  
  return(iterCall.df)
}

pingFinalize.otherLoci <- function( iterCall.df, copyCall.df ){
  
  lociVect <- colnames(iterCall.df)
  lociVect <- lociVect[ !lociVect %in% c('KIR3DL1S1','KIR2DL23','KIR2DS35','KIR3DL1','KIR3DS1','KIR2DS3','KIR2DS5','KIR2DL2','KIR2DL3')]
  
  for( currentLocus in lociVect ){
    locus.call.list <- sapply( rownames( iterCall.df ), function(x){
      locus.call <- iterCall.df[x,currentLocus]
      locus.copy <- as.integer( copyCall.df[x,currentLocus] )
      
      locus.call.bool <- !is.na(locus.call) & locus.call != 'failed'
      locus.call.unresolved <- grepl('$',locus.call,fixed=T)
      
      if( locus.call.unresolved ){
        locus.call <- paste0(currentLocus,'*unresolved')
      }
      
      if( locus.call.bool ){
        locus.het.bool <- grepl('+',locus.call,fixed=T)
      }
      
      if( locus.call.bool ){
        if( !locus.het.bool ){
          locus.allele1.vect <- strsplit(locus.call,' ',fixed=T)[[1]]
          
          if( locus.copy > 1 ){
            locus.call.mat <- combinations(length(locus.allele1.vect), locus.copy, locus.allele1.vect,repeats.allowed = T)
            locus.call.vect <- apply( locus.call.mat, 1, paste0, collapse='+')
          }else{
            locus.call.vect <- apply( expand.grid( locus.allele1.vect, paste0(currentLocus,'*null') ), 1, paste0, collapse='+' )
            #locus.call.vect <- locus.allele1.vect
          }
        }else{
          locus.hetCall.vect <- tstrsplit(locus.call,' ',fixed=T)
          #locus.hetCall.list <- tstrsplit(locus.hetCall.vect,'+',fixed=T)
          #locus.call.vect <- apply( unique( expand.grid(locus.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
          
          locus.call.vect <- unlist( locus.hetCall.vect )
        }
      }else if(locus.copy == 0){
        locus.call.vect <- paste0(currentLocus,'*null+',currentLocus,'*null')
      }else{
        locus.call.vect <- 'failed'
      }
      
      if( length(locus.call.vect) > 64 ){
        locus.call.str <- 'failed'
      }else{
        locus.call.str <- paste0(locus.call.vect, collapse=' ')
      }
      
      return(locus.call.str)
    })
    iterCall.df[names(locus.call.list),currentLocus] <- unlist( locus.call.list )
  }
  
  return(iterCall.df)
}

# PING2 allele finalization
pingFinalize.format_calls <- function( resultsDirectory ){
  iterCall.df <- read.csv(file.path(resultsDirectory,'iterAlleleCalls.csv'),
                          stringsAsFactors=F,check.names=F,row.names=1)
  copyCall.df <- read.csv(file.path(resultsDirectory,'manualCopyNumberFrame.csv'),
                          stringsAsFactors=F,check.names=F,row.names=1)
  
  goodSampleID.vect <- rownames(iterCall.df)[ !apply( iterCall.df, 1, function(x) all( is.na(x) ) ) ]
  iterCall.df <- iterCall.df[goodSampleID.vect,]
  copyCall.df <- copyCall.df[goodSampleID.vect,]
  
  mod.locusVect <- c("KIR3DP1","KIR2DP1","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL2","KIR2DS4",
                     "KIR2DL1","KIR2DS1","KIR2DL5","KIR2DL23","KIR2DS35","KIR3DL1S1")
  iterCall.df <- pingFinalize.combineL23( iterCall.df )
  iterCall.df <- pingFinalize.combineS35( iterCall.df, copyCall.df )
  iterCall.df <- pingFinalize.combineL1S1( iterCall.df, copyCall.df )
  iterCall.df <- pingFinalize.otherLoci( iterCall.df, copyCall.df )
  
  iterCall.df <- iterCall.df[,mod.locusVect]
  
  write.csv( iterCall.df, file.path( resultsDirectory, 'finalAlleleCalls.csv') )
  return( paste0( resultsDirectory, 'finalAlleleCalls.csv') )
}
