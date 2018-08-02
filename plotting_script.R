library(reshape2)

run.generate_copy_number_graphs <- function(countRatioDF, kffDF){
  
  ## Iterate over all KIR loci, create a plot for each one
  for(currentLocus in kirLocusList){
    
    ## Determine the rank of the ratios (x-axis order from lowest to highest)
    countRatioDF$ratioRank <- rank(countRatioDF[,currentLocus], ties.method = 'first')
    
    ## Special dual locus graphing for KIR loci in tight LD
    if(currentLocus == 'KIR2DL1'){
      
      ## Graph 2DL1 and 2DP1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DL1', 'KIR2DP1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
        geom_point() + 
        scale_color_manual(values=c('KIR2DL1' = '#000000', 'KIR2DP1' = '#9b9b9b')) + 
        scale_shape_manual(values=c('KIR2DL1'=20,'KIR2DP1'=20)) +
        scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
        scale_x_continuous(name='Sample rank') +
        ggtitle(currentLocus) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR2DP1'){
      
      ## Graph 2DP1 and 2DL1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DP1', 'KIR2DL1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
              geom_point() + 
              scale_color_manual(values=c('KIR2DP1' = '#000000', 'KIR2DL1' = '#9b9b9b')) + 
              scale_shape_manual(values=c('KIR2DP1'=20,'KIR2DL1'=20)) +
              scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
              scale_x_continuous(name='Sample rank') +
              ggtitle(currentLocus) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR2DL2'){
      
      ## Graph 2DL2 and 2DL3 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DL2', 'KIR2DL3', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
        geom_point() + 
        scale_color_manual(values=c('KIR2DL2' = '#000000', 'KIR2DL3' = '#9b9b9b')) + 
        scale_shape_manual(values=c('KIR2DL2'=20,'KIR2DL3'=20)) +
        scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
        scale_x_continuous(name='Sample rank') +
        ggtitle(currentLocus) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR2DL3'){
      
      ## Graph 2DL3 and 2DL2 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR2DL3', 'KIR2DL2', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
        geom_point() + 
        scale_color_manual(values=c('KIR2DL3' = '#000000', 'KIR2DL2' = '#9b9b9b')) + 
        scale_shape_manual(values=c('KIR2DL3'=20,'KIR2DL2'=20)) +
        scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
        scale_x_continuous(name='Sample rank') +
        ggtitle(currentLocus) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR3DL1'){
      
      ## Graph 3DL1 and 3DS1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR3DL1', 'KIR3DS1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
        geom_point() + 
        scale_color_manual(values=c('KIR3DL1' = '#000000', 'KIR3DS1' = '#9b9b9b')) + 
        scale_shape_manual(values=c('KIR3DL1'=20,'KIR3DS1'=20)) +
        scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
        scale_x_continuous(name='Sample rank') +
        ggtitle(currentLocus) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)))
    }else if(currentLocus == 'KIR3DS1'){
      
      ## Graph 3DS1 and 3DL1 on the same graph
      countRatioDF.melt <- melt(countRatioDF[,c('KIR3DS1', 'KIR3DL1', 'ratioRank')], id.vars='ratioRank', variable.name='locus', value.name='ratio')
      
      print(ggplot(countRatioDF.melt, aes(x=ratioRank, y=ratio, shape=locus, color=locus)) +
        geom_point() + 
        scale_color_manual(values=c('KIR3DS1' = '#000000', 'KIR3DL1' = '#9b9b9b')) + 
        scale_shape_manual(values=c('KIR3DS1'=20,'KIR3DL1'=20)) +
        scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF.melt$ratio))) +
        scale_x_continuous(name='Sample rank') +
        ggtitle(currentLocus) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)))
    }else{
      
      ## Graph all other loci by themselves
      print(ggplot(countRatioDF, aes(x=ratioRank, y=countRatioDF[[currentLocus]])) +
        geom_point() +
        scale_color_manual(values=c(currentLocus = '#000000')) + 
        scale_shape_manual(values=c(currentLocus = 20)) +
        scale_y_continuous(name=paste(currentLocus, ' / KIR3DL3 read ratio'), limits=c(0, max(countRatioDF[[currentLocus]]))) +
        scale_x_continuous(name='Sample rank') +
        ggtitle(currentLocus) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)))
    }
    
    ## Save each plot
    ggsave(file.path(resultsDirectory, paste0(currentLocus, '_copy_number_plot.pdf')), width=10, height=10)
  }
}
