copiedMsfDirectory <- 'Resources/ipdkir_resources/copied_msf/'
referenceResourceDirectory <- 'Resources/ipdkir_resources/reference_resources/'
kirLocusVect <- kir.locus.vect

locusRefList <- general.initialize_locus_ref_object()
locusRefList <- initLocusRef.read_raw_msf(locusRefList, copiedMsfDirectory)
locusRefList <- initLocusRef.create_bed(locusRefList, referenceResourceDirectory, kirLocusFeatureNameList, writeBed=F)

initialize_genotype_workflow <- function(){
  
}

KIR2DL1.filter_contam_reads_from_sam_file <- function(sam_file, sampleName){
  # Custom probes to identify contaminating read pairs responsible for incorrect base at CDS position 13
  probe1 <- "CACCGGCAGCACCATGTCGCTCACGGT"
  probe2 <- "CACAGACAGCACCATGTCGCTCATGGTC"
  probe3 <- "TCCGGCACCACCATGTCGCTCATGGTCA"
  probe4 <- "CACCGGCAGCACCATGTTGCTCATGGTC"
  probe5 <- "CACCGGCAGCACCATGTCGCTCATGGTC"
  probe6 <- "CACCGGCAGCACCATGTCGCTCACGGTC"
  probe7 <- "CACGGTCGTCAGCGTGGC"
  
  # Create data table of SAM file
  #sam.df <- read.table(sam_file, skip = 3, col.names = 1:20,header = F, check.names=F, fill=T,sep = "\t", stringsAsFactors = F)
  sam.df <- read.table(sam_file, skip = 3,header = F, check.names=F, fill=T,sep = "\t", stringsAsFactors = F)
  
  output.samTable <- as.data.table(sam.df)
  colnames(output.samTable)[1]  <- 'read_name'
  colnames(output.samTable)[3]  <- 'reference_name'
  colnames(output.samTable)[4]  <- 'ref_pos'
  colnames(output.samTable)[10] <- 'read_seq'
  
  output.samTable.subset <- output.samTable[grepl(probe1,output.samTable$read_seq)  |
                                              grepl(probe2,output.samTable$read_seq)|
                                              grepl(probe3,output.samTable$read_seq)|
                                              grepl(probe4,output.samTable$read_seq)|
                                              grepl(probe5,output.samTable$read_seq)|
                                              grepl(probe6,output.samTable$read_seq)|
                                              grepl(probe7,output.samTable$read_seq),]
  
  contam_read_ids <- output.samTable.subset$read_name
  
  # Filter contaminating read pairs out of SAM file
  new_sam_file <- paste0(sampleName,"_2DL1.filt.sam")
  sink(new_sam_file)
  con = file(sam_file, "r")
  while (TRUE){
    line = readLines(con, n=1)
    if (length(line) == 0){
      break
    }
    read.id <- unlist(strsplit(line,"\t"))[1]
    if (! read.id %in% contam_read_ids){cat(paste0(line, "\n"))}
  } 
  
  close(con)
  sink()
  
  # Replace SAM file with filtered SAM file
  file.remove(sam_file)
  system2("mv", c(new_sam_file,sam_file))
  return(NULL)
}
