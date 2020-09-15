copiedMsfDirectory <- 'Resources/ipdkir_resources/copied_msf/'
referenceResourceDirectory <- 'Resources/ipdkir_resources/reference_resources/'
kirLocusVect <- kir.locus.vect

# ----- Generate aligned SNP directories -----
# Create directory for storing files relating to allele calling
alleleFileDirectory <- file.path(resultsDirectory,'alleleFiles')
if(!file.exists(alleleFileDirectory)){
  dir.create(alleleFileDirectory)
}

filterAlleleFileDirectory <- file.path(alleleFileDirectory,'filterFiles')
if(!file.exists(filterAlleleFileDirectory)){
  dir.create(filterAlleleFileDirectory)
}

alignmentFileDirectory <- file.path(resultsDirectory,'alignmentFiles')
if(!file.exists(alignmentFileDirectory)){
  dir.create(alignmentFileDirectory)
}

# ----- Generating reference object list for each locus -----
locusRefList <- general.initialize_locus_ref_object()
locusRefList <- initLocusRef.read_raw_msf(locusRefList, copiedMsfDirectory)
locusRefList <- initLocusRef.create_bed(locusRefList, referenceResourceDirectory, kirLocusFeatureNameList, writeBed=F)

# Read in reference allele dataframe
referenceAlleleDF <- read.csv('Resources/genotype_resources/master_haplo_iteration_testing_v10.csv',row.names=1,stringsAsFactors = F)

# Read in reference allele 5'UTR and 3'UTR extensions
UTRextList <- general.read_fasta('Resources/genotype_resources/KIR_UTR_ext.fasta')

# Generate df's to store allele calls for both workflows covering all samples and loci
alleleDFPathList <- allele.setup_results_df( locusRefList, filterLocusConv, resultsDirectory, sampleList )



# ----- Extra support stuff for allele calling -----

# Generate known SNP df's for allele calling
knownSnpDFList <- allele.create_allele_resources(locusRefList, alleleFileDirectory)

# Add KIR2DL23 combined SNP df to the knownSnpDFList [filter specific]
knownSnpDFList <- allele.generate_L23_SNP_df( kirLocusFeatureNameList, knownSnpDFList )

# Add KIR2DS35 combined SNP df to the knownSnpDFList [filter specific]
knownSnpDFList <- allele.generate_S35_SNP_df( kirLocusFeatureNameList, knownSnpDFList )

maskedPosList <- list('KIR2DS4'=paste0('E5_',84:105),'KIR2DL4'='E7_105')

# ----- Support functions for cleaning up allele call results -----
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

## New function for removing masked SNP positions
filter_masked_snps <- function(current.locus, KIR_sample_snps, snpMaskDF){
  ## If the current locus is in the snp mask DF
  if(current.locus %in% row.names(snpMaskDF)){
    ## Pull out the positions that should be masked
    maskPosVect <- tstrsplit(snpMaskDF[current.locus,'hPos.pPos'],'.',fixed=T)[[2]]
    
    ## Remove those positions from allele calling
    KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$position %in% maskPosVect,]
  }
  
  return(KIR_sample_snps)
}

# Remove Problematic SNP positions that are a product of contamination from a different locus
#  - Verified through sanger sequencing experiments
filter_contam_snps <- function(currentSample, current.locus, sampleSnpDF,PctFilter2DL1=0.20, PctFilter2DL5=0.20){
  
  if( current.locus == 'KIR2DL1' ){
    
    cat('\nADJUSTING KIR2DL1 SNPs')
    vcfPath <- currentSample$filterVCFList[[ 'KIR2DL1' ]]
    
    cat('\n\tReading VCF')
    vcfDT <- general.read_VCF(vcfPath)
    
    # 07/16/2020 CONDITION FOR REMOVING pos 14933, 15179 (E8_32, E9_127) 2DL2 contamination
    if( as.numeric(currentSample$geneContent$KIR2DL2) > 0 ){
      
      # If heterozygous, nullify this position
      if( length(unique(sampleSnpDF[,'E9_127'])) > 1 ){
        cat('\n\t\tRemoving E9_127')
        sampleSnpDF[,'E9_127'] <- NA
      }
      
      # if heterozygous, nullify this position
      if( length(unique(sampleSnpDF[,'E8_32'])) > 1 ){
        cat('\n\t\tRemoving E8_32')
        sampleSnpDF[,'E8_32'] <- NA
      }
      
    }
    
    # 14933 & 15179 is 2DL2 contamination that arises because of 2DL2 reads aligning to 2DL1*005 (allele that is not real), Pauls bowtie reference has 005
    #   - These two positions were validated by Danillo after resequencing several individuals
    contam_snp_pos_vect <- c("4788","4822","4950","4991","5009","5011","6614","14933","15179")
    
    contamSnpDT <- vcfDT[ POS %in% contam_snp_pos_vect ]
    
    # Set up data frame to assess the presence of a block of 2DS1 contamination 
    #block_contam_snp_pos <- x[x$V2 %in% c("6614","6733","6759","6760","6786","6818","10072","14418"),]
    block_pos_vect      <- c("6733","6759","6760","6786","6818","10072","14418")
    block_feat_vect <- c("E5_153","E5_179","E5_180","E5_206","E5_238","E6_45","E7_81")
    blockContamSnpDT <- vcfDT[ POS %in% block_pos_vect ]
    
    refR1_vec <- c()
    refR2_vec <- c()
    altR1_vec <- c()
    altR2_vec <- c()
    
    block_dp4 <- tstrsplit(tstrsplit(blockContamSnpDT$INFO,"DP4=")[[2]],";")[[1]]
    for(i in 1:length(block_dp4)){
      dp4_vec <- unlist(strsplit(block_dp4[i],","))
      refR1_vec <- c(refR1_vec, as.numeric(dp4_vec[1]))
      refR2_vec <- c(refR2_vec,as.numeric(dp4_vec[2]))
      altR1_vec <- c(altR1_vec,as.numeric(dp4_vec[3]))
      altR2_vec <- c(altR2_vec,as.numeric(dp4_vec[4]))
    }
    
    blockContamSnpDT$RefR1 <- refR1_vec
    blockContamSnpDT$RefR2 <- refR2_vec
    blockContamSnpDT$AltR1 <- altR1_vec
    blockContamSnpDT$AltR2 <- altR2_vec
    
    ## Additional Conditions
    # 1. Except when: 4991=C or 5011=T, in those cases do not remove position 5009
    # 2. Position 4788: if this position is present with a 'G', 
    #                   then substract the read depth of 5009 from 4788, and remove the 'G' call from pos 4788
    #                   if depth(5009) > depth(4788), then make depth(4788) = 0
    
    alt_call_4950 <- contamSnpDT[POS == "4950"]$ALT
    alt_call_4991 <- contamSnpDT[POS == "4991"]$ALT
    alt_call_5011 <- contamSnpDT[POS == "5011"]$ALT
    alt_call_6614 <- contamSnpDT[POS == "6614"]$ALT
    presence_2DS1 <- FALSE
    
    # Catch for position not in VCF
    if( length(alt_call_4950) != 0 ){
      bool.4950 = T
      if( alt_call_4950 == "A" | alt_call_4950 == "G" ){presence_2DS1 <- TRUE}
    }else{
      bool.4950 = F
    }
    
    if( nrow(contamSnpDT) > 0 && any(is_nuc(sampleSnpDF[,'E4_261'])) ){ # 5009 == E4_261, 4950 == E4_202
      if( presence_2DS1 ){
        if( any(is_nuc(sampleSnpDF[,'E4_202'])) ) { 
          cat('\n\t\tRemoving E4_202')
          sampleSnpDF[,'E4_202'] <- NA 
          }
        
        if( any(is_nuc(sampleSnpDF[,'E4_261'])) ) {
          
          if( length(alt_call_4991) == 0 | length(alt_call_5011) == 0 ){
            bool.4991 = F
            bool.5011 = F
            cat('\n\t\tRemoving E4_261')
            sampleSnpDF[,'E4_261'] <- NA
          }else if( alt_call_4991 == "C" | alt_call_5011 == "T" ){
            bool.4991 = T
            bool.5011 = T
            # If 4991=C (Ref = A) or 5011=T (Ref = G), keep 5009
            sampleSnpDF <- sampleSnpDF
          }else{
            bool.4991 = T
            bool.5011 = T
            cat('\n\t\tRemoving E4_261')
            sampleSnpDF[,'E4_261'] <- NA
          }
        }
      }
    }
    
    # Additional condition
    # - Remove contamination from position 4788 (G) and 4822 (A)
    
    dp4_4950 <- as.character(contamSnpDT[ POS == '4950' ]$INFO)
    if( length(dp4_4950) != 0 ){
      dp4_4950 <- unlist(strsplit(dp4_4950, "DP4="))[2]
      dp4_4950 <- unlist(strsplit(dp4_4950, ";"))[1]
      dp4_4950 <- as.integer(unlist(strsplit(dp4_4950,",")))
      dp4_4950_ref <- dp4_4950[1] + dp4_4950[2]
      dp4_4950_alt <- dp4_4950[3] + dp4_4950[4]
    }
    
    # 4788
    ref_call_4788 <- contamSnpDT[ POS == '4788' ]$REF
    alt_call_4788 <- contamSnpDT[ POS == '4788' ]$ALT
    if( length(ref_call_4788) != 0 & length(alt_call_4788) != 0){
      bool.4788 = T
      if( alt_call_4788 == "." ){ alt_call_4788 <- ref_call_4788 }
      dp4_4788 <- as.character(contamSnpDT[ POS == '4788' ]$INFO)
      dp4_4788 <- unlist(strsplit(dp4_4788, "DP4="))[2]
      dp4_4788 <- unlist(strsplit(dp4_4788, ";"))[1]
      dp4_4788 <- as.integer(unlist(strsplit(dp4_4788,",")))
      dp4_4788_ref <- dp4_4788[1] + dp4_4788[2]
      dp4_4788_alt <- dp4_4788[3] + dp4_4788[4]
    }else{
      bool.4788 = F
    }
    
    # 4822
    ref_call_4822 <- contamSnpDT[ POS == '4822' ]$REF
    alt_call_4822 <- contamSnpDT[ POS == '4822' ]$ALT
    if( length(ref_call_4822) != 0 & length(alt_call_4822) != 0 ){
      bool.4822 = T
      dp4_4822 <- as.character(contamSnpDT[ POS == '4822' ]$INFO)
      dp4_4822 <- unlist(strsplit(dp4_4822, "DP4="))[2]
      dp4_4822 <- unlist(strsplit(dp4_4822, ";"))[1]
      dp4_4822 <- as.integer(unlist(strsplit(dp4_4822,",")))
      dp4_4822_ref <- dp4_4822[1] + dp4_4822[2]
      dp4_4822_alt <- dp4_4822[3] + dp4_4822[4]
    }else{
      bool.4822 = F
    }
    
    dp4_5009 <- as.character(contamSnpDT[ POS == '5009' ]$INFO)
    if( length(dp4_5009) != 0 ){
      bool.5009 = T
      dp4_5009 <- unlist(strsplit(dp4_5009, "DP4="))[2]
      dp4_5009 <- unlist(strsplit(dp4_5009, ";"))[1]
      dp4_5009 <- as.integer(unlist(strsplit(dp4_5009,",")))
      dp4_5009_ref <- dp4_5009[1] + dp4_5009[2]
    }else{
      bool.5009 = F
    }
    
    # depth values needed to resolve position 6614
    ref_call_6614 <- contamSnpDT[ POS == '6614' ]$REF
    alt_call_6614 <- contamSnpDT[ POS == '6614' ]$ALT
    if( length(ref_call_6614) != 0 & length(alt_call_6614) != 0 ){
      bool.6614 = T
      dp4_6614 <- as.character(contamSnpDT[ POS == '6614' ]$INFO)
      dp4_6614 <- unlist(strsplit(dp4_6614, "DP4="))[2]
      dp4_6614 <- unlist(strsplit(dp4_6614, ";"))[1]
      dp4_6614 <- as.integer(unlist(strsplit(dp4_6614,",")))
      dp4_6614_ref <- dp4_6614[1] + dp4_6614[2]
      dp4_6614_alt <- dp4_6614[3] + dp4_6614[4]
    }else{
      bool.6614 = F
    }
    
    # Intronic position 6497 
    ref_call_6497 <- vcfDT[ POS == '6497' ]$REF
    alt_call_6497 <- vcfDT[ POS == '6497' ]$ALT
    if( length(ref_call_6497) != 0 & length(alt_call_6497) != 0 ){
      bool.6497 = T
      dp4_6497 <- as.character(vcfDT[ POS == '6497' ]$INFO)
      dp4_6497 <- unlist(strsplit(dp4_6497, "DP4="))[2]
      dp4_6497 <- unlist(strsplit(dp4_6497, ";"))[1]
      dp4_6497 <- as.integer(unlist(strsplit(dp4_6497,",")))
      dp4_6497_ref <- dp4_6497[1] + dp4_6497[2]
      dp4_6497_alt <- dp4_6497[3] + dp4_6497[4]
    }else{
      bool.6497 = F
    }
    
    # Remove contamination "G" from position 4788,4822
    # subtracted DP4 of 4950 alt call (A or G) from 4788 alt call (C) 
    # If less that 20% of DP4 remains for 4788 alt call after subtraction, force reference call
    ## FIX THIS CONDITION based on the above statements: (Ask Danillo, removing contam G should involve substracting from ref call not alt call?)
    
    if( bool.4788 & bool.4950 ){
      if( presence_2DS1 && ref_call_4788 == "G" ){
        dp4_thresh_4788 <- PctFilter2DL1 * dp4_4788_ref
        diff_4950_vs_4788 <- 0
        if (dp4_4950_alt >= dp4_4788_ref){
          diff_4950_vs_4788 <- 0
        }else {
          diff_4950_vs_4788 <- dp4_4788_ref - dp4_4950_alt
        }
        
        # Force reference call for 4788 - 4950 dp4 is less the the theshold (for old PING the reference is wrong compared to IPD-KIR at this position, use alt call)
        if( diff_4950_vs_4788 < dp4_thresh_4788 ){
          
          if( vcfDT[ POS == '4788' ]$ALT == '.' ){
            repNuc <- NA
          }else{
            repNuc <- vcfDT[ POS == '4788' ]$ALT
          }
          
          cat('\n\t\tForcing ALT E4_40')
          sampleSnpDF[,'E4_40'] <- repNuc
        }
      }
    }
    
    if( bool.4822 & bool.4950 ){
      ## Remove contaminating alternate call from position 4822 if 2DS1 reads are present 
      if( presence_2DS1 && alt_call_4822 == "A" ){
        dp4_thresh_4822 <- PctFilter2DL1 * dp4_4822_alt
        diff_4950_vs_4822 <- 0
        if (dp4_4950_alt >= dp4_4822_alt){
          diff_4950_vs_4822 <- 0
        }else {
          diff_4950_vs_4822 <- dp4_4822_alt - dp4_4950_alt
        }
        
        # Force reference call for 4822 if the difference between the DP4 of 4950 and 4822 small
        if( diff_4950_vs_4822 < dp4_thresh_4822 ){
          cat('\n\t\tForcing REF E4_74')
          sampleSnpDF[,'E4_74'] <- vcfDT[ POS == '4822' ]$REF
        }
        
      }
    }
    
    ## Remove contaminating alternate call for position 6614 if based on subtraction from intronic position 6497 indicative of 2DS1 contamination
    #### New condition: if position 2DS1 is present and 6614 is heterozygous, remove this position from consideration in calling (not reliable)
    if( bool.6614 ){
      if( presence_2DS1 && alt_call_6614 == "C" ){
        cat('\n\t\tRemoving E5_34')
        sampleSnpDF[,'E5_34'] <- NA
      }
    }
    
    ## Handling Block of contaminated SNP positions introduced by 2DS1
    # Block conditions:
    # - if every positon in the block has an alt call (6 total), remove these positions from consideration
    # - if one or more (but not all) of the block has an alt call force reference calls
    if( presence_2DS1 ){
      block_presence <- TRUE
      make_block_ref <- FALSE
      remove_block   <- TRUE
      for( alt_bp in blockContamSnpDT$ALT ){
        if( !is_nodel_nuc(alt_bp) ){
          make_block_ref <- TRUE
          remove_block   <- FALSE
          block_presence <- FALSE
        }
      }
      
      # If at least one SNP is missing in the block, it is 2DS1 contamination, 
      # so the reference calls can be force, otherwise remove all block positions from calling
      if( make_block_ref ){
        cat('\n\t\tForcing REF for snp block')
        # If block is present, force reference calls, since 4950 is A/G, 2DS1 is present
        sampleSnpDF[1,block_feat_vect] <- vcfDT[ POS %in% block_pos_vect ]$REF
        sampleSnpDF[2,block_feat_vect] <- vcfDT[ POS %in% block_pos_vect ]$REF
      }else if( block_presence && remove_block ){
        cat('\n\t\tRemoving snp block')
        sampleSnpDF[1,block_feat_vect] <- NA
        sampleSnpDF[2,block_feat_vect] <- NA
      }
    }
  }
  
  # 2DL5: filter all Exon 1 variable positions. In general 2DL5 suffers from 3DP1 contamination. 
  #       Exon 1 has strong sequence similarity across many of the KIR loci, so PING filters out a lot of reads
  #       in this particular exon region. The low depth makes pos 283 unreliable
  if (current.locus == "KIR2DL5"){
    
    cat('\nADJUSTING KIR2DL5 SNPs')
    
    vcfPath <- currentSample$filterVCFList[[ 'KIR2DL5' ]]
    
    cat('\n\tReading VCF from:',vcfPath)
    vcfDT <- general.read_VCF(vcfPath)
    
    # Contaminated positions (PING Vcf):
    # - 283 -> contamination from most KIR loci, this why it is being omitted from calling
    # - 1081 -> contamination from multiple loci, in paritcular we can tag the 2DS2 and 2DP1 contamination
    contam_snps_to_omit <- c("283")
    cat('\n\t\tRemoving E1_16')
    sampleSnpDF[,'E1_16'] <- NA
    
    # Handle contamination at position 1081
    # - dp4_vec[1] = Ref DP4, dp4_vec[2] = Alt DP4
    
    bool.1081 <- 1081 %in% vcfDT$POS
    bool.1112 <- 1112 %in% vcfDT$POS
    bool.1115 <- 1115 %in% vcfDT$POS
    
    if( bool.1081 ){
      dp4_1081_vec <- get_DP4_select_snp(pos = "1081", vcfDT = vcfDT)
    }
    if( bool.1112 ){
      dp4_1112_vec <- get_DP4_select_snp(pos = "1112", vcfDT = vcfDT)
    }
    if( bool.1115 ){
      dp4_1115_vec <- get_DP4_select_snp(pos = "1115", vcfDT = vcfDT)
    }
    
    if( bool.1081 & bool.1112 & bool.1115 ){
      # Get base calls for positions of interest and their tags
      vcfDT[ POS == '1115' ]$ALT
      alt_1081 <- as.character(vcfDT[ POS == '1081' ]$ALT)
      alt_1112 <- as.character(vcfDT[ POS == '1112' ]$ALT)
      alt_1115 <- as.character(vcfDT[ POS == '1115' ]$ALT)
      
      if (alt_1081 == "C" && dp4_1081_vec[2] > 0){
        total_contam_dp <- 0
        
        alt_dp4_1081 <- as.numeric(dp4_1081_vec[2])
        alt_dp4_1112 <- as.numeric(dp4_1112_vec[2])
        alt_dp4_1115 <- as.numeric(dp4_1115_vec[2])
        
        if (alt_1112 == "A"){total_contam_dp <- total_contam_dp + alt_dp4_1112}
        if (alt_1115 == "T"){total_contam_dp <- total_contam_dp + alt_dp4_1115}
        
        force_ref <- FALSE
        dp4_thresh_1081 <- PctFilter2DL5 * alt_dp4_1081
        if (total_contam_dp >= alt_dp4_1081){
          diff_dp4 <- 0
          force_ref <- TRUE
        } else{
          diff_dp4 <- alt_dp4_1081 - total_contam_dp
          if (diff_dp4 < dp4_thresh_1081){
            force_ref <- TRUE
          }
        }
        
        if (force_ref){
          cat('\n\t\tForcing REF for E2_27')
          sampleSnpDF[,'E2_27'] <- vcfDT[ POS == '1081' ]$REF
        }
      }
    }
  }
  
  return(sampleSnpDF)
}

get_DP4_select_snp <- function(pos, vcfDT){
  
  dp4 <- as.character( vcfDT[ POS == pos ]$INFO )
  dp4 <- unlist(strsplit(dp4, "DP4="))[2]
  dp4 <- unlist(strsplit(dp4, ";"))[1]
  dp4 <- as.integer(unlist(strsplit(dp4,",")))
  dp4_ref <- dp4[1] + dp4[2]
  dp4_alt <- dp4[3] + dp4[4]
  
  dp4_vec <- c(dp4_ref,dp4_alt)
  
  return(dp4_vec)
}

filter_2DL5_allele_str <- function(allele_str, currentSample, DPthresh,PctFilter2DL5=0.20){
  
  vcfPath <- currentSample$filterVCFList[[ 'KIR2DL5' ]]
  
  cat('\n\tReading VCF')
  vcfDT <- general.read_VCF(vcfPath)
  
  presence2DL5A = TRUE
  presence2DL5B = TRUE
  
  ### Extract position 1320
  #    - differentiates 2DL5A from 2DL5B, specifically 2DL5A*00101 vs. 2DL5B*00801
  target_pos = "1320"
  bool.1320 <- 1320 %in% vcfDT$POS
  
  if( bool.1320 ){
    ref_call_1320 <- as.character(vcfDT[ POS == '1320' ]$REF)
    alt_call_1320 <- as.character(vcfDT[ POS == '1320' ]$ALT)
    
    dp4_1320 <- as.character( vcfDT[ POS == '1320' ]$INFO )
    dp4_1320 <- unlist(strsplit(dp4_1320, "DP4="))[2]
    dp4_1320 <- unlist(strsplit(dp4_1320, ";"))[1]
    dp4_1320 <- as.integer(unlist(strsplit(dp4_1320,",")))
    dp4_1320_ref <- dp4_1320[1] + dp4_1320[2]
    dp4_1320_alt <- dp4_1320[3] + dp4_1320[4]
    
    dp4_threshold <- PctFilter2DL5 * dp4_1320_ref
    if (ref_call_1320 == "G" && alt_call_1320 == "."){
      if (dp4_1320_alt < dp4_threshold){
        presence2DL5B <- FALSE
      }
    }
    
    # filter allele string based on 2DL5A/B presence absence
    allele_list <- unlist(strsplit(allele_str,"/"))
    if (presence2DL5B == FALSE && presence2DL5A == TRUE){
      allele_list <- allele_list[grep("KIR2DL5A",allele_list)]
    } else if (presence2DL5B == TRUE && presence2DL5A == FALSE){
      allele_list <- allele_list[grep("KIR2DL5B",allele_list)]
    }
    allele_str <- paste(allele_list,collapse = '/')
  }
  return(allele_str)
}

# Filter 3DL2 genotype calls based on the presence/absence of select alleles
#  - 3DL2*00201 presence: C at poition 193 (1167)
#  - 3DL1*00103 presence: G at pos 285(1259) 
#  - 3DL2*010 presence: G at pos 285(1259) and T at pos 1044(2018)
filter_3DL2_genos <- function(genos.df, sample,vcf.location.gen, DPthresh,PctFilter3DL2=0.20){
  # Format genos data frame into one string
  genos_str_temp <- c()
  for(i in row.names(genos.df)){
    g <- paste(c(genos.df[i,1],"+",genos.df[i,2]),collapse = "")
    genos_str_temp <- c(genos_str_temp,g)
  }
  genos_str_temp <- paste(genos_str_temp,collapse = ":")
  
  # Only use 3DL2 conditions on these specific ambiguities
  target_ambiguities <- c(
    "KIR3DL2_00201+KIR3DL2_00701:KIR3DL2_01001+KIR3DL2_015:KIR3DL2_10701+KIR3DL2_00601",
    "KIR3DL2_00103+KIR3DL2_00701:KIR3DL2_01001+KIR3DL2_00601",
    "KIR3DL2_00101+KIR3DL2_00201:KIR3DL2_00202+KIR3DL2_00103",
    "KIR3DL2_00103+KIR3DL2_10701:KIR3DL2_01001+KIR3DL2_00201")
  
  ### Execute conditions
  genos.df.fin <- genos.df
  if (genos_str_temp %in% target_ambiguities){
    vcf.genomic.df    <- read.table(paste(c(vcf.location.gen,sample,"_3DL2nuc.vcf"),collapse = ""),header = F)
    vcf.df.target.pos <- vcf.genomic.df[vcf.genomic.df$V2 %in% c("1167","1259","2018"),]
    dp_str <- tstrsplit(vcf.df.target.pos$V8,"DP4=")[[2]]
    dp_str <- tstrsplit(dp_str,";")[[1]]
    
    # Extract raw DP 
    raw_dp_str <- tstrsplit(vcf.df.target.pos$V8,"DP=")[[2]]
    raw_dp_str <- tstrsplit(raw_dp_str,";")[[1]]
    vcf.df.target.pos$DP <- raw_dp_str
    
    # Extract total depth(DP4) for both the reference and alternate call at these select positions
    ref_dp4_vec <- c()
    alt_dp4_vec <- c()
    thresh_vec  <- c()
    for (dp4 in dp_str) {
      all_dp4 = tstrsplit(dp4,",")
      ref_dp4 <- as.numeric(all_dp4[1]) + as.numeric(all_dp4[2])
      alt_dp4 <- as.numeric(all_dp4[3]) + as.numeric(all_dp4[4])
      dp4_thresh <- PctFilter3DL2 * ref_dp4
      
      ref_dp4_vec <- c(ref_dp4_vec,ref_dp4)
      alt_dp4_vec <- c(alt_dp4_vec,alt_dp4)
      thresh_vec  <- c(thresh_vec,dp4_thresh)
    }
    
    vcf.df.target.pos$refDP4    <- ref_dp4_vec
    vcf.df.target.pos$altDP4    <- alt_dp4_vec
    vcf.df.target.pos$threshDP4 <- thresh_vec
    
    # Determine the presence/absence of select alleles
    presence010   <- FALSE
    presence00201 <- FALSE
    presence00103 <- FALSE
    
    # 3DL2*010 presence/absence
    entry_285  <- vcf.df.target.pos[vcf.df.target.pos$V2 == "1259",]
    entry_1044 <- vcf.df.target.pos[vcf.df.target.pos$V2 == "2018",]
    if (entry_285$V5 == "G" && as.numeric(entry_285$altDP4) > as.numeric(entry_285$threshDP4)){
      if (as.numeric(entry_285$DP) > DPthresh) {presence00103 <- TRUE}
      if (entry_1044$V5 == "T" && as.numeric(entry_1044$altDP4) > as.numeric(entry_1044$threshDP4)){
        if (as.numeric(entry_285$DP) > DPthresh && as.numeric(entry_1044$DP) > DPthresh){
          presence010 <- TRUE
        }
      }
    }
    
    # 3DL2*00201 presence/absence
    entry_193  <- vcf.df.target.pos[vcf.df.target.pos$V2 == "1167",]
    if (entry_193$V5 == "C" && as.numeric(entry_193$altDP4) > as.numeric(entry_193$threshDP4)){
      if (as.numeric(entry_193$DP) > DPthresh){presence00201 <- TRUE}
    }
    
    # remove genotypes ambiuities that contain alleles determined to be absent by this function
    rows_to_rm <- c()
    for (i in row.names(genos.df)) {
      if(any(grepl("3DL2_010",genos.df[i,])) && presence010 == FALSE){
        rows_to_rm <- c(rows_to_rm,i)
      }
      
      if(any(grepl("3DL2_00201",genos.df[i,])) && presence00201 == FALSE){
        rows_to_rm <- c(rows_to_rm,i)
      }
      
      if(any(grepl("3DL2_00103",genos.df[i,])) && presence00103 == FALSE){
        rows_to_rm <- c(rows_to_rm,i)
      }
    }
    
    genos.df.fin <- genos.df[! row.names(genos.df) %in% rows_to_rm,]
    if (nrow(genos.df.sub) == 0){genos.df.fin <- genos.df}
  }
  
  return(genos.df.fin)
} 

# Implemented KIR2DL1*010 and *01201 filter (check for 2DS1 contaminatin and alters allele calls)
allele.custom_2DL1_allele_filter <- function(currentSample,
                                             formattedAlleleVect){
  
  old.formattedAlleleVect <- formattedAlleleVect
  
  vcfPath <- currentSample$filterVCFList[[ 'KIR2DL1' ]]
  
  cat('\n\tReading VCF')
  vcfDT <- general.read_VCF(vcfPath)
  
  null_allele <- "KIR2DL1*null"
  
  ## 1. get target allele defining positions from Genomic Vcf file
  # Conversion (IPD-KIR to Vcf)
  # 2DL1*010:       2469=3735,2497=3763
  # 2DL1*01201:     1311=2577
  # 2DL1*01201/010: 1321=2587
  unique_pos_2DL1_010   <- c(3735,3763) # SNPs specific to 2DL1*010
  unique_pos_2DL1_01201 <- c(2577)      # SNP specific to 2DL1*01201
  unique_pos_both       <- c(2587)      # SNP specfic to both 2DL1*010 and 2DL1*01201
  target_pos_vec        <- c(unique_pos_2DL1_010,unique_pos_2DL1_01201, unique_pos_both)
  
  #vcf.genomic.df <- read.table(paste(c(vcf.location.gen,sample,"_2DL1nuc.vcf"),collapse = ""),header = F)
  
  #vcf.genomic.df_subset <- vcf.genomic.df[vcf.genomic.df$V2 %in% target_pos_vec,]
  vcfDT.subset <- vcfDT[ POS %in% target_pos_vec ]
  
  # If all target positions are not defined, return original allele calls
  allPos.bool <- all( target_pos_vec %in% vcfDT$POS )
  if( !allPos.bool ){
    return(formattedAlleleVect)
  }
  
  ## 2. determine status of 2DL1*010 and 2DL1*01201
  status_010   <- FALSE
  status_01201 <- FALSE
  # Target SNPs of interest
  allele_filter_pos_frame <- read.table("Resources/genotype_resources/allele_filter_2DL1.genomic.csv",sep = ",",header = T,stringsAsFactors = F,row.names = 1)
  snps_2DL1_00101 <- paste(allele_filter_pos_frame["KIR2DL1*00101",], collapse = "")
  # Reference calls from Genomic Vcf
  subset_ref_calls <- vcfDT.subset$REF
  subset_ref_calls_str <- paste(subset_ref_calls,collapse = "")
  
  
  subset_alt_calls <- c()
  alt_calls_bool <- is_nodel_nuc(vcfDT.subset$ALT)
  if( all(alt_calls_bool) == FALSE && snps_2DL1_00101 == subset_ref_calls_str ){
    subset_alt_calls <- NULL
    status_010   <- FALSE
    status_01201 <- FALSE
  }else{
    subset_alt_calls <- vcfDT.subset$ALT
    names(subset_alt_calls) <- vcfDT.subset$POS
    
    # positive filter
    if( allele_filter_pos_frame["KIR2DL1*010","X3735"] == as.character(subset_alt_calls["3735"]) &&
        allele_filter_pos_frame["KIR2DL1*010","X3763"] == as.character(subset_alt_calls["3763"]) ){
      status_010 <- TRUE
    }
    if( allele_filter_pos_frame["KIR2DL1*01201","X2577"] == as.character(subset_alt_calls["2577"]) ){
      status_01201 <- TRUE
    } 
    # negative filter
    if (as.character(subset_alt_calls["2587"]) != "T"){
      status_010   <- FALSE
      status_01201 <- FALSE
    }
  }
  
  ## 3. Based on status keep or remove 2DL1 alleles from genotype dataframe
  
  if(!status_010){ 
    formattedAlleleVect <- grep('KIR2DL1*010', formattedAlleleVect, fixed=T, value=T, invert=T)
  }
  
  if(!status_01201){
    formattedAlleleVect <- grep('KIR2DL1*01201', formattedAlleleVect, fixed=T, value=T, invert=T)
  }
  
  if( length(formattedAlleleVect) == 0 ){
    formattedAlleleVect <- old.formattedAlleleVect
  }
  
  return(formattedAlleleVect)
}


# ----- Container functions to coordinate FILTER and ITER alignments

# function to run the ITER alignments
ping_iter.run_alignments <- function( currentSample, threads ){
  
  cat('\nLoading ref DF')
  currentSample <- sampleObj.loadRefDF(currentSample, referenceAlleleDF) # Subset reference allele dataframe by present loci, save to sample object
  
  cat('\nWriting reference files')
  currentSample <- sampleObj.writeRefFastaBed(currentSample, locusRefList, alignmentFileDirectory) # Write fasta reference file for sample object based on refDF
  
  currentSample <- sampleObj.iterBowtie2Index(currentSample, bowtie2Build, threads) # Converts fasta file from previous line into a bowtie2 index
  
  currentSample <- sampleObj.iterBowtie2Align(currentSample, bowtie2, threads, deleteSam=T) # Align sample to bowtie2 index
  
  currentSample <- sampleObj.iterVCFGen(currentSample, samtools, bcftools, threads) # Convert SAM file into VCF
  
  return( currentSample )
  
}

# Iter SNP consolidation workflow + allele calling
ping_iter.allele <- function( currentSample ){
  
  # Initialze data frames for consolidating SNPs across all samples for each gene
  snpDFList <- allele.initialize_SNP_tables(alleleFileDirectory, locusRefList, referenceAlleleDF, workflow='iter')
  
  currentSample <- allele.iter_alignments_to_snp_dfs(currentSample, locusRefList, referenceAlleleDF, minDP, kirLocusFeatureNameList)
  
  currentSample <- allele.combine_iter_snps(currentSample, snpDFList)
  
  if( 'KIR2DL2' %in% rownames(currentSample$refAlleleDF) & 'KIR2DL3' %in% rownames(currentSample$refAlleleDF) ){
    currentSample <- allele.iter_combine_KIR2DL23( currentSample, knownSnpDFList, alleleFileDirectory, snpDFList )
  }
  
  currentSample <- allele.setup_iter_allele_call_df(currentSample)
  
  cat('\n\nFinding allele matches for',currentSample$name)
  for( currentLocus in colnames(currentSample[['iterAlleleCallDF']]) ){
    currentSample <- allele.call_allele(currentSample, currentLocus, alleleFileDirectory, knownSnpDFList, alleleDFPathList$iter$newAllelePath, filterLocusConv, workflow = 'iter')
  }
  
  cat('\nWriting allele matches to', alleleDFPathList$iter$alleleCallPath )
  currentSample <- allele.save_call( currentSample, alleleDFPathList$iter$alleleCallPath, 'iter' )
  
  return( currentSample )
}

# Filter align workflow
ping_filter.run_alignments <- function( currentSample, threads ){
  
  currentSample <- sampleObj.filterAlign.setup(currentSample, alignmentFileDirectory)
  
  if(currentSample[['filterRefDirectory']] == 'failed'){
    return(currentSample)
  }
  
  #locusPresenceList <- sapply(intersect(names(currentSample$geneContent), names(currentSample$copyNumber)), function(locusName){
  #  as.numeric(currentSample$copyNumber[[locusName]]) > 0 | as.numeric(currentSample$geneContent[[locusName]]) > 0
  #})
  
  locusPresenceList <- sapply(kirLocusVect, function(locusName){ # GC / CN flexible solution
    
    cnBool <- F
    if( length(currentSample$copyNumber) > 0 ){
      cnBool <- as.numeric(currentSample$copyNumber[[locusName]]) > 0
    }
    
    gcBool <- F
    if( length(currentSample$geneContent) > 0 ){
      gcBool <- as.numeric(currentSample$geneContent[[locusName]]) > 0
    }
    
    return( cnBool | gcBool )
  })
  
  currentSample[[ 'samplePresentLocusVect' ]] <- c()
  
  currentSample <- sampleObj.filterAlign.KIR3DL3(currentSample, bowtie2, samtools, bcftools, threads)
  currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR3DL3' )
  
  currentSample <- sampleObj.filterAlign.KIR3DL2(currentSample, bowtie2, samtools, bcftools, threads)
  currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR3DL2' )
  
  if( locusPresenceList[['KIR3DL1']] & locusPresenceList[['KIR3DS1']] ){
    currentSample <- sampleObj.filterAlign.KIR3DL1S1(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR3DL1het', 'KIR3DS1het' )
  }
  
  if( locusPresenceList[['KIR3DL1']] & !locusPresenceList[['KIR3DS1']] ){
    currentSample <- sampleObj.filterAlign.KIR3DL1(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR3DL1' )
  }
  
  if( locusPresenceList[['KIR3DS1']] & !locusPresenceList[['KIR3DL1']] ){
    currentSample <- sampleObj.filterAlign.KIR3DS1(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR3DS1' )
  }
  
  # If 2DS5 present regardless of 2DS3
  if( locusPresenceList[['KIR2DS5']] ){
    currentSample <- sampleObj.filterAlign.KIR2DS35(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DS35' )
  }
  
  # If 2DS3 present and not 2DS5
  if( locusPresenceList[['KIR2DS3']] & !locusPresenceList[['KIR2DS5']] ){
    currentSample <- sampleObj.filterAlign.KIR2DS3(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DS3' )
  }
  
  if( locusPresenceList[['KIR2DS4']] ){
    currentSample <- sampleObj.filterAlign.KIR2DS4(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DS4' )
  }
  
  if( locusPresenceList[['KIR2DP1']] ){
    currentSample <- sampleObj.filterAlign.KIR2DP1(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DP1' )
  }
  
  currentSample <- sampleObj.filterAlign.KIR2DL23(currentSample, bowtie2, samtools, bcftools, threads)
  currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DL23', 'KIR2DL2', 'KIR2DL3' )
  
  if( locusPresenceList[['KIR2DL1']] ){
    currentSample <- sampleObj.filterAlign.KIR2DL1(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DL1' )
  }
  
  if( locusPresenceList[['KIR2DL5']] ){
    currentSample <- sampleObj.filterAlign.KIR2DL5(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DL5' )
  }
  
  if( locusPresenceList[['KIR2DL4']] ){
    currentSample <- sampleObj.filterAlign.KIR2DL4(currentSample, bowtie2, samtools, bcftools, threads)
    currentSample[[ 'samplePresentLocusVect' ]] <- c( currentSample[[ 'samplePresentLocusVect' ]], 'KIR2DL4' )
  }
  
  return(currentSample)
  
}

# Filter allele calling
ping_filter.allele <- function( currentSample ){
  filterSnpDFList <- allele.initialize_filter_SNP_tables(filterAlleleFileDirectory, locusRefList, names(filterLocusConv), filterLocusConv)
  
  # Filter SNP processing workflow + allele calling
  cat('\n\nProcessing filtration alignments for',currentSample$name)
  currentSample <- allele.filter_alignments_to_snp_dfs( currentSample, locusRefList, minDP, kirLocusFeatureNameList, filterRefFastaList, filterLocusConv, filterSnpDFList, knownSnpDFList )
  
  currentSample <- allele.setup_filter_allele_call_df( currentSample )
  
  if( currentSample[['filterRefDirectory']] == 'failed' ){
    return( currentSample )
  }
  
  for( currentLocus in currentSample$samplePresentLocusVect ){
    currentSample <- allele.call_allele( currentSample, currentLocus, filterAlleleFileDirectory, knownSnpDFList, alleleDFPathList$filter$newAllelePath, filterLocusConv, 'filter' )
  }
  
  currentSample <- allele.save_call( currentSample, alleleDFPathList$filter$alleleCallPath, 'filter' )
  
  return(currentSample)
}

