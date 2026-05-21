# This script
#   1) Solve ambiguous genotype using allele frequency
#   2) Perform special processing for KIR3DL2 and KIR3DL3 genotyping using PHASE

library(dplyr)
library(stringr)
library(purrr)
library(tidyverse)
library(tidyr)

### ALLELE FREQUENCY TO SOLVE AMBIGUOUS GENOTYPE
args <- commandArgs(trailingOnly = TRUE)
ping_dir <- getwd()
setwd(args) #Set this to your PING results directory
# setwd('/home/rsuseno/PING_input_output/subset_LRC_output_cli')
# setwd('/home/rsuseno/PING_OUTPUTS/IND130_132-133_results')
# setwd('/home/rsuseno/PING_input_output/KIR3DL2_rerun_Feb2025_results')
# setwd('/home/rsuseno/PING_OUTPUTS/OneLambdaNewProbe')
dir.create("phase_output", showWarnings = FALSE)
final_allele <- read.csv('finalAlleleCalls.csv')

# Cut all alleles call to three digits
process_genotypes <- function(genotypes) {
  revisedGenotypes <- c()  # Initialize an empty vector
  
  # Create a regex pattern for all known KIR genes
  all_kir_genes <- c("KIR3DP1", "KIR2DP1", "KIR2DS2", "KIR2DL4", "KIR3DL3",
                     "KIR3DL2", "KIR2DS4", "KIR2DL1", "KIR2DS1", 
                     "KIR2DL5", "KIR2DL5A", "KIR2DL5B",
                     "KIR2DL23", "KIR2DL2", "KIR2DL3",
                     "KIR2DS35", "KIR2DS3", "KIR2DS5", 
                     "KIR3DL1S1","KIR3DL1", "KIR3DS1")
  gene_regex <- paste0("(", paste(all_kir_genes, collapse = "|"), ")")
  
  for (g in genotypes) {
    # Split genotypes on spaces (since multiple genotype pairs are space-separated)
    genotype_list <- unlist(strsplit(g, " "))
    
    # Process each genotype in the list
    processed_list <- sapply(genotype_list, function(pair) {
      # Split individual alleles within a pair (assuming alleles are separated by "+")
      alleles <- unlist(strsplit(pair, "\\+"))
      
      # Process each allele
      processed_alleles <- sapply(alleles, function(a) {
        if (grepl("null|unresolved", a)) {
          return(a)  # Keep "null" or "unresolved" unchanged
        } else {
          # Generalized regex to match any known KIR gene and trim to 3 digits
          return(gsub(paste0(gene_regex, "\\*(\\d{3})\\d*"), "\\1*\\2", a, perl = TRUE))
        }
      })
      
      return(paste(processed_alleles, collapse = "+"))  # Recombine into a genotype pair
    })
    
    # Recombine the processed genotypes into a single string
    revisedGenotypes <- c(revisedGenotypes, paste(processed_list, collapse = " "))
  }
  
  return(revisedGenotypes)
}
cut_allele <- final_allele
for (i in 2:length(colnames(final_allele))){
  gene = colnames(final_allele[i])
  result <- process_genotypes(final_allele[[gene]])
  cut_allele[[gene]] <- result
}

# Get allele frequency and fix the calls accordingly
get_frequency <- function(gene, alleleFile) {
  freq <- list()
  count_problematic <- 0
  for (i in alleleFile[[gene]]) {
    if (grepl(" |unresolved|failed", i)) {
      count_problematic <- count_problematic + 1
      next
    }
    
    genotypes <- unlist(strsplit(i, "\\+"))
    
    for (allele in genotypes) {
      if (allele %in% names(freq)) {
        freq[[allele]] <- freq[[allele]] + 1
      } else {
        freq[[allele]] <- 1
      }
    }
  }
  count_ok <- nrow(final_allele) - count_problematic
  allele_names <- (names(freq))
  freq <- as.numeric(freq)
  allele_frequency <- freq / (2 * count_ok)
  freq_list <- setNames(as.list(allele_frequency), allele_names)
  return(freq_list)
}
get_revised_calls <- function(gene, final_allele, freq){
  revisedCalls <- list()  # Initialize as an empty list
  for (idx in seq_along(final_allele[[gene]])) {
    i <- final_allele[[gene]][idx]  # Access element by index
    genotype <- unique(unlist(strsplit(i, " ")))
    
    # Case 1: Unresolved or failed samples, leave as is
    if (grepl("unresolved|failed", i)) {
      revisedCalls[[idx]] <- genotype
    }
    # Case 2: Normal genotype (no space), leave as is
    # else if (!grepl(" ", i)) {
    else if (length(genotype) == 1){
      revisedCalls[[idx]] <- genotype
    }
    
    # Case 3: Ambiguous genotype (contains space)
    # else if (grepl(" ", i)) {
    else if (length(genotype > 1)) {
      candidates <- genotype
      valid_candidates <- c()  # Temporary vector to store valid candidates
      
      for (c in candidates) {
        alleles <- unlist(strsplit(c, "\\+"))  # Split candidate into alleles
        
        # Ensure all alleles exist in the frequency table
        if (all(alleles %in% names(freq))) {
          valid_candidates <- c(valid_candidates, c)  # Store valid genotype
        }
      }
      
      # Store as a list of possible calls (or NA if no valid candidates were found)
      if (length(valid_candidates) > 0) {
        revisedCalls[[idx]] <- valid_candidates
      } else {
        revisedCalls[[idx]] <- genotype  # in the case that none of the candidates have all their alleles in the freq table, leave it as ambiguous
      }
    }
  }
  
  return(revisedCalls)
}

freq_allele <- cut_allele
for (i in 2:length(colnames(final_allele))){
  gene = colnames(final_allele[i])
  freq <- get_frequency(gene,cut_allele)
  revised <- get_revised_calls(gene, cut_allele, freq)
  freq2 <- Filter(function(x) x >= 0.01, freq)
  revised2 <- list()
  for (idx in seq_along(revised)){
    genotype = revised[idx][[1]]
    # print(genotype)
    if (length(genotype) == 1){
      revised2[[idx]] <- genotype
    }
    else{
      candidates <- genotype
      valid_candidates <- c()  # Temporary vector to store valid candidates
      
      for (c in candidates) {
        alleles <- unlist(strsplit(c, "\\+"))  # Split candidate into alleles
        
        # Ensure all alleles exist in the frequency table
        if (all(alleles %in% names(freq2))) {
          valid_candidates <- c(valid_candidates, c)  # Store valid genotype
        }
      }
      
      # Store as a list of possible calls (or NA if no valid candidates were found)
      if (length(valid_candidates) > 0) {
        revised2[[idx]] <- valid_candidates
      } else {
        revised2[[idx]] <- genotype  # in the case that none of the candidates have all their alleles in the freq table, leave it as ambiguous
      }
    }
  }
  freq_allele[[gene]] <- revised2
  freq_allele[[gene]] <- sapply(freq_allele[[gene]], paste, collapse = " ")
}

### KIR3DL2/3 GENOTYPING SCRIPT 1 ###
prep_PHASE_input <- function(gene){
  allele_snps <- NULL
  if (gene == 'KIR3DL2'){
    # Read excel file (3DL2 specific)
    # KIR3DL2_alleleSNPs3.csv <- readxl::read_excel("/home/rsuseno/PING/Resources/PHASE_resources/KIR3DL2_alleleSNPs_copy.xlsx")[-1, ]
    KIR3DL2_alleleSNPs3.csv <- readxl::read_excel(file.path(ping_dir, "Resources", "PHASE_resources", "KIR3DL2_alleleSNPs_copy.xlsx"))[-1, ]
    mapping <- setNames(c("1", "3", "2", "2"), c("A", "T", "C", "G"))
    KIR3DL2_filtered <- KIR3DL2_alleleSNPs3.csv %>% 
      select(where(~ length(unique(.x[! .x %in% c("*", ".")])) > 1)) %>% 
      mutate(across(everything(), ~ {
        uniq <- unique(.x[! .x %in% c("*", ".")])
        if(length(uniq) > 2) str_replace_all(.x, mapping) else .x
      }))
    colnames(KIR3DL2_filtered)[1] <- "Allele"
    allele_snps <- KIR3DL2_filtered
  }
  else if (gene == 'KIR3DL3'){
    # Read excel file (3DL3 specific)
    # KIR3DL3_alleles_snps <- read_table("/home/rsuseno/PING/Resources/PHASE_resources/3DL3_44SNP_POS_ping2.txt", col_names = FALSE)
    KIR3DL3_alleles_snps <- read_table(file.path(ping_dir, "Resources", "PHASE_resources", "3DL3_44SNP_POS_ping2.txt"), col_names = FALSE)
    names(KIR3DL3_alleles_snps) <- as.character(unlist(KIR3DL3_alleles_snps[1, ])); KIR3DL3_alleles_snps <- KIR3DL3_alleles_snps[-1, ]
    mapping <- setNames(c("1", "3", "2", "2"), c("A", "T", "C", "G"))
    KIR3DL3_filtered <- KIR3DL3_alleles_snps %>% 
      select(where(~ length(unique(.x[! .x %in% c("*", ".")])) > 1)) %>% 
      mutate(across(everything(), ~ {
        uniq <- unique(.x[! .x %in% c("*", ".")])
        if(length(uniq) > 2) str_replace_all(.x, mapping) else .x
      }))
    colnames(KIR3DL3_filtered)[1] <- "Allele"
    allele_snps <- KIR3DL3_filtered
  }
  
  # Task 2: Convert the PING Output into a PHASE input # run this on console for more accurate result
  # Filter for KIR3DL2 and DP in the file names of the target samples 
  # "snp_output" is name of folder with samples in directory
  specific_folder <- "snp_output"
  
  # List all CSV in target folder
  specific_files <- list.files(path = specific_folder, pattern = "\\.csv$", full.names = TRUE)
  
  # start list to store data frames
  target_files <- list()
  
  # Repeat for each CSV file
  for (csv_file in specific_files) {
    # Extract the file name
    sample_name <- tools::file_path_sans_ext(basename(csv_file))
    
    # Check if the file name contains both "KIR3DL2" or "KIR3DL3 and "DP"
    if (grepl(gene, sample_name) && grepl("DP", sample_name)) {
      # Load the CSV file as a data frame
      csv_data <- read.csv(csv_file)
      
      # Store the data frame in the list
      target_files[[sample_name]] <- csv_data
    }
  }
  
  #Task: retain only positions with variable positions  
  # isolate names in filtered KIR3DL2
  isolated_columns <- colnames(allele_snps)
  
  # Use map to iterate over target_3DL2_files and select only the columns present in KIR3DL2_filtered
  variable_list <- map(target_files, ~ .x %>% select(all_of(intersect(isolated_columns, colnames(.x)))))
  
  # Applying S and M PHASE Criteria to genotypes
  # Letters to numbers: A -> "1", T -> "3", C and G -> "2"
  mapping <- c("A" = "1", "T" = "3", "C" = "2", "G" = "2")
  letter_to_number<- function(x) str_replace_all(x, mapping)
  
  # Identify M (multi-allelic) marked columns in KIR3DL3_filtered (skip Column 1 because its Allele)
  M_columns <- names(allele_snps)[-1][
    sapply(allele_snps[-1], function(col) {
      suppressWarnings(all(!is.na(as.numeric(as.character(col)))))
    })
  ]
  
  # Apply depth coverage of > 20 and a depth ratio of 25%, Change M marked positions to Multi-allelic condition (numeric),  set the S marked positions as letters 
  threshold_samples <- imap_dfr(variable_list, function(dataset, samples_name) {
    dataset %>% 
      summarise(across(where(is.numeric), ~{
        sum_of_rows <- sum(.x[1:6], na.rm = TRUE)
        target_rows <- .x[1:4] > 0.25 * sum_of_rows
        set_labels <- c("A", "T", "C", "G")[target_rows]
        col_type <- if (cur_column() %in% M_columns) "M" else "S"
        if (sum_of_rows < 20) return(if (col_type == "M") "-1-1" else "??")
        if (sum(target_rows) > 2) return(NA_character_)
        res <- if (length(set_labels) == 1) paste0(set_labels, set_labels) else paste0(set_labels, collapse = "")
        if (col_type == "M") letter_to_number(res) else res
      }, .names = "{.col}_result")) %>% 
      pivot_longer(everything(), names_to = "Position", values_to = "Nucleotides") %>% 
      mutate(Sample_ID = samples_name,
             Position = str_remove(Position, "_result$"))
  }) %>% 
    filter(!is.na(Nucleotides))
  
  threshold_samples <- threshold_samples %>%
    filter(!grepl("^I", Position) & Position != "X5UTR_24")
  
  # --- Additional Code: Remove samples with >50% unknown results ---
  unknown_samples <- threshold_samples %>%
    group_by(Sample_ID) %>%
    summarise(
      total = n(),
      unknown = sum(Nucleotides %in% c("??", "-1-1"))
    ) %>%
    filter(unknown / total > 0.5) %>%
    pull(Sample_ID)
  
  if (length(unknown_samples) > 0) {
    message("Samples removed due to >50% unknown results: ", paste(unknown_samples, collapse = ", "))
  }
  
  threshold_samples <- threshold_samples %>%
    filter(!(Sample_ID %in% unknown_samples))
  # --- End Additional Code ---
  
  # Phase Prep to Get Locus Types 
  # Changing the threshold_samples frame into PHASE conditions 
  # PHASE_prep <- threshold_samples %>%
  #   distinct(Position, .keep_all = TRUE) %>%  # Keep only unique samples
  #   mutate(Nucleotides = case_when(
  #     grepl("^[A|G|T|C]{2}$", Nucleotides) ~ "S",  # if it is  2 letters -> A, G, T, C
  #     grepl("^[1-4]{2}$", Nucleotides) ~ "M",       # if it is  2 numbers -> 1 to 4
  #     Nucleotides == "??" ~ "?",
  #     Nucleotides == "-1-1" ~ "-1",
  #     TRUE ~ Nucleotides                           
  #   ))
  PHASE_prep <- tibble(Position = names(allele_snps)[-1]) %>%
    mutate(Locus_type = if_else(Position %in% M_columns, "M", "S"))
  
  ## Creates dataframe that takes that summarizes the nucleotides in each position for each of the 
  samples_summarized <- threshold_samples %>%
    group_by(Sample_ID) %>%
    summarise(
      ID_Sample = first(Sample_ID),  
      Haplotype = paste0(Nucleotides, collapse = "")
    ) %>%
    ungroup() %>%
    select(ID_Sample, Haplotype)
  
  # For Haplotypes to have spaces between each value
  #Check if the string length in samples_summarized is x2 of the Locus type length in PHASE_prep (to ensure that the genotype length matches the locus type length)
  samples_summarized <- samples_summarized %>%
    filter(nchar(Haplotype) == 2 * nrow(PHASE_prep))
  
  samples_summarized <- samples_summarized %>%
    mutate(Haplotype = sapply(Haplotype, function(x) {
      spaced <- paste(unlist(str_split(x, "")), collapse = " ")
      str_replace_all(spaced, "-\\s+1", "-1")
    }))
  
  
  ##Final txt file to be inputted into PHASE v 2.1 Software 
  #paste the sample ID and their haplotype 
  sample_info <- paste(paste(samples_summarized$ID_Sample, samples_summarized$Haplotype, sep = "\n"), collapse = "\n")
  
  #paste the locus type in PHASE format
  Locus_type_string <- paste(PHASE_prep$Locus_type, collapse = "")
  
  #List the number of samples and locus types in pipeline
  number_of_samples <- length(samples_summarized$ID_Sample)
  number_of_locus_types <- length(PHASE_prep$Locus_type)
  
  final_output <- paste(number_of_samples, number_of_locus_types,  Locus_type_string, sample_info, sep = "\n")
  
  if (gene == 'KIR3DL2'){
    writeLines(final_output, "phase_output/project_PHASE_3DL2.txt")
  }
  else if (gene == 'KIR3DL3'){
    writeLines(final_output, "phase_output/project_PHASE_3DL3.txt")
  }
  
  return(final_output)
}
kir3dl2 <- prep_PHASE_input('KIR3DL2')
kir3dl3 <- prep_PHASE_input('KIR3DL3')

### CALL PHASE ###
# INPUT: project_PHASE_3DL*.txt
# OUTPUT: PHASE_3DL*.out, PHASE_3DL*.out_pairs
system2('PHASE',
        args = c('-f1', 
                 'phase_output/project_PHASE_3DL2.txt', 
                 'phase_output/PHASE_3DL2.out'))
system2('PHASE',
        args = c('-f1', 
                 'phase_output/project_PHASE_3DL3.txt', 
                 'phase_output/PHASE_3DL3.out'))

### KIR3DL2/3 GENOTYPING SCRIPT 2 ###
process_PHASE_output <- function(gene){
  allele_haplotypes <- NULL
  test_out_pairs <- NULL
  test.out <- NULL
  if (gene == 'KIR3DL2'){
    # allele_haplotypes <- read.table("/home/rsuseno/PING/Resources/PHASE_resources/KIR3dl2_dictionary_output.txt", header = FALSE, sep = "", strip.white = TRUE, fill = TRUE, col.names = c("Allele", "Haplotypes"))
    allele_haplotypes <- read.table(file.path(ping_dir, "Resources", "PHASE_resources", "KIR3dl2_dictionary_output.txt"), header = FALSE, sep = "", strip.white = TRUE, fill = TRUE, col.names = c("Allele", "Haplotypes"))
    test_out_pairs <- readLines("phase_output/PHASE_3DL2.out_pairs")
    test.out <- read.table("phase_output/PHASE_3DL2.out",
                           header = FALSE, 
                           sep = "\t", 
                           strip.white = FALSE,
                           stringsAsFactors = FALSE,
                           colClasses = "character")
  }
  else if (gene == 'KIR3DL3'){
    # allele_haplotypes <- read.table("/home/rsuseno/PING/Resources/PHASE_resources/3DL3_alleles_haplotypes.txt", header = FALSE, sep = "", strip.white = TRUE, fill = TRUE, col.names = c("Allele", "Haplotypes"))
    allele_haplotypes <- read.table(file.path(ping_dir, "Resources", "PHASE_resources", "3DL3_alleles_haplotypes.txt"), header = FALSE, sep = "", strip.white = TRUE, fill = TRUE, col.names = c("Allele", "Haplotypes"))
    test_out_pairs <- readLines("phase_output/PHASE_3DL3.out_pairs")
    test.out <- read.table("phase_output/PHASE_3DL3.out",
                           header = FALSE, 
                           sep = "\t", 
                           strip.white = FALSE,
                           stringsAsFactors = FALSE,
                           colClasses = "character")
  }
  
  # Locate the markers BEGINLIST_SUMMARY and ENDLIST_SUMMARY in test.out and extract the haplotype lines 
  begin_location <- which(trimws(test.out$V1) == "BEGIN LIST_SUMMARY")
  end_location   <- which(trimws(test.out$V1) == "END LIST_SUMMARY")
  
  
  if(length(begin_location) == 0 || length(end_location) == 0) {
    stop("Cannot find begin and end markers in test.out")
  }
  
  # Extract only the haplotype string by removing the start and end numbers in each line [assuming there are spaces between the leading number, haplotype string, and ending number]
  extracted_lines <- test.out$V1[(begin_location + 1):(end_location - 1)]
  
  extracted <- sapply(extracted_lines, function(line) {
    # remove spaces to match dictionary  
    line <- trimws(line)
    
    target <- sub("^[0-9]+\\s+(.*?)\\s+[0-9]+(\\.[0-9]+)?$", "\\1", line)
    
    gsub("\\s+", "", target)
    
  })
  
  # Update the 3DL3 dictionary with any new haplotypes (matched as new_allele_)
  new <- extracted[!extracted %in% allele_haplotypes$Haplotypes]
  if (length(new)) {
    new_alleles <- paste0("new_allele_", seq_along(new))
    allele_haplotypes <- rbind(allele_haplotypes, data.frame(Allele = new_alleles,
                                                             Haplotypes = new,
                                                             stringsAsFactors = FALSE))
    
    # maintain row names as numbers 
    rownames(allele_haplotypes) <- NULL

    # Write new alleles to a text file
    output <- data.frame(Allele = new_alleles,
                       Haplotypes = new,
                       stringsAsFactors = FALSE)

    write.table(output,
              file <- paste0("newAlleles_", gene, ".txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE)  # <-- important if running multiple times
  } 
  
  #Task 2 
  ## Convert test.out_pairs txt to df with IND, Haplotype1, Haplotype2, and Probability columns
  # Remove spaces and locate lines with "IND:"
  lines <- trimws(test_out_pairs)
  
  group <- cumsum(grepl("^IND:", lines))
  
  # Create test_out_pairs dataframe focusing on each IND group (with their associated Haplotype 1, Haplotype 2, and Probability)
  test_out_pairs <- do.call(rbind, lapply(unique(group), function(g) {
    subvec <- lines[group == g]
    
    id <- sub("^IND:\\s*", "", subvec[1])
    
    # Lines with Haplotype 1, Haplotype 2, and Probability as variation entries
    variations <- subvec[-1]
    if (length(variations) == 0) return(NULL)
    
    variation_rows <- do.call(rbind, lapply(variations, function(x) {
      parts <- gsub("\\s+", "", trimws(unlist(strsplit(x, ","))))
      # If there are fewer than 3 parts, fill in with NA.
      if (length(parts) < 3) parts <- c(parts, rep(NA, 3 - length(parts)))
      parts[1:3]
    }))
    # Create the dataframe 
    data.frame(IND = id,
               Haplotype_1       = variation_rows[, 1],
               Haplotype_2       = variation_rows[, 2],
               Probability  = variation_rows[, 3],
               stringsAsFactors = FALSE)
  }))
  
  # For each sample, match the haplotypes with the alleles in the 3DL3 dictinary. 
  test_out_pairs$Allele_1 <- allele_haplotypes$Allele[match(test_out_pairs$Haplotype_1, allele_haplotypes$Haplotypes)]
  test_out_pairs$Allele_2 <- allele_haplotypes$Allele[match(test_out_pairs$Haplotype_2, allele_haplotypes$Haplotypes)]
  
  #Considering test.out only has haplotypes with higher probabilities, some haplotypes may not be found in either test.out or the 3DL3 dictionary. These haplotypes often have lower probabilities. 
  # Option: If a sample has more than two pairs of haplotypes, only retain pairs with higher probabilities. 
  
  test_out_pairs <- test_out_pairs %>%
    group_by(IND) %>%
    slice_max(order_by = Probability, n = 2, with_ties = FALSE) %>%  
    ungroup()
  
  
  # Option: further filtering of dataframe to focus on samples with just two pairs of haplotypes. 
  # If these haplotypes are not in test.out, identified as NA under the Allele columns, then only retain haplotype pairs with the higher probability (only one of the pairs of haplotypes remains) 
  test_out_pairs <- test_out_pairs %>%
    group_by(IND) %>%
    group_modify(~ if(nrow(.x) == 2 && any(is.na(.x$Allele_1) | is.na(.x$Allele_2))) {
      slice_max(.x, Probability, n = 1, with_ties = FALSE)
    } else {
      .x
    }) %>%
    ungroup()
  
  # Cut to three digits and merge alleles into one column
  phase_df <- NULL
  if (gene == 'KIR3DL2'){
    phase_df <- test_out_pairs %>%
      mutate(
        Allele_1 = ifelse(str_starts(Allele_1, "KIR"), sub("(KIR[^*]*\\*\\d{3})\\d*", "\\1", Allele_1), Allele_1),
        Allele_2 = ifelse(str_starts(Allele_2, "KIR"), sub("(KIR[^*]*\\*\\d{3})\\d*", "\\1", Allele_2), Allele_2),
        Alleles = paste(Allele_1, Allele_2, sep = "+")
      )
  }
  if (gene == 'KIR3DL3'){
    phase_df <- test_out_pairs %>%
      mutate(
        Allele_1 = ifelse(str_starts(Allele_1, "KIR"), sub("(KIR[^_]*_\\d{3})\\d*", "\\1", Allele_1), Allele_1),
        Allele_2 = ifelse(str_starts(Allele_2, "KIR"), sub("(KIR[^_]*_\\d{3})\\d*", "\\1", Allele_2), Allele_2),
        Alleles = paste(Allele_1, Allele_2, sep = "+")
      )
  }
  
  # If there's multiple entries, pick the one with prob > 0.6
  phase_df <- phase_df %>%
    group_by(IND) %>%
    filter(n() == 1 | Probability >= 0.6) %>%
    ungroup() %>%
    select(IND, Alleles)
  
  return(phase_df)
}
phase_kir3dl2 <- process_PHASE_output('KIR3DL2')
phase_kir3dl3 <- process_PHASE_output('KIR3DL3')
rm(list = setdiff(ls(), c("freq_allele", "phase_kir3dl2", "phase_kir3dl3")))

### MERGE ALLELE FREQUENCY AND PHASE ###
# Extract ID from X (first section before '_')
freq_allele <- freq_allele %>%
  mutate(IND_ID = str_extract(X, "^[^_]+"))  # Extracts the first section before '_'

# Extract ID from IND (third section after two '_')
phase_kir3dl2 <- phase_kir3dl2 %>%
  mutate(IND_ID = word(IND, 3, sep = "_"))  # Extracts the third section after two underscores
phase_kir3dl3 <- phase_kir3dl3 %>%
  mutate(IND_ID = word(IND, 3, sep = "_"))  # Extracts the third section after two underscores

merged_df <- left_join(freq_allele, phase_kir3dl2, by = "IND_ID") %>%
  mutate(KIR3DL2 = coalesce(Alleles, KIR3DL2)) %>%
  select(-Alleles, -IND)  # Remove the redundant column if necessary

merged_df <- left_join(merged_df, phase_kir3dl3, by = "IND_ID") %>%
  mutate(KIR3DL3 = coalesce(Alleles, KIR3DL3)) %>%
  select(-Alleles, -IND)  # Remove the redundant column if necessary


write.csv(merged_df, "kirAllelesForAnalysis.csv")

#### KIR3DL1/KIR3DL2 FUSION DETECTION (Deprecated)
# print('Running KIR3DL1/KIR3DL2 fusion detection')
# predictedCopyNumber <- read.csv("predictedCopyNumberFrame.csv")
# kffCountFrame <- read.csv("kffCountFrame.csv")
# # Extract the relevant columns
# kff_selected <- kffCountFrame[, c("X", "X.KIR3DL12.v", "X.KIR3DL12.v_rc")]
# pred_selected <- predictedCopyNumber[, c("X", "KIR3DL1", "KIR3DL2", "KIR3DS1")]
# # Merge the two dataframes using 'X' as the key
# merged_df <- merge(kff_selected, pred_selected, by = "X", all = FALSE)
# merged_df$KIR3DL1S1 <- rowSums(merged_df[, c("KIR3DL1", "KIR3DS1")], na.rm = TRUE)
# # Filter for fusion
# filtered_df <- merged_df[
#   merged_df$`X.KIR3DL12.v` != 0 | 
#   merged_df$`X.KIR3DL12.v_rc` != 0 & 
#   merged_df$KIR3DL2 == 1 & 
#   merged_df$KIR3DL1S1 == 1, 
# ]
# # print(paste('dimension filtered_df:',dim(filtered_df)))
# write.csv(filtered_df, "fusionSamples.csv")
