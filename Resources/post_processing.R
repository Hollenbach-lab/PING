# This script
#   1) Selects only the first of the allele calls in each column (the pair before a space)
#      to resolve ambiguity in calls. The pair that remains should be the one with generally 
#      lower allele group number. Check that this is how the allele pairs are being ordered, 
#      otherwise the ambiguity resolution should be done manually
#   2) Formats finalAlleleCalls.csv into copy 1 and 2 for each gene of each sample
#   3) Incorporates calls from iterAlleleCalls.csv where there are unresolved calls

library(dplyr)

### ALLELE FREQUENCY TO SOLVE AMBIGUOUS GENOTYPE
args <- commandArgs(trailingOnly = TRUE)
print(args)
setwd(args) #Set this to your PING results directory
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
  
  
  
  # ambiguous_count_before = 0
  # for (i in final_allele[[gene]]){
  #   if (grepl(" ", i)) {
  #     ambiguous_count_before = ambiguous_count_before + 1
  #   }  
  # }
  
  # ambiguous_count_after = 0
  # for (i in revisedCalls){
  #   if (length(i)>1) {
  #     ambiguous_count_after = ambiguous_count_after + 1
  #   }
  # }
  # print(paste("Ambiguous", gene,"call before allele frequency = ",ambiguous_count_before))
  # print(paste("Ambiguous", gene,"call after allele frequency = ",ambiguous_count_after))
  # return(list(revisedCalls = revisedCalls, ambiguous_count_before=ambiguous_count_before, ambiguous_count_after=ambiguous_count_after))
  return(revisedCalls)
}

fixed_allele <- cut_allele
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
  fixed_allele[[gene]] <- revised2
  fixed_allele[[gene]] <- sapply(fixed_allele[[gene]], paste, collapse = " ")
}

write.csv(fixed_allele, "kirAllelesForAnalysis.csv")


#### KIR3DL1/KIR3DL2 FUSION DETECTION

print('Running KIR3DL1/KIR3DL2 fusion detection')
predictedCopyNumber <- read.csv("predictedCopyNumberFrame.csv")
kffCountFrame <- read.csv("kffCountFrame.csv")

# Extract the relevant columns
kff_selected <- kffCountFrame[, c("X", "X.KIR3DL12.v", "X.KIR3DL12.v_rc")]
# print(paste('dimension kff_selected:',dim(kff_selected)))
pred_selected <- predictedCopyNumber[, c("X", "KIR3DL1", "KIR3DL2", "KIR3DS1")]
# print(paste('dimension pred_selected:',dim(pred_selected)))
# Merge the two dataframes using 'X' as the key
merged_df <- merge(kff_selected, pred_selected, by = "X", all = FALSE)
merged_df$KIR3DL1S1 <- rowSums(merged_df[, c("KIR3DL1", "KIR3DS1")], na.rm = TRUE)
# print(paste('dimension merged_df:',dim(merged_df)))
# print(merged_df)
# Filter for fusion
filtered_df <- merged_df[
  merged_df$`X.KIR3DL12.v` != 0 | 
  merged_df$`X.KIR3DL12.v_rc` != 0 & 
  merged_df$KIR3DL2 == 1 & 
  merged_df$KIR3DL1S1 == 1, 
]
# print(paste('dimension filtered_df:',dim(filtered_df)))
write.csv(filtered_df, "fusionSamples.csv")
