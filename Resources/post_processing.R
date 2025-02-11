# This script
#   1) Selects only the first of the allele calls in each column (the pair before a space)
#      to resolve ambiguity in calls. The pair that remains should be the one with generally 
#      lower allele group number. Check that this is how the allele pairs are being ordered, 
#      otherwise the ambiguity resolution should be done manually
#   2) Formats finalAlleleCalls.csv into copy 1 and 2 for each gene of each sample
#   3) Incorporates calls from iterAlleleCalls.csv where there are unresolved calls

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
print(args)

setwd(args) #Set this to your PING results directory
finalAlleleCallsProcessed <- read.csv("finalAlleleCalls.csv")


finalAlleleCallsProcessed$KIR3DP1 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR3DP1)
finalAlleleCallsProcessed$KIR2DP1 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DP1)
finalAlleleCallsProcessed$KIR2DS2 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DS2)
finalAlleleCallsProcessed$KIR2DL4 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DL4)
finalAlleleCallsProcessed$KIR3DL3 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR3DL3)
finalAlleleCallsProcessed$KIR3DL2 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR3DL2)
finalAlleleCallsProcessed$KIR2DS4 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DS4)
finalAlleleCallsProcessed$KIR2DL1 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DL1)
finalAlleleCallsProcessed$KIR2DS1 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DS1)
finalAlleleCallsProcessed$KIR2DL5 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DL5)
finalAlleleCallsProcessed$KIR2DL23 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DL23)
finalAlleleCallsProcessed$KIR2DS35 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR2DS35)
finalAlleleCallsProcessed$KIR3DL1S1 <-  sub(" .*", "", finalAlleleCallsProcessed$KIR3DL1S1)

finalAlleleCallsProcessed$KIR3DP1_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DP1), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR3DP1_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DP1), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DP1_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DP1), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DP1_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DP1), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DS2_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS2), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DS2_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS2), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DL4_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL4), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DL4_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL4), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR3DL3_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DL3), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR3DL3_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DL3), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR3DL2_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DL2), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR3DL2_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DL2), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DS4_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS4), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DS4_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS4), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DL1_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL1), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DL1_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL1), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DS1_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS1), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DS1_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS1), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DL5_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL5), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DL5_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL5), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DL23_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL23), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DL23_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DL23), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR2DS35_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS35), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR2DS35_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR2DS35), "\\+"), "[", 1)

finalAlleleCallsProcessed$KIR3DL1S1_1 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DL1S1), "\\+"), "[", 2)
finalAlleleCallsProcessed$KIR3DL1S1_2 <- lapply(strsplit(as.character(finalAlleleCallsProcessed$KIR3DL1S1), "\\+"), "[", 1)


finalAlleleCallsProcessed$KIR3DP1 <-  NULL
finalAlleleCallsProcessed$KIR2DP1 <-  NULL
finalAlleleCallsProcessed$KIR2DS2 <-  NULL
finalAlleleCallsProcessed$KIR2DL4 <-  NULL
finalAlleleCallsProcessed$KIR3DL3 <-  NULL
finalAlleleCallsProcessed$KIR3DL2 <-  NULL
finalAlleleCallsProcessed$KIR2DS4 <-  NULL
finalAlleleCallsProcessed$KIR2DL1 <-  NULL
finalAlleleCallsProcessed$KIR2DS1 <-  NULL
finalAlleleCallsProcessed$KIR2DL5 <-  NULL
finalAlleleCallsProcessed$KIR2DL23 <-  NULL
finalAlleleCallsProcessed$KIR2DS35 <-  NULL
finalAlleleCallsProcessed$KIR3DL1S1 <-  NULL



iteralleleprocess <- read.csv("iterAlleleCalls.csv", row.names = NULL)
#iteralleleprocess$KIR3DP1_resolved <- iteralleleprocess$KIR2DS5
iteralleleprocess$X <- iteralleleprocess$row.names
iteralleleprocess$row.names <- NULL

processedanditer <- merge(finalAlleleCallsProcessed, iteralleleprocess, by = "X", all.x = TRUE)
processedanditer <- as.data.frame(processedanditer)

processedanditer <- apply(processedanditer,2,as.character)
# write.csv(processedanditer, "processedanditer.csv")


finalAlleleCallsProcessed <- as.data.frame(finalAlleleCallsProcessed)
# write.csv(finalAlleleCallsProcessed, "finalAlleleCallsProcessed.csv")

processedanditer_test <- processedanditer

processedanditer_test <- as.data.frame(processedanditer_test)
processedanditer_test[] <- lapply(processedanditer_test, as.character)


processedanditer_test[is.na(processedanditer_test)] <- 0

processedanditer_test$KIR2DS5_resolved <- processedanditer_test$KIR2DS5
processedanditer_test$KIR2DL3_resolved <- processedanditer_test$KIR2DL3
processedanditer_test$KIR2DP1_resolved <- processedanditer_test$KIR2DP1
processedanditer_test$KIR2DS3_resolved <- processedanditer_test$KIR2DS3
processedanditer_test$KIR2DS2_resolved <- processedanditer_test$KIR2DS2
processedanditer_test$KIR2DL4_resolved <- processedanditer_test$KIR2DL4
processedanditer_test$KIR3DL3_resolved <- processedanditer_test$KIR3DL3
processedanditer_test$KIR2DL2_resolved <- processedanditer_test$KIR2DL2
processedanditer_test$KIR3DL2_resolved <- processedanditer_test$KIR3DL2
processedanditer_test$KIR2DS4_resolved <- processedanditer_test$KIR2DS4
processedanditer_test$KIR2DL1_resolved <- processedanditer_test$KIR2DL1
processedanditer_test$KIR2DS1_resolved <- processedanditer_test$KIR2DS1
processedanditer_test$KIR2DL5_resolved <- processedanditer_test$KIR2DL5
processedanditer_test$KIR3DL1_resolved <- processedanditer_test$KIR3DL1
processedanditer_test$KIR3DS1_resolved <- processedanditer_test$KIR3DS1
processedanditer_test$KIR3DP1_resolved <- processedanditer_test$KIR3DP1 # RS JW ADD 3/26

processedanditer_test$KIR3DP1_1 <- ifelse(processedanditer_test$KIR3DP1_1 == "KIR3DP1*unresolved", processedanditer_test$KIR3DP1_resolved, processedanditer_test$KIR3DP1_1)
processedanditer_test$KIR3DP1_2 <- ifelse(processedanditer_test$KIR3DP1_2 == "KIR3DP1*unresolved", processedanditer_test$KIR3DP1_resolved, processedanditer_test$KIR3DP1_2)

# RS JW COMMENT 3/26
# processedanditer_test$KIR2DS5_resolved <- processedanditer_test$KIR2DL3
# processedanditer_test$KIR2DL3_resolved <- processedanditer_test$KIR2DP1
# processedanditer_test$KIR2DP1_resolved <- processedanditer_test$KIR2DS3
# processedanditer_test$KIR2DS3_resolved <- processedanditer_test$KIR2DS2
# processedanditer_test$KIR2DS2_resolved <- processedanditer_test$KIR2DL4
# processedanditer_test$KIR2DL4_resolved <- processedanditer_test$KIR3DL3
# processedanditer_test$KIR3DL3_resolved <- processedanditer_test$KIR3DL1
# processedanditer_test$KIR2DL2_resolved <- processedanditer_test$KIR3DL2
# processedanditer_test$KIR3DL2_resolved <- processedanditer_test$KIR2DS4
# processedanditer_test$KIR2DS4_resolved <- processedanditer_test$KIR2DL1
# processedanditer_test$KIR2DL1_resolved <- processedanditer_test$KIR2DS1
# processedanditer_test$KIR2DS1_resolved <- processedanditer_test$KIR2DL5
# processedanditer_test$KIR2DL5_resolved <- processedanditer_test$KIR2DL23
# processedanditer_test$KIR3DL1_resolved <- processedanditer_test$KIR3DS1
# processedanditer_test$KIR3DS1_resolved <- processedanditer_test$KIR2DL2





processedanditer_test$KIR2DS35_1 <- ifelse(processedanditer_test$KIR2DS35_1 == "KIR2DS5*unresolved", processedanditer_test$KIR2DS5_resolved, processedanditer_test$KIR2DS35_1)
processedanditer_test$KIR2DS35_2 <- ifelse(processedanditer_test$KIR2DS35_2 == "KIR2DS5*unresolved", processedanditer_test$KIR2DS5_resolved, processedanditer_test$KIR2DS35_2)
processedanditer_test$KIR2DS35_1 <- ifelse(processedanditer_test$KIR2DS35_1 == "KIR2DS3*unresolved", processedanditer_test$KIR2DS3_resolved, processedanditer_test$KIR2DS35_1)
processedanditer_test$KIR2DS35_2 <- ifelse(processedanditer_test$KIR2DS35_2 == "KIR2DS3*unresolved", processedanditer_test$KIR2DS3_resolved, processedanditer_test$KIR2DS35_2)


processedanditer_test$KIR2DL23_1 <- ifelse(processedanditer_test$KIR2DL23_1 == "KIR2DL3*unresolved", processedanditer_test$KIR2DL3_resolved, processedanditer_test$KIR2DL23_1)
processedanditer_test$KIR2DL23_2 <- ifelse(processedanditer_test$KIR2DL23_2 == "KIR2DL3*unresolved", processedanditer_test$KIR2DL3_resolved, processedanditer_test$KIR2DL23_2)
processedanditer_test$KIR2DL23_1 <- ifelse(processedanditer_test$KIR2DL23_1 == "KIR2DL2*unresolved", processedanditer_test$KIR2DL2_resolved, processedanditer_test$KIR2DL23_1)
processedanditer_test$KIR2DL23_2 <- ifelse(processedanditer_test$KIR2DL23_2 == "KIR2DL2*unresolved", processedanditer_test$KIR2DL2_resolved, processedanditer_test$KIR2DL23_2)


processedanditer_test$KIR3DL1S1_1 <- ifelse(processedanditer_test$KIR3DL1S1_1 == "KIR3DL1*unresolved", processedanditer_test$KIR3DL1_resolved, processedanditer_test$KIR3DL1S1_1)
processedanditer_test$KIR3DL1S1_2 <- ifelse(processedanditer_test$KIR3DL1S1_2 == "KIR3DL1*unresolved", processedanditer_test$KIR3DL1_resolved, processedanditer_test$KIR3DL1S1_2)
processedanditer_test$KIR3DL1S1_1 <- ifelse(processedanditer_test$KIR3DL1S1_1 == "KIR3DS1*unresolved", processedanditer_test$KIR3DS1_resolved, processedanditer_test$KIR3DL1S1_1)
processedanditer_test$KIR3DL1S1_2 <- ifelse(processedanditer_test$KIR3DL1S1_2 == "KIR3DS1*unresolved", processedanditer_test$KIR3DS1_resolved, processedanditer_test$KIR3DL1S1_2)


processedanditer_test$KIR2DP1_1 <- ifelse(processedanditer_test$KIR2DP1_1 == "KIR2DP1*unresolved", processedanditer_test$KIR2DP1_resolved, processedanditer_test$KIR2DP1_1)
processedanditer_test$KIR2DP1_2 <- ifelse(processedanditer_test$KIR2DP1_2 == "KIR2DP1*unresolved", processedanditer_test$KIR2DP1_resolved, processedanditer_test$KIR2DP1_2)


processedanditer_test$KIR2DS2_1 <- ifelse(processedanditer_test$KIR2DS2_1 == "KIR2DS2*unresolved", processedanditer_test$KIR2DS2_resolved, processedanditer_test$KIR2DS2_1)
processedanditer_test$KIR2DS2_2 <- ifelse(processedanditer_test$KIR2DS2_2 == "KIR2DS2*unresolved", processedanditer_test$KIR2DS2_resolved, processedanditer_test$KIR2DS2_2)


processedanditer_test$KIR2DL4_1 <- ifelse(processedanditer_test$KIR2DL4_1 == "KIR2DL4*unresolved", processedanditer_test$KIR2DL4_resolved, processedanditer_test$KIR2DL4_1)
processedanditer_test$KIR2DL4_2 <- ifelse(processedanditer_test$KIR2DL4_2 == "KIR2DL4*unresolved", processedanditer_test$KIR2DL4_resolved, processedanditer_test$KIR2DL4_2)


processedanditer_test$KIR3DL3_1 <- ifelse(processedanditer_test$KIR3DL3_1 == "KIR3DL3*unresolved", processedanditer_test$KIR3DL3_resolved, processedanditer_test$KIR3DL3_1)
processedanditer_test$KIR3DL3_2 <- ifelse(processedanditer_test$KIR3DL3_2 == "KIR3DL3*unresolved", processedanditer_test$KIR3DL3_resolved, processedanditer_test$KIR3DL3_2)



processedanditer_test$KIR3DL2_1 <- ifelse(processedanditer_test$KIR3DL2_1 == "KIR3DL2*unresolved", processedanditer_test$KIR3DL2_resolved, processedanditer_test$KIR3DL2_1)
processedanditer_test$KIR3DL2_2 <- ifelse(processedanditer_test$KIR3DL2_2 == "KIR3DL2*unresolved", processedanditer_test$KIR3DL2_resolved, processedanditer_test$KIR3DL2_2)


processedanditer_test$KIR2DS4_1 <- ifelse(processedanditer_test$KIR2DS4_1 == "KIR2DS4*unresolved", processedanditer_test$KIR2DS4_resolved, processedanditer_test$KIR2DS4_1)
processedanditer_test$KIR2DS4_2 <- ifelse(processedanditer_test$KIR2DS4_2 == "KIR2DS4*unresolved", processedanditer_test$KIR2DS4_resolved, processedanditer_test$KIR2DS4_2)

processedanditer_test$KIR2DL1_1 <- ifelse(processedanditer_test$KIR2DL1_1 == "KIR2DL1*unresolved", processedanditer_test$KIR2DL1_resolved, processedanditer_test$KIR2DL1_1)
processedanditer_test$KIR2DL1_2 <- ifelse(processedanditer_test$KIR2DL1_2 == "KIR2DL1*unresolved", processedanditer_test$KIR2DL1_resolved, processedanditer_test$KIR2DL1_2)


processedanditer_test$KIR2DS1_1 <- ifelse(processedanditer_test$KIR2DS1_1 == "KIR2DS1*unresolved", processedanditer_test$KIR2DS1_resolved, processedanditer_test$KIR2DS1_1)
processedanditer_test$KIR2DS1_2 <- ifelse(processedanditer_test$KIR2DS1_2 == "KIR2DS1*unresolved", processedanditer_test$KIR2DS1_resolved, processedanditer_test$KIR2DS1_2)


processedanditer_test$KIR2DL5_1 <- ifelse(processedanditer_test$KIR2DL5_1 == "KIR2DL5*unresolved", processedanditer_test$KIR2DL5_resolved, processedanditer_test$KIR2DL5_1)
processedanditer_test$KIR2DL5_2 <- ifelse(processedanditer_test$KIR2DL5_2 == "KIR2DL5*unresolved", processedanditer_test$KIR2DL5_resolved, processedanditer_test$KIR2DL5_2)


final_alleles_completely_processed <- select(processedanditer_test, c("X"        ,       "KIR3DP1_1"   ,     "KIR3DP1_2"    ,    "KIR2DP1_1"      ,  "KIR2DP1_2"    ,    "KIR2DS2_1"    ,    "KIR2DS2_2"    ,    "KIR2DL4_1"    ,    "KIR2DL4_2"  ,     
                                                                      "KIR3DL3_1"    ,    "KIR3DL3_2"   ,     "KIR3DL2_1"    ,    "KIR3DL2_2"    ,    "KIR2DS4_1"    ,    "KIR2DS4_2"    ,    "KIR2DL1_1"    ,    "KIR2DL1_2"   ,     "KIR2DS1_1"  ,     
                                                                      "KIR2DS1_2"   ,     "KIR2DL5_1"    ,    "KIR2DL5_2"     ,   "KIR2DL23_1"    ,   "KIR2DL23_2"    ,   "KIR2DS35_1"   ,    "KIR2DS35_2"    ,   "KIR3DL1S1_1"    ,  "KIR3DL1S1_2"))


stop_pos = 13

final_alleles_completely_processed$KIR3DP1_1 <- substr(final_alleles_completely_processed$KIR3DP1_1, start = 1, stop = 13)

final_alleles_completely_processed_2 <- final_alleles_completely_processed

final_alleles_completely_processed_2$KIR3DP1_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR3DP1_2)
final_alleles_completely_processed_2$KIR3DP1_2 <- substr(final_alleles_completely_processed_2$KIR3DP1_2, start = 1, stop = stop_pos)


final_alleles_completely_processed_2$KIR2DS35_1 <- substr(final_alleles_completely_processed_2$KIR2DS35_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DS35_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DS35_2)
final_alleles_completely_processed_2$KIR2DS35_2 <- substr(final_alleles_completely_processed_2$KIR2DS35_2, start = 1, stop = stop_pos)


final_alleles_completely_processed_2$KIR2DL23_1 <- substr(final_alleles_completely_processed_2$KIR2DL23_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DL23_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DL23_2)
final_alleles_completely_processed_2$KIR2DL23_2 <- substr(final_alleles_completely_processed_2$KIR2DL23_2, start = 1, stop = stop_pos)


final_alleles_completely_processed_2$KIR3DL1S1_1 <- substr(final_alleles_completely_processed_2$KIR3DL1S1_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR3DL1S1_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR3DL1S1_2)
final_alleles_completely_processed_2$KIR3DL1S1_2 <- substr(final_alleles_completely_processed_2$KIR3DL1S1_2, start = 1, stop = stop_pos)


final_alleles_completely_processed_2$KIR2DP1_1 <- substr(final_alleles_completely_processed_2$KIR2DP1_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DP1_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DP1_2)
final_alleles_completely_processed_2$KIR2DP1_2 <- substr(final_alleles_completely_processed_2$KIR2DP1_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR2DS2_1 <- substr(final_alleles_completely_processed_2$KIR2DS2_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DS2_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DS2_2)
final_alleles_completely_processed_2$KIR2DS2_2 <- substr(final_alleles_completely_processed_2$KIR2DS2_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR2DL4_1 <- substr(final_alleles_completely_processed_2$KIR2DL4_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DL4_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DL4_2)
final_alleles_completely_processed_2$KIR2DL4_2 <- substr(final_alleles_completely_processed_2$KIR2DL4_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR3DL3_1 <- substr(final_alleles_completely_processed_2$KIR3DL3_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR3DL3_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR3DL3_2)
final_alleles_completely_processed_2$KIR3DL3_2 <- substr(final_alleles_completely_processed_2$KIR3DL3_2, start = 1, stop = stop_pos)




final_alleles_completely_processed_2$KIR3DL2_1 <- substr(final_alleles_completely_processed_2$KIR3DL2_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR3DL2_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR3DL2_2)
final_alleles_completely_processed_2$KIR3DL2_2 <- substr(final_alleles_completely_processed_2$KIR3DL2_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR2DS4_1 <- substr(final_alleles_completely_processed_2$KIR2DS4_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DS4_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DS4_2)
final_alleles_completely_processed_2$KIR2DS4_2 <- substr(final_alleles_completely_processed_2$KIR2DS4_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR2DL1_1 <- substr(final_alleles_completely_processed_2$KIR2DL1_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DL1_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DL1_2)
final_alleles_completely_processed_2$KIR2DL1_2 <- substr(final_alleles_completely_processed_2$KIR2DL1_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR2DS1_1 <- substr(final_alleles_completely_processed_2$KIR2DS1_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DS1_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DS1_2)
final_alleles_completely_processed_2$KIR2DS1_2 <- substr(final_alleles_completely_processed_2$KIR2DS1_2, start = 1, stop = stop_pos)



final_alleles_completely_processed_2$KIR2DL5_1 <- substr(final_alleles_completely_processed_2$KIR2DL5_1, start = 1, stop = 13)
final_alleles_completely_processed_2$KIR2DL5_2 <- sub(".*\\+", "", final_alleles_completely_processed_2$KIR2DL5_2)
final_alleles_completely_processed_2$KIR2DL5_2 <- substr(final_alleles_completely_processed_2$KIR2DL5_2, start = 1, stop = stop_pos)


write.csv(final_alleles_completely_processed_2, "kirAllelesForAnalysis.csv")


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
