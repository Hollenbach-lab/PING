#!/usr/bin/env Rscript

# This script is written to update the IPD-KIR resources that's released every year

msf_dir = 'raw_msf_2025'

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(methods)
})

cat("\nPING IPD-KIR resource rebuild\n\n")

## ------------------------------------------------------------------
## Locate repo root
## ------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
fileArg <- grep("^--file=", args, value = TRUE)

if (length(fileArg) > 0) {
  scriptPath <- normalizePath(
    sub("^--file=", "", fileArg[1]),
    mustWork = TRUE
  )

  repoRoot <- normalizePath(
    file.path(dirname(scriptPath), ".."),
    mustWork = TRUE
  )
} else {
  repoRoot <- normalizePath(".", mustWork = TRUE)
}

cat("Repo root:", repoRoot, "\n")

## ------------------------------------------------------------------
## Source required files
## ------------------------------------------------------------------

source(file.path(repoRoot, "Resources", "general_functions.R"))
source(file.path(repoRoot, "Resources", "ping_allele.R"))
source(file.path(repoRoot, "Resources", "build_ipdkir_functions.R"))

## genotype_alignment_functions.R normally gets is_nuc() indirectly
## from other runtime files, so define it here for standalone use.
is_nuc <- function(chr) {
  as.character(chr) %in% c("A", "T", "C", "G", ".")
}

## ------------------------------------------------------------------
## Settings and paths
## ------------------------------------------------------------------

threads <- as.integer(Sys.getenv("PING_THREADS", "8"))

bowtie2Build <- Sys.which("bowtie2-build")
if (bowtie2Build == "") {
  stop("bowtie2-build not found in PATH")
}

kirLocusList <- kir.locus.vect

copiedMsfDirectory <- file.path(
  repoRoot,
  "Resources",
  "ipdkir_resources",
  msf_dif
)

referenceResourceDirectory <- file.path(
  repoRoot,
  "Resources",
  "ipdkir_resources",
  "reference_resources"
)

snpDFDirectory <- file.path(
  repoRoot,
  "Resources",
  "genotype_resources",
  "SNP_files"
)

annotatedAlleleDirectory <- file.path(
  repoRoot,
  "Resources",
  "genotype_resources",
  "extended_SNP_files"
)

gcResourceDirectory <- file.path(
  repoRoot,
  "Resources",
  "gc_resources"
)

filledKirReferenceDirectory <- file.path(
  gcResourceDirectory,
  "filled_kir_reference"
)

UTRextPath <- file.path(
  repoRoot,
  "Resources",
  "genotype_resources",
  "KIR_UTR_ext.fasta"
)

dir.create(referenceResourceDirectory, recursive = TRUE, showWarnings = FALSE)
dir.create(snpDFDirectory, recursive = TRUE, showWarnings = FALSE)
dir.create(annotatedAlleleDirectory, recursive = TRUE, showWarnings = FALSE)
dir.create(filledKirReferenceDirectory, recursive = TRUE, showWarnings = FALSE)

UTRextList <- general.read_fasta(UTRextPath)

cat("MSF directory:", copiedMsfDirectory, "\n")
cat("SNP directory:", snpDFDirectory, "\n")
cat("Threads:", threads, "\n\n")

## ------------------------------------------------------------------
## Sanity check MSF files
## ------------------------------------------------------------------

expectedRawFiles <- paste0(kirLocusList, "_raw.msf")
actualRawFiles <- list.files(
  copiedMsfDirectory,
  pattern = "_raw\\.msf$"
)

missingRawFiles <- setdiff(expectedRawFiles, actualRawFiles)

if (length(missingRawFiles) > 0) {
  stop(
    "Missing expected raw MSF files:\n",
    paste(missingRawFiles, collapse = "\n")
  )
}

problemRawFiles <- intersect(
  actualRawFiles,
  c("KIR2DL5A_raw.msf", "KIR2DL5B_raw.msf")
)

if (length(problemRawFiles) > 0) {
  stop(
    "Remove these files from new_msf before running because PING expects one KIR2DL5 file:\n",
    paste(problemRawFiles, collapse = "\n")
  )
}

## ------------------------------------------------------------------
## Rebuild resources
## ------------------------------------------------------------------

cat("\n----- Building reference object from raw MSF files -----\n")

old.locusRefList <- general.initialize_locus_ref_object()

old.locusRefList <- initLocusRef.read_raw_msf(
  old.locusRefList,
  copiedMsfDirectory
)


old.locusRefList <- initLocusRef.create_bed(
  old.locusRefList,
  referenceResourceDirectory,
  kirLocusFeatureNameList,
  writeBed = TRUE
)

cat("\n----- Creating allele SNP resources -----\n")

alleleDFPathList <- allele.create_allele_resources(
  old.locusRefList,
  snpDFDirectory
)

## Initialize a fresh locusRefList for the final object
locusRefList <- general.initialize_locus_ref_object()

filled.snpDFList <- new.initLocusRef.read_snp_df(
  locusRefList,
  snpDFDirectory
)

filled.snpDFList <- new.initLocusRef.extend_5UTR(
  filled.snpDFList,
  UTRextList
)

filled.snpDFList <- new.initLocusRef.extend_3UTR(
  filled.snpDFList,
  UTRextList
)

remove(old.locusRefList)

locusRefList <- new.initLocusRef.snpDFtoLocusRefAlleleSeq(
  filled.snpDFList,
  locusRefList
)

locusRefList <- new.initLocusRef.snpDFtoLocusRefBed(
  filled.snpDFList,
  locusRefList,
  kirLocusFeatureNameList
)

cat("\n----- Saving RDS resources -----\n")

saveRDS(
  filled.snpDFList,
  file = file.path(repoRoot, "Resources", "filled.snpDFList.rds")
)

saveRDS(
  locusRefList,
  file = file.path(repoRoot, "Resources", "locusRefList.rds")
)

utils::zip(
  zipfile = file.path(repoRoot, "Resources", "locusRefList.rds.zip"),
  files = file.path(repoRoot, "Resources", "locusRefList.rds")
)

cat("Saved Resources/filled.snpDFList.rds\n")
cat("Saved Resources/locusRefList.rds\n")

## Update gc_allele_reference.csv if naming change
## Example: the csv previously had 2DP1*004, but the new allele has 2DP1*0040101 and 2DP4*0040102
## This section updates gc_allele_reference.csv from 2DP1*004 to 2DP1*0040101

gcRefPath <- file.path(
  repoRoot,
  "Resources",
  "genotype_resources",
  "gc_allele_reference.csv"
)

update_gc_reference_csv(
  referenceCSVPath = gcRefPath,
  filled.snpDFList = filled.snpDFList,
  outputCSVPath = gcRefPath,
  logPath = file.path(
    repoRoot,
    "Resources",
    "genotype_resources",
    "gc_allele_reference_update.log"
  ),
  strict = TRUE
)

cat("\n----- Writing extended allele SNP CSV files -----\n")

for (locus in kirLocusList) {
  write.csv(
    filled.snpDFList[[locus]],
    file.path(
      annotatedAlleleDirectory,
      paste0(locus, "_alleleSNPs.csv")
    )
  )
}

## ------------------------------------------------------------------
## Build full filled FASTA + Bowtie2 index
## ------------------------------------------------------------------

cat("\n----- Writing KIR_gen_onelines_filled.fasta -----\n")

kirReferenceFasta <- file.path(
  filledKirReferenceDirectory,
  "KIR_gen_onelines_filled.fasta"
)

kirReferenceIndex <- file.path(
  filledKirReferenceDirectory,
  "KIR_gen_onelines_filled"
)

fastaCon <- file(kirReferenceFasta, "w")

for (locus in names(filled.snpDFList)) {
  for (alleleName in rownames(filled.snpDFList[[locus]])) {
    alleleStr <- paste0(
      filled.snpDFList[[locus]][alleleName, ],
      collapse = ""
    )

    alleleStr <- gsub(".", "", alleleStr, fixed = TRUE)
    alleleStr <- gsub("*", "N", alleleStr, fixed = TRUE)

    general.write_fasta(
      fastaCon,
      alleleName,
      alleleStr
    )
  }
}

close(fastaCon)

system2(
  bowtie2Build,
  c(
    kirReferenceFasta,
    kirReferenceIndex,
    "--quiet",
    "--threads",
    as.character(threads)
  )
)

cat("Built Bowtie2 index:", kirReferenceIndex, "\n")

## ------------------------------------------------------------------
## Build compact filled FASTA + Bowtie2 index
## ------------------------------------------------------------------

cat("\n----- Writing KIR_compact_filled.fasta -----\n")

kirCompactFasta <- file.path(
  filledKirReferenceDirectory,
  "KIR_compact_filled.fasta"
)

kirCompactIndex <- file.path(
  filledKirReferenceDirectory,
  "KIR_compact_filled"
)

## The original PING code uses the existing compact FASTA as a template
## to decide which alleles belong in the compact reference.
if (file.exists(kirCompactFasta)) {
  test.fa <- general.read_fasta(kirCompactFasta)

  alleleVect <- unlist(
    lapply(names(filled.snpDFList), function(x) {
      unlist(
        sapply(names(test.fa), function(y) {
          grep(
            y,
            rownames(filled.snpDFList[[x]]),
            fixed = TRUE,
            value = TRUE
          )
        })
      )
    }),
    use.names = FALSE
  )
} else {
  warning(
    "Existing KIR_compact_filled.fasta not found. ",
    "Using all alleles for compact FASTA."
  )

  alleleVect <- unlist(
    lapply(filled.snpDFList, rownames),
    use.names = FALSE
  )
}

fastaCon <- file(kirCompactFasta, "w")

for (locus in names(filled.snpDFList)) {
  for (alleleName in rownames(filled.snpDFList[[locus]])) {
    if (alleleName %in% alleleVect) {
      alleleStr <- paste0(
        filled.snpDFList[[locus]][alleleName, ],
        collapse = ""
      )

      alleleStr <- gsub(".", "", alleleStr, fixed = TRUE)
      alleleStr <- gsub("*", "N", alleleStr, fixed = TRUE)

      general.write_fasta(
        fastaCon,
        alleleName,
        alleleStr
      )
    }
  }
}

close(fastaCon)

system2(
  bowtie2Build,
  c(
    kirCompactFasta,
    kirCompactIndex,
    "--quiet",
    "--threads",
    as.character(threads)
  )
)

cat("Built Bowtie2 index:", kirCompactIndex, "\n")

cat("\nIPD-KIR resource rebuild complete.\n")