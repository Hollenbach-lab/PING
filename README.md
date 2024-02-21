# PING (Pushing Immunogenetics to the Next Generation)
An R-based bioinformatic pipeline to determine killer-cell immunoglobulin-like receptor (KIR) copy number and high-resolution genotypes from short-read sequencing data.

## Data compatibility
Paired-end KIR targeted sequencing data

## Setting up pipeline
Download the pipeline code with the following line:
```
git clone https://github.com/Hollenbach-lab/PING.git
cd ./PING
```
## Setting up container
To ensure a reliable run, we have containerized our image in Singularity. We have tested it on Singularity version 3.11.4. The Singularity recipe file, MHConstructor/container/mhconstructor.def, can be built using one of the following commands:
```
cd ./container
sudo singularity build ping.sif ping.def
cd ../
```
If you do not have sudo privilege, you can utilize the `fakeroot` option by Singularity which will let you build the container still, by using the command below instead: 
```
cd ./container
singularity build --fakeroot ping.sif ping.def
cd ../
```

## Running PING
Ensure that you are within the `PING` directory and you can run the entirety of the pipeline using the following command:

```
singularity exec ./container/ping.sif Rscript PING_run.R --fqDirectory <fastq_location> --fastqPattern <fastq_pattern> --resultsDirectory <output_location>
```

#### Environment setup
Change lines 33-39 to fit your data/environment, this does not need to be changed to run the included example dataset.
  - `33 rawFastqDirectory <- 'test_sequence/'` Set to raw sequence directory or extracted fastq directory if extraction has already been performed
  - `34 fastqPattern <- 'fastq'` Use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
  - `35 threads <- 4` Number of threads to use during bowtie2 alignments
  - `36 resultsDirectory <- '3_test_sequence_results/'` Set the results directory, one will be created if it does not already exist (all pipeline output will be recorded here)
  - `37 shortNameDelim <- ''` Set a delimiter to shorten sample ID's (ID will be characters before delim)
  - `38 setup.minDP <- 8` Minimum depth for calling variants used for genotype-matched references (set lower if using low-depth data, the default of 8 should work for most data)
  - `39 final.minDP <- 20` Minimum depth for calling variants used for final genotype determination (set lower if using low-depth data, the default of 20 should work for most data)

## Running included test data
We have included 5 test sequences to run through the pipeline, they are located in the test_sequence/ directory. These samples were picked to cover a range of KIR haplotypes.

The default input settings (lines 33-39) are setup to run on this data without modification.

The copy thresholding functionality will not work well for such a small cohort, but we included preset thresholds for the example dataset.

## Running your own data
Update line 33 `rawFastqDirectory <- 'test_sequence/'` to point to your own data directory.

Update line 34 `fastqPattern <- 'fastq'` to match your data naming if necessary. For example, if your sequencing data is named [SAMPLE_ID]\_R1_fq.gz, you would change line 34 to `fastqPattern <- 'fq'`.

Change line 84 option `use.threshFile=T` to `use.threshFile=F`, this will cause PING to prompt for copy number thresholding. If using Rstudio, copy number plots will be displayed in the 'Plots' panel, these plots can also be found in \[resultsDirectory\]/copyPlots/\[GENE\]\_copy\_number\_plot.html as html files to be opened with a web browser. 

It is recommended to only run the script to this line first if you are using your own data, then set the copy number thresholds, then run the rest of the script. 

## PING output
Copy number output can be found at `[resultsDirectory]/manualCopyNumberFrame.csv`

Genotype output can be found at `[resultsDirectory]/finalAlleleCalls.csv`

Aligned SNP tables can be found in `[resultsDirectory]/alignmentFiles/[sampleID]/iterAlign/`

Copy number graphs can be found in `[resultsDirectory]/copyPlots/`

#### Unresolved genotypes
If PING is unable to perfectly match aligned SNPs to known KIR allele sequences an unresolved call will be produced.

Unresolved genotype information can be found in `[resultsDirectory]/iterAlleleCalls.csv`, where the closest allele match is recorded along with the mismatched SNP information in the following format:
`[closest_matched_allele]$[exon]_[position].[nucleotide]`

Where closest matched allele is the allele genotyping that best matches the aligned SNPs, nucleotide denotes the mismatched nucleotide located at the indicated exon and position within the exon. Multiple mismatched SNPs are connected with the `^` symbol.

## Troubleshooting
Please save a copy of your R Console output and contact me through github or email at wesley.marin@ucsf.edu

# Citations
Please cite:


## PING (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008904)

Marin WM, Dandekar R, Augusto DG, Yusufali T, Heyn B, et al. (2021) High-throughput Interpretation of Killer-cell Immunoglobulin-like Receptor Short-read Sequencing Data with PING. PLOS Computational Biology 17(8): e1008904. https://doi.org/10.1371/journal.pcbi.1008904


## IPD-KIR (https://www.ebi.ac.uk/ipd/kir/)

Robinson J, Waller MJ, Stoehr P, Marsh SGE. IPD-the Immuno Polymorphism Database. Nucleic Acids Research (2005), 331:D523-526


## Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
