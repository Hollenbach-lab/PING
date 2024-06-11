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
To ensure a reliable run, we have containerized our image in [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html). Please install Singularity on your system before proceeding. We have tested PING to run on Singularity version 3.11.4. 

Once Singularity has been installed, obtain the image using one of these following commands:

1. Build with sudo
```
sudo singularity build ping.sif ping.def
```
2. Build with `--fakeroot`, if sudo is not possible
```
singularity build --fakeroot ping.sif ping.def
```
3. Pull directly from Sylabs
```
singularity pull ping.sif library://rsuseno/rsuseno/ping:latest
```

## Running PING
Ensure that you are within the `PING` directory and you can run the entirety of the pipeline using the following command:

```
singularity exec ping.sif Rscript PING_run.R 
  --fqDirectory <fastq_location> 
  --resultsDirectory <output_location> 
  --fastqPattern <fastq_pattern> 
  --threads <number_of_threads>
```

Listed below are the arguments needed to run PING:
  - `--fqDirectory` Set to raw sequence directory or extracted fastq directory if extraction has already been performed
  - `--resultsDirectory` Set the results directory, one will be created if it does not already exist (all pipeline output will be recorded here)
  - (OPTIONAL) `--fastqPattern` Specify a pattern on the `fqDirectory` to only process  specific samples. For example, if your sequencing data is named [SAMPLE_ID]\_R1_fq.gz, you would change it to '_fq_'. Additionally, you can use '_KIR_' to find already extracted files. (default = '_fastq_')
  - (OPTIONAL) `--threads` Number of threads to use during bowtie2 alignments (default = 4)


## Running included test data
We have included 10 test sequences to run through the pipeline, they are located in the test_sequence/ directory. These samples are meant to test that all the installations were done properly. You can run the following code to execute the test:
```
singularity exec ping.sif Rscript PING_run.R 
  --fqDirectory test_sequence
  --resultsDirectory test_sequence_output 
```


## PING output
Copy number output can be found at `[resultsDirectory]/predictedCopyNumberFrame.csv`

Genotype output can be found at `[resultsDirectory]/finalAlleleCalls.csv`

Aligned SNP tables can be found in `[resultsDirectory]/alignmentFiles/[sampleID]/iterAlign/`

Copy number graphs can be found in `[resultsDirectory]/copyPlots/`

#### Unresolved genotypes
If PING is unable to perfectly match aligned SNPs to known KIR allele sequences an unresolved call will be produced.

Unresolved genotype information can be found in `[resultsDirectory]/iterAlleleCalls.csv`, where the closest allele match is recorded along with the mismatched SNP information in the following format:
`[closest_matched_allele]$[exon]_[position].[nucleotide]`

Where closest matched allele is the allele genotyping that best matches the aligned SNPs, nucleotide denotes the mismatched nucleotide located at the indicated exon and position within the exon. Multiple mismatched SNPs are connected with the `^` symbol.

## Troubleshooting
Please save a copy of your R Console output and contact me through Github or email at wesley.marin@ucsf.edu or rayo.suseno@ucsf.edu

# Citations
Please cite:


## PING (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008904)

Marin WM, Dandekar R, Augusto DG, Yusufali T, Heyn B, et al. (2021) High-throughput Interpretation of Killer-cell Immunoglobulin-like Receptor Short-read Sequencing Data with PING. PLOS Computational Biology 17(8): e1008904. https://doi.org/10.1371/journal.pcbi.1008904


## IPD-KIR (https://www.ebi.ac.uk/ipd/kir/)

Robinson J, Waller MJ, Stoehr P, Marsh SGE. IPD-the Immuno Polymorphism Database. Nucleic Acids Research (2005), 331:D523-526


## Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
