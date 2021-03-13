import csv
from functools import lru_cache

@lru_cache(None)
def get_tgp_urls():
    with open('../data/tgp_30x.tsv') as csvfile:
        csvr = csv.reader(csvfile,delimiter='\t')
        next(csvr)
        return {x[5]:x[0] for x in csvr}

def get_tgp_url(sample):
    return get_tgp_urls()[sample]


rule retrieve_fastq:
    output:
        "../data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell:
        """wget -O {output} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"""



rule retrieve_cram:
    input:
      "../data/kir.bed"
    params:
        lambda wildcards: get_tgp_url(wildcards.sample)
    output:
        cramf="../input/{sample}_kir.cram",
    conda:
        "../../envs/ping.yaml"                
    shell:
        "samtools view -L {input} -o {output.cramf} {params}"
