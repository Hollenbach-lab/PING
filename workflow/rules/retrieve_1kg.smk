import csv
from functools import lru_cache
from urllib.parse import urlparse
from urllib.parse import ParseResult
from pathlib import Path


def make_s3(path):
  s3r = Path("1000genomes") / path
  return ParseResult(scheme='http',netloc='s3.amazonaws.com',path=str(s3r),params='',query='',fragment='').geturl()

@lru_cache(None)
def get_tgp_urls():
    with open('../data/tgp_30x.tsv') as csvfile:
        csvr = csv.reader(csvfile,delimiter='\t')
        return {x[0]:make_s3(x[1]) for x in csvr}

def get_tgp_url(sample):
    return get_tgp_urls()[sample]



rule retrieve_fastq:
  """Download GRCh38 from the 1000 genomes project"""
    output:
        "../data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell:
        """wget -O {output} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"""


rule retrieve_cram:
    """Download the subset of the 1000 genomes CRAM file corresponding to the kir regions"""
    input:
        bedf="../data/kir_regions.bed",
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
    params:
        lambda wildcards: get_tgp_url(wildcards.sample)
    output:
        cramf="../input/tgp/cram/{sample}_kir.cram",
    conda:
        "../../envs/ping.yaml"
    resources:
        mem_mb=8000,
        time=str(timedelta(hours=3))
    threads:
        4
    shell:
        "samtools view --reference {input.ref_fq} --threads {threads} -M -L {input.bedf} -o {output.cramf} {params}"
