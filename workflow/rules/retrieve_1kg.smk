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


def ftp2s3(sample):
    ftp_url = urlparse(get_tgp_urls()[sample])
    ftp_url_path = Path(ftp_url.path)
    rel_p = ftp_url_path.relative_to('/vol1/run')
    rel_pp= rel_p.parents[1]
    new_rel_p=rel_p.relative_to(rel_pp)





rule retrieve_fastq:
    output:
        "../data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell:
        """wget -O {output} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"""



rule retrieve_cram:
    input:
      "../data/kir.bed"
    params:
        lambda wildcards: ftp2s3(wildcards.sample)
    output:
        cramf="../input/{sample}_kir.cram",
    conda:
        "../../envs/ping.yaml"
    resources:
        mem_mb=8000,
        time=str(timedelta(hours=3))
    threads:
        4
    shell:
        "samtools view --threads {threads} -M -L {input} -o {output.cramf} {params}"
