rule bam_extract:
    """Download the subset of the 1000 genomes CRAM file corresponding to the kir regions"""
    input:
        bedf="../data/kir_regions.bed",
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
    output:
        bamf="../input/tgp/bam/{sample}_kir.bam",
    conda:
        "../../envs/ping.yaml"
    resources:
        mem_mb=8000,
        time=str(timedelta(hours=3))
    threads:
        4
    shell:
        "samtools view --reference {input.ref_fq} --threads {threads} -M -L {input.bedf} -o {output.bamf} ../input/{pref}/bam/{sample}.bam"

rule index_cram:
  """index the subset CRAM file"""
    input:
        "../input/{pref}/cram/{sample}_kir.cram"
    output:
        "../input/{pref}/cram/{sample}_kir.cram.crai"
    conda:
        "../../envs/ping.yaml"
    resources:
        mem_mb=2000,
        time=str(timedelta(minutes=30))
    threads:
        4
    shell:
        """samtools index -@ {threads} {input}"""


rule index_bam:
  """index the subset BAM file"""
    input:
        "../input/{pref}/bam/{sample}_kir.bam"
    output:
        "../input/{pref}/bam/{sample}_kir.bam.bai"
    conda:
        "../../envs/ping.yaml"
    resources:
        mem_mb=2000,
        time=str(timedelta(minutes=30))
    threads:
        4
    shell:
        """samtools index -@ {threads} {input}"""



rule cram2fastq:
    """convert a cram to fastq suitable for PING"""
    input:
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cramf="../input/{pref}/cram/{sample}_kir.cram",
        crami="../input/{pref}/cram/{sample}_kir.cram.crai",
    output:
        fq1="../input/{pref}/fastq/{sample}_KIR_1.fastq.gz",
        fq2="../input/{pref}/fastq/{sample}_KIR_2.fastq.gz"
    threads: 4
    conda:
        "../../envs/ping.yaml"        
    shell:
        """bazam -Dsamjdk.reference_fasta={input.ref_fq} -Xmx4g  -bam {input.cramf} -r1 {output.fq1} -r2 {output.fq2}"""


rule bam2fastq:
    """convert a bam to fastq suitable for PING"""
    input:
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cramf="../input/{pref}/bam/{sample}_kir.bam",
        crami="../input/{pref}/bam/{sample}_kir.bam.bai",
    output:
        fq1="../input/{pref}/fastq/{sample}_KIR_1.fastq.gz",
        fq2="../input/{pref}/fastq/{sample}_KIR_2.fastq.gz"
    threads: 4
    conda:
        "../../envs/ping.yaml"        
    shell:
        """bazam -Dsamjdk.reference_fasta={input.ref_fq} -Xmx4g  -bam {input.cramf} -r1 {output.fq1} -r2 {output.fq2}"""




rule run_PING:
  """run PING"""
    input:
        fq1="../input/{pref}/fastq/{sample}_KIR_1.fastq.gz",
        fq2="../input/{pref}/fastq/{sample}_KIR_2.fastq.gz"
    output:
        "../output/{pref}/{sample}.tar.gz",
    resources:
        mem_mb=33000,
        time=str(timedelta(hours=3))
    params:
        workingDirectory="..",
        resultsDirectory="output/{pref}/{sample}",
        samplename="{sample}",
        shortNameDelim="_",
        setup_hetRatio="0.25",
        final_hetRatio="0.25",
        setup_minDP="6",
        final_minDP="6",
        copy_readBoost="T",
        setup_readBoost="T",
        final_readBoost="F",
        readBoost_thresh="2",
        allele_fullAlign="F",
        copy_fullAlign="F",
    conda:
        "../../envs/ping.yaml"
    log:
        "../output/{pref}/{sample}_log.txt"        
    threads:
        8
    script:
        """../scripts/PING_run.R"""
