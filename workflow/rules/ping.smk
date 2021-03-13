rule cram2fastq:
    input:
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cramf="../input/{sample}_kir.cram",
    output:
        fq1="../input/{sample}/fastq_1.fastq.gz",
        fq2="../input/{sample}/fastq_2.fastq.gz"
    threads: 12
    container:
        "docker://registry.code.roche.com/knoblauch.nicholas/ping/ping"
    conda:
        "../envs/ping.yaml"        
    shell:
       """samtools sort {input.cramf} --threads {threads} | samtools fastq --reference {input.ref_fq} -1 {output.fq1} -2 {output.fq2} --threads {threads}"""

rule run_PING:
    input:
        fq1="../input/{sample}/fastq_1.fastq.gz",
        fq2="../input/{sample}/fastq_2.fastq.gz"
    output:
        "../output/{sample}.tar.gz"
    params:
        workingDirectory="..",
        rawFastqDirectory="input/{sample}",
        fastqPattern="fastq",
        resultsDirectory="output/{sample}",
        shortNameDelim="_",
        setup_hetRatio="0.25",
        final_hetRatio="0.25",
        setup_minDP="8",
        final_minDP="10",
        copy_readBoost="T",
        setup_readBoost="T",
        final_readBoost="F",
        readBoost_thresh="2",
        allele_fullAlign="F",
        copy_fullAlign="F",
    conda:
        "../envs/ping.yaml"
    container:
        "docker://registry.code.roche.com/knoblauch.nicholas/ping/ping"
    threads:
        26
    script:
        """../scripts/PING_run.R"""
