## Snakefile - FASTQC
##
# configfile: "config/config.yaml"


paired_samples = glob_wildcards("FASTQ/rna_paired/merged/{sample_paired}_mergedF.fastq").sample_paired
single_end_samples = glob_wildcards("FASTQ/single_end/{sample}.fastq").sample


## paired_end_fastq_to_fastqc: quality control checks on unmerged paired end sequence data
## .fastq files must be located in a directory named resources/rna_paired/unmerged
## output will be directed to a directory named /FASTQC/{sample}/{sample}


rule paired_end_fastq_to_fastqc:
    input:
        "FASTQ/rna_paired/unmerged/{sample}_1_1.fastq",
        "FASTQ/rna_paired/unmerged/{sample}_1_2.fastq",
        "FASTQ/rna_paired/unmerged/{sample}_2_1.fastq",
        "FASTQ/rna_paired/unmerged/{sample}_2_2.fastq",
    output:
        "FASTQC/{sample}/{sample}_1_1_fastqc.html",
        "FASTQC/{sample}/{sample}_1_1_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/paired_end_fastq_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir = lambda wildcards: "FASTQC/{}".format(wildcards.sample),
        logdir = "logs/paired_end_fastq_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """

## single_end_fastq_to_fastqc: quality control checks on single end sequence data
## .fastq files must be located in a directory named FASTQ/sinle_end/
## output will be directed to a directory named /FASTQC/single_end/{sample}

rule single_end_fastq_to_fastqc:
    input:
        "FASTQ/single_end/{sample}.fastq",
    output:
        "FASTQC/single_end/{sample}/{sample}_fastqc.html",
        "FASTQC/single_end/{sample}/{sample}_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/single_end_fastq_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir = lambda wildcards: "FASTQC/single_end/{}".format(wildcards.sample),
        logdir = "logs/single_end_fastq_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """


## filtered_fwd_rna_to_fastqc: quality control checks on filtered forward sequence data
## .fastq files must be located in a directory named /sortmerna_files and end with *_rRNA_filtered.fastq
## output will be directed to a directory named /FASTQC/sortmerna/filtered/{sample}


rule filtered_fwd_rna_to_fastqc:
    input:
        "results/sortmerna_files/rRNAf/{sample}_fwd.fq"
    output:
        "FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fwd_fastqc.html",
        "FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fwd_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/filtered_fwd_rna_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir=lambda wildcards: "FASTQC/sortmerna/filtered/{}_fwd".format(wildcards.sample),
        logdir = "logs/filtered_fwd_rna_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """

## filtered_rev_rna_to_fastqc: quality control checks on filtered reverse sequence data
## .fastq files must be located in a directory named /sortmerna_files/rRNAf and end with _rev.fq
## output will be directed to a directory named /FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc


rule filtered_rev_rna_to_fastqc:
    input:
        "results/sortmerna_files/rRNAf/{sample}_rev.fq"
    output:
        "FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc.html",
        "FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/filtered_rev_rna_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir=lambda wildcards: "FASTQC/sortmerna/filtered/{}_rev".format(wildcards.sample),
        logdir = "logs/filtered_rev_rna_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """


## filtered_unpaired_rna_to_fastqc: quality control checks on unpaired sorted sequence data
## .fastq files must be located in a directory named /unpaired/rRNAf/ and end with .fq
## output will be directed to a directory named /unpaired/filtered/{sample}/{sample}_fastqc


rule filtered_unpaired_rna_to_fastqc:
    input:
        "results/sortmerna_files/unpaired/rRNAf/{sample}.fq"
    output:
        "FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.html",
        "FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/filtered_unpaired_rna_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir=lambda wildcards: "FASTQC/sortmerna/unpaired/filtered/{}".format(wildcards.sample),
        logdir = "logs/filtered_unpaired_rna_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """


## trimmomatic_paired_rna_to_fastqc: quality control checks on trimmed paired sequence data
## .fq.gz files must be located in a directory named trimmomatic_files and end with *.fq.gz
## output will be directed to a directory named /FASTQC/trimmomatic/{sample}


rule trimmomatic_paired_rna_to_fastqc:
    input:
        "results/trimmomatic_files/{sample}_fwd_p.fq.gz",
        "results/trimmomatic_files/{sample}_fwd_up.fq.gz",
        "results/trimmomatic_files/{sample}_rev_p.fq.gz",
        "results/trimmomatic_files/{sample}_rev_up.fq.gz"
    output:
        "FASTQC/trimmomatic/{sample}/{sample}_fwd_p_fastqc.html",
        "FASTQC/trimmomatic/{sample}/{sample}_fwd_p_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/trimmomatic_paired_rna_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir=lambda wildcards: "FASTQC/trimmomatic/{}".format(wildcards.sample),
        logdir = "logs/trimmomatic_paired_rna_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """


## trimmomatic_unpaired_rna_to_fastqc: quality control checks on trimmed unpaired sequence data
## .fq.gz files must be located in a directory named trimmomatic_files/unpaired and end with *.fq.gz
## output will be directed to a directory named /trimmomatic/unpaired/{sample}/{sample}_fwd


rule trimmomatic_unpaired_rna_to_fastqc:
    input:
        "results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz"
    output:
        "FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.html",
        "FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/trimmomatic_unpaired_rna_to_fastqc/{sample}.log"
    envmodules:
        "/apps/modules/modulefiles/tools/java/12.0.2.lua"
        "/apps/modules/modulefiles/tools/fastqc/0.11.9"
    params:
        outdir=lambda wildcards: "FASTQC/trimmomatic/unpaired/{}".format(wildcards.sample),
        logdir = "logs/trimmomatic_unpaired_rna_to_fastqc/",
        threads=24,
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        fastqc -o {params.outdir} -t {params.threads} {input} 2> {log}
        """