## Snakefile - FASTQC
##
# configfile: "config/config.yaml"


fastq_samples = glob_wildcards("resources/rna_paired/unmerged/{sample}.fastq").sample
single_end_fastq_samples = glob_wildcards("./FASTQ/single_end/{sample}.fastq").sample
filtered_rna_samples_fwd = glob_wildcards("./results/sortmerna_files/rRNAf/{sample}_fwd.fq").sample
filtered_rna_samples_rev = glob_wildcards("./results/sortmerna_files/rRNAf/{sample}_rev.fq").sample
unpaired_rna_samples = glob_wildcards("./results/sortmerna_files/unpaired/rRNAf/{sample}.fq").sample
trimmed_fwd_paired_rna_samples = glob_wildcards("./results/trimmomatic_files/{sample}_fwd_p.fq.gz").sample
trimmed_unpaired_rna_samples = glob_wildcards ("results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz").sample

## fastq_to_fastqc: quality control checks on unmerged paired end sequence data
## .fastq files must be located in a directory named resources/rna_paired/unmerged
## output will be directed to a directory named /FASTQC/{sample}/{sample}


rule fastq_to_fastqc:
    input:
        "resources/rna_paired/unmerged/{sample}.fastq",
    output:
        "FASTQC/{sample}/{sample}_fastqc.html",
        "FASTQC/{sample}/{sample}_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    params:
        outdir = lambda wildcards: "FASTQC/{}".format(wildcards.sample),
        threads=24,
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """

## single_end_fastq_to_fastqc: quality control checks on single end sequence data
## .fastq files must be located in a directory named FASTQ/sinle_end/
## output will be directed to a directory named /FASTQC/single_end/{sample}

rule single_end_fastq_to_fastqc:
    input:
        "FASTQ/single_end/{sample}.fastq",
    output:
        "FASTQC/single_end/{sample}_fastqc.html",
        "FASTQC/single_end/{sample}_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    params:
        outdir = lambda wildcards: "FASTQC/single_end/{}".format(wildcards.sample),
        threads=24,
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """


## filtered__fwd_rna_to_fastqc: quality control checks on filtered forward sequence data
## .fastq files must be located in a directory named /sortmerna_files and end with *_rRNA_filtered.fastq
## output will be directed to a directory named /FASTQC/sortmerna/filtered/{sample}


rule filtered_fwd_rna_to_fastqc:
    input:
        "results/sortmerna_files/rRNAf/{sample}_fwd.fq"
    output:
        "FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fastqc.html",
        "FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fastqc.zip",
    params:
        outdir=lambda wildcards: "FASTQC/sortmerna/filtered/{}_fwd".format(wildcards.sample),
        threads=24,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """

## filtered_rev_rna_to_fastqc: quality control checks on filtered forward sequence data
## .fastq files must be located in a directory named /sortmerna_files and end with *_rRNA_filtered.fastq
## output will be directed to a directory named /FASTQC/sortmerna/filtered/{sample}


rule filtered_rev_rna_to_fastqc:
    input:
        "results/sortmerna_files/rRNAf/{sample}_rev.fq"
    output:
        "FASTQC/sortmerna/filtered/{sample}_rev/{sample}_fastqc.html",
        "FASTQC/sortmerna/filtered/{sample}_rev/{sample}_fastqc.zip",
    params:
        outdir=lambda wildcards: "FASTQC/sortmerna/filtered/{}_rev".format(wildcards.sample),
        threads=24,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """


## unpaired_rna_to_fastqc: quality control checks on unpaired sequence data
## .fastq files must be located in a directory named /sortmerna_files and end with *_rRNA_filtered.fastq
## output will be directed to a directory named /FASTQC/sortmerna/filtered/{sample}


rule unpaired_rna_to_fastqc:
    input:
        "results/sortmerna_files/unpaired/rRNAf/{sample}.fq"
    output:
        "FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.html",
        "FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.zip",
    params:
        outdir=lambda wildcards: "FASTQC/sortmerna/unpaired/filtered/{}".format(wildcards.sample),
        threads=24,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """


## trimmomatic_paired_rna_to_fastqc: quality control checks on trimmed sequence data
## .fq.gz files must be located in a directory named Trimmomatic_files and end with *.fq.gz
## output will be directed to a directory named /FASTQC/trimmomatic/{sample}


rule trimmomatic_paired_rna_to_fastqc:
    input:
        "results/trimmomatic_files/{sample}_fwd_p.fq.gz",
        "results/trimmomatic_files/{sample}_fwd_up.fq.gz",
        "results/trimmomatic_files/{sample}_rev_p.fq.gz",
        "results/trimmomatic_files/{sample}_rev_up.fq.gz"
    output:
        "FASTQC/trimmomatic/{sample}/{sample}_fastqc.html",
        "FASTQC/trimmomatic/{sample}/{sample}_fastqc.zip"
    params:
        outdir=lambda wildcards: "FASTQC/trimmomatic/{}".format(wildcards.sample),
        threads=24,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """


## trimmomatic_unpaired_rna_to_fastqc: quality control checks on trimmed sequence data
## .fq.gz files must be located in a directory named Trimmomatic_files and end with *.fq.gz
## output will be directed to a directory named /FASTQC/trimmomatic/{sample}


rule trimmomatic_unpaired_rna_to_fastqc:
    input:
        "results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz"
    output:
        "FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.html",
        "FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.zip"
    params:
        outdir=lambda wildcards: "FASTQC/trimmomatic/unpaired/{}".format(wildcards.sample),
        threads=24,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {params.threads} {input}
        """