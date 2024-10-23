## Snakefile - MULTIQC
##
configfile: "config/config.yaml"

rule paired_end_multiqc:
    input:
        "FASTQC/{sample}/{sample}_1_1_fastqc.html",
        "FASTQC/{sample}/{sample}_1_2_fastqc.html",
        "FASTQC/{sample}/{sample}_2_1_fastqc.html",
        "FASTQC/{sample}/{sample}_2_2_fastqc.html"
    output:
        "MULTIQC/{sample}/{sample}_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/paired_end_multiqc/{sample}.log"
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    params:
        outdir = lambda wildcards: "MULTIQC/{}".format(wildcards.sample),
        logdir = "logs/paired_end_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc FASTQC/{wildcards.sample}/ -o {params.outdir} -n {wildcards.sample}_multiqc_report.html 2> {log}
        """

rule single_end_multiqc:
    input:
        "FASTQC/single_end/{sample}/{sample}_fastqc.html"
    output:
        "MULTIQC/single_end/{sample}/{sample}_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/single_end_multiqc/{sample}.log"
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    params:
        outdir = lambda wildcards: "MULTIQC/single_end/{}".format(wildcards.sample),
        logdir = "logs/single_end_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc FASTQC/single_end/{wildcards.sample}/ -o {params.outdir} -n {wildcards.sample}_multiqc_report.html 2> {log}
        """

rule sortmerna_paired_multiqc:
    input:
        "FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fwd_fastqc.html",
        "FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc.html"
    output:
        "MULTIQC/sortmerna/{sample}_fwd_paired_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/sortmerna_paired_multiqc/{sample}.log"
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    params:
        outdir = lambda wildcards: "MULTIQC/sortmerna/{}".format(wildcards.sample),
        logdir = "logs/sortmerna_paired_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc FASTQC/sortmerna/filtered/{wildcards.sample}_fwd/ -o {params.outdir} -n {wildcards.sample}_fwd_paired_multiqc_report.html 2> {log}
        multiqc FASTQC/sortmerna/filtered/{wildcards.sample}_rev/ -o {params.outdir} -n {wildcards.sample}_rev_paired_multiqc_report.html 2> {log}
        """

rule sortmerna_unpaired_multiqc:
    input:
        "FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.html"
    output:
        "MULTIQC/sortmerna/{sample}/{sample}_unpaired_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/sortmerna_unpaired_multiqc/{sample}.log"
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    params:
        outdir = lambda wildcards: "MULTIQC/sortmerna/{}".format(wildcards.sample),
        logdir = "logs/sortmerna_unpaired_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc FASTQC/sortmerna/unpaired/filtered/{wildcards.sample}/ -o {params.outdir} -n {wildcards.sample}_unpaired_multiqc_report.html 2> {log}
        """

rule trimmomatic_paired_multiqc:
    input:
        "FASTQC/trimmomatic/{sample}/{sample}_fwd_p_fastqc.html",
        "FASTQC/trimmomatic/{sample}/{sample}_fwd_up_fastqc.html",
        "FASTQC/trimmomatic/{sample}/{sample}_rev_p_fastqc.html",
        "FASTQC/trimmomatic/{sample}/{sample}_rev_up_fastqc.html"
    output:
        "MULTIQC/trimmomatic/{sample}/{sample}_paired_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/trimmomatic_paired_multiqc/{sample}.log"
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    params:
        outdir = lambda wildcards: "MULTIQC/trimmomatic/{}".format(wildcards.sample),
        logdir = "logs/trimmomatic_paired_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc FASTQC/trimmomatic/{wildcards.sample}/ -o {params.outdir} -n {wildcards.sample}_paired_multiqc_report.html 2> {log}
        """

rule trimmomatic_unpaired_multiqc:
    input:
        "FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.html",
    output:
        "MULTIQC/trimmomatic/unpaired/{sample}/{sample}_unpaired_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/trimmomatic_unpaired_multiqc/{sample}.log"
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    params:
        outdir = lambda wildcards: "MULTIQC/trimmomatic/unpaired/{}".format(wildcards.sample),
        logdir = "logs/trimmomatic_unpaired_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc FASTQC/trimmomatic/unpaired/{wildcards.sample}/ -o {params.outdir} -n {wildcards.sample}_unpaired_multiqc_report.html 2> {log}
        """