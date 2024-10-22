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
        # TO DO #
        ,
    params:
        outdir = lambda wildcards: "MULTIQC/{}".format(wildcards.sample),
        logdir = "logs/paired_end_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc {input} -o {params.outdir} -n {wildcards.sample}_multiqc_report.html 2> {log}
        """

rule single_end_multiqc:
    input:
        "FASTQC/single_end/{sample}/{sample}_fastqc.html"
    output:
        "MULTIQC/single_end/{sample}_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/single_end_multiqc/{sample}.log"
    envmodules:
        # TO DO #
        ,
    params:
        outdir = lambda wildcards: "MULTIQC/single_end/{}".format(wildcards.sample),
        logdir = "logs/single_end_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc {input} -o {params.outdir} -n {wildcards.sample}_multiqc_report.html 2> {log}
        """

rule sortmerna_paired_multiqc:
    input:
        "FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fwd_fastqc.html",
        "FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc.html"
    output:
        "MULTIQC/sortmerna/{sample}_paired_multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/sortmerna_paired_multiqc/{sample}.log"
    envmodules:
        # TO DO #
        ,
    params:
        outdir = lambda wildcards: "MULTIQC/sortmerna/{}".format(wildcards.sample),
        logdir = "logs/sortmerna_paired_multiqc/"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc {input} -o {params.outdir} -n {wildcards.sample}_multiqc_report.html 2> {log}
        """