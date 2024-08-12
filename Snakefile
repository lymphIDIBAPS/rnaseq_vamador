include: "workflow/rules/fastqc.smk"
include: "workflow/rules/rnaseq_paired.smk"
include: "workflow/rules/rnaseq_not_paired.smk"


# configfile: "config/config.yaml"


rule all:
    input:
        # expand("results/sortmerna_files/unpaired/rRNA/{sample}.log", sample = single_end_sample_name),
        # expand("results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz", sample = single_end_sample_name),
        # 
        # Below I have the paired results
        # 
        # expand("resources/rna_paired/merged/{sample}_mergedF.fastq", sample=paired_end_sample_name),
        # expand("results/sortmerna_files/rRNA/{sample}_rev.fq", sample=paired_end_sample_name),
        # expand("results/trimmomatic_files/{sample}_fwd_p.fq.gz", sample=paired_end_sample_name),
        # ("resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"),
        # 
        # 
        # The lines below will perform fastqc in the raw RNA data
        #
        expand("FASTQC/{sample}/{sample}_1_1_fastqc.html", sample = paired_samples),
        expand("FASTQC/single_end/{sample}/{sample}_fastqc.html", sample = single_end_samples),
        #
        # The lines below will perform fastqc in the filtered and trimmed RNA
        #
        expand("FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fwd_fastqc.html", sample = paired_samples),
        expand("FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc.html", sample = paired_samples),
        expand("FASTQC/trimmomatic/{sample}/{sample}_fwd_p_fastqc.html", sample = paired_samples),
        expand("FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.html", sample = single_end_samples),
        expand("FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.html", sample = single_end_samples),
        # 
        # The lines below run the full quantification pipeline, for unpaired and paired samples
        #
        expand("results/kallisto_files/unpaired/{sample}/abundance.h5", sample = single_end_sample_name),
        expand("results/kallisto_files/{sample}/abundance.h5", sample = paired_end_sample_name)
    shell:
        'echo "I just run subrules!"'