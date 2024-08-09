include: "workflow/rules/fastqc.smk"
include: "workflow/rules/rnaseq_paired.smk"
include: "workflow/rules/rnaseq_not_paired.smk"


# configfile: "config/config.yaml"


rule all:
    input:
        # expand("results/sortmerna_files/unpaired/rRNA/{sample}.log", sample = sample_name),
        # expand("results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz", sample = sample_name),
        # 
        # Below I have the paired results
        # 
        # expand("resources/rna_paired/merged/{sample}_mergedF.fastq", sample=L1_F),
        # expand("results/sortmerna_files/rRNA/{sample}_rev.fq", sample=L1_F),
        # expand("results/trimmomatic_files/{sample}_fwd_p.fq.gz", sample=L1_F),
        # ("resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"),
        # 
        # 
        # The lines below will perform fastqc in the raw RNA data
        #
        # expand("FASTQC/{sample}/{sample}_fastqc.html", sample=fastq_samples),
        # expand("FASTQC/single_end/{sample}_fastqc.html", sample=single_end_fastq_samples),
        #
        # The lines below will perform fastqc in the filtered and trimmed RNA
        #
        expand("FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fastqc.html", sample = fastq_samples),
        expand("FASTQC/sortmerna/filtered/{sample}_rev/{sample}_fastqc.html", sample = fastq_samples),
        expand("FASTQC/trimmomatic/{sample}/{sample}_fwd_p_fastqc.html", sample = fastq_samples),
        expand("FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.html", sample = fastq_samples),
        expand("FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.html", sample = fastq_samples),

        # The lines below run the full quantification pipeline, for unpaired and paired samples
        #
        # expand("results/kallisto_files/unpaired/{sample}/abundance.h5", sample = sample_name),
        # expand("results/kallisto_files/{sample}/abundance.h5", sample = L1_F),
    shell:
        'echo "I just run subrules!"'