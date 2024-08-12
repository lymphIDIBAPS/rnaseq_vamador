## Snakefile - RNAseq Paired
##
# configfile: "config/config.yaml"


# Define the directory containing the FASTQ files
fastq_dir = "FASTQ/rna_paired/unmerged"

# Define patterns to match specific files
paired_end_sample_name = glob_wildcards(fastq_dir + "/{sample}_1_1.fastq").sample


# In the original scipt, we have the following code snippet, that i dont
# know why is present. Will research.
# dir=$(basename "${1/_1_sortmerna/}")


## merge_technical_replicates: merge data from L1_F and L2_F
## .fastq files must be located in a directory named FASTQ/rna_paired
## output will be directed to a directory named FASTQ/rna_paired/merged


rule merge_technical_replicates:
    input:
        L1_F = "FASTQ/rna_paired/unmerged/{sample}_1_1.fastq",
        L2_F = "FASTQ/rna_paired/unmerged/{sample}_2_1.fastq",
        L1_R = "FASTQ/rna_paired/unmerged/{sample}_1_2.fastq",
        L2_R = "FASTQ/rna_paired/unmerged/{sample}_2_2.fastq",
    output:
        merged_F = "FASTQ/rna_paired/merged/{sample}_mergedF.fastq",
        merged_R = "FASTQ/rna_paired/merged/{sample}_mergedR.fastq",
    shell:
        """
        mkdir -p FASTQ/rna_paired/merged
        cat {input.L1_F} {input.L2_F} > {output.merged_F}
        cat {input.L1_R} {input.L2_R} > {output.merged_R}
        """


## rna_filtering: filter RNA
## .fastq files must be located in a directory named FASTQ/rna_paired/merged
## output will be directed to a directory named results/sortmerna_files/rRNA or /rRNAf, 
## the first are aligned reads and the second rejected reads 


rule rna_filtering:
    input:
        merged_F = "FASTQ/rna_paired/merged/{sample}_mergedF.fastq",
        merged_R = "FASTQ/rna_paired/merged/{sample}_mergedR.fastq",
    output:
        aligned = "results/sortmerna_files/rRNA/{sample}_rev.fq",
        forward = "results/sortmerna_files/rRNAf/{sample}_fwd.fq",
        other = "results/sortmerna_files/rRNAf/{sample}_rev.fq",
    params:
        aligned = "results/sortmerna_files/rRNA/{sample}",
        other = "results/sortmerna_files/rRNAf/{sample}",
        threads = 24
    conda:
        "../envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/sortmerna_files/rRNA
        mkdir -p results/sortmerna_files/rRNAf
        sortmerna --ref resources/rRNA_databases_v4/smr_v4.3_default_db.fasta --reads {input.merged_F} --reads {input.merged_R} --aligned {params.aligned} --other {params.other} --workdir /home/oscar/rnaseq --fastx --paired_in -threads {params.threads} -out2 -v --idx-dir ./idx
        rm -r ./kvdb/
        """

# IDX Dir is the index directory, which will be built the first time we run sortmeRNA
# --idx-dir ./idx

#  The correct --ref databases must be this one if we want to copy exactly the pipeline form Marta Sureda:
#  --ref resources/rRNA_databases_v4/silva-bac-16s-id90.fasta --ref resources/rRNA_databases_v4/silva-bac-23s-id98.fasta --ref resources/rRNA_databases_v4/silva-arc-16s-id95.fasta  --ref resources/rRNA_databases_v4/silva-arc-23s-id98.fasta --ref resources/rRNA_databases_v4/silva-euk-18s-id95.fasta  --ref resources/rRNA_databases_v4/silva-euk-28s-id98.fasta --ref resources/rRNA_databases_v4/rfam-5s-database-id98.fasta  --ref resources/rRNA_databases_v4/rfam-5.8s-database-id98.fasta
#  However, I have found that the creator of the package recommends using only the following database:
#  smr_v4.3_default_db.fasta , from:(https://github.com/sortmerna/sortmerna/issues/292)



## rna_trimming: trim and crop filtered data
## trim and crop Illumina (FASTQ) data and remove adapters.
## .fq files will be located in a directory named results/sortmerna_files/rRNAf
## output will be directed to a directory named results/trimmomatic_files/{sample}, 


rule rna_trimming:
    input:
        forward = "results/sortmerna_files/rRNAf/{sample}_fwd.fq",
        rev = "results/sortmerna_files/rRNAf/{sample}_rev.fq"
    output:
        forward_paired = "results/trimmomatic_files/{sample}_fwd_p.fq.gz",
        forward_unpaired = "results/trimmomatic_files/{sample}_fwd_up.fq.gz",
        rev_paired = "results/trimmomatic_files/{sample}_rev_p.fq.gz",
        rev_unpaired = "results/trimmomatic_files/{sample}_rev_up.fq.gz"
    params:
        threads = 24
    conda:
        "../envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/trimmomatic_files
        trimmomatic PE -threads {params.threads} -phred33 {input.forward} {input.rev} \
        {output.forward_paired} {output.forward_unpaired} {output.rev_paired} {output.rev_unpaired} \
        ILLUMINACLIP:TruSeq3-PE:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50
        """

# This configuration of the trimmomatic run is taken from the script of Marta Sureda.
# This will perform the following in this order:

# Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look 
# for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped 
# if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single 
# ended reads a score of 10, (about 17 bases). Scan the read with a 5-base wide sliding window, cutting when 
# the average quality per base drops below 20. Drop reads which are less than 50 bases long after these steps

# Paired-end mode requires 2 input files (for forward and reverse reads) and 4 output files (for
# forward paired, forward unpaired, reverse paired and reverse unpaired reads).

# 2 input files: FASTQ/sortmerna_files/$rRNAf_F FASTQ/sortmerna_files/$rRNAf_R

# 4 output files: 

# FASTQ/Trimmomatic_files/$trimmed_F
# FASTQ/Trimmomatic_files/$trimmed_F2
# FASTQ/Trimmomatic_files/$trimmed_R
# FASTQ/Trimmomatic_files/$trimmed_R2

# (for forward paired, forward unpaired, reverse paired and reverse unpaired reads).



## kallisto_index: builds an index
## from a FASTA formatted file of target sequences. Compute intensive rule

rule kallisto_index:
    input:
        index_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        index_out_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"
    params:
        threads = 24
    conda:
        "../envs/rnaseq.yaml"
    shell:
        """
        kallisto index -i {output.index_out_path} --threads={params.threads} {input.index_path} 
        """

# kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-100.idx Homo_sapiens.GRCh38.cdna.all.fa.gz


## kallisto_quant: runs the quantification algorithm
## outputs three files: abundance.h5 (read by sleuth), abundance.tsv (plaintext) and run_info.json (log file)

rule kallisto_quant:
    input:
        index_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx",
        rev_paired = "results/trimmomatic_files/{sample}_rev_p.fq.gz",
        for_paired = "results/trimmomatic_files/{sample}_fwd_p.fq.gz",
    output:
        abundance = "results/kallisto_files/{sample}/abundance.h5",
        abundance_tsv = "results/kallisto_files/{sample}/abundance.tsv",
    conda:
        "../envs/rnaseq.yaml"
    params:
        threads = 24,
        output_dir = "results/kallisto_files/{sample}/"
    shell:
        """
        mkdir -p {params.output_dir}
        kallisto quant -i {input.index_path} -o {params.output_dir} -t {params.threads} {input.for_paired} {input.rev_paired}
        """