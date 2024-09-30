## Snakefile - RNAseq not Paired
##
configfile: "config/config.yaml"

import os

# Define the directory containing the FASTQ files
single_end_dir = config["fastq_dir_single_end"]

# Define patterns to match specific files
single_end_sample_name = glob_wildcards(f"{single_end_dir}/{{sample}}.fastq").sample


## sortmerna_not_paired: filter RNA from non paired reads
## .fastq files must be located in a directory named FASTQ/single_end
## output will be directed to a directory named results/sortmerna_files/unpaired/rRNA or /rRNAf, 
## the first are aligned reads and the second rejected reads 


rule sortmerna_not_paired:
    input:
        ref = ["resources/rRNA_databases_v4/smr_v4.3_default_db.fasta"],
        reads = "FASTQ/single_end/{sample}.fastq",
    output:
        aligned = "results/sortmerna_files/unpaired/rRNA/{sample}.fq",
        other = "results/sortmerna_files/unpaired/rRNAf/{sample}.fq",
        stats = "results/sortmerna_files/unpaired/rRNA/{sample}.log",
    params:
        extra = "--idx-dir ./idx",
    threads: 6
    resources:
        mem_mb=4096,  # amount of memory for building the index
    conda:
        "../envs/sortmerna.yaml"
    log:
        "logs/sortmerna_not_paired/{sample}.log",
    envmodules:
        "sortmerna/4.3.6"
    wrapper:
        # Use first line if we are running in an offline cluster, 
        # the second one if we are running it in an online machine
        "file:/slgpfs/projects/rcli/rcli186814/snakemake_wrappers/bio/sortmerna/"
        # "v3.14.0/bio/sortmerna"


## rna_trimming_not_paired: trim and crop filtered data
## trim and crop Illumina (FASTQ) data and remove adapters.
## input files will be located in a directory named results/sortmerna_files/unpaired/rRNAf
## output will be directed to a directory named results/trimmomatic_files/unpaired/{sample}

rule rna_trimming_not_paired:
    input:
        read = "results/sortmerna_files/unpaired/rRNAf/{sample}.fq"
    output:
        file = "results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz"
    params:
        threads = 6,
        log_dir = "logs/rna_trimming_not_paired/"
    conda:
        "../envs/rnaseq.yaml"
    log:
        "logs/rna_trimming_not_paired/{sample}.log"
    envmodules:
        "java/12.0.2"
    shell:
        """
        mkdir -p results/trimmomatic_files/unpaired {params.log_dir}
        java -jar /apps/TRIMMOMATIC/0.39/trimmomatic-0.39.jar SE -threads {params.threads} -phred33 -trimlog {log} {input.read} \
        {output.file} \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50
        """
# If we want to run in our PC, we need to remove "java -jar /apps/TRIMMOMATIC/0.39/trimmomatic-0.39.jar" and leave only "trimmomatic" 


# This will perform the following in this order:

# Remove Illumina adapters provided in the TruSeq3-SE.fa file (provided). Initially Trimmomatic will look 
# for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped 
# if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single 
# ended reads a score of 10, (about 17 bases). Scan the read with a 5-base wide sliding window, cutting when 
# the average quality per base drops below 20. Drop reads which are less than 50 bases long after these steps
# Single-end mode requires 1 input files and outputs 1 file.



## kallisto_quant: runs the quantification algorithm
## outputs three files: abundance.h5 (read by sleuth), abundance.tsv (plaintext) and run_info.json (log file)
## The index is constructed with the rule kallisto_index in the rnaseq_paired.smk file


rule kallisto_quant_not_paired:
    input:
        index_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx",
        sample_not_paired = "results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz",
    output:
        abundance = "results/kallisto_files/unpaired/{sample}/abundance.h5"
    params:
        threads = 6,
        output_dir = lambda wildcards, output: os.path.dirname(output.abundance),
        # "results/kallisto_files/unpaired/{sample}/"
        log_dir = "logs/kallisto_quant_not_paired/"
    conda:
        "../envs/rnaseq.yaml"
    log:
        "logs/kallisto_quant_not_paired/{sample}.log"
    envmodules:
        "intel/2018.3",
        "impi/2018.3",
        "zlib/1.2.11",
        "gcc/12.2.0",
        "hdf5/1.10.2",
        "szip/2.1.1",
        "kallisto/0.46.1"
    shell:
        """
        mkdir -p {params.output_dir} {params.log_dir}
        kallisto quant -i {input.index_path} -o {params.output_dir} --single -l 260 -s 20 -t {params.threads} {input.sample_not_paired} >{log}
        """