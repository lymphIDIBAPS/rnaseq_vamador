## Snakefile - RNAseq Single End
##

# Define the directory containing the FASTQ files
single_end_dir = "FASTQ/single_end"

# Define patterns to match specific files
sample_name = glob_wildcards(single_end_dir + "/{sample}.fastq").sample


## rna_filtering: filter RNA from non paired reads
## .fastq files must be located in a directory named FASTQ/single_end
## output will be directed to a directory named results/sortmerna_files/unpaired/rRNA or /rRNAf, 
## the first are aligned reads and the second rejected reads 

# rule make_directory:
#     input:
#         reads = "FASTQ/single_end/"
#     output:
#         dir = directory("results/sortmerna_files/unpaired/rRNA")
#     shell:
#         """
#         mkdir results/sortmerna_files/unpaired/rRNA
#         mkdir results/sortmerna_files/unpaired/rRNAf
#         """


rule rna_filtering_not_paired:
    input:
        reads = "FASTQ/single_end/{sample}.fastq",
    output:
        aligned = "results/sortmerna_files/unpaired/rRNA/{sample}.log"
    params:
        aligned = "results/sortmerna_files/unpaired/rRNA/{sample}",
        other = "results/sortmerna_files/unpaired/rRNAf/{sample}",
        threads = 24
    shadow: "minimal"
    conda:
        "../envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/sortmerna_files/unpaired/rRNA results/sortmerna_files/unpaired/rRNAf
        sortmerna --ref /home/oscar/rnaseq/resources/rRNA_databases_v4/smr_v4.3_default_db.fasta --reads {input.reads} --aligned {params.aligned} --other {params.other} --workdir /home/oscar/rnaseq --fastx -threads {params.threads} -v --idx-dir ./idx
        """

# --kvdb {params.kvdb}

# mkdir -p results/sortmerna_files/unpaired/rRNA results/sortmerna_files/unpaired/rRNAf
# sortmerna --ref /home/oscar/rnaseq/resources/rRNA_databases_v4/smr_v4.3_default_db.fasta --reads FASTQ/single_end/Zcr2_02161AAC_CAGATC.fastq --aligned results/sortmerna_files/unpaired/rRNA/Zcr2_02161AAC_CAGATC --other results/sortmerna_files/unpaired/rRNAf/Zcr2_02161AAC_CAGATC --workdir /home/oscar/rnaseq --fastx -threads 24 -v --idx-dir ./idx
# rm -r ./kvdb/

# mkdir -p results/sortmerna_files/unpaired/rRNA results/sortmerna_files/unpaired/rRNAf
# sortmerna --ref /home/oscar/rnaseq/resources/rRNA_databases_v4/smr_v4.3_default_db.fasta --reads FASTQ/single_end/Zcr3_02163AAC_AGTTCC.fastq --aligned results/sortmerna_files/unpaired/rRNA/Zcr3_02163AAC_AGTTCC --other results/sortmerna_files/unpaired/rRNAf/Zcr3_02163AAC_AGTTCC --workdir /home/oscar/rnaseq --fastx -threads 24 -v --idx-dir ./idx
# rm -r ./kvdb/

#Sample Z138WT1

# sortmerna --ref ./sortmernaDB/silva-bac-16s-id90.fasta,
# ./sortmernaDB/index/silva-bac-16s-db:./sortmernaDB/silva-bac-23s-id98.fasta,
# ./sortmernaDB/index/silva-bac-23s-db:./sortmernaDB/silva-arc-16s-id95.fasta,
# ./sortmernaDB/index/silva-arc-16s-db:./sortmernaDB/silva-arc-23s-id98.fasta,
# ./sortmernaDB/index/silva-arc-23s-db:./sortmernaDB/silva-euk-18s-id95.fasta,
# ./sortmernaDB/index/silva-euk-18s-db:./sortmernaDB/silva-euk-28s-id98.fasta,
# ./sortmernaDB/index/silva-euk-28s:./sortmernaDB/rfam-5s-database-id98.fasta,
# ./sortmernaDB/index/rfam-5s-db:./sortmernaDB/rfam-5.8s-database-id98.fasta,
# ./sortmernaDB/index/rfam-5.8s-db 
# --reads ./FASTQ/Zwt1_02158AAC_ATCACG.fastq 
# --aligned ./sortmerna_files/ZWT1_rRNA 
# --other ./sortmerna_files/ZWT1_rRNAf 
# --fastx --log -a 16 -v


## Trimming reads

rule rna_trimming_not_paired:
    input:
        read = "results/sortmerna_files/unpaired/rRNAf/{sample}.fq"
    output:
        file = "results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz"
    params:
        threads = 24
    conda:
        "../envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/trimmomatic_files/unpaired
        trimmomatic SE -threads {params.threads} -phred33 {input.read} \
        {output.file} \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50
        """

# This configuration of the trimmomatic run is taken from the script of Marta Sureda.
# This will perform the following in this order:

# Remove Illumina adapters provided in the TruSeq3-SE.fa file (provided). Initially Trimmomatic will look 
# for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped 
# if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single 
# ended reads a score of 10, (about 17 bases). Scan the read with a 5-base wide sliding window, cutting when 
# the average quality per base drops below 20. Drop reads which are less than 50 bases long after these steps

# Single-end mode requires 1 input files and outputs 1 file.


rule kallisto_quant_not_paired:
    input:
        index_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx",
        sample_not_paired = expand("results/trimmomatic_files/unpaired/{sample}_fwd.fq.gz", sample = sample_name)
    output:
        abundance = "results/kallisto_files/unpaired/abundance.h5"
    params:
        threads = 24,
        output_dir = "results/kallisto_files/unpaired"
    conda:
        "../envs/rnaseq.yaml"
    shell:
        """
        mkdir -p {params.output_dir}
        kallisto quant -i {input.index_path} -o {params.output_dir} --single -l 260 -s 20 -t {params.threads} {input.sample_not_paired}
        """

# kallisto quant -i ./Homo_sapiens.GRCh38.cdna.all.release-100.idx 
# -o ./Kallisto_files/ZWT1_aligned 
# --single -l 260 -s 20 --threads=16 
# ./Trimmomatic_files/ZWT1_trimmed.fq.gz