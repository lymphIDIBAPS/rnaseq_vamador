containerized: "docker://Dockerfile"

include: "workflow/rules/fastqc.smk"
include: "workflow/rules/rnaseq_paired.smk"
include: "workflow/rules/rnaseq_not_paired.smk"
include: "workflow/rules/multiqc.smk"


configfile: "config/config.yaml"

# CONDITIONAL RULES #
# FASTQC if specified
run_fastqc = config["run_fastqc"] == "yes"

# MULTIQC if specified
run_multiqc = config["run_multiqc"] == "yes"


rule all:
    input:
        ## The lines below will perform fastqc in the raw RNA data
        #
        # expand("FASTQC/{sample}/{sample}_1_1_fastqc.html", sample = paired_samples) if run_fastqc else [],
        # expand("FASTQC/single_end/{sample}/{sample}_fastqc.html", sample = single_end_samples) if run_fastqc else [],
        #
        ## The lines below will perform fastqc in the filtered and trimmed RNA
        #
        # expand("FASTQC/sortmerna/filtered/{sample}_fwd/{sample}_fwd_fastqc.html", sample = paired_samples) if run_fastqc else [],
        # expand("FASTQC/sortmerna/filtered/{sample}_rev/{sample}_rev_fastqc.html", sample = paired_samples) if run_fastqc else [],
        # expand("FASTQC/trimmomatic/{sample}/{sample}_fwd_p_fastqc.html", sample = paired_samples) if run_fastqc else [],
        # expand("FASTQC/sortmerna/unpaired/filtered/{sample}/{sample}_fastqc.html", sample = single_end_samples) if run_fastqc else [],
        # expand("FASTQC/trimmomatic/unpaired/{sample}/{sample}_fwd_fastqc.html", sample = single_end_samples) if run_fastqc else [],
        #
        ## The lines below will perform multiQC in the raw RNA data
        #
        expand("MULTIQC/single_end/{sample}/{sample}_multiqc_report.html", sample = single_end_samples) if run_multiqc else [],
        expand("MULTIQC/{sample}/{sample}_multiqc_report.html", sample = paired_samples) if run_multiqc else [],
        #
        ## The lines below will perform multiQC in the filtered and trimmed RNA
        # 
        expand("MULTIQC/sortmerna/{sample}/{sample}_unpaired_multiqc_report.html", sample = single_end_samples) if run_multiqc else [],
        expand("MULTIQC/sortmerna/{sample}/{sample}_fwd_paired_multiqc_report.html", sample = paired_samples) if run_multiqc else [],
        expand("MULTIQC/trimmomatic/unpaired/{sample}/{sample}_unpaired_multiqc_report.html", sample = single_end_samples) if run_multiqc else [],
        expand("MULTIQC/trimmomatic/{sample}/{sample}_paired_multiqc_report.html", sample = paired_samples) if run_multiqc else [],
        # 
        ## The lines below run the full quantification pipeline, for unpaired and paired samples
        #
        expand("results/kallisto_files/unpaired/{sample}/abundance.h5", sample = single_end_sample_name),
        expand("results/kallisto_files/{sample}/abundance.h5", sample = paired_end_sample_name)
    shell:
        'echo "I just run subrules!"'