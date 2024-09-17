FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="59717d7f889a90d879bf63b2efcbf4d426de31d536af2c51732242af61cf2ecc"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.14.0/bio/sortmerna/environment.yaml
#   prefix: /conda-envs/fd6a53f2993f45fb1ee437ad48c7bd56
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - sortmerna =4.3.6
RUN mkdir -p /conda-envs/fd6a53f2993f45fb1ee437ad48c7bd56
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.14.0/bio/sortmerna/environment.yaml /conda-envs/fd6a53f2993f45fb1ee437ad48c7bd56/environment.yaml

# Conda environment:
#   source: workflow/envs/rnaseq.yaml
#   prefix: /conda-envs/c35889471b04dbd77e47c3c6ecebec68
#   name: rnaseq_snakemake
#   channels:
#     - conda-forge
#     - defaults
#     - bioconda
#   dependencies:
#     - python
#     # should be sortmerna=4.3.2, but conda cant solve env with 4.3.2
#     - sortmerna=4.3.2
#     - trimmomatic=0.39
#     - kallisto=0.46.1
#     - gsea
#   
#     ## R (Bioconductor)
#     - bioconductor-tximport=1.12.3
#     - bioconductor-deseq2=1.24.0
#   
#     ## R (CRAN)
#     - r-base
RUN mkdir -p /conda-envs/c35889471b04dbd77e47c3c6ecebec68
COPY workflow/envs/rnaseq.yaml /conda-envs/c35889471b04dbd77e47c3c6ecebec68/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/fd6a53f2993f45fb1ee437ad48c7bd56 --file /conda-envs/fd6a53f2993f45fb1ee437ad48c7bd56/environment.yaml && \
    mamba env create --prefix /conda-envs/c35889471b04dbd77e47c3c6ecebec68 --file /conda-envs/c35889471b04dbd77e47c3c6ecebec68/environment.yaml && \
    mamba clean --all -y
