# RNA Seq Analysis Pipeline for V. Amador

This is a pipeline written in Python and bash, and with Snakemake as a workflow manager, that will output *kallisto* files (abundance.h5) from RNA seq samples.

The paired samples must be placed in the directory named ***FASTQ/rna_paired***, and depending if they are merged or not, in the ***/merged*** or ***/unmerged*** directories.

## Installation

Use the package manager [miniconda](https://docs.anaconda.com/miniconda/) to install miniconda3.

### Windows
Run these three commands with ***Windows Power Shell*** to quickly and quietly download the latest 64-bit Windows installer, rename it to a shorter file name, silently install, and then delete the installer.
```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe
Start-Process -FilePath ".\miniconda.exe" -ArgumentList "/S" -Wait
del miniconda.exe
```
After installing, from the Start menu, open the “Anaconda Powershell Prompt (miniconda3)”.

You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```

Now i will download Mamba to create snakemake
I have trouble with SSL cert
You can usually get a copy by clicking on the padlock icon in your browser when visiting any https site, then click around to view certificate, and download in PEM format.
Then I will point conda to it in our system. 
```bash
conda config --set ssl_verify <pathToYourFile>.pem
```

Next we will install mamba in our base environment with:
```bash
conda install -n base -c conda-forge mamba
```
To check:
```bash
mamba --version
# mamba 1.X.X
# conda 24.X.X
```
Create a new environment with snakemake:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```
Activate the snakemake environment:
```bash
conda activate snakemake
```

## Usage

```bash
# For a test run of the pipeline
snakemake --use-conda -j 24 -np

# For a real test of the pipeline
snakemake --use-conda -j 24
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)