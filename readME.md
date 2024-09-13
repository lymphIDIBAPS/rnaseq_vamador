# RNA Seq Analysis Pipeline for V. Amador
![GitHub last commit](https://img.shields.io/github/last-commit/Programa-de-neoplasias-linfoides/rnaseq_virginia)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/Programa-de-neoplasias-linfoides/rnaseq_virginia)
----
This is a pipeline written in Python and bash, and with Snakemake as a workflow manager, that will output *kallisto* files (abundance.h5) from RNA seq samples.

The paired samples must be placed in the directory named ***FASTQ/rna_paired***, and depending if they are merged or not, in the ***/merged*** or ***/unmerged*** directories.

# Table of Contents
1. [Installation](#installation)
2. [Windows](#windows)
3. [macOS](#macos)
4. [Snakemake Environment](#snakemake-environment)
5. [Snakemake Usage](#snakemake-usage)

## Example2
## Third Example
## [Fourth Example](http://www.fourthexample.com) 

## Installation

<!-- Use the package manager [miniconda](https://docs.anaconda.com/miniconda/) to install miniconda3. -->

### Windows
Since Windows does not have access to the majority of packages we nedd in the pipeline, we need to install *Linux on Windows*, also known as *WSL*. In the ***Windows Power Shell***:
```powershell
wsl --install
# This command will install the Ubuntu distribution of Linux.
```
If you run into an issue during the installation process, please check [the installation section of the troubleshooting guide.](https://learn.microsoft.com/en-us/windows/wsl/troubleshooting#installation-issues).

#### Set up your Linux user info

Once you have installed WSL, you will need to create a user account and password for your newly installed Linux distribution. See the Best practices for [setting up a WSL development environment guide to learn more](https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password).

#### Install Miniconda 3
Once you have a working shell, run these three commands to quickly and quietly download the latest 64-bit Linux miniconda 3 installer, rename it to a shorter file name, silently install, and then delete the installer.
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:
```bash
~/miniconda3/bin/conda init bash
```
```zsh
~/miniconda3/bin/conda init zsh
```
### macOS
These four commands download the latest M1 version of the MacOS installer, rename it to a shorter file name, silently install, and then delete the installer:
```bash
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:
```bash
~/miniconda3/bin/conda init bash
```
```zsh
~/miniconda3/bin/conda init zsh
```

You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```
<!-- #### Install mamba
Now we will install ***mamba***, that is a fast-ligthweigth package manager aking to conda.

In case you are downloading from the Clinic network, you will have trouble with the SLL certificate. TO solve any problems, do the following:
1. You can usually get a copy by clicking on the padlock icon in your browser when visiting any https site, then click around to view certificate, and download in PEM format.
Then we will point conda to it in our system. 
```bash
conda config --set ssl_verify <pathToYourFile>.pem
```

2. Next we will install mamba in our base environment with:
```bash
conda install -n base -c conda-forge mamba
```
To check:
```bash
mamba --version
# mamba 1.X.X
# conda 24.X.X
``` -->

---
#### Snakemake-Environment
Create a new environment with snakemake installed:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```
Activate the snakemake environment with:
```bash
conda activate snakemake
```

## Snakemake-Usage

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