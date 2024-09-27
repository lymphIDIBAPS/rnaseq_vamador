# RNA Seq Analysis Pipeline for V. Amador
![GitHub last commit](https://img.shields.io/github/last-commit/Programa-de-neoplasias-linfoides/rnaseq_virginia)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/Programa-de-neoplasias-linfoides/rnaseq_virginia)
----
This is a pipeline written in Python and bash, and with Snakemake as a workflow manager, that will output *kallisto* files (abundance.h5) from RNA seq samples. These samples can be in .fastq or compressed format. 

The samples can be placed in any directory, but the path must be specified in the ***config/config.yaml*** file. 

# Table of Contents
1. [Installation](#installation)
2. [Linux](#linux)
3. [Windows](#windows)
4. [macOS](#macos)
5. [Snakemake Environment](#snakemake-environment)
6. [Clone the repository](#clone-the-repository)
7. [Snakemake Usage](#snakemake-usage)

## Installation

<!-- Use the package manager [miniconda](https://docs.anaconda.com/miniconda/) to install miniconda3. -->
### Linux
#### Install Miniconda 3
Open a Linux shell, then run these three commands to quickly and quietly download the latest 64-bit Linux miniconda 3 installer, rename it to a shorter file name, silently install, and then delete the installer.
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
You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```
---
### Windows
<details><summary>Windows Installation</summary>

Since Windows does not have access to the majority of packages we need in the pipeline, we need to install *Linux on Windows*, also known as *WSL*. In the ***Windows Power Shell***:
```powershell
wsl --install
# This command will install the Ubuntu distribution of Linux.
```
If you run into an issue during the installation process, please check [the installation section of the troubleshooting guide](https://learn.microsoft.com/en-us/windows/wsl/troubleshooting#installation-issues).

#### Set up your Linux user info

Once you have installed WSL, you will need to create a user account and password for your newly installed Linux distribution. See the Best practices for [setting up a WSL development environment guide to learn more](https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password).

#### Install Miniconda 3
Once you have a working shell in your *WSL*, run these three commands to quickly and quietly download the latest 64-bit Linux miniconda 3 installer, rename it to a shorter file name, silently install, and then delete the installer.
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
You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```
</details>

----
### macOS
<details><summary>macOS Installation</summary>
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
</details>

---
### Snakemake Environment
Now, with miniconda installed in our machine, we can create a new environment with snakemake installed:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

In case you are downloading from the Clinic network, you will have trouble with the SLL certificate. To solve any problems, do the following:
1. You can usually get a copy by clicking on the padlock icon in your browser when visiting any https site, then click around to view certificate, and download in PEM format.

2. Then we will point conda to it in our system. 
```bash
conda config --set ssl_verify <pathToYourFile>.pem
```

Once installed, we must activate and move into the snakemake environment with:
```bash
conda activate snakemake
snakemake --version
# 8.16.X
```
If at any time we want to exit the environment, we can with```conda deactivate```, and to get back in with ```conda activate snakemake```.
To see the packages we have currently installed in the environment, we can with```conda list```.
### Clone the repository

1. Above the list of files, click Code.

![image_code](image.png)

2. Copy the URL for the repository. To clone the repository using HTTPS, under "HTTPS", copy the link provided.
3. Open a Terminal.
4. Change the current working directory to the location where you want the cloned directory. For example, ```cd rna_seq_vamador```. Make sure that the directory exists before you move into it.
5. Type ```git clone https://github.com/Programa-de-neoplasias-linfoides/rnaseq_virginia.git```.
6. Press Enter to create your local clone.
```bash
git clone https://github.com/Programa-de-neoplasias-linfoides/rnaseq_virginia.git
> Cloning into `rna_seq_vamador`...
> remote: Counting objects: 10, done.
> remote: Compressing objects: 100% (8/8), done.
> remove: Total 10 (delta 1), reused 10 (delta 1)
> Unpacking objects: 100% (10/10), done.
```
## Snakemake Usage
When we have the cloned repository, we can proced and add our sample data to the FASTQ directory. This is not mandatory, as in ***config/config.yaml*** file we can edit and set any path to our sample data. 

In the same file we can edit the number of threads our computer has, so it will run adapted to the current resources we have available. 

```bash
# For a test run of the pipeline
snakemake --use-conda -np

# For a real run of the pipeline
snakemake --use-conda
```

## Run the pipeline in a HPC
If we have many samples and our computer does not have enough computational power, we can run the pipeline in a cluster. This pipeline has been prepared to run in the [StarLife](https://www.bsc.es/supportkc/docs/StarLife/intro) cluster, in the [BSC](https://www.bsc.es/).

1. Make a new directory named ***/slgpfs/*** in your computer and mount it to the same directory in StarLife:
```bash
mkdir /home/user/slgpfs
sshfs -o allow_other your_bsc_user@sl1.bsc.es:/slgpfs/ /home/user/slgpfs/
```
This will allow you to see and work on the custer from your computer system directly.

2. On your computer, navigate to the directory: ***/home/user/slgpfs/projects/group_folder***

3. Download and extract the following file to the directory, in which we have a full conda environment ready to run snakemake:
[Snakemake Conda Environment](https://drive.google.com/file/d/1ZQB02jZhhr5G_DlbhDKFHiqZ7AVuCxqx/view?usp=sharing)

4. Clone this repository in the directory, following the steps from [Clone the repository](#clone-the-repository)

5. Now, connect to the cluster:
```bash
ssh your_username@sl1.bsc.es # or
ssh your_username@sl2.bsc.es
```
6. In the cluster, navigate to the cloned repository: ***/slgpfs/projects/group_folder/rna_seq_vamador***

7. Now, activate the snakemake_bsc environment:
```bash
source ../snakemake_bsc/bin/activate
```
In your terminal, you should now see something like: ```(snakemake_bsc) your_username@sllogin1```

8. At this point, in your local machine, you can move your samples to the directory ***/slgpfs/projects/group_folder/rna_seq_vamador/FASTQ***. 

9. Now, you can run the pipeline from the cluster with the command:
```bash
# For a test run of the pipeline
snakemake --profile config/slurm/ --use-envmodules -np

# For a real run of the pipeline
snakemake --profile config/slurm/ --use-envmodules
``` 
Remember to check the files in ***/config/slurm/config.yaml*** for the cluster configuration and the 
## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)