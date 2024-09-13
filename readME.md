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

Check for a good installation with 
```bash
conda --version
# conda 24.7.*

conda list
# outputs a list of packages installed in the current environment (base)
```

## Usage

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)