# fDOG - Feature-aware Directed OrtholoG search
[![published in: MBE](https://img.shields.io/badge/published%20in-MBE-ff69b4)](https://doi.org/10.1093/molbev/msaf120)
[![PyPI version](https://badge.fury.io/py/fdog.svg)](https://pypi.org/project/fdog/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Github Build](https://github.com/BIONF/fDOG/workflows/build/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17250793.svg)](https://doi.org/10.5281/zenodo.17250793)


# Table of Contents
* [How to install](#how-to-install)
     * [Install the fDOG package](#install-the-fdog-package)
     * [Setup fDOG](#setup-fdog)
* [Usage](#usage)
* [fDOG data set](#fdog-data-set)
     * [Adding a new gene set into fDOG](#adding-a-new-gene-set-into-fdog)
     * [Adding a list of gene sets into fDOG](#adding-a-list-of-gene-sets-into-fdog)
* [fDOG-Assembly](#fdog-assembly)
* [Bugs](#bugs)
* [How to cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*fDOG* tool is distributed as a python package called *fdog*. It is compatible with [Python ≥ v3.12](https://www.python.org/downloads/).

## Install the fDOG package

**_RECOMMENDATION:_** Install fDOG in a fresh conda environment to ensure compatibility and avoid conflicts with other packages.

### Using a Conda Environment

1. Follow [this instruction](https://mamba.readthedocs.io/en/latest/) to install Mamba or Micromamba (faster package manager for conda)

2. Create a new environment

```
mamba create -n fdog python -y
```

3. Activate the environment

```
mamba activate fdog
```

4. Install *fdog* using `pip`:

```
python3 -m pip install fdog
```

### Without conda

1. Install *fdog* globally (requires admin rights)

```
python3 -m pip install fdog
```

2. Install *fdog* for a single user (no admin rights needed)

```
python3 -m pip install --user fdog
```

3. Add local bin to PATH (if using `--user`)

Append this line to the end of your **~/.bashrc** or **~/.bash_profile** file

```
export PATH=$HOME/.local/bin:$PATH
```

Then, reload the current terminal to apply the change (or run `source ~/.bashrc`)

## Setup fDOG

After installing *fdog*, you need to setup *fdog* to get its dependencies and pre-calculated data.

**IMPORTANT NOTE**: if you haven't installed [greedyFAS](https://github.com/BIONF/FAS), it will be automatically installed during the *fDOG* setup. After installation, you must run [setupFAS](https://github.com/BIONF/FAS/wiki/setupFAS) before using *fDOG* with *FAS*! This step is required to configure *FAS* correctly. You can run *fDOG* without *FAS* by adding the `--fasOff` option to your command. However, it is recommended to use *FAS* to access all the features of *fDOG*.

You can setup fDOG by running this command
```
fdog.setup -d /output/path/for/fdog/data
```

[Pre-calculated data set](https://github.com/BIONF/fDOG/wiki/Input-and-Output-Files#data-structure) of fdog will be saved in `/output/path/for/fdog/data`. After the setup run successfully, you can start using *fdog*. **Please make sure to check if you need to run [setupFAS](https://github.com/BIONF/FAS/wiki/setupFAS) first.**

You will get a warning if any of the dependencies are not ready to use, please solve those issues and rerun `fdog.setup`.

*For debugging the setup, please create a log file by running the setup as e.g. `fdog.setup | tee log.txt` and send us that log file, so that we can trouble shoot the issues. Most of the problems can be solved by just re-running the setup.*

# Usage

Once *fdog* is installed and set up correctly, it can be run using the provided sample input file 'infile.fa'

## Running fDOG with FAS
```
fdog.run --seqFile infile.fa --jobName test --refspec HUMAN@9606@qfo24_02
```

## Running fDOG without FAS

If FAS has not been set up, add the `--fasOff` option

```
fdog.run --seqFile infile.fa --jobName test --refspec HUMAN@9606@qfo24_02 --fasOff
```

## Output

All output files will be saved in your **current working directory** with the prefix specified in `--jobName` (e.g. `test`).

## Viewing all options

You can have an overview about all available options with the command
```
fdog.run -h
```

Please find more information in [our wiki](https://github.com/BIONF/fDOG/wiki) to learn about the [input and outputs files](https://github.com/BIONF/fDOG/wiki/Input-and-Output-Files) of *fdog*.

# fDOG data set

Within the data package we provide a set of [81 reference taxa](https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/QfO_release_2024_02.tar.gz). They will be automatically downloaded during the setup. This data comes "ready to use" with the *fdog* framework. Species data must be present in the three directories listed below:

* searchTaxa_dir (Contains sub-directories for proteome fasta files for each species)
* coreTaxa_dir (Contains sub-directories for BLAST databases made with `makeblastdb` out of your proteomes)
* annotation_dir (Contains feature annotation files for each proteome)

For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version])

*fdog* is not limited to those 81 reference taxa. If needed the user can manually add further gene sets (multiple fasta format) using provided functions.

## Adding a new gene set into fDOG
For adding **one gene set**, please use the `fdog.addTaxon` function:
```
fdog.addTaxon -f newTaxon.fa -i tax_id [-o /output/directory] [-n abbr_tax_name] [-c] [-v protein_version] [-a]
```

in which, the first 3 arguments are required including `newTaxon.fa` is the gene set that need to be added, `tax_id` is its NCBI taxonomy ID, `/output/directory` is where the sub-directories can be found (*genome_dir*, *blast_dir* and *weight_dir*). If not given, new taxon will be added into the same directory of pre-calculated data. Other arguments are optional, which are `-n` for specify your own taxon name (if not given, an abbriviate name will be suggested based on the NCBI taxon name of the input `tax_id`), `-c` for calculating the BLAST DB (only needed if you need to include your new taxon into the list of taxa for compilating the core set), `-v` for identifying the genome/proteome version (default will be the current date <YYMMDD>), and `-a` for turning off the annotation step (*not recommended*).

## Adding a list of gene sets into fDOG
For adding **more than one gene set**, please use the `fdog.addTaxa` script:
```
fdog.addTaxa -i /path/to/newtaxa/fasta -m mapping_file [-o /output/directory] [-c]
```
in which, `/path/to/taxa/fasta` is a folder where the FASTA files of all new taxa can be found. `mapping_file` is a tab-delimited text file, where you provide the taxonomy IDs that stick with the FASTA files:

```
#filename	tax_id	abbr_tax_name	version
filename1.fa	12345678
filename2.faa	9606
filename3.fasta	4932	my_fungi
...
```

The header line (started with #) is a Must. The values of the last 2 columns (abbr. taxon name and genome version) are, however, optional. If you want to specify a new version for a genome, you need to define also the abbr. taxon name, so that the genome version is always at the 4th column in the mapping file.

_**NOTE:** After adding new taxa into *fdog*, you should [check for the validity of the new data](https://github.com/BIONF/fDOG/wiki/Check-data-validity) before running fdog._

# fDOG-Assembly

*fDOG-Assembly* is an extension of *fDog* that enables searching for orthologs directly within unannotated genome assemblies. For more details about *fDOG-Assembly*, please refer to our [wiki page](https://github.com/BIONF/fDOG/wiki/fDOG-Assembly).

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/fDOG/issues/new) or be in touch via email.

# How to cite
Tran V, Langschied F, Muelbaier H, Dosch J, Arthen F, Balint M, Ebersberger I. 2025. Feature architecture-aware ortholog search with fDOG reveals the distribution of plant cell wall-degrading enzymes across life. Molecular Biology and Evolution:msaf120. https://doi.org/10.1093/molbev/msaf120

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Vinh Tran](https://github.com/trvinh)
- [Hannah Muelbaier](https://github.com/HannahBioI)
- [Holger Bergmann](https://github.com/holgerbgm)

# Contact
For further support or bug reports please contact: ebersberger@bio.uni-frankfurt.de
