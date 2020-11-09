# fDOG - Feature-aware Directed OrtholoG search
[![PyPI version](https://badge.fury.io/py/fdog.svg)](https://pypi.org/project/fdog/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/BIONF/fDOG.svg?branch=master)](https://travis-ci.com/BIONF/fDOG)

# Table of Contents
* [How to install](#how-to-install)
     * [Install the fDOG package](#install-the-fdog-package)
     * [Setup fDOG](#setup-fdog)
* [Usage](#usage)
* [fDOG data set](#fdog-data-set)
     * [Adding a new gene set into fDOG](#adding-a-new-gene-set-into-fdog)
     * [Adding a list of gene sets into fDOG](#adding-a-list-of-gene-sets-into-fdog)
* [Bugs](#bugs)
* [How to cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*fDOG* tool is distributed as a python package called *fdog*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

## Install the fDOG package
You can install *fdog* using `pip`:
```
python3 -m pip install fdog
```

or, in case you do not have admin rights, and don't use package systems like Anaconda to manage environments you need to use the `--user` option:
```
python3 -m pip install --user fdog
```

and then add the following line to the end of your **~/.bashrc** or **~/.bash_profile** file, restart the current terminal to apply the change (or type `source ~/.bashrc`):

```
export PATH=$HOME/.local/bin:$PATH
```

## Setup fDOG

After installing *fdog*, you need to setup *fdog* to get its dependencies and pre-calculated data.

**NOTE**: in case you haven't installed [greedyFAS](https://github.com/BIONF/FAS) before, it will be installed automatically within *fDOG* setup. However, you need to run [setupFAS](https://github.com/BIONF/FAS/wiki/setupFAS) after *fDOG* setup finished before actually using *fDOG*! 

You can setup fDOG by running this command
```
fdog.setup -o /output/path/for/fdog/data
```
or, in case you are using Anaconda
```
fdog.setup -o /output/path/for/fdog/data --conda
```

*You should have the sudo password ready, otherwise some missing dependencies cannot be installed. See [dependency list](#dependencies) for more info. If you do not have root privileges, ask your admin to install those dependencies using `fdog.setup --lib` command.*

[Pre-calculated data set](https://github.com/BIONF/fDOG/wiki/Input-and-Output-Files#data-structure) of fdog will be saved in `/output/path/for/fdog/data`. After the setup run successfully, you can start using *fdog*. **Please make sure to check if you need to run [setupFAS](https://github.com/BIONF/FAS/wiki/setupFAS) first.**

You will get a warning if any of the dependencies are not ready to use, please solve those issues and rerun `fdog.setup`.

*For debugging the setup, please create a log file by running the setup as e.g. `fdog.setup | tee log.txt` for Linux/MacOS or `fdog.setup --conda | tee log.txt` for Anaconda and send us that log file, so that we can trouble shoot the issues. Most of the problems can be solved by just re-running the setup.*

# Usage
*fdog* will run smoothly with the provided sample input file 'infile.fa' if everything is set correctly.

```
fdog.run --seqFile infile.fa --seqName test --refspec HUMAN@9606@3
```
The output files with the prefix `test` will be saved at your current working directory.
You can have an overview about all available options with the command
```
fdog.run -h
```

Please find more information in [our wiki](https://github.com/BIONF/fDOG/wiki) to learn about the [input and outputs files](https://github.com/BIONF/fDOG/wiki/Input-and-Output-Files) of *fdog*.

# fDOG data set

Within the data package we provide a set of 78 reference taxa. They can be automatically downloaded during the setup. This data comes "ready to use" with the *fdog* framework. Species data must be present in the three directories listed below:

* genome_dir (Contains sub-directories for proteome fasta files for each species)
* blast_dir (Contains sub-directories for BLAST databases made with `makeblastdb` out of your proteomes)
* weight_dir (Contains feature annotation files for each proteome)

For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version])

*fdog* is not limited to those 78 taxa. If needed the user can manually add further gene sets (multiple fasta format) using provided functions.

## Adding a new gene set into fDOG
For adding **one gene set**, please use the `fdog.addTaxon` function:
```
fdog.addTaxon -f newTaxon.fa -i tax_id [-o /output/directory] [-n abbr_tax_name] [-c] [-v protein_version] [-a]
```

in which, the first 3 arguments are required including `newTaxon.fa` is the gene set that need to be added, `tax_id` is its NCBI taxonomy ID, `/output/directory` is where the sub-directories can be found (*genome_dir*, *blast_dir* and *weight_dir*). If not given, new taxon will be added into the same directory of pre-calculated data. Other arguments are optional, which are `-n` for specify your own taxon name (if not given, an abbriviate name will be suggested based on the NCBI taxon name of the input `tax_id`), `-c` for calculating the BLAST DB (only needed if you need to include your new taxon into the list of taxa for compilating the core set), `-v` for identifying the genome/proteome version (default will be 1), and `-a` for turning off the annotation step (*not recommended*).

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

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/fDOG/issues/new) or be in touch via email.

# How to cite
Ebersberger, I., Strauss, S. & von Haeseler, A. HaMStR: Profile hidden markov model based search for orthologs in ESTs. BMC Evol Biol 9, 157 (2009), [doi:10.1186/1471-2148-9-157](https://doi.org/10.1186/1471-2148-9-157)

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Vinh Tran](https://github.com/trvinh)
- [Holger Bergmann](https://github.com/holgerbgm)

# Contact
For further support or bug reports please contact: ebersberger@bio.uni-frankfurt.de
