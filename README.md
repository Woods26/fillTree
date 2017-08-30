# fillTree

Simple script to fill out fasta files with missing species by substituting the closest avalible reletive using a newick tree for reference.

## Install

This script was written in python 3.6 (could work with earlier, havn't tested) and requires [biopython](http://biopython.org/).

An environment.yml file is included for installation as a conda env (probably overkill, but whatever).

 - Install [Anaconda](https://www.anaconda.com/download/)
 - download this git repo, unzip (or equivalent), navigate to the root directory of this repo.
 - run: `conda env create -f environment.yml`
 - run: `source activate fillTree`
 
## Usage
 
```
usage: fillTree.py [-h] [-f FASTA_PATH] [-o OUTPUT_PATH] [-t TREE_FN]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA_PATH, --fasta_path FASTA_PATH
                        path to input fasta directory
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        path to output fasta directory (be careful, will
                        delete contents
  -t TREE_FN, --tree_fn TREE_FN
                        path to newick tree file
```
