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
## Logs

fillTree produces logs of the form:
```
INFO:fillTree:species_list=['Aobl', 'Afra', 'Asus', 'Ccap', 'Bcur', 'Bmin', 'Bole', 'Bjar', 'Btry', 'Blat', 'Bdor', 'Bcor', 'Bzon']
INFO:fillTree:file=orth10136_1834-2121.padded.fasta
INFO:fillTree:missing=Aobl - substituting=Afra - distance=14.227426000000001 - next_nearest=Asus - next_distance=14.228821
INFO:fillTree:missing=Bmin - substituting=Bole - distance=8.895795999999999 - next_nearest=Bcur - next_distance=14.586389999999998
INFO:fillTree:missing=Bjar - substituting=Bdor - distance=5.712655 - next_nearest=Bzon - next_distance=5.715592999999999
INFO:fillTree:missing=Btry - substituting=Bdor - distance=5.711662 - next_nearest=Bzon - next_distance=5.714599999999999
INFO:fillTree:missing=Blat - substituting=Bdor - distance=0.048722999999999995 - next_nearest=Bzon - next_distance=0.051661
INFO:fillTree:file=orth2374_1109-1760.padded.fasta
INFO:fillTree:file=orth10425_453-875.padded.fasta
INFO:fillTree:missing=Aobl - substituting=Ccap - distance=62.569479 - next_nearest=Bcur - next_distance=84.552573
INFO:fillTree:missing=Afra - substituting=Ccap - distance=76.787487 - next_nearest=Bcur - next_distance=98.77058099999999
INFO:fillTree:missing=Asus - substituting=Ccap - distance=76.788882 - next_nearest=Bcur - next_distance=98.771976
INFO:fillTree:missing=Bmin - substituting=Bole - distance=8.895795999999999 - next_nearest=Bcur - next_distance=14.586389999999998
INFO:fillTree:missing=Bjar - substituting=Btry - distance=0.019292999999999998 - next_nearest=Blat - next_distance=5.707153999999999
```
