# phylo-fasta-fudger

This tool "fills in" or "fudges" individual-specific missing data in phylogenomic datasets by substituting missing individuals' data with that of their closest relative (determined from a reference phylogeny). This can facilitate analyses that require no missing data (e.g. using gene trees to estimate a species tree where all gene trees must have the same taxa).

The input for **phylo-fasta-fudger** is:
1. A directory of nucleotide/peptide fasta files (one per gene, and they must end with `.fasta` suffix). Although not required, using nucleotide/peptide alignments as input facilitates downstream analyses.
2. A reference phylogeny in newick format that contains at least all taxa represented in the gene-specific alignments (taxon names must match between the tree and fasta files). Most easily this could be a well-supported tree from concatenated analysis (supermatrix) of the same dataset (thus containing all of the same taxa). A reference tree derived from different data will also work, however if the tree contains more taxa than the fasta files, all taxa from the tree will be filled in with the data of their closest relative based on the relationships in the tree.

Closest relatives are determined by summing the branch lengths between pairs of taxa using [distance calculations](https://github.com/biopython/biopython/blob/master/Bio/Phylo/BaseTree.py#L482) in the [Phylo](http://biopython.org/wiki/Phylo) package of BioPython. If a taxon has multiple closest relatives of equal distance based on the tree, a random relative is chosen. The random seed can be specified as an argument to recreate results of previous runs.

## Install

This script was written in python 3.6 (could work with earlier, haven't tested) and requires [biopython](http://biopython.org/).

An environment.yml file is included for installation as a conda env (probably overkill, but whatever).

 - Install [Anaconda](https://www.anaconda.com/download/)
 - download this git repo, unzip (or equivalent), navigate to the root directory of this repo.
 - run: `conda env create -f environment.yml`
 - run: `source activate phylo-fasta-fudger`

## Usage
 
```
usage: fudger.py [-h] [-f FASTA_PATH] [-o OUTPUT_PATH] [-t TREE_FN] [-s SEED]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA_PATH, --fasta_path FASTA_PATH
                        path to input fasta directory
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        path to output fasta directory
  -t TREE_FN, --tree_fn TREE_FN
                        path to newick tree file
  -s SEED, --seed SEED  optionally set seed to make random repeatable
```

Note: if OUTPUT_PATH directory is empty or non-existant it will be used/created, otherwise an adjacent directory will be created with an automatically incrementing suffix such as `_01`.

Example input data from Dupuis *et al.* (*in review*). See also [here](https://github.com/popphylotools/HiMAP).

## Logs

phylo-fasta-fudger produces logs of the form:
```
INFO:phylo-fasta-fudger:random_seed=4321145145430681311
INFO:phylo-fasta-fudger:species_list=['Aobl', 'Afra', 'Asus', 'Ccap', 'Bcur', 'Bmin', 'Bole', 'Bjar', 'Btry', 'Blat', 'Bdor', 'Bcor', 'Bzon']
INFO:phylo-fasta-fudger:file=orth10136_1834-2121.padded.fasta
INFO:phylo-fasta-fudger:missing=Aobl - substituting=Afra - distance=14.227426000000001 - next_nearest=Asus - next_distance=14.228821
INFO:phylo-fasta-fudger:missing=Bmin - substituting=Bole - distance=8.895795999999999 - next_nearest=Bcur - next_distance=14.586389999999998
INFO:phylo-fasta-fudger:missing=Bjar - substituting=Bdor - distance=5.712655 - next_nearest=Bzon - next_distance=5.715592999999999
INFO:phylo-fasta-fudger:missing=Btry - substituting=Bdor - distance=5.711662 - next_nearest=Bzon - next_distance=5.714599999999999
INFO:phylo-fasta-fudger:missing=Blat - substituting=Bdor - distance=0.048722999999999995 - next_nearest=Bzon - next_distance=0.051661
INFO:phylo-fasta-fudger:file=orth2374_1109-1760.padded.fasta
INFO:phylo-fasta-fudger:file=orth10425_453-875.padded.fasta
INFO:phylo-fasta-fudger:missing=Aobl - substituting=Ccap - distance=62.569479 - next_nearest=Bcur - next_distance=84.552573
INFO:phylo-fasta-fudger:missing=Afra - substituting=Ccap - distance=76.787487 - next_nearest=Bcur - next_distance=98.77058099999999
INFO:phylo-fasta-fudger:missing=Asus - substituting=Ccap - distance=76.788882 - next_nearest=Bcur - next_distance=98.771976
INFO:phylo-fasta-fudger:missing=Bmin - substituting=Bole - distance=8.895795999999999 - next_nearest=Bcur - next_distance=14.586389999999998
INFO:phylo-fasta-fudger:missing=Bjar - substituting=Btry - distance=0.019292999999999998 - next_nearest=Blat - next_distance=5.707153999999999
```
