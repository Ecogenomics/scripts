# treecompare

Comparison of phylogenetic tree topology and taxonomic differences between trees.

This package is not recommended for public use. It was developed for some of the
analyses in the following manuscript: 
[*A rank-normalized archaeal taxonomy based on genome phylogeny resolves widespread incomplete and uneven classifications*](https://www.biorxiv.org/content/10.1101/2020.03.01.972265v1).

## Usage
After installing treecompare using either pip or conda, run the following command:

 `treecompare --trees example/input_trees.tsv --out_dir /tmp/treecompare --cpus 1`
 
 Note that the `trees` input file contains two columns (tab delimited). The first
 is the unique id for the tree, and the second is the absolute path to the **rooted tree**.


 
 
 ## Releases
 #### 0.1.0
 * Initial release.
