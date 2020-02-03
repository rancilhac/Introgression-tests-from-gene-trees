# Introgression tests from gene trees

# Introduction
#
This script performs introgression tests in 4-taxa oriented trees, in a similar fashion as ABBA/BABA tests, but using gene trees as input data instead of SNPs. Thus, it is usefull to investigate introgression signals in large phylogenomic datasets.

# Installation and usage
#
This script is written in the R language, and based on the ape package [1], which must be installed. No other installation is needed, just source the script in R.

The function can be called as follow : `ABBA_GT(treeset, taxa, Nbs, bplot)`

* `treeset` : phylogenetic trees contained in a multiphylo object (see ape documentation);
* `taxa` : a vector containing the labels of the focal taxa, in the following order : (outgroup, (H1, (H2, H3)));
* `Nbs` : number of bootstrap replicates;
* `bplot` : logical, whether to plot the pattern counts or not (default=T).

