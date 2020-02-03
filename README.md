# Introgression tests from gene trees

# Introduction
#
This script performs introgression tests in 4-taxa oriented trees, in a similar fashion as ABBA/BABA tests, but using gene trees as input data instead of SNPs. Thus, it is usefull to investigate introgression signals in large phylogenomic datasets. See the "Rationale" part for more informations on how it works.

# Installation and usage
#
This script is written in the R language, and based on the ape package [1], which must be installed. No other installation is needed, just source the script in R.

The function can be called as follow : `ABBA_GT(treeset, taxa, Nbs, bplot)`

* `treeset` : phylogenetic trees contained in a multiphylo object (see ape documentation);
* `taxa` : a vector containing the labels of the focal taxa, in the following order : (outgroup, (H1, (H2, H3)));
* `Nbs` : number of bootstrap replicates;
* `bplot` : logical, whether to plot the pattern counts or not (default=T).

The output is a list of several items

 * `$pruned_trees` contains the phylogenetic trees retained for the test
 * `$patterns_counts` contains the counts for the three alternative branching patterns
 * `$D.stat` contains the D statistic (cf. below)
 * `$bootstrap_results$replicates` contains the D statistic for each of the bootstrap replicates
 * `$bootstrap_results$p-value` contains the p-value that Dstat = 0.
 
 Additionnaly, if `bplot=T`, the script returns a plot displaying the distribution of the D-stat for all the bootstrap replicates (left pannel), and a barplot of the counts of the three branching patterns.   
