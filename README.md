# Introgression tests from gene trees

# Introduction
#
This script performs introgression tests in 4-taxa oriented trees, in a similar fashion as ABBA/BABA tests [1,2], but using gene trees as input data instead of SNPs. Thus, it is usefull to investigate introgression signals in large phylogenomic datasets. See the "Rationale" part for more informations on how it works.

# Installation and usage
#
This script is written in the R language, and based on the ape package [3], which must be installed. No other installation is needed, just source the script in R.

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

# Rationale
#
This test is an adaptation of the Patterson's D statistic test, or ABBA/BABA test, to be used on gene trees from unlinked loci rather than bi-allelic markers. 
We consider a 4-taxa oriented species tree noted (outgroup, (H1, (H2, H3))). When considering single genes genealogies, two branching patterns are possible in addition to the one concordant with the species tree : either D1=(H2, (H1,H3)) or D2=(H3, (H1, H2)). For unlinked loci, under stochastic processes alone (i.e. ILS), both discordant patterns are expected to be present in the same proportion in the gene trees pool. Thus, an excess of one of these patterns indicates introgression. In a similar way as the ABBA/BABA test, we here calculate a statistic to describe the difference of proportions between both patterns, as follow `(D1-D2)/(D1+D2)`. In the abscence of introgression, this statistic is expected to be equal to zero.

# References
#
[1]
[2]
[3] Ape
