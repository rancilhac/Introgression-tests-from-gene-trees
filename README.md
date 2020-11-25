# Introgression tests from gene trees

# Introduction
#
This script performs introgression tests in 4-taxa oriented trees, in a similar fashion as ABBA/BABA tests [1,2], but using gene trees as input data instead of SNPs. Thus, it is usefull to investigate introgression signals in large phylogenomic datasets. See the "Rationale" part for more informations on how it works.

When using this script, please cite: Rancilhac, L., Irisarri, I., Angelini, C., Arntzen, J. W., Babik, W., Bossuyt, F., Künzel, S., Lüddecke, T., Pasmans, F., Sanchez, E., Weisrock, D., Veith, M., Wielstra, B., Steinfarz, S., Hofreiter, M., Philippe, H., & Vences. (2020). Phylotranscriptomic evidence for pervasive ancient hybridization among Old World salamanders. Molecular Phylogenetics and Evolution, 106967. https://doi.org/10.1016/j.ympev.2020.106967

# Installation and usage
#
This script is written in R language, and based on the ape package [3], which must be installed. No other installation is needed, just source the script in R.

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

We consider a 4-taxa oriented species tree noted (outgroup, (H1, (H2, H3))). When considering single genes genealogies, two branching patterns are possible in addition to the one concordant with the species tree : either D1=(H3, (H1,H2)) or D2=(H2, (H1, H3)). For unlinked loci, under stochastic processes alone (i.e. ILS), both discordant patterns are expected to be present in the same proportion in the gene trees pool. Thus, an excess of one of these patterns indicates introgression. 

In a similar way as the ABBA/BABA test, we here calculate a statistic to describe the difference of proportions between both patterns, as follow `(D1-D2)/(D1+D2)`. In the abscence of introgression, this statistic is expected to be equal to zero. If the statistic is greater than zero, it indicates an excess of D1, and thus intregression between H1 and H2. On the contrary, if it is negative, this indicates an excess of D2 and introgression between H2 and H3. The significance of the statistic is assessed with bootstraping, using a user-specified number of replicates.

Since phylogenetic trees can give misleading signal due to various artifacts, we recommand to collapse the poorly supported nodes of the gene trees. The function will automaticaly remove the trees for which the focal nodes are not resolved.

# References
#
[1] Green, R.E., Krause, J., Briggs, A.W., Maricic, T., Stenzel, U., Kircher, M.,[...], Pääbo, S., 2010. A draft sequence of the Neandertal genome. Science 328, 710–722. https://dx.doi.org/10.1126/science.1188021 

[2] Durand, E.Y., Patterson, N., Reich, D., Slatkin, M., 2011. Testing for ancient admixture between closely related populations. Mol. Biol. Evol. 28, 2239–2252. https://doi.org/10.1093/molbev/msr048 

[3] Paradis, E., Claude, J., Strimmer, K., 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20, 289–290. https://doi.org/10.1093/bioinformatics/btg412
