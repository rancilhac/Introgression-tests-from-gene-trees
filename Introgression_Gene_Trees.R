############### Perform introgression test from sets of gene trees ###############
## Loïs Rancilhac, 04/2019
## Dependencies : ape

library(ape)

## Define function to calculate dstat
# This function counts the patterns and calculte the Dstat from a table where the lines are the individuals and the columns are SNPs
Dstat <- function(table){
  ABBA <- 0
  BABA <- 0
  for(j in seq(1, ncol(table))){
    if(table[2,j] == "B" & table[3,j] == "B") { ABBA <- ABBA+1 }
    else if(table[1,j] == "B" & table[3,j] == "B"){BABA <- BABA+1 }
  }
  return((ABBA - BABA)/(ABBA + BABA))
}

## Bootstrapping function
# This function generate N bootstrap replicates from a SNP table and assess the significativity of the observed Dstat
BT_sign <- function(table.bs, N, obs){
  D.stat <- c()
  for(i in seq(1, N)){
    bootab <- table.bs[, sample(ncol(table.bs), replace=T)]
    cur.D.stat <- Dstat(bootab)
    D.stat <- c(D.stat, cur.D.stat)
  }
  BT_output <- list(D.stat, abs(obs/sd(D.stat)), 2*(1-pnorm(abs(obs/sd(D.stat)))))
  names(BT_output) <- c("replicates", "z", "p-value")
  return(BT_output)
}

## Full test function
# This function performs the introgression test from a set of gene trees and 4 taxa, and performs bootstrap to assess the signigicativity of the Dstat

ABBA_GT <- function(treeset, taxa, Nbs, bplot=T) {
  cat("Outgroup = ", taxa[1], "\n", sep='')
  cat("H1 = ", taxa[2],"\n", sep='')
  cat("H2 = ", taxa[3],"\n", sep='')
  cat("H3 = ", taxa[4],"\n", sep='')
  #select only the gene trees that contains the 4 taxa of interest
  complete_trees <- list()
  selected_trees <- c()
  it=0
  for(i in seq(1, length(treeset))){
    if(isTRUE(length(grep(taxa[1], treeset[[i]]$tip.label)) > 0 ) & 
       isTRUE(length(grep(taxa[2], treeset[[i]]$tip.label)) > 0) &
       isTRUE(length(grep(taxa[3], treeset[[i]]$tip.label)) > 0) &
       isTRUE(length(grep(taxa[4], treeset[[i]]$tip.label)) > 0)) {
      it <- it+1
      complete_trees[[it]] <- treeset[[i]]
      selected_trees <- c(selected_trees, i)
    }
  }
  class(complete_trees) <- "multiPhylo"
  # extract the 4 taxa of interrest from the selected gene trees
  pruned_trees <- list()
  for(i in seq(1, length(complete_trees))){
    cur_ptree <- keep.tip(complete_trees[[i]], taxa)
    pruned_trees[[i]] <- cur_ptree
  }
  class(pruned_trees) <- "multiPhylo"
  # root the gene trees with the outgroup
  pruned_trees <- root(pruned_trees, taxa[1], resolve.root=T)
  # check that the trees are resolved, remove the ones that are not
  unresolved <- c()
  for(i in seq(1, length(pruned_trees))){
    if(pruned_trees[[i]]$Nnode < 3) { unresolved <- c(unresolved, i) }
  }
  pruned_trees <- pruned_trees[-unresolved]
  class(pruned_trees) <- "multiPhylo" 
  print(pruned_trees)                          
  # initiate variables for patterns counts
  Out <- c()
  H1 <- c()
  H2 <- c()
  H3 <- c()
  AABB_count <- 0
  ABBA_count <- 0
  BABA_count <- 0
  # translation from trees to patterns and counts
  for(i in seq(1,length(pruned_trees))) {
    if(isTRUE(is.monophyletic(pruned_trees[[i]], taxa[c(3,4)]))) {
      Out <- c(Out, "A")
      H1 <- c(H1, "A")
      H2 <- c(H2, "B")
      H3 <- c(H3, "B")
      AABB_count <- AABB_count+1
    }
    else if(isTRUE(is.monophyletic(pruned_trees[[i]], taxa[c(2,3)]))) {
      Out <- c(Out, "A")
      H1 <- c(H1, "B")
      H2 <- c(H2, "B")
      H3 <- c(H3, "A")
      ABBA_count <- ABBA_count+1
    }
    else if(isTRUE(is.monophyletic(pruned_trees[[i]], taxa[c(2,4)]))) {
      Out <- c(Out, "B")
      H1 <- c(H1, "A")
      H2 <- c(H2, "B")
      H3 <- c(H3, "A")
      BABA_count <- BABA_count+1
    }
  }
  # store the patterns in a table
  ABBA_table <- rbind(Out, H1, H2, H3)
  count_table <- cbind(AABB_count, BABA_count, ABBA_count)
  colnames(count_table) <- c("AABB", "BABA", "ABBA")
  # calculate the D statistic
  Dstat <- function(table){
    ABBA <- 0
    BABA <- 0
    for(j in seq(1, ncol(table))){
      if(table[2,j] == "B" & table[3,j] == "B") { ABBA <- ABBA+1 }
      else if(table[1,j] == "B" & table[3,j] == "B"){BABA <- BABA+1 }
    }
    return((ABBA - BABA)/(ABBA + BABA))
  }
  D.stat <- (ABBA_count - BABA_count)/(ABBA_count + BABA_count)
  # perform bootstrap
  D.stat.bs <- c()
  for(i in seq(1, Nbs)){
    ABBA <- 0
    BABA <- 0
    bootab <- ABBA_table[, sample(ncol(ABBA_table), replace=T)]
    cur.D.stat <- Dstat(bootab)
    D.stat.bs <- c(D.stat.bs, cur.D.stat) 
  }
  BT_output <- list(D.stat.bs, 2*(1-pnorm(abs(D.stat/sd(D.stat.bs)))))
  names(BT_output) <- c("replicates", "p-value")
  
  if(bplot == T) {
    par(mfrow=c(1,2))
    xdens <- paste("Dstat (p-value = ", BT_output[[2]], ")", sep='')
    maindens <- paste("density distribution of the D statistic", "\n", "for ", Nbs, " bootstrap replicates", sep='')
    plot(density(BT_output[[1]]), main=maindens, xlab=xdens)
    abline(v=D.stat)
    barplot(count_table, main="patterns counts ", 
            names.arg=c(paste(taxa[3], taxa[4], sep="\n"), paste(taxa[2], taxa[4], sep="\n"), paste(taxa[2], taxa[3], sep="\n")))
  }
  ABBA_list <- list(pruned_trees, count_table, D.stat, BT_output)
  names(ABBA_list) <- c("pruned_trees", "patterns_counts", "D.stat", "bootstrap_results")
  return(ABBA_list)
  
}
