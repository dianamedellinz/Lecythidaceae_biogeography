#!/usr/bin/env Rscript

# Rscript reroot_genetrees.R [genetrees.treefile] [genes_tree_order.txt] [outgroup_list.txt] [output prefix] [output suffix]
# for Astral tree:
# Rscript reroot_genetrees.R [astral.treefile] [file with one line with name you want] [outgroup_list.txt]
library(ape)
library(phytools)

#parse arguments
argv <- commandArgs(trailingOnly = TRUE)
geneTreesFile <- argv[1]
geneListFile <- argv[2]
outgroupFile <- argv[3]
outputPrefix <- argv[4]
outputSuffix <- argv[5]

#read in files
geneTrees <- read.tree(geneTreesFile)
locus.df <- read.table(geneListFile)
loci <- locus.df[,1]
outgroup.df <- read.table(outgroupFile)
outgroups <- outgroup.df[,1]

reroot_and_ladderize_tree <- function(inputTree, rootNode) {
  if (is.rooted(inputTree) == TRUE) {
    inputTree <- unroot(inputTree)
  }
  rootedTree <- reroot(inputTree, rootNode, position=0.5 * inputTree$edge.length[which(inputTree$edge[,2] == rootNode)])
  ladderizedTree <- ladderize(rootedTree, right=TRUE)
  return(ladderizedTree)
}

if (class(geneTrees) == "phylo") {
  gt <- geneTrees
  foundOGs <- gt$tip.label[!is.na(match(gt$tip.label,outgroups))]
  if (length(foundOGs) <= 1) {
    rrNode <- which(gt$tip.label == foundOGs[1])
  } else {
    rrNode <- getMRCA(gt, tip=foundOGs)
  }
  gt.rerooted <- reroot_and_ladderize_tree(gt, rrNode)
  fname<-paste0(loci[1])
  write.tree(gt.rerooted,file=fname)
} else {
#reroot gene trees and output each as a separate treefile
  for (g in 1:length(geneTrees)) {
    gt <- geneTrees[[g]]
    foundOGs <- gt$tip.label[!is.na(match(gt$tip.label,outgroups))]
    if (length(foundOGs) == 1) {
      rrNode <- which(gt$tip.label == foundOGs[1])
    } else if (length(foundOGs) == 0) {
      print(paste("Gene",loci[g],"is missing all outgroups. Try rerooting manually.",sep=" "))
      next
    } else {
      rrNode <- getMRCA(gt, tip=foundOGs)
    }

    gt.rerooted <- reroot_and_ladderize_tree(gt, rrNode)
    fname<-paste0(outputPrefix,loci[g],outputSuffix)
    write.tree(gt.rerooted,file=fname)
 }
}
