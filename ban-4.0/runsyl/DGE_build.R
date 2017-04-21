#!/usr/bin/env Rscript

#Setting working directory
#setwd("/home/montalvo/scripts/ban-4.0/Rscripts")

rm(list=ls())

args <- commandArgs(trailingOnly=T)

### Loading edgeR & locfit libraries
library(edgeR)
library(locfit)

### Loading tissues gene expression data
tissues <- read.delim(args[1], header=T, row.names="ID")
head(tissues)

### Building DGE Object
groups.tiss <- factor(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
design.tiss <- model.matrix(~0+groups.tiss)
treat.tiss <- c("brain", "eye", "heart", "intestine", "kidney", "liver", "lung", "muscle", "ovary", "pancreas", "skin", "spleen", "stomach", "testis")
rownames(design.tiss) <- colnames(tissues)
colnames(design.tiss) <- treat.tiss
design.tiss
tissues.DGE <- DGEList(counts=tissues, group=groups.tiss)
tissues.DGE

### Filtering data
keep <- rowSums(cpm(tissues) > 1) >= 2
tissues <- tissues[keep,]

### Calculate Normalized Factors and estimate disperssion
tissues.DGE <- calcNormFactors(tissues.DGE)
tissues.DGE <- estimateDisp(tissues.DGE, design.tiss)

### Fit data to model Quasi-likellihood F distribution
tissues.fit <- glmQLFit(tissues.DGE, design.tiss)

save(list=ls(all.names=T), file=paste(".", args[1], ".RData", sep=""), envir=.GlobalEnv)
