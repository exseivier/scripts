#!/usr/bin/env Rscript
library(edgeR)
library(locfit)

rm(list=ls())
args <- commandArgs(trailingOnly=T)
selected.tissues <- eval(parse(text=args[2]))
load(paste(".", args[1], ".RData", sep=""))
selected.tissues
selected.tissues.len <- length(selected.tissues)
selected.tissues.len
indexes <- match(selected.tissues, treat.tiss)
indexes
remaining.tissues <- treat.tiss[-indexes]
remaining.tissues
remaining.tissues.len <- length(remaining.tissues)
paste(selected.tissues, collapse="+")
#class(eval(parse(text=paste(selected.tissues, collapse="+"))))
my.ctst <- makeContrasts(cts=((eval(parse(text=paste(selected.tissues, collapse="+")))/selected.tissues.len)-(eval(parse(text=paste(remaining.tissues, collapse="+")))/remaining.tissues.len)), levels=design.tiss)
my.ctst

tissues.qlf.cts <- glmQLFTest(tissues.fit, contrast=my.ctst[,"cts"])

tissues.qlf.cts$table$PValue_log <- abs(log10(tissues.qlf.cts$table$PValue))
tissues.qlf.cts$table$PValue_log_sig <- ifelse(tissues.qlf.cts$table$logFC < 0, tissues.qlf.cts$table$PValue_log*-1, tissues.qlf.cts$table$PValue_log)
tissues.qlf.cts$table <- tissues.qlf.cts$table[order(tissues.qlf.cts$table$logFC),]
tissues.qlf.cts$OOtable <- tissues.qlf.cts$table[order(tissues.qlf.cts$table$PValue_log_sig),]

name.splitted <- strsplit(args[1], "\\.")
name.splitted <- name.splitted[[1]]

write.table(tissues.qlf.cts$OOtable, file=paste(name.splitted[1], paste(selected.tissues, collapse="_"), paste(name.splitted[3:length(name.splitted)], collapse="."), sep="."))


