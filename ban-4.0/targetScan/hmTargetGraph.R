#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly=T)

data.tgs <- read.table(args[1], header=T, row.names=1)
head(data.tgs)
pdf(paste(args[1], ".pdf", sep=""))
plot(data.tgs$Percent_on_L, data.tgs$Percent_on_S)
dev.off()
