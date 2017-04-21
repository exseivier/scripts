#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly=T)

data.tgs <- read.table(args[1], header=T, row.names=1)
output.suffix <- args[2]
head(data.tgs)
data.tgs$Percent_on_L <- data.tgs$X/data.tgs$M
data.tgs$Percent_on_S <- data.tgs$Y/data.tgs$N
pdf(paste(args[1], "_", args[2], ".pdf", sep=""))
plot(data.tgs$Percent_on_L, data.tgs$Percent_on_S)
dev.off()
