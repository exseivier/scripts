#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
counts.file <- args[1]
output.file <- args[2]
width <- as.integer(args[3])
height <- as.integer(args[4])
main <- args[5]
xlab <- args[6]
ylab <- args[7]
cummulative_freq <- function(v){
	v.breaks <- seq(min(v), max(v), by=50)
	v.cut <- cut(v, v.breaks, right=F)
	v.freq <- table(v.cut)
	v.cumfreq <- c(0, cumsum(v.freq/length(v)))
	return(list(v.cumfreq, v.breaks))
}

df <- read.table(counts.file, header=T, row.names=1)
names <- colnames(df)
colors <- sample(colors(), size=length(df), replace=F)
output <- cummulative_freq(df[,1])
pdf(output.file, width=width, height=height)
plot(output[[2]], output[[1]], main=main, xlab=xlab, ylab=ylab, type="l", col=colors[1])
for(i in 2:length(df)) {
	output <- cummulative_freq(df[,i])
	lines(output[[2]], output[[1]], col=colors[i])
}
legend("bottomright", names, col=colors, lwd=1)
dev.off()
