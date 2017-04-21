#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
olddir <- getwd()
data.path <- args[1]
n <- as.integer(args[2])
setwd(data.path)
print(paste("Getting in in ", getwd(), sep=""))
lf <- list.files(pattern=".data$")
columns <- c(5)
total <- length(lf)
rows <- total/columns
pdf("targetScan.pdf", width=columns*4, height=rows*3)
par(mfrow=c(rows, columns))
par(mar=c(3,5,3,5))
for (file in lf) {
	print(paste("Accessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene", "cs++")
	rownames(table) <- 1:length(table[,1])
	array_of_genes <- as.character(table$gene)
	N <- length(array_of_genes)
	window <- n <- N/15
	k <- N/2
	x_axis <- c()
	y_axis <- c()
	y2_axis <- c()
	for (i in seq(1, N, by=n)) {
		subarray_of_genes <- array_of_genes[i:window]
		print(subarray_of_genes)
		print(length(subarray_of_genes))
		x <- length(grep(".S", subarray_of_genes))
		prb <- phyper(x, k, N-k, n)
		x_axis <- c(x_axis, window -  (n/2))
		y_axis <- c(y_axis, prb)
		y2_axis <- c(y2_axis, (x/n)*100)
		window = window + n
	}
	#print(x_axis)
	#print(y_axis)
	q1 <- (25*N)/100
	mode <- (50*N)/100
	q2 <- (75*N)/100
	plot(x_axis, y_axis, main=file, ylim=c(0.0, 1.0), xlab="", ylab="", type="l", axes=F)
	axis(2, ylim=c(0.0, 1.0), col="black", las=1)
	mtext("P(X <= x)", side=2, line=2.5)
	box()
	par(new=T)
	plot(x_axis, y2_axis, col="blue", ylim=c(0, 100), xlab="", ylab="", axes=F)
	axis(4, ylim=c(0,100), col="blue", las=1)
	mtext("Percent of L genes", side=4, line=2.5)
	axis(1, pretty(1:N))
	mtext("Genes", side=1, line=2.5, col="black")
	abline(v=q1, col="red")
	abline(v=mode, col="red")
	abline(v=q2, col="red")
}
dev.off()

.silenced <- function (){
pdf("targetScan_boxplot.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for (file in lf) {
	print(paste("Accessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene", "cs++")
	rownames(table) <- 1:length(table[,1])
	array_of_genes <- as.character(table$gene)
	for (i in c(0, 25, 50, 75)) {
		subarray_of_genes <- array_of_genes[i:i+25]
		L <- length(grep(".L", subarray_of_genes))
		S <- length(grep(".S", subarray_of_genes))

	}
	#print(x_axis)
	#print(y_axis)
	q1 <- (25*N)/100
	mode <- (50*N)/100
	q2 <- (75*N)/100
	plot(x_axis, y_axis, main=file, xlab="Genes", ylab="P(X <= x)", type="l")
	abline(v=q1, col="red")
	abline(v=mode, col="red")
	abline(v=q2, col="red")
}
dev.off()
} #END .silenced


lf <- list.files(pattern=".sorted$")
pdf("targetScan_paralogs.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for (file in lf) {
	print(paste("Accsessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene_name", "geneL", "csL", "geneS", "csS")
	rownames(table) <- 1:length(table[,1])
	table <- table[order(table$csL),]
	x_axis <- 1:length(table[,1])
	y_axis <- table$csS
	plot(x_axis, y_axis, main=file, xlab="L genes", ylab="CS++ for S genes")
	abline(h=median(y_axis), col="red")
	abline(h=median(table$csL), col="blue")
}
dev.off()

lf <- list.files(pattern=".sorted$")
pdf("targetScan_paralogs_boxplot.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for (file in lf) {
	print(paste("Accsessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene_name", "geneL", "csL", "geneS", "csS")
	rownames(table) <- 1:length(table[,1])
	table <- table[order(table$csL),]
	N <- length(table[,1])
	q1 <- (25*N)/100
	median <- (50*N)/100
	q2 <- (75*N)/100
	q1.table <- table[0:q1,]
	median.table <- table[q1:median,]
	q2.table <- table[median:q2,]
	last.table <- table[q2:N,]
	y_min <- min(min(table$csL), min(table$csS))
	y_max <- max(max(table$csL), max(table$csS))
	boxplot(q1.table$csL, main=file, at=1, xlim=c(0,8), ylim=c(y_min - 1, y_max + 1), col="red")
	boxplot(q2.table$csS, at=2, add=T, col="red")
	boxplot(median.table$csL, at=3, add=T, col="blue")
	boxplot(median.table$csS, at=4, add=T, col="blue")
	boxplot(q2.table$csL, at=5, add=T, col="violet")
	boxplot(q2.table$csS, at=6, add=T, col="violet")
	boxplot(last.table$csL, at=7, add=T, col="green")
	boxplot(last.table$csS, at=8, add=T, col="green")
}
dev.off()

setwd(olddir)
print(paste("Returning to ", getwd(), sep=""))


