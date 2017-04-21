#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
xlae.rdc <- read.table(args[1], header=F, row.names=1)
Lxtro.rdc <- read.table(args[2], header=F, row.names=1)
Sxtro.rdc <- read.table(args[3], header=F, row.names=1)


w.lsvlx <- wilcox.test(xlae.rdc[,1], Lxtro.rdc[,1], alternative=c("two.sided"))
w.lsvsx <- wilcox.test(xlae.rdc[,1], Sxtro.rdc[,1], alternative=c("two.sided"))
w.lxvsx <- wilcox.test(Lxtro.rdc[,1], Sxtro.rdc[,1], alternative=c("two.sided"))

pdf("boxplot_rdc.pdf", width=4, height=4)
par(mfrow=c(1,1))
data <- cbind(xlae.rdc[,1], Lxtro.rdc[,1], Sxtro.rdc[,1])
boxplot(data, names=c("L Vs S", "L Vs Xtro", "S Vs Xtro"))
mtext(paste("w.lsvlx = ", w.lsvlx$p.value, "\n", "w.lsvsx = ", w.lsvsx$p.value, "\n", "w.lxvsx = ", w.lxvsx$p.value, sep=""), side=1, line=4, at=1)
dev.off()

pdf("hist_rdc.pdf", width=3.5, height=10)
par(mfrow=c(3,1))
hist(xlae.rdc[,1], main="L Vs S", breaks=50)
abline(v=0.5, col="red")
abline(v=median(xlae.rdc[,1]), col="black", lwd=2)
hist(Lxtro.rdc[,1], main="L Vs Xtro", breaks=50)
abline(v=0.5, col="red")
abline(v=median(Lxtro.rdc[,1]), col="black", lwd=2)
hist(Sxtro.rdc[,1], main="S Vs Xtro", breaks=50)
abline(v=0.5, col="red")
abline(v=median(Sxtro.rdc[,1]), col="black", lwd=2)

dev.off()
