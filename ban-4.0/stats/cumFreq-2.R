#!/usr/bin/env Rscript

### UTRs Analysis

# Parameters settings
library(Hmisc)
par(mfrow=c(1,1))

args <- commandArgs(trailingOnly=T)

threeUTR.xlae.L.File <- args[1]
threeUTR.xlae.S.File <- args[2]
threeUTR.xtro.File <- args[3]

### Load Data

# 3UTR Xlaevis L
threeUTR.xlae.L.table <- read.table(threeUTR.xlae.L.File, header=T, row.names="Key")
# 3UTR Xlaevis S
threeUTR.xlae.S.table <- read.table(threeUTR.xlae.S.File, header=T, row.names="Key")
# 3UTR Xtropicalis
threeUTR.xtro.table <- read.table(threeUTR.xtro.File, header=T, row.names="Key")

### Getting Data to Plotting

L.breaks <- seq(min(threeUTR.xlae.L.table$mRNA_len), 5, by=0.005)
S.breaks <- seq(min(threeUTR.xlae.S.table$mRNA_len), 5, by=0.005)
xtro.breaks <- seq(min(threeUTR.xtro.table$mRNA_len), 5, by=0.005)

L.cut <- cut(threeUTR.xlae.L.table$mRNA_len, L.breaks)
S.cut <- cut(threeUTR.xlae.S.table$mRNA_len, S.breaks)
xtro.cut <- cut(threeUTR.xtro.table$mRNA_len, xtro.breaks)

wilcoxLS <- wilcox.test(threeUTR.xlae.L.table$mRNA_len, threeUTR.xlae.S.table$mRNA_len)
wilcoxLS
wilcoxLX <- wilcox.test(threeUTR.xlae.L.table$mRNA_len, threeUTR.xtro.table$mRNA_len)
wilcoxLX
wilcoxSX <- wilcox.test(threeUTR.xlae.S.table$mRNA_len, threeUTR.xtro.table$mRNA_len)
wilcoxSX

L.freq <- table(L.cut)
S.freq <- table(S.cut)
xtro.freq <- table(xtro.cut)

L.cumfreq <- c(0, cumsum(L.freq / length(threeUTR.xlae.L.table$mRNA_len)))
S.cumfreq <- c(0, cumsum(S.freq / length(threeUTR.xlae.S.table$mRNA_len)))
xtro.cumfreq <- c(0, cumsum(xtro.freq / length(threeUTR.xtro.table$mRNA_len)))
# UTR's length cummulative frequency plot
pdf("cumFreq.pdf", width = 8, height = 8)
plot(L.breaks, L.cumfreq, type="l", col="blue", main="3' UTRs length", ylab="Frequency", xlab="Kb", xlim=c(0, 5), lwd=1)
lines(S.breaks, S.cumfreq, col="red", lwd=1)
lines(xtro.breaks, xtro.cumfreq, col="green", lwd=1)
legend("bottomright", col=c("blue", "red", "green"), c("Xlae.L", "Xlae.S", "Xtro"), lwd=2)
text(3,0.5, paste("LS_W ", wilcoxLS$p.value, "\n", "LX_W ", wilcoxLX$p.value, "\n", "SX_W ", wilcoxSX$p.value, sep=""))
dev.off()

# Seed enrichment cummulative frequency plot
seed.L.breaks <- seq(min(threeUTR.xlae.L.table$SPKM_by_mRNA), max(threeUTR.xlae.L.table$SPKM_by_mRNA), by=5)
seed.L.cut <- cut(threeUTR.xlae.L.table$SPKM_by_mRNA, seed.L.breaks)
seed.L.freq <- table(seed.L.cut)
seed.S.breaks <- seq(min(threeUTR.xlae.S.table$SPKM_by_mRNA), max(threeUTR.xlae.S.table$SPKM_by_mRNA), by=5)
seed.S.cut <- cut(threeUTR.xlae.S.table$SPKM_by_mRNA, seed.S.breaks)
seed.S.freq <- table(seed.S.cut)
seed.X.breaks <- seq(min(threeUTR.xtro.table$SPKM_by_mRNA), max(threeUTR.xtro.table$SPKM_by_mRNA), by=5)
seed.X.cut <- cut(threeUTR.xtro.table$SPKM_by_mRNA, seed.X.breaks)
seed.X.freq <- table(seed.X.cut)
x_max = max(max(threeUTR.xlae.L.table$SPKM_by_mRNA), max(threeUTR.xlae.S.table$SPKM_by_mRNA), max(threeUTR.xtro.table$SPKM_by_mRNA))
seed.L.cumfreq <- c(0, cumsum(seed.L.freq / length(threeUTR.xlae.L.table$SPKM_by_mRNA)))
seed.S.cumfreq <- c(0, cumsum(seed.S.freq / length(threeUTR.xlae.S.table$SPKM_by_mRNA)))
seed.X.cumfreq <- c(0, cumsum(seed.X.freq / length(threeUTR.xtro.table$SPKM_by_mRNA)))


wilcoxLS.sE <- wilcox.test(threeUTR.xlae.L.table$SPKM_by_mRNA, threeUTR.xlae.S.table$SPKM_by_mRNA, alternative=c("less"))
wilcoxLS.sE
wilcoxLX.sE <- wilcox.test(threeUTR.xlae.L.table$SPKM_by_mRNA, threeUTR.xtro.table$SPKM_by_mRNA, alternative=c("less"))
wilcoxLX.sE
wilcoxSX.sE <- wilcox.test(threeUTR.xlae.S.table$SPKM_by_mRNA, threeUTR.xtro.table$SPKM_by_mRNA, alternative=c("less"))
wilcoxSX.sE



pdf("seedEnrichment_cumFreq.pdf", width=8, height=8)
plot(seed.L.breaks, seed.L.cumfreq, type="l", col="blue", main="Seed enrichment", ylab="Frequency", xlab="Seed Counts", xlim=c(10, 30), lwd=1)
lines(seed.S.breaks, seed.S.cumfreq, col="red", lwd=1)
lines(seed.X.breaks, seed.X.cumfreq, col="green", lwd=1)
text(25,0.5, paste("LS_W ", wilcoxLS.sE$p.value, "\n", "LX_W ", wilcoxLX.sE$p.value, "\n", "SX_W ", wilcoxSX.sE$p.value, sep=""))
legend("bottomright", col=c("blue", "red", "green"), c("Xlae.L", "Xlae.S", "Xtro"), lwd=2)
dev.off()


.silenced_load <- function(){
LEN_3utr <- read.delim("Xla.v91.1.8.3.2_3utr.len", header=FALSE)
LEN_5utr <- read.delim("Xla.v91.1.8.3.2_5utr.len", header=FALSE)
LEN_3utr <- LEN_3utr[!duplicated(LEN_3utr[,1]),]
LEN_5utr <- LEN_5utr[!duplicated(LEN_5utr[,1]),]
rownames(LEN_3utr) <- LEN_3utr[,1]
rownames(LEN_5utr) <- LEN_5utr[,1]
LEN_3utr <- LEN_3utr[,2:3]
LEN_5utr <- LEN_5utr[,2:3]
cnames <- c("len", "sg_type")
colnames(LEN_3utr) <- cnames
colnames(LEN_5utr) <- cnames
}

.silenced_plotting <- function(){
### Plotting Cummulative frequency graphs 5' Vs 3' UTRs
LEN_3utr.Breaks <- seq(min(LEN_3utr$len), max(LEN_3utr$len), by=5)
LEN_5utr.Breaks <- seq(min(LEN_5utr$len), max(LEN_5utr$len), by=5)
LEN_3utr.cut <- cut(LEN_3utr$len, LEN_3utr.Breaks)
LEN_5utr.cut <- cut(LEN_5utr$len, LEN_5utr.Breaks)
LEN_3utr.freq <- table(LEN_3utr.cut)
LEN_5utr.freq <- table(LEN_5utr.cut)
LEN_3utr.cumfreq <- c(0, cumsum(LEN_3utr.freq / length(LEN_3utr$len)))
LEN_5utr.cumfreq <- c(0, cumsum(LEN_5utr.freq / length(LEN_5utr$len)))
plot(LEN_3utr.Breaks, LEN_3utr.cumfreq, type="l", col="blue", main = "5' & 3' UTR length distirbution", ylab="Frequency", xlab="UTR length", xlim=c(0,20000))
lines(LEN_5utr.Breaks, LEN_5utr.cumfreq, col="red")
legend(12000, 0.2, legend=c("5' UTR", "3' UTR"), col=c("red", "blue"), lty=1, cex=0.7, pt.cex=0.7, y.intersp=1.5, x.intersp=1.5)
mtext(side=3, adj=5, at=-1, line=2, "A")

### Plotting Cumulative frequency graphs 5' and 3' UTRs (L Vs S subgenome)
LEN_3utr.L <- LEN_3utr[LEN_3utr$sg_type == "L",]
LEN_3utr.S <- LEN_3utr[LEN_3utr$sg_type == "S",]
LEN_5utr.L <- LEN_5utr[LEN_5utr$sg_type == "L",]
LEN_5utr.S <- LEN_5utr[LEN_5utr$sg_type == "S",]
# Plotting 3' UTRs
LEN_3utr.L.breaks <- seq(min(LEN_3utr.L$len), max(LEN_3utr.L$len), by=5)
LEN_3utr.S.breaks <- seq(min(LEN_3utr.S$len), max(LEN_3utr.S$len), by=5)
LEN_3utr.L.cut <- cut(LEN_3utr.L$len, LEN_3utr.L.breaks)
LEN_3utr.S.cut <- cut(LEN_3utr.S$len, LEN_3utr.S.breaks)
LEN_3utr.L.freq <- table(LEN_3utr.L.cut)
LEN_3utr.S.freq <- table(LEN_3utr.S.cut)
LEN_3utr.L.cumfreq <- c(0, cumsum(LEN_3utr.L.freq / length(LEN_3utr.L$len)))
LEN_3utr.S.cumfreq <- c(0, cumsum(LEN_3utr.S.freq / length(LEN_3utr.S$len)))
#plot(LEN_3utr.L.breaks, LEN_3utr.L.cumfreq, type="l", col="blue")
#lines(LEN_3utr.S.breaks, LEN_3utr.S.cumfreq, col="red")

# Plotting 5' UTRs
LEN_5utr.L.breaks <- seq(min(LEN_5utr.L$len), max(LEN_5utr.L$len), by=5)
LEN_5utr.S.breaks <- seq(min(LEN_5utr.S$len), max(LEN_5utr.S$len), by=5)
LEN_5utr.L.cut <- cut(LEN_5utr.L$len, LEN_5utr.L.breaks)
LEN_5utr.S.cut <- cut(LEN_5utr.S$len, LEN_5utr.S.breaks)
LEN_5utr.L.freq <- table(LEN_5utr.L.cut)
LEN_5utr.S.freq <- table(LEN_5utr.S.cut)
LEN_5utr.L.cumfreq <- c(0, cumsum(LEN_5utr.L.freq / length(LEN_5utr.L$len)))
LEN_5utr.S.cumfreq <- c(0, cumsum(LEN_5utr.S.freq / length(LEN_5utr.S$len)))
plot(LEN_3utr.L.breaks, LEN_3utr.L.cumfreq, type="l", col="blue", main = "L & S UTR length distribution", ylab="Frequency", xlab="UTR length", xlim=c(0, 20000))
lines(LEN_3utr.S.breaks, LEN_3utr.S.cumfreq, col="blue", lty=2)
lines(LEN_5utr.L.breaks, LEN_5utr.L.cumfreq, col="red")
lines(LEN_5utr.S.breaks, LEN_5utr.S.cumfreq, col="red", lty=2)
legend(11000, 0.2, legend = c("L 5 UTR", "S 5 UTR", "L 3 UTR", "S 3 UTR"), col=c("red", "red", "blue", "blue"), lty=1:2, cex=0.7, pt.cex=0.7, y.intersp=1.5, x.intersp=1.5)
subplot(plot(LEN_3utr.L.breaks, LEN_3utr.L.cumfreq, type="l", col="blue", cex.lab=0.8, ylab="Frequency", xlab="UTR length", xlim=c(0,5000), ylim=c(0.6,1.0)), 13000, 0.7, size=c(2,2.7))
subplot(plot(LEN_3utr.S.breaks, LEN_3utr.S.cumfreq, type="l", col="blue", lty=2, xlab="", ylab="", xlim=c(0,5000), ylim=c(0.6,1.0)), 13000, 0.7, size=c(2,2.7))
subplot(plot(LEN_5utr.L.breaks, LEN_5utr.L.cumfreq, type="l", col="red", xlab="", ylab="", xlim=c(0,5000), ylim=c(0.6,1.0)), 13000, 0.7, size=c(2,2.7))
subplot(plot(LEN_5utr.S.breaks, LEN_5utr.S.cumfreq, type="l", col="red", lty=2, xlab="", ylab="", xlim=c(0,5000), ylim=c(0.6,1.0)), 13000, 0.7, size=c(2,2.7))
mtext(side=3, adj=5, at=-1, line=2, "B")
}

# Final plotting
