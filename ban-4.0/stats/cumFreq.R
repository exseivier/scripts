### UTRs Analysis

# Parameters settings
pdf("cumFreq.pdf", width = 10, height = 10)
library(Hmisc)
par(mfrow=c(1,2))


### Load Data
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

dev.off()
# Final plotting
