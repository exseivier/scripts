#!/usr/bin/env Rscript

### Analysing seddMap results of 3' and 5' UTRs; and mRNA and Promoter sequences
args <- commandArgs(trailingOnly=T)
# Load table
DATA_3utr <- read.delim(args[1], header=F)
DATA_5utr <- read.delim(args[2], header=F)
DATA_mRNA <- read.delim(args[3], header=F)

# row names settings
rownames(DATA_3utr) <- DATA_3utr[,1]
rownames(DATA_5utr) <- DATA_5utr[,1]
rownames(DATA_mRNA) <- DATA_mRNA[,1]

# Deleting first column
DATA_3utr <- DATA_3utr[,2:9]
DATA_5utr <- DATA_5utr[,2:9]
DATA_mRNA <- DATA_mRNA[,2:9]

# Colnames settings
cnames <- c("ID.L", "len_kb.L", "counts.L", "SPKM.L", "ID.S", "len_kb.S", "counts.S", "SPKM.S")
colnames(DATA_3utr) <- cnames
colnames(DATA_5utr) <- cnames
colnames(DATA_mRNA) <- cnames

# Deleting NA rows
DATA_3utr <- DATA_3utr[rowSums(is.na(DATA_3utr)) < 1,]
DATA_5utr <- DATA_5utr[rowSums(is.na(DATA_5utr)) < 1,]
DATA_mRNA <- DATA_mRNA[rowSums(is.na(DATA_mRNA)) < 1,]

# Complete cases
DATA_3utr <- DATA_3utr[complete.cases(DATA_3utr),]
DATA_5utr <- DATA_5utr[complete.cases(DATA_5utr),]
DATA_mRNA <- DATA_mRNA[complete.cases(DATA_mRNA),]

# Checkin structure
str(DATA_3utr)
str(DATA_5utr)
str(DATA_mRNA)

#DATA_3utr$SPKM.L <- as.numeric(DATA_3utr$SPKM.L
#DATA_3utr$SPKM.S <- as.numeric(DATA_3utr$SPKM.S
#DATA_5utr$SPKM.L <- as.numeric(DATA_5utr$SPKM.S

# Correlation coefficient
cor(DATA_3utr$SPKM.L, DATA_3utr$SPKM.S)
cor(DATA_5utr$SPKM.L, DATA_5utr$SPKM.S)
cor(DATA_mRNA$SPKM.L, DATA_mRNA$SPKM.S)

# Lineal model Fitting
DATA_3utr.fit <- glm(DATA_3utr$SPKM.S ~ DATA_3utr$SPKM.L)
DATA_5utr.fit <- glm(DATA_5utr$SPKM.S ~ DATA_5utr$SPKM.L)
DATA_mRNA.fit <- glm(DATA_mRNA$SPKM.S ~ DATA_mRNA$SPKM.L)

# Ploting
max_3utr <- max(max(DATA_3utr$SPKM.L), max(DATA_3utr$SPKM.S))
max_5utr <- max(max(DATA_5utr$SPKM.L), max(DATA_5utr$SPKM.S))
max_mRNA <- max(max(DATA_mRNA$SPKM.L), max(DATA_mRNA$SPKM.S))
pdf("seedMap.pdf", width = 3.5, height = 10)
par(mfrow=c(3,1))
plot(DATA_3utr$SPKM.L, DATA_3utr$SPKM.S, ylim=c(0, max_3utr), xlim=c(0, max_3utr), main = "Correlation of miRNA seed enrichment\non L & S 3'UTR", ylab="S 3'UTR", xlab="L 3'UTR")
abline(DATA_3utr.fit)
text(70, 10, labels=paste(expression(R), " = ", round(cor(DATA_3utr$SPKM.L, DATA_3utr$SPKM.S), 4), "\n", "Slope = ", round(DATA_3utr.fit$coefficients[[2]], 4),  sep=""))
mtext(side=3, adj=5, at=-1.3, line=2, "A")
plot(DATA_5utr$SPKM.L, DATA_5utr$SPKM.S, ylim=c(0, max_5utr), xlim=c(0, max_5utr), main = "Correlation of miRNA seed enrichment\non X. laevis L genes & X. tropicalis ortologs", ylab="Xtro 3'UTR", xlab="L 3'UTR")
abline(DATA_5utr.fit)
text(70, 10, labels=paste(expression(R^2), " = ",  round(cor(DATA_5utr$SPKM.L, DATA_5utr$SPKM.S), 4), "\n", "Slope = ", round(DATA_5utr.fit$coefficients[[2]], 4),  sep=""))
mtext(side=3, adj=5, at=-1.8, line=2, "B")
plot(DATA_mRNA$SPKM.L, DATA_mRNA$SPKM.S, ylim=c(0, max_mRNA), xlim=c(0, max_mRNA), main = "Correlation of miRNA seed enrichment\non X. laevis S genes & X. tropicalis ortologs", ylab="Xtro 3'UTR", xlab="S 3'UTR")
abline(DATA_mRNA.fit)
text(70, 10, labels=paste(expression(R^2), " = ", round(cor(DATA_mRNA$SPKM.L, DATA_mRNA$SPKM.S), 4), "\n", "Slope = ", round(DATA_mRNA.fit$coefficients[[2]], 4), sep=""))
mtext(side=3, adj=5, at=-1.3, line=2, "C")
dev.off()

# Plotting seed enrichment lines
pdf("Enrichment_seedmap_lines.pdf", width=8, height=8)

plot(1:length(rownames(DATA_3utr)), DATA_3utr$SPKM.L, type="l", col="blue", main="Seed Enrichment Lines Plot", ylab="Seed Counts", xlab="Genes", lwd=1)
lines(1:length(rownames(DATA_3utr)), DATA_3utr$SPKM.S, col="red")
lines(1:length(rownames(DATA_5utr)), DATA_5utr$SPKM.S, col="green")
legend("bottomright", col=c("blue", "red", "green"), c("Xlae.L", "Xlae.S", "Xtro"), lwd=2)

dev.off()
