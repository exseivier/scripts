#!/usr/bin/env Rscript

### Analysing seddMap results of 3' and 5' UTRs; and mRNA and Promoter sequences

# Load table
DATA_3utr <- read.delim("Prom_Xla_total_ori.smp.sorted", header=F)
DATA_5utr <- read.delim("Prom_Xla_total_ana1.smp.sorted", header=F)
DATA_mRNA <- read.delim("Prom_Xla_total_ana2.smp.sorted", header=F)
DATA_mRNA2 <- read.delim("Prom_Xla_total_ana3.smp.sorted", header=F)

# row names settings
rownames(DATA_3utr) <- DATA_3utr[,1]
rownames(DATA_5utr) <- DATA_5utr[,1]
rownames(DATA_mRNA) <- DATA_mRNA[,1]
rownames(DATA_mRNA2) <- DATA_mRNA2[,1]

# Deleting first column
DATA_3utr <- DATA_3utr[,2:9]
DATA_5utr <- DATA_5utr[,2:9]
DATA_mRNA <- DATA_mRNA[,2:9]
DATA_mRNA2 <- DATA_mRNA2[,2:9]

# Colnames settings
cnames <- c("ID.L", "len_kb.L", "counts.L", "SPKM.L", "ID.S", "len_kb.S", "counts.S", "SPKM.S")
colnames(DATA_3utr) <- cnames
colnames(DATA_5utr) <- cnames
colnames(DATA_mRNA) <- cnames
colnames(DATA_mRNA2) <- cnames

# Deleting NA rows
DATA_3utr <- DATA_3utr[rowSums(is.na(DATA_3utr)) < 1,]
DATA_5utr <- DATA_5utr[rowSums(is.na(DATA_5utr)) < 1,]
DATA_mRNA <- DATA_mRNA[rowSums(is.na(DATA_mRNA)) < 1,]
DATA_mRNA2 <- DATA_mRNA2[rowSums(is.na(DATA_mRNA2)) < 1,]

# Complete cases
DATA_3utr <- DATA_3utr[complete.cases(DATA_3utr),]
DATA_5utr <- DATA_5utr[complete.cases(DATA_5utr),]
DATA_mRNA <- DATA_mRNA[complete.cases(DATA_mRNA),]
DATA_mRNA2 <- DATA_mRNA2[complete.cases(DATA_mRNA2),]

# Checkin structure
str(DATA_3utr)
str(DATA_5utr)
str(DATA_mRNA)
str(DATA_mRNA2)

# Correlation coefficient
cor(DATA_3utr$SPKM.L, DATA_3utr$SPKM.S)
cor(DATA_5utr$SPKM.L, DATA_5utr$SPKM.S)
cor(DATA_mRNA$SPKM.L, DATA_mRNA$SPKM.S)
cor(DATA_mRNA2$SPKM.L, DATA_mRNA2$SPKM.S)

# Lineal model Fitting
DATA_3utr.fit <- glm(DATA_3utr$SPKM.S ~ DATA_3utr$SPKM.L)
DATA_5utr.fit <- glm(DATA_5utr$SPKM.S ~ DATA_5utr$SPKM.L)
DATA_mRNA.fit <- glm(DATA_mRNA$SPKM.S ~ DATA_mRNA$SPKM.L)
DATA_mRNA2.fit <- glm(DATA_mRNA2$SPKM.S ~ DATA_mRNA2$SPKM.L)

# Ploting
max_3utr <- max(max(DATA_3utr$SPKM.L), max(DATA_3utr$SPKM.S))
max_5utr <- max(max(DATA_5utr$SPKM.L), max(DATA_5utr$SPKM.S))
max_mRNA <- max(max(DATA_mRNA$SPKM.L), max(DATA_mRNA$SPKM.S))
max_mRNA2 <- max(max(DATA_mRNA2$SPKM.L), max(DATA_mRNA2$SPKM.S))

pdf("seedMap.pdf", width = 3.5, height = 10)
par(mfrow=c(3,1))
plot(DATA_3utr$SPKM.L, DATA_3utr$SPKM.S, ylim=c(0, max_3utr), xlim=c(0, max_3utr), main = "Correlation of miRNA seed enrichment\non L & S 3'UTR", ylab="S 3'UTR", xlab="L 3'UTR")
abline(DATA_3utr.fit)
text(800, 800, labels=paste(expression(R^2), round(cor(DATA_3utr$SPKM.L, DATA_3utr$SPKM.S), 4), sep=" = "))
mtext(side=3, adj=5, at=-1.3, line=2, "A")
plot(DATA_5utr$SPKM.L, DATA_5utr$SPKM.S, ylim=c(0, max_5utr), xlim=c(0, max_5utr), main = "Correlation of miRNA seed enrichment\non X. laevis L genes & X. tropicalis ortologs", ylab="Xtro 3'UTR", xlab="L 3'UTR")
abline(DATA_5utr.fit)
text(800, 800, labels=paste(expression(R^2), round(cor(DATA_5utr$SPKM.L, DATA_5utr$SPKM.S), 4), sep=" = "))
mtext(side=3, adj=5, at=-1.8, line=2, "B")
plot(DATA_mRNA$SPKM.L, DATA_mRNA$SPKM.S, ylim=c(0, max_mRNA), xlim=c(0, max_mRNA), main = "Correlation of miRNA seed enrichment\non X. laevis S genes & X. tropicalis ortologs", ylab="Xtro 3'UTR", xlab="S 3'UTR")
abline(DATA_mRNA.fit)
text(500, 500, labels=paste(expression(R^2), round(cor(DATA_mRNA$SPKM.L, DATA_mRNA$SPKM.S), 4), sep=" = "))
mtext(side=3, adj=5, at=-1.3, line=2, "C")
plot(DATA_mRNA2$SPKM.L, DATA_mRNA2$SPKM.S, ylim=c(0, max_mRNA2), xlim=c(0, max_mRNA2), main = "Correlation of miRNA seed enrichment\non X. laevis S genes & X. tropicalis ortologs", ylab="Xtro 3'UTR", xlab="S 3'UTR")
abline(DATA_mRNA2.fit)
text(500, 500, labels=paste(expression(R^2), round(cor(DATA_mRNA2$SPKM.L, DATA_mRNA2$SPKM.S), 4), sep=" = "))
mtext(side=3, adj=5, at=-1.3, line=2, "C")

dev.off()
