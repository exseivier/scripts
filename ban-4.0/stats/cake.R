### Exploratory analysis

# Data loading
genes.total <- 45099
genes.L <- 16835
genes.S <- 12697
genes.O <- 15567
genes.total.var <- 1617
genes.L.var <- 878
genes.S.var <- 598
genes.O.var <- 141
genes.total.NoVar <- genes.total - genes.total.var
genes.L.NoVar <- genes.L - genes.L.var
genes.S.NoVar <- genes.S - genes.S.var
genes.O.NoVar <- genes.O - genes.O.var
percent.genes.L.var <- (genes.L.var / genes.L) * 100
percent.genes.S.var <- (genes.S.var / genes.S) * 100
percent.genes.O.var <- (genes.O.var / genes.O) * 100
percent.genes.L.NoVar <- (genes.L.NoVar / genes.L) * 100
percent.genes.S.NoVar <- (genes.S.NoVar / genes.S) * 100
percent.genes.O.Novar <- (genes.O.NoVar / genes.O) * 100


# Parameter settings
pdf("cakes.pdf")
par(mfrow=c(2,2))
par(mar=c(4, 5, 4, 4))

# LSO-Plotting
slices <- c(genes.L, genes.S, genes.O)
lbls <- c("L genes", "S genes", "Orphans")
pct <- round(slices/sum(slices) * 100)
lbls <- paste(lbls, pct)
lbls <- paste(lbls, "%", sep="")
pie(slices, labels=lbls, col=rainbow(length(lbls)), main="Percentage of Orphans\nand L or S genes")
mtext(side=3, adj=1, at=-1.3, line=2, "A")

# Variant transcripts plotting
slices <- c(genes.total.NoVar, genes.total.var)
lbls <- c("\t1\nTranscript", "\t+1\nTranscript")
pct <- round(slices/sum(slices) * 100)
lbls <- paste(lbls, pct)
lbls <- paste(lbls, "%", sep="")
pie(slices, labels=lbls, col=rainbow(length(lbls)), main="Percentage of genes\nwith transcript variants")
mtext(side=3, adj=1, at=-1.8, line=2, "B")

# LSO-variant transcripts plotting
slices <- c(genes.L.var, genes.S.var, genes.O.var)
lbls <- c("L genes", "S genes", "Orphans")
pct <- round(slices/sum(slices) * 100)
lbls <- paste(lbls, pct)
lbls <- paste(lbls, "%", sep="")
pie(slices, labels=lbls, col=rainbow(length(lbls)), main="Percentage of Orphans,\nL or S genes\nwith trasncript variants")
mtext(side=3, adj=1, at=-1.3, line=2, "C")

# Normalized LSO-variant plotting
#par(mar=c(5,5,5,5))
NoVar <- c(percent.genes.L.NoVar, percent.genes.S.NoVar, percent.genes.O.Novar)
var <- c(percent.genes.L.var, percent.genes.S.var, percent.genes.O.var)
header <- c("L genes", "S genes", "O genes")
data <- rbind(NoVar, var)
colnames(data) <- header
x <- barplot(data, col=rainbow(length(data[,1])), main="Normalized percentage of Orphans,\n L or S genes with transcript variants", ylim=c(0,110), ylab = "%")
text(x, y=70, label=c(paste(round(percent.genes.L.NoVar, 2), "%", sep=" "), paste(round(percent.genes.S.NoVar, 2), "%", sep=" "), paste(round(percent.genes.O.Novar, 2), "%", sep=" ")), pos=3, cex=0.8, col="blue")
text(x, y=98, label=c(paste(round(percent.genes.L.var, 2), "%", sep=" "), paste(round(percent.genes.S.var, 2), "%", sep=" "), paste(round(percent.genes.O.var, 2), "%", sep=" ")), pos=3, cex=0.8, col="black")
mtext(side=3, adj=1, at=-1, line=2, "D")

dev.off()