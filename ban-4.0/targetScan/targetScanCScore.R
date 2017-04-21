#!/usr/bin/env Rscript

files <- list.files(".", pattern="^frog.")
files
headers <- c()
ws <- c()
for(i in 1:length(files)){
	print(paste("Processing ", files[i], " ...", sep=""))
	table <- read.table(files[i], header=T, row.names="Transcript")
	table.L <- table[grep(".L", rownames(table)),]
	table.S <- table[grep(".S", rownames(table)),]
	#print(head(table.L, n=25))
	#print(head(table.S, n=25))
	score.L <- table.L[,1]
	score.S <- table.S[,1]
	score.L <- score.L * -1
	score.S <- score.S * -1
	w <- wilcox.test(score.L, score.S)$p.value
	headers <- c(headers, files[i])
	ws <- c(ws, w)
	#####################SILENCED##########################################
	.slenced <- function(){
	L.breaks <- seq(min(score.L), max(score.L), by=0.01)
	S.breaks <- seq(min(score.S), max(score.S), by=0.01)
	L.cut <- cut(score.L, L.breaks)
	S.cut <- cut(score.S, S.breaks)
	L.freq <- table(L.cut)
	S.freq <- table(S.cut)
	L.cumfreq <- c(0, cumsum(L.freq/length(score.L)))
	S.cumfreq <- c(0, cumsum(S.freq/length(score.S)))
	plot(L.breaks, L.cumfreq, type="l", col="blue", main="L Vs S genes", ylab="Frequency", xlab="Context Score")
	lines(S.breaks, S.cumfreq, col="red")
	}
	######################################################################
}

wilcox.table <- cbind(headers, ws)
wilcox.table
write.table(wilcox.table[as.numeric(wilcox.table[,2]) <= 0.05,], file="wilcox.output.txt")

finalists <- wilcox.table[as.numeric(wilcox.table[,2]) <= 0.05, 1]
finalists
columns <- as.integer(length(finalists)/4)
rows <- as.integer(length(finalists)/columns) + 1
pdf("CScore_cumfreq.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for(item in finalists){
	print(paste("Plotting ", item, " ...", sep=""))
	table <- read.table(item, header=T, row.names="Transcript")
	table.L <- table[grep(".L", rownames(table)),]
	table.S <- table[grep(".S", rownames(table)),]
	#print(head(table.L, n=25))
	#print(head(table.S, n=25))
	score.L <- table.L[,1]
	score.S <- table.S[,1]
	score.L <- score.L * -1
	score.S <- score.S * -1
	L.breaks <- seq(min(score.L), max(score.L), by=0.01)
	S.breaks <- seq(min(score.S), max(score.S), by=0.01)
	L.cut <- cut(score.L, L.breaks)
	S.cut <- cut(score.S, S.breaks)
	L.freq <- table(L.cut)
	S.freq <- table(S.cut)
	L.cumfreq <- c(0, cumsum(L.freq/length(score.L)))
	S.cumfreq <- c(0, cumsum(S.freq/length(score.S)))
	w <- wilcox.test(score.L, score.S)$p.value
	plot(L.breaks, L.cumfreq, type="l", col="blue", main=paste(item, "\nL Vs S genes", sep=""), ylab="Frequency", xlab="Context Score")
	lines(S.breaks, S.cumfreq, col="red")
	text(max(max(score.L), max(score.S))-0.3,0.6, paste("W", format(round(w, 4), nsmall=4), sep="="))
	legend("bottomright", c("L genes", "S genes"), col=c("blue", "red"), lty=1)
}

dev.off()
