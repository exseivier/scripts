#!/usr/bin/env Rscript
# This script reads data obtained from sylamer analysis results and plots them.
# Author: Javier Montalvo-Arredondo
# Contact: javier.montalvo@cinvestav.mx
# Version 0.2

print("Hello world!")
args = commandArgs(trailingOnly=TRUE)

max_num <- function(DATA){
	# Finds the maximum value of DATA and the corresponding miRNA.
	p <- 0
	for(i in 2:length(DATA[,1])){
		maximum <- max(as.numeric(as.vector(as.matrix(DATA[i, 2:length(DATA[1,])]))))
		if(p < maximum){
			p <- maximum
			mirna <- as.vector(as.matrix(DATA[i,1]))
			print(p)
			print(mirna)
			print(i)
		}
	}
	return(c(p, mirna))
}

min_num <- function(DATA){
	# Finds the minimum value of DATA.
	p <- 0
	for(i in 2:length(DATA[,1])){
		minimum <- min(as.numeric(as.vector(as.matrix(DATA[i, 2:length(DATA[1,])]))))
		if(p > minimum){
			p <- minimum
			print(p)
			print(i)
		}
	}
	return(p)
}


elploteo_01 <- function(infile){
	# Input file -> infile

	# Reading data.
	# DATA structure contains the x axis labels at the first row of the data
	# and after that every row contains the values of enrichment for each kmer
	# observed at the left hand side of the sequences compared to the rigth side
	# at each determined point.
	print(paste("Loading DATA...", infile, sep=" "))
	DATA <- read.delim(infile, header=F)
	#DATA <- read.delim(paste("sylByOnlyGrow/", infile, sep=""), header=TRUE, check.names=FALSE, row.names=1)

	# Getting labels for x axis.
	x <- as.vector(as.matrix(DATA[1, 2:length(DATA[1, ])]))

	# Getting data from first row.
	yp <- as.numeric(as.vector(as.matrix(DATA[2, 2:length(DATA[1,])])))

	# Searching for maximum number in DATA and the corresponding miRNA.
	print("Getting maximum number...")
	results <- max_num(DATA)
	print(paste("maxnum", results[1], sep=" "))
	print(paste("mirna", results[2], sep=" "))
	#maxnum <- 50

	# Getting the minimum number from DATA.
	# It was set -50 by default.
	print("Getting minimum number...")
	#minnum <- min_num(DATA)
	minnum <- -50

	# Plotting data from the first row of DATA.
	print("Plotting DATA...")
	pdf(paste(infile, ".pdf", sep=""))
	plot(x, yp, type="l", ylim=c(minnum, as.numeric(results[1])), main="3' UTRs", col="gray")

	# Plotting the data from row 3 to the end of DATA.
	for(i in 3:length(DATA[,1])){
		lines(x, as.numeric(as.vector(as.matrix(DATA[i, 2:length(DATA[1,])]))), col=ifelse(DATA[i, 1] == results[2], "red", "gray"))
	}

	dev.off()
}

#########################################################################################

shift_stack <- function(x, value, N) {
	# Shifts the stack of topN miRNAs seed sequences to the rigth
	# if the over-representation of that miRNA is greater than the past
	# miRNA observed.
	# Hypothetical example
	# ["ATCGTAG", "AGGCTAG", "AGCCTGA"] -> "GGTCATT"
	# "ATCGTAG" -> ["AGGCTAG", "AGCCTGA", "GGTCATT"]
	# And the first miRNA that was in the stack will be removed
	if(length(x) < N){
		stack <- c(x, value)
		return(stack)
	} else {
		stack <- c(x, value)
		stack <- stack[-1]
		return(stack)
	}
}

max_num2 <- function(DATA){
	# Finds the maximum value of DATA and the corresponding miRNA.
	p <- 0
	for(i in 1:length(DATA[,1])){
		maximum <- max(DATA[i,])
		#print(head(DATA[i, row_len]))
		if(p < maximum){
			p <- maximum
			print(p)
			mirna <- rownames(DATA)[i]
			print(mirna)
		}
	}
	return(c(p, mirna))
}

min_num2 <- function(DATA){
	# Finds the minimum value of DATA.
	p <- 0
	for(i in length(DATA[,1])){
		minimum <- min(DATA[i,])
		if(p > minimum){
			p <- minimum
			mirna <- rownames(DATA)[i]
		}
	}
	return(c(p, mirna))
}

topN <- function(DATA, N){
	# Returns an array of top N of overrepresented miRNAs. N = {1, 2, 3, ..., N}
	row_names <- rownames(DATA)
	row_len <- length(DATA[1,]) %/% 3
	maximum_nums <- c()
	for(i in 1:length(DATA[,1])){
		maximum_nums <- c(maximum_nums, max(DATA[i,1:row_len]))
	}
	maximum_nums <- as.matrix(maximum_nums)
	rownames(maximum_nums) <- row_names
	ordered_row_names <- rownames(as.matrix(maximum_nums[order(-maximum_nums[,1]),]))
	N_ordered_row_names <- ordered_row_names[1:N]
	return(N_ordered_row_names)
}


elploteo_02 <- function(infile, N){
	# Input file -> infile

	# Reading data.
	# DATA structure contains the x axis labels at the first row of the data
	# and after that every row contains the values of enrichment for each kmer
	# observed at the left hand side of the sequences compared to the rigth side
	# at each determined point.
	print(paste("Loading DATA...", infile, sep=" "))
	
	#DATA <- read.delim(infile, header=F)
	DATA <- read.delim(infile, header=TRUE, check.names=FALSE, row.names=1, comment.char="#")
	DATA <- as.matrix(DATA)
	DATA <- cbind("0"=0, DATA)
	
	METADATA <- readLines(infile, n=10)
	METADATA <- METADATA[grep("#", METADATA)]
	METADATA <- gsub("#", "", METADATA)
	print(METADATA)

	#print(head(DATA))
	
	# Getting labels for x axis.
	x_string <- colnames(DATA)
	x <- as.numeric(x_string)

	# Getting kmers names
	kmers <- rownames(DATA)

	# Getting the topN miRNAs
	print("Selecting top N miRNAs...")
	topNmirnas <- topN(DATA, N)
	print(topNmirnas)
	
	# Getting data from first row.
	yp <- as.numeric(as.vector(DATA[1,]))

	# Sellecting colors
	#total_colors <- colors(distinct=T)
	#set.seed(15888)
	#my_colors <- sample(total_colors, N)
	my_colors <- c("red", "blue", "violet", "green", "black")
	#colour_counter <- ifelse(kmers[1] %in% topNmirnas, 2, 1)
	col1 <- ifelse(kmers[1] %in% topNmirnas, my_colors[match(kmers[i], topNmirnas)], "gray")

	# Searching for maximum number in DATA and the corresponding miRNA.
	print("Getting maximum number...")
	max_results <- max_num2(DATA)
	print(paste("maxnum", max_results[1], sep=" "))
	print(paste("mirna", max_results[2], sep=" "))
	#maxnum <- 50

	# Getting the minimum number from DATA.
	# It was set -50 by default.
	print("Getting minimum number...")
	min_results <- min_num2(DATA)
	print(paste("minnum", min_results[1], sep=" "))
	print(paste("mirna", min_results[2], sep=" "))
	#minnum <- -50

	# Plotting data from the first row of DATA.
	print("Plotting DATA...")
	plot(x, yp, type="l", ylim=c(as.numeric(min_results[1]) - 5, as.numeric(max_results[1]) + 20), ylab="P(X)", xlab=METADATA[2], main=METADATA[1], col=col1, cex.main=1)
	mtext(LETTERS[j], side=3, adj=-0.05, line=1.3, cex=1.5)

	# Plotting the data from row 2 to the end of DATA.
	for(i in 2:length(DATA[,1])){
		lines(x, as.vector(DATA[i, ]), col=ifelse(kmers[i] %in% topNmirnas, my_colors[match(kmers[i], topNmirnas)], "gray"))
		if(kmers[i] %in% topNmirnas){
			print(paste("ploting", kmers[i], "with color", my_colors[match(kmers[i], topNmirnas)], sep=" "))
		}
	#	colour_counter <- ifelse(kmers[i] %in% topNmirnas, colour_counter+1, colour_counter)
	}
	legend("topright", topNmirnas, col=my_colors[1:N], lty=1)
}

path <- args[1]
N <- args[2]
setwd(paste(path, "sylResults", sep="/"))
getwd()
lsdir <- list.files()
for(i in 1:length(lsdir)){
	print("Processing...")
	setwd(paste(path, lsdir[i], sep="/"))
	print(getwd())
	lsf <- list.files(pattern=".syl$")
	lsf_len <- length(lsf)
	rows <- (lsf_len %/% 2) + (lsf_len %% 2)
	print(paste("Rows: ", rows, sep = ""))
	columns <- ifelse(lsf_len > 1, 2, 1)
	print(paste("Columns: ", columns, sep=""))
	pdf(paste(lsdir[i], ".pdf", sep=""), width=4*columns, height=4*rows)
	par(mfrow=c(rows, columns))
	par(oma=c(0,0,2,0))
	for(j in 1:length(lsf)){
		elploteo_02(lsf[j], N)
	}
	title(main=lsdir[i], outer=T, cex=0.5)
	dev.off()
	setwd("../")
}

.silenced <- function(){
path <- args[1]
setwd(paste(path, "sylResults", sep="/"))
getwd()
lsdir <- list.files()
lsdir
lsdir_len <- length(lsdir)
for(i in 1:length(lsdir)){
	print("Processing...")
	setwd(paste(path, lsdir[i], sep="/"))
	print(getwd())
	lsf <- list.files(pattern=".syl$")
	for(j in 1:length(lsf)){
		elploteo_01(lsf[j])
	}
	dev.off()
	setwd("../")
}
}



