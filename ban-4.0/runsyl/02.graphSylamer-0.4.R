#!/usr/bin/env Rscript
# This script reads data obtained from sylamer analysis results and plots them.
# Author: Javier Montalvo-Arredondo
# Contact: javier.montalvo@cinvestav.mx
# Version 0.2

args = commandArgs(trailingOnly=T)

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
	print("I am inside min_num2")
	for(i in length(DATA[,1])){
		minimum <- min(DATA[i,])
		if(p > minimum){
			p <- minimum
			print(p)
			mirna <- rownames(DATA)[i]
			print(mirna)
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


elploteo_02 <- function(infile, N, miRNAs.table){
	# Input file -> infile

	# Reading data.
	# DATA structure contains the x axis labels at the first row of the data
	# and after that every row contains the values of enrichment for each kmer
	# observed at the left hand side of the sequences compared to the rigth side
	# at each determined point.

	#Loading DATASET
	print(paste("\nLoading DATA...", infile, sep=" "))
	DATA <- read.delim(infile, header=TRUE, check.names=FALSE, row.names=1, comment.char="#")
	DATA <- as.matrix(DATA)
	DATA <- cbind("0"=0, DATA)

	#Loading METADATA
	print("Loading Metadata...")
	METADATA <- readLines(infile)
	METADATA <- METADATA[grep("#", METADATA)]
	METADATA <- gsub("#", "", METADATA)
	print(METADATA)
	
	#Getting selected miRNAs
	selected.mirna <- c()
	if (length(METADATA) >= 3){
		for (a in 3:length(METADATA)) {
			entry <- strsplit(METADATA[a], "\t")
			entry <- cbind(entry[[1]][1], entry[[1]][2], entry[[1]][3])
			selected.mirna <- rbind(selected.mirna, entry)
		}
		selected.mirna <- as.data.frame(selected.mirna)
		colnames(selected.mirna) <- c("sequence", "name", "color")
		selected.mirna$sequence <- as.character(selected.mirna$sequence)
		selected.mirna$name <- as.character(selected.mirna$name)
		selected.mirna$color <- as.character(selected.mirna$color)
		print ("Selected miRNAs...")
		str(selected.mirna)
	}
	else {
	print("No selected miRNAs!")
	}

	print("Filtering data...")
	DATA <- DATA[rownames(DATA) %in% miRNAs.table$sequence, ]
	
	# Getting labels for x axis.
	x_string <- colnames(DATA)
	x <- as.numeric(x_string)

	# Getting kmers names
	kmers <- rownames(DATA)
		
	# Getting data from first row.
	yp <- as.numeric(as.vector(DATA[1,]))
	
	#Plotting DATA
	print("Plotting DATA...")
	if (length(selected.mirna[,1]) > 0){
		plot(x, yp, type="l", ylim=c(-5, 20), ylab="P(X)", xlab=METADATA[2], main=METADATA[1], col=ifelse(rownames(DATA)[1] %in% selected.mirna$sequence, selected.mirna[selected.mirna$sequence == rownames(DATA)[1], "color"], "gray") , cex.main=1)
		mtext(LETTERS[j], side=3, adj=-0.05, line=1.3, cex=1.5)

		# Plotting the data from row 2 to the end of DATA.
		for(i in 2:length(DATA[,1])){
			lines(x, as.vector(DATA[i, ]), col=ifelse(rownames(DATA)[i] %in% selected.mirna$sequence, selected.mirna[selected.mirna$sequence == rownames(DATA)[i], "color"], "gray"))
		}
		legend("topright", legend=selected.mirna$name, col=selected.mirna$color, lty=1)
	}
	else {
		plot(x, yp, type="l", ylim=c(-5, 20), ylab="P(X)", xlab=METADATA[2], main=METADATA[1], col = "gray" , cex.main=1)
		mtext(LETTERS[j], side=3, adj=-0.05, line=1.3, cex=1.5)

		# Plotting the data from row 2 to the end of DATA.
		for(i in 2:length(DATA[,1])){
			lines(x, as.vector(DATA[i, ]), col = "gray")
		}
	}
}

path <- args[1]
setwd(paste(path, "sylResults", sep="/"))
getwd()
lsdir <- list.files()
for(i in 1:length(lsdir)){
	print("Processing...")
	
	#Setting working directory
	setwd(paste(path, lsdir[i], sep="/"))
	print(getwd())
	
	#Getting the file name of *.syl files list
	lsf <- list.files(pattern=".syl$")
	lsf_len <- length(lsf)
	
	#Getting the miRNA table file name
	miRNAs.file <- list.files(pattern=".dat$")
	miRNAs.table <- read.table(miRNAs.file, header=T, row.names=1)
	print("miRNAs table uploaded!")
	miRNAs.table$color <- rainbow(length(miRNAs.table[,1]))
	
	#Setting canvas to plot data
	rows <- (lsf_len %/% 2) + (lsf_len %% 2)
	print(paste("Rows: ", rows, sep = ""))
	columns <- ifelse(lsf_len > 1, 2, 1)
	print(paste("Columns: ", columns, sep=""))
	pdf(paste(lsdir[i], ".pdf", sep=""), width=4*columns, height=4*rows)
	par(mfrow=c(rows, columns))
	par(oma=c(0,0,2,0))
	
	#Plotting data
	for(j in 1:length(lsf)){
		elploteo_02(lsf[j], N, miRNAs.table)
	}
	
	title(main=lsdir[i], outer=T, cex=0.5)
	dev.off()
	setwd("../")
}

