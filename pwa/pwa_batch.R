# Combining 

library(Biostrings)

args <- commandArgs(trailingOnly=T)

# Making combinations "number of files" choose combinations of two of them
# If you have 3 files; it takes them and makes the posible combinations as follow
# {1,2}; {1,3}; {2,3}
# and then iterates over those combination and makes pairwise alignment with the sequences
# of the files which positions in the passed arguments are represented in those combinations.

print("Making combinations...")
first.gen <- args[length(args)-3]
second.gen <- args[length(args)-2]
method <- args[length(args)-1]
output <- args[length(args)]
combinations <- c()
for ( i in 1:(length(args) - 5) ) {
	a = i+1
	for ( j in a:(length(args) - 4) ) {
		combinations <- rbind(combinations, c(i,j))
	}
}
combinations

# Loading sequences...
print("Loading sequences...")
sequences <- c()

for ( i in 1:length(combinations[,1]) ) {
	DNAobject <- readDNAStringSet(args[i], format="fasta")
	sequences <- c(sequences, DNAobject)
}

# Checking sequences loading
print("Loaded sequences are...")
for ( i in 1:length(sequences) ) {
	print(paste("Sequences from file ", args[i], sep=""))
	#print(names(sequences[[i]])[as.integer(first.gen):as.integer(second.gen)])
}

names <- names(sequences[[1]])[as.integer(first.gen):as.integer(second.gen)]

# Performing PWA for every combination
alignments <- c()
column.names <- c()
for ( i in 1:length(combinations[,1]) ) {
	picked.seq1 <- combinations[i,1]
	picked.seq2 <- combinations[i,2]
	print(paste("Performing BATCH PWA over ", args[picked.seq1], " and ", args[picked.seq2], sep=""))
	alignments <- c(alignments, pairwiseAlignment(sequences[[picked.seq1]][as.integer(first.gen):as.integer(second.gen)], sequences[[picked.seq2]][as.integer(first.gen):as.integer(second.gen)], type=method))
	column.names <- c(column.names, paste(args[picked.seq1], " Vs ", args[picked.seq2], sep=""))
}

pids <- c()
for ( i in 1:length(alignments) ) {
	pids <- cbind(pids, pid(alignments[[i]]))
}

pids
pids <- as.data.frame(pids)
pids
rownames(pids) <- names
colnames(pids) <- column.names
pids

write.table(pids, output, sep="\t")



