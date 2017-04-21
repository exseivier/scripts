#!/usr/bin/env Rscript

print("\n\n\nR analysis begins\n********************\n")

args <- commandArgs(trailingOnly=T)
library("Biostrings")
input1 <- args[1]
input2 <- args[2]

input1.tr <- readDNAStringSet(input1, format="fasta", use.names=T)
#input1.tr <- sample(input1.tr)
input1.tr <- input1.tr[order(names(input1.tr)),]
input2.tr <- readDNAStringSet(input2, format="fasta", use.names=T)
#input2.tr <- sample(input2.tr)
input2.tr <- input2.tr[order(names(input2.tr)),]
count_TRUE <- 0
count_FALSE <- 0
input1.tr.names <- sort(names(input1.tr))
head(input1.tr.names, n=20)
input2.tr.names <- sort(names(input2.tr))
head(input2.tr.names, n=20)

#//
#	Checking if input1.tr.names and input2.tr.names have the same quantity of elements
#//
if (length(input1.tr.names) == length(input2.tr.names)) {
	print("Equal number of elements")
} else {
	print("Different number of elements")
}

#//
#	Selecting unique elements of input1.tr and input2.tr
#//
input1.unique <- c()
pivot <- (1/200) * length(input1.tr.names)
for (i in seq(length(input1.tr.names))) {
	found <- "FALSE"
	for (j in seq(length(input2.tr.names))){
		if(i > pivot) {
			print(paste("Input1: ", i, ". Input2: ", j, sep=""))
			pivot = pivot + ((1/200) * length(input1.tr.names))
		}
		if (input1.tr.names[i] == input2.tr.names[j]) {
			input2.tr.names <- input2.tr.names[j * -1]
			found <- "TRUE"
			break
		}
	}
	# HERE AFTER BREAK...
	if (found == "FALSE") {
		input1.unique <- c(input1.unique, input1.tr.names[i])
	} else {
		#	I think that nothing happens!
	}
}
writeLines(input1.unique, "input1.unique.txt", sep="\n")
print(paste("Unique elements number of input1 is: ", length(input1.unique)))
writeLines(input2.tr.names, "input2.unique.txt", sep="\n")
print(paste("Unique elements number of input2 is: ", length(input2.tr.names)))

#//
#	Checking sequences one by one.
#//
true <- c()
false <- c()
for (i in seq(length(input1.tr))) {
	count <- 0
	for (j in seq(length(input2.tr))) {
		#print(j)
		if (input1.tr[i] == input2.tr[j]) {
			count <- count + 1
			input2.tr <- input2.tr[j * -1]
			break
		}
		if (j > 500) {
			break
		}
	}
	# HERE AFTER BREAK...
	if (count == 0) {
		count_FALSE <- count_FALSE + 1
		print(paste("Checking ", names(input1.tr[i]), " FALSE", sep=""))
		false <- c(false, paste("Checking ", names(input1.tr[i]), " FALSE", sep=""))
	}
	else {
		count_TRUE <- count_TRUE + 1
		print(paste("Checking ", names(input1.tr[i]), " TRUE", sep=""))
		true <- c(true, paste("Checking ", names(input1.tr[i]), " TRUE", sep=""))
	}
	
}
writeLines(false, "false.txt", sep="\n")
writeLines(true, "true.txt", sep="\n")

