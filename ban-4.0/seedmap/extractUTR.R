#!/usr/bin/env Rscript

print("Hello world!")
library(rtracklayer)
library(Biostrings)

# Paths
# Javo modified it, getting input path by command line in bash
args <- commandArgs(trailingOnly=TRUE)
genome.Path <- args[1]
outPath <- "3_UTR_OUT"

system(paste("mkdir", outPath, sep=" "))

transcriptome.Dir <- dirname(args[2])
transcriptome.File <- basename(args[2])
genome.Dir <- dirname(args[1])
genome.File <- basename(args[1])

# min length utr seq
minLen <- 50

# Import gff
gff <- import.gff3(args[2])

# Get 3'utrs gff
threeUTR <- gff[grepl("three_prime_UTR", gff$type)]
threeUTR$type <- sub("^trans.+:", "", as.character(threeUTR$Parent))
threeUTR <- threeUTR[!is.na(threeUTR$type)]
three_UTR_gff <- paste(outPath, "/three_UTR.gff", sep = "")
export.gff(threeUTR, con = paste(outPath, "/three_UTR.gff", sep = ""))

# Get 3'utrs exon sequences
system(paste("bedtools getfasta -s -name -fi",  paste(genome.Dir, genome.File, sep="/"), "-bed", three_UTR_gff, "-fo", paste(outPath, "three_UTR_exons.fa", sep="/"), sep = " "))
fa <- readDNAStringSet(paste(outPath, "three_UTR_exons.fa", sep="/"))

# Cei modified, since now bedtools outputs positions added to names
names(fa) <- sub("::.+","",names(fa))

# Concat 3'utrs chunks
faTab <- aggregate(as.character(fa), by = list(names(fa)), function(x) paste(x, collapse = ""))
faVector <- (faTab$x)
names(faVector) <- faTab$Group.1
faThreeUTR <- DNAStringSet(faVector, use.names = TRUE)

# Cei modified to put the minLen filter here, instead of at the exon level
faThreeUTR <- faThreeUTR[width(faThreeUTR) >= minLen]

# Save 3'utrs in fasta format
writeXStringSet(faThreeUTR, filepath = (paste(outPath, paste("three_UTR", genome.File, sep="_"), sep = "/")))



