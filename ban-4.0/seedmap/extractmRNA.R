#!/usr/bin/env Rscript

print("Hello world!")
library(rtracklayer)
library(Biostrings)

# Paths
# Javo modified it, getting input path by command line in bash
args <- commandArgs(trailingOnly=TRUE)
genome.Path <- args[1]
outPath <- "mRNA_OUT"

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
mRNA <- gff[grepl("exon", gff$type)]
mRNA$type <- sub("^trans.+:", "", as.character(mRNA$Parent))
mRNA <- mRNA[!is.na(mRNA$type)]
mRNA_gff <- paste(outPath, "/mRNA.gff", sep = "")
export.gff(mRNA, con = paste(outPath, "/mRNA.gff", sep = ""))

# Get 3'utrs exon sequences
system(paste("bedtools getfasta -s -name -fi",  paste(genome.Dir, genome.File, sep="/"), "-bed", mRNA_gff, "-fo", paste(outPath, "mRNA_exons.fa", sep="/"), sep = " "))
fa <- readDNAStringSet(paste(outPath, "mRNA_exons.fa", sep="/"))

# Cei modified, since now bedtools outputs positions added to names
names(fa) <- sub("::.+","",names(fa))

# Concat 3'utrs chunks
faTab <- aggregate(as.character(fa), by = list(names(fa)), function(x) paste(x, collapse = ""))
faVector <- (faTab$x)
names(faVector) <- faTab$Group.1
famRNA <- DNAStringSet(faVector, use.names = TRUE)

# Cei modified to put the minLen filter here, instead of at the exon level
famRNA <- famRNA[width(famRNA) >= minLen]

# Save 3'utrs in fasta format
writeXStringSet(famRNA, filepath = (paste(outPath, paste("mRNA", genome.File, sep="_"), sep = "/")))



