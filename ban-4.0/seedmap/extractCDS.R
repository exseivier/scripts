#!/usr/bin/env Rscript

print("Extracting CDS...")
library(rtracklayer)
library(Biostrings)

# Paths
# Javo modified it, getting input path by command line in bash
args <- commandArgs(trailingOnly=TRUE)
genome.Path <- args[1]
outPath <- "CDS_OUT"

system(paste("mkdir", outPath, sep=" "))

transcriptome.Dir <- dirname(args[2])
transcriptome.File <- basename(args[2])
genome.Dir <- dirname(args[1])
genome.File <- basename(args[1])

# min length utr seq
minLen <- 100

# Import gff
gff <- import.gff3(args[2])

# Get CDS gff
CDS <- gff[grepl("three_prime_UTR", gff$type)]
CDS$type <- sub("^trans.+:", "", as.character(CDS$Parent))
CDS <- CDS[!is.na(CDS$type)]
CDS_gff <- paste(outPath, "/CDSs.gff", sep = "")
export.gff(CDS, con = paste(outPath, "/CDSs.gff", sep = ""))

# Get CDS exon sequences
system(paste("bedtools getfasta -s -name -fi",  paste(genome.Dir, genome.File, sep="/"), "-bed", CDS_gff, "-fo", paste(outPath, "CDS_exons.fa", sep="/"), sep = " "))
fa <- readDNAStringSet(paste(outPath, "CDS_exons.fa", sep="/"))

# Cei modified, since now bedtools outputs positions added to names
names(fa) <- sub("::.+","",names(fa))

# Concat 3'utrs chunks
faTab <- aggregate(as.character(fa), by = list(names(fa)), function(x) paste(x, collapse = ""))
faVector <- (faTab$x)
names(faVector) <- faTab$Group.1
faCDS <- DNAStringSet(faVector, use.names = TRUE)

# Cei modified to put the minLen filter here, instead of at the exon level
faCDS <- faCDS[width(faCDS) >= minLen]

# Save 3'utrs in fasta format
writeXStringSet(faCDS, filepath = (paste(outPath, paste("CDS", genome.File, sep="_"), sep = "/")))

