#!/usr/bin/env Rscript

print("Hello world!")
library(rtracklayer)
library(Biostrings)

#//
#	Function definition section
#//

la_reversa_con_platano <- function(x) {
	rc_strs <- c()
	for (item in x) {
		item <- DNAString(item)
		item <- reverseComplement(item)
		item <- as.character(item)
		rc_strs <- c(rc_strs, item)
	}
	return(rc_strs)
}


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
gff <- gff[grepl("exon", gff$type)]
gff$type <- sub("^trans.+:", "", as.character(gff$Parent))
gff <- gff[!is.na(gff$type)]
gff_plus <- gff[strand(gff) == "+"]
gff_minus <- gff[strand(gff) == "-"]
gff_plus <- gff_plus[order(as.character(gff_plus$type))]
gff_minus <- gff_minus[order(as.character(gff_minus$type))]
mRNA_plus_gff <- paste(outPath, "/mRNA_plus.gff", sep = "")
mRNA_minus_gff <- paste(outPath, "/mRNA_minus.gff", sep="")
export.gff(gff_plus, con = mRNA_plus_gff)
export.gff(gff_minus, con = mRNA_minus_gff)

system(paste("bedtools getfasta -name -fi",  paste(genome.Dir, genome.File, sep="/"), "-bed", mRNA_plus_gff, "-fo", paste(outPath, "mRNA_plus_exons.fa", sep="/"), sep = " "))
fa_plus <- readDNAStringSet(paste(outPath, "mRNA_plus_exons.fa", sep="/"))
system(paste("bedtools getfasta -name -fi",  paste(genome.Dir, genome.File, sep="/"), "-bed", mRNA_minus_gff, "-fo", paste(outPath, "mRNA_minus_exons.fa", sep="/"), sep = " "))
fa_minus <-  readDNAStringSet(paste(outPath, "mRNA_minus_exons.fa", sep="/"))

names(fa_plus) <- sub("::.+","",names(fa_plus))
names(fa_minus) <- sub("::.+","",names(fa_minus))

faTab_plus <- aggregate(as.character(fa_plus), by = list(names(fa_plus)), function(x) paste(x, collapse = ""))
faVec_plus <- (faTab_plus$x)
names(faVec_plus) <- faTab_plus$Group.1
famRNA_plus <- DNAStringSet(faVec_plus, use.names = TRUE)

faTab_minus <- aggregate(as.character(fa_minus), by = list(names(fa_minus)), function(x) paste(x, collapse = ""))
faVec_minus <- (la_reversa_con_platano(faTab_minus$x))
names(faVec_minus) <- faTab_minus$Group.1
famRNA_minus <- DNAStringSet(faVec_minus, use.names = TRUE)

writeXStringSet(famRNA_plus, filepath = (paste(outPath, paste("mRNA", genome.File, sep="_"), sep = "/")))
writeXStringSet(famRNA_minus, filepath = (paste(outPath, paste("mRNA", genome.File, sep="_"), sep = "/")), append=T)

.silenced <- function() {
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

}

