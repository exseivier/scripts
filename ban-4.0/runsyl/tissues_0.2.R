#!/usr/bin/env Rscript

#Setting working directory
#setwd("/home/montalvo/scripts/ban-4.0/Rscripts")

args <- commandArgs(trailingOnly=T)

### Loading edgeR & locfit libraries
library(edgeR)
library(locfit)

### Loading tissues gene expression data
tissues <- read.delim(args[1], header=T, row.names="ID")
head(tissues)

### Building DGE Object
groups.tiss <- factor(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
design.tiss <- model.matrix(~0+groups.tiss)
treat.tiss <- c("brain", "eye", "heart", "intestine", "kidney", "liver", "lung", "muscle", "ovary", "pancreas", "skin", "spleen", "stomach", "testis")
rownames(design.tiss) <- colnames(tissues)
colnames(design.tiss) <- treat.tiss
design.tiss
tissues.DGE <- DGEList(counts=tissues, group=groups.tiss)
tissues.DGE

### Filtering data
keep <- rowSums(cpm(tissues) > 1) >= 2
tissues <- tissues[keep,]

### Calculate Normalized Factors and estimate disperssion
tissues.DGE <- calcNormFactors(tissues.DGE)
tissues.DGE <- estimateDisp(tissues.DGE, design.tiss)

### Fit data to model Quasi-likellihood F distribution
tissues.fit <- glmQLFit(tissues.DGE, design.tiss)


selected.tissues <- eval(parse(text=args[2]))
#selected.tissues
selected.tissues.len <- length(selected.tissues)
#selected.tissues.len
indexes <- match(selected.tissues, treat.tiss)
#indexes
remaining.tissues <- treat.tiss[-indexes]
#remaining.tissues
remaining.tissues.len <- length(remaining.tissues)
#paste(selected.tissues, collapse="+")
#class(eval(parse(text=paste(selected.tissues, collapse="+"))))
my.ctst <- makeContrasts(cts=((eval(parse(text=paste(selected.tissues, collapse="+")))/selected.tissues.len)-(eval(parse(text=paste(remaining.tissues, collapse="+")))/remaining.tissues.len)), levels=design.tiss)
my.ctst

tissues.qlf.cts <- glmQLFTest(tissues.fit, contrast=my.ctst[,"cts"])

tissues.qlf.cts$table$PValue_log <- abs(log10(tissues.qlf.cts$table$PValue))
tissues.qlf.cts$table$PValue_log_sig <- ifelse(tissues.qlf.cts$table$logFC < 0, tissues.qlf.cts$table$PValue_log*-1, tissues.qlf.cts$table$PValue_log)
tissues.qlf.cts$table <- tissues.qlf.cts$table[order(tissues.qlf.cts$table$logFC),]
tissues.qlf.cts$OOtable <- tissues.qlf.cts$table[order(tissues.qlf.cts$table$PValue_log_sig),]

write.table(tissues.qlf.cts$OOtable, file=paste("tissues.", paste(selected.tissues, collapse="_"), ".OO.txt", sep=""))

.silenced <- function(){
### Making contrasts
my.ctst.tiss <- makeContrasts(Brain=(brain-((eye+heart+intestine+kidney+liver+lung+muscle+ovary+pancreas+skin+spleen+stomach+testis)/13)), levels=design.tiss)

### QLF Testing
tissues.qlf.brain <- glmQLFTest(tissues.fit, contrast=my.ctst.tiss[,"Brain"])
head(tissues.qlf.brain$table)

### Sorting tissues.qlf.brain$table
tissues.qlf.brain$Otable <- tissues.qlf.brain$table[tissues.qlf.brain$table$PValue <= 1e-1,]
tissues.qlf.brain$Otable <- tissues.qlf.brain$Otable[order(tissues.qlf.brain$Otable[,1]),]
dim(tissues.qlf.brain$Otable[tissues.qlf.brain$Otable$logFC <= 0,])

## Sorting tissues.qlf.brain$table by the value of log10(tissues.qlf.brain$table$PValue)
## with the same sign of the corresponding logFC
tissues.qlf.brain$table$PValue_log <- abs(log10(tiussues.qlf.brain$table$PValue))
tissues.qlf.brain$table$PValue_log_sig <- ifelse(tissues.qlf.brain$table$logFC < 0, tissues.qlf.brain$table$PValue_log*-1, tissues.qlf.brain$table$PValue_log)
tissues.qlf.brain$table <- tissues.qlf.brain$table[order(tissues.qlf.brain$table$logFC),]
tissues.qlf.brain$OOtable <- tissues.qlf.brain$table[order(tissues.qlf.brain$table$PValue_log_sig),]

##writing to file
write.table(tissues.qlf.brain$table, file="tissues.updown.txt")
write.table(tissues.qlf.brain$Otable, file="tissues.updown.O.txt")
write.table(tissues.qlf.brain$OOtable, file="tissues.updown.OO.txt")

# Saving workspace
save(list = ls(all.names = TRUE), file = "tissues.RData", envir = .GlobalEnv)
}
