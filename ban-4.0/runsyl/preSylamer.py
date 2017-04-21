#!/usr/bin/env python
from sys import argv
from Utils import formatseq, headsArray

"""
preSylamer Usage: preSylamer.py <UTR-fasta-file> <edgeR-result-table>

This script compares the header of each UTR sequence with the ID or
gene symbol from the edgeR analysis and select only the entries of
the edgeR result table which ID is present in the header of the UTR
sequences.

"""


# Check if seq is formated
try:
	fh = open("%s.out" % argv[1], "r")
	fh.close()
except IOError:
	formatseq(argv[1], None)
	print "Formating sequences file..."

# Create the headers array
heads = headsArray("%s.out" % argv[1])

# Select metched entries
IN = open(argv[2], "r")
OUT = open("%s.filt" % argv[2], "w+")
for entry in IN:
	entry = entry.strip()
	if entry.split("\t")[0] in heads:
		OUT.write("%s\n" % entry)

IN.close()
OUT.close()
