#!/usr/bin/env python

# This script chages the header of a fasta sequence entry
# based on a database called parentDB. This database must
# contain almost 2 columns, the entry of the first column
# must be equal to the current header of the sequence, and
# the entry of the second column must be the header replacement.
# The info of the second column was extracted from the gff file
# of the transcriptome using bash tools such as grep, cut...

from sys import argv, exit
from Utils import PDStruct, replace_header

def main():
	# Input fasta file.
	inputfile = argv[1]
	# Output fasta file.
	outputfile = argv[2]
	# parentDB database.
	parentDB = argv[3]
	# HASH data structure constructed with data from parentDB.
	pds = PDStruct(parentDB)
	# For each fasta sequence, the current header is replaced by
	# the entry of the second column in parentDB.
	if (replace_header(inputfile, outputfile, pds)):
		print "Success!"
	else:
		print "Something went wrong!"
		exit(0)

if __name__ == "__main__":
	main()
