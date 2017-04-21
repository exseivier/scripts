#!/usr/bin/env python

from Utils import geneModel, load_Ome, rev_comp
from sys import argv, exit
import pickle

message = """
Usage:
	eProm.py <genome_file> <gff_file> <output_file>
"""

genomefile = argv[1]
if genomefile == "--help" or genomefile == "-h":
	print message
	exit(0)
gff_file = argv[2]
out_file = argv[3]
OUT = open(out_file, "w+")

try:
	IN = open("%s.dat" % genomefile, "rb")
	genome = pickle.load(IN)
	IN.close()
	print "Try success"
except Exception:
	genome = load_Ome(genomefile)
	print "except Exception sucesss"

gene_models = geneModel(gff_file)


for key, value in gene_models.iteritems():
	if value["ori"] == "+":
		if value["start"]-1500 < 0:
			this_start = 0
		else:
			this_start = value["start"]-1500
		sequence = genome[">%s" % value["locus"]][this_start:value["start"]+3]
	elif value["ori"] == "-":
		sequence = rev_comp(genome[">%s" % value["locus"]][value["end"]-3:value["end"]+1500])
	
	OUT.write(">%s|%s\n%s\n" % (value["name"], key.split(".")[0], sequence))


OUT.close()
