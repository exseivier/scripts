#!/usr/bin/env python
from sys import argv, path
path.append('/home/montalvo/scripts/ban-4.0/utils/')
import pickle
from Utils import load_Ome, count_matches, rev_comp, rev_trans
transcriptome = {} # Declaring empty hash
miRNAome = {} # Declaring empty hash
mRNAcounts = {} # Declaring empty hash
total_out = argv[3]
bymiRNA_out = argv[4]
try:
	transF = open("%s.dat" % argv[1], "rb")
	transcriptome = pickle.load(transF)
	transF.close()
except IOError:
	transcriptome = load_Ome(argv[1])
try:
	miRNAF = open("%s.dat" % argv[2], "rb")
	miRNAome = pickle.load(miRNAF)
	miRNAF.close()
except IOError:
	miRNAome = load_Ome(argv[2])

header = ""
for miHeader in miRNAome.iterkeys():
	header += "%s\t" % miHeader.split(" ")[0][1:]


whole_counts = 0
for transHeader, sequence in transcriptome.iteritems():
	total_counts = 0
	mRNAcounts[transHeader] = {}
	mRNAcounts[transHeader]["mRNA_len"] = float(len(sequence)) / 1000.0
	for miHeader, seed in miRNAome.iteritems():
		mRNAcounts[transHeader][miHeader] = {}
		mRNAcounts[transHeader][miHeader]["counts"] = count_matches(rev_trans(seed[2:8]), sequence)
		total_counts += mRNAcounts[transHeader][miHeader]["counts"]
	mRNAcounts[transHeader]["total_counts"] = total_counts
	whole_counts += total_counts

F1M = float(whole_counts)/1000000.0
for transHeader in mRNAcounts.iterkeys():
	mRNAcounts[transHeader]["SPKM_by_mRNA"] = \
	float(mRNAcounts[transHeader]["total_counts"]) / (mRNAcounts[transHeader]["mRNA_len"] * F1M)
	for miHeader in miRNAome.iterkeys():
		mRNAcounts[transHeader][miHeader]["SPKM_by_miRNA"] = \
		float(mRNAcounts[transHeader][miHeader]["counts"]) / (mRNAcounts[transHeader]["mRNA_len"] * F1M)


mRNAOUT = open(total_out, "w+")
miRNAOUT = open(bymiRNA_out, "w+")


mRNAOUT.write("Key\tmRNA_len\ttotal_counts\tSPKM_by_mRNA\n")
miRNAOUT.write("GeneName\t%s\n" % header)

for key, values in mRNAcounts.iteritems():
	print key
	mRNAOUT.write("%s\t%f\t%f\t%f\n" % (key, values["mRNA_len"], values["total_counts"], values["SPKM_by_mRNA"]))
	print values["mRNA_len"]
	print values["total_counts"]
	print values["SPKM_by_mRNA"]
	str_out = ""
	for miHeader in miRNAome.iterkeys():
		str_out += "%s\t" %mRNAcounts[key][miHeader]["SPKM_by_miRNA"]
	miRNAOUT.write("%s\t%s\n" % (key, str_out))

mRNAOUT.write("\n\nF1M: %f\tWhole counts: %f\n" % (F1M, whole_counts))
print ""
print F1M
print whole_counts
mRNAOUT.close()
miRNAOUT.close()
