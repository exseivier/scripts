#!/usr/bin/env python
"""
Concatenates the mapped reads counts of all tratments obtained from analysis
with programms like HTSeq-count, Kallisto, etc for every gene name
Requieres a filename list with the counts or abundance text files, and the
OUTPUT path
Example:
	./joinCounts.py <filelist.txt> <outputfile.txt> <column number>

The outputfile will be used in posterior analysis with software like edgeR

"""


from sys import argv, exit

LISTFILE = argv[1]
OUTPUT = argv[2]
column = int(argv[3])
IN = open(LISTFILE, "r")
header = ""
data_structure = {}
for line in IN:
	line = line.strip()
	header += "\t%s" % line.split("/")[-2]
	PROCESSING_FILE = open(line, "r")
	waste = PROCESSING_FILE.readline().strip()
	for entry in PROCESSING_FILE:
		entry = entry.strip()
		fields = entry.split("\t")
		if len(fields) < 2:
			continue
		try:
			data_structure[fields[0]].append(round(float(fields[column-1])))
		except Exception:
			data_structure[fields[0]] = [round(float(fields[column-1]))]
	PROCESSING_FILE.close()
OUT = open(OUTPUT, "w+")
OUT.write("ID%s\n" % header)
for key, value in data_structure.iteritems():
	OUT.write("%s" % key)
	for item in value:
		OUT.write("\t%.4f" % item)
	OUT.write("\n")
OUT.close()
