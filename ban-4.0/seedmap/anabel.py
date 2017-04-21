#!/usr/bin/env python

from sys import argv, exit
from random import randrange

def anagram(seed):
	"""(STR) -> STR
	Creates an anagram of a k-mer
	"""
	anaSeed = []
	indexes = []
	i = 0
	while i < len(seed):
		randnum = randrange(0, len(seed))
		if randnum in indexes:
			pass
		else:
			indexes.append(randnum)
			i += 1
	for j in xrange(len(seed)):
		anaSeed.append(seed[indexes.pop()])
	return "".join(anaSeed)
# Reads a fasta file of kmers and for each kmer creates an anagram for everyone
# and write the results in fasta file
anagrams = {}
IN = open(argv[1], "r")
OUT = open("ana_%s" % argv[1], "w+")
for line in IN:
	line = line.strip()
	if line[0] == ">":
		header = line
		anagrams[header] = ""
	else:
		anagrams[header] = anagram(line)

for key, value in anagrams.iteritems():
	OUT.write("%s\n%s\n" % (key, value))

IN.close()
OUT.close()
