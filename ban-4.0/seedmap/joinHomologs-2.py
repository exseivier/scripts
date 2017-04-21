#!/usr/bin/env python

"""
Sorts the data obtained from seedMap.py.
It only requires the output file from seedMap.py
"""

from sys import argv

filename = argv[1]
delimiter = argv[2]
print delimiter
DATA = {}
IN = open(filename, "r")
OUT = open("%s.sorted" % filename, "w+")
for line in IN:
	line = line.strip()
	sline = line.split(delimiter)[0].split(".")
	print sline
	if sline[-1] == "L" or sline[-1] == "S":
		shortID = "".join(sline[:-1])
	else:
		shortID = "".join(sline)
		print shortID
	try:
		DATA[shortID].append(line)
	except Exception:
		DATA[shortID] = [line]

for key, value in DATA.iteritems():
	if len(value) == 2:
		value = sorted(value)
		OUT.write("%s\t%s\t%s\n" % (key, value[0], value[1]))

IN.close()
OUT.close()
