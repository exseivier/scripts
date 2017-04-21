#! /usr/bin/env python

# This script measures the coefficient of the reciprocal 'in cis' elements
# difference between two paralogs.
from sys import argv, exit

def summing(array):
	"""
	Returns the sum of the all values from the array
	"""
	summed = 0.0
	for item in array:
		summed += float(item)
	return summed

def percentage(array, summed):
	"""
	Returns an array of values. each value is the percentage, taking into account 
	the sum of all values as one hundred percent (100%)
	"""
	arrayOfPercentages = []
	coeff = 0.0
	for item in array:
		if summed > 0:
			coeff = float(item)/float(summed)
			coeff *= 100.0
		else:
			coeff = 0
		arrayOfPercentages.append(coeff)
	return arrayOfPercentages

def absoluteValue(value):
	"""
	Returns the absolute value
	"""
	if value < 0:
		return value * -1
	else:
		return value

def reciprocalDifferenceCoefficient(arrays):
	"""
	Returns the reciprocal difference coefficient
	"""
	sum_a0 = summing(arrays[0])
	sum_a1 = summing(arrays[1])
	if sum_a0 > 0 and sum_a1 > 0:
		summed = 0.0
		for i in xrange(len(arrays[0])):
			summed += absoluteValue(arrays[0][i] - arrays[1][i])
		return float(summed)/200.0
	elif sum_a0 <= 0 and sum_a1 > 0:
		return 1.0
	elif sum_a0 > 0 and sum_a1 <= 0:
		return 1.0
	elif sum_a0 <= 0 and sum_a1 <= 0:
		return 0.0
	else:
		print "WTF!!!"


laDataStructure = {}
# filename of data file
filename = argv[1]
output = "%s.rdc" % filename
# Open file. Declaring a file handler
DATA = open(filename, "r")
# Generating the basic structure
for line in DATA:
	line.strip()
	sliced_line = line.split(">")
	# Adding the raw data for L subgenome
	laDataStructure[sliced_line[1]] = {sliced_line[2].split("\t")[0] : [sliced_line[2].split("\t")[1:-1]]}
	# Appending the sum of the values fo L subgenome
	laDataStructure[sliced_line[1]][sliced_line[2].split("\t")[0]].append(summing(sliced_line[2].split("\t")[1:-1]))
	# Calculating percentages of the raw values for L subgenome
	laDataStructure[sliced_line[1]][sliced_line[2].split("\t")[0]].append(percentage(sliced_line[2].split("\t")[1:-1], laDataStructure[sliced_line[1]][sliced_line[2].split("\t")[0]][1]))
	# Adding the raw data for S subgenome
	laDataStructure[sliced_line[1]].update({sliced_line[3].split("\t")[0] : [sliced_line[3].split("\t")[1:]]})
	# Appending the sum of the values of S subgenome
	laDataStructure[sliced_line[1]][sliced_line[3].split("\t")[0]].append(summing(sliced_line[3].split("\t")[1:]))
	# Calculating percentages of the raw values for S subgenome
	laDataStructure[sliced_line[1]][sliced_line[3].split("\t")[0]].append(percentage(sliced_line[3].split("\t")[1:], laDataStructure[sliced_line[1]][sliced_line[3].split("\t")[0]][1]))

""" laDataStructure was checked... [OK]
count = 0
for key, value in laDataStructure.iteritems():
	print key
	print value
	count += 1
	if count == 5:
		break
"""
# Calculating the reciprocal difference coefficient
OUT = open(output, "w+")
for key, values in laDataStructure.iteritems():
	arrays = []
	for IDsg, array in values.iteritems():
		arrays.append(array[2])
	rdc = reciprocalDifferenceCoefficient(arrays)
	OUT.write("%s\t%f\n" %(key, rdc))

OUT.close()


DATA.close()
