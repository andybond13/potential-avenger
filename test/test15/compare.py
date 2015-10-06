#/usr/bin/python

import numpy

def main(output,gold):
	status = 1

	#read results & compare
	data = numpy.loadtxt(output,skiprows=1)
	data_gold = numpy.loadtxt(gold,skiprows=1)
	assert(len(data) == len(data_gold))

	for i in range(0,len(data)):
		if (data[i] != data_gold[i]):
			return 0

	return 1

