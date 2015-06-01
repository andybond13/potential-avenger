#!/usr/bin/python
from os import listdir
from os.path import isfile, join
import contextlib
from contextlib import nested
import math

def transform1(block):
	L = block[-1][0]
	#print L
	r = L/(2.0* math.pi)
	outblock = []
	xout = []
	yout = []
	xorig = []

	for line in block:
		x = line[0]
		xorig.append(x)
		line[0] = r*math.cos(x/r)
		xout.append(line[0])
		line[1] = r*math.sin(x/r)
		yout.append(line[1])
		line[2] = 0.0
		str1 = ' '.join(str(e) for e in line)+' \n'
		outblock.append(str1)
	return outblock,L,xorig,xout,yout

def transform2(block,L,x0,x1,y1):
	Uend = block[-1][0]
	#print Uend
	rp = (L+Uend)/(2.0* math.pi)
	outblock = []

	index = 0
	for line in block:
		#print line
		u = line[0]
		line[0] = rp*math.cos((x0[index]+u)/rp) - x1[index]
		line[1] = rp*math.sin((x0[index]+u)/rp) - y1[index]
		line[2] = 0.0
		str1 = ' '.join(str(e) for e in line)+' \n'
		outblock.append(str1)
		index += 1
	return outblock

def main(infiles):
	for i in infiles:
		if i[-4:] != ".vtk":
			continue
	
		if "circle" in i:
			continue

		output = i[:-4]+"_circle.vtk"
		print i,"->",output
		status = ""
		count = 0
		lineno = 0
#		outblock = [][]

		#print i
		with contextlib.nested(open(i),open(output,"w")) as (fIn, fOut):
			for line in fIn:
				#print line.strip()
				if (status == "points"):
					lineno += 1
					lss = line.strip().split()
					if (len(lss) >= 3):
						#print lss,len(lss)
						outblock[lineno-1][0] = float(lss[0])
						outblock[lineno-1][1] = float(lss[1])
						outblock[lineno-1][2] = float(lss[2])

				if (status == "displacement"):
					lineno += 1
					lss = line.strip().split()
					if (len(lss) >= 3):
						#print lss,len(lss)
						outblock[lineno-1][0] = float(lss[0])
						outblock[lineno-1][1] = float(lss[1])
						outblock[lineno-1][2] = float(lss[2])

				if (status == ""):
					fOut.write(line)

				if (len(line.strip().split()) >= 3):
					if (line.strip().split()[0] == "POINTS"):
						status = "points"
						count = int(line.strip().split()[1])
						outblock = [[0 for x in range(3)] for x in range(count)] 	
						#print i,count

					if (line.strip().split()[1] == "displacements"):
						status = "displacement"
						outblock = [[0 for x in range(3)] for x in range(count)] 	

				if (lineno == count and count > 0):
#					print "end points",lineno,count,status
					lineno = 0
					#transform outblock
					if (status == "points"):
						outblock,L,x0,x1,y1 = transform1(outblock)
						assert(len(x1) == count)
					if (status == "displacement"):
						assert(len(x1) == count)
						outblock = transform2(outblock,L,x0,x1,y1)
					status = ""
					for outline in outblock:
						fOut.write(outline)



path = "results/vtkFiles/"
list = [path+f for f in listdir(path) if isfile(join(path,f)) ]
main(list)
