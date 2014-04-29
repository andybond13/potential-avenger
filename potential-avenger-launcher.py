#!/usr/bin/python
import os

#clear results
command = "rm *.vtk"
os.system(command)

#define input arguments
ts_refine = 4
end_t = 0.7
Nelt = 100
lc = 0.09
intOrder = 1
strain_rate = 1.0
printVTK = 1
oneAtATime = 1
minOpenDist = 0.05
alpha = 0.5		#hardening coefficient bounded by 0 & 1: 0-very brittle, 1-very plastic
localOnly = 1	#1 for local only, 0 for local/non-local

#run program
command = "./potential-avenger.exe %f %f %f %u %f %u %u %u %f %f %u" % (strain_rate, ts_refine, end_t, Nelt, lc, intOrder, printVTK, oneAtATime, minOpenDist, alpha, localOnly)
print command
os.system(command)
