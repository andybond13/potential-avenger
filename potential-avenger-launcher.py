#!/usr/bin/python
import os

#clear results
command = "rm results/vtkFiles/*.vtk"
os.system(command)

#define input arguments
ts_refine = 4
end_t = 2.5e-5*100
Nelt = 14#20
lc = 0.072
startWithLoad = 0 
strain_rate = 1 #100.0
printVTK = 1
oneAtATime = 1
minOpenDist = 0.05
alpha = 0.5		#hardening coefficient bounded by 0 & 1: 0-very brittle, 1-very plastic
localOnly = 1	#1 for local only, 0 for local/non-local
visualizeCracks = 1 #1 to visualize cracks (elements disappear when d = 1)
fullCompression = 1 #0: s = E*e*(1-d) always, 1: s=E*e*(1-d) in tension, s=E*e in compression

#run program
command = "./potential-avenger.exe %f %f %f %u %f %u %u %u %f %f %u %u %u" % (strain_rate, ts_refine, end_t, Nelt, lc, startWithLoad, printVTK, oneAtATime, minOpenDist, alpha, localOnly, visualizeCracks,fullCompression)
print command
os.system(command)
