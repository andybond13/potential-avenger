#!/usr/bin/python
import os,time

#clear results
command = "rm results/vtkFiles/*.vtk"
os.system(command)

prefix = ""
suffix = ""

#define input arguments
ts_refine = 1
Gc = 8.313e-5
Yc = 8.1967e-4
startWithLoad = 1
strain_rate = 1
lc = 0.02 
end_t = 4e-2 * 0.25 / strain_rate /100 /2#* 10 
print "end_t=",end_t
Nelt = 1/lc * 8 
printVTK = 0
oneAtATime = 0
minOpenDist = 0
alpha = 2*Yc*lc/Gc
print alpha
assert(alpha < 0.5)
sm = "LIN"         #SQRT - match sqrt cohesive TSL or LIN - match linear cohesive TSL
localOnly = 0   #1 for local only, 0 for local/non-local
visualizeCracks = 0 #1 to visualize cracks (elements disappear when d = 1)
fullCompression = 0 #0: s = E*e*(1-d) always, 1: s=E*e*(1-d) in tension, s=E*e in compression
elemDeath = 1 #0: element death not on, 1: element death on

#for Nelt in list:
print lc
#run program
command = prefix+"./potential-avenger.exe %f %f %f %u %f %u %u %u %f %f %u %u %u %s %u" % (strain_rate, ts_refine, end_t, Nelt, lc, startWithLoad, printVTK, oneAtATime, minOpenDist, alpha, localOnly, visualizeCracks,fullCompression,sm,elemDeath)+suffix
print command
os.system(command)
