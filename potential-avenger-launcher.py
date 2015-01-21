#!/usr/bin/python
import os

#clear results
command = "rm results/vtkFiles/*.vtk"
os.system(command)

#define input arguments
ts_refine = 4
end_t = 8.0e-5#2.5e-5*100
list = [1000]
Gc = 5.902e-5
Yc = 8.1967e-4
lc = 0.75 * Gc/Yc#0.072
startWithLoad = 1 
strain_rate = 0.25 #1 #100.0
printVTK = 1
oneAtATime = 0
minOpenDist = 0.00
#(SQRT) hardening coefficient bounded by 0 & 1: 0-very brittle, 1-very plastic
#(LIN) lambda = 2*sigmac*lc/(E*wc); wc = 2 * Gc/sigmac; sigmac^2 = 2*E*Yc; lambda = 2*Yc*lc/Gc
alpha = 0.95 	#2*Yc*lc/Gc
alpha = 2*Yc*lc/Gc
sm = "LIN"			#SQRT - match sqrt cohesive TSL or LIN - match linear cohesive TSL
localOnly = 0	#1 for local only, 0 for local/non-local
visualizeCracks = 0 #1 to visualize cracks (elements disappear when d = 1)
fullCompression = 1 #0: s = E*e*(1-d) always, 1: s=E*e*(1-d) in tension, s=E*e in compression
elemDeath = 1 #0: element death not on, 1: element death on

#for Nelt in list:
#run program
command = prefix+"./potential-avenger.exe %f %f %f %u %f %u %u %u %f %f %u %u %u %s %u" % (strain_rate, ts_refine, end_t, Nelt, lc, startWithLoad, printVTK, oneAtATime, minOpenDist, alpha, localOnly, visualizeCracks,fullCompression,sm,elemDeath)+suffix
print command
os.system(command)
