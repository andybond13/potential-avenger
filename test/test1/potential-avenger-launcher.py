#!/usr/bin/python
import os

#clear results
command = "rm results/vtkFiles/*.vtk"
os.system(command)

#define input arguments
ts_refine = 4
end_t = 8.0e-5#2.5e-5*100
Nelt = 100#20
Gc = 5.902e-5
Yc = 8.1967e-4
lc = 0.75 * Gc/Yc#0.072
startWithLoad = 1 
strain_rate = 0.25 #1 #100.0
printVTK = 0
oneAtATime = 0
minOpenDist = 0.00
alpha = 0.25	#(SQRT) hardening coefficient bounded by 0 & 1: 0-very brittle, 1-very plastic; (LIN) 
CV = 0.0

#model options
sm = "SQRT"          #SQRT - match sqrt cohesive TSL or LIN - match linear cohesive TSL
oneAtATime = 0      #0 or 1. If 1, cracks open 1 per timestep. If 0, all can open at once
minOpenDist = 0     #minimum distance between two non-local insertions
TLSoption = 1       #0, 1, or 2. 0 for local only. 1 for local/non-local. 2 for non-local-only.
visualizeCracks = 0 #1 to visualize cracks (elements disappear when d = 1)
fullCompression = 0 #0 or 1. If 0, s = E*e*(1-d) always, if 1: s=E*e*(1-d) in tension, s=E*e in compression
elemDeath = 0       #0 or 1. If 0, element death not on, if 1: element death on
frontExtension = 1  #0, 1, or 2. If 0, explicit delay-damage front extension. If 1, implicit-type NR solve. If 2, explicit NR-style solve (use with maxIteration = 1 for example) 
maxIteration = 20   #>= 1. If 1, single iteration.

#run program
command = "./potential-avenger.exe %f %f %f %u %f %f %f %u %u %s %u %f %u %u %u %u %u %u" % (strain_rate, ts_refine, end_t, Nelt, lc, alpha, CV, startWithLoad, printVTK, sm, oneAtATime, minOpenDist, TLSoption, visualizeCracks, fullCompression, elemDeath, frontExtension, maxIteration)

print command
os.system(command)
