#!/usr/bin/python
import os,time

#clear results
command = "rm results/vtkFiles/*.vtk"
os.system(command)

prefix = ""
suffix = ""

#define input arguments
ts_refine = 1
#end_t = 8.0e-5#2.5e-5*100
Gc = 8.313e-5
Yc = 8.1967e-4
startWithLoad = 1
SR = [0.025,0.075,0.25,0.75,2.5,7.5,25,75,250,750,2500,7500,25000,75000,250000,750000]
LC = [5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,2.483e-2,1.889e-2,1.182e-2,7.077e-3,4.43e-3,2.651e-3,1.660e-3]
strain_rate = 0.000100
lc = 0.02
#lc = 3./4.*Gc/Yc
print "*** lc = ",lc
end_t = 1.0e-4 #4e-2 * 0.25 / strain_rate  
print "end_t=",end_t
#Nelt = max(1000, 20*strain_rate)
#Nelt = 1/lc * 8 
Nelt = 1000
printVTK = 0
oneAtATime = 0
minOpenDist = 0 
#(SQRT) hardening coefficient bounded by 0 & 1: 0-very brittle, 1-very plastic
#(LIN) lambda = 2*sigmac*lc/(E*wc); wc = 2 * Gc/sigmac; sigmac^2 = 2*E*Yc; lambda = 2*Yc*lc/Gc
##alpha = 1- 0.2*1*Yc/(8*DE_CZM[i])   #2*Yc*lc/Gc
##alpha = max(0,1- NF_CZM[i]*Yc*1*0.2*lc*4/3/(DE_CZM[i]*1) )  #2*Yc*lc/Gc
alpha = 2*Yc*lc/Gc
print alpha
assert(alpha <= 0.5)
#alpha = 0.0;
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
