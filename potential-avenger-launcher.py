#!/usr/bin/python
import os,time,sys

#clear results
command = "rm results/vtkFiles/*.vtk"
os.system(command)

prefix = ""
suffix = ""

#define openMP & Cilk theads
numThreads = 1
##os.environ["OMP_NUM_THREADS"] = str(numThreads)
##os.system("echo $OMP_NUM_THREADS")
os.environ["CILK_NWORKERS"] = str(numThreads)
os.system("echo $CILK_NWORKERS")

#define input arguments
ts_refine = 1
#end_t = 8.0e-5#2.5e-5*100
Gc = 8.313e-5
Yc = 8.1967e-4
startWithLoad = 1
SR = [0.025,0.075,0.25,0.75,2.5,7.5,25,75,250,750,2500,7500,25000,75000,250000,750000]
LC = [5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,5.071e-2,2.483e-2,1.889e-2,1.182e-2,7.077e-3,4.43e-3,2.651e-3,1.660e-3]
strain_rate = 7500
#lcmult = float(sys.argv[1])
#lc = 0.01
lc = LC[13] * 0.5
#lc = 3./4.*Gc/Yc
print "*** lc = ",lc
end_t = 4e-2 * 0.25 / strain_rate * 10.0 #* 0.005 #50 #2.5 #0.01 for SR=0.25 
print "end_t=",end_t
#Nelt = max(1000, 20*strain_rate)
#Nelt = 1/lc * 8 
#Nelt = 10000
Nelt = int(sys.argv[1])
#Nelt = int(8/lc) 
#print sys.argv[1] 
CV = 0.0


#(SQRT) hardening coefficient bounded by 0 & 1: 0-very brittle, 1-very plastic
#(LIN) lambda = 2*sigmac*lc/(E*wc); wc = 2 * Gc/sigmac; sigmac^2 = 2*E*Yc; lambda = 2*Yc*lc/Gc
alpha = 2*Yc*lc/Gc
print alpha
assert(alpha <= 0.5)
#model options
sm = "LIN"         	#SQRT - match sqrt cohesive TSL or LIN - match linear cohesive TSL
oneAtATime = 0		#0 or 1. If 1, cracks open 1 per timestep. If 0, all can open at once
minOpenDist = 0		#minimum distance between two non-local insertions
TLSoption = 1   	#0, 1, or 2. 0 for local only. 1 for local/non-local. 2 for non-local-only.		note: for non-local-only, be careful with nucleation amount (e.g. 1.0*h) vs. lc... if h > lc, then you're putting in d=1. (avoid h>lc anyway)
visualizeCracks = 0 #1 to visualize cracks (elements disappear when d = 1)
fullCompression = 0 #0 or 1. If 0, s = E*e*(1-d) always, if 1: s=E*e*(1-d) in tension, s=E*e in compression
elemDeath = 0 		#0 or 1. If 0, element death not on, if 1: element death on
frontExtension = 2	#0, 1, or 2. If 0, explicit delay-damage front extension. If 1, implicit-type NR solve. If 2, explicit NR-style solve (use with maxIteration = 1 for example) 
maxIteration = 5	#>= 1. If 1, single iteration.

#output option
printVTK = 0 #100*(Nelt/1000)
#for Nelt in list:
print lc
#run program
#os.system("echo $OMP_NUM_THREADS")
command = prefix+"./potential-avenger.exe %f %f %f %u %f %f %f %u %u %s %u %f %u %u %u %u %u %u" % (strain_rate, ts_refine, end_t, Nelt, lc, alpha, CV, startWithLoad, printVTK, sm, oneAtATime, minOpenDist, TLSoption, visualizeCracks, fullCompression, elemDeath, frontExtension, maxIteration)+suffix
print command
os.system(command)
