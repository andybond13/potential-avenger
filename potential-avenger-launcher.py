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

#run program
command = "./potential-avenger.exe %f %f %f %u %f %u %u %u" % (strain_rate, ts_refine, end_t, Nelt, lc, intOrder, printVTK, oneAtATime)
print command
os.system(command)
