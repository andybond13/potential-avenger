#!/usr/bin/python
import os

#clear results
command = "rm *.vtk"
os.system(command)

#define input arguments
ts_refine = 2
end_t = 2
Nelt = 20
lc = 0.02
intOrder = 2
strain_rate = 2.5
printVTK = 0

#run program
command = "./potential-avenger.exe %f %f %f %u %f %u %u" % (strain_rate, ts_refine, end_t, Nelt, lc, intOrder, printVTK)
print command
os.system(command)
