#!/usr/bin/python
import os

#define input arguments
ts_refine = 2
end_t = 2
Nelt = 20
lc = 0.02
intOrder = 2
strain_rate = 2.5

#run program
command = "./potential-avenger.exe %f %f %f %u %f %u" % (strain_rate, ts_refine, end_t, Nelt, lc, intOrder)
print command
os.system(command)
