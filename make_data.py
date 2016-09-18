import os
import sys
import numpy as np

cmdarg = False
try:
    cmdarg = sys.argv[1]
except:
    None

if cmdarg == "new_time":
    #delete time-datafile and start a new
    first_line = "#data for CPU time in c++ program"
    time_filename = "data/dderiv_time_c++.dat"
    #os.system("rm " + time_filename)
    outfile = open(time_filename, "w")
    outfile.write(first_line + "\n")
    outfile.close()

#make time data for LU-decomposition
#run LU-decomposition for log10(n) = 1,2,3
n_range_LU = np.logspace(1,3,num=3)
for n in n_range_LU:
    os.system("./../build-project1-Desktop_Qt_5_7_0_GCC_64bit-Release/project1 %d 1"%n)
#run tridiag for log10(n) = 1,2,3,4,5,6
n_range_tridiag = np.logspace(1,6,num=6)
for n in n_range_tridiag:
    os.system("./../build-project1-Desktop_Qt_5_7_0_GCC_64bit-Release/project1 %d "%int(n))

#make u(x) data for tridiagonal methods
#run tridiag methods for n = logspace(10,10000,some number of variable)
n_range_tridiag = np.logspace(1, 4, num=4)
for n in n_range_tridiag:
    os.system("./../build-project1-Desktop_Qt_5_7_0_GCC_64bit-Release/project1 %d "%int(n))



