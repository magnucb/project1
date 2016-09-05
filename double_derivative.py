import numpy as np
import astropy as ap
import matplotlib.pylab as pl
import sys
import os

"""
Solve (d/dx)^2 u = f for known f.
Double differentials using Euler forward and backward methods:
->(d/dx)^2 u = (u[i+1] - 2 u[i] + u[i-1])/h^2 = f
-> u[i+1] - 2 u[i] + u[i-1] = h^2*f = y[i]
-> generalize -> a[i]u[i-1] + b[i]u[i] + c[i]u[i+1] = y[i]
"""

curdir = os.getcwd()
try:
    version = int(sys.argv[2])
except IndexError:
    version = "0"
except ValueError:
    sys.exit("Commandline-argument must be integer")

def write2file(outstring,
               filename=curdir+"/data/testfile_v%d.dat",
               append=True):
    """
    If 'append' is True:
    -open 'filename'.
    -append 'outstring' to end of file
    -close file
    If 'append' is False:
    -open a new file 'filename' (deleting the old one)
    -write 'outstring' to file
    -close file
    """
    outstring = str(outstring)
    if append:
        with open(filename,'a') as outfile:
            outfile.write(outstring)
    else:
        with open(filename,'w') as outfile:
            outfile.write(outstring)
    return outstring

def forwardalg():
    None
def backwardalg():
    None

n=10; x0=0.0; x1=1.0; h=(x1-x0)/(n-1.0)
#vectors of tridiagonal matrix A
a = np.ones(n)
b = -2*np.ones(n)
c = np.ones(n)
#vectors u, y and x
u = np.zeros(n)
x = np.linspace(x0,x1,n) 
y = h*h*100*np.exp(-10*x)

