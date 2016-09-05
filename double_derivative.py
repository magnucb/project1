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
#NB program defaults version to 0, if succesful run as diff. version.
curdir = os.getcwd()
try:
    n = int(sys.argv[1])
except IndexError:
    n = 10
except ValueError:
    sys.exit("Commandline-argument must be integer \nFirst number of gridpoints, then version number")

try:
    version = int(sys.argv[2])
except IndexError:
    version = "0"
except ValueError:
    sys.exit("Commandline-argument must be integer \nFirst number of gridpoints, then version number")
    
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

def forwardalg_vec(tri_bottom, tri_mid, tri_top, vert):
    for i in range(1,len(vert)):
        k = tri_bottom[i]/float(tri_mid[i-1])
        tri_mid[i] -= k*tri_top[i-1]
        vert[i] -= k*vert[i-1]
    return tri_mid, vert

def forwardalg_i(ai, bi, yi, bim, cim, yim):
    k = ai/float(bim) # 1 flop
    bi_new = bi - k*cim # 2 flops
    yi_new = yi - k*yim # 2 flops
    return bi_new, yi_new

def backwardalg_vec(tri_mid, tri_top, vert):
    i = len(vert) -2
    while i > 0:
        k = tri_top[i]/float(tri_mid[i+1])
        vert[i] -= k*vert[i+1]
        i -= 1
    return vert

def backwardalg_i(ci, yi, bip, yip):
    k = ci/float(bip) # 1 flop
    yi_new = yi - k*yip # 2 flops
    return yi_new

x0=0.0; x1=1.0; h=(x1-x0)/(n-1.0)
#vectors of tridiagonal matrix A
a = np.ones(n)
b = -2*np.ones(n)
c = np.ones(n)
A = np.diag(a[1:], -1) + np.diag(b, 0) + np.diag(c[:-1], 1)
#vectors y and x
x = np.linspace(x0,x1,n) 
y = h*h*100*np.exp(-10*x)
u_exact = 1-(1-np.exp(-10))*x-np.exp(-10*x)

#forward substitution
b, y = forwardalg_vec(a,b,c,y)
print "finished forward substitution"
#backward substitution
y = backwardalg_vec(b,c,y)
print "finished backward substitution"

#solve diagonal matrix equation
u = y/b.astype(float)

outfile = curdir+"/data/dderiv_python_n%d_v%s.dat"%(n,version)
write2file("x, f, u", filename=outfile, append=False)
for i in range(len(x)):
    write2file("%1.20e, %1.20e, %1.20e"%(x[i], y[i]*h*h, u[i]),
               filename=outfile,
               append=True)

#find relative error
#push u_exact and u up by one to correct for u_exact[i] = 0
u = u - 1
u_exact = u_exact - 1
eps = np.abs((u - u_exact)/u_exact)
