import numpy as np
import astropy as ap
import matplotlib.pylab as pl
import sys
import os
import scipy.linalg.lu as LU

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

def general_tridiag(tri_bottom, tri_mid, tri_top, vert):
    #args: arrays for tridiagonal matrix (below, on and above), array for vertical solution
    for i in range(1,len(vert)):
        k = tri_bottom[i]/float(tri_mid[i-1])
        tri_mid[i] -= k*tri_top[i-1]
        vert[i] -= k*vert[i-1]
    i = len(vert) - 2
    while i > 0:
        k = tri_top[i]/float(tri_mid[i+1])
        vert[i] -= k*vert[i+1]
        i -= 1
    return vert

def specific_tridiag(vert):
    #arg: vertical array of solution
    n = len(vert) #size of matrix/arrays
    u = np.zeros(n) #solution u of poisson equation
    d = np.zeros(n); d[0] = -2.#array of diagonals
    #forward subst.
    for i range(1, n):
        k = -(i)/float(d[i-1]) # 1flop
        vert[i] -= k*vert[i-1] # 2flop
        d[i] = -(i+1)/float(i) #(2.5flops)
    #solve upper triangular equation s
    u[-1] = vert[-1]/float(d[-1]) # 1flop
    for i in range(1,n):
        u[-1-i] = (vert[-1-i] - u[-i])/float(d[-1-i]) # 2flops
    return u

def general_LU_decomp(A_matrix, vert):
    #arg: matrix of linear equation, array of solution
    n = len(vert)
    P, L, U = LU(a = A_matrix, overwrite_a = False, check_finite = True)#scipy-function
    # solve L*w = y
    w = np.zeros(n); w[:] = vert[:] #make array w equal to vert
    for i in range(1,n):
        for j in range(0,i):
            w[i] -= L[i,j]*w[j] #modify w according to 'w_i = y_i - sum(l_ij*w_j)'
    u= np.zeros(n); u[:] = w[:] #make array y equal to w
    #TODO DEFINITIVELY BUGS IN THE INDEXECTION BELOW
    for i in range(2,n):
        for j in range(i,n):
            u[-1*i] -= U[-1*i,j]*u[j]/U[-1*i,-1*i] 
    return u

def test_diag(d):
    #d must be 1-D array
    n = len(d)
    exact = np.zeros(n)
    for i in range(n):
        exact[i] = - float(i+2)/float(i+1)
    return abs((d-exact)/exact)

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

#double check diagonal of matrix
diag_error = test_diag(b)
max_diag_error = max(diag_error)
if max_test_diag >= 1e-14:
    print "maximum error in diagonal is: e=", max_test_diag

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
u = u + 1
u_exact = u_exact + 1
eps = np.abs((u - u_exact)/u_exact)

outfile = curdir+"/data/epsilon_python_n%d_v%s.dat"%(n,version)
if os.path.isfile(outfile):
    write2file("n, eps, h", filename=outfile, append=False)
write2file("%0.2e, %0.20e, %0.20e"%(e,eps,h),
           filename=outfile,
           append=True)
