import pylab as pyl
import os
import sys

def plot_generator(version, n):
    """
    plot generator of generated data
    """
    datafile = open(curdir+"/data/dderiv_u_python_v%s_n%d.dat"%(version,n))
    data = []

    for line in datafile:
        linesplit = [item.replace(",","") for item in line.split()]
        data.append(linesplit)

    columns = data[0]
    data    = pyl.array(data[1:]).astype(pyl.float64)
    for i in xrange(len(columns)-1):
        pyl.figure() # comment out this line to unify the plots ... when their dimensions correlate
        pyl.plot(data[:,0], data[:,i+1], label=r"%s" % columns[i+1])
        pyl.xlabel("x")
        pyl.ylabel(r"%s" % columns[i+1])
        pyl.title(r"Plot\ of\ %s\ over\ x" % columns[i+1])
        pyl.legend(loc='best')
    
    # pyl.savefig("evil_plot.png", dpi=400)
    pyl.show()
    datafile.close()

curdir = os.getcwd()
version = int(sys.argv[1])
n       = int(sys.argv[2])
plot_generator(version, n)