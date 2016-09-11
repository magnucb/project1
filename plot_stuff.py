import pylab as pyl
import os
import sys

curdir = os.getcwd()


def plot_generator(version, n):
    """
    creates a plot of the generated data
    """
    datafile = open(curdir+"/data/dderiv_u_python_v%s_n%d.dat"%(version,n))
    biglist = []

    for line in datafile:
        linesplit = line.split
        biglist.append(linesplit)