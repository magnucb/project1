import pylab as pyl
import numpy as np
import os
import sys

curdir = os.getcwd()
version = 3
data_dict = {}
column_dict = {}
n_range = [10,100,1000]

for n in n_range:
    #loop through different n's
    with open(curdir+"/data/dderiv_u_python_v%s_n%d.dat"%(version,n), 'r') as infile:
        full_file = infile.read() #read entire file into text
        lines = full_file.split('\n') #separate by EOL-characters
        lines = lines[:-1] #remove last line (empty line)
        keys = lines.pop(0).split(', ') #use top line as keys for dict.
        dict_of_content = {}
        for i,key in zip(range(len(keys)), keys):
            #loop over keys: h, f2c, f3c
            dict_of_content[key] = []
            for j in range(len(lines)):
                # loop over all lines
                line = lines[j].split(', ')
                word = line[i]
                try:
                    word = float(word)
                except ValueError: #word cannot be turned to number
                    print word
                    sys.exit("There is something wrong with your data-file \n'%s' cannot be turned to numbers"%word)
                dict_of_content[key].append(word)
        data_dict["n=%d"%n] = dict_of_content

def u_exact(x):
    u =  1.0 - (1.0 - pyl.exp(-10.0))*x - pyl.exp(-10.0*x)
    return u

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

def harry_plotter():
    #plot pre-game
    pyl.figure()
    pyl.grid(True)
    pyl.title("function u for different steplengths")
    pyl.ylabel("u(x)")
    pyl.xlabel("x")
    for n in n_range:
        x = pyl.array(data_dict["n=%d"%n]["x"])
        u_gen = pyl.array(data_dict["n=%d"%n]["u_gen"])
        u_spec = pyl.array(data_dict["n=%d"%n]["u_spec"])
        u_LU = pyl.array(data_dict["n=%d"%n]["u_LU"])
        pyl.plot(x,u_gen, 'g-', label="general tridiag, n=%d"%n)
        pyl.plot(x,u_spec, 'r-', label="specific tridiag, n=%d"%n)
        pyl.plot(x,u_LU, 'b-', label="LU-dekomp, n=%d"%n)
        
    u_ex = u_exact(x)
    pyl.plot(x, u_ex, 'k-', label="exact, n=%d"%len(x))
    #pyl.legend(loc="best",prop={"size":8})
    pyl.show()
    return None

def compare_methods(n):
    """
    For a specific length 'n' compare all three methods 
    with the exact function.
    """
    x = pyl.array(data_dict["n=%d"%n]["x"])
    gen = pyl.array(data_dict["n=%d"%n]["u_gen"])
    spec = pyl.array(data_dict["n=%d"%n]["u_spec"])
    LU = pyl.array(data_dict["n=%d"%n]["u_LU"])
    exact = u_exact(x)
    pyl.figure("compare methods")
    pyl.grid(True)
    pyl.hold(True)
    pyl.xlabel("x")
    pyl.ylabel("u(x)")
    pyl.title("function u for three different methods (n=%d)"%n)

    pyl.plot(x, exact, 'k-', label="exact")
    pyl.plot(x, gen, 'b-', label="general tridiagonal")
    pyl.plot(x, spec, 'g-', label="specific tridiagonal")
    pyl.plot(x, LU, 'r-', label="LU-decomp.")
    pyl.legend(loc='best', prop={'size':9})
    pyl.savefig(curdir+"/img/compare_methods_n%d.png"%n)
               
def compare_approx_n(n_range=[10,100,1000], approx_string="general"):
    """
    For all n's available, plot the general approximation and 
    exact solution
    """
    if approx_string == "general":
        approx_key = "u_gen"
    elif approx_string == "specific":
        approx_key = "u_spec"
    else:
        sys.exit("In function 'compare_approx_n', wrong argument 'approx_string'")
    pyl.figure("compare %s"%approx_string)
    pyl.grid(True)
    pyl.hold(True)
    pyl.xlabel("x")
    pyl.ylabel("u(x)")
    pyl.title("approximation by %s tridiagonal method"%approx_string)

    for n in n_range:
        x = pyl.array(data_dict["n=%d"%n]["x"])
        u_approx = pyl.array(data_dict["n=%d"%n][approx_key])
        pyl.plot(x, u_approx, '--', label="n=%1.1e"%n)

    x = pyl.linspace(0,1,1001)
    exact = u_exact(x)
    pyl.plot(x, exact, '-', label="exact")
    pyl.legend(loc='best', prop={'size':9})
    pyl.savefig(curdir+"/img/compare_%s_n_n%d.png"%(approx_string,n))
    
def epsilon_plots(n_range=[10,100,1000]):
    eps_max = pyl.zeros(len(n_range))
    h = pyl.zeros(len(n_range))
    for i, n in enumerate(n_range):
        x = pyl.array(data_dict["n=%d"%n]["x"])
        u = u_exact(x) #+1.0 #add one to avoid divide-by-zero
        v = pyl.array(data_dict["n=%d"%n]["u_gen"]) #+1.0 #add one to avoid divide-by-zero
        #print max(v), min(v)
        #print max(v), min(v)
        #print abs((max(v)-max(u))/max(u)), abs((min(v)-min(u))/min(u))
        eps = pyl.zeros(len(v)) #eps = pyl.log10(abs((v-u)/u))
        for j in range(len(eps)):
            eps[j] = pyl.log10(abs((v[j]-u[j])/u[j]))
            print abs((v[j]-u[j])/u[j])
            print u[j]                      
            print eps[j]
            print " "
        eps_max[i] = max(eps)
        sys.exit()
        h[i] = pyl.log10(1.0/(n+1))
    pyl.figure("epsilon")
    pyl.grid(True)
    pyl.hold(True)
    pyl.xlabel(r"\log_{10}(h)")
    pyl.ylabel(r"$\epsilon = \log_{10}\left(\frac{u_{approx}-e_{exact}}{u_{exact}}\right)$")
    pyl.title("log-plot of epsilon against step-length h")
    pyl.plot(h, eps_max, 'ko')
    print h
    print eps_max
    pyl.legend(loc='best')
    pyl.savefig(curdir+"/img/epsilon.png")

#make plots
compare_methods(n=1000)
#pyl.show()
compare_approx_n(approx_string="general")
compare_approx_n(approx_string="specific")
epsilon_plots()
#pyl.show()
