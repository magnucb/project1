import pylab as pyl
import os
import sys

curdir = os.getcwd()
data_dict = {} #dictionary of files
n_range_LU = pyl.logspace(1,3,num=3)
n_range_tridiag = pyl.logspace(1, 4, num=4)

for n in n_range_tridiag:
    #loop through different n's
    with open(curdir+"/data/dderiv_u_c++_n%d_tridiag.dat"%(int(n)), 'r') as infile:
        full_file = infile.read() #read entire file into text
        lines = full_file.split('\n') #separate by EOL-characters
        lines = lines[:-1] #remove last line (empty line)
        keys = lines.pop(0).split(', ') #use top line as keys for dictionary
        dict_of_content = {}
        for i,key in zip(range(len(keys)), keys): #loop over keys, create approp. arrays
            dict_of_content[key] = [] #empty list
            for j in range(len(lines)): #loop though the remaining lines of the data-set
                line = lines[j].split(', ') #split the line into string-lists
                word = line[i] # append the correct value to the correct list with the correct key
                try:
                    word = float(word) #check if value can be float
                except ValueError: #word cannot be turned to number
                    print word
                    sys.exit("There is something wrong with your data-file \n'%s' cannot be turned to numbers"%word)
                dict_of_content[key].append(word) 
        data_dict["n=%d"%n] = dict_of_content #add complete dictionary to dictionary of files

def u_exact(x):
    u =  1.0 - (1.0 - pyl.exp(-10.0))*x - pyl.exp(-10.0*x)
    return u

def compare_methods(n):
    """
    For a specific length 'n' compare both methods 
    with the exact function.
    """
    x = pyl.array(data_dict["n=%d"%n]["x"])
    gen = pyl.array(data_dict["n=%d"%n]["u_gen"])
    spec = pyl.array(data_dict["n=%d"%n]["u_spec"])
    exact = u_exact(x)
    pyl.figure("compare methods")
    pyl.grid(True)
    pyl.hold(True)
    pyl.xlabel("x")
    pyl.ylabel("u(x)")
    pyl.title("function u for three different methods (n=%d)"%n)

    pyl.plot(x, exact, 'k-', label="exact")
    pyl.plot(x, gen, 'b--', label="general tridiagonal")
    pyl.plot(x, spec, 'g-', label="specific tridiagonal")
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
        n = int(n)
        x = pyl.array(data_dict["n=%d"%n]["x"])
        u_approx = pyl.array(data_dict["n=%d"%n][approx_key])
        pyl.plot(x, u_approx, '--', label="n=%1.1e"%n)

    x = pyl.linspace(0,1,1001)
    exact = u_exact(x)
    pyl.plot(x, exact, '-', label="exact")
    pyl.legend(loc='best', prop={'size':9})
    pyl.savefig(curdir+"/img/compare_%s_exact.png"%(approx_string))
    
def epsilon_plots(n_range=[10,100,1000]):
    eps_max = pyl.zeros(len(n_range))
    h = pyl.zeros(len(n_range))
    for i, n in enumerate(n_range):
        x = pyl.array(data_dict["n=%d"%n]["x"])
        u = u_exact(x) 
        v = pyl.array(data_dict["n=%d"%n]["u_gen"])
        #calculate eps_max by finding max of |v_i - u_i|
        max_diff_uv = 0; jmax = 0;
        for j in range(int(n)):
            diff_uv = abs(v[j]-u[j])
            if diff_uv > max_diff_uv:
                max_diff_uv = diff_uv
                jmax = j
        if jmax == 0 or jmax == n-1:
            sys.exit("There is an error in calculating the max_epsilon")
        eps_max[i] = pyl.log10(max_diff_uv/float(abs(u[jmax])))
        h[i] = pyl.log10(1.0/(n+1))
        print "%d & %1.4f & %1.4f \\\\"%(n, h[i], eps_max[i])
    pyl.figure("epsilon")
    pyl.grid(True)
    pyl.hold(True)
    pyl.xlabel(r"\log_{10}(h)")
    pyl.ylabel(r"$\epsilon = \log_{10}\left(\frac{u_{approx}-e_{exact}}{u_{exact}}\right)$")
    pyl.title("log-plot of epsilon against step-length h")
    pyl.plot(h, eps_max, 'k-')
    pyl.legend(loc='best')
    pyl.savefig(curdir+"/img/epsilon.png")

#make plots
#compare_methods(n=10)
compare_approx_n(approx_string="general") #exercise b
epsilon_plots(n_range=pyl.logspace(1, 4, num=4)) #exercise d
pyl.show()
