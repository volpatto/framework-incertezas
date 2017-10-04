# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib import rc
import sys

def uniform(x,a,b):
    if (x<=b and x>=a):
        return 1.0/(b-a)
    else:
        return 0.0

filename1 = sys.argv[1]
varname = sys.argv[2]
ext = ".dat"
data1 = np.genfromtxt(filename1)
V1 = data1[:]

plt.rc('text',usetex=True)
Va = min(V1)
Vb = max(V1)
nSturgis = int(1 + 3.3*np.log10(len(V1)))
n, bins, patches = plt.hist(V1,nSturgis,normed=1)
v_uniform = np.vectorize(uniform)
plt.plot(bins, v_uniform(bins,Va,Vb), 'r--', linewidth=2,label=("$\\displaystyle\\frac{1}{b-a}$ = %.2e \n a = %.2e, b = %.2e" % (1.0/(Vb-Va),Va,Vb)))
plt.legend(loc='best',borderpad=0.5)
plt.ylabel(r'$\rho$ (densidade de probabilidade)')
plt.xlabel(varname)
plt.savefig('plotHistLastU%s.eps' % varname)
plt.show()

