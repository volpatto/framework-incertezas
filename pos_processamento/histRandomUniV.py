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
ifi = int(sys.argv[2])
ext = ".dat"
Vlast = np.array([])

plt.rc('text',usetex=True)
for i in range(1,ifi+1):
	name = "%s%s%s" % (filename1, str(i), ext)
	data1 = np.genfromtxt(name)
	t1, V1 = data1[:,0], data1[:,1]
	Vlast = np.append(Vlast, V1[-1])

Va = min(Vlast)
Vb = max(Vlast)
n, bins, patches = plt.hist(Vlast,50,normed=1)
v_uniform = np.vectorize(uniform)
plt.plot(bins, v_uniform(bins,Va,Vb), 'r--', linewidth=2,label=("$\\displaystyle\\frac{1}{b-a}$ = %.2e \n a = %.2e, b = %.2e" % (1.0/(Vb-Va),Va,Vb)))
plt.legend(loc='best',borderpad=0.5)
plt.ylabel(r'$\rho$ (densidade de probabilidade)')
plt.xlabel(r'$V$ (volts)')
plt.savefig('plotHistLastU.eps')
plt.show()

