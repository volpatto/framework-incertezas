# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib import rc
import sys

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

nSturgis = int(1 + 3.3*np.log10(len(Vlast)))
n, bins, patches = plt.hist(Vlast,nSturgis,normed=1)
plt.ylabel(r'$\rho$ (Densidade de probabilidade)')
plt.xlabel('V')
plt.savefig('plotHistLast.eps')
mu = np.mean(Vlast)
sigma = np.std(Vlast)
xx = np.linspace(min(Vlast),max(Vlast),len(Vlast))
y1 = mlab.normpdf(xx, mu, sigma)
plt.plot(xx, y1, 'r--', linewidth=2,label=("$\mu = $%.3e\n$\sigma = $%.3e" % (mu,sigma)))
plt.legend(loc='best',borderpad=0.5)
plt.savefig('plotHistPDFLast.eps')

