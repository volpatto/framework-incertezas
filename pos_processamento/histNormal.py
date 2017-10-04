# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib import rc
import sys

filename1 = sys.argv[1]
varname = sys.argv[2]

data1 = np.genfromtxt(filename1)
V1 = data1[:]

plt.rc('text',usetex=True)
nSturgis = int(1 + 3.3*np.log10(len(V1)))
n, bins, patches = plt.hist(V1,nSturgis,normed=1)
plt.ylabel(r'$\rho$ (Densidade de probabilidade)')
plt.xlabel(varname)
plt.savefig('plotHist%s.eps' % (varname))
#plt.clf()

mu = np.mean(V1)
sigma = np.std(V1)
xx = np.linspace(min(V1),max(V1),len(V1))
y1 = mlab.normpdf(xx, mu, sigma)
plt.plot(xx, y1, 'r--', linewidth=2,label=("$\mu = $%.3e\n$\sigma = $%.3e" % (mu,sigma)))
plt.legend(loc='best',borderpad=0.5)
plt.savefig('plotHistPDF%s.eps' % (varname))

