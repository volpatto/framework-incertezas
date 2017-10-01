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

for i in range(1,ifi+1):
	name = "%s%s%s" % (filename1, str(i), ext)
	data1 = np.genfromtxt(name)
	t1, V1 = data1[:,0], data1[:,1]
	Vlast = np.append(Vlast, V1[-1])

n, bins, patches = plt.hist(Vlast,50)
mu = np.mean(Vlast)
sigma = np.std(Vlast)
plt.clf()

# setup the normalization and the colormap
normalize = mcolors.Normalize(vmin=min(mlab.normpdf(bins, mu, sigma)), vmax=max(mlab.normpdf(bins, mu, sigma)))
colormap = cm.jet

plt.rc('text',usetex=True)
for i in range(1,ifi+1):
	name = "%s%s%s" % (filename1, str(i), ext)
	data1 = np.genfromtxt(name)
	t1, V1 = data1[:,0], data1[:,1]
	y1 = mlab.normpdf(V1[-1], mu, sigma)
	plt.plot(t1, V1, color=colormap(normalize(y1)))


# setup the colorbar
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
#scalarmappaple.set_array(alphavec)
scalarmappaple.set_array(mlab.normpdf(bins, mu, sigma))
clb = plt.colorbar(scalarmappaple)
clb.ax.set_title(r'$\rho$')

# show the figure
plt.xlabel(r't (segundos)')
plt.ylabel(r'$V$ (volts)')
plt.title("$\mu = $ %.3e, $\sigma = $ %.3e" % (mu,sigma))
plt.grid(True)
plt.savefig('plotFillNorm.eps')
plt.show()

