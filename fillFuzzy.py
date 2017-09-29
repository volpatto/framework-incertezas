# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib import rc
import sys

filename1 = sys.argv[1]
ifi = int(sys.argv[2])
ext = ".dat"
alphavec = np.array([])

# setup the normalization and the colormap
normalize = mcolors.Normalize(vmin=0.0, vmax=1.0)
colormap = cm.jet

plt.rc('text',usetex=True)
for i in range(1,ifi+1):
	name = "%s%s%s" % (filename1, str(i), ext)
	data1 = np.genfromtxt(name)
	t1, alpha = data1[:,0], data1[:,1]
	Vl, Vr = data1[:,4], data1[:,7]
	alphavec = np.append(alphavec, alpha)
	# plot
	plt.plot(t1, Vl, color=colormap(normalize(alpha[0])))
	plt.plot(t1, Vr, color=colormap(normalize(alpha[0])))


# setup the colorbar
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(alphavec)
clb = plt.colorbar(scalarmappaple)
clb.ax.set_title(r'$\alpha$')

# show the figure
plt.xlabel(r't (segundos)')
plt.ylabel(r'$V$ (volts)')
plt.grid(True)
plt.savefig('plotFillFuzzy.eps')
plt.show()

