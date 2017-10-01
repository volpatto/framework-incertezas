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
for i in range(ifi,ifi+1):
        name = "%s%s%s" % (filename1, str(i), ext)
        data1 = np.genfromtxt(name)
        t1, alpha = data1[:,0], data1[:,1]
        Vl, Vr = data1[:,4], data1[:,7]
        alphavec = np.append(alphavec, alpha)
        # plot
        plt.plot(t1, Vl, color=colormap(normalize(alpha[0])))
        plt.plot(t1, Vr, color=colormap(normalize(alpha[0])))
        plt.fill_between(t1, Vl, Vr, color=colormap(normalize(alpha[0])))


# setup the colorbar
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(alphavec)
clb = plt.colorbar(scalarmappaple)
clb.ax.set_title(r'$\alpha$')

# show the figure
plt.xlabel(r't (segundos)')
plt.ylabel(r'$V$ (volts)')
plt.title(r'$\alpha$-cut = %s' % (alpha[0]))
plt.grid(True)
plt.ylim(ymin=0)
#plt.legend(loc='best',borderpad=0.5)
plt.savefig('%s%s%s'  % ('plotFillFuzzy', str(int(10.0*alpha[0])),'.eps'))
plt.show()

