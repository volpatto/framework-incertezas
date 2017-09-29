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
Vlast = np.array([])

plt.rc('text',usetex=True)
for i in range(1,ifi+1):
	name = "%s%s%s" % (filename1, str(i), ext)
	data1 = np.genfromtxt(name)
	t1, V1 = data1[:,0], data1[:,1]
	Vlast = np.append(Vlast, V1[-1])

plt.hist(Vlast,100,normed=True,stacked=True)
plt.ylabel(r'Frequ\^encia')
plt.xlabel(r'$V$ (volts)')
plt.savefig('plotHistLast.eps')
plt.show()

