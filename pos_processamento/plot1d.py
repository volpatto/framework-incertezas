# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import sys

filename = sys.argv[1]
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,1]

plt.rc('text',usetex=True)
#plt.rc('font',family='helvica')
plt.plot(x,y,'.')
#plt.ylim(ymax=1.1)
plt.grid(True)
plt.xlabel(r't (segundos)')
plt.ylabel(r'$V$ (volts)')
plt.savefig('plot1d.eps')
#plt.show()
