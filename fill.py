# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import sys

filename1 = sys.argv[1]
filename2 = sys.argv[2]
data1 = np.genfromtxt(filename1)
data2 = np.genfromtxt(filename2)
t1, V1 = data1[:,0], data1[:,1]
t2, V2 = data2[:,0], data2[:,1]

plt.rc('text',usetex=True)
#plt.rc('font',family='helvica')
plt.plot(t1,V1,'b')
plt.plot(t2,V2,'b')
plt.fill_between(t1,V1,V2,color='blue')
#plt.ylim(ymax=1.1)
plt.grid(True)
plt.xlabel(r't (segundos)')
plt.ylabel(r'$V$ (volts)')
plt.savefig('plotFill.eps')
#plt.show()
