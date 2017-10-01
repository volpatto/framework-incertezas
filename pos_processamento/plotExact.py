# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import sys

def exata(t,y0,R,C):
    return y0*np.exp(-t/(R*C))

trange = np.linspace(0.0, 10.0, 200)
R = 1.0
C = 3.0
V0 = 1.0

filename = sys.argv[1]
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,1]

plt.rc('text',usetex=True)
#plt.rc('font',family='helvica')
plt.plot(x,y,'.',label=r'Runge-Kutta')
plt.plot(trange,exata(trange,V0,R,C),'r',label=r"Anal\'itica")
#plt.ylim(ymax=1.1)
plt.grid(True)
plt.xlabel(r't (segundos)')
plt.ylabel(r'$V$ (volts)')
plt.legend(loc='best',borderpad=0.5)
plt.savefig('plot1d.eps')
#plt.show()
