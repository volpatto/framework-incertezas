# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
#rc('text',usetex=True)

def alpha(x,r1,r1m,r2p,r2):
    if (x>=r1 and x<r1m):
        a = 1./4.37e3
        #b = -4.3257e5/4.43
        b = -4.3263e2/4.37
        return a*x+b
    elif (x>r2p and x<=r2):
        a = -1./4.43e3
        #b = 4.4737e5/4.37
        b = 4.4743e2/4.43
        return a*x+b
    else:
        return 1.0

r1 = 4.3263e5
r1m = 4.37e5
r2p = 4.43e5
r2 = 4.4743e5
x = np.linspace(r1,r2,200)
valpha = np.vectorize(alpha)

plt.rc('text',usetex=True)
#plt.rc('font',family='helvica')
plt.plot(x,valpha(x,r1,r1m,r2p,r2))
plt.ylim(ymax=1.1)
plt.grid(True)
plt.xlabel(r'R ($\Omega$)')
plt.ylabel(r'$\alpha$')
plt.show()
