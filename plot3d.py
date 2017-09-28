# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys

filename = sys.argv[1]
data = np.genfromtxt(filename)
x, y, z = data[:,0], data[:,1], data[:,2]

XX, YY, ZZ = np.meshgrid(x,y,z)
#ZZ = z.reshape(XX.shape)
#ZZ = np.meshgrid(z)

fig = plt.figure()
ax = fig.gca(projection='3d')
#surf = ax.plot_surface(XX, YY, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
print ZZ
plt.imshow(ZZ)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
