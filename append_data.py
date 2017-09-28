# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
name_new = sys.argv[3]
data1 = np.genfromtxt(file1)
data2 = np.genfromtxt(file2)

np.savetxt(name_new,np.transpose([data1[:,0],data2[:,0],data1[:,1]]))
