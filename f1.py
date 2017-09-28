import matplotlib.pyplot as plt
import numpy as np
x = np.linspace(0, 100, 50) 
y1 = -(x-50)**2 + 2500
y2 = (200*x-50)
plt.plot(x,y1)
plt.plot(x,y2)
plt.fill_between(x,y1, y2, color='green')
plt.show()
