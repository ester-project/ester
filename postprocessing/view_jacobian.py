import matplotlib.pyplot as plt
import numpy as np

mat = np.loadtxt('J-19.txt')
x = mat[:, 0]
y = mat[:, 1]
plt.xlim(xmin=0.,xmax=150.)
plt.ylim(ymin=-150.,ymax=0.)

plt.scatter(x, -y) 
plt.show()
