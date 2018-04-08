import sys
import matplotlib.pyplot as plt
import numpy as np


#plt.close()
a = [[], []] # A list with 2 empty lists
f = open('norm.txt', 'r')
k=0
for line in f:
      data = line.split()
      for i, value in enumerate(line.split()):
           a[i].append(float(data[i]))
      #print data[0],data[1]
plt.plot(a[0][:],a[1][:])
plt.yscale('log')
plt.xlabel('Iterations')
plt.ylabel('Norm2')
plt.show()
