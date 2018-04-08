import sys
import matplotlib.pyplot as plt
import numpy as np


#plt.close()
a = [[], []] # A list with 2 empty lists
b = [[], []] # A list with 2 empty lists
c = [[], []] # A list with 2 empty lists
d = [[], []] # A list with 2 empty lists
f = open('norm.txt', 'r')
k=0
for line in f:
      data = line.split()
      if (k % 4 == 0):
       for i, value in enumerate(line.split()):
           a[i].append(float(data[i]))
      if (k % 4 == 1):
       for i, value in enumerate(line.split()):
           b[i].append(float(data[i]))
      if (k % 4 == 2):
       for i, value in enumerate(line.split()):
           c[i].append(float(data[i]))
      if (k % 4 == 3):
       for i, value in enumerate(line.split()):
           d[i].append(float(data[i]))
      #print data[0],data[1]
      k=k+1
pa, =plt.plot(a[0][:],a[1][:])
pb, =plt.plot(b[0][:],b[1][:])
pc, =plt.plot(c[0][:],c[1][:])
pd, =plt.plot(d[0][:],d[1][:])
plt.yscale('log')
plt.xlabel('Iterations')
plt.ylabel('Norm2')
plt.legend([pa,pb,pc,pd],["Phi","Pres","Temp","Map"],loc=3)

plt.show()
