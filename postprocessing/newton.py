# Test de la methode de Newton sur un cas degenere...

import matplotlib.pyplot as plt
import numpy as np

nx=100
n=10

x=0.5
fp=-n*(1-x)**(n-1)
f=(1-x)**n
a=[]
for i in range(nx):
   dx=-f/fp
   x=x+dx
   fp=-n*(1.-x)**(n-1)
   f=(1.-x)**n
   print x
   a.append(x)

b=1.-np.asarray(a)
plt.yscale('log')
plt.plot(b)
plt.show()
