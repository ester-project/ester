# Hints for plotting 1D quantities with domains interfaces

import sys
from ester import *
import matplotlib.pyplot as plt
import numpy as np

a=star1d('M5.h5')
plt.close()
#plt.plot(a.r[:],np.log10(a.rho[:]))
#gradt=np.dot(a.D,a.T)

#plt.plot(a.r[240:],a.s[240:]/a.s[0])
plt.plot(a.r,a.s)
for i in range(a.ndomains+1) : 
    plt.axvline(x=a.xif[i],ls=':')

plt.legend()
#x=[]
#x.append(a.r[30])
#x.append(a.r[30])
#x=x/1e11
#y=[0,1]
#plt.plot(x,y)
plt.show()
