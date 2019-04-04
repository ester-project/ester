# Plot of the Brunt-Vaisala frequency in a meridional plane
# convection zones appear as whiie places.

from ester import *
import matplotlib.pyplot as plt
import numpy as np

a=star2d('/home/rieutord/Ester/local/runs/Master/M3O7.h5')

Req=a.Re
Rpol=a.Rp
print Req/R_SUN,Rpol/R_SUN

x=a.r*np.sin(a.th)
y=a.r*np.cos(a.th)
x=np.c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
y=np.c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
z=np.log10(np.c_[a.N2,a.N2[:,-2::-1],a.N2[:,1:],a.N2[:,-2::-1]])
plt.figure(1,figsize=(7,7))
plt.contourf(x,y,z,60,cmap=plt.cm.RdBu)
plt.axis([-1.1,1.1,-1.1,1.1])

plt.show()
