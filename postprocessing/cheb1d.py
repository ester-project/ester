# Compute the 1D Chebyshev spectrum of the density

from ester import *
import matplotlib.pyplot as plt
import numpy as np


a=star1d('M3.h5') # open a 1D solution

print('Compute the Chebyshev spectrum of the density')

# a.P is the projection matrix of the Gauss-Lobatto grid

x=a.r**2
sp1d=abs(np.dot(a.P,a.rho))
jfirst=np.zeros(a.ndomains,dtype=np.int)

for i in np.arange(a.ndomains-1)+1:
   jfirst[i]=jfirst[i-1]+a.npts[i]

for i in range(a.ndomains):
    sp0=sp1d[jfirst[i]] # normalization by the first coefficient
    sp1d[jfirst[i]:jfirst[i]+a.npts[i],0]=sp1d[jfirst[i]:jfirst[i]+a.npts[i],0]/sp0
    plt.plot(np.arange(a.npts[i])+jfirst[i],sp1d[jfirst[i]:jfirst[i]+a.npts[i]])

ym=min(sp1d)
vertical=np.array([1e-20,1e1])
hori=np.array([0,0])
for i in range(a.ndomains):
    plt.plot(hori+jfirst[i],vertical,'k--')

plt.xlim(xmin=0,xmax=a.nr)
plt.ylim(ymin=ym,ymax=1)
plt.yscale('log')
plt.xlabel(r'grid point of the n$^{\rm th}$ domain + n')
plt.ylabel(r'|a$_n$|')
plt.title('1D Chebyshev spectra')

plt.show()
