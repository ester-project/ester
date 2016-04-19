# Comoute the 1D Chebyshev spectrum of the density

from ester import *
import matplotlib.pyplot as plt
import numpy as np


a=star1d('M5paper_OPAL') # open a 1D solution

print 'Compute the Chebyshev spectrum of the density'

# a.P is the projection matrix of the Gauss-Lobatto grid

sp1d=abs(np.dot(a.P,a.rho))
jfirst=np.zeros(a.ndomains)

for i in np.arange(a.ndomains-1)+1:
   jfirst[i]=jfirst[i-1]+a.npts[i]

for i in range(a.ndomains):
    sp0=sp1d[jfirst[i]] # normalization by the first coefficient
    plt.plot(np.arange(a.npts[i])+jfirst[i],sp1d[jfirst[i]:jfirst[i]+a.npts[i]]/sp0)

vert=np.array([1e-20,1e1])
hori=np.array([0,0])
for i in range(a.ndomains):
    plt.plot(hori+jfirst[i],vert,'k--')

plt.xlim(xmin=0,xmax=a.nr)
plt.ylim(ymin=1e-15,ymax=1)
plt.yscale('log')
plt.xlabel('grid point of nth domain and n')
plt.ylabel(r'a$_n$')
plt.title('1D Chebyshev spectra')

plt.show()
