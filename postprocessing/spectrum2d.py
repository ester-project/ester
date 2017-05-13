# Compute the 2D Chebyshev Spherical Harmonic spectrum of the density

from ester import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

filename='M5_O95_nth30'
a=star2d(filename) # open a 2D solution

print 'Compute the 2D spectrum of the density'

# a.P is the projection matrix of the Gauss-Lobatto grid
# a.P_00 is the projection matrix of the Gauss-Legendre grid for a field
# symmetric with respect to equator


sp_leg=np.zeros((a.nr,a.nth))
print 'check the shape',sp_leg.shape

# Gauss Legendre spectrum
# we eliminate the pole/equator points that are respectively at a.nth+1
# and 0

for i in range(a.nr):
    sp_leg[i,:]=np.dot(a.P_00,a.rho[i,1:a.nth+1]) 

for k in range(a.nth):
    sp_leg[:,k]=abs(np.dot(a.P,sp_leg[:,k]))

# Normalisation by the coefficient (0,0) in each domain
jfirst=np.zeros(a.ndomains,dtype=np.int)

for i in np.arange(a.ndomains-1)+1:
   jfirst[i]=jfirst[i-1]+a.npts[i]

for i in range(a.ndomains):
    sp0=sp_leg[jfirst[i],0] # normalization by the (n=0,l=0) coefficient
    sp_leg[jfirst[i]:jfirst[i]+a.npts[i],:]=sp_leg[jfirst[i]:jfirst[i]+a.npts[i],:]/sp0

x=np.arange(a.nth)
y=np.arange(a.nr)

rap=float(a.nth)/a.nr
cax = plt.matshow(np.log10(abs(sp_leg)),origin='lower',aspect=rap)
cbar = plt.colorbar(cax,fraction=0.05,pad=0.1)
# fraction = 5 % and pad = white padding between bar and plot in inches!

#plt.grid(True)
#plt.xlabel('2D spectrum')
plt.title('Index of spherical harmonic')
plt.ylabel(r'Chebyshev index + domain rank')

li=[str(2*x[i]) for i in range(a.nth)]
plt.xticks(x,li)
plt.yticks(jfirst)


fig = plt.gcf()
fig.set_size_inches(12., 10.,forward=True)
plt.savefig(filename+'_sp2D.png',dpi=100, bbox_inches=0)

plt.show()
