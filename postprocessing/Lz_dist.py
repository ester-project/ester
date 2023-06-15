# visualize the distribution of angular momentum among the domains

import sys
from ester import *
import matplotlib.pyplot as plt
from numpy import *
import os

path0=os.path.expanduser( '~/Ester/')

filename=path0+'build/runs/M12_O5_evol.h5_0014'

a=star2d(filename)
jfirst=zeros(a.ndomains,dtype=int)
for i in arange(a.ndomains-1)+1:
       jfirst[i]=jfirst[i-1]+a.npts[i-1]

Omega_scale=sqrt(a.pc/a.rhoc)/a.Rp

Lz=(a.r*sin(a.th))**2*a.w*a.R**5*Omega_scale*a.rho*a.rhoc*a.rz*a.r**2
# Integrateur radial: I(1,nr)
Icheb=a.I.ravel() # 2D array with 1 column ==> 1D array

# We first integrate over theta on every radial grid point
dlz=[]
for  i in range(a.nr):
      dlz.append(dot(a.It,Lz[i,:]))

dlza=2*pi*array(dlz)
# Then we integrate over 'r' on the first ncore grid points
Lz_tot1=dot(dlza,Icheb)
print('Lz total = %6.3e'%Lz_tot1)

# examen Lz par domaine
Lz_dom=[]
r_dom=[]
for i in range(a.ndomains):
#    print(jfirst[i],jfirst[i]+a.npts[i])
    i1=jfirst[i]
    i2=jfirst[i]+a.npts[i]
    dlz=[]
    for j in range(i1,i2):
      dlz.append(dot(a.It,Lz[j,:]))
    dlza=2*pi*array(dlz)

    Icheb_dom=Icheb[i1:i2]
    Lz_dom.append(dot(dlza,Icheb_dom))
    r_dom.append((a.r[i1,0]+a.r[i2-1,0])/2)

Lz_tot2=sum(Lz_dom)
print('Lz total bis = %6.3e'%Lz_tot2)

plt.plot(r_dom,Lz_dom)
plt.yscale('log')

plt.savefig('Lz_dom.png', bbox_inches=0)

plt.show()
