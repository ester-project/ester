# visualize rho(r) for various co-latitudes theta
# include the path of the python tools in your environment 
# variable something like
# setenv PYTHONPATH $HOME/Ester/local/lib/python2.7/site-packages/ester

import sys
from ester import *
import matplotlib.pyplot as plt
import numpy as np

#a=star1d('M2ideal')
#a=star1d('M199s.h5')
a=star1d('toto.h5')
UMA=1.67353249e-24
RGP=K_BOL/UMA
print 'Teff = ',a.Teff
Teff=(a.L/4/np.pi/a.R**2/SIG_SB)**0.25
print 'Teff calc = ',Teff
plt.close()

Pe=np.reshape(0*a.r,a.nr,1)
cp=np.reshape(a.cp,a.nr,1)

# computes the jfirst grid points
jfirst=np.zeros(a.ndomains,dtype=np.int)

for i in np.arange(a.ndomains-1)+1:
       jfirst[i]=jfirst[i-1]+a.npts[i-1]

try:
 Peclet=a.Peclet
except:
 Peclet=0.

for i in range(a.ndomains-1) : 
  j0=jfirst[i]
  j1=jfirst[i+1]-1
  if (i==0) :
        Pe[j0:j1]=100.*Peclet #*cp[j0:j1]/RGP
  if (i==9) :
        Pe[j0:j1]=Peclet #*cp[j0:j1]/RGP

r=np.reshape(a.r,a.nr,1)
xi=np.reshape(a.conduct,a.nr,1)
del_ad=np.reshape(a.del_ad,a.nr,1)
temp=a.Tc*np.reshape(a.T,a.nr,1)
T=np.reshape(a.T,a.nr,1)
xi_e=xi*(1.+Pe)
xi_t=xi*Pe
gradt=a.Tc*np.reshape(np.dot(a.D,a.T),a.nr,1)/a.R
gradT=np.reshape(np.dot(a.D,a.T),a.nr,1)
grads=np.reshape(np.dot(a.D,a.s),a.nr,1)
gradlnp=np.reshape(np.dot(a.D,a.p)/a.p,a.nr,1)/a.R
#Flux=-4*np.pi*(xi_e*gradt-xi_t*temp*del_ad*gradlnp)*r**2*a.R**2/a.L
Flux=-(Pe*T*grads+gradT)*xi*r**2/xi[300]
Fc=-(Pe*T*grads)*xi*r**2/xi[300]
Frad=-(gradT)*xi*r**2/xi[300]
#Flux=-4*np.pi*xi*gradt*r**2*a.R**2/a.L

nst=250
nf=360
plt.plot(r[nst:nf],Flux[nst:nf])
plt.plot(r[nst:nf],Frad[nst:nf])
plt.plot(r[nst:nf],Fc[nst:nf])
idom=range(a.ndomains+1)
#for i in range(a.ndomains+1) : 
for i in idom[9:10] : 
    plt.axvline(x=a.xif[i],ls=':')

plt.show()
