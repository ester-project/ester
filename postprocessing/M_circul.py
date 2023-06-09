# visualize streamline of the meridional circulation in (r,theta)
# include the path of the python tools in your environment 
# variable something like
# setenv PYTHONPATH $HOME/Ester/local/lib/python2.7/site-packages/ester

import sys
from ester import *
from numpy import *
import matplotlib.pyplot as plt

filename='M20O3nth40.h5'
a=star2d(filename)


nth=a.nth
Req=a.Re/R_SUN
Rpol=a.Rp/R_SUN
print('check Req = %6.3f'%Req,'and Rpol = %6.3f'%Rpol,'in Rsun')
Omega_scale=sqrt(a.pc/a.rhoc)/a.Rp
omega=Omega_scale*a.w
print('omega pole = %6.5e'%omega[-1,nth+1])
print('omega equator = %6.5e'%omega[-1,0])

# the cartesian grid
x=a.r*sin(a.th)
y=a.r*cos(a.th)
# the stellar surface
xb=x[-1,:]
yb=y[-1,:]
x=c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
y=c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
z=c_[a.G,a.G[:,-2::-1],a.G[:,1:],a.G[:,-2::-1]]
# G is not the true stream function but  R*sin(theta)*G is a true
# meridional stream function so
Psi=x*z

M=a.M/M_SUN
param='M = %4.1f'%M+' omega = %4.2f'%a.Omega_bk
plt.figure(1,figsize=(7,7))
plt.contour(x,y,Psi,30,colors='r')
plt.title('meridional streamlines '+param)
plt.plot(xb,yb,'-',color='k')
plt.plot(xb,-yb,'-',color='k')
plt.plot(-xb,yb,'-',color='k')
plt.plot(-xb,-yb,'-',color='k')

plt.axis([-5.5,5.5,-5.5,5.5])
plt.axis('scaled')
plt.savefig(filename+'_Psi.png', bbox_inches=0)

plt.show()
