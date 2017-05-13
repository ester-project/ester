# visualize omega(r,theta)
# include the path of the python tools in your environment 
# variable something like
# setenv PYTHONPATH $HOME/Ester/local/lib/python2.7/site-packages/ester

import sys
from ester import *
import matplotlib.pyplot as plt
import numpy as np

grav = 6.67384e-8;

filename='M5_O95'
a=star2d(filename)


nth=a.nth
Req=a.Re/R_SUN
Rpol=a.Rp/R_SUN
print 'check Req and Rpol',Req,Rpol,'in Rsun'
Omega_scale=np.sqrt(a.pc/a.rhoc)/a.Rp
omega=Omega_scale*a.w
print 'omega pole',omega[-1,nth+1]
print 'omega equator',omega[-1,0]

x=a.r*np.sin(a.th)/R_SUN
y=a.r*np.cos(a.th)/R_SUN
x=np.c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
y=np.c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
z=np.c_[omega,omega[:,-2::-1],omega[:,1:],omega[:,-2::-1]]
plt.figure(1,figsize=(7,7))
plt.contourf(x,y,z,60,cmap=plt.cm.RdBu)
plt.title('Rotation rate (rad/s)')
plt.colorbar()
plt.axis([-5.5,5.5,-5.5,5.5])
plt.axis('scaled')
plt.savefig(filename+'_w.png', bbox_inches=0)

plt.show()
