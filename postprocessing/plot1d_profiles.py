# Script giving an example for plotting profiles of physical quantities
# from 1d models

from ester import *
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os

path=os.path.expanduser( '~/Ester/runs/Master/')
a=star1d(path+'M5.h5') # open a 1D solution

Rayon=a.R/R_SUN
print(Rayon) # print the radius

lnrho=np.log(a.rho)         # ln(rho)
dlnrho=np.dot(a.D,lnrho)    # radial derivative of ln(rho)

lnp=np.log(a.p)
dlnp=np.dot(a.D,lnp)
Pr_rad=a.cp*a.T/5/C_LIGHT**2

phi=a.phi		   # gravitational potential
dphi=np.dot(a.D,phi)	   # gravity

plt.plot(a.r/a.R,dphi,'ro') # plot red dots
#plt.yscale('log')        # if log scale needed
#plt.xscale('log')


plt.show()