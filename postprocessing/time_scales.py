#
# visualize the distribution of angular momentum among the domains
# and now its evolution through a series of models


import os
import subprocess as cmd

import matplotlib
#matplotlib.use('Tkagg')

import sys
from ester import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import integrate
from scipy.special import sph_harm
import glob, h5py 
fontsize = 14

matplotlib.rc('xtick', labelsize=fontsize)
matplotlib.rc('ytick', labelsize=fontsize)


def volume_integral(a,f):

    ncore=a.npts[0] # we assume 1 domain in the core

    # Integrateur radial: I(1,nr)
    b=a.I.ravel() # 2D array with 1 column ==> 1D array
    Icheb=b[ncore:]

    dml=[]
    # We first integrate over theta on every radial grid point
    for  i in range(ncore,a.nr): 
        dm=np.dot(a.It,f[i,:])
        dml.append(dm)

    dma=2*np.pi*np.array(dml)
    # Then we integrate over 'r' on the first ncore grid points
    M=np.dot(dma,Icheb)
    return M


#ls =cmd.getoutput('ls -v M3_O8nth24_ev*').split()
#ls =cmd.getoutput('ls -v M12_O5_evol.h5_*').split()
ls =cmd.getoutput('ls -v M12_O15ev*').split()

mH  = 1.007825
mHe = 4.0026033
MYR=86400*365.25*1e6
Tratio   = []
star_age = []
ic=0
for fichier in ls:
    print()
    print(fichier)
    a=star2d(fichier)
    rho=a.rho*a.rhoc

    Omega_scale=np.sqrt(a.pc/a.rhoc)/a.Rp
    Km = 16*SIG_SB*(a.T*a.Tc)**3/(3*a.opacity*rho**2*a.cp)
    N2 = a.N2
    N2[N2 < 0] = 0
    # We compute <N^2/Omega^2*K_diff>
    f=N2/(a.w*Omega_scale)**2/Km*a.rz**(a.r)**2
    v=a.rz*(a.r)**2
    mean_NOK = volume_integral(a, f)
    vol      = volume_integral(a, v)
    mean_NOK = mean_NOK/vol

    # We compute <N^2/Omega^2> weighted by mass
    f=N2/(a.w*Omega_scale)**2*a.rz**(a.r)**2*rho
    dm=a.rz*(a.r)**2*rho
    mean_NO = volume_integral(a, f)
    mass    = volume_integral(a, dm)
    mean_NO = mean_NO/mass
    eta=a.r[a.npts[0],0] # scaled radius of the core at equator

    if ic == 0:
       Kmm=Km*rho*a.Rp**3
       mean_kappa=volume_integral(a,Km)/vol
       f=rho*a.r**2*a.rz*a.Rp**3
       Mass_env=volume_integral(a,f)
       m_kappa=volume_integral(a,Kmm)/Mass_env
       print()
       print('Volume averaged and mass averaged heat diffusivity:')
       print('kappa vol = %6.3e'%mean_kappa,'kappa mass = %6.3e'%m_kappa)
       print('R_core/R_star = %6.3f'%eta,' at equator')

    T = mean_NOK*a.R**2/MYR / (np.pi**2 + 4)/np.pi**2*(1-eta)**2
    Tmass = mean_NO*a.R**2/MYR/m_kappa / (np.pi**2 + 4)/np.pi**2*(1-eta)**2

    Q = (4*mH-mHe)*C_LIGHT**2
    dXdt = 4*mH*MYR/Q*a.eps[0,0]
    Tx = a.X[0,0]/dXdt
    print('Xc = %6.3f'%a.Xc,'Tbaro mass = %5.2e Myr'%Tmass,'Tbaro vol = %5.2e Myr'%T, \
              ' Tnuc = %5.2e Myr'%Tx, \
             'Tnuc/Tbaro_vol = %6.3f'%(Tx/T))
    Tratio.append(Tx/T)
    star_age.append(a.Xc)
    T_KH=GRAV*a.M**2/a.R/a.L/MYR # Kelvin-Helmholtz
    Omega_eq=a.w[-1,0]*Omega_scale
    T_ES=T_KH*GRAV*a.M/a.R**3/Omega_eq**2 # Eddington-Sweet
    T_diff=a.R**2/m_kappa/MYR # diffusion mass-weighted
    print('T ES = %5.2e Myr'%T_ES,' T_KH = %5.2e Myr'%T_KH, \
    ' T_diff = %5.2e Myr'%T_diff,'Tnuc/T_ES = %6.3f'%(Tx/T_ES))
    ic+=1

plt.figure()
plt.plot(star_age, Tratio)
plt.show()
