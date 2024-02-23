# coding: utf-8
# Compute the volume and the mass of a model

from ester import *
from numpy import *

a=star2d('M3O5.h5')
dv = a.r**2 * a.rz
dm = a.rho*a.r**2 * a.rz

Icheb=a.I.ravel()

# We first integrate over theta on every radial grid point
# remove pole and equator added by ester
Vrad=[]
Mrad=[]
for i in range(a.nr):
    Vrad.append(dot(a.It[1:-1],dv[i,1:-1]))
    Mrad.append(dot(a.It[1:-1],dm[i,1:-1]))

V = 2*pi*dot(Vrad, Icheb)*a.R**3
Mass=2*pi*dot(Mrad, Icheb)*a.rhoc*a.R**3/M_SUN

print('Mass %6.3f'%Mass+' Msun')
print('Vol %6.3e'%V+' cm^3')


