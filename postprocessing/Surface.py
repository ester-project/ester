# coding: utf-8
# computes the non-dimensional surface of the star. For spherical star S=4*pi
# The star is axisymmetric

from ester import *
from numpy import *

file_of_star_model='Master/M3_O1.h5'
a=star2d(file_of_star_model)

# We first define fields without equator and poles 
r=a.r[-1,1:-1]
rt=a.rt[-1,1:-1]
weights=a.It[1:-1] # actually we eliminate 0 introduced at pole and equator
# dS is the surface element

dS=2*pi*r**2*sqrt(1+rt**2/r**2)
Surf=dot(weights,dS)

print('The non-dimensional Rpol=1 is %6.3f'%Surf,' S/S_sphere = %6.3f'%(Surf/4/pi))
