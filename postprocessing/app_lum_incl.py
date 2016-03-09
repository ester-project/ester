# Script that plots the apparent luminosity of a rotating star as a
# function of the inclination of the line-of-sight with respect to the
# rotation axis. Units are incl. angle in degree, apparent luminosity in
# solar units.
# Note: app_lum=4*pi*d^2*flux_received_at distance 'd' along the
# line-of-sight characterized by the angle i, limb darkening is not
# taken into account. Details in star2d_extra.cpp
# Author : M. Rieutord (9/3/2016)

from ester import *
import matplotlib.pyplot as plt

np=20
lum_app=[0.0 for i in range(np)]
angle=[float(i)/float(np-1)*90 for i in range(np)]
a=star2d('Achernar2d')
print 'Luminosity of the star (solar unit) =',a.L/L_SUN
print ' '

k=0
for i in angle:
	app_l= a.apparent_luminosity(i)
	print i, 'angle in degrees, app lum =',app_l
	lum_app[k]=app_l
	k=k+1


plt.plot(angle,lum_app,'b-')

plt.xlabel(r'Line of sight angle (degree)')     #labeling the plot
plt.ylabel(r'Apparent luminosity (Lsun)')

plt.show()
