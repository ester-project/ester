# Python script to view the stellar surface from a polar orbit
# 360 images are produced to make a film of the views from a polar
# orbit around the star.
# Use the script "make_the_film" to generate the avi file.

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
import numpy as np
import string
from ester import *

a=star2d('Achernar2d')

nth=a.nth+2 # we take care that the plotted fields as read by "a" have
            # nth+2 points in latitude
u = np.linspace(0, 2 * np.pi, 100) # angle phi
v = a.th[-1,:] # angle theta
w=v+np.pi/2
v=np.concatenate([w[0:nth-1],v]) # angle theta pi (South pole) to 0

r_sud=a.r[-1,:]
r_nord=a.r[-1,nth:0:-1]
r=np.concatenate([r_nord,r_sud])/a.Re

teff_rsh=np.reshape(a.Teff,nth)
teff_ext=np.concatenate([teff_rsh[nth:0:-1],teff_rsh])


x = np.outer(np.cos(u), r*np.sin(v))
y = np.outer(np.sin(u), r*np.sin(v))
z = np.outer(np.ones(np.size(u)), r*np.cos(v))
Teff = np.outer(np.ones(np.size(u)), teff_ext)

N=1-(Teff/Teff.max())**4 # so that white of "Blues" shows the max of Teff
m=cm.ScalarMappable(cmap=cm.Blues, norm = pltc.Normalize(vmin=0,vmax=1.,clip=False))

m.set_array(N) 

for i in range(0,360):
	print(i,'out of 360 degrees')
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=m.to_rgba(N),linewidth=0, antialiased=True, shade=True, alpha=0.9) # alpha = opacity of the surface
	ax.set_xlim(-1,1)
	ax.set_ylim(-1,1)
	ax.set_zlim(-1,1)
	ax._axis3don=False
	ax.patch.set_facecolor('black')
	ax.auto_scale_xyz([-1,1],[-1,1],[-1,1]) # true axis ratio
	incl=90-i
	ax.view_init(incl,30)  # view angle
	plt.axis('image')
	srtnum='{0:03}'.format(i)
	plt.savefig('film/anim_orbit'+srtnum+'.png', bbox_inches=0)
        plt.close()
