# Make a series of .png files viewing (here) the differential rotation
# (a.w field) from a series of ester models
# files written in "film/" directory
# preliminary step to make a film with other tools
# the script make_the_film, which uses mencoder, is a possible solution

import sys
from ester import *
import matplotlib.pyplot as plt
import math
import commands
import string
import numpy as np

Lsun=3.8396e33;
sigma=5.67e-5; # Stefan constant, cgs
R_SUN=6.95508e10;

liste=commands.getoutput('ls $HOME/uqbar/ESTER/Z0.014/M5_film/M5Xc1sol6565_film_0*')

lis=liste.split()
print('nb of images ',len(lis))

te=[0.0 for i in range(len(lis))]

i=0
for filename in lis:
	print(filename)
	a=star2d(filename)
	print('Xc',i,a.Xc)
	Req=a.Re/R_SUN
	Rpol=a.Rp/R_SUN
        print(Req,Rpol)
       
	x=a.r*np.sin(a.th)/R_SUN # en km
        y=a.r*np.cos(a.th)/R_SUN # en km
        x=np.c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
        y=np.c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
        pi=math.acos(-1)/3600.
        z=2*pi/np.c_[a.w,a.w[:,-2::-1],a.w[:,1:],a.w[:,-2::-1]]
        plt.figure(1,figsize=(7,7))
        plt.contourf(x,y,z,60,cmap=plt.cm.RdBu)
        plt.axis([-5.5,5.5,-5.5,5.5])
	srtnum='{0:03}'.format(i)
	plt.savefig('film/anim_M5_'+srtnum+'.png', bbox_inches=0)
        plt.close()
        i=i+1
