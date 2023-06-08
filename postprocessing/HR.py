# This script plot the HR tracks of a series of models using
# 1/ true luminosity and average effective temperature
# 2/ apparent polar luminosity and associated apparent effective temperature
# 3/ apparent equatorial luminosity and associated apparent effective temperature

from ester import *
import matplotlib.pyplot as plt
from numpy import *
import math
import subprocess as cmd
import os

path=os.path.expanduser( '~/Ester/')


sigma=SIG_SB
path=path+'local/runs/Evol/'
lpath=len(path)
liste=cmd.getoutput('ls '+path+'M5_ev_000*')

lis=liste.split()
n=len(lis)

te=zeros(n)
t_app_pole=zeros(n)
t_app_eq=zeros(n)
lum=zeros(n)
lum_app_pole=zeros(n)
lum_app_eq=zeros(n)


i=0
for filename in lis:
	print(filename[lpath:])
	a=star2d(filename)
	nth=a.nth+2
	integrand=a.r[a.nr-1,:]*sqrt(a.r[a.nr-1,:]**2+a.rt[a.nr-1,:]**2)
	integrand=a.R**2*2*pi*integrand
	surf=dot(integrand,a.It) # surface calculee
	te[i]=math.log10((a.L/surf/sigma)**(0.25))
	lum[i]=math.log10(a.L/L_SUN)
	print('Temp eff =',(a.L/surf/sigma)**0.25,'lum = ',a.L/L_SUN,' Req = ',a.Re/R_SUN)
	app_l= a.apparent_luminosity(0) # polar app. luminosity
	#print('app lum pole=',app_l)
	t_app_pole[i]=log10((L_SUN*app_l/surf/sigma)**(0.25))
	lum_app_pole[i]=log10(app_l)
	app_l= a.apparent_luminosity(90) # equat. app. luminosity
	#print('app lum eq=',app_l)
	t_app_eq[i]=math.log10((L_SUN*app_l/surf/sigma)**(0.25))
	lum_app_eq[i]=math.log10(app_l)
	i=i+1


xm=min(concatenate([te,t_app_eq]))-0.05
xmm=max(concatenate([te,t_app_eq]))+0.05
ym=min(concatenate([lum_app_eq,lum]))-0.05
ymm=max(concatenate([lum_app_pole,lum]))+0.05

plt.xlim(xmin=xmm,xmax=xm)      # set a userdefined x-range
plt.ylim(ymin= ym,ymax=ymm)      # set a userdefined y-range
plt.xlabel(r'$\log T_{\mathrm{eff}}$')     #labeling the plot
plt.ylabel(r'$\log\,L/L_\odot$')

plt.plot(te,lum,'b-')
plt.plot(t_app_pole,lum_app_pole,'r--')
plt.plot(t_app_eq,lum_app_eq,'g--')

plt.show()

