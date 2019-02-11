# This script plot the HR tracks of a series of models using
# 1/ true luminosity and average effective temperature
# 2/ apparent polar luminosity and associated apparent effective temperature
# 3/ apparent equatorial luminosity and associated apparent effective temperature

from ester import *
import matplotlib.pyplot as plt
import numpy as np
import math
import commands

plt.switch_backend('TKAgg')
sigma=SIG_SB
#liste=commands.getoutput('ls *ev_000*.h5')
liste=commands.getoutput('ls *ev_0*.h5')
MYRS=365.25*86400*1e6

lis=liste.split()

te=[0.0 for i in range(len(lis))]
lum=[0.0 for i in range(len(lis))]
age=[0.0 for i in range(len(lis))]

i=0
for filename in lis:
	#print filename
	a=star2d(filename)
	nth=a.nth
	#print 'Age = ',a.age,' Xc = ',a.X[0,0],' Tc = ',a.Tc,' rho_c=',a.rhoc,'opa(surf) = ',a.opacity[-1,0]
	lu=a.L/L_SUN
	te[i]=a.Teff[0,0]
	age[i]=a.age
	lum[i]=math.log10(a.L/L_SUN)
	if (i==0) : 
		E=10.**lum[i]*L_SUN*age[i]*MYRS
		#print 'Energy emitted (ergs) = %.3e'%E
	if (i>0) : 
		L0=10.**lum[i-1]*L_SUN
		L1=10.**lum[i]*L_SUN
		L=(L0+L1)/2
		E=E+L*(age[i]-age[i-1])*MYRS
		#print 'Energy emitted (ergs) = %.3e'%E
	print i,'Age = %.3f'%a.age,' Xc = %.2e'%a.X[0,0],' Tc = %.3e'%a.Tc,' rho_c = %.2f'%a.rhoc,'Lum = %.2f'%lu, \
                'E emitted = %.3e'%E
        i=i+1

print 'nombre de modeles =',i

xm=min(te)-1000
xmm=max(te)+1000
ym=min(lum)-0.05
ymm=max(lum)+0.05
plt.figure(1,figsize=(14,14))

plt.xlim(xmin=xmm,xmax=xm)      # set a userdefined x-range
plt.ylim(ymin= ym,ymax=ymm)      # set a userdefined y-range
plt.xticks(fontsize=18) #, rotation=90)
plt.yticks(fontsize=18) #, rotation=90)

plt.xlabel(r'$T_{\mathrm{eff}}$',fontsize=18)     #labeling the plot
plt.ylabel(r'$\log\,L/L_\odot$',fontsize=18)

plt.plot(te,lum,'b.')
print te[5],lum[5],age[5]

for k in range(i):
	if (k%10 == 0):
		text=plt.annotate(' age = %.2f'%age[k],(te[k],lum[k]),fontsize=18) #,xytext=(-2,1))
		dy=(lum[k+1]-lum[k])/(ymm-ym)
		dx=(te[k+1]-te[k])/(xmm-xm)
		angle=-180./np.pi*np.arctan(dy/dx)-90.
		#print 'angle=',angle
		text.set_rotation(angle)

plt.show()
