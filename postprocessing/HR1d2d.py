# This script plots the HR tracks of a series of models using
# 1/ true luminosity and average effective temperature
# 2/ apparent polar luminosity and associated apparent effective temperature
# 3/ apparent equatorial luminosity and associated apparent effective temperature
# and a 1D evolution model for reference

from ester import *
import matplotlib.pyplot as plt
from numpy import *
import subprocess as cmd
import os

path0=os.path.expanduser( '~/Ester/')


plt.close()
#plt.switch_backend('TKAgg')
sigma=SIG_SB
#liste=cmd.getoutput('ls *ev_000*.h5')
#---------------------- 1D models-----------------------------------
path=path0+'local/runs/Achernar/'
lpath=len(path)
liste=cmd.getoutput('ls '+path+'*_ev_*.h5')
MYRS=365.25*86400*1e6

lis=liste.split()

nf=len(lis)
te=zeros(nf)
lum=zeros(nf)
age=zeros(nf)

i=0

for filename in lis[0:10]:
	a=star1d(filename)
	lu=a.L/L_SUN
	te[i]=a.Teff[0,0]
	age[i]=float(i) # fake
	lum[i]=a.L/L_SUN
	if (i==0) : 
		E=lum[i]*L_SUN*age[i]*MYRS
		#print 'Energy emitted (ergs) = %.3e'%E
	if (i>0) : 
		L0=lum[i-1]*L_SUN
		L1=lum[i]*L_SUN
		L=(L0+L1)/2
		E=E+L*(age[i]-age[i-1])*MYRS
		#print 'Energy emitted (ergs) = %.3e'%E
	Esun=E/C_LIGHT**2/M_SUN
	print(filename[lpath:],'Age = %.3f'%age[i],' Xc = %.2e'%a.X[0,0],' Tc = %.3e'%a.Tc,\
         ' rho_c = %.2f'%a.rhoc,'Lum = %.2f'%lu) #,'E emitted = %.3e'%Esun
	i=i+1
i1d=i
print('nombre de modeles 1D =',i1d)

#---------------------- 2D models-----------------------------------
# we add 2D models
path=path0+'local/runs/olympe/Achernar2D_O4/'
lpath=len(path)
liste=cmd.getoutput('ls '+path+'Ach2D_ev_0*.h5')

lis=liste.split()
nf=len(lis)

te2D=zeros(nf)
t_app_eq=zeros(nf)
t_app_pole=zeros(nf)
lum2D=zeros(nf)
lum_app_pole=zeros(nf)
lum_app_eq=zeros(nf)
age2D=zeros(nf)

i=0
for filename in lis:
        a=star2d(filename)
        age2D[i]=a.age
        integrand=a.r[a.nr-1,:]*sqrt(a.r[a.nr-1,:]**2+a.rt[a.nr-1,:]**2)
        integrand=a.R**2*2*pi*integrand
        surf=dot(integrand,a.It) # surface calculee
        te2D[i]=(a.L/surf/sigma)**0.25
        lum2D[i]=a.L/L_SUN
        print(filename[lpath:],' Temp eff =',te2D[i],'lum = ',lum2D[i],' Req = ',a.Re/R_SUN)
        app_l= a.apparent_luminosity(0) # polar app. luminosity
        t_app_pole[i]=(L_SUN*app_l/surf/sigma)**0.25
        lum_app_pole[i]=app_l
        app_l= a.apparent_luminosity(90) # equat. app. luminosity
        t_app_eq[i]=(L_SUN*app_l/surf/sigma)**0.25
        lum_app_eq[i]=app_l
        i=i+1

xm=12000 #min(te)-1000
xmm=22000 #max(te)+1000
ym=900. #min(lum)-0.05
ymm=3000. #max(lum)+0.05
#plt.figure(1,figsize=(14,14))
fts=18

plt.xlim(xmin=xmm,xmax=xm)      # set a userdefined x-range
plt.ylim(ymin= ym,ymax=ymm)      # set a userdefined y-range
plt.xticks(1000*(12+2*arange(6)),('12','14','16','18','20','22'),fontsize=fts)
plt.yticks(1000*(1+arange(3)),('1000','2000','3000'),fontsize=fts)

plt.xlabel(r'$T_{\mathrm{eff}}$',fontsize=fts)     #labeling the plot
plt.ylabel(r'$L/L_\odot$',fontsize=fts)
plt.title('Evolution of a 6.22 M$_\odot$ with $\Omega=0$, $\Omega_0/\Omega_c=0.4$\
 and $\Omega_0/\Omega_c=0.5$',fontsize=fts)

plt.plot(te,lum,'k.')
plt.plot(te2D,lum2D,'ro')
plt.plot(t_app_pole,lum_app_pole,'r--')
plt.plot(t_app_eq,lum_app_eq,'r--')
#plt.yscale('log')

for k in range(i1d-1):
	if (age[k] < 50):
		kmod=5
	else:
		kmod=30
	if (k%kmod == 0):
		if age[k] < 50:
			ha='right'
			va='bottom'
		else:
			ha='left'
			va='top'
		if age[k] < 50 or age[k] > 53.9:
			text=plt.annotate(' age = %.2f'%age[k],(te[k],lum[k]),fontsize=fts, \
	           horizontalalignment=ha, verticalalignment=va) #,xytext=(-2,1))
			dy=(lum[k+1]-lum[k])/(ymm-ym)
			dx=(te[k+1]-te[k])/(xmm-xm)
			angle=-180./pi*arctan(dy/dx)-90.
			#print 'angle=',angle
			text.set_rotation(angle)
for k in range(i-1):
	if (age2D[k] < 50):
		kmod=5
	else:
		kmod=30
	if (k%kmod == 0):
		if age2D[k] < 50:
			ha='left'
			va='top'
		else:
			ha='left'
			va='top'
		#if age2D[k] < 50 or age2D[k] > 53.9:
		if age2D[k] > 48:
			text=plt.annotate(' age = %.2f'%age2D[k],(te2D[k],lum2D[k]),fontsize=fts, \
	           horizontalalignment=ha, verticalalignment=va) #,xytext=(-2,1))
			dy=(lum2D[k+1]-lum2D[k])/(ymm-ym)
			dx=(te2D[k+1]-te2D[k])/(xmm-xm)
			angle=-180./pi*arctan(dy/dx)-90.
			text.set_rotation(angle)

# Let's add another 2D track
path=path0+'local/runs/olympe/Achernar2D_O5/'
lpath=len(path)
liste=cmd.getoutput('ls '+path+'Ach2DO5_ev_0*.h5')

lis=liste.split()
nf=len(lis)

te2D=zeros(nf)
t_app_pole=zeros(nf)
t_app_eq=zeros(nf)
lum2D=zeros(nf)
lum_app_pole=zeros(nf)
lum_app_eq=zeros(nf)
age2D=zeros(nf)

i=0
for filename in lis:
        a=star2d(filename)
        age2D[i]=a.age
        integrand=a.r[a.nr-1,:]*sqrt(a.r[a.nr-1,:]**2+a.rt[a.nr-1,:]**2)
        integrand=a.R**2*2*pi*integrand
        surf=dot(integrand,a.It) # surface calculee
        te2D[i]=(a.L/surf/sigma)**0.25
        lum2D[i]=a.L/L_SUN
        print(filename[lpath:],' Temp eff =',te2D[i],'lum = ',lum2D[i],' Req = ',a.Re/R_SUN)
        app_l= a.apparent_luminosity(0) # polar app. luminosity
        t_app_pole[i]=(L_SUN*app_l/surf/sigma)**0.25
        lum_app_pole[i]=app_l
        app_l= a.apparent_luminosity(90) # equat. app. luminosity
        t_app_eq[i]=(L_SUN*app_l/surf/sigma)**0.25
        lum_app_eq[i]=app_l
        i=i+1

print(i-1)
plt.plot(te2D,lum2D,'go')
plt.plot(t_app_pole,lum_app_pole,'g--')
plt.plot(t_app_eq,lum_app_eq,'g--')


for k in range(i-1):
	if (age2D[k] < 50):
		kmod=i-2
	else:
		kmod=30
	if (k%kmod == 0):
		if age2D[k] < 50:
			ha='left'
			va='top'
		else:
			ha='left'
			va='top'
		if age2D[k] < 50 or age2D[k] > 53.9:
			text=plt.annotate(' age = %.2f'%age2D[k],(te2D[k],lum2D[k]),fontsize=fts, \
	           horizontalalignment=ha, verticalalignment=va)
			dy=(lum2D[k+1]-lum2D[k])/(ymm-ym)
			dx=(te2D[k+1]-te2D[k])/(xmm-xm)
			angle=-180./pi*arctan(dy/dx)-90.
			text.set_rotation(angle)

plt.show()

