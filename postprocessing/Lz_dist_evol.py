# visualize the distribution of angular momentum among the domains
# and now its evolution through a series of models

import sys
from ester import *
import matplotlib.pyplot as plt
from numpy import *
import os
import subprocess as cmd



def eval_lz_domains(filename):

    a=star2d(filename)
    jfirst=zeros(a.ndomains,dtype=int)
    for i in arange(a.ndomains-1)+1:
           jfirst[i]=jfirst[i-1]+a.npts[i-1]
    
    Omega_scale=sqrt(a.pc/a.rhoc)/a.Rp
    
    Lz=(a.r*sin(a.th))**2*a.w*a.R**5*Omega_scale*a.rho*a.rhoc*a.rz*a.r**2
    # Integrateur radial: I(1,nr)
    Icheb=a.I.ravel() # 2D array with 1 column ==> 1D array
    
    # We first integrate over theta on every radial grid point
    dlz=[]
    for  i in range(a.nr):
          dlz.append(dot(a.It,Lz[i,:]))
    
    dlza=2*pi*array(dlz)
    # Then we integrate over 'r' on the first ncore grid points
    Lz_tot1=dot(dlza,Icheb)
    print('Lz total = %6.3e'%Lz_tot1)
    
    # examen Lz par domaine
    Lz_dom=[]
    r_dom=[]
    for i in range(a.ndomains):
        i1=jfirst[i]
        i2=jfirst[i]+a.npts[i]
        dlz=[]
        for j in range(i1,i2):
          dlz.append(dot(a.It,Lz[j,:]))
        dlza=2*pi*array(dlz)
    
        Icheb_dom=Icheb[i1:i2]
        Lz_dom.append(dot(dlza,Icheb_dom))
        r_dom.append((a.r[i1,0]+a.r[i2-1,0])/2)
    
    Lz_tot2=sum(Lz_dom)
    #print('Lz total bis = %6.3e'%Lz_tot2)
    return Lz_dom,r_dom

path0=os.path.expanduser( '~/Ester/')

ls=cmd.getoutput('ls -v *.h5_*').split()

lz=[]
r=[]
plt.figure(1,figsize=(12,12))
plt.rcParams.update({'font.size': 18})
for fichier in ls:
    #print(fichier)
    lf,rf= eval_lz_domains(fichier)
    lz.append(lf)
    r.append(rf)
    plt.plot(rf,lf)
    plt.yscale('log')
    

nfont=18
plt.xticks(fontsize=nfont)
plt.yticks(fontsize=nfont)
plt.xlabel('r/Rp',fontsize=nfont)
plt.ylabel('Lz / Domain',fontsize=nfont)
plt.title('Evolution of the angular momentum per domain')
plt.savefig('Lz_dom.png', bbox_inches=0)
    
plt.figure(2,figsize=(12,12))
nt=len(lz)
data=zeros((nt,12))
for it in range(nt):
   for idom in range(12):
      data[it,idom]=lz[it][idom] 

for idom in range(12):
   plt.plot(data[:,idom],label=str(idom))
plt.yscale('log')
plt.legend()
plt.show()
