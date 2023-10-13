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

ls=cmd.getoutput('ls -v M10O75evol_*.h5').split()

lz=[]
r=[]
plt.figure(1,figsize=(12,12))
plt.rcParams.update({'font.size': 18})
ic=0
for fichier in ls:
    #print(fichier)
    lf,rf= eval_lz_domains(fichier)
    if ic==0:
       lf0=array(lf)
    lz.append(lf)
    r.append(rf)
    alz=array(lf)
    plt.plot(lf/lf0)
    #plt.yscale('log')
    ic+=1
    

nfont=18
tick_loc=arange(len(lf))
strr=[str(i+1) for i in tick_loc]
plt.xticks(ticks=tick_loc,labels=strr,fontsize=nfont)
plt.yticks(fontsize=nfont)
plt.xlabel('Domain number',fontsize=nfont)
plt.ylabel('Lz / Lz(t=0)',fontsize=nfont)
plt.title('Evolution of the angular momentum per domain')
plt.savefig('Lz_dom.png', bbox_inches=0)
    
plt.show()
