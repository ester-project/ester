# Plot of the Brunt-Vaisala frequency in a meridional plane
# with a zoom on the surface layers
# and a plot of the isotherm T=50000K where kappa-mechanism is active

from ester import *
from numpy import *
import matplotlib.pyplot as plt

a=star2d('../runs/Master/M3O7.h5')

# Isotherme à 50000K ----------------------------------------------
temp=a.Tc*a.T
T_min=-temp.max() # minimum autorisé pour le graphique
temp[temp <T_min]=T_min
z=c_[temp,temp[:,-2::-1],temp[:,1:],temp[:,-2::-1]]

n_depth=-250
th_ester=c_[a.th+pi/2,a.th[:,1:]] # from pi to 0
r_ester=c_[a.r[:,-1::-1],a.r[:,1:]] # we avoid double vue at pi/2
temp_n=c_[temp[:,-1::-1],temp[:,1:]] 
cm=plt.cm.autumn
cs=plt.contour(th_ester[n_depth:-1,:]/pi,r_ester[n_depth:-1,:]-r_ester[-1,:],temp_n[n_depth:-1,:],levels=[50000])

# retrieve data from the isocontour plot
num_levels = len(cs.allsegs)
num_element = len(cs.allsegs[0])  # in level 0
num_vertices = len(cs.allsegs[0][0])  # of element 0, in level 0
num_coord = len(cs.allsegs[0][0][0])


print('nb of pts in the isotherm =',num_vertices)
x_curve=[]
y_curve=[]
for i in range(num_vertices):
  x_curve.append(cs.allsegs[0][0][i][0])
  y_curve.append(cs.allsegs[0][0][i][1])
plt.close()
# Isotherm at 50000K -- Data are collected ------------------------

Req=a.Re
Rpol=a.Rp
print('Radii (Rsun) =',Req/R_SUN,Rpol/R_SUN)

x=a.r*sin(a.th)
y=a.r*cos(a.th)

# numpy.c_ = <numpy.lib.index_tricks.CClass object>
#    Translates slice objects to concatenation along the second axis.
x=c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
y=c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]

N2min=-a.N2.max() # minimum autorisé pour le graphique
a.N2[a.N2 <N2min]=N2min
bv=sign(a.N2)*sqrt(abs(a.N2))/2/pi*86400
z=c_[bv,bv[:,-2::-1],bv[:,1:],bv[:,-2::-1]]

fig=plt.figure(1,figsize=(7,7))
n=5
nxny=(n,n)
plt.subplot2grid(nxny, (0, 0), rowspan=n-1,colspan=n)

cm=plt.cm.RdBu_r
graph=plt.contourf(x,y,z,60,cmap=cm)
size=1.3

plt.axis('off')
plt.axis('equal')
plt.colorbar(graph,orientation='horizontal') 
plt.title('Brunt-Vaisala frequency')

plt.subplot2grid(nxny, (n-1, 0),colspan=n)
n_depth=-250
theta=c_[a.th+pi/2,a.th[:,1:]] # va de pi à 0
rnew=c_[a.r[:,-1::-1],a.r[:,1:]] # en évitant le doublon à pi/2
bvn=c_[bv[:,-1::-1],bv[:,1:]] # en évitant le doublon à pi/2
plt.contourf(theta[n_depth:-1,:]/pi,rnew[n_depth:-1,:]-rnew[-1,:],bvn[n_depth:-1,:],60,cmap=cm)
plt.plot(x_curve,y_curve,'k--') # isotherm at 50000K
plt.ylim(ymin=-0.035,ymax=0)
plt.xlabel(r'co-latitude (fraction of $\pi$)')
plt.ylabel('Depth')


plt.show()
plt.savefig('BV.png')
