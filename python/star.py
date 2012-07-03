from numpy import *
from matplotlib.pyplot import *
import os
import struct

rcParams['patch.antialiased']=False

class star2d:
	def __init__(self,file):
		names=['R/R_SUN','M/M_SUN','L/L_SUN','M','L','Xr','r','z','th','opa.k','w','p','rho','eos.G1','eos.del_ad','T','D','opa.xi','Teff','gsup','rex','map.R','R','nuc.eps','N2','eos.del_ad','Omega','phi','phiex','Dex','Dt','Dtodd','Dt2','Omega_bk','eos.cp','eos.G3_1','eos.cv','eos.d','eos.chi_T','eos.chi_rho','opa.dlnxi_lnrho','opa.dlnxi_lnT','nuc.pp','nuc.cno','vr','vt','stream','G','psi','virial','energy_test']		
		fp=open('/tmp/python_star_template','w')
		fp.write('\\conf{equator=1}\n\\conf{pole=1}\n\\conf{dim=1}\n')
		fp.write('${nr}${nth}${nex}${ndomains}${npts}${conv}')
		for x in names:
			fp.write('${'+x+'}')
		fp.close()
		status=os.system('gen_output '+file+' < /tmp/python_star_template > /tmp/python_star_out')
		if status:
			self.nr=0
			return
		fp=open('/tmp/python_star_out','rb')	
		self.nr=struct.unpack('i',fp.read(4))[0]
		self.nth=struct.unpack('i',fp.read(4))[0]+2
		self.nex=struct.unpack('i',fp.read(4))[0]
		self.ndom=struct.unpack('i',fp.read(4))[0]
		self.npts=fromfile(fp,'i',self.ndom)
		self.conv=struct.unpack('i',fp.read(4))[0]
		for x in names:
			x=x.replace('/','_')
			x=x.replace('.','_')
			if x=='eps' or x=='eps_c' or x=='M' or x=='L' or x=='Omega_bk' or x=='R' or x=='Omega' or x=='Omegac' or x=='R_R_SUN' or x=='Omega_bk' or x=='M_M_SUN' or x=='L_L_SUN' or x=='virial' or x=='energy_test':
				s="self."+x+"=struct.unpack('d',fp.read(8))[0]"
			elif x=='Teff' or x=='gsup':
				s="self."+x+"=fromfile(fp,'d',"+str(self.nth)+")"
			elif x=='th':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nth)+"),[1,"+str(self.nth)+"],'F')"
			elif x=='z':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr)+"),["+str(self.nr)+",1],'F')"
			elif x=='D':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr*self.nr)+"),["+str(self.nr)+","+str(self.nr)+"],'F')"
			elif x=='Dt' or x=='Dt2' or x=='Dtodd':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nth*self.nth)+"),["+str(self.nth)+","+str(self.nth)+"],'F')"
			elif x=='map_R':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.ndom*self.nth)+"),["+str(self.ndom)+","+str(self.nth)+"],'F')"
			elif x=='rex' or x=='phiex':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nex*self.nth)+"),["+str(self.nex)+","+str(self.nth)+"],'F')"
			elif x=='Dex':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nex*self.nex)+"),["+str(self.nex)+","+str(self.nex)+"],'F')"
			else:
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr*self.nth)+"),["+str(self.nr)+","+str(self.nth)+"],'F')"
			exec(s)
		fp.close()
		os.system('rm /tmp/python_star_template')
		os.system('rm /tmp/python_star_out')
		self.th=dot(ones((self.nr,1)),self.th)
		self.z=dot(self.z,ones((1,self.nth)))
		self.rz=dot(self.D,self.r)
		self.rzz=dot(self.D,self.rz)
		self.rt=dot(self.r,self.Dt)
		self.rtt=dot(self.r,self.Dt2)
		self.rzt=dot(self.rz,self.Dt)
		self.gzz=1./self.rz**2*(1.+self.rt**2/self.r**2)
		self.gzz[0,:]=1./self.rz[0,:]**2
		self.gzt=-self.rt/self.r**2/self.rz
		self.gzt[0,:]=0
		self.gtt=1./self.r**2
	
	def draw(self,z):
		x=self.r*sin(self.th)
		y=self.r*cos(self.th)
		x=c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
		y=c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
		z=c_[z,z[:,-2::-1],z[:,1:],z[:,-2::-1]]
		pcolor(x,y,z,rasterized=True)
		axis('image')
		
	def draw11(self,z):
		x=self.r*sin(self.th)
		y=self.r*cos(self.th)
		x=c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
		y=c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
		z=c_[z,-z[:,-2::-1],z[:,1:],-z[:,-2::-1]]
		pcolor(x,y,z)
		axis('image')
		
	def drawc(self,z,N=7,**kwargs):
		x=self.r*sin(self.th)
		y=self.r*cos(self.th)
		x=c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
		y=c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
		z=c_[z,z[:,-2::-1],z[:,1:],z[:,-2::-1]]
		hh=ishold()
		xs=self.r[-1]*sin(self.th[-1])
		ys=self.r[-1]*cos(self.th[-1])
		xs=r_[xs,-xs[-2::-1],-xs[1:],xs[-2::-1]]
		ys=r_[ys,ys[-2::-1],-ys[1:],-ys[-2::-1]]
		plot(xs,ys,'k')
		if not hh:
			hold(True)
		contour(x,y,z,N,**kwargs)
		axis('image')
		hold(hh)
	
	def drawc11(self,z,N=7,**kwargs):
		x=self.r*sin(self.th)
		y=self.r*cos(self.th)
		x=c_[x,-x[:,-2::-1],-x[:,1:],x[:,-2::-1]]
		y=c_[y,y[:,-2::-1],-y[:,1:],-y[:,-2::-1]]
		z=c_[z,-z[:,-2::-1],z[:,1:],-z[:,-2::-1]]
		hh=ishold()
		xs=self.r[-1]*sin(self.th[-1])
		ys=self.r[-1]*cos(self.th[-1])
		xs=r_[xs,-xs[-2::-1],-xs[1:],xs[-2::-1]]
		ys=r_[ys,ys[-2::-1],-ys[1:],-ys[-2::-1]]
		plot(xs,ys,'k')
		if not hh:
			hold(True)
		contour(x,y,z,N,**kwargs)
		axis('image')
		hold(hh)
		
		
		
	
	
		
