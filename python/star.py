from numpy import *
from matplotlib.pyplot import *
import os
import struct
import tempfile

rcParams['patch.antialiased']=False

class star2d:
	def __init__(self,file):
		names=['th','z','D','r','Dt','Dt2','xif','surff','Omega','Omega_bk','Omegac','X','Z','Xc',
				'rhoc','Tc','pc','M','R','Rp','Re','L','R/R_SUN','Rp/R_SUN','Re/R_SUN',
				'M/M_SUN','L/L_SUN','rex','phi','phiex','rho','p','T','w','G','Xr',
				'N2','opa.k','opa.xi','opa.dlnxi_lnT','opa.dlnxi_lnrho','eos.G1','eos.cp',
				'eos.del_ad','eos.G3_1','eos.cv','eos.prad','eos.chi_T','eos.chi_rho',
				'eos.d','nuc.eps','nuc.pp','nuc.cno','Teff','gsup','I','Dex','Dt','Dtodd','Dt2',
				'It','map.R','Dtodd','vr','vt','virial','energy_test','eos.s','opa','eos']
		fd,template_file=tempfile.mkstemp(prefix='star_template_',suffix='.tmp')
		fp=os.fdopen(fd,'w')
		fp.write('\\conf{equator=1}\n\\conf{pole=1}\n\\conf{dim=1}\n')
		fp.write('${nr}${nth}${nex}${ndomains}${npts}${conv}')
		for x in names:
			fp.write('${'+x+'}')
		fp.close()
		fd,out_file=tempfile.mkstemp(prefix='star_out_',suffix='.tmp')
		status=os.system('gen_output '+file+' < '+template_file+' > '+out_file)
		if status:
			self.nr=0
			return
		fp=os.fdopen(fd,'rb')
		self.nr=struct.unpack('i',fp.read(4))[0]
		self.nth=struct.unpack('i',fp.read(4))[0]+2
		self.nex=struct.unpack('i',fp.read(4))[0]
		self.ndomains=struct.unpack('i',fp.read(4))[0]
		self.npts=fromfile(fp,'i',self.ndomains)
		self.conv=struct.unpack('i',fp.read(4))[0]
		for x in names:
			x=x.replace('/','_')
			x=x.replace('.','_')
			if x=='eps' or x=='eps_c' or x=='M' or x=='L' or x=='Omega_bk' or x=='R' or x=='Omega' \
				 or x=='Omegac' or x=='R_R_SUN' or x=='Omega_bk' or x=='M_M_SUN' or x=='L_L_SUN' \
				 or x=='virial' or x=='energy_test' or x=='surff' or x=='X' or x=='Z' or x=='Xc' \
				 or x=='rhoc' or x=='pc' or x=='Tc' or x=='Rp' or x=='Re' or x=='Rp_R_SUN' or x=='Re_R_SUN':
				s="self."+x+"=struct.unpack('d',fp.read(8))[0]"
			elif x=='opa' or x=='eos':
				c=list();
				cc=fp.read(1);
				while not cc=='\x00':
					c.append(cc);
					cc=fp.read(1);
				s=str();
				s=s.join(c);
				s="self."+x+"='"+s+"'"
			elif x=='Teff' or x=='gsup':
				s="self."+x+"=fromfile(fp,'d',"+str(self.nth)+")"
			elif x=='I':
				s="self."+x+"=fromfile(fp,'d',"+str(self.nr)+")"
			elif x=='It':
				s="self."+x+"=fromfile(fp,'d',"+str(self.nth)+")"
			elif x=='th':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nth)+"),[1,"+str(self.nth)+"],'F')"
			elif x=='z':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr)+"),["+str(self.nr)+",1],'F')"
			elif x=='xif':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.ndomains+1)+"),["+str(self.ndomains+1)+",1],'F')"
			elif x=='D':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr*self.nr)+"),["+str(self.nr)+","+str(self.nr)+"],'F')"
			elif x=='Dt' or x=='Dt2' or x=='Dtodd':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nth*self.nth)+"),["+str(self.nth)+","+str(self.nth)+"],'F')"
			elif x=='map_R':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.ndomains*self.nth)+"),["+str(self.ndomains)+","+str(self.nth)+"],'F')"
			elif x=='rex' or x=='phiex':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nex*self.nth)+"),["+str(self.nex)+","+str(self.nth)+"],'F')"
			elif x=='Dex':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nex*self.nex)+"),["+str(self.nex)+","+str(self.nex)+"],'F')"
			else:
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr*self.nth)+"),["+str(self.nr)+","+str(self.nth)+"],'F')"
			exec(s)
		fp.close()
		os.system('rm '+template_file)
		os.system('rm '+out_file)
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
		
		
		
	
	
		
