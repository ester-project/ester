from numpy import *
from matplotlib.pyplot import *
import os
import struct
import tempfile

rcParams['patch.antialiased']=False

SIG_SB=5.670400e-5;
K_BOL=1.3806503e-16;
HYDROGEN_MASS=1.67353249e-24;
A_RAD=7.565767e-15;
GRAV=6.67384e-8;
C_LIGHT=2.99792458e10;

M_SUN=1.9891e33;
R_SUN=6.95508e10;
L_SUN=3.8396e33;

class star2d:
	def __init__(self,file):
		names=['th','z','D','r','Dt','Dt2','xif','surff','Omega','Omega_bk','Omegac','X','Z','Xc',
				'rhoc','Tc','pc','M','R','Rp','Re','L','R/R_SUN','Rp/R_SUN','Re/R_SUN',
				'M/M_SUN','L/L_SUN','rex','phi','phiex','rho','p','T','w','G',
				'N2','opa.k','opa.xi','opa.dlnxi_lnT','opa.dlnxi_lnrho','eos.G1','eos.cp',
				'eos.del_ad','eos.G3_1','eos.cv','eos.prad','eos.chi_T','eos.chi_rho',
				'eos.d','nuc.eps','nuc.pp','nuc.cno','Teff','gsup','I','Dex','Dtodd',
				'It','map.R','vr','vt','virial','energy_test','eos.s','opa','eos',
				'Xr','Yr','Zr','X_H','X_He3','X_He4','X_C12','X_C13','X_N14','X_N15','X_O16',
				'X_O17','Mcore','Lz','Lzcore']
		fd,template_file=tempfile.mkstemp(prefix='star_template_',suffix='.tmp')
		fp=os.fdopen(fd,'w')
		fp.write('\\conf{equator=1}\n\\conf{pole=1}\n\\conf{dim=1}\n')
		fp.write('${nr}${nth}${nex}${ndomains}${npts}${conv}')
		for x in names:
			fp.write('${'+x+'}')
		fp.close()
		fd,out_file=tempfile.mkstemp(prefix='star_out_',suffix='.tmp')
		status=os.system('ester output '+file+' < '+template_file+' > '+out_file)
		if status:
			raise ValueError("Error reading file")
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
				 or x=='rhoc' or x=='pc' or x=='Tc' or x=='Rp' or x=='Re' or x=='Rp_R_SUN' or x=='Re_R_SUN' \
				 or x=='Mcore' or x=='Lz' or x=='Lzcore':
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
				s="self."+x+"=reshape(fromfile(fp,'d',"+str((self.ndomains+1)*self.nth)+"),["+str(self.ndomains+1)+","+str(self.nth)+"],'F')"
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
		
		
class star1d:
	def __init__(self,file):
		names=['D','r','xif','surff','X','Z','Xc',
				'rhoc','Tc','pc','M','R','L','R/R_SUN',
				'M/M_SUN','L/L_SUN','phi','rho','p','T',
				'N2','opa.k','opa.xi','opa.dlnxi_lnT','opa.dlnxi_lnrho','eos.G1','eos.cp',
				'eos.del_ad','eos.G3_1','eos.cv','eos.prad','eos.chi_T','eos.chi_rho',
				'eos.d','nuc.eps','nuc.pp','nuc.cno','Teff','gsup','I',
				'eos.s','opa','eos',
				'Xr','Yr','Zr','X_H','X_He3','X_He4','X_C12','X_C13','X_N14','X_N15','X_O16'
				,'X_O17','Mcore']
		fd,template_file=tempfile.mkstemp(prefix='star_template_',suffix='.tmp')
		fp=os.fdopen(fd,'w')
		fp.write('\\conf{dim=1}\n')
		fp.write('${nr}${ndomains}${npts}${conv}')
		for x in names:
			fp.write('${'+x+'}')
		fp.close()
		fd,out_file=tempfile.mkstemp(prefix='star_out_',suffix='.tmp')
		status=os.system('ester output '+file+' < '+template_file+' > '+out_file)
		if status:
			raise ValueError("Error reading file")
		fp=os.fdopen(fd,'rb')
		self.nr=struct.unpack('i',fp.read(4))[0]
		self.ndomains=struct.unpack('i',fp.read(4))[0]
		self.npts=fromfile(fp,'i',self.ndomains)
		self.conv=struct.unpack('i',fp.read(4))[0]
		self.nth=1
		for x in names:
			x=x.replace('/','_')
			x=x.replace('.','_')
			if x=='eps' or x=='eps_c' or x=='M' or x=='L' or x=='R' \
				 or x=='R_R_SUN' or x=='M_M_SUN' or x=='L_L_SUN' \
				 or x=='virial' or x=='energy_test' or x=='surff' or x=='X' or x=='Z' or x=='Xc' \
				 or x=='rhoc' or x=='pc' or x=='Tc' or x=='Mcore':
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
			elif x=='xif':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.ndomains+1)+"),["+str(self.ndomains+1)+",1],'F')"
			elif x=='D':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr*self.nr)+"),["+str(self.nr)+","+str(self.nr)+"],'F')"
			elif x=='map_R':
				s="self."+x+"=reshape(fromfile(fp,'d',"+str((self.ndomains+1)*self.nth)+"),["+str(self.ndomains+1)+","+str(self.nth)+"],'F')"
			else:
				s="self."+x+"=reshape(fromfile(fp,'d',"+str(self.nr*self.nth)+"),["+str(self.nr)+","+str(self.nth)+"],'F')"
			exec(s)
		fp.close()
		os.system('rm '+template_file)
		os.system('rm '+out_file)
		self.z=self.r
		self.rz=dot(self.D,self.r)
		self.rzz=dot(self.D,self.rz)
		self.rt=0*self.r
		self.rtt=0*self.r
		self.rzt=0*self.r
		self.gzz=1./self.rz**2*(1.+self.rt**2/self.r**2)
		self.gzz[0,:]=1./self.rz[0,:]**2
		self.gzt=-self.rt/self.r**2/self.rz
		self.gzt[0,:]=0
		self.gtt=1./self.r**2
		
def star_evol(file):
	
	n=0
	A=list();
	while(True):
		try:
			B=star2d(file+'_'+str(n).zfill(4));
		except:
			return A;
		print(file+'_'+str(n).zfill(4));
		A.append(B)
		n+=1;
		
	
	
	
	
	
	

