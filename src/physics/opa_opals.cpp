#include<math.h>
#include"matrix.h"
#include"constants.h"
#include"physics.h"
#include"numdiff.h"
#include<string.h>
#include LAPACK

int opa_opals(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa) {

	static int init=0;
	static diff_gl glt,glr;
	static matrix opac,dopacT,dopacR;
	int i,j,n,N;
	int nt=50,nr=15;
	double a=0.01,b=2;

	if(!init) {
		char str[512];
		FILE *fp;
		matrix r,t,k,ki,Tr,Tt,A;
		
		sprintf(str,"%s/tables/opal/GN93hz",ESTER_ROOT);
		fp=fopen(str,"rt");
		while(strncmp(str,"TABLE # 73",10)) fgets(str,512,fp);
		for(i=0;i<3;i++) fgets(str,512,fp);
		r=zeros(1,19);
		fscanf(fp,"%*s");
		for(i=0;i<19;i++) fscanf(fp,"%lf",r.data()+i);
		for(i=0;i<2;i++) fgets(str,512,fp);
		k=9.999*ones(70,19);
		t=zeros(70,1);
		for(i=0;i<70;i++) {
			n=20;
			if(i>56) n=19;
			if(i>57) n=18;
			if(i>59) n=17;
			if(i>63) n=16;
			if(i>68) n=15;
			ki.read(1,n,fp);
			k.setblock(i,i,0,n-2,ki.block(0,0,1,n-1));
			t(i)=ki(0,0);
		}
		fclose(fp);
		
		glt.set_ndomains(1);
		glt.set_npts(nt);
		glt.set_xif(3.75,8.70);
		glt.init();
		glr.set_ndomains(1);
		glr.set_npts(nr);
		glr.set_xif(-8.,1.);
		glr.init();
	
		glr.eval(glr.x,r.transpose(),Tr);
		glt.eval(glt.x,t,Tt);
		Tr=Tr.transpose();
		
		n=0;
		for(i=0;i<70;i++) 
			for(j=0;j<19;j++) n+=(k(i,j)!=9.999);
		n+=glr.npts[0]*glt.npts[0];
		A=zeros(n,nr*nt);
		opac=zeros(n,1);
		n=0;
		for(j=0;j<19;j++) 
			for(i=0;i<70;i++) {
				if(k(i,j)==9.999) continue;
				ki=(Tt.row(i)*Tr.col(j)).transpose();
				ki.redim(1,ki.nrows()*ki.ncols());
				A.setrow(n,ki);
				opac(n)=k(i,j);
				n++;
			}
		for(j=0;j<glr.npts[0];j++) 
			for(i=0;i<glt.npts[0];i++) {
				ki=(glt.P.block(0).row(i)*glr.P.block(0).transpose().col(j)).transpose();
				ki.redim(1,ki.nrows()*ki.ncols());
				ki=ki*a*pow(10,b*((double) i/glt.npts[0]+(double) j/glr.npts[0]));
				ki=ki.transpose();
				ki.redim(1,ki.nrows()*ki.ncols());
				A.setrow(n,ki);
				opac(n)=0;
				n++;
			}
		
		char trans='N';
		int N,M,nrhs=1,info,lwork;
		double *work;
	
		N=A.nrows();
		M=glt.npts[0]*glr.npts[0];
		work=new double[1];lwork=-1;
		dgels_(&trans,&N,&M,&nrhs,A.data(),&N,opac.data(),&N,work,&lwork,&info);
		lwork=work[0];
		delete [] work;
		work=new double[lwork];
		dgels_(&trans,&N,&M,&nrhs,A.data(),&N,opac.data(),&N,work,&lwork,&info);
		delete [] work;
		opac=opac.block(0,M-1,0,0);
		opac.redim(glt.npts[0],glr.npts[0]);
		dopacT=(glt.D,opac);
		dopacR=(opac,glr.D.transpose());
		init=1;
	}

	int error=0;
	matrix logR,logT,logk,dlogkT,dlogkR,q,qT,qR;
	int nf,nc;
	
	nf=T.nrows();nc=T.ncols();
	n=T.nrows()*T.ncols();
	logT=log10(T);
	logT.redim(n,1);
	logR=log10(rho/pow(T/1e6,3));
	logR.redim(n,1);
	logk=zeros(n,1);
	dlogkT=zeros(n,1);
	dlogkR=zeros(n,1);
	q=glt.eval(opac,logT);
	qT=glt.eval(dopacT,logT);
	qR=glt.eval(dopacR,logT);
	q=q.transpose();
	qT=qT.transpose();
	qR=qR.transpose();
	for(i=0;i<n;i++) {
		logk(i)=glr.eval(q.col(i),logR(i))(0);
		dlogkT(i)=glr.eval(qT.col(i),logR(i))(0);
		dlogkR(i)=glr.eval(qR.col(i),logR(i))(0);
	}	
	logk.redim(nf,nc);
	dlogkT.redim(nf,nc);
	dlogkR.redim(nf,nc);
	
	opa.k=pow(10,logk);
	opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
	opa.dlnxi_lnrho=-1-dlogkR;
    opa.dlnxi_lnT=3-dlogkT+3*dlogkR;
	
	return error;
		
}
		
