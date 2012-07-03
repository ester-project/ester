#include"symbolic.h"
#include<string.h>
#include<math.h>
#include<stdlib.h>

void symbolic::create(int n_sc,int n_vec,int max_der) {

	int i,j,k,n,l,m;
	
	maxder=max_der;
	nsc=n_sc;
	nvec=n_vec;
	jnum=0;
	j=1;
	jr=new int*[maxder+2];
	for(i=0;i<maxder+2;i++) jr[i]=new int[maxder+2-i];
	for(i=0;i<maxder+2;i++) 
		for(k=0;k<=i;k++) jr[i-k][k]=j++;
	jsint=j++;
	jcost=j++;
	jsc=new int***[nsc];
	for(n=0;n<nsc;n++) {
		jsc[n]=new int**[maxder+1];
		for(i=0;i<maxder+1;i++) {
			jsc[n][i]=new int*[maxder+1-i];
			for(k=0;k<maxder+1-i;k++) jsc[n][i][k]=new int[maxder+1-i-k];
		}
		for(i=0;i<maxder+1;i++) 
			for(k=0;k<=i;k++)
				for(l=0;l<=k;l++) jsc[n][i-k][k-l][l]=j++;
	}
	jvec=new int****[nvec];
	for(n=0;n<nvec;n++) {
		jvec[n]=new int***[3];
		for(m=0;m<3;m++) {
			jvec[n][m]=new int**[maxder+1];
			for(i=0;i<maxder+1;i++) {
				jvec[n][m][i]=new int*[maxder+1-i];
				for(k=0;k<maxder+1-i;k++) jvec[n][m][i][k]=new int[maxder+1-i-k];
			}
			for(i=0;i<maxder+1;i++) 
				for(k=0;k<=i;k++)
					for(l=0;l<=k;l++) jvec[n][m][i-k][k-l][l]=j++;
		}
	}
	jvec_=new int****[nvec];
	for(n=0;n<nvec;n++) {
		jvec_[n]=new int***[3];
		for(m=0;m<3;m++) {
			jvec_[n][m]=new int**[maxder+1];
			for(i=0;i<maxder+1;i++) {
				jvec_[n][m][i]=new int*[maxder+1-i];
				for(k=0;k<maxder+1-i;k++) jvec_[n][m][i][k]=new int[maxder+1-i-k];
			}
			for(i=0;i<maxder+1;i++) 
				for(k=0;k<=i;k++)
					for(l=0;l<=k;l++) jvec_[n][m][i-k][k-l][l]=j++;
		}
	}	
	char str[20];
	scname=new char*[nsc];
	for(i=0;i<nsc;i++) {
		scname[i]=new char[20];
		sprintf(str,"a%d",i);
		strcpy(scname[i],str);
	}
	vecname=new char*[nvec];
	for(i=0;i<nvec;i++) {
		vecname[i]=new char[20];
		sprintf(str,"V%d",i);
		strcpy(vecname[i],str);
	}
	sccodename=new char*[nsc];
	for(i=0;i<nsc;i++) {
		sccodename[i]=new char[20];
		sprintf(str,"a%d",i);
		strcpy(sccodename[i],str);
	}
	veccodename=new char*[nvec];
	for(i=0;i<nvec;i++) {
		veccodename[i]=new char[20];
		sprintf(str,"V%d",i);
		strcpy(veccodename[i],str);
	}
	scparity=new int*[nsc];
	for(i=0;i<nsc;i++) {
		scparity[i]=new int[2];
		scparity[i][0]=0;scparity[i][1]=0;
	}
	vecparity=new int**[nvec];
	for(i=0;i<nvec;i++) {
		vecparity[i]=new int*[3];
		for(k=0;k<3;k++) vecparity[i][k]=new int[2];
		vecparity[i][0][0]=0;vecparity[i][0][1]=0;
		vecparity[i][1][0]=1;vecparity[i][1][1]=1;
		vecparity[i][2][0]=0;vecparity[i][2][1]=0;
	}
	N=j;
	initd=1;

}

void symbolic::destroy() {

	if(!initd) return;
	
	int i,n,k,m;
	
	for(i=0;i<maxder+2;i++) delete [] jr[i];
	delete [] jr;

	for(n=0;n<nsc;n++) {
		for(i=0;i<maxder+1;i++) {
			for(k=0;k<maxder+1-i;k++) delete [] jsc[n][i][k];
			delete [] jsc[n][i];
		}
		delete [] jsc[n];
	}
	delete [] jsc;
	
	for(n=0;n<nvec;n++) {
		for(m=0;m<3;m++) {
			for(i=0;i<maxder+1;i++) {
				for(k=0;k<maxder+1-i;k++) delete [] jvec[n][m][i][k];
				delete [] jvec[n][m][i];
			}
			delete [] jvec[n][m];
		}
		delete [] jvec[n];
	}
	delete [] jvec;
	
	for(n=0;n<nvec;n++) {
		for(m=0;m<3;m++) {
			for(i=0;i<maxder+1;i++) {
				for(k=0;k<maxder+1-i;k++) delete [] jvec_[n][m][i][k];
				delete [] jvec_[n][m][i];
			}
			delete [] jvec_[n][m];
		}
		delete [] jvec_[n];
	}
	delete [] jvec_;
	
	for(n=0;n<nsc;n++) delete [] scname[n];
	delete [] scname;
	for(n=0;n<nvec;n++) delete [] vecname[n];
	delete [] vecname;
	for(n=0;n<nsc;n++) delete [] sccodename[n];
	delete [] sccodename;
	for(n=0;n<nvec;n++) delete [] veccodename[n];
	delete [] veccodename;
	for(n=0;n<nsc;n++) delete [] scparity[n];
	delete [] scparity;
	for(n=0;n<nvec;n++) {
		for(m=0;m<3;m++) delete [] vecparity[n][m];
		delete [] vecparity[n];
	}
	delete [] vecparity;
	
	initd=0;
	
}

void symbolic::set_scname(int isc,const char *name) {

	delete [] scname[isc];
	scname[isc]=new char[strlen(name)+1];
	strcpy(scname[isc],name);
}

void symbolic::set_vecname(int ivec,const char *name) {

	delete [] vecname[ivec];
	vecname[ivec]=new char[strlen(name)+1];
	strcpy(vecname[ivec],name);
}

void symbolic::set_sccodename(int isc,const char *name) {

	delete [] sccodename[isc];
	sccodename[isc]=new char[strlen(name)+1];
	strcpy(sccodename[isc],name);
}

void symbolic::set_veccodename(int ivec,const char *name) {

	delete [] veccodename[ivec];
	veccodename[ivec]=new char[strlen(name)+1];
	strcpy(veccodename[ivec],name);
}

void symbolic::set_scparity(int isc,int par0,int par1) {

	scparity[isc][0]=par0;
	scparity[isc][1]=par1;
}

void symbolic::set_vecparity(int ivec,int icomp,int par0,int par1) {

	vecparity[ivec][icomp][0]=par0;
	vecparity[ivec][icomp][1]=par1;
}

void symbolic::init(matrix &g,int n) {

	g=zeros(N,n);

}

matrix symbolic::g(int i,int j) {

	matrix gij;

	if(i==0&&j==0) {
		init(gij,2);
		gij(jnum,0)=1;
		gij(jr[1][0],0)=-2;
		gij(jnum,1)=1;
		gij(jr[0][0],1)=-2;
		gij(jr[1][0],1)=-2;
		gij(jr[0][1],1)=2;
	} else if(i==1&&j==0||i==0&&j==1) {
		init(gij,1);
		gij(jnum,0)=-1;
		gij(jr[0][1],0)=1;
		gij(jr[0][0],0)=-2;
		gij(jr[1][0],0)=-1;
	} else if(i==1&&j==1) {
		init(gij,1);
		gij(jnum,0)=1;
		gij(jr[0][0],0)=-2;
	} else if(i==2&&j==2) {
		init(gij,1);
		gij(jnum,0)=1;
		gij(jr[0][0],0)=-2;
		gij(jsint,0)=-2;
	} else {
		init(gij,1);
		gij(jnum,0)=0;
	}
	
	return gij;
}

matrix symbolic::g_(int i,int j) {

	matrix gij;

	if(i==0&&j==0) {
		init(gij,1);
		gij(jnum,0)=1;
		gij(jr[1][0],0)=2;
	} else if(i==1&&j==0||i==0&&j==1) {
		init(gij,1);
		gij(jnum,0)=1;
		gij(jr[0][1],0)=1;
		gij(jr[1][0],0)=1;
	} else if(i==1&&j==1) {
		init(gij,2);
		gij(jnum,0)=1;
		gij(jr[0][0],0)=2;
		gij(jnum,1)=1;
		gij(jr[0][1],1)=2;
	} else if(i==2&&j==2) {
		init(gij,1);
		gij(jnum,0)=1;
		gij(jr[0][0],0)=2;
		gij(jsint,0)=2;
	} else {
		init(gij,1);
		gij(jnum,0)=0;
	}
	
	return gij;
}

matrix symbolic::r(int exp) {

	matrix A;
	init(A,1);
	
	A(jnum,0)=1;
	A(jr[0][0],0)=exp;
	return A;
}

matrix symbolic::sint(int exp) {

	matrix A;
	init(A,1);
	
	A(jnum,0)=1;
	A(jsint,0)=exp;
	return A;
}

matrix symbolic::cost(int exp) {

	matrix A;
	init(A,1);
	
	A(jnum,0)=1;
	A(jcost,0)=exp;
	return A;
}

matrix symbolic::sc(int isc,int exp) {

	matrix A;
	init(A,1);
	
	A(jnum,0)=1;
	A(jsc[isc][0][0][0],0)=exp;
	return A;
}

matrix symbolic::vec(int ivec,int icomp) {

	matrix A;
	init(A,1);
	
	A(jnum,0)=1;
	A(jvec[ivec][icomp][0][0][0],0)=1;
	return A;
}

matrix symbolic::vec_(int ivec,int icomp) {

	matrix A;
	init(A,1);
	
	A(jnum,0)=1;
	A(jvec_[ivec][icomp][0][0][0],0)=1;
	return A;
}

matrix symbolic::perm(int i,int j,int k) {
	
	matrix e;

	init(e,1);

	if(i==0&&j==1&&k==2||i==1&&j==2&&k==0||i==2&&j==0&&k==1) e(jnum,0)=1;
	else if(i==2&&j==1&&k==0||i==0&&j==2&&k==1||i==1&&j==0&&k==2) e(jnum,0)=-1;
	else e(jnum,0)=0;
	
	e(jr[0][0],0)=-2;
	e(jr[1][0],0)=-1;
	e(jsint,0)=-1;
	
	return e;
}

matrix symbolic::perm_(int i,int j,int k) {
	
	matrix e;

	init(e,1);

	if(i==0&&j==1&&j==2||i==1&&j==2&&j==0||i==2&&j==0&&j==1) e(jnum,0)=1;
	else if(i==2&&j==1&&j==0||i==0&&j==2&&j==1||i==1&&j==0&&j==2) e(jnum,0)=-1;
	else e(jnum,0)=0;
	
	e(jr[0][0],0)=2;
	e(jr[1][0],0)=1;
	e(jsint,0)=1;
	
	return e;
}

matrix symbolic::simplify(matrix a) {

	matrix ind,b;
	int i,j,k;

	b=zeros(a.nrows(),a.ncols());
	ind=zeros(1,a.ncols());
	k=0;
	for(i=0;i<a.ncols();i++) {
		if(ind(i)) continue;
		for(j=i+1;j<a.ncols();j++) {
			if(isequal(a.block(1,N-1,i,i),a.block(1,N-1,j,j))) {
				a(jnum,i)+=a(jnum,j);
				ind(j)=1;
			}
		}
		if(fabs(a(jnum,i))>1e-12) {
			b.setcol(k,a.col(i));
			k++;
		}
		ind(i)=1;
	}
	
	if(k) b=b.block(0,N-1,0,k-1);
	else init(b,1);

	return b;

}

matrix symbolic::simplify_trig(matrix a) {

	int mod=1,i;

	a=simplify(a);
	while(mod) {
		mod=0;
		mod=simplify_trig_1(a);
		a=simplify(a);
		mod+=simplify_trig_2(a);
		a=simplify(a);
	}
	return a;
	
}

int symbolic::simplify_trig_1(matrix &a) {

	matrix ind,b,c;
	int i,j,mod=0;
	
	for(i=0;i<a.ncols();i++) {
		c=prod(a.col(i),cost(-2));
		for(j=i+1;j<a.ncols();j++) {
			if(isequal(c.block(1,N-1,0,0),prod(a.col(j),sint(-2)).block(1,N-1,0,0))) {
				a.setcol(i,c);
				a(jnum,j)-=a(jnum,i);
				mod=1;
				break;
			}
		}
	}

	return mod;

}

int symbolic::simplify_trig_2(matrix &a) {

	matrix ind,b,c;
	int i,j,mod=0;
	
	for(i=0;i<a.ncols();i++) {
		c=prod(a.col(i),sint(-2));
		for(j=i+1;j<a.ncols();j++) {
			if(isequal(c.block(1,N-1,0,0),prod(a.col(j),cost(-2)).block(1,N-1,0,0))) {
				a.setcol(i,c);
				a(jnum,j)-=a(jnum,i);
				mod=1;
				break;
			}
		}
	}

	return mod;

}

matrix symbolic::add(matrix a,matrix b) {

	matrix c;
	init(c,a.ncols()+b.ncols());

	c.setblock(0,N-1,0,a.ncols()-1,a);
	c.setblock(0,N-1,a.ncols(),a.ncols()+b.ncols()-1,b);

	return c;

}

matrix symbolic::neg(matrix a) {

	matrix c(a);

	c.setblock(jnum,jnum,0,a.ncols()-1,-a.block(jnum,jnum,0,a.ncols()-1));

	return c;

}

matrix symbolic::prod(matrix a,matrix b) {

	matrix c;
	int i,j;

	init(c,a.ncols()*b.ncols());

	for(i=0;i<a.ncols();i++) {
		for(j=0;j<b.ncols();j++) {
			c.setcol(i*b.ncols()+j,a.col(i)+b.col(j));
			c(jnum,i*b.ncols()+j)=a(jnum,i)*b(jnum,j);
		}
	}

	return c;

}

void symbolic::write2(const char* var,int ex) {

	printf("%s",var);
	if(ex!=1) printf("^{%d}",ex);
}

void symbolic::write3(const char* var,int ex) {

	if(ex!=1) printf("\\left(");
	printf("%s",var);
	if(ex!=1) printf("\\right)^{%d}",ex);
}

void symbolic::write_frac(double f) {

	int n,a=0,b;
	
	for(n=1;n<1000;n++) {
		if(fabs(round(f*n)-f*n)<1e-8) {
			a=round(f*n);
			b=n;
			break;
		}
	}
	if(a==0&&f!=0) printf("%8.8f",f);
	else if(b==1) printf("%d",a);
	else {
		if(a<0) printf("-");
		printf("\\frac{%d}{%d}",abs(a),b);
	}

}

void symbolic::write(matrix a) {

	int i,j,k,l,frac,n,m;
	char str[4096];

	for(i=0;i<a.ncols();i++) {
		frac=0;
		if(!exist(a.block(1,N-1,i,i)!=0)) {
			if(i>0&&a(jnum,i)>=0) printf("+");
			write_frac(a(jnum,i));
			continue;
		}
		if(a(jnum,i)>=0&&i>0) printf("+");
		if(a(jnum,i)!=1&&a(jnum,i)!=-1) write_frac(a(jnum,i));
		if(a(jnum,i)==-1) printf("-");
		if(exist(a.block(1,N-1,i,i)<0)) frac=1;
		if(frac) printf("\\frac{");
		for(j=0;j<=frac;j++) {
			if(!exist(round(a.block(jr[0][0],jcost,i,i))>0)&&!j&&frac) printf("1");
			else if(!exist(round(a.block(1,N-1,i,i))>0)&&!j) printf("1");
			for(k=0;k<maxder+2;k++) {
				for(l=0;l<=k;l++) {
					strcpy(str,"r");
					if(k>0) strcat(str,"_{");
					for(n=0;n<k-l;n++) strcat(str,"\\zeta");
					for(n=0;n<l;n++) strcat(str,"\\theta");
					if(k>0) strcat(str,"}");
					if(a(jr[k-l][l],i)>0&&!j||a(jr[k-l][l],i)<0&&j) write2(str,int(fabs(a(jr[k-l][l],i))));
				}
			}
			if(a(jsint,i)>0&&!j||a(jsint,i)<0&&j) {write2("\\sin",int(fabs(a(jsint,i))));printf("\\theta ");}
			if(a(jcost,i)>0&&!j||a(jcost,i)<0&&j) {write2("\\cos",int(fabs(a(jcost,i))));printf("\\theta ");}
			if(frac&&!j) printf("}{");
		}
		if(frac) printf("}");
		for(n=0;n<nsc;n++) {
			for(j=0;j<maxder+1;j++) 
				for(k=0;k<=j;k++)
					for(l=0;l<=k;l++) {
						if(!j) {
							sprintf(str," %s",scname[n]);
							if(a(jsc[n][j-k][k-l][l],i)!=0) write2(str,a(jsc[n][j-k][k-l][l],i));
						} else {
							strcpy(str,"\\frac{\\partial");
							if(j>1) sprintf(str,"%s^%d",str,j);
							sprintf(str,"%s %s}{",str,scname[n]);
							if(j-k) strcat(str,"\\partial\\zeta");
							if(j-k>1) sprintf(str,"%s^%d",str,j-k);
							if(k-l) strcat(str,"\\partial\\theta");
							if(k-l>1) sprintf(str,"%s^%d",str,k-l);
							if(l) strcat(str,"\\partial\\varphi");
							if(l>1) sprintf(str,"%s^%d",str,l);
							strcat(str,"}");
							if(a(jsc[n][j-k][k-l][l],i)!=0) write3(str,a(jsc[n][j-k][k-l][l],i));
						}
					}
		}
		for(n=0;n<nvec;n++) {
			for(m=0;m<3;m++) {
				for(j=0;j<maxder+1;j++) 
					for(k=0;k<=j;k++)
						for(l=0;l<=k;l++) {
							if(!j) {
								sprintf(str," %s",vecname[n]);
								if(m==0) strcat(str,"^\\zeta");
								if(m==1) strcat(str,"^\\theta");
								if(m==2) strcat(str,"^\\varphi");
							} else {
								strcpy(str,"\\frac{\\partial");
								if(j>1) sprintf(str,"%s^%d",str,j);
								sprintf(str,"%s %s",str,vecname[n]);
								if(m==0) strcat(str,"^\\zeta");
								if(m==1) strcat(str,"^\\theta");
								if(m==2) strcat(str,"^\\varphi");
								strcat(str,"}{");
								if(j-k) strcat(str,"\\partial\\zeta");
								if(j-k>1) sprintf(str,"%s^%d",str,j-k);
								if(k-l) strcat(str,"\\partial\\theta");
								if(k-l>1) sprintf(str,"%s^%d",str,k-l);
								if(l) strcat(str,"\\partial\\varphi");
								if(l>1) sprintf(str,"%s^%d",str,l);
								strcat(str,"}");
							}
							if(a(jvec[n][m][j-k][k-l][l],i)!=0) write3(str,a(jvec[n][m][j-k][k-l][l],i));
						}
			}
		}
		for(n=0;n<nvec;n++) {
			for(m=0;m<3;m++) {
				for(j=0;j<maxder+1;j++) 
					for(k=0;k<=j;k++)
						for(l=0;l<=k;l++) {
							if(!j) {
								sprintf(str," %s",vecname[n]);
								if(m==0) strcat(str,"_\\zeta");
								if(m==1) strcat(str,"_\\theta");
								if(m==2) strcat(str,"_\\varphi");
							} else {
								strcpy(str,"\\frac{\\partial");
								if(j>1) sprintf(str,"%s^%d",str,j);
								sprintf(str,"%s %s",str,vecname[n]);
								if(m==0) strcat(str,"_\\zeta");
								if(m==1) strcat(str,"_\\theta");
								if(m==2) strcat(str,"_\\varphi");
								strcat(str,"}{");
								if(j-k) strcat(str,"\\partial\\zeta");
								if(j-k>1) sprintf(str,"%s^%d",str,j-k);
								if(k-l) strcat(str,"\\partial\\theta");
								if(k-l>1) sprintf(str,"%s^%d",str,k-l);
								if(l) strcat(str,"\\partial\\varphi");
								if(l>1) sprintf(str,"%s^%d",str,l);
								strcat(str,"}");
							}
							if(a(jvec_[n][m][j-k][k-l][l],i)!=0) write3(str,a(jvec_[n][m][j-k][k-l][l],i));
						}
			}
		}		
	}
}

void symbolic::writecode2(char *str,int derz,int dert,int par0,int par1) {
	
	char str2[512];
	int i,dt2,dt;
	
	if(!derz&&!dert) return;
	strcpy(str2,str);
	strcpy(str,"(");
	for(i=0;i<derz;i++) strcat(str,"D,");
	strcat(str,str2);
	dt2=dert/2;
	dt=dert%2;
	for(i=0;i<dt2;i++) {
		if(!par0&&!par1) strcat(str,",Dt2");
		else sprintf(str,"%s,map.leg.D2_%d%d",str,par0,par1);
	}
	for(i=0;i<dt;i++) {
		if(!par0&&!par1) strcat(str,",Dt");
		else sprintf(str,"%s,map.leg.D_%d%d",str,par0,par1);
	}
	strcat(str,")");
	
}

void symbolic::writecode3(char *str,int n,int *start) {

	char str2[512];
	int i;

	if(n==0) {str[0]='\0';return;}
	strcpy(str2,str);
	if(*start) {
		if(n>0) sprintf(str,"*%s",str2);
		else sprintf(str,"/%s",str2);
	} else {
		if(n>0) sprintf(str,"%s",str2);
		else sprintf(str,"1./%s",str2);
	}
	*start=1;
	for(i=1;i<abs(n);i++) {
		if(n>0) sprintf(str,"%s*%s",str,str2);
		else sprintf(str,"%s/%s",str,str2);
	}

}

void symbolic::writecode_frac(double f) {

	int n,a=0,b;
	
	for(n=1;n<1000;n++) {
		if(fabs(round(f*n)-f*n)<1e-8) {
			a=round(f*n);
			b=n;
			break;
		}
	}
	if(a==0&&f!=0) printf("%e",f);
	else if(b==1) printf("%d.",a);
	else {
		printf("%d./%d.",abs(a),b);
	}

}

void symbolic::writecode(matrix a) {

	int i,j,k,l,n,m,start;
	char str[1024];

	for(i=0;i<a.ncols();i++) {
		start=0;
		if(a(jnum,i)<0) printf("-");
		else printf("+");
		if(a(jnum,i)!=1||!exist(a.block(1,N-1,i,i)!=0)) {
			writecode_frac(fabs(a(jnum,i)));
			start=1;
		}
		for(k=0;k<maxder+2;k++) {
			for(l=0;l<=k;l++) {
				if(a(jr[k-l][l],i)!=0) {
					if(k==0) strcpy(str,"r");
					else if(k==1) {
						if(l==0) strcpy(str,"map.rz");
						else strcpy(str,"map.rt");
					} 
					else if(k==2) {
						if(l==0) strcpy(str,"map.rzz");
						else if(l==1) strcpy(str,"map.rzt");
						else strcpy(str,"map.rtt");
					}
					else if(k-l>=2) {
						strcpy(str,"map.rzz");
						writecode2(str,k-l-2,l,0,0);
					}
					else {
						strcpy(str,"map.rtt");
						writecode2(str,k-l,l-2,0,0);
					}
					writecode3(str,a(jr[k-l][l],i),&start);
					printf("%s",str);
				}
			}
		}
		if(a(jsint,i)!=0) {
			strcpy(str,"sin(th)");
			writecode3(str,a(jsint,i),&start);
			printf("%s",str);
		}
		if(a(jcost,i)!=0) {
			strcpy(str,"cos(th)");
			writecode3(str,a(jcost,i),&start);
			printf("%s",str);
		}
		for(n=0;n<nsc;n++) {
			for(j=0;j<maxder+1;j++) 
				for(k=0;k<=j;k++) {
					l=0;
					if(a(jsc[n][j-k][k-l][l],i)!=0) {
						strcpy(str,sccodename[n]);
						writecode2(str,j-k,k,scparity[n][0],scparity[n][1]);
						writecode3(str,a(jsc[n][j-k][k-l][l],i),&start);
						printf("%s",str);
					}
				}
		}
		for(n=0;n<nvec;n++) {
			for(m=0;m<3;m++) {
				for(j=0;j<maxder+1;j++) 
					for(k=0;k<=j;k++) {
						l=0;
						if(a(jvec[n][m][j-k][k-l][l],i)!=0) {
							sprintf(str,"%s[%d]",veccodename[n],m);
							writecode2(str,j-k,k,vecparity[n][m][0],vecparity[n][m][1]);
							writecode3(str,a(jvec[n][m][j-k][k-l][l],i),&start);
							printf("%s",str);
						}
					}
			}
		}
		for(n=0;n<nvec;n++) {
			for(m=0;m<3;m++) {
				for(j=0;j<maxder+1;j++) 
					for(k=0;k<=j;k++) {
						l=0;
						if(a(jvec_[n][m][j-k][k-l][l],i)!=0) {
							sprintf(str,"%s_[%d]",veccodename[n],m);
							writecode2(str,j-k,k,vecparity[n][m][0],vecparity[n][m][1]);
							writecode3(str,a(jvec_[n][m][j-k][k-l][l],i),&start);
							printf("%s",str);
						}
					}
			}
		}		
	}
}

void symbolic::write_group_sc(int isc,matrix a) {

	matrix ind;
	matrix b;
	int i,i1,i2,j,var,k,k1=0;
	ind=zeros(1,a.ncols());
	char str[1024];
	
	b=a;
	
	for(i=maxder;i>=0;i--) for(i1=0;i1<=i;i1++) for(i2=0;i2<=i1;i2++) {
		var=jsc[isc][i-i1][i1-i2][i2];
		k=0;
		for(j=0;j<a.ncols();j++) {
			if(a(var,j)&&!ind(j)) {
				b.setcol(k,a.col(j));
				b(var,k)--;
				k++;
				ind(j)=1;
			}
		}
		if(k) {
			if((k!=1||b(jnum,0)>0)&&k1) printf("+");
			k1=1;
			if(k>1) {
				printf("\\left(");
				if(exist(b.block(0,N-1,0,k-1))) write(b.block(0,N-1,0,k-1));
				printf("\\right)");
			} else 
				if(exist(b.block(1,N-1,0,0))||(b(jnum,0)!=1&&b(jnum,0)!=-1)) write(b.block(0,N-1,0,k-1));
				else if(b(jnum,0)==-1) printf("-");
			if(!i) {
				sprintf(str," %s ",scname[isc]);
				write2(str,1);
			} else {
				strcpy(str,"\\frac{\\partial");
				if(i>1) sprintf(str,"%s^%d",str,i);
				sprintf(str,"%s %s}{",str,scname[isc]);
				if(i-i1) strcat(str,"\\partial\\zeta");
				if(i-i1>1) sprintf(str,"%s^%d",str,i-i1);
				if(i1-i2) strcat(str,"\\partial\\theta");
				if(i1-i2>1) sprintf(str,"%s^%d",str,i1-i2);
				if(i2) strcat(str,"\\partial\\varphi");
				if(i2>1) sprintf(str,"%s^%d",str,i2);
				strcat(str,"}");
				write3(str,1);
			}
		}
		
	}
	k=0;
	for(j=0;j<a.ncols();j++) {
		if(ind(j)==0) {
			b.setcol(k,a.col(j));
			k++;
			ind(j)=1;
		}
	}
	if(k) {
		if(b(jnum,0)>=0) printf("+");
	 	write(b.block(0,N-1,0,k-1));
	 }
}

void symbolic::write_group_sc_u(int isc,matrix a,int start) {

	matrix ind;
	matrix b;
	int i,i1,i2,j,var,k,k1=0;
	ind=zeros(1,a.ncols());
	char str[1024];
	
	b=a;
	if(start==0) k1=1;
	for(i=maxder;i>=0;i--) for(i1=0;i1<=i;i1++) for(i2=0;i2<=i1;i2++) {
		var=jsc[isc][i-i1][i1-i2][i2];
		k=0;
		for(j=0;j<a.ncols();j++) {
			if(a(var,j)&&!ind(j)) {
				b.setcol(k,a.col(j));
				b(var,k)--;
				k++;
				ind(j)=1;
			}
		}
		if(k) {
			if((k!=1||b(jnum,0)>0)&&k1) printf("+");
			k1=1;
			if(k>1) {
				printf("\\left(");
				if(exist(b.block(0,N-1,0,k-1))) write(b.block(0,N-1,0,k-1));
				printf("\\right)");
			} else 
				if(exist(b.block(1,N-1,0,0))||(b(jnum,0)!=1&&b(jnum,0)!=-1)) write(b.block(0,N-1,0,k-1));
				else if(b(jnum,0)==-1) printf("-");
			if(!i) {
				sprintf(str," %s ",scname[isc]);
				write2(str,1);
			} else {
				strcpy(str,"\\frac{\\partial");
				if(i>1) sprintf(str,"%s^%d",str,i);
				sprintf(str,"%s %s}{",str,scname[isc]);
				if(i-i1) strcat(str,"\\partial\\zeta");
				if(i-i1>1) sprintf(str,"%s^%d",str,i-i1);
				if(i1-i2) strcat(str,"\\partial\\theta");
				if(i1-i2>1) sprintf(str,"%s^%d",str,i1-i2);
				if(i2) strcat(str,"\\partial\\varphi");
				if(i2>1) sprintf(str,"%s^%d",str,i2);
				strcat(str,"}");
				write3(str,1);
			}
		}
		
	}
	
}

void symbolic::write_group_vec(int ivec,matrix a) {

	matrix ind;
	matrix b;
	int i,i1,i2,m,j,var,k,k1=0;
	ind=zeros(1,a.ncols());
	char str[1024];
	
	b=a;
	
	for(m=0;m<3;m++) for(i=maxder;i>=0;i--) for(i1=0;i1<=i;i1++) for(i2=0;i2<=i1;i2++) {
		var=jvec[ivec][m][i-i1][i1-i2][i2];
		k=0;
		for(j=0;j<a.ncols();j++) {
			if(a(var,j)&&!ind(j)) {
				b.setcol(k,a.col(j));
				b(var,k)--;
				k++;
				ind(j)=1;
			}
		}
		if(k) {
			if((k!=1||b(jnum,0)>0)&&k1) printf("+");
			k1=1;
			if(k>1) {
				printf("\\left(");
				if(exist(b.block(0,N-1,0,k-1))) write(b.block(0,N-1,0,k-1));
				printf("\\right)");
			} else 
				if(exist(b.block(1,N-1,0,0))||(b(jnum,0)!=1&&b(jnum,0)!=-1)) write(b.block(0,N-1,0,k-1));
				else if(b(jnum,0)==-1) printf("-");
			if(!i) {
				sprintf(str," %s",vecname[ivec]);
				if(m==0) strcat(str,"^\\zeta");
				if(m==1) strcat(str,"^\\theta");
				if(m==2) strcat(str,"^\\varphi");
				write2(str,1);
			} else {
				strcpy(str,"\\frac{\\partial");
				if(i>1) sprintf(str,"%s^%d",str,i);
				sprintf(str,"%s %s",str,vecname[ivec]);
				if(m==0) strcat(str,"^\\zeta");
				if(m==1) strcat(str,"^\\theta");
				if(m==2) strcat(str,"^\\varphi");
				strcat(str,"}{");
				if(i-i1) strcat(str,"\\partial\\zeta");
				if(i-i1>1) sprintf(str,"%s^%d",str,i-i1);
				if(i1-i2) strcat(str,"\\partial\\theta");
				if(i1-i2>1) sprintf(str,"%s^%d",str,i1-i2);
				if(i2) strcat(str,"\\partial\\varphi");
				if(i2>1) sprintf(str,"%s^%d",str,i2);
				strcat(str,"}");
				write3(str,1);
			}
		}
		
	}
	k=0;
	for(j=0;j<a.ncols();j++) {
		if(ind(j)==0) {
			b.setcol(k,a.col(j));
			k++;
			ind(j)=1;
		}
	}
	if(k) {
		if(b(jnum,0)>=0) printf("+");
	 	write(b.block(0,N-1,0,k-1));
	 }
}

void symbolic::write_group_vec_(int ivec,matrix a) {

	matrix ind;
	matrix b;
	int i,i1,i2,m,j,var,k,k1=0;
	ind=zeros(1,a.ncols());
	char str[1024];
	
	b=a;
	
	for(m=0;m<3;m++) for(i=maxder;i>=0;i--) for(i1=0;i1<=i;i1++) for(i2=0;i2<=i1;i2++) {
		var=jvec_[ivec][m][i-i1][i1-i2][i2];
		k=0;
		for(j=0;j<a.ncols();j++) {
			if(a(var,j)&&!ind(j)) {
				b.setcol(k,a.col(j));
				b(var,k)--;
				k++;
				ind(j)=1;
			}
		}
		if(k) {
			if((k!=1||b(jnum,0)>0)&&k1) printf("+");
			k1=1;
			if(k>1) {
				printf("\\left(");
				if(exist(b.block(0,N-1,0,k-1))) write(b.block(0,N-1,0,k-1));
				printf("\\right)");
			} else 
				if(exist(b.block(1,N-1,0,0))||(b(jnum,0)!=1&&b(jnum,0)!=-1)) write(b.block(0,N-1,0,k-1));
				else if(b(jnum,0)==-1) printf("-");
			if(!i) {
				sprintf(str," %s",vecname[ivec]);
				if(m==0) strcat(str,"^\\zeta");
				if(m==1) strcat(str,"^\\theta");
				if(m==2) strcat(str,"^\\varphi");
				write2(str,1);
			} else {
				strcpy(str,"\\frac{\\partial");
				if(i>1) sprintf(str,"%s^%d",str,i);
				sprintf(str,"%s %s",str,vecname[ivec]);
				if(m==0) strcat(str,"^\\zeta");
				if(m==1) strcat(str,"^\\theta");
				if(m==2) strcat(str,"^\\varphi");
				strcat(str,"}{");
				if(i-i1) strcat(str,"\\partial\\zeta");
				if(i-i1>1) sprintf(str,"%s^%d",str,i-i1);
				if(i1-i2) strcat(str,"\\partial\\theta");
				if(i1-i2>1) sprintf(str,"%s^%d",str,i1-i2);
				if(i2) strcat(str,"\\partial\\varphi");
				if(i2>1) sprintf(str,"%s^%d",str,i2);
				strcat(str,"}");
				write3(str,1);
			}
		}
		
	}
	k=0;
	for(j=0;j<a.ncols();j++) {
		if(ind(j)==0) {
			b.setcol(k,a.col(j));
			k++;
			ind(j)=1;
		}
	}
	if(k) {
		if(b(jnum,0)>=0) printf("+");
	 	write(b.block(0,N-1,0,k-1));
	 }
}


int symbolic::isnull(matrix a) {

	if(a.ncols()==0) return 1;
	if(!exist(abs(a.row(jnum))>1e-12)) return 1;

	return 0;
}

matrix symbolic::derive(matrix a,int var) {

	matrix b,c;
	int i,j,j1,j2,m,n,k;
	
	b=a.block(1,N-1,0,a.ncols()-1);
	b=(b!=0);
	k=sum(b);
	
	init(b,k);
	k=0;
	for(i=0;i<a.ncols();i++) {
		for(j=0;j<maxder+2;j++) for(j1=0;j1<=j;j1++) {
			if(a(jr[j-j1][j1],i)) {
				c=a.col(i);
				c(jnum)*=c(jr[j-j1][j1]);
				c(jr[j-j1][j1])--;
				if(j==maxder+1) c(jnum=0);
				else {
					if(var==0) c(jr[j-j1+1][j1])++;
					if(var==1) c(jr[j-j1][j1+1])++;
					if(var==2) c(jnum)=0;
				}
				b.setcol(k++,c);
			}	
		}
		if(a(jsint,i)) {
			c=a.col(i);
			c(jnum)*=c(jsint);
			c(jsint)--;
			if(var==1) c(jcost)++;
			else c(jnum)=0;
			b.setcol(k++,c);
		}
		if(a(jcost,i)) {
			c=a.col(i);
			c(jnum)*=c(jcost);
			c(jcost)--;
			if(var==1) {c(jsint)++;c(jnum)*=-1;}
			else c(jnum)=0;
			b.setcol(k++,c);
		}
		for(n=0;n<nsc;n++) for(j=0;j<maxder+1;j++) for(j1=0;j1<=j;j1++) for(j2=0;j2<=j1;j2++){
			if(a(jsc[n][j-j1][j1-j2][j2],i)) {
				c=a.col(i);
				c(jnum)*=c(jsc[n][j-j1][j1-j2][j2]);
				c(jsc[n][j-j1][j1-j2][j2])--;
				if(j==maxder) c(jnum=0);
				else {
					if(var==0) c(jsc[n][j-j1+1][j1-j2][j2])++;
					if(var==1) c(jsc[n][j-j1][j1-j2+1][j2])++;
					if(var==2) c(jsc[n][j-j1][j1-j2][j2+1])++;
				}
				b.setcol(k++,c);
			}	
		}
		for(n=0;n<nvec;n++) for(m=0;m<3;m++) for(j=0;j<maxder+1;j++) for(j1=0;j1<=j;j1++) for(j2=0;j2<=j1;j2++){
			if(a(jvec[n][m][j-j1][j1-j2][j2],i)) {
				c=a.col(i);
				c(jnum)*=c(jvec[n][m][j-j1][j1-j2][j2]);
				c(jvec[n][m][j-j1][j1-j2][j2])--;
				if(j==maxder) c(jnum=0);
				else {
					if(var==0) c(jvec[n][m][j-j1+1][j1-j2][j2])++;
					if(var==1) c(jvec[n][m][j-j1][j1-j2+1][j2])++;
					if(var==2) c(jvec[n][m][j-j1][j1-j2][j2+1])++;
				}
				b.setcol(k++,c);
			}	
		}
		for(n=0;n<nvec;n++) for(m=0;m<3;m++) for(j=0;j<maxder+1;j++) for(j1=0;j1<=j;j1++) for(j2=0;j2<=j1;j2++){
			if(a(jvec_[n][m][j-j1][j1-j2][j2],i)) {
				c=a.col(i);
				c(jnum)*=c(jvec_[n][m][j-j1][j1-j2][j2]);
				c(jvec_[n][m][j-j1][j1-j2][j2])--;
				if(j==maxder) c(jnum=0);
				else {
					if(var==0) c(jvec_[n][m][j-j1+1][j1-j2][j2])++;
					if(var==1) c(jvec_[n][m][j-j1][j1-j2+1][j2])++;
					if(var==2) c(jvec_[n][m][j-j1][j1-j2][j2+1])++;
				}
				b.setcol(k++,c);
			}	
		}
	}
	
	return b;

}

matrix symbolic::subst_sc(int isc,matrix a,matrix x) {

	matrix b,res,c;
	int i,j,k,l;
	int i1,j1,k1;

	for(l=0;l<2;l++) {
		k=0;
		for(i=0;i<a.ncols();i++) {
			b=a.col(i);
			for(i1=0;i1<maxder+1;i1++) for(j1=0;j1<=i1;j1++) for(k1=0;k1<=j1;k1++) {
				if(a(jsc[isc][i1-j1][j1-k1][k1],i)>0) {
					for(j=0;j<b.ncols();j++) b(jsc[isc][i1-j1][j1-k1][k1],j)-=a(jsc[isc][i1-j1][j1-k1][k1],i);
					c=x;
					for(j=0;j<i1-j1;j++) c=derive(c,0);
					for(j=0;j<j1-k1;j++) c=derive(c,1);
					for(j=0;j<k1;j++) c=derive(c,2);
					for(j=0;j<a(jsc[isc][i1-j1][j1-k1][k1],i);j++) b=prod(b,c);
				}
			}
			if(l) res.setblock(0,N-1,k,k+b.ncols()-1,b);
			k+=b.ncols();
		}
		if(!l) init(res,k);
	}
	return simplify(res);
}

matrix symbolic::subst_vec(int ivec,matrix a,matrix xz,matrix xt,matrix xp) {

	matrix b,res,c,x[3];
	int i,j,k,l;
	int i1,j1,k1,m;

	x[0]=xz;x[1]=xt;x[2]=xp;
	for(l=0;l<2;l++) {
		k=0;
		for(i=0;i<a.ncols();i++) {
			b=a.col(i);
			for(m=0;m<3;m++) {
				for(i1=0;i1<maxder+1;i1++) for(j1=0;j1<=i1;j1++) for(k1=0;k1<=j1;k1++) {
					if(a(jvec[ivec][m][i1-j1][j1-k1][k1],i)>0) {
						for(j=0;j<b.ncols();j++) b(jvec[ivec][m][i1-j1][j1-k1][k1],j)-=a(jvec[ivec][m][i1-j1][j1-k1][k1],i);
						c=x[m];
						for(j=0;j<i1-j1;j++) c=derive(c,0);
						for(j=0;j<j1-k1;j++) c=derive(c,1);
						for(j=0;j<k1;j++) c=derive(c,2);
						for(j=0;j<a(jvec[ivec][m][i1-j1][j1-k1][k1],i);j++) b=prod(b,c);
					}
				}
			}
			if(l) res.setblock(0,N-1,k,k+b.ncols()-1,b);
			k+=b.ncols();
		}
		if(!l) init(res,k);
	}
	return simplify(res);
}

matrix symbolic::subst_vec_(int ivec,matrix a,matrix xz,matrix xt,matrix xp) {

	matrix b,res,c,x[3];
	int i,j,k,l;
	int i1,j1,k1,m;
	
	x[0]=xz;x[1]=xt;x[2]=xp;
	for(l=0;l<2;l++) {
		k=0;
		for(i=0;i<a.ncols();i++) {
			b=a.col(i);
			for(m=0;m<3;m++) {
				for(i1=0;i1<maxder+1;i1++) for(j1=0;j1<=i1;j1++) for(k1=0;k1<=j1;k1++) {
					if(a(jvec_[ivec][m][i1-j1][j1-k1][k1],i)>0) {
						for(j=0;j<b.ncols();j++) b(jvec_[ivec][m][i1-j1][j1-k1][k1],j)-=a(jvec_[ivec][m][i1-j1][j1-k1][k1],i);
						c=x[m];
						for(j=0;j<i1-j1;j++) c=derive(c,0);
						for(j=0;j<j1-k1;j++) c=derive(c,1);
						for(j=0;j<k1;j++) c=derive(c,2);
						for(j=0;j<a(jvec_[ivec][m][i1-j1][j1-k1][k1],i);j++) b=prod(b,c);
					}
				}
			}
			if(l) res.setblock(0,N-1,k,k+b.ncols()-1,b);
			k+=b.ncols();
		}
		if(!l) init(res,k);
	}
	return simplify(res);
}

void symbolic::covariant(matrix &Az,matrix &At,matrix &Ap) {

	matrix Bz,Bt,Bp;
	int i;
			
	Bz=prod(g_(0,0),Az);
	Bz=add(Bz,prod(g_(0,1),At));
	Bz=add(Bz,prod(g_(0,2),Ap));
	
	Bt=prod(g_(1,0),Az);
	Bt=add(Bt,prod(g_(1,1),At));
	Bt=add(Bt,prod(g_(1,2),Ap));
	
	Bp=prod(g_(2,0),Az);
	Bp=add(Bp,prod(g_(2,1),At));
	Bp=add(Bp,prod(g_(2,2),Ap));
	
	Az=simplify(Bz);
	At=simplify(Bt);
	Ap=simplify(Bp);
}

void symbolic::contravariant(matrix &Az,matrix &At,matrix &Ap) {

	matrix Bz,Bt,Bp;
	int i;
			
	Bz=prod(g(0,0),Az);
	Bz=add(Bz,prod(g(0,1),At));
	Bz=add(Bz,prod(g(0,2),Ap));
	
	Bt=prod(g(1,0),Az);
	Bt=add(Bt,prod(g(1,1),At));
	Bt=add(Bt,prod(g(1,2),Ap));
	
	Bp=prod(g(2,0),Az);
	Bp=add(Bp,prod(g(2,1),At));
	Bp=add(Bp,prod(g(2,2),Ap));
	
	Az=simplify(Bz);
	At=simplify(Bt);
	Ap=simplify(Bp);
}

matrix symbolic::subst_covariant(int ivec,matrix a) {

	matrix Az,At,Ap;
	matrix res;
	
	Az=vec_(ivec,0);At=vec_(ivec,0);Ap=vec_(ivec,0);
	
	contravariant(Az,At,Ap);
	
	res=subst_vec(ivec,a,Az,At,Ap);
	return res;
}

matrix symbolic::subst_contravariant(int ivec,matrix a) {

	matrix Az,At,Ap;
	matrix res;
	
	Az=vec(ivec,0);At=vec(ivec,1);Ap=vec(ivec,2);
	
	covariant(Az,At,Ap);
	
	res=subst_vec_(ivec,a,Az,At,Ap);
	return res;
}

matrix symbolic::christoffel(int i,int j,int k) {

	matrix G,half;
	int l;

	init(half,1);
	init(G,1);
	half(jnum)=0.5;
	for(l=0;l<3;l++) {
		G=add(G,prod( g(i,l),derive(g_(l,j),k)  ));
		G=add(G,prod( g(i,l),derive(g_(l,k),j)  ));
		G=add(G,prod( g(i,l),neg(derive(g_(j,k),l))  ));
	}

	G=prod(half,G);

	return simplify(G);

}

matrix symbolic::gradient(matrix a,int i) {

	matrix b;
	b=derive(a,i);
	return b;
}

matrix symbolic::laplacian(matrix B) {

	int i,j,l;
	matrix A;

	init(A,1);

	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			A=add(A,derive(prod(g(i,j),derive(B,j)),i));
			for(l=0;l<3;l++) {
				A=add(A,prod(christoffel(i,l,i),prod(g(l,j),derive(B,j))));
			}
		}
	}
	return simplify(A);
}

matrix symbolic::curl(matrix v_z,matrix v_t,matrix v_p,int i) {

	matrix a,v_[3];
	int j,k;
	v_[0]=v_z;v_[1]=v_t;v_[2]=v_p;
	
	init(a,1);
	
	for(j=0;j<3;j++) {
		for(k=0;k<3;k++) {
			a=add(a,prod(perm(i,j,k),derive(v_[k],j)));
		}
	}
	
	return simplify(a);
}

matrix symbolic::vect_laplacian(matrix vz,matrix vt,matrix vp,int i) {
	
	int j,k,l,m;
	matrix A,v[3];
	
	v[0]=vz;v[1]=vt;v[2]=vp;
	
	init(A,1);
	
	for(j=0;j<3;j++) {
		for(k=0;k<3;k++) {
			A=add(A,prod(g(j,k), derive( derive(v[i],j),k) ));
			for(l=0;l<3;l++) {
				A=add(A,prod(g(j,k), derive( prod(christoffel(i,l,j),v[l]) ,k)  ));
				A=add(A,prod(g(j,k), prod(christoffel(i,l,k),derive(v[l],j))  ));
				A=add(A,prod(g(j,k), neg( prod(christoffel(l,j,k),derive(v[i],l)) ) ));
				for(m=0;m<3;m++) {
					A=add(A,prod(g(j,k), prod(christoffel(i,l,k), prod( christoffel(l,m,j) , v[m]  )  )  ));
					A=add(A,prod(g(j,k), prod(neg(christoffel(l,j,k)), prod( christoffel(i,m,l) , v[m]  )  )  ));
				}
			}
		}
	}	

	return simplify(A);
}

matrix symbolic::div(matrix vz,matrix vt,matrix vp) {

	matrix A,v[3];
	int i,l;
	
	v[0]=vz;v[1]=vt;v[2]=vp;
	init(A,1);
	for(i=0;i<3;i++) {
		A=add(A,derive(v[i],i));
		for(l=0;l<3;l++) {
			A=add(A,prod(christoffel(i,l,i),v[l]));
		}
	}
	return simplify(A);
}

matrix symbolic::div_tens(matrix T[3][3],int i) {

	matrix A;
	int j,l;
	
	init(A,1);
	for(j=0;j<3;j++) {
		A=add(A,derive(T[i][j],j));
		for(l=0;l<3;l++) {
			A=add(A,prod(christoffel(i,l,j),T[l][j]));
			A=add(A,prod(christoffel(j,l,j),T[i][l]));
		}
	}
	return simplify(A);
}

matrix symbolic::advec(matrix wz,matrix wt,matrix wp,matrix vz,matrix vt,matrix vp,int i) {

	matrix A,w[3],v[3];
	int j,l;
	
	w[0]=wz;w[1]=wt;w[2]=wp;
	v[0]=vz;v[1]=vt;v[2]=vp;
	init(A,1);
	
	for(j=0;j<3;j++) {
		A=add(A,prod(w[j],derive(v[i],j)));
		for(l=0;l<3;l++) {
			A=add(A,prod(w[j],prod(christoffel(i,l,j),v[l])));
		}
	}
	return simplify(A);

}

matrix symbolic::dot_prod(matrix wz,matrix wt,matrix wp,matrix vz,matrix vt,matrix vp) {

	matrix A,w[3],v[3];
	int j;
	
	w[0]=wz;w[1]=wt;w[2]=wp;
	v[0]=vz;v[1]=vt;v[2]=vp;
	init(A,1);
	for(j=0;j<3;j++) A=add(A,prod(w[j],v[j]));
	return simplify(A);
}

matrix symbolic::cross_prod(matrix w_z,matrix w_t,matrix w_p,matrix v_z,matrix v_t,matrix v_p,int i) {

	matrix a,v_[3],w_[3];
	int j,k;
	v_[0]=v_z;v_[1]=v_t;v_[2]=v_p;
	w_[0]=w_z;w_[1]=w_t;w_[2]=w_p;
	
	init(a,1);
	
	for(j=0;j<3;j++) {
		for(k=0;k<3;k++) {
			a=add(a,prod(perm(i,j,k),prod(w_[k],v_[j])));
		}
	}
	
	return simplify(a);
}

matrix symbolic::axisymmetric(matrix a) {

	int n,m,i,j,k,l;
	
	for(i=0;i<a.ncols();i++) {
		for(n=0;n<nsc;n++) for(j=1;j<maxder+1;j++) for(k=1;k<=j;k++) for(l=1;l<=k;l++) 
			if(a(jsc[n][j-k][k-l][l],i)) a(jnum,i)=0;
		for(n=0;n<nvec;n++) for(m=0;m<3;m++) for(j=1;j<maxder+1;j++) for(k=1;k<=j;k++) for(l=1;l<=k;l++) 
			if(a(jvec[n][m][j-k][k-l][l],i)) a(jnum,i)=0;
		for(n=0;n<nvec;n++) for(m=0;m<3;m++) for(j=1;j<maxder+1;j++) for(k=1;k<=j;k++) for(l=1;l<=k;l++) 
			if(a(jvec_[n][m][j-k][k-l][l],i)) a(jnum,i)=0;
	}
	return simplify(a);

}

matrix symbolic::spherical(matrix a) {

	int i,j,k,n;
	
	for(i=0;i<a.ncols();i++) {
		for(n=0;n<nsc;n++) for(j=1;j<maxder+1;j++) for(k=0;k<=j;k++)
			if(j==1&&k==0) a(jr[j-k][k],i)=0;
			else if(a(jr[j-k][k],i)) a(jnum,i)=0;
	}
	return simplify(a);
}

matrix symbolic::factor_sc(matrix a,int isc,int derz,int dert,int derp) {

	matrix b;
	int i,j;
	
	init(b,1);
	if(derz+dert+derp>maxder) return a;
	j=jsc[isc][derz][dert][derp];
	for(i=0;i<a.ncols();i++) {
		if(a(j,i)) {
			a(j,i)--;
			b=add(b,a.col(i));
		}
	}
	
	return simplify(b);

}

matrix symbolic::factor_vec(matrix a,int ivec,int icomp,int derz,int dert,int derp) {

	matrix b;
	int i,j;
	
	init(b,1);
	if(derz+dert+derp>maxder) return a;
	j=jvec[ivec][icomp][derz][dert][derp];
	for(i=0;i<a.ncols();i++) {
		if(a(j,i)) {
			a(j,i)--;
			b=add(b,a.col(i));
		}
	}
	
	return simplify(b);

}

matrix symbolic::factor_vec_(matrix a,int ivec,int icomp,int derz,int dert,int derp) {

	matrix b;
	int i,j;
	
	init(b,1);
	if(derz+dert+derp>maxder) return a;
	j=jvec_[ivec][icomp][derz][dert][derp];
	for(i=0;i<a.ncols();i++) {
		if(a(j,i)) {
			a(j,i)--;
			b=add(b,a.col(i));
		}
	}
	
	return simplify(b);

}

matrix symbolic::stress(matrix vz,matrix vt,matrix vp,int i,int j) {

	matrix A,v[3];
	int k,l;
	
	v[0]=vz;v[1]=vt;v[2]=vp;
	init(A,1);
	A(jnum,0)=2./3;
	A=prod(prod(A,div(vz,vt,vp)),g(i,j));
	A=neg(A);
	for(k=0;k<3;k++) {
		A=add(A, prod(g(i,k), derive(v[j],k) ) );
		A=add(A, prod(g(j,k), derive(v[i],k) ) );
		for(l=0;l<3;l++) {
			A=add(A, prod(g(i,k), prod(christoffel(j,l,k),v[l]) ));
			A=add(A, prod(g(j,k), prod(christoffel(i,l,k),v[l]) ));
		}
	}
	
	return simplify(A);

}


