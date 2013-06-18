#include"symbolic.h"
#include<stdlib.h>
#include<string.h>

symterm::symterm(const symbolic *context) {

	nvar=context->nvar;
	maxder=context->maxder;
	next=NULL;
	num=0;
	sint=0;
	cost=0;
	r=new int*[maxder+2];
	for(int i=0;i<maxder+2;i++) {
		r[i]=new int[maxder+2-i];
		for(int j=0;j<maxder+2-i;j++) r[i][j]=0;
	}
	var=new int**[nvar];
	for(int n=0;n<nvar;n++) {
		var[n]=new int*[maxder+1];
		for(int i=0;i<maxder+1;i++) {
			var[n][i]=new int[maxder+1-i];
			for(int j=0;j<maxder+1-i;j++) var[n][i][j]=0;
		}
	}
}

symterm::symterm(const symterm &t) {
	
	nvar=t.nvar;
	maxder=t.maxder;
	next=NULL;
	num=t.num;
	sint=t.sint;
	cost=t.cost;
	r=new int*[maxder+2];
	for(int i=0;i<maxder+2;i++) {
		r[i]=new int[maxder+2-i];
		for(int j=0;j<maxder+2-i;j++) r[i][j]=t.r[i][j];
	}
	var=new int**[nvar];
	for(int n=0;n<nvar;n++) {
		var[n]=new int*[maxder+1];
		for(int i=0;i<maxder+1;i++) {
			var[n][i]=new int[maxder+1-i];
			for(int j=0;j<maxder+1-i;j++) var[n][i][j]=t.var[n][i][j];
		}
	}
}

symterm::~symterm() {

	for(int i=0;i<maxder+2;i++) delete [] r[i];
	delete [] r;
	for(int n=0;n<nvar;n++) {
		for(int i=0;i<maxder+1;i++) delete [] var[n][i];
		delete [] var[n];
	}
	delete [] var;

}

bool symterm::operator==(const symterm &t) const {

	for(int i=0;i<maxder+2;i++) 
		for(int j=0;j<maxder+2-i;j++) 
			if(r[i][j]!=t.r[i][j]) return false;
	if(sint!=t.sint) return false;
	if(cost!=t.cost) return false;
	for(int n=0;n<nvar;n++) 
		for(int i=0;i<maxder+1;i++) 
			for(int j=0;j<maxder+1-i;j++) 
				if(var[n][i][j]!=t.var[n][i][j]) return false;

	return true;
}

symterm symterm::operator*(const symterm &t) const {

	symterm tnew(*this);
	
	tnew.num=num*t.num;
	for(int i=0;i<maxder+2;i++) 
		for(int j=0;j<maxder+2-i;j++) 
			tnew.r[i][j]=r[i][j]+t.r[i][j];
	tnew.sint=sint+t.sint;
	tnew.cost=cost+t.cost;
	for(int n=0;n<nvar;n++) 
		for(int i=0;i<maxder+1;i++) 
			for(int j=0;j<maxder+1-i;j++) 
				tnew.var[n][i][j]=var[n][i][j]+t.var[n][i][j];
	
	return tnew;

}

symterm symterm::operator/(const symterm &t)  const {

	symterm tnew(*this);
	
	tnew.num=num/t.num;
	for(int i=0;i<maxder+2;i++) 
		for(int j=0;j<maxder+2-i;j++) 
			tnew.r[i][j]=r[i][j]-t.r[i][j];
	tnew.sint=sint-t.sint;
	tnew.cost=cost-t.cost;
	for(int n=0;n<nvar;n++) 
		for(int i=0;i<maxder+1;i++) 
			for(int j=0;j<maxder+1-i;j++) 
				tnew.var[n][i][j]=var[n][i][j]-t.var[n][i][j];
	
	return tnew;

}

void symterm::dump() const {

	printf("\tnum\t: %f\n",num);
	printf("\tr\t: ");
	for(int i=0;i<maxder+2;i++) 
		for(int j=0;j<maxder+2-i;j++) printf("[%d%d] %d ",i,j,r[i][j]);
	printf("\n");
	printf("\tsint\t: %d\n",sint);
	printf("\tcost\t: %d\n",cost);
	for(int n=0;n<nvar;n++) {
		printf("\t(%d)\t: ",n);
		for(int i=0;i<maxder+1;i++) 
			for(int j=0;j<maxder+1-i;j++) printf("[%d%d] %d ",i,j,var[n][i][j]);
		printf("\n");
	}
	
}

symbolic::symbolic(int n_var,int max_der): r(_r),sint(_sint),cost(_cost),one(_one),
				_g(CONTRAVARIANT,CONTRAVARIANT), _g_(COVARIANT,COVARIANT), g(_g),g_(_g_),sqrt_g(_sqrt_g) {

	if(!n_var) {
		fprintf(stderr,"Symbolic: nvar must be greater than zero\n");
		exit(1);
	}

	nvar=n_var;
	maxder=max_der;
	tol=1e-14;
	
	var_name=new char*[nvar];
	var_value=new matrix[nvar];
	var_par=new int[nvar];
	is_const=new bool[nvar];
	for(int i=0;i<nvar;i++) {
		var_name[i]=new char[32];
		var_name[i][0]='\0';
		var_par[i]=0;
		is_const[i]=false;
	}
	
	symterm *t;
	
	_r.context=this;
	t=_r.addterm();
	t->r[0][0]=1;
	t->num=1;
	
	_sint.context=this;
	t=_sint.addterm();
	t->sint=1;
	t->num=1;
	
	_cost.context=this;
	t=_cost.addterm();
	t->cost=1;
	t->num=1;
	
	_one.context=this;
	t=_one.addterm();
	t->num=1;
	
	_g.set_context(this);
	_g(0,0)=(r*r+Dt(r)*Dt(r))/r/r/Dz(r)/Dz(r);
	_g(0,1)=-Dt(r)/r/r/Dz(r);
	_g(1,0)=g(0,1);
	_g(1,1)=1/r/r;
	_g(2,2)=1/r/r/sint/sint;
	
	_g_.set_context(this);
	_g_(0,0)=Dz(r)*Dz(r);
	_g_(0,1)=Dz(r)*Dt(r);
	_g_(1,0)=g_(0,1);
	_g_(1,1)=r*r+Dt(r)*Dt(r);
	_g_(2,2)=r*r*sint*sint;
	
	_sqrt_g=r*r*Dz(r)*sint;
	
	
}

symbolic::~symbolic() {
	
	for(int i=0;i<nvar;i++) 
		delete [] var_name[i];
	delete [] var_name;
	delete [] var_value;
	delete [] var_par;
	delete [] is_const;

}

sym symbolic::regvar(const char *name) {

	if(name[0]=='\0') {
		fprintf(stderr,"Symbolic: Invalid name\n");
		exit(1);
	}
	if(!strcmp(name,"r")) {
		fprintf(stderr,"Symbolic: Variable name \"r\" is reserved\n");
		exit(1);
	}
	int i=0;
	while(var_name[i][0]!='\0') {
		if(!strcmp(var_name[i],name)) {
			fprintf(stderr,"Symbolic: Can't register variable (already registered)\n");
			exit(1);
		}
		i++;
		if(i==nvar) {
			fprintf(stderr,"Symbolic: Can't register variable (increase nvar)\n");
			exit(1);
		}
	}
	strncpy(var_name[i],name,31);
	
	return var(name);

}

sym symbolic::regconst(const char *name) {

	sym s;
	s=regvar(name);
	is_const[id(name)]=true;
	
	return s;

}

int symbolic::id(const char *name) const {

	if(name[0]=='\0') {
		fprintf(stderr,"Symbolic: Invalid name\n");
		exit(1);
	}
	
	if(!strcmp(name,"r")) return -1;	
	
	int i=0;
	while(strcmp(var_name[i],name)) {
		i++;
		if(i==nvar) {
			fprintf(stderr,"Symbolic: Unknown variable \"%s\"\n",name);
			exit(1);
		}
	}
	return i;

}

sym symbolic::var(const char *name) const {

	int i=id(name);
	
	sym s;
	
	if(i==-1) {
		s=r;
		return s;
	}
	
	symterm *t;
	
	s.context=this;
	t=s.addterm();
	t->var[i][0][0]=1;
	t->num=1;

	return s;

}

void symbolic::set_map(const mapping &map_set) {

	map=map_set;

}

void symbolic::set_value(const char *name,const matrix &value,int parity) {

	int i=id(name);
	if(i==-1) {
		fprintf(stderr,"Symbolic: Can't assign value to variable \"r\", use set_map instead\n");
		exit(1);
	}
	var_value[i]=value;
	if(parity!=0&&parity!=1&&parity!=10&&parity!=11) {
		fprintf(stderr,"Symbolic: Unknown parity %d\n",parity);
		exit(1);
	}
	var_par[i]=parity;

}

void symbolic::write_term(const symterm &t,FILE *fp) const {
	
	char str[64];
	
	fprintf(fp,"%+f",t.num);
	for(int i=0;i<maxder+2;i++) {
		for(int j=0;j<maxder+2-i;j++) {
			if(t.r[i][j]) {
				if(i==0&&j==0) strcpy(str,"r");
				else if(i==1&&j==0) strcpy(str,"rz");
				else if(i==2&&j==0) strcpy(str,"rzz");
				else if(i==0&&j==1) strcpy(str,"rt");
				else if(i==0&&j==2) strcpy(str,"rtt");
				else if(i==1&&j==1) strcpy(str,"rzt");
				else {
					int i2,j2;
					if(i==0) { strcpy(str,"rtt");i2=0;j2=j-2; }
					if(i==1) { strcpy(str,"rz");i2=0;j2=j; }
					if(i>=2) { strcpy(str,"rzz");i2=i-2;j2=j; }
					format_deriv(str,i2,j2,00);
				}
				write_var(fp,str,t.r[i][j]);
			}
		}
	}
	if(t.sint) write_var(fp,"sin(th)",t.sint);
	if(t.cost) write_var(fp,"cos(th)",t.cost);
	for(int n=0;n<nvar;n++) {
		for(int i=0;i<maxder+1;i++) {
			for(int j=0;j<maxder+1-i;j++) {
				if(t.var[n][i][j]) {
					strcpy(str,var_name[n]);
					format_deriv(str,i,j,var_par[n]);
					write_var(fp,str,t.var[n][i][j]);
				}
			}
		}
	}
	
}

void symbolic::format_deriv(char *str,int dz,int dt,int par) const {

	if(dz==0&&dt==0) return;
	char prefix[64],suffix[32];
	
	strcpy(prefix,"(");
	for(int i=0;i<dz;i++) strcat(prefix,"D,");
	
	char Dt2[16],Dt[16];
	if(par==00) { strcpy(Dt2,"Dt2");strcpy(Dt,"Dt"); }
	if(par==01) { strcpy(Dt2,"leg.D2_01");strcpy(Dt,"leg.D_01"); }
	if(par==10) { strcpy(Dt2,"leg.D2_10");strcpy(Dt,"leg.D_10"); }
	if(par==11) { strcpy(Dt2,"leg.D2_11");strcpy(Dt,"leg.D_11"); }
	suffix[0]='\0';
	for(int i=0;i<dt/2;i++) { strcat(suffix,",");strcat(suffix,Dt2);}
	if(dt%2) { strcat(suffix,",");strcat(suffix,Dt);}
	strcat(suffix,")");

	strcat(prefix,str);
	strcpy(str,prefix);
	strcat(str,suffix);
	
}

void symbolic::write_var(FILE *fp,const char *str,int n) const {

	if(n>3) fprintf(fp,"*pow(%s,%d)",str,n);
	else if(n<-3) fprintf(fp,"/pow(%s,%d)",str,-n);
	else if(n>0) 
		for(int i=0;i<n;i++) fprintf(fp,"*%s",str);
	else if(n<0) 
		for(int i=0;i<-n;i++) fprintf(fp,"/%s",str);

}

matrix symbolic::eval(const symterm &t) const {
	
	matrix res,m;

	res=t.num*ones(1,1);
	
	for(int i=0;i<maxder+2;i++) {
		for(int j=0;j<maxder+2-i;j++) {
			if(t.r[i][j]) {
				if(i==0&&j==0) m=map.r;
				else if(i==1&&j==0) m=map.rz;
				else if(i==2&&j==0) m=map.rzz;
				else if(i==0&&j==1) m=map.rt;
				else if(i==0&&j==2) m=map.rtt;
				else if(i==1&&j==1) m=map.rzt;
				else {
					int i2,j2;
					if(i==0) { m=map.rtt;i2=0;j2=j-2; }
					if(i==1) { m=map.rz;i2=0;j2=j; }
					if(i>=2) { m=map.rzz;i2=i-2;j2=j; }
					eval_deriv(m,i2,j2,00);
				}
				res*=eval_var(m,t.r[i][j]);
			}
		}
	}
	if(t.sint) res*=eval_var(sin(map.th)*ones(map.gl.N,1),t.sint);
	if(t.cost) res*=eval_var(cos(map.th)*ones(map.gl.N,1),t.cost);
	for(int n=0;n<nvar;n++) {
		for(int i=0;i<maxder+1;i++) {
			for(int j=0;j<maxder+1-i;j++) {
				if(t.var[n][i][j]) {
					m=var_value[n];
					eval_deriv(m,i,j,var_par[n]);
					res*=eval_var(m,t.var[n][i][j]);
				}
			}
		}
	}
	
	return res;
}

void symbolic::eval_deriv(matrix &m,int dz,int dt,int par) const {

	if(dz==0&&dt==0) return;
	matrix_block_diag pre;
	matrix post;
	
	pre=eye(map.D);
	for(int i=0;i<dz;i++) pre=(map.D,pre);
	
	matrix Dt2,Dt;
	if(par==00) { Dt2=map.leg.D2_00;Dt=map.leg.D_00; }
	if(par==01) { Dt2=map.leg.D2_01;Dt=map.leg.D_01; }
	if(par==10) { Dt2=map.leg.D2_10;Dt=map.leg.D_10; }
	if(par==11) { Dt2=map.leg.D2_11;Dt=map.leg.D_11; }
	post=eye(map.leg.npts);
	for(int i=0;i<dt/2;i++) post=(post,Dt2);
	if(dt%2) post=(post,Dt);
	
	m=(pre,m,post);
	
}

matrix symbolic::eval_var(const matrix &m,int n) const {

	matrix res;
	
	res=ones(m.nrows(),m.ncols());

	if(n>3||n<-3) res=pow(m,n);
	else if(n>0) 
		for(int i=0;i<n;i++) res*=m;
	else if(n<0) 
		for(int i=0;i<-n;i++) res/=m;
		
	return res;

}

sym symbolic::Dz(const sym &s) const {

	sym snew;
	snew.context=s.context;
	
	int maxder=s.context->maxder;
	int nvar=s.context->nvar;
	
	symterm *t,*tnew;
	
	t=s.terms;
	while(t!=NULL) {
		for(int i=0;i<maxder+2;i++) {
			for(int j=0;j<maxder+2-i;j++) {
				if(t->r[i][j]!=0) {
					if(i==maxder+1||j==maxder+1-i) {
						fprintf(stderr,"WARNING: Symbolic: Maximum order of derivative exceeded (increase maxder)\n");
						continue;
					}
					tnew=snew.addterm(*t);
					tnew->num*=tnew->r[i][j];
					tnew->r[i][j]--;
					tnew->r[i+1][j]++;
				}
			}
		}
		for(int n=0;n<nvar;n++) {
			for(int i=0;i<maxder+1;i++) {
				for(int j=0;j<maxder-i+1;j++) {
					if(t->var[n][i][j]!=0&&!is_const[n]) {
						if(i==maxder||j==maxder-i) {
							fprintf(stderr,"WARNING: Symbolic: Maximum order of derivative exceeded (increase maxder)\n");
							continue;
						}	
						tnew=snew.addterm(*t);
						tnew->num*=tnew->var[n][i][j];
						tnew->var[n][i][j]--;
						tnew->var[n][i+1][j]++;
					}
				}
			}
		}
		t=t->next;
	}
	
	return snew.simplify();

}

sym symbolic::Dt(const sym &s) const {

	sym snew;
	snew.context=s.context;
	
	int maxder=s.context->maxder;
	int nvar=s.context->nvar;
	
	symterm *t,*tnew;
	
	t=s.terms;
	while(t!=NULL) {
		for(int i=0;i<maxder+2;i++) {
			for(int j=0;j<maxder+2-i;j++) {
				if(t->r[i][j]!=0) {
					if(i==maxder+1||j==maxder+1-i) {
						fprintf(stderr,"WARNING: Symbolic: Maximum order of derivative exceeded (increase maxder)\n");
						continue;
					}
					tnew=snew.addterm(*t);
					tnew->num*=tnew->r[i][j];
					tnew->r[i][j]--;
					tnew->r[i][j+1]++;
				}
			}
		}
		if(t->sint!=0) {
			tnew=snew.addterm(*t);
			tnew->num*=tnew->sint;
			tnew->sint--;
			tnew->cost++;
		}
		if(t->cost!=0) {
			tnew=snew.addterm(*t);
			tnew->num*=-tnew->cost;
			tnew->cost--;
			tnew->sint++;
		}
		for(int n=0;n<nvar;n++) {
			for(int i=0;i<maxder+1;i++) {
				for(int j=0;j<maxder+1-i;j++) {
					if(t->var[n][i][j]!=0&&!is_const[n]) {
						if(i==maxder||j==maxder-i) {
							fprintf(stderr,"WARNING: Symbolic: Maximum order of derivative exceeded (increase maxder)\n");
							continue;
						}
						tnew=snew.addterm(*t);
						tnew->num*=tnew->var[n][i][j];
						tnew->var[n][i][j]--;
						tnew->var[n][i][j+1]++;
					}
				}
			}
		}
		t=t->next;
	}

	return snew.simplify();

}


sym DzDt(const sym &s) {
	
	return DzDt(s,1,1);

}

sym Dz(const sym &s,int n) {

	return DzDt(s,n,0);

}
sym Dt(const sym &s,int n) {

	return DzDt(s,0,n);

}

sym DzDt(const sym &s,int nz,int nt) {

	sym snew(s);

	for(int i=0;i<nz;i++) snew=Dz(snew);
	for(int i=0;i<nt;i++) snew=Dt(snew);
	
	return snew;
			
}

sym_vec symbolic::contravariant(const sym_vec &s) const {

	if(s.type==CONTRAVARIANT) return s;
	
	sym_vec snew(CONTRAVARIANT);
	snew.set_context(this);
	
	for(int i=0;i<3;i++) snew(i)=g(i,0)*s(0)+g(i,1)*s(1)+g(i,2)*s(2);

	return snew;

}

sym_vec symbolic::covariant(const sym_vec &s) const {

	if(s.type==COVARIANT) return s;
	
	sym_vec snew(COVARIANT);
	snew.set_context(this);
	
	for(int i=0;i<3;i++) snew(i)=g_(i,0)*s(0)+g_(i,1)*s(1)+g_(i,2)*s(2);

	return snew;

}

sym_tens symbolic::contravariant_contravariant(const sym_tens &s) const {

	sym_tens snew(s),s2(s);
	snew.set_type(CONTRAVARIANT,CONTRAVARIANT);

	if(s.type[0]!=CONTRAVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g(i,0)*s(0,j)+g(i,1)*s(1,j)+g(i,2)*s(2,j);	
	}
	if(s.type[1]!=CONTRAVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g(0,j)+s2(i,1)*g(1,j)+s2(i,2)*g(2,j);	
	}
	
	return snew;

}

sym_tens symbolic::contravariant_covariant(const sym_tens &s) const {

	sym_tens snew(s),s2(s);
	snew.set_type(CONTRAVARIANT,COVARIANT);
	
	if(s.type[0]!=CONTRAVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g(i,0)*s(0,j)+g(i,1)*s(1,j)+g(i,2)*s(2,j);	
	}
	if(s.type[1]!=COVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g_(0,j)+s2(i,1)*g_(1,j)+s2(i,2)*g_(2,j);	
	}
	
	return snew;

}

sym_tens symbolic::covariant_contravariant(const sym_tens &s) const {

	sym_tens snew(s),s2(s);
	snew.set_type(COVARIANT,CONTRAVARIANT);

	if(s.type[0]!=COVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g_(i,0)*s(0,j)+g_(i,1)*s(1,j)+g_(i,2)*s(2,j);	
	}
	if(s.type[1]!=CONTRAVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g(0,j)+s2(i,1)*g(1,j)+s2(i,2)*g(2,j);	
	}
	
	return snew;

}

sym_tens symbolic::covariant_covariant(const sym_tens &s) const {

	sym_tens snew(s),s2(s);
	snew.set_type(COVARIANT,COVARIANT);

	if(s.type[0]!=COVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g_(i,0)*s(0,j)+g_(i,1)*s(1,j)+g_(i,2)*s(2,j);	
	}
	if(s.type[1]!=COVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g_(0,j)+s2(i,1)*g_(1,j)+s2(i,2)*g_(2,j);	
	}
	
	return snew;

}

sym symbolic::spherical(const sym &s) const {

	sym snew(s);
	symterm *t;
	

	for(t=snew.terms;t!=NULL;t=t->next) {
		for(int i=0;i<maxder+2;i++) 
			for(int j=0;j<maxder+2-i;j++) {
				if(i||j) {
					if(i==1&&j==0) t->r[i][j]=0;
					else {
						if(t->r[i][j]) t->num=0;
					}
				}
			}
	}
	
	return snew.simplify();

}

sym_vec symbolic::spherical(const sym_vec &s) const {

	sym_vec snew(s);
	
	for(int i=0;i<3;i++)
		snew(i)=spherical(s(i));

	return snew;
}

sym_tens symbolic::spherical(const sym_tens &s) const {

	sym_tens snew(s);
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			snew(i,j)=spherical(s(i,j));

	return snew;

}


sym symbolic::d(const sym &s,int i) const {

	switch(i) {
		case 0: return Dz(s);
		case 1: return Dt(s);
	}
	
	return 0*one;

}

sym symbolic::christoffel(int i,int j,int k) const {
	
	sym G;
	G.context=this;
	
	for(int l=0;l<3;l++) 
		G=G+g(i,l)*( d(g_(l,j),k) + d(g_(l,k),j) - d(g_(j,k),l) );
			
	return G*0.5;	
}

sym symbolic::covderiv(const sym_vec &v,int i,int k) const {

	sym dv(d(v(i),k));

	if(v.type==CONTRAVARIANT) {
		for(int l=0;l<3;l++) dv=dv+G(i,l,k)*v(l);
	} else {
		for(int l=0;l<3;l++) dv=dv-G(l,i,k)*v(l);
	}
	
	return dv;

}

sym_vec symbolic::gradient(const sym &s) const {

	sym_vec grad(COVARIANT);
	grad.set_context(this);
	
	grad(0)=Dz(s);
	grad(1)=Dt(s);
	grad(2)=0*one;

	return grad;
}

sym_tens symbolic::gradient(const sym_vec &v) const {

	sym_tens grad(v.type,COVARIANT);
	grad.set_context(this);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			grad(i,j)=covderiv(v,i,j);

	return grad;
}

sym symbolic::divergence(const sym_vec &v) const {

	sym s;
	sym_vec V;
	s.context=this;
	
	V=contravariant(v);
	
	for(int i=0;i<3;i++) 
		s=s+d(V(i),i)+d(sqrt_g,i)/sqrt_g*V(i);
		//s=s+covderiv(V,i,i);
		
		
	return s;

}

sym_vec symbolic::divergence(const sym_tens &t) const {

	sym_vec v(CONTRAVARIANT);
	sym_tens T;
	v.set_context(this);
	
	T=contravariant_contravariant(t);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) {
			v(i)=v(i)+d(T(i,j),j);
			for(int k=0;k<3;k++)
				v(i)=v(i)+G(i,k,j)*T(k,j)+G(j,k,j)*T(i,k);
		}
	return v.set_variance(t.type[0]);

}

sym symbolic::laplacian(const sym &s) const {

	return divergence(gradient(s));

}

double symbolic::perm(int i,int j,int k) const {

	if(j==(i+1)%3&&k==(j+1)%3) return 1;
	if(j==(i+2)%3&&k==(j+2)%3) return -1;
	return 0;

}

sym_vec symbolic::curl(const sym_vec &v) const {

	sym_vec res(CONTRAVARIANT),V;
	res.set_context(this);
	
	V=covariant(v);
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				res(i)=res(i)+perm(i,j,k)*covderiv(V,k,j)/sqrt_g;
	
	return res.set_variance(v);
	

}


sym_vec symbolic::laplacian(const sym_vec &v) const {

	return (div(grad(v))).set_variance(v);

}

sym_tens symbolic::stress(const sym_vec &v) const {

	sym_tens t;
	
	t.set_context(this);
	
	t=grad(v)+grad(v).T()-2./3.*g*div(v);
				
	return t;
}






