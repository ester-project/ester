#include"symbolic.h"
#include"solver.h"
#include<stdlib.h>
#include<cmath>
#include<string.h>

sym::sym() {
	terms=NULL;
	context=NULL;
}

sym::~sym() {
	
	destroy();
}

void sym::destroy() {

	symterm *t1,*t2;
	
	t1=terms;
	while(t1!=NULL) {
		t2=t1;
		t1=t1->next;
		delete t2;
	}
	terms=NULL;
}

sym::sym(const sym &s) {

	copy(s);
}

void sym::copy(const sym &s) {

	context=s.context;
	
	if(s.terms==NULL) {
		terms=NULL;
		return;
	}
	
	symterm *t1,*t2;
	terms=new symterm(*s.terms);
	t1=terms;
	t2=s.terms->next;
	while(t2!=NULL) {
		t1->next=new symterm(*t2);
		t1=t1->next;
		t2=t2->next;
	}

}

sym &sym::operator=(const sym &s) {

	if(&s==this) return *this;
	destroy();
	copy(s);
	return *this;

}

const symbolic *sym::check_context() const {

	if(context==NULL) {
		fprintf(stderr,"Symbolic: sym object not initialized\n");
		exit(1);
	}
	return context;
}

const symbolic *sym::check_context(const sym &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym::check_context(const sym_vec &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym::check_context(const sym_tens &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym_vec::check_context() const {

	const symbolic *context;

	context=s[0].check_context();
	for(int i=1;i<3;i++) 
		if(context!=s[i].check_context()) {
			fprintf(stderr,"Symbolic: sym_vec components have different contexts\n");
			exit(1);
		}

	return context;
}

const symbolic *sym_vec::check_context(const sym &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym_vec::check_context(const sym_vec &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym_vec::check_context(const sym_tens &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym_tens::check_context() const {

	const symbolic *context;

	context=s[0][0].check_context();
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) 
			if(context!=s[i][j].check_context()) {
				fprintf(stderr,"Symbolic: sym_tens components have different contexts\n");
				exit(1);
			}

	return context;
}

const symbolic *sym_tens::check_context(const sym &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym_tens::check_context(const sym_vec &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

const symbolic *sym_tens::check_context(const sym_tens &s) const {
	
	const symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		fprintf(stderr,"Symbolic: Wrong context\n");
		exit(1);
	}
	return c;
}

void sym_tens::set_context(const symbolic *context) {

	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) 
			s[i][j].context=context;

}

void sym_vec::set_context(const symbolic *context) {

	for(int i=0;i<3;i++) 
		s[i].context=context;

}
symterm *sym::addterm() {

	symterm **t;
	
	t=&terms;
	while(*t!=NULL) 
		t=&((*t)->next);
	*t=new symterm(context);
	
	return *t;

}

symterm *sym::addterm(const symterm &st) {

	symterm **t;
	
	t=&terms;
	while(*t!=NULL) 
		t=&((*t)->next);
	*t=new symterm(st);
	
	return *t;

}


void sym::dump() const {

	symterm *t1;
	int n=0;
	
	printf("Context(%p)\n",context);
	t1=terms;
	while(t1!=NULL) {
		printf("(%d): \n",n);
		t1->dump();
		t1=t1->next;
		n++;
	}

}

sym &sym::simplify() {

	if(terms==NULL) return *this;
	bool modified=true;
	while(modified) {
		simplify0();
		modified=trig_simplify();
	}
	trig_simplify2();
		
	return *this;
}

void sym::simplify0() {

	symterm **t1,**t2,*t;
	t1=&terms;
	while(*t1!=NULL) {
		t2=&((*t1)->next);
		while(*t2!=NULL) {
			if(**t1==**t2) {
				(*t1)->num+=(*t2)->num;
				t=(*t2)->next;
				delete *t2;
				*t2=t;
			} else t2=&((*t2)->next);
		}
		if(fabs((*t1)->num)<=context->tol) {
			t=(*t1)->next;
			delete *t1;
			*t1=t;
		} else t1=&((*t1)->next);
	}

}

bool sym::trig_simplify() {

	symterm *t1,**t2,*t;
	bool modified=false;
	t1=terms;
	while(t1!=NULL) {
		t2=&terms;
		while(*t2!=NULL) {
			if(!(t1==*t2)) {
				t1->sint-=2;
				(*t2)->cost-=2;
				if(*t1==**t2&&fabs(t1->num-(*t2)->num)<=context->tol) {
					modified=true;
					t=(*t2)->next;
					delete *t2;
					*t2=t;
				} else {
					t1->sint+=2;
					(*t2)->cost+=2;
					t2=&((*t2)->next);
				}
			} else t2=&((*t2)->next);
		}
		t1=t1->next;
	}
	return modified;

}

void sym::trig_simplify2() {

	symterm *t1,**t2,*t;
	
	t1=terms;
	while(t1!=NULL) {
		t2=&terms;
		while(*t2!=NULL) {
			if(!(t1==*t2)) {
				(*t2)->sint-=2;
				if(*t1==**t2&&fabs(t1->num+(*t2)->num)<=context->tol) {
					t1->cost+=2;
					t=(*t2)->next;
					delete *t2;
					*t2=t;
				} else {
					(*t2)->sint+=2;
					(*t2)->cost-=2;
					if(*t1==**t2&&fabs(t1->num+(*t2)->num)<=context->tol) {
						t1->sint+=2;
						t=(*t2)->next;
						delete *t2;
						*t2=t;
					} else {
						(*t2)->cost+=2;
						t2=&((*t2)->next);
					}
				}
			} else t2=&((*t2)->next);
		}
		t1=t1->next;
	}

}

sym sym::operator+(const sym &s) const {

	check_context(s);
	sym snew(*this);
	symterm *t1,*t2;
	t1=snew.terms;
	snew.copy(s);
	t2=snew.terms;
	if(t1!=NULL) {
		snew.terms=t1;
		while(t1->next!=NULL) t1=t1->next;
		t1->next=t2;
	} else snew.terms=t2;
	
	return snew.simplify();

}

sym sym::operator+(const double x) const {

	check_context();
	sym snum(context->one);
	
	snum.terms->num=x;
	
	return *this+snum;

}

sym operator+(const double x,const sym &s) {

	return s+x;
}

sym operator+(const sym &s) {

	s.check_context();
	return s;
}

sym sym::operator-(const sym &s) const {

	return *this+(-s);

}

sym sym::operator-(const double x) const {

	return *this+(-x);

}

sym operator-(const double x,const sym &s) {

	return (-s)+x;
}

sym operator-(const sym &s) {

	s.check_context();
	sym snew(s);
	symterm *t;
	t=snew.terms;
	while(t!=NULL) {
		t->num=-t->num;
		t=t->next;
	}
	
	return snew;
}

sym sym::operator*(const sym &s) const {

	check_context(s);
	sym snew;
	snew.context=context;
	
	symterm *t1,*t2;
	t1=terms;
	while(t1!=NULL) {
		t2=s.terms;
		while(t2!=NULL) {
			snew.addterm((*t1)*(*t2));
			t2=t2->next;
		}
		t1=t1->next;
	}

	return snew.simplify();
}

sym sym::operator*(const double x) const {

	check_context();
	sym snew(*this);
	symterm *t;
	t=snew.terms;
	while(t!=NULL) {
		t->num*=x;
		t=t->next;
	}
	
	return snew;
}

sym operator*(const double x,const sym &s) {

	return s*x;

}

sym sym::operator/(const sym &s) const  {

	check_context(s);
	if(s.terms==NULL) {
		fprintf(stderr,"Symbolic: Division by zero\n");
		exit(1);
	}
	if(s.terms->next!=NULL) {
		fprintf(stderr,"Symbolic: Sorry, but division by a sym object with more than one term is not allowed by the implementation\n");
		exit(1);
	}
	sym snew;snew.context=context;
	symterm *t;
	t=terms;
	while(t!=NULL) {
		snew.addterm((*t)/(*s.terms));
		t=t->next;
	}
	return snew.simplify();

}

sym sym::operator/(const double x) const {
	
	return (*this)*(1./x);

}

sym operator/(const double x,const sym &s) {
	
	return x*s.context->one/s;

}



void sym::write(FILE *fp) const {

	check_context();

	symterm *t;

	t=terms;
	if(t==NULL) fprintf(fp,"0");
	while(t!=NULL) {
		context->write_term(*t,fp);
		t=t->next;
	}
	fprintf(fp,"\n");

}

matrix sym::eval() const {

	check_context();

	symterm *t;
	matrix res;

	t=terms;
	res=zeros(1,1);
	while(t!=NULL) {
		res+=context->eval(*t);
		t=t->next;
	}
	return res*ones(context->map.gl.N,context->map.leg.npts);

}

sym_vec::sym_vec(sym_vec_type type_set) {
	
	type=type_set;

}

sym_tens::sym_tens(sym_vec_type type_set0,sym_vec_type type_set1) {
	
	type[0]=type_set0;
	type[1]=type_set1;

}

sym &sym_vec::operator()(int i) {
	
	if(i<0||i>2) {
		fprintf(stderr,"Symbolic: sym_vec index out of range\n");
		exit(1);
	}
	return s[i];
	
}

sym &sym_tens::operator()(int i,int j) {
	
	if(i<0||i>2||j<0||j>2) {
		fprintf(stderr,"Symbolic: sym_tens index out of range\n");
		exit(1);
	}
	return s[i][j];
	
}

sym sym_vec::operator()(int i) const {
	
	if(i<0||i>2) {
		fprintf(stderr,"Symbolic: sym_vec index out of range\n");
		exit(1);
	}
	return s[i];
	
}

sym sym_tens::operator()(int i,int j) const {
	
	if(i<0||i>2||j<0||j>2) {
		fprintf(stderr,"Symbolic: sym_tens index out of range\n");
		exit(1);
	}
	return s[i][j];
	
}

bool sym_vec::is_covariant() {return type==COVARIANT;}
bool sym_vec::is_contravariant() {return type==CONTRAVARIANT;}
bool sym_tens::is_covariant(int i) {return type[i%2]==COVARIANT;}
bool sym_tens::is_contravariant(int i) {return type[i%2]==CONTRAVARIANT;}
void sym_vec::set_type(sym_vec_type type_set) {type=type_set;}
void sym_tens::set_type(sym_vec_type type_set0,sym_vec_type type_set1) {type[0]=type_set0;type[1]=type_set1;}

sym_vec sym_vec::set_variance(const sym_vec &v) const {

	return set_variance(v.type);
}
sym_vec sym_vec::set_variance(sym_vec_type new_type) const {

	sym_vec vnew;
	
	switch(new_type) {
		case COVARIANT:
			vnew=covariant(*this);
			break;
		case CONTRAVARIANT:
			vnew=contravariant(*this);
	}
	return vnew;

}


sym_vec sym_vec::operator+(const sym_vec &v) const {

	sym_vec vnew,v2;
	vnew.set_context(check_context(v));
	vnew.type=type;
	
	v2=v.set_variance(*this);
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]+v2.s[i];

	return vnew;

}

sym_vec sym_vec::operator-(const sym_vec &v) const {

	return (*this)+(-v);

}

sym_vec operator+(const sym_vec &v) {

	v.check_context();
	return v;
}
sym_vec operator-(const sym_vec &v) {

	sym_vec vnew;
	vnew.set_context(v.check_context());
	vnew.set_type(v.type);
	
	for(int i=0;i<3;i++) 
		vnew(i)=-v(i);
	return vnew;

}

sym_vec operator*(const sym &s,const sym_vec &v) {return v*s;}
sym_vec operator*(double x,const sym_vec &v) {return v*x;}

sym_vec sym_vec::operator*(const sym &q) const {
	
	sym_vec vnew;
	vnew.set_context(check_context(q));
	vnew.type=type;
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]*q;

	return vnew;

}

sym_vec sym_vec::operator*(double x) const {

	sym_vec vnew;
	vnew.set_context(check_context());
	vnew.type=type;
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]*x;

	return vnew;

}

sym_vec sym_vec::operator/(const sym &q) const {
	
	sym_vec vnew;
	vnew.set_context(check_context(q));
	vnew.type=type;
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]/q;

	return vnew;

}

sym_vec sym_vec::operator/(double x) const {

	sym_vec vnew;
	vnew.set_context(check_context());
	vnew.type=type;
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]/x;

	return vnew;

}

sym_vec cross(const sym_vec &v1,const sym_vec &v2) {

	sym_vec vnew,a,b;
	const symbolic *C;
	
	C=v1.check_context(v2);
	vnew.set_context(C);
	sym_vec_type q;
	
	q=v1.type;
	a=covariant(v1);
	b=covariant(v2);
	
	vnew.type=CONTRAVARIANT;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				vnew(i)=vnew(i)+C->perm(i,j,k)/C->sqrt_g*a(j)*b(k);
	vnew=vnew.set_variance(q);
	
	return vnew;

}

sym_tens tensor(const sym_vec &v1,const sym_vec &v2) {

	sym_tens tnew(v1.type,v2.type);

	tnew.set_context(v1.check_context(v2));
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
			tnew(i,j)=v1(i)*v2(j);
		}
	}
	return tnew;

}

sym_tens sym_tens::set_variance(const sym_tens &t) const {

	return set_variance(t.type[0],t.type[1]);

}

sym_tens sym_tens::set_variance(sym_vec_type t1,sym_vec_type t2) const {

	if(t1==COVARIANT) {
		if(t2==COVARIANT) return covariant_covariant(*this);
		else return covariant_contravariant(*this);
	} else {
		if(t2==COVARIANT) return contravariant_covariant(*this);
		else return contravariant_contravariant(*this);
	}
}


sym_tens sym_tens::operator+(const sym_tens &t) const {

	sym_tens tnew,t2;
	tnew.set_context(check_context(t));
	tnew.set_type(type[0],type[1]);
	
	t2=t.set_variance(*this);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) 
			tnew.s[i][j]=s[i][j]+t2.s[i][j];

	return tnew;

}

sym_tens sym_tens::operator-(const sym_tens &t) const {

	return (*this)+(-t);

}

sym_tens operator+(const sym_tens &t) {

	t.check_context();
	return t;
}
sym_tens operator-(const sym_tens &t) {

	sym_tens tnew;
	tnew.set_context(t.check_context());
	tnew.set_type(t.type[0],t.type[1]);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew(i,j)=-t(i,j);
	return tnew;

}

sym_tens operator*(const sym &s,const sym_tens &t) {return t*s;}
sym_tens operator*(double x,const sym_tens &t) {return t*x;}

sym_tens sym_tens::operator*(const sym &q) const {
	
	sym_tens tnew;
	tnew.set_context(check_context(q));
	tnew.type[0]=type[0];
	tnew.type[1]=type[1];
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew.s[i][j]=s[i][j]*q;

	return tnew;

}

sym_tens sym_tens::operator*(double x) const {
	
	sym_tens tnew;
	tnew.set_context(check_context());
	tnew.type[0]=type[0];
	tnew.type[1]=type[1];
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew.s[i][j]=s[i][j]*x;

	return tnew;

}


sym_tens sym_tens::operator/(const sym &q) const {
	
	sym_tens tnew;
	tnew.set_context(check_context(q));
	tnew.type[0]=type[0];
	tnew.type[1]=type[1];
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew.s[i][j]=s[i][j]/q;

	return tnew;

}

sym_tens sym_tens::operator/(double x) const {
	
	sym_tens tnew;
	tnew.set_context(check_context());
	tnew.type[0]=type[0];
	tnew.type[1]=type[1];
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew.s[i][j]=s[i][j]/x;

	return tnew;

}

sym sym_vec::operator,(const sym_vec &v) const {

	sym snew;
	snew.context=check_context(v);
	
	sym_vec v1,v2;
	v1=contravariant(*this);
	v2=covariant(v);
	
	for(int i=0;i<3;i++) 
		snew=snew+v1.s[i]*v2.s[i];

	return snew;

}

sym sym_tens::operator%(const sym_tens &t) const {

	sym snew;
	snew.context=check_context(t);
	
	sym_tens t1,t2;
	t1=contravariant_contravariant(*this);
	t2=covariant_covariant(t);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			snew=snew+t1.s[i][j]*t2.s[i][j];

	return snew;

}

sym_tens sym_tens::operator,(const sym_tens &t) const {

	sym_tens tnew;
	tnew.set_context(check_context(t));
	
	sym_vec_type q1,q2;
	q1=type[0];q2=t.type[1];
	tnew.set_type(q1,q2);
	
	sym_tens t1;
	q1=(type[1]==COVARIANT?CONTRAVARIANT:COVARIANT);
	t1=t.set_variance(q1,q2);
		
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				tnew.s[i][j]=tnew.s[i][j]+s[i][k]*t1.s[k][j];

	return tnew;

}

sym_vec sym_tens::operator,(const sym_vec &v) const {

	sym_vec vnew;
	vnew.set_context(check_context(v));
	
	vnew.type=type[0];
	
	sym_vec v1;
	v1=v.set_variance(type[1]==COVARIANT?CONTRAVARIANT:COVARIANT);
		
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			vnew.s[i]=vnew.s[i]+s[i][j]*v1.s[j];

	return vnew;

}


sym_vec sym_vec::operator,(const sym_tens &t) const {

	sym_vec vnew;
	vnew.set_context(check_context(t));
	
	vnew.type=t.type[1];
	
	sym_vec v1;
	v1=set_variance(t.type[0]==COVARIANT?CONTRAVARIANT:COVARIANT);
		
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			vnew.s[i]=vnew.s[i]+t.s[j][i]*v1.s[j];

	return vnew;

}

sym sym_vec::D(const sym &q) const {

	return (*this,grad(q));

}

sym_vec sym_vec::D(const sym_vec &v) const {

	sym_vec vnew(CONTRAVARIANT);
	const symbolic *C;
	C=check_context(v);
	vnew.set_context(C);
	
	sym_vec a,b;
	
	a=contravariant(*this);
	b=contravariant(v);
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			vnew(i)=vnew(i)+a(j)*C->covderiv(b,i,j);
			
	return vnew.set_variance(v);
	
}

sym sym::jacobian(const char *name,int dz,int dt) const {

	const symbolic *C;
	
	C=check_context();
	int ivar=C->id(name);
	
	if(ivar==-1) return jacobian_r(dz,dt);
	
	sym J;
	J.context=C;
	if(dz+dt>C->maxder) return J;

	
	for(symterm *t=terms;t!=NULL;t=t->next) {
		if(t->var[ivar][dz][dt]!=0) {
			symterm *tnew=J.addterm(*t);
			tnew->num*=tnew->var[ivar][dz][dt];
			tnew->var[ivar][dz][dt]--;
		}
	}
	
	return J;

}

sym sym::jacobian_r(int dz,int dt) const {

	const symbolic *C;
	sym J;
	C=check_context();
	J.context=C;
	if(dz+dt>C->maxder+1) return J;

	for(symterm *t=terms;t!=NULL;t=t->next) {
		if(t->r[dz][dt]!=0) {
			symterm *tnew=J.addterm(*t);
			tnew->num*=tnew->r[dz][dt];
			tnew->r[dz][dt]--;
		}
	}
	
	return J;

}

void sym::add(solver *op,const char *eq_name,const char *var_name) const {

	add(op,eq_name,var_name,ones(1,1));

}

void sym::add(solver *op,const char *eq_name,const char *var_name,const matrix &d) const {

	const symbolic *C;
	C=check_context();
	
	int maxder=C->maxder;
	int parity;
	if (!strcmp(var_name,"r")) {
		maxder++;
		parity=0;
	} else {
		parity=C->var_par[C->id(var_name)];
	}
	
	sym q;
	matrix Dt,Dt2;
	switch(parity) {
		case 00:
			Dt=C->map.leg.D_00;
			Dt2=C->map.leg.D2_00;
			break;
		case 01:
			Dt=C->map.leg.D_01;
			Dt2=C->map.leg.D2_01;
			break;
		case 10:
			Dt=C->map.leg.D_10;
			Dt2=C->map.leg.D2_10;
			break;
		case 11:
			Dt=C->map.leg.D_11;
			Dt2=C->map.leg.D2_11;
	}
	
	matrix_block_diag pre;
	matrix post;
	for(int i=0;i<maxder+1;i++)
		for(int j=0;j<maxder+1;j++) {
			q=jacobian(var_name,i,j);
			if(q.terms!=NULL) {
				pre=eye(C->map.D);
				post=eye(C->map.leg.npts);
				for(int k=0;k<i;k++) pre=(C->map.D,pre);
				for(int k=0;k<j/2;k++) post=(post,Dt2);
				if(j%2) post=(post,Dt);
				if(i==0&&j==0) op->add_d(eq_name,var_name,d*q.eval());
				else if(j==0) op->add_l(eq_name,var_name,d*q.eval(),pre);
				else if(i==0) op->add_r(eq_name,var_name,d*q.eval(),post);
				else op->add_lr(eq_name,var_name,d*q.eval(),pre,post);
			}
		}

}

void sym::add_bc(solver *op,int n,const char *type,const char *eq_name,const char *var_name,const matrix &d) const {

	const symbolic *C;
	C=check_context();
	
	int maxder=C->maxder;
	int parity;
	if (!strcmp(var_name,"r")) {
		maxder++;
		parity=0;
	} else {
		parity=C->var_par[C->id(var_name)];
	}
	
	sym q;
	matrix Dt,Dt2;
	switch(parity) {
		case 00:
			Dt=C->map.leg.D_00;
			Dt2=C->map.leg.D2_00;
			break;
		case 01:
			Dt=C->map.leg.D_01;
			Dt2=C->map.leg.D2_01;
			break;
		case 10:
			Dt=C->map.leg.D_10;
			Dt2=C->map.leg.D2_10;
			break;
		case 11:
			Dt=C->map.leg.D_11;
			Dt2=C->map.leg.D2_11;
	}
	
	int j0=0;
	for(int i=0;i<n;i++) {
		j0+=C->map.gl.npts[i];
	}
	if(!strcmp(type,"top1")) j0+=C->map.gl.npts[n]-1;
	if(!strcmp(type,"top2")) j0+=C->map.gl.npts[n];
	if(!strcmp(type,"bot1")) j0--;
	
	void (solver::*add_d) (int,const char *,const char *,const matrix &);
	void (solver::*add_l) (int,const char *,const char *,const matrix &,const matrix &);
	void (solver::*add_r) (int,const char *,const char *,const matrix &,const matrix &);
	void (solver::*add_lr) (int,const char *,const char *,const matrix &,const matrix &,const matrix &);
	
	if(!strcmp(type,"top1")) {
		add_d=&solver::bc_top1_add_d;
		add_l=&solver::bc_top1_add_l;
		add_r=&solver::bc_top1_add_r;
		add_lr=&solver::bc_top1_add_lr;
	}
	if(!strcmp(type,"top2")) {
		add_d=&solver::bc_top2_add_d;
		add_l=&solver::bc_top2_add_l;
		add_r=&solver::bc_top2_add_r;
		add_lr=&solver::bc_top2_add_lr;
	}
	if(!strcmp(type,"bot1")) {
		add_d=&solver::bc_bot1_add_d;
		add_l=&solver::bc_bot1_add_l;
		add_r=&solver::bc_bot1_add_r;
		add_lr=&solver::bc_bot1_add_lr;
	}
	if(!strcmp(type,"bot2")) {
		add_d=&solver::bc_bot2_add_d;
		add_l=&solver::bc_bot2_add_l;
		add_r=&solver::bc_bot2_add_r;
		add_lr=&solver::bc_bot2_add_lr;
	}
	
	matrix pre,D;
	if(!strcmp(type,"top1")) D=C->map.D.block(n);
	if(!strcmp(type,"top2")) D=C->map.D.block(n+1);
	if(!strcmp(type,"bot1")) D=C->map.D.block(n-1);
	if(!strcmp(type,"bot2")) D=C->map.D.block(n);
				
	matrix post;
	for(int i=0;i<maxder+1;i++)
		for(int j=0;j<maxder+1;j++) {
			q=jacobian(var_name,i,j);
			if(q.terms!=NULL) {
				pre=eye(C->map.gl.npts[n]);
				for(int k=0;k<i;k++) pre=(D,pre);
				if(!strcmp(type,"top1")) pre=pre.row(-1);
				if(!strcmp(type,"top2")) pre=pre.row(0);
				if(!strcmp(type,"bot1")) pre=pre.row(-1);
				if(!strcmp(type,"bot2")) pre=pre.row(0);
				post=eye(C->map.leg.npts);
				for(int k=0;k<j/2;k++) post=(post,Dt2);
				if(j%2) post=(post,Dt);
				if(i==0&&j==0) (op->*add_d)(n,eq_name,var_name,d*q.eval().row(j0));
				else if(j==0) (op->*add_l)(n,eq_name,var_name,d*q.eval().row(j0),pre);
				else if(i==0) (op->*add_r)(n,eq_name,var_name,d*q.eval().row(j0),post);
				else (op->*add_lr)(n,eq_name,var_name,d*q.eval().row(j0),pre,post);
			}
		}


}

