#include "ester-config.h"
#include "symbolic.h"

#include <cstdlib>
#include <typeinfo>
#include <iostream>

using namespace std;

sym_vec::sym_vec(sym_vec_type type_set) {
	
	type=type_set;

}

void sym_vec::set_context(symbolic *context) {

	for(int i=0;i<3;i++) 
		s[i].context=context;

}

symbolic *sym_vec::check_context() const {

	symbolic *context;

	context=s[0].check_context();
	for(int i=1;i<3;i++) 
		if(context!=s[i].check_context()) {
			cerr<<"Symbolic: sym_vec components have different contexts\n";
			exit(1);
		}

	return context;
}

symbolic *sym_vec::check_context(const sym &s) const {
	
	symbolic *c;
	
	bool allow_null=typeid(*(s.expr))==typeid(sym::sym_num);

	if(s.context==NULL&&allow_null) {
		return check_context();
	}
	
	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context"<<endl;
		exit(1);
	}
	return c;
}

symbolic *sym_vec::check_context(const sym_vec &s) const {
	
	symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context\n";
		exit(1);
	}
	return c;
}

symbolic *sym_vec::check_context(const sym_tens &s) const {
	
	symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context\n";
		exit(1);
	}
	return c;
}

bool sym_vec::is_covariant() {return type==COVARIANT;}
bool sym_vec::is_contravariant() {return type==CONTRAVARIANT;}
void sym_vec::set_type(sym_vec_type type_set) {type=type_set;}

sym_vec sym_vec::set_variance(const sym_vec &v) const {return set_variance(v.type);}

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

sym &sym_vec::operator()(int i) {
	
	if(i<0||i>2) {
		cerr<<"Symbolic: sym_vec index out of range\n";
		exit(1);
	}
	return s[i];
	
}


sym sym_vec::operator()(int i) const {
	
	if(i<0||i>2) {
		cerr<<"Symbolic: sym_vec index out of range\n";
		exit(1);
	}
	return s[i];
	
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

sym_vec sym_vec::operator+(const sym_vec &v) const {

	sym_vec vnew,v2;
	vnew.set_context(check_context(v));
	vnew.type=type;
	
	v2=v.set_variance(*this);
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]+v2.s[i];

	return vnew;

}

sym_vec operator-(const sym_vec &v) {

	sym_vec vnew;
	vnew.set_context(v.check_context());
	vnew.set_type(v.type);
	
	for(int i=0;i<3;i++) 
		vnew(i)=-v(i);
	return vnew;

}

sym_vec sym_vec::operator*(const sym &q) const {
	
	sym_vec vnew;
	vnew.set_context(check_context(q));
	vnew.type=type;
	
	for(int i=0;i<3;i++) 
		vnew.s[i]=s[i]*q;

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

sym_vec cross(const sym_vec &v1,const sym_vec &v2) {

	sym_vec vnew,a,b;
	symbolic *C;
	
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

sym sym_vec::D(const sym &q) const {

	return (*this,grad(q));

}

sym_vec sym_vec::D(const sym_vec &v) const {

	return (grad(v),*this);
	
}

std::ostream& operator<<(std::ostream& os, const sym_vec &v) {

	os<<v(0)<<endl;
	os<<v(1)<<endl;
	os<<v(2)<<endl;

	return os;

}

