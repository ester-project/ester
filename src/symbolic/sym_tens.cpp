#include "ester-config.h"
#include "symbolic.h"

#include <cstdlib>
#include <typeinfo>
#include <iostream>

using namespace std;

sym_tens::sym_tens(sym_vec_type type_set0,sym_vec_type type_set1) {
	
	type[0]=type_set0;
	type[1]=type_set1;

}

void sym_tens::set_context(symbolic *context) {

	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) 
			s[i][j].context=context;

}

symbolic *sym_tens::check_context() const {

	symbolic *context;

	context=s[0][0].check_context();
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) 
			if(context!=s[i][j].check_context()) {
				cerr<<"Symbolic: sym_tens components have different contexts\n";
				exit(1);
			}

	return context;
}

symbolic *sym_tens::check_context(const sym &s) const {
	
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

symbolic *sym_tens::check_context(const sym_vec &s) const {
	
	symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context\n";
		exit(1);
	}
	return c;
}

symbolic *sym_tens::check_context(const sym_tens &s) const {
	
	symbolic *c;
	
	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context";
		exit(1);
	}
	return c;
}

bool sym_tens::is_covariant(int i) {return type[i%2]==COVARIANT;}
bool sym_tens::is_contravariant(int i) {return type[i%2]==CONTRAVARIANT;}
void sym_tens::set_type(sym_vec_type type_set0,sym_vec_type type_set1) {type[0]=type_set0;type[1]=type_set1;}

sym_tens sym_tens::set_variance(const sym_tens &t) const {return set_variance(t.type[0],t.type[1]);}

sym_tens sym_tens::set_variance(sym_vec_type t1,sym_vec_type t2) const {

	if(t1==COVARIANT) {
		if(t2==COVARIANT) return covariant_covariant(*this);
		else return covariant_contravariant(*this);
	} else {
		if(t2==COVARIANT) return contravariant_covariant(*this);
		else return contravariant_contravariant(*this);
	}
}


sym &sym_tens::operator()(int i,int j) {
	
	if(i<0||i>2||j<0||j>2) {
		cerr<<"Symbolic: sym_tens index out of range\n";
		exit(1);
	}
	return s[i][j];
	
}


sym sym_tens::operator()(int i,int j) const {
	
	if(i<0||i>2||j<0||j>2) {
		cerr<<"Symbolic: sym_tens index out of range\n";
		exit(1);
	}
	return s[i][j];
	
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

sym_tens operator-(const sym_tens &t) {

	sym_tens tnew;
	tnew.set_context(t.check_context());
	tnew.set_type(t.type[0],t.type[1]);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew(i,j)=-t(i,j);
	return tnew;

}

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

sym_tens sym_tens::T() const {

	sym_tens tnew(this->type[1],this->type[0]);
	
	tnew.set_context(check_context(*this));
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			tnew(i,j)=(*this)(j,i);
	
	return tnew;

}

std::ostream& operator<<(std::ostream& os, const sym_tens &v) {

	os<<v(0,0)<<" \t"<<v(0,1)<<" \t"<<v(0,2)<<endl;
	os<<v(1,0)<<" \t"<<v(1,1)<<" \t"<<v(1,2)<<endl;
	os<<v(2,0)<<" \t"<<v(2,1)<<" \t"<<v(2,2)<<endl;

	return os;

}



