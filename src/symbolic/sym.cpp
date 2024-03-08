#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "symbolic.h"

#include <cstdlib>
#include <typeinfo>

using namespace std;

sym::sym() {
	expr=new sym_num(0);
	context=NULL;
}
sym::~sym() {
	if(expr!=NULL) delete expr;
}

sym::sym(const sym &s) {
	context=s.context;
	expr=s.expr->clone();

}

sym::sym(const double &q) {
	context=NULL;
	expr=new sym::sym_num(q);
}

sym::sym(const int &q) {
	context=NULL;
	expr=new sym::sym_num(q);
}

sym::sym(const rational &q) {
	context=NULL;
	expr=new sym::sym_num(q.eval());
}

sym &sym::operator=(const sym &s) {
	context=s.context;
	delete expr;
	expr=s.expr->clone();
	return *this;
}

ostream& operator<<(ostream& os, const sym &s) {

	return os<<*s.expr;

}

sym trig_simplify(const sym &s) {
	sym snew(s);
	bool temp=symbolic::trig_simplify;
	symbolic::trig_simplify=true;
	snew.expr=snew.expr->reduce();
	symbolic::trig_simplify=temp;
	return snew;
}
sym_vec trig_simplify(const sym_vec &v) {
	sym_vec w;
	for(int i=0;i<3;i++) w(i)=trig_simplify(v(i));
	return w;
}
sym_tens trig_simplify(const sym_tens &v) {
	sym_tens w;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) w(i,j)=trig_simplify(v(i,j));
	return w;
}

symbolic *sym::check_context() const {

	if(context==NULL) {
		cerr<<"Symbolic: sym object not initialized"<<endl;
		exit(1);
	}
	return context;
}

symbolic *sym::check_context(const sym &s) const {

	symbolic *c;

	bool allow_null_1=typeid(*expr)==typeid(sym::sym_num);
	bool allow_null_2=typeid(*(s.expr))==typeid(sym::sym_num);

	if(this->context==NULL&&s.context==NULL&&allow_null_1&&allow_null_2) {
		return NULL;
	}
	if(this->context==NULL&&allow_null_1) {
		return s.check_context();
	}
	if(s.context==NULL&&allow_null_2) {
		return check_context();
	}

	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context"<<endl;
		exit(1);
	}
	return c;
}

symbolic *sym::check_context(const sym_vec &s) const {

	symbolic *c;

	bool allow_null=typeid(*expr)==typeid(sym::sym_num);

	if(context==NULL&&allow_null) {
		return s.check_context();
	}

	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context"<<endl;
		exit(1);
	}
	return c;
}

symbolic *sym::check_context(const sym_tens &s) const {

	symbolic *c;

	bool allow_null=typeid(*expr)==typeid(sym::sym_num);

	if(context==NULL&&allow_null) {
		return s.check_context();
	}

	c=check_context();

	if(c!=s.check_context()) {
		cerr<<"Symbolic: Wrong context"<<endl;
		exit(1);
	}
	return c;
}

sym operator+(const sym &s1,const sym &s2) {

	sym snew;

	snew.context=s1.check_context(s2);
	delete snew.expr;
	snew.expr=s1.expr->clone();
	snew.expr=snew.expr->add(*s2.expr)->reduce();

	return snew;
}

sym operator*(const sym &s1,const sym &s2) {

	sym snew;

	snew.context=s1.check_context(s2);
	delete snew.expr;
	snew.expr=s1.expr->clone();
	snew.expr=snew.expr->mult(*s2.expr)->reduce();

	return snew;
}

sym pow(const sym &s,const rational &q) {

	if(typeid(*s.expr)==typeid(sym::sym_num))
		return sym(pow(((sym::sym_num *)s.expr)->value,q.eval()));

	sym snew;

	snew.context=s.check_context();
	delete snew.expr;
	snew.expr=s.expr->clone();
	snew.expr=snew.expr->pow(q)->reduce();

	return snew;
}

sym sin(const sym &s) {

	sym snew;

	snew.context=s.check_context();
	delete snew.expr;
	snew.expr=s.expr->clone();
	snew.expr=snew.expr->sin()->reduce();

	return snew;
}

sym cos(const sym &s) {

	sym snew;

	snew.context=s.check_context();
	delete snew.expr;
	snew.expr=s.expr->clone();
	snew.expr=snew.expr->cos()->reduce();

	return snew;
}

sym tan(const sym &s) {

	return sin(s)/cos(s);
}

sym exp(const sym &s) {

	sym snew;

	snew.context=s.check_context();
	delete snew.expr;
	snew.expr=s.expr->clone();
	snew.expr=snew.expr->exp()->reduce();

	return snew;
}

sym log(const sym &s) {

	sym snew;

	snew.context=s.check_context();
	delete snew.expr;
	snew.expr=s.expr->clone();
	snew.expr=snew.expr->log()->reduce();

	return snew;
}

sym pow(const sym &s,const sym &ex) {

	return exp(ex*log(s));

}

sym pow(const double &n,const sym &ex) {
	return exp(ex*log(n));
}

sym diff(const sym &f,const sym &x) {

	if(typeid(*x.expr)==typeid(sym::symbol)) {
		if(((sym::symbol *)x.expr)->is_indep) {
			return jacobian(f,x);
		}
	}

	ester_critical("(symbolic) Can derive only respect to independent symbols");
}
sym jacobian(const sym &f,const sym &a) {

	if(typeid(*a.expr)!=typeid(sym::symbol)&&typeid(*a.expr)!=typeid(sym::sym_deriv)) {
		ester_critical("(symbolic) Can calculate jacobian only with respect to a symbol or a derivative of a symbol");
	}

	sym snew;
	snew.context=f.check_context(a);
	delete snew.expr;
	snew.expr=f.expr->clone();
	snew.expr=snew.expr->derive(*a.expr)->reduce();

	return snew;
}

bool sym::operator==(const sym &s) const {return *expr==*s.expr;}

matrix sym::eval() const {
	return expr->eval();
}

void sym::add(solver *op,std::string eq_name,std::string var_name) const
		{check_context()->add(*this,op,0,"in",eq_name,var_name,ones(1,1));}
void sym::add(solver *op,std::string eq_name,std::string var_name,const matrix &d) const
		{check_context()->add(*this,op,0,"in",eq_name,var_name,d);}
void sym::add_ex(solver *op,int n,std::string eq_name,std::string var_name) const
		{check_context()->add(*this,op,n,"ex",eq_name,var_name,ones(1,1));}
void sym::add_ex(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const
		{check_context()->add(*this,op,n,"ex",eq_name,var_name,d);}

void sym::bc_top1_add(solver *op,int n,std::string eq_name,std::string var_name) const
		{check_context()->add_bc(*this,op,n,"top1",eq_name,var_name,ones(1,1));}
void sym::bc_top2_add(solver *op,int n,std::string eq_name,std::string var_name) const
		{check_context()->add_bc(*this,op,n,"top2",eq_name,var_name,ones(1,1));}
void sym::bc_bot1_add(solver *op,int n,std::string eq_name,std::string var_name) const
		{check_context()->add_bc(*this,op,n,"bot1",eq_name,var_name,ones(1,1));}
void sym::bc_bot2_add(solver *op,int n,std::string eq_name,std::string var_name) const
		{check_context()->add_bc(*this,op,n,"bot2",eq_name,var_name,ones(1,1));}
void sym::bc_top1_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const
		{check_context()->add_bc(*this,op,n,"top1",eq_name,var_name,d);}
void sym::bc_top2_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const
		{check_context()->add_bc(*this,op,n,"top2",eq_name,var_name,d);}
void sym::bc_bot1_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const
		{check_context()->add_bc(*this,op,n,"bot1",eq_name,var_name,d);}
void sym::bc_bot2_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const
		{check_context()->add_bc(*this,op,n,"bot2",eq_name,var_name,d);}

std::ostream& sym::print(std::ostream& os) const {
    this->expr->print(os);
    return os;
}
