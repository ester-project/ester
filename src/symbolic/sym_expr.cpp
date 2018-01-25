#include "ester-config.h"
#include"symbolic.h"
#include<iostream>
#include<typeinfo>
#include <algorithm> 
#include <cmath> 

using namespace std;

#define COMPARE(A,B) ((A)>(B)?1:((A)==(B)?0:-1))

sym::sym_flags::sym_flags() {
	collect = false;
	trig_simplify = false;
}

bool sort_pair_d(const pair<sym::sym_expr *,double> &a,const pair<sym::sym_expr *,double> &b) {

	int c;
	c=a.first->comp(*(b.first));
	if(c==-1) return true;
	if(c==1) return false;
	return a.second < b.second;

}

bool sort_pair_r(const pair<sym::sym_expr *,rational> &a,const pair<sym::sym_expr *,rational> &b) {

	int c;
	c=a.first->comp(*(b.first));
	if(c==-1) return true;
	if(c==1) return false;
	return a.second < b.second;

}

std::ostream& operator<<(std::ostream& os, const sym::sym_expr &s) {

	return s.print(os);
}

sym::sym_expr *sym::sym_expr::add(const sym_expr &s) {
	return sym_add::create(this,s.clone());
}

sym::sym_expr *sym::sym_expr::mult(const sym_expr &s) {
	return sym_prod::create(this,s.clone());
}

sym::sym_expr *sym::sym_expr::pow(const rational &q) {
	return sym_prod::create_pow(this,q);
}

sym::sym_expr *sym::sym_expr::sin() {
	return sym_sin::create(this);
}

sym::sym_expr *sym::sym_expr::cos() {
	return sym_cos::create(this);
}

sym::sym_expr *sym::sym_expr::exp() {
	return sym_exp::create(this);
}

sym::sym_expr *sym::sym_expr::log() {
	return sym_log::create(this);
}

////////////////////////// sym_num ////////////////////////////////

sym::sym_num::sym_num(const double &q) {
	value=q;
}

sym::sym_num *sym::sym_num::clone() const {
	return new sym_num(*this);
}

int sym::sym_num::nodeCount() const {
	return 1;
}

int sym::sym_num::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_num *q;
	q=(sym_num *) &s;
	
	return COMPARE(value,q->value);

}

sym::sym_expr *sym::sym_num::reduce(sym_flags flags) {

	value=symbolic::round_to_tol(value);
	
	return this;

}

sym::sym_expr *sym::sym_num::derive(const sym_expr &) {
	value=0;
	return this;
}

matrix sym::sym_num::eval() const {
	return value*ones(1,1);
} 

ostream &sym::sym_num::print(ostream &os) const {

	if(value<0) os<<"(";
	os<<value;
	if(value<0) os<<")";
	return os;
}


//////////////////////////////// symbol /////////////////////////////////

sym::symbol *sym::symbol::clone() const {
	return new symbol(*this);
}

int sym::symbol::nodeCount() const {
	return 1;
}

int sym::symbol::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	symbol *q;
	q=(symbol *) &s;
	
	return COMPARE(name,q->name);

}

sym::sym_expr *sym::symbol::derive(const sym_expr &s) {
	
	if(is_indep || is_const) {
		if(s==*this) {delete this;return new sym_num(1);}
		else {delete this;return new sym_num(0);}
	}
	
	if(typeid(s)==typeid(symbol)) {
		symbol *symb;
		symb=(symbol *)(&s);
		if(symb->is_indep) {
			sym_expr *s_spec;
			s_spec=context->derive(*this,*symb);
			if(s_spec!=NULL) {
				delete this;
				return s_spec;
			}
			sym_deriv *sderiv=sym_deriv::create(this,*symb);
			sderiv->context=context;
			return sderiv;
		}
	}
	
	if(s==*this) {delete this;return new sym_num(1);}
	else {delete this;return new sym_num(0);}
	
}

matrix sym::symbol::eval() const {
	return context->get_value(*this);
} 

ostream &sym::symbol::print(ostream &os) const {
	return os<<name;
}

//////////////////////////////// sym_deriv /////////////////////////////////

sym::sym_deriv::~sym_deriv() {
	delete oper;
}

sym::sym_deriv::sym_deriv(const sym_deriv &s) {
	context=s.context;
	oper=s.oper->clone();
	var=s.var;
}

sym::sym_deriv *sym::sym_deriv::clone() const {
	return new sym_deriv(*this);
}

int sym::sym_deriv::nodeCount() const {
	int cnt = 1;
	cnt += oper->nodeCount();
	return cnt;
}

int sym::sym_deriv::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_deriv *q;
	q=(sym_deriv *) &s;
	
	int c=oper->comp(*(q->oper));
	if(c!=0) return c;
	
	return var.comp(q->var);

}

sym::sym_expr *sym::sym_deriv::reduce(sym_flags flags) {

	oper->reduce(flags);
	if(typeid(*oper)==typeid(sym_deriv)) {
		sym_deriv *sderiv;
		sderiv=(sym_deriv *)oper;
		if(sderiv->var.name<var.name) {
			symbol stemp;
			stemp=var;
			var=sderiv->var;
			sderiv->var=stemp;
			sderiv->reduce(flags);
		}
	}
	return this;

}

sym::sym_deriv *sym::sym_deriv::create(sym_expr *s1,const symbol &symb) {

	sym_deriv *snew;
	
	snew=new sym_deriv();
	snew->oper=s1;
	snew->var=symb;
	return snew;
}


sym::sym_expr *sym::sym_deriv::derive(const sym_expr &s) {
	
	if(typeid(s)==typeid(symbol)) {
		symbol *symb;
		symb=(symbol *)(&s);
		if(symb->is_indep) {
			sym_expr *s_spec;
			s_spec=context->derive(*this,*symb);
			if(s_spec!=NULL) {
				delete this;
				return s_spec;
			}
			sym_deriv *sderiv=sym_deriv::create(this,*symb);
			sderiv->context=context;
			return sderiv;
		}
	}
	
	if(s==*this) {delete this;return new sym_num(1);}
	else {delete this;return new sym_num(0);}
	
}

matrix sym::sym_deriv::eval() const {
	return context->get_value(*this);
} 

ostream &sym::sym_deriv::print(ostream &os) const {

	os<<"D("<<*oper<<","<<var<<")";	
	
	return os;
}


///////////////////////// sym_add //////////////////////////////////

sym::sym_add::~sym_add() {

	for(unsigned int i=0;i<oper.size();i++) {
		delete oper[i].first;
	}

}

sym::sym_add::sym_add(const sym_add &s) {

	oper=s.oper;
	for(unsigned int i=0;i<oper.size();i++) 
		oper[i].first=s.oper[i].first->clone();

}

sym::sym_add *sym::sym_add::clone() const {
	return new sym_add(*this);
}

int sym::sym_add::nodeCount() const {
	int cnt = 1;
	for(unsigned int i=0; i<oper.size(); i++)
		cnt += oper[i].first->nodeCount();
	return cnt;
}

int sym::sym_add::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_add *q;
	q=(sym_add *) &s;
	
	int c;
	c=COMPARE(oper.size(),(q->oper).size());
	if(c!=0) return c;

	for(unsigned int i=0;i<oper.size();i++) {
		c=oper[i].first->comp( *(q->oper[i].first) );
		if(c!=0) return c;
		c=COMPARE(oper[i].second,q->oper[i].second);
		if(c!=0) return c;
	}
	return c;
}

sym::sym_expr *sym::sym_add::reduce(sym_flags flags) {
	
	sym_add s_old(*this);

	for(unsigned int i=0;i<oper.size();i++)
        oper[i].first=oper[i].first->reduce(flags);
	for(unsigned int i=0;i<oper.size();i++)
        oper[i].second=symbolic::round_to_tol(oper[i].second);
	
	
// Find children nodes of type sym_add and merge them in current node
//  a + [b + c] --> a + b + c

	int n=oper.size();
	for(int	i=0;i<n;i++) {
		if(typeid(*(oper[i].first))==typeid(sym_add)) {
			sym_add *s;
			double q;
			s=(sym_add *)oper[i].first;
			q=oper[i].second;
			oper[i].first=s->oper[0].first;
			oper[i].second=q*s->oper[0].second;
			s->oper[0].first=NULL;
			for(unsigned int j=1;j < s->oper.size();j++) {
				oper.push_back(make_pair(s->oper[j].first,q*s->oper[j].second));
				s->oper[j].first=NULL;
			}
			delete s;
		}
	}

// Find children nodes of type sym_prod with nodes of type sym_num that can be
// merged into the coefficient within current node
//  a + 1*[2*b*c] --> a + 2*[1*b*c]

	n=oper.size();
	for(int	i=0;i<n;i++) {
		if(typeid(*(oper[i].first))==typeid(sym_prod)) {
			sym_prod *s;
			s=(sym_prod *)oper[i].first;
			for(unsigned int j=0;j<s->oper.size();j++) {
				if(typeid(*(s->oper[j].first))==typeid(sym_num)) {
					sym_num *snum;
					snum=(sym_num *)s->oper[j].first;
					oper[i].second*=snum->value; // The exponent s->oper[j].second should be 1 if the sym_prod
													// node has been reduced
					snum->value=1;		// It will be deleted by the next reduction of the child node
					break;    // After reduction of the sym_prod node ther will be only 1 sym_num node
				}
			}			
		}
	}

// sym_num nodes should have coefficient=1
	
	n=oper.size();
	for(int	i=0;i<n;i++) {
		if(typeid(*(oper[i].first))==typeid(sym_num)) 
			if(oper[i].second!=1) {
				sym_num *s;
				s=(sym_num *)oper[i].first;
				s->value*=oper[i].second;
				oper[i].second=1;
			}	
	}

	if(flags.trig_simplify) {
		n=oper.size();
		for(int i=0;i<n;i++) {
			// sin(x)^2 or cos(x)^2 should be in a sym_prod node
			if(typeid(*(oper[i].first))==typeid(sym_prod)&&oper[i].second!=0) {	
				sym_prod *sprod;
				sprod=(sym_prod *)oper[i].first;
				// Try to find a sin or cos in the sym_prod node with exponent at least 2
				for(unsigned int j=0;j<sprod->oper.size();j++) {
					bool found=false;
					sym_expr *test,*snew;
					if(typeid(*(sprod->oper[j].first))==typeid(sym_sin)&&sprod->oper[j].second>=2) {
						found=true;
						test=sprod->clone();
						((sym_prod *)test)->oper[j].second-=2;
						test=test->mult(sym_num(oper[i].second));
						snew=test->clone();
						sym_expr* scos;
						scos=((sym_sin *)sprod->oper[j].first)->oper->clone(); //argument of sin
						scos=sym_cos::create(scos);
						scos=sym_prod::create_pow(scos,2);
						test=test->mult(*scos)->reduce(flags);
						delete scos;
					}
					else if(typeid(*(sprod->oper[j].first))==typeid(sym_cos)&&sprod->oper[j].second>=2) {
						found=true;
						test=sprod->clone();
						((sym_prod *)test)->oper[j].second-=2;
						test=test->mult(sym_num(oper[i].second));
						snew=test->clone();
						sym_expr* ssin;
						ssin=((sym_cos *)sprod->oper[j].first)->oper->clone(); //argument of cos
						ssin=sym_sin::create(ssin);
						ssin=sym_prod::create_pow(ssin,2);
						test=test->mult(*ssin)->reduce(flags);
						delete ssin;
					}
					if(found) {
						found=false;
						// Check if another node in the original sym_add node is = to test
						for(int k=0;k<n;k++) {
							if(k==i) continue;
							if(oper[k].second==0) continue;
							sym_expr *test2;
							test2=oper[k].first->clone();
							test2=test2->mult(sym_num(oper[k].second))->reduce(flags);
							if(*test2==*test) {
								found=true;
								delete oper[i].first;
								oper[i].first=snew->clone()->reduce(flags);
								oper[i].second=1;
								oper[k].second=0;
								delete test2;
								break;
							}
							delete test2;
						}
						delete test;
						delete snew;
						if(found) break;
					}
				}
			}
		} 
	}


// Sort
	sort(oper.begin(),oper.end(),sort_pair_d);

// Add equal terms

	n=oper.size();
	for(int	i=n-1;i>0;i--) {
		if(typeid(*(oper[i].first))==typeid(sym_num)) {
			if(typeid(*(oper[i-1].first))==typeid(sym_num)) { // Add numbers
				sym_num *s,*s0;
				s=(sym_num *)oper[i].first;
				s0=(sym_num *)oper[i-1].first;
				s0->value+=s->value;
				s->value=0;
			}
		} else {
			if(*(oper[i].first)==*(oper[i-1].first)) { // Add other expressions
				oper[i-1].second+=oper[i].second;
				oper[i].second=0;
			}
		}
		
	}

// Remove zeros
	for(unsigned int i=0;i<oper.size();i++) {
		if(oper[i].second==0) {
			delete oper[i].first;
			oper.erase(oper.begin()+i--);
		} else if(typeid(*(oper[i].first))==typeid(sym_num)) {
			if( ((sym_num *)oper[i].first)->value==0) {
				delete oper[i].first;
				oper.erase(oper.begin()+i--);
			} 		
		}
	}

// Check the new size and change to appropriate node type
	if(oper.size()==0) {
		delete this;
		return new sym_num(0);
	}
	if(oper.size()==1) {
		sym_expr *s;
		if(oper[0].second==1)
			s=oper[0].first;
		else {
			sym_num *f;
			f=new sym_num(oper[0].second);
			s=sym_prod::create(f,oper[0].first);
		}
		oper[0].first=NULL;
		delete this;
		return s->reduce(flags);
	}

// Collect terms with common factors
	if(flags.collect) {
		std::vector<sym_expr *> factors;
		// Collect possible factors in factors[]
		for(unsigned int i=0; i<oper.size(); i++) {
			if(typeid(*(oper[i].first))==typeid(sym_prod)) {
				sym_prod *sprod = (sym_prod *) oper[i].first;
				for(unsigned int j=0; j<sprod->oper.size(); j++) {
					bool already = false;
					for(unsigned int k=0; k<factors.size(); k++) {
						if (*(sprod->oper[j].first) == *(factors[k])) {
							already = true;
							break;
						}
					}
					if(!already)
						factors.push_back(sprod->oper[j].first->clone());
				}
			}
			else {
				bool already = false;
				for(unsigned int k=0; k<factors.size(); k++) {
					if (*(oper[i].first) == *(factors[k])) {
						already = true;
						break;
					}
					if(!already)
						factors.push_back(oper[i].first->clone());
				}
			}
		}
		// Create common factor
		// The factor is chosen if the minimum exponent among all the terms is != 0
		// Factors with negative exponent in at least 1 term are always chosen
		sym_prod *fact = new sym_prod();
		sym_prod *fact_inv = new sym_prod();
		for(unsigned int i=0; i<factors.size(); i++) {
			rational minexp = 999999999;
			for(unsigned int j=0; j<oper.size(); j++) {
				if(typeid(*(oper[j].first))==typeid(sym_prod)) {
					rational exponent = ((sym_prod *) oper[j].first)->get_exponent(*factors[i]);
					minexp = exponent < minexp ? exponent : minexp;
				}
				else if(*(oper[j].first) == *factors[i]) {
					minexp = minexp > 1 ? 1 : minexp;
				}
				else {
					minexp = minexp > 0 ? 0 : minexp;
				}
			}
			if(minexp != 0) {
				fact->oper.push_back(make_pair(factors[i]->clone(), minexp));
				fact_inv->oper.push_back(make_pair(factors[i]->clone(), -minexp));
			}
		}
		for(unsigned int i=0; i<factors.size(); i++)
			delete factors[i];
		if(fact->oper.size()) { // Extract global factor (fact)
			// Divide all terms by fact_inv
			for(unsigned int i=0; i<oper.size(); i++) {
				sym_prod *snew = sym_prod::create(fact_inv->clone(), oper[i].first);
				oper[i].first = snew->reduce(flags);
			}
			// Multiply globally bu fact
			sym_prod *snew = sym_prod::create(fact, this);
			return snew->reduce(flags);
		}
	}



	if(comp(s_old)) return reduce(flags);
	else return this;
	
}

sym::sym_expr *sym::sym_add::derive(const sym_expr &s) {

	for(unsigned int i=0;i<oper.size();i++)
		oper[i].first=oper[i].first->derive(s);

	return this;
}

matrix sym::sym_add::eval() const {
	
	matrix m=zeros(1,1);
	
	for(unsigned int i=0;i<oper.size();i++) {
		m+=oper[i].second*oper[i].first->eval();
	}
	return m;
} 

ostream &sym::sym_add::print(ostream &os) const {

	os<<"(";
	for(unsigned int i=0;i<oper.size();i++) {
		if(oper[i].second>=0&&i) os<<"+";
		if(oper[i].second!=1) {
			if(oper[i].second==-1) 
				os<<"-";
			else
				os<<oper[i].second<<"*";
		}
		os<<*(oper[i].first);
	}
	os<<")";
	return os;
}

sym::sym_add *sym::sym_add::create(sym_expr *s1,sym_expr *s2) {

	sym_add *snew;
	
	snew=new sym_add();
	snew->oper.push_back(make_pair(s1,1));
	snew->oper.push_back(make_pair(s2,1));

	return snew;
}

sym::sym_expr *sym::sym_add::multiply(const sym_add &s) {

	sym_add *snew;
	snew=new sym_add();
	for(unsigned int i=0;i<oper.size();i++) {
		for(unsigned int j=0;j<s.oper.size();j++) {
			snew->oper.push_back(make_pair(
				sym_prod::create(oper[i].first->clone(),s.oper[j].first->clone()),
				oper[i].second*s.oper[j].second));
		}
	}
	return snew;

}

sym::sym_expr *sym::sym_add::power(int n, sym_flags flags) {

	sym_add *sa,*sb,*stemp;

	sb=new sym_add();
	sb->oper.push_back(make_pair(new sym_num(1),1));
	sa=clone();
	while(n>1) {
		if(n&1) { // If n odd
			stemp=(sym_add *)sb->multiply(*sa);
			delete sb;
			sb=stemp;
		}
		stemp=(sym_add *)sa->multiply(*sa);
		delete sa;
		sa=stemp;
		
		n/=2;
	}
	
	sym_add *snew=(sym_add *)sa->multiply(*sb);
	delete sa;
	delete sb;
	
	return snew->reduce(flags);
	
}

/////////////////////////////// sym_prod ///////////////////////////

sym::sym_prod::~sym_prod() {

	for(unsigned int i=0;i<oper.size();i++) {
		delete oper[i].first;
	}

}

sym::sym_prod::sym_prod(const sym_prod &s) {

	oper=s.oper;
	for(unsigned int i=0;i<oper.size();i++) 
		oper[i].first=s.oper[i].first->clone();

}

sym::sym_prod *sym::sym_prod::clone() const {
	return new sym_prod(*this);
}

int sym::sym_prod::nodeCount() const {
	int cnt = 1;
	for(unsigned int i=0; i<oper.size(); i++)
		cnt += oper[i].first->nodeCount();
	return cnt;
}

int sym::sym_prod::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_prod *q;
	q=(sym_prod *) &s;
	
	int c;
	c=COMPARE(oper.size(),(q->oper).size());
	if(c!=0) return c;

	for(unsigned int i=0;i<oper.size();i++) {
		c=oper[i].first->comp( *(q->oper[i].first) );
		if(c!=0) return c;
		c=COMPARE(oper[i].second,q->oper[i].second);
		if(c!=0) return c;
	}
	return c;
}

sym::sym_expr *sym::sym_prod::reduce(sym_flags flags) {

	sym_prod s_old(*this);

	for(unsigned int i=0;i<oper.size();i++) oper[i].first=oper[i].first->reduce(flags);

// Find children nodes of type sym_prod and merge them in current node
//  a * [b * c] --> a * b * c

	unsigned int n=oper.size();
	for(unsigned int i=0;i<n;i++) {
		if(typeid(*(oper[i].first))==typeid(sym_prod)) {;
			sym_prod *s;
			rational q;
			s=(sym_prod *)oper[i].first;
			q=oper[i].second;
			oper[i].first=s->oper[0].first;
			oper[i].second=q*s->oper[0].second;
			s->oper[0].first=NULL;
			for(unsigned int j=1; j<s->oper.size();j++) {
				oper.push_back(make_pair(s->oper[j].first,q*s->oper[j].second));
				s->oper[j].first=NULL;
			}
			delete s;
		}
	}

// sym_num nodes should have exponent=1

	n=oper.size();
	for(unsigned int i=0;i<n;i++) {
		if(typeid(*(oper[i].first))==typeid(sym_num)) 
			if(oper[i].second!=1) {
				sym_num *s;
				s=(sym_num *)oper[i].first;
				s->value=std::pow(s->value,oper[i].second.eval());
				oper[i].second=1;
			}	
	}

// sym_exp nodes should have exponent=1

	n=oper.size();
	for(unsigned int i=0;i<n;i++) {
		if(typeid(*(oper[i].first))==typeid(sym_exp)) 
			if(oper[i].second!=1&&oper[i].second!=0) {
				double ex=oper[i].second.eval();
				oper[i].second=1;
				sym_exp *sexp;
				sexp=(sym_exp *)oper[i].first;
				sexp->oper=sym_prod::create(sexp->oper,new sym_num(ex));
			}	
	}
	
// Sort
	sort(oper.begin(),oper.end(),sort_pair_r);

// Multiply equal terms

	n=oper.size();
	for(int	i=n-1;i>0;i--) {
		if(typeid(*(oper[i].first))==typeid(sym_num)) {
			if(typeid(*(oper[i-1].first))==typeid(sym_num)) { // Multiply numbers
				sym_num *s,*s0;
				s=(sym_num *)oper[i].first;
				s0=(sym_num *)oper[i-1].first;
				s0->value*=s->value;
				s->value=1;
			}
		} else {
			if(*(oper[i].first)==*(oper[i-1].first)) { // Multiply other expressions
				oper[i-1].second+=oper[i].second;
				oper[i].second=0;
			}
		}
		
	}
	
// Expand powers of sym_add terms
	if(!flags.collect) {
		for(unsigned int i=0;i<oper.size();i++) {
			if(typeid(*(oper[i].first))!=typeid(sym_add)) continue;
			if(oper[i].second==0) continue;
			if(oper[i].second.den()!=1) {
				// Split integer and fractional part of the exponent
				int ex=(int)floor(oper[i].second.eval());
				if(ex!=0) {
					oper[i].second-=ex;
					oper.push_back(make_pair(oper[i].first->clone(),ex));
				}
			}
			if(oper[i].second.num()!=1) {
				sym_expr *s;
				int ex=abs(oper[i].second.num());
				s=((sym_add *)oper[i].first)->power(ex, flags);
				delete oper[i].first;
				oper[i].first=s;
				oper[i].second=rational(oper[i].second.num()/ex,oper[i].second.den());
			}	
		}
// Sort
	sort(oper.begin(),oper.end(),sort_pair_r);
	
	}
// Remove ones
	for(unsigned int i=0;i<oper.size();i++) {
		if(oper[i].second==0) {
			delete oper[i].first;
			oper.erase(oper.begin()+i--);
		} else if(typeid(*(oper[i].first))==typeid(sym_num)) {
			if( ((sym_num *)oper[i].first)->value==1) {
				delete oper[i].first;
				oper.erase(oper.begin()+i--);
			} 		
		}
	}
	
// Find zeros
	for(unsigned int i=0;i<oper.size();i++) {
		if(typeid(*(oper[i].first))==typeid(sym_num)) {
			if( ((sym_num *)oper[i].first)->value==0) {
				delete this;
				return new sym_num(0);
			}
		}
	}

	
// Check the new size and change to appropriate node type
	if(oper.size()==0) {
		delete this;
		return new sym_num(1);
	}
	if(oper.size()==1) {
		sym_expr *s;
		if(oper[0].second==1) {
			s=oper[0].first;
			oper[0].first=NULL;
			delete this;
			return s->reduce(flags);
		}
	}

// Multiply sym_add terms

	if(!flags.collect) {
		n=oper.size();
		for(unsigned int i=0;i<n;i++) {
			if(typeid(*(oper[i].first))!=typeid(sym_add)) continue;
			if(oper[i].second!=1) continue;
			sym_add *sadd,*snew;
			sadd=(sym_add *)oper[i].first;
			snew=new sym_add();
			for(unsigned int j=0;j<sadd->oper.size();j++) {
				sym_prod *sprod;
				sprod=clone();
				delete sprod->oper[i].first;
				sprod->oper[i].first=sadd->oper[j].first->clone();
				snew->oper.push_back(make_pair(sprod,sadd->oper[j].second));
			}
			delete this;
			return snew->reduce(flags);
		}
	}

// First coefficient of sym_add terms should be one
// otherwise (2*x+2*y)/(x+y) won't be simplified

	for(unsigned int i=0; i<oper.size(); i++) {
		if(typeid(*(oper[i].first))!=typeid(sym_add)) continue;
		sym_add *sadd = (sym_add *)oper[i].first;
		if(sadd->oper[0].second == 1) continue;
		double fact = sadd->oper[0].second;
		sadd->oper[0].second = 1;
		for(unsigned int j=1; j<sadd->oper.size(); j++) {
			sadd->oper[j].second /= fact;
		}
		sym_num *snum = new sym_num();
		snum->value = std::pow(fact, oper[i].second.eval());
		oper.push_back(make_pair(snum, rational(1)));
	}

	if(comp(s_old)) return reduce(flags);
	else return this;
	
}

sym::sym_expr *sym::sym_prod::derive(const sym_expr &s) {

	sym_add *snew;
	snew=new sym_add();

	for(unsigned int i=0;i<oper.size();i++) {
		sym_prod *sprod;
		sprod=clone();
		sprod->oper.push_back(make_pair(oper[i].first->clone()->derive(s),1));
		sprod->oper[i].second-=1;
		snew->oper.push_back(make_pair(sprod,(oper[i].second).eval()));
	}
	
	delete this;

	return snew;
}

rational sym::sym_prod::get_exponent(const sym_expr &factor) {
	// Return exponent of factor
	for (unsigned int i=0; i<oper.size(); i++) {
		if (*oper[i].first == factor) {
			return oper[i].second;
		}
	}

	return 0;
}

matrix sym::sym_prod::eval() const {
	
	matrix m=ones(1,1);
	
	for(unsigned int i=0;i<oper.size();i++) {
		m*=::pow(oper[i].first->eval(),oper[i].second.eval());
	}
	return m;
} 

ostream &sym::sym_prod::print(ostream &os) const {

	for(unsigned int i=0;i<oper.size();i++) {
        if(oper[i].second>=0&&i) os<<"*";
        else if(oper[i].second<0) {
            if(i==0) os<<"1/";
            else os<<"/";
        }
		os<<*(oper[i].first);
		if(abs(oper[i].second)!=1) {
			os<<"^";
			if(oper[i].second.den()!=1) os<<"(";
			os<<abs(oper[i].second);
			if(oper[i].second.den()!=1) os<<")";
		}
	}

	return os;
}

sym::sym_prod *sym::sym_prod::create(sym_expr *s1,sym_expr *s2) {

	sym_prod *snew;
	
	snew=new sym_prod();
	snew->oper.push_back(make_pair(s1,1));
	snew->oper.push_back(make_pair(s2,1));

	return snew;
}

sym::sym_prod *sym::sym_prod::create_pow(sym_expr *s,const rational &q) {

	sym_prod *snew;
	
	snew=new sym_prod();
	snew->oper.push_back(make_pair(s,q));

	return snew;
}

///////////////////////// sym_sin //////////////////////////////

sym::sym_sin::~sym_sin() {
	delete oper;
}

sym::sym_sin::sym_sin(const sym_sin &s) {
	oper=s.oper->clone();
}


sym::sym_sin *sym::sym_sin::clone() const {
	return new sym_sin(*this);
}

int sym::sym_sin::nodeCount() const {
	int cnt = 1;
	cnt += oper->nodeCount();
	return cnt;
}

int sym::sym_sin::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_sin *q;
	q=(sym_sin *) &s;
	
	return oper->comp(*(q->oper));

}

sym::sym_expr *sym::sym_sin::reduce(sym_flags flags) {
	
	oper=oper->reduce(flags);
	if(typeid(*oper)==typeid(sym_num)) {
		double val=((sym_num *)oper)->value;
		delete this;
		return new sym_num(::sin(val));
	}
	return this;

}

sym::sym_expr *sym::sym_sin::derive(const sym_expr &s) {
	
	sym_expr *arg;
	sym_cos *scos;
	scos=sym_cos::create(oper);
	arg=oper->clone();
	oper=NULL;
	delete this;
	
	return sym_prod::create(scos,arg->derive(s));
	
}

matrix sym::sym_sin::eval() const {
	return ::sin(oper->eval());
} 

ostream &sym::sym_sin::print(ostream &os) const {
	return os<<"sin("<<*oper<<")";
}

sym::sym_sin *sym::sym_sin::create(sym_expr *s) {

	sym_sin *snew;
	
	snew=new sym_sin();
	snew->oper=s;
	
	return snew;
}

///////////////////////// sym_cos //////////////////////////////

sym::sym_cos::~sym_cos() {
	delete oper;
}

sym::sym_cos::sym_cos(const sym_cos &s) {
	oper=s.oper->clone();
}


sym::sym_cos *sym::sym_cos::clone() const {
	return new sym_cos(*this);
}

int sym::sym_cos::nodeCount() const {
	int cnt = 1;
	cnt += oper->nodeCount();
	return cnt;
}

int sym::sym_cos::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_cos *q;
	q=(sym_cos *) &s;
	
	return oper->comp(*(q->oper));

}

sym::sym_expr *sym::sym_cos::reduce(sym_flags flags) {
	
	oper=oper->reduce(flags);
	if(typeid(*oper)==typeid(sym_num)) {
		double val=((sym_num *)oper)->value;
		delete this;
		return new sym_num(::cos(val));
	}
	return this;

}

sym::sym_expr *sym::sym_cos::derive(const sym_expr &s) {
	
	sym_expr *arg;
	sym_sin *ssin;
	ssin=sym_sin::create(oper);
	arg=oper->clone();
	oper=NULL;
	delete this;
	
	return sym_prod::create(ssin->mult(sym_num(-1)),arg->derive(s));
	
}

matrix sym::sym_cos::eval() const {
	return ::cos(oper->eval());
} 

ostream &sym::sym_cos::print(ostream &os) const {
	return os<<"cos("<<*oper<<")";
}

sym::sym_cos *sym::sym_cos::create(sym_expr *s) {

	sym_cos *snew;
	
	snew=new sym_cos();
	snew->oper=s;
	
	return snew;
}

///////////////////////// sym_exp //////////////////////////////

sym::sym_exp::~sym_exp() {
	delete oper;
}

sym::sym_exp::sym_exp(const sym_exp &s) {
	oper=s.oper->clone();
}


sym::sym_exp *sym::sym_exp::clone() const {
	return new sym_exp(*this);
}

int sym::sym_exp::nodeCount() const {
	int cnt = 1;
	cnt += oper->nodeCount();
	return cnt;
}

int sym::sym_exp::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_exp *q;
	q=(sym_exp *) &s;
	
	return oper->comp(*(q->oper));

}

sym::sym_expr *sym::sym_exp::reduce(sym_flags flags) {
	
	oper=oper->reduce(flags);
	if(typeid(*oper)==typeid(sym_num)) {
		double val=((sym_num *)oper)->value;
		delete this;
		return new sym_num(::exp(val));
	}
	if(typeid(*oper)==typeid(sym_log)) {
		sym_expr *arg;
		arg=((sym_log *)oper)->oper;
		((sym_log *)oper)->oper=NULL;
		delete this;
		return arg;
	}
	if(typeid(*oper)==typeid(sym_add)) {
		sym_add *sadd;
		sadd=(sym_add *)oper;
		for(unsigned int i=0;i<sadd->oper.size();i++) {
			if(typeid(*(sadd->oper[i].first))==typeid(sym_log)&&sadd->oper[i].second!=0) {
				double dex=sadd->oper[i].second;
				rational qex=0;
				for(int k=1;k<=10;k++) { //Check if the exponent can be repr. by a small fraction
					if(round(k*dex)==k*dex) {
						qex=rational((int)round(k*dex),k);
						break;
					}
				}
				if(qex==0) continue;
				sym_expr *arg;
				arg=((sym_log *)sadd->oper[i].first)->oper->clone();
				sadd->oper[i].second=0;
				return sym_prod::create( sym_prod::create_pow(arg,qex), this)->reduce(flags);
			}
		}
	}
	return this;

}

sym::sym_expr *sym::sym_exp::derive(const sym_expr &s) {
	
	sym_expr *arg;
	arg=oper->clone();
	
	return sym_prod::create(this,arg->derive(s));
	
}

matrix sym::sym_exp::eval() const {
	return ::exp(oper->eval());
} 

ostream &sym::sym_exp::print(ostream &os) const {
	return os<<"exp("<<*oper<<")";
}

sym::sym_exp *sym::sym_exp::create(sym_expr *s) {

	sym_exp *snew;
	
	snew=new sym_exp();
	snew->oper=s;
	
	return snew;
}

///////////////////////// sym_log //////////////////////////////

sym::sym_log::~sym_log() {
	delete oper;
}

sym::sym_log::sym_log(const sym_log &s) {
	oper=s.oper->clone();
}


sym::sym_log *sym::sym_log::clone() const {
	return new sym_log(*this);
}

int sym::sym_log::nodeCount() const {
	int cnt = 1;
	cnt += oper->nodeCount();
	return cnt;
}

int sym::sym_log::comp(const sym_expr &s) const {

	if(typeid(*this)!=typeid(s)) return COMPARE(order(),s.order());
	
	sym_log *q;
	q=(sym_log *) &s;
	
	return oper->comp(*(q->oper));

}

sym::sym_expr *sym::sym_log::reduce(sym_flags flags) {
	
	oper=oper->reduce(flags);
	if(typeid(*oper)==typeid(sym_num)) {
		double val=((sym_num *)oper)->value;
		delete this;
		return new sym_num(::log(val));
	}
	if(typeid(*oper)==typeid(sym_exp)) {
		sym_expr *arg;
		arg=((sym_exp *)oper)->oper;
		((sym_exp *)oper)->oper=NULL;
		delete this;
		return arg;
	}
	if(typeid(*oper)==typeid(sym_prod)) {
		sym_prod *sprod;
		sprod=(sym_prod *)oper;
		if(sprod->oper.size()==1) {
			oper=sprod->oper[0].first;
			sprod->oper[0].first=NULL;
			sym_expr *snew;
			snew=sym_prod::create(new sym_num(sprod->oper[0].second.eval()),this);
			delete sprod;
			return snew->reduce(flags);
		}
		for(unsigned int i=0;i<sprod->oper.size();i++) {
			if(typeid(*(sprod->oper[i].first))==typeid(sym_exp)&&sprod->oper[i].second==1) {				
				sym_expr *arg;
				arg=((sym_exp *)sprod->oper[i].first)->oper->clone();
				sprod->oper[i].second=0;
				return sym_add::create( arg, this)->reduce(flags);
			}
		}
	}
	
	
	return this;

}

sym::sym_expr *sym::sym_log::derive(const sym_expr &s) {
	
	sym_expr *arg;
	arg=oper;
	oper=NULL;
	delete this;
	sym_expr *s1=sym_prod::create_pow(arg->clone(),-1);
	return sym_prod::create( s1,arg->derive(s));
	
}

matrix sym::sym_log::eval() const {
	return ::log(oper->eval());
} 

ostream &sym::sym_log::print(ostream &os) const {
	return os<<"log("<<*oper<<")";
}

sym::sym_log *sym::sym_log::create(sym_expr *s) {

	sym_log *snew;
	
	snew=new sym_log();
	snew->oper=s;
	
	return snew;
}


