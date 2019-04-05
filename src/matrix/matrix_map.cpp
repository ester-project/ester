#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "matrix.h"

using namespace std;

matrix_map_elem::operator double_map() {
	double_map res;
	matrix_map_elem::iterator it;
	for(it=begin();it!=end();it++) {
		res[it->first]=*(it->second);
	}
	return res;
}

matrix_map_elem &matrix_map_elem::operator=(const double_map &a) {
	matrix_map_elem::iterator it;
	for(it=begin();it!=end();it++) {
		if (a.count(it->first))
			*(it->second)=a[it->first];
	}
	return *this;
}


matrix_map_elem &matrix_map_elem::operator=(double a) {
	matrix_map_elem::iterator it;
	for(it=begin();it!=end();it++) {
		*(it->second)=a;
	}
	return *this;
}

matrix_map_elem &matrix_map_elem::operator=(matrix_map_elem &a) {
	
	double_map b(a);
	
	return *this=b;
}


matrix_map_elem matrix_map::operator()(int nfil, int ncol) {
	
	matrix_map_elem elem;
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) {
		elem[it->first]=&(it->second)(nfil,ncol);
	}

	return elem;
}

const double_map matrix_map::operator()(int nfil, int ncol) const {
	
	double_map elem;
	matrix_map::const_iterator it;
	for(it=begin();it!=end();it++) {
		elem[it->first]=(it->second)(nfil,ncol);
	}

	return elem;
}

const matrix_map matrix_map::row(int irow) const {

	matrix_map res;
	matrix_map::const_iterator it;
	for(it=begin();it!=end();it++) 
		res[it->first]=(it->second).row(irow);
	return res;
}

const matrix_map matrix_map::col(int icol) const {

	matrix_map res;
	matrix_map::const_iterator it;
	for(it=begin();it!=end();it++) 
		res[it->first]=(it->second).col(icol);
	return res;
}

const matrix_map matrix_map::block(int irow1,int irow2,int icol1,int icol2) const {

	matrix_map res;
	matrix_map::const_iterator it;
	for(it=begin();it!=end();it++) 
		res[it->first]=(it->second).block(irow1,irow2,icol1,icol2);
	return res;
}

const matrix_map matrix_map::block_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol) const {

	matrix_map res;
	matrix_map::const_iterator it;
	for(it=begin();it!=end();it++) 
		res[it->first]=(it->second).block_step(irow1,irow2,drow,icol1,icol2,dcol);
	return res;
}

matrix_map &matrix_map::setrow(int irow,const matrix_map &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		if (a.count(it->first))
			(it->second).setrow(irow,a[it->first]);
	return *this;
}

matrix_map &matrix_map::setcol(int icol,const matrix_map &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		if (a.count(it->first))
			(it->second).setcol(icol,a[it->first]);
	return *this;
}

matrix_map &matrix_map::setblock(int irow1,int irow2,int icol1,int icol2,const matrix_map &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		if (a.count(it->first))
			(it->second).setblock(irow1,irow2,icol1,icol2,a[it->first]);
	return *this;
}

matrix_map &matrix_map::setblock_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol,const matrix_map &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		if (a.count(it->first))
			(it->second).setblock_step(irow1,irow2,drow,icol1,icol2,dcol,a[it->first]);
	return *this;
}

matrix_map &matrix_map::setrow(int irow,const matrix &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		(it->second).setrow(irow,a);
	return *this;
}

matrix_map &matrix_map::setcol(int icol,const matrix &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		(it->second).setcol(icol,a);
	return *this;
}

matrix_map &matrix_map::setblock(int irow1,int irow2,int icol1,int icol2,const matrix &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		(it->second).setblock(irow1,irow2,icol1,icol2,a);
	return *this;
}

matrix_map &matrix_map::setblock_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol,const matrix &a) {
	matrix_map::iterator it;
	for(it=begin();it!=end();it++) 
		(it->second).setblock_step(irow1,irow2,drow,icol1,icol2,dcol,a);
	return *this;
}

matrix matrix_map::sum() {

	matrix res(zeros(1,1));
	matrix_map::iterator it;
	for(it=begin();it!=end();it++)
		res+=it->second;
	return res;

}

double double_map::sum() {

	double res=0;
	double_map::iterator it;
	for(it=begin();it!=end();it++)
		res+=it->second;
	return res;

}

matrix_map matrix_map::operator*(double d) const {

	matrix_map::const_iterator it;
	matrix_map A;
	
	for(it=begin();it!=end();it++) 
		A[it->first]=(it->second)*d;
	
	return A;
}

matrix_map matrix_map::operator*(const matrix &m) const {

	matrix_map::const_iterator it;
	matrix_map A;
	
	for(it=begin();it!=end();it++) 
		A[it->first]=(it->second)*m;
	
	return A;
}

matrix_map operator*(const double_map &d,const matrix &m) {

	double_map::const_iterator it;
	matrix_map A;
	
	for(it=d.begin();it!=d.end();it++) 
		A[it->first]=(it->second)*m;
	
	return A;

}

matrix_map &matrix_map::operator*=(const double_map &d) {
	
	double_map::const_iterator it;
	
	for(it=d.begin();it!=d.end();it++) {
		if(count(it->first)) {
			(*this)[it->first]*=(it->second);
		}
	}
	return *this;

}

matrix_map &matrix_map::operator/=(const double_map &d) {
	
	double_map::const_iterator it;
	
	for(it=d.begin();it!=d.end();it++) {
		if(count(it->first)) {
			(*this)[it->first]/=(it->second);
		}
	}
	return *this;

}

matrix_map &matrix_map::operator+=(const matrix_map &m) {
	
	matrix_map::const_iterator it;
	
	for(it=m.begin();it!=m.end();it++) {
		if(count(it->first)) {
			(*this)[it->first]+=(it->second);
		}
	}
	return *this;

}

matrix_map &matrix_map::operator-=(const matrix_map &m) {
	
	matrix_map::const_iterator it;
	
	for(it=m.begin();it!=m.end();it++) {
		if(count(it->first)) {
			(*this)[it->first]-=(it->second);
		}
	}
	return *this;

}

matrix_map &matrix_map::operator*=(const matrix_map &m) {
	
	matrix_map::const_iterator it;
	
	for(it=m.begin();it!=m.end();it++) {
		if(count(it->first)) {
			(*this)[it->first]*=(it->second);
		}
	}
	return *this;

}

matrix_map &matrix_map::operator/=(const matrix_map &m) {
	
	matrix_map::const_iterator it;
	
	for(it=m.begin();it!=m.end();it++) {
		if(count(it->first)) {
			(*this)[it->first]/=(it->second);
		}
	}
	return *this;

}

