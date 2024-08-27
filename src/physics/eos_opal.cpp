#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "physics.h"
#include "constants.h"

#include <iostream>
#include <cmath>

extern"C" {
	void zfs_interp_eos5_(double *z);
	void eos5_xtrin_(double *x,double *ztab,double *t6,double *p,double *r);
	extern struct{
		double esact,eos[10];
	} eeos_;
	extern struct{
		int itime;
	} lreadco_;
}

int eos_opal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos) {
    
    matrix t6,p_mb;
    int i,N,error=0;
    
    static double Z_table=-99;

    t6=T*1e-6;
    p_mb=p*1e-12;

    if(Z!=Z_table) {
    	lreadco_.itime=0;
    	zfs_interp_eos5_(&Z);
	    Z_table=Z;
    }
    
    if(error) {
    	printf("Can't initialize OPAL EOS table\n");
    	return error;
    }

    rho.dim(T.nrows(),T.ncols());
    eos.s.dim(T.nrows(),T.ncols());
    eos.G1.dim(T.nrows(),T.ncols());
    eos.del_ad.dim(T.nrows(),T.ncols());
    eos.G3_1.dim(T.nrows(),T.ncols());
    eos.d.dim(T.nrows(),T.ncols());
    eos.cp.dim(T.nrows(),T.ncols());
    eos.cv.dim(T.nrows(),T.ncols());
    eos.chi_rho.dim(T.nrows(),T.ncols());
    eos.chi_T.dim(T.nrows(),T.ncols());
        
    eos.prad=A_RAD/3*pow(T,4);
    
    N=T.nrows()*T.ncols();

	double Xi,Zi,t6i,p_mbi,rhoi;
    for(i=0;i<N;i++) {
    	Xi=X(i);Zi=Z;t6i=t6(i);p_mbi=p_mb(i);
    	eos5_xtrin_(&Xi,&Zi,&t6i,&p_mbi,&rhoi);
    	rho(i)=rhoi;
    	eos.s(i)=1e6*(*(eeos_.eos+2));
		eos.G1(i)=*(eeos_.eos+7);
		eos.del_ad(i)=1/(*(eeos_.eos+8));
    	eos.G3_1(i)=*(eeos_.eos+7)/(*(eeos_.eos+8));
    	eos.d(i)=(*(eeos_.eos+6))/(*(eeos_.eos+5));
    	eos.cp(i)=1e6*(*(eeos_.eos+7))*(*(eeos_.eos+4))/(*(eeos_.eos+5));
    	eos.cv(i)=1e6*(*(eeos_.eos+4));
    	eos.chi_rho(i)=*(eeos_.eos+5);
    	eos.chi_T(i)=*(eeos_.eos+6);
        if (fabs(rhoi - (-9e99)) < 1e-10) {
            ester_critical(
                    "Values outside OPAL eos table:\n"
                    "  X = %e\n"
                    "  Z = %e\n"
                    "  T = %e\n"
                    "  p = %e", Xi, Zi, t6i, p_mbi);
        }
    }
    if(exist(rho==-9e99)) {
        ester_critical("Values outside OPAL eos table");
        return 1; // Useless because ester_critical already exit
    }

    return 0;
	
}

