#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include<cmath>
#include"physics.h"

int nuc_simple(const composition_map &comp,const matrix &T,const matrix &rho,nuc_struct &nuc) {

	matrix T9;
	double C,N,XCNO;
	double Z=comp.Z()(-1); //(-1) as in the last value? or negative of it...
	matrix X(comp.X());
	
	T9=T/1e9;
	C=0.173285;
	N=0.053152;
	// O=0.482273;
	XCNO=Z*(C+N);

// nuc.pp from the CESAM code (cited in Espinosa Lara Rieutord 2007)
// Actually used in routine pp1.f of CESAM code citing GONG

	nuc.pp=8.2398e4*X*X*pow(T9,-2./3.)*exp(-3.6*pow(T9,-1./3.));
// nuc.cno from Kippenhahn & Weigert (1991)
	nuc.cno=8.67e25*XCNO*X*pow(T9,-2./3.)*exp(-15.228*pow(T9,-1./3.))*(1+0.027*pow(T9,1./3.)-0.778*pow(T9,2./3.)-0.149*T9);
		
// nuc.talpha combining Angulo+1999, revision of l2 term by Tsumura+2021 and Lang+1998
        
        // eps_talpha * rho 
        
        double f_talpha;
        matrix X4 = comp["He4"]; // the abundance of Helium 4 
        
        f_talpha = 1;
        nuc.talpha1 = 5.09e8*f_talpha*pow(rho *  X4,3)*pow(T9,-3)*exp(-4.4027/T9); // erg/s/cm^3 --> formula was simplified, f_talpha is 1? 
        
        // eps_talpha 
        
        double Q_talpha_erg = 0.0000116558398575; // 7.275 MeV in erg 
        double Na =  6.02214076e23; // Avagadro's number 
        double m_alpha = 4.0330 * 1/Na ; // 1 AMU = 1/Na, mass of helium nucleus is 2 protons (1p = 1.0078) and 2 neutrons (1n = 1.0087)
        
        nuc.talpha2 = Q_talpha_erg/6 * pow(rho*X4/m_alpha,3) *
         (2.05e-8*pow(T9,-3)*exp(-4.405/T9)+5.315e-7*pow(T9,-3)*exp(-27.4/T9))/pow(Na,2); // eps_talpha = Q*r_talpha = Q * (rho*X4/m_alpha)**3 / 6 * lamba_talpha
        nuc.talpha2*=rho; // here we're assuming in the notes eps is volumetric energy, multiplying by density changes to mass energy
        
        // third formula for error checking, we're getting different values
        
        /*
        matrix l1 = (2.43e9/pow(T9,2/3) * exp(-13.49/pow(T9,1/3)-pow(T9/0.15,2))*(1+74.5*T9)+
        		6.09e5/pow(T9,3/2)*exp(-1.054/T9));
        matrix l2 = 3.055e-10/pow(T9,2/3) * exp(-23.135/pow(T9,1/3)-pow(T9/0.4,2))*(1+187.12*T9 + 4.294e3*pow(T9,2))+
        		4.909e-14/pow(T9,3/2)*exp(-3.35/T9)+
        		9.551e-12/pow(T9,3/2)*exp(-26.84/T9);
        nuc.talpha3 = (l1 * l2)/pow(Na,2) /6 * pow(rho*X4/m_alpha,3) * Q_talpha_erg;
        nuc.talpha3*=rho;
        */
        
        matrix l1=(2.43e9/pow(T9,2./3)*exp(-13.49/pow(T9,1./3)-pow(T9/0.15,2))*(1+74.5*T9)+6.09e5/pow(T9,1.5)*exp(-1.054/T9));
        matrix l2=3.055e-10/pow(T9,2./3)*exp(-23.135/pow(T9,1./3)-pow(T9/0.4,2))*(1+187.12*T9+4.294e3*pow(T9,2))+4.909e-14/pow(T9,1.5)*exp(-3.35/T9)
        +9.551e-12/pow(T9,1.5)*exp(-26.84/T9);   
        
        nuc.talpha3 = l1*l2/6*pow(rho*X4/m_alpha,3)/pow(Na,2) * Q_talpha_erg;
        nuc.talpha3*=rho;
        
	nuc.pp*=rho;
	nuc.cno*=rho;
	
	//printf("\n nuc.pp: %e \n", nuc.pp(0));
	//printf("nuc.cno: %e \n" , nuc.cno(0));
	//printf("nuc.talpha1: %e \n" , nuc.talpha1(0));
	//printf("nuc.T9: %e\n" , T9(0));
	//printf("nuc.l1*l2: %e, \n" , l1(0)*l2(0));
	//printf("nuc.talpha2: %e \n" , nuc.talpha2(0));
	//printf("nuc.talpha3: %e \n" , nuc.talpha3(0));
		
	nuc.eps=nuc.pp+nuc.cno;
	nuc.dlneps_lnrho=ones(T.nrows(),T.ncols());
	nuc.dlneps_lnT=nuc.pp/nuc.eps*(-2./3.+1.2*pow(T9,-1./3.))+
		       nuc.cno/nuc.eps*(-2./3.+15.228/3.*pow(T9,-1./3.)+
		       (0.027/3.*pow(T9,1./3.)-2.*0.778/3.*pow(T9,2./3.)-0.149*T9)/(1+0.027*pow(T9,1./3.)-0.778*pow(T9,2./3.)-0.149*T9));

	for(int i=0; i<nuc.eps.nrows(); i++)
		for(int j=0; j<nuc.eps.ncols(); j++)
			if (X(i,j)==0) nuc.dlneps_lnT(i,j) = 0;

	return 0;
}
