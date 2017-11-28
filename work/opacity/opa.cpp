/* Calculates the opacity for a given X,Z, rho & T using OPAL tables
 for T >= pow(10,3.75) and Houdek for T < pow(10,3.75). 
The mixture is that of GN93

Z 	The metallicity
X 	The hydrogen mass fraction
t6 	The temperature in millions of degrees Kelvin 

Usage : opa  X Z rho T 
Example: opa 0.7064571423 0.02 1.38e-7 2800.
	 gives:
	 Opacity= 5.368516e-02
*/

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>


extern"C" {
	void opacgn93_(double *Z,double *X,double *t6,double *r);
	extern struct{
		double opact,dopact,dopacr,dopactd;
	} e_;
}

double opa_opal_i(double X,double Z,double T,double rho,double &dlogkt,double &dlogkr) {
    
    double t6,r,Xs,Zs;
    
    t6=T*1e-6;
    r=rho/t6/t6/t6;
	Xs= X;
	Zs= Z;
    
    opacgn93_(&Zs,&Xs,&t6,&r);
    dlogkt=e_.dopact;
    dlogkr=e_.dopacr;
    
	return (double) e_.opact;
}

extern"C" {
	void maceps_(double *eps);
	void opinit_(double *eps,int *iorder,char *tabnam,int *imode,int tabnamlen);
	void opints_(double *x,double *z,double *tlg,double *rlg,double *opalg,
			double *opr,double *opt,double *opx,double *opz,int *iexp,int *ier);
}

double opa_houdek_i(double X,double Z,double T,double rho,double &dlogkt,double &dlogkr) {
    
    double logT,t6,logR,logk,dlogkx,dlogkz,Xs,Zs;
    int iexp,ier;
        
    logT=log10(T);
    t6=T*1e-6;
    logR=log10(rho/t6/t6/t6);
    Xs=(float) X;
    Zs=(float) Z;	
    
    opints_(&Xs,&Zs,&logT,&logR,&logk,&dlogkr,&dlogkt,&dlogkx,&dlogkz,&iexp,&ier);
    
	return (double) logk;
}
int main(int argc, char *argv[]) {

//	std::ofstream fhrd ("output_temp.txt", std::ios::in | std::ios::out | std::ios::ate);
	double opa,dlogkt,dlogkr;
	double rho,X,Z,T;

	
	
	if (argc != 5){		
	printf("erreur\n");
	}

	X=atof(argv[1]); 
	Z=atof(argv[2]);
	rho=atof(argv[3]);
	T=atof(argv[4]);

	int i;
	int N=80;
	char tabnam[80];
	int iorder, imode;
	double eps;
	static int init=0;

	if (T < pow(10,3.75)) {
   		printf("Using Houdek\n");
		
		if(!init) {

			iorder=4;
			imode=2;
			maceps_(&eps);
			sprintf(tabnam,"/home/rieutord/Ester/local/share/ester/tables/houdek/v9/OPINTPATH_AX"); 
			for(i=strlen(tabnam);i<N;i++) tabnam[i]=' ';
			opinit_(&eps,&iorder,tabnam,&imode,N);
			init=1;
		}		
	

		opa=opa_houdek_i(X,Z,T,rho,dlogkt,dlogkr);
		opa=pow(10,opa);

	}
	else {
   		printf("Using OPAL\n");
		opa=opa_opal_i(X,Z,T,rho,dlogkt,dlogkr);
		opa=pow(10,opa);
	}

	printf("Opacity= %e\n", opa);

//        fhrd  << opa << std::endl;
	return 0;

}
