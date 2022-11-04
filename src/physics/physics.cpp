#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "physics.h"
#include <string.h>

int opa_calc(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa) {

	int error=0;

	if(!strcmp(opa.name,"opal")) {
		error=opa_opal(X,Z,T,rho,opa);
	} else if(!strcmp(opa.name,"houdek")) {
    	error=opa_houdek(X,Z,T,rho,opa);
	} else if(!strcmp(opa.name,"kramer")) {
        error=opa_kramer(T,rho,opa);
	} else if(!strcmp(opa.name,"cesam")) {
		error=opa_cesam(X, Z, T, rho, opa);
    } else if(!strcmp(opa.name,"opmono")) {
		error=opa_opmono(X, Z, T, rho, opa);  
     } else if(!strcmp(opa.name,"opmesa")) {
		error=opa_opmesa(X, Z, T, rho, opa);          
    } else {
        ester_err("Unknown opacity method: %s",opa.name);
    	return 1;
    }

	return error;
}

int eos_calc(const matrix &X,double Z,const matrix &T,const matrix &p,
        matrix &rho,eos_struct &eos) {

    int error=0;

    if(!strcmp(eos.name,"ideal"))
        error=eos_ideal(X,Z,T,p,rho,eos);
    else if(!strcmp(eos.name,"ideal+rad"))
        error=eos_idealrad(X,Z,T,p,rho,eos);
    else if(!strcmp(eos.name,"opal"))
        error=eos_opal(X,Z,T,p,rho,eos);
    else if(!strcmp(eos.name,"freeeos"))
        error = eos_freeeos(X, Z, T, p, rho, eos);
    else {
        ester_err("Unknown equation of state: %s",eos.name);
        return 1;
    }

    return error;

}

int nuc_calc(const matrix_map &X,const matrix &T,const matrix &rho,
		nuc_struct &nuc) {

	int error=0;

	if(!strcmp(nuc.name,"simple")) {
		error=nuc_simple(X,T,rho,nuc);
	} else if(!strcmp(nuc.name,"cesam")) {
		error=nuc_cesam(X,T,rho,nuc);
    } else {
        ester_err("Unknown nuc. reac. type: %s",nuc.name);
    	return 1;
    }

	return error;

}

int atm_calc(const matrix &X,double Z,const matrix &g,const matrix &Teff,
		const char *eos_name,const char *opa_name,atm_struct &atm) {

	int error=0;

	if(!strcmp(atm.name,"onelayer")) {
		error=atm_onelayer(X,Z,g,Teff,eos_name,opa_name,atm);
    } else {
        ester_err("Unknown atmosphere type: %s", atm.name);
    	return 1;
    }

	return error;


}




