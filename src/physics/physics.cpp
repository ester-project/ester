#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "physics.h"
#include <string.h>

#include <ostream>
#include <iostream>


/*
// Singleton instance of abundance manager
AbundanceManager& AbundanceManager::getInstance() {
    static AbundanceManager instance;
    return instance;
}

void AbundanceManager::setAbundances(const std::map<std::string, float>& data) {
    abundances = data;
}

const std::map<std::string, float>& AbundanceManager::getAbundances() const {
    return abundances;
}
*/ 

//int opa_calc(const matrix &X,double Z,const matrix &T,const matrix &rho,
//		opa_struct &opa) {

// Define and initialize the global variable

//AbundanceMap global_abundance_map;  // Default initialization


int opa_calc(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa,const double Xsol, const double Ysol, const double Zsol){//,
//		const AbundanceMap& abundance_map) {
		
	// Use the initialized global abundance map
    	//AbundanceMap& abundance_map = global_abundance_map;

	/*
        for (const auto& pair : abundance_map.comp_abund) {
        	const std::string& key = pair.first;
        	float comp_abund_value = pair.second;
        	
        	std::cout << key << ": " << "comp_abund=" << comp_abund_value << std::endl; 
        }
        
        for (const auto& pair : abundance_map.A_weights) {
        	const std::string& key = pair.first;
        	double A_weight_value = pair.second;
        	
        	std::cout << key << ": " << "A_weights=" << A_weight_value << std::endl;
    	}*/

	int error=0;

	if(!strcmp(opa.name,"opal")) {
		error=opa_opal(X,Z,T,rho,opa);
	} else if(!strcmp(opa.name,"houdek")) {
    	error=opa_houdek(X,Z,T,rho,opa);
	} else if(!strcmp(opa.name,"kramer")) {
        error=opa_kramer(T,rho,opa);
	} else if(!strcmp(opa.name,"cesam")) {
		error=opa_cesam(X, Z, T, rho, opa);    
    } else if (!strcmp(opa.name,"mono")) { 
        //error=opa_opmesa(X, Z, T, rho, opa);
        //error=opa_opmesa(X, Z, T, rho, opa,Xsol,Ysol,Zsol,abundance_map);
        error=opa_opmesa(X, Z, T, rho, opa,Xsol,Ysol,Zsol);
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

	if(!strcmp(nuc.name,"simple_ON") || !strcmp(nuc.name,"simple_CN") || !strcmp(nuc.name,"simple")) {
		error=nuc_simple(X,T,rho,nuc);
	} else if(!strcmp(nuc.name,"cesam")) {
		error=nuc_cesam(X,T,rho,nuc);
    } else {
        ester_err("Unknown nuc. reac. type: %s",nuc.name);
    	return 1;
    }

	return error;

}

//int atm_calc(const matrix &X,double Z,const matrix &g,const matrix &Teff,
//		const char *eos_name,const char *opa_name,atm_struct &atm,const double Xsol, const double Ysol, const double Zsol,
//		const AbundanceMap& abundance_map) {

int atm_calc(const matrix &X,double Z,const matrix &g,const matrix &Teff,
		const char *eos_name,const char *opa_name,atm_struct &atm,const double Xsol, const double Ysol, const double Zsol) {

	int error=0;

	if(!strcmp(atm.name,"onelayer")) {
		//error=atm_onelayer(X,Z,g,Teff,eos_name,opa_name,atm,Xsol,Ysol,Zsol,abundance_map);
		error=atm_onelayer(X,Z,g,Teff,eos_name,opa_name,atm,Xsol,Ysol,Zsol);
    } else {
        ester_err("Unknown atmosphere type: %s", atm.name);
    	return 1;
    }

	return error;


}




