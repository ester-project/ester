#include"star.h"
#include<string.h>
#include<stdio.h>

int star2d::opacity() {

	int error=0;

	if(!strcmp(opa.name,"opal")) {
		error=opa_opal(Xr,Z,Tc*T,rhoc*rho,opa);
	} else if(!strcmp(opa.name,"houdek")) {
    	error=opa_houdek(Xr,Z,T*Tc,rho*rhoc,opa);
	} else if(!strcmp(opa.name,"kramer")) {
    	error=opa_kramer(T*Tc,rho*rhoc,opa);
    } else {
    	printf("Unknown opacity method: %s\n",opa.name);
    	return 1;
    }
	
	return error;
}

int star2d::nuclear() {

	int error=0;
	
	if(!strcmp(nuc.name,"simple")) {
		error=nuc_simple(Xr,Z,T*Tc,rho*rhoc,nuc);
    } else {
    	printf("Unknown nuc. reac. type: %s\n",nuc.name);
    	return 1;
    }
	
	return error;

}

int star2d::eq_state() {

	int error;
	
	if(!strcmp(eos.name,"ideal"))
		error=eos_ideal(Xr,Z,T*Tc,p*pc,rho,eos);
	else if(!strcmp(eos.name,"ideal+rad"))
		error=eos_idealrad(Xr,Z,T*Tc,p*pc,rho,eos);
	else if(!strcmp(eos.name,"opal")) 
		error=eos_opal(Xr,Z,T*Tc,p*pc,rho,eos);
	else {
		printf("Unknown equation of state: %s\n",eos.name);
    	return 1;
    }
	
	rhoc=rho(0);
	rho=rho/rhoc;
	
	return error;

}


