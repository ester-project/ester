#include"star.h"
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

void star2d::opacity() {

	int error;

	error=opa_calc(Xr,Z,Tc*T,rhoc*rho,opa);
	
	if(error) exit(1);
	
}

void star2d::nuclear() {

	int error;

	error=nuc_calc(Xr,Z,T*Tc,rho*rhoc,nuc);
	
	if(error) exit(1);

}

void star2d::eq_state() {

	int error;
	
	error=eos_calc(Xr,Z,T*Tc,p*pc,rho,eos);	

	rhoc=rho(0);
	rho=rho/rhoc;
	
	if(error) exit(1);

}

void star2d::atmosphere() {

	int error;
	
	if(!strcmp(atm.name,"simple")) {
		ps=surff*gsup()/opa.k.row(-1)/pc;
		matrix q;
		q=surff*p.row(-1)/ps;
		Ts=pow(q,0.25)*Teff()/Tc;
		return;
	}
	
	
	error=atm_calc(Xr,Z,gsup(),Teff(),eos.name,opa.name,atm);
	
	if(error) exit(1);
	
	matrix q;
	ps=surff*atm.ps/pc;
	q=surff*p.row(-1)/ps;//q=ones(1,nth);
	Ts=pow(q,0.25)*atm.Ts/Tc;

}




