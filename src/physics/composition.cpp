#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "physics.h"

double_map initial_composition(double X, double Z) {
	double_map comp;
	
	comp["H"]=X;
	comp["He3"]=3.15247417638132e-04*(1.-X-Z);
	comp["He4"]=(1.-X-Z)-comp["He3"];
	comp["C12"]=Z*1.71243418737847e-01;
	comp["C13"]=Z*2.06388003380057e-03;
	comp["N14"]=Z*5.29501630871695e-02;
	comp["N15"]=Z*2.08372414940812e-04;
	comp["O16"]=Z*4.82006487350336e-01;
	comp["O17"]=Z*1.95126448826986e-04;
	
	double tot=comp.sum();
	
	comp["Ex"] =1-tot;
	
	return comp;

}

double_map initial_composition_cno_cycle_core(double X, double Z) {
	double_map comp;
	// The function sets the initial composition to account for the almost instanteous conversion of C12 to N14 through the CNO cycle.
	// Afterwards, the abundance of N14 is dominated by the reaction of O16 to F17, which decays to O17.
	// O17 is then coverted to N14 via O17 + p --> N14 + He4.  
	comp["H"]=X;
	comp["He3"]=3.15247417638132e-04*(1.-X-Z);
	comp["He4"]=(1.-X-Z)-comp["He3"];
	comp["C12"]=0; // Z*1.71243418737847e-01;
	comp["C13"]=Z*2.06388003380057e-03;
	comp["N14"]=Z*5.29501630871695e-02 + 14.003074*Z*(1.71243418737847e-01/12.011 + 1.95126448826986e-04/16.9991315); //Assumes Z = Zini throughout the evolution!
	comp["N15"]=Z*2.08372414940812e-04;
	comp["O16"]=Z*4.82006487350336e-01;
	comp["O17"]=0;  // Z*1.95126448826986e-04;
	
	double tot=comp.sum();
	
	comp["Ex"] =1-tot;
	
	return comp;

}

matrix composition_map::X() const {

	return (*this)["H"];

}

matrix composition_map::Y() const {

	return (*this)["He3"]+(*this)["He4"];

}

matrix composition_map::Z() const {

	return 1.-X()-Y();

}

