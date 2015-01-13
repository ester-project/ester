#include "ester-config.h"
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

matrix composition_map::X() const {

	return (*this)["H"];

}

matrix composition_map::Y() const {

	return (*this)["He3"]+(*this)["He4"];

}

matrix composition_map::Z() const {

	return 1.-X()-Y();

}
