#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include <set>
#include <string>
#include "utils.h"
#include "physics.h"
#include "parser.h"

bool _init_metals = true;

bool first = true;

double_map initial_composition(double X, double Z) {
	if(!first){
		//ester_err("initial_composition must be called only once");
	}
	first = false;

	double_map comp;
	file_parser fp;
	
	double Y = 1. - (X + Z);

	comp["H"] = X;
	// TODO: should this one be hardcoded too?
	comp["He3"] = 3.15247417638132e-04 * Y;
	comp["He4"] = Y - comp["He3"];
	// Note that the previous way of setting He3 He4 ensure that He3+He4 = Y

	// TODO: change this hardcoded
	char file[] = "metal-mix.cfg";

	char* arg = NULL;
	char* val = NULL;
	std::set<std::string> metals = {"C12","C13","N14","N15","O16","O17","Ne20","Ne22","Na23",
									"Mg24","Mg25","Mg26","Al27","Si28","P31","S32","Cl35","Cl37",
									"A40","Ca40","Ti","Cr","Mn55","Fe","Ni"};
	// Initialization of comp
	for(std::string metal: metals){
		comp[metal] = .0;
	}
	double metals_fraction = .0;

	if(!fp.open(file)){
		printf("Can't open configuration file %s\n", file);
		perror("Error:");
		exit(1);
	} else {
		int line;
		while(line = fp.get(arg,val)) {
			if(val == NULL){
				printf("Syntax error in configuration file %s, line %d\n", file, line);
				exit(1);
			}
			if(metals.find(arg) == metals.end()){
				printf("%s is unknown, possible metals composition are: ", val);
				for(std::string metal: metals){
					printf("%s ", metal);
				}
				puts("\n");
				exit(1);
			}
			metals_fraction += atof(val);
			comp[arg] = Z * atof(val);
		}
	}
	fp.close();
	comp["Ex"] = 1 - comp.sum();

	if(_init_metals){
		printf("A config file %s has been used to config metal composition\n", file);
		printf("The sum of metals fractions (of Z sum of metal mass fraction) is %f, and should be as close to 1 as possible.\n", metals_fraction);
		_init_metals = false;
	}

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

