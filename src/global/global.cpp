#include "ester-config.h"
#include "matrix.h"

#include <cmath>
#include <map>
#include <string>

extern const double PI=acos(-1);
extern const double SIG_SB=5.670400e-5;
extern const double K_BOL=1.3806503e-16;
extern const double UMA=1.660538922e-24; // Atomic Mass Unit cgs
extern const double RGP=K_BOL/UMA;
extern const double HYDROGEN_MASS=1.67353249e-24;
extern const double A_RAD=7.565767e-15;
extern const double GRAV=6.67384e-8;
extern const double C_LIGHT=2.99792458e10;
extern const double MYR=3.15576e13; // one Myear in seconds

extern const double M_SUN=1.9891e33;
extern const double R_SUN=6.95508e10;
extern const double L_SUN=3.8396e33;

#include "matrix.h"

extern "C"
const char *ester_version() {
    return VERSION;
}

double_map AMASS = create_double_map
		("n",1.008665)
		("p",1.00727647)
		("H",1.007825)
		("H2",2.0141018)
		("He3",3.0160293)
		("He4",4.0026033)
		("Li6",6.015121)
		("Li7",7.0160040)
		("Be7",7.0169292)
		("Be9",9.0121821)
		("B11",11.0093055)
		("C12",12.)
		("C13",13.0033548)
		("N13",13.0057386)
		("N14",14.003074)
		("N15",15.001089)
		("O16",15.9949146)
		("O17",16.9991315)
		("O18",17.9991604)
		("F18",18.0009377)
		("F19",18.9984032)
		("Ne20",19.9924402)
		("Ne21",20.9938467)
		("Ne22",21.9913855)
		("Na23",22.9897697)
		("Mg23",22.9941249)
		("Mg24",23.9850419)
		("Mg25",24.985837)
		("Mg26",25.982593)
		("Si28",28.0855)
		("Al27",26.9854)
		("S32",31.972070)
		("P31",30.973762)
		("Fe56",55.847);
