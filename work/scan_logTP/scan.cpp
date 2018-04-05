//#include "ester-config.h"
#include "utils.h"
#include "physics.h"
#include "constants.h"

#include <iostream>
#include <cmath>

extern "C" {
    void free_eos_(int *ifoption, int *ifmodified, int *ifion, int *kif_in,
            double *eps, int *neps, double *match_variable, double *tl, double*
            fl, double *t, double *rho, double *rl, double *p, double *pl,
            double *cf, double *cp, double *qf, double *qp, double *sf, double
            *st, double *grada, double *rtp, double *rmue, double *fh2, double
            *fhe2, double *fhe3, double *xmu1, double *xmu3, double *eta, double
            *gamma1, double *gamma2, double *gamma3, double *h2rat, double
            *h2plusrat, double *lambda, double *gamma_e, double *degeneracy,
            double *pressure, double *density, double *energy, double *enthalpy,
            double *entropy, int *iteration_count);
}

//int main(const matrix &X, double Z, const matrix &T, const matrix &p,
        //matrix &rho, eos_struct &eos) {

int main(int argc, char *argv[]) {

	double Z=0.02;
	double X=0.7;
    int ifoption = 1;
    int ifmodified = 2;
    int ifion = 0;
    int kif_in = 1;
    int neps = 20;
    double *eps = new double[neps];
    double tl;
    double fl;
    double rl;
    double pl;
    double cf;
    double cp;
    double qf;
    double qp;
    double sf;
    double st;
    double grada;
    double rtp;

    double rmue;
    double fh2;
    double fhe2;
    double fhe3;
    double xmu1;
    double xmu3;
    double eta;

    double gamma1;
    double gamma2;
    double gamma3;
    double h2rat;
    double h2plusrat;
    double lambda;
    double gamma_e;

    double degeneracy[3];
    double pressure[3];
    double density[3];
    double energy[3];
    double enthalpy[3];
    double entropy[3];

    int iteration_count;

        double_map comp = initial_composition(X, Z);

        eps[0] = comp["H"] / 1.008e0;                       // H
        eps[1] = (comp["He3"] + comp["He4"]) / 4.0026e0;   // He3 + He4
        eps[2] = (comp["C12"] + comp["C13"]) / 12.0111e0; // C12 + C13
        eps[3] = (comp["N14"] + comp["N15"]) / 14.0067e0; // N14 + N15
        eps[4] = (comp["O16"] + comp["O17"]) / 15.9994e0; // O16 + O17
        eps[5] = 0.0; // Ne
        eps[6] = 0.0; // Na
        eps[7] = 0.0; // Mg
        eps[8] = 0.0; // AL
        eps[9] = 0.0; // Si
        eps[10] = .0; // P
        eps[11] = .0; // S
        eps[12] = .0; // Cl
        eps[13] = .0; // A
        eps[14] = .0; // Ca
        eps[15] = .0; // Ti
        eps[16] = .0; // Cr
        eps[17] = .0; // Mn
        eps[18] = .0; // Fe
        eps[19] = .0; // Ni

	double Tmin=5000.;
	double Tmax=1e5;
	double pmin=2500.;
	double pmax=4e7;
	int ntemp=10;
	int npres=10;
	int i,j;
	double match_variable,pi,rhoi;
	double t,wl,grada0,dlnw,Dnab_Dlnw;
	dlnw=1e-3;

	FILE *NA;
	NA=fopen("dnabla","w");
	for (i=0;i<ntemp;i++) {
	  tl=(3.9+i*(5.-3.9)/(ntemp-1))*log(10.);
	  t=pow(10.,tl/log(10.));
          printf("t=%e  ",t);  

	  for (j=0;j<npres;j++) {

// on est rectangulaire dans la variable log(w)=log(p)-4log(T) d'apres MESA
	    wl=(-12.4+j*(-11.4+12.4)/(npres-1))*log(10.);
	    pl=wl+4*tl;
	    pi=pow(10.,pl/log(10.));
            //printf("p=%e pl=%e ",pi,pl);  
            match_variable = log(pi);

            free_eos_(&ifoption, &ifmodified,
                &ifion, &kif_in, eps, &neps, &match_variable, &tl, &fl,
                &t, &rhoi, &rl, &pi, &pl, &cf, &cp, &qf, &qp, &sf, &st, &grada, &rtp,
                &rmue, &fh2, &fhe2, &fhe3, &xmu1, &xmu3, &eta,
                &gamma1, &gamma2, &gamma3, &h2rat, &h2plusrat, &lambda, &gamma_e,
                degeneracy, pressure, density, energy, enthalpy, entropy,
                &iteration_count);

	    grada0=grada;

	    wl=wl+dlnw;
	    pl=wl+4*tl;
            pi=pow(10.,pl/log(10.));
            match_variable = log(pi);
            
            free_eos_(&ifoption, &ifmodified,
                &ifion, &kif_in, eps, &neps, &match_variable, &tl, &fl, 
                &t, &rhoi, &rl, &pi, &pl, &cf, &cp, &qf, &qp, &sf, &st, &grada, &rtp,
                &rmue, &fh2, &fhe2, &fhe3, &xmu1, &xmu3, &eta, 
                &gamma1, &gamma2, &gamma3, &h2rat, &h2plusrat, &lambda, &gamma_e,
                degeneracy, pressure, density, energy, enthalpy, entropy,
                &iteration_count);
	    Dnab_Dlnw=(grada-grada0)/dlnw;

	    fprintf(NA,"%e  %e  %e\n",tl,wl-dlnw/2,Dnab_Dlnw);

			}
		}
    delete[] eps;
    fclose(NA);
    return 0;
}
