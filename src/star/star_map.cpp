#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "star.h"

extern "C" {
#include <stdlib.h>
}

std::vector<double> star2d::init_domain_weight(std::vector<int> &domain_type_new) {
// Control the pressure drop in each domain
	int ndom = domain_type_new.size();
	std::vector<double> domain_weight_new(ndom, 1.);
	int n = 0;
	while(domain_type_new[n] == CORE) n++;
	if (n < ndom) domain_weight_new[n] = 8.;
	n++;
	if (n < ndom) domain_weight_new[n] = 4.;
	n++;
	if (n < ndom) domain_weight_new[n] = 2.;
	//domain_weight_new[ndomains-1] = 8.;
	//domain_weight_new[ndomains-2] = 4.;
	//domain_weight_new[ndomains-3] = 2.;

	return domain_weight_new;
}

void star2d::update_domain_weights() {

}

// take a model and modify the resolution (nth, ndomains, number of

void star2d::remap(int ndom_in,int *npts_in,int nth_in,int nex_in) {
    remapper red(map);  // declaration object of class remapper

    red.set_ndomains(ndom_in);
    red.set_npts(npts_in);
    red.set_nt(nth_in);
    red.set_nex(nex_in);

    if(ndom_in!=ndomains) 
        remap_domains(ndom_in,red); // the new R_i are now known

    map=red.get_map(); // update the mapping
    interp(&red); // interpolate the variable on the new update

}


void star2d::check_map() {

    update_domain_weights();

    matrix Rcc;
    double pcc;

    int convective_core = check_convec(pcc, Rcc);
    if (conv && !convective_core) {  // Convective core disappeared
    	if (remove_convective_core()) return;
    }
    else if(!conv && convective_core) { // Convective core appeared
    	add_convective_core(pcc, Rcc);
    	return;
    }

    remapper red(map);
    if(!remap_domains(ndomains, red)) return;
    if(config.verbose) {printf("Remapping...");fflush(stdout);}
    map=red.get_map();
    interp(&red);
    if(config.verbose) printf("Done\n");

}

bool star2d::remap_domains(int ndom, remapper &red) {

	std::vector<int> index;
	int nzones = count_zones(index);
	matrix pzone;
	matrix Rzone = get_zone_itfs(pzone);
	std::vector<int> index_new;
	std::vector<int> ndom_zone = distribute_domains(ndom, pzone, index_new);
	if (index_new == index) return false;

	std::vector<int> domain_type_new = resample_domain_type(ndom, index_new, index);
	std::vector<double> domain_weight_new = init_domain_weight(domain_type_new);
	matrix Rnew = get_new_boundaries(pzone, Rzone, ndom_zone, domain_weight_new);
	red.set_R(Rnew);
	for(int n=0; n<=nzones; n++)
		red.set_fixed(index_new[n], index[n]);

	domain_type = domain_type_new;
	domain_weight = domain_weight_new;
	conv=0;
	int n=0;
	while(domain_type[n++]==CORE) conv++;

	return true;
}

int star2d::remove_convective_core() {
	for(int n=0; n<conv; n++) domain_type[n] = RADIATIVE;
	conv = 0;
	return 0;
}

void star2d::add_convective_core(double pcc, matrix Rcc) {

	if (ndomains < 2) {
		fprintf(stderr, "Error: At least 2 domains are needed to deal with core convection.\n");
		exit(1);
	}
	if (domain_type[1] != RADIATIVE) {
		fprintf(stderr, "Error: At least 1 radiative domain needed above convective core. Try increasing number of domains\n");
	}
	domain_type[0] = CORE;
	std::vector<int> index;
	int nzones = count_zones(index);
	matrix pzone;
	matrix Rzone = get_zone_itfs(pzone);
	Rzone.setrow(1, Rcc);
	pzone(1) = pcc;
	std::vector<int> index_new;
	std::vector<int> ndom_zone = distribute_domains(ndomains, pzone, index_new);
	std::vector<int> domain_type_new = resample_domain_type(ndomains, index_new, index);
	std::vector<double> domain_weight_new = init_domain_weight(domain_type_new);
	matrix Rnew = get_new_boundaries(pzone, Rzone, ndom_zone, domain_weight_new);

	domain_type = domain_type_new;
	domain_weight = domain_weight_new;
	conv=0;
	int n=0;
	while(domain_type[n++]==CORE) conv++;

	remapper red(map);
	red.set_R(Rnew);
	for(int n=0; n<=nzones; n++)
		if (n != 1) red.set_fixed(index_new[n], index[n]); // Not for new interface CC-RZ
	map = red.get_map();
	interp(&red);

}

int star2d::count_zones(std::vector<int> &index) {
	int nzones = 1;
	index.clear();
	index.push_back(0);
	int zone_type = domain_type[0];
	for(int n=1; n<ndomains; n++) {
		if (domain_type[n] != zone_type) {
			zone_type = domain_type[n];
			nzones++;
			index.push_back(n);
		}
	}
	index.push_back(ndomains);
	return nzones;
}

matrix star2d::get_zone_itfs(matrix &pzone) {
	std::vector<int> index;
	int nzones = count_zones(index);
	matrix Rzone(nzones+1, nth);
	pzone = zeros(nzones+1, 1);
	for(int n=0; n<=nzones; n++) {
		Rzone.setrow(n, map.R.row(index[n]));
		int j0 = 0;
		for(int i=0; i<index[n]; i++) j0 += map.npts[i];
		if (j0 == nr) j0 = nr-1;
		pzone(n) = log(map.leg.eval_00(PRES.row(j0), 0)(0));
	}
	return Rzone;

}

// Distribute domains between zones
std::vector<int> star2d::distribute_domains(int ndom, matrix pzone, std::vector<int> &index_new) {

	int nzones = pzone.nrows() - 1;
	if (nzones > ndom) {
		fprintf(stderr,"Error: At least %d domains are needed for this model\n",nzones);
		exit(1);
	}

	std::vector<int> ndom_zone(nzones, 1);
	for(int n=nzones; n<ndom; n++) {
		double dp = 0;
		int k = 0;
		for(int i=0; i<nzones; i++) {
			double dpzone = (pzone(i) - pzone(i+1)) / ndom_zone[i];
			if ( dpzone  > dp ) {
				k = i;
				dp = dpzone;
			}
		}
		ndom_zone[k]++;
	}

	index_new.resize(nzones + 1);
	index_new[0] = 0;
	for(int n=0; n<nzones; n++) {
		index_new[n+1] = index_new[n] + ndom_zone[n];
	}

	return ndom_zone;

}

std::vector<int> star2d::resample_domain_type(int ndom_new, std::vector<int> index_old, std::vector<int> index_new) {
	std::vector<int> domain_type_new(ndom_new);
	int izone = 0;
	int type = 0;
	for (int n=0; n<ndom_new; n++) {
		if (n == index_new[izone]) {
			type = domain_type[index_old[izone]];
			izone++;
		}
		domain_type_new[n] = type;
	}
	return domain_type_new;
}

matrix star2d::get_new_boundaries(matrix pzone, matrix Rzone, std::vector<int> ndom_zone, std::vector<double> weight) {

	int ndom = 0;
	int nzones = ndom_zone.size();
	for (int n=0; n<nzones; n++) ndom += ndom_zone[n];
	matrix Rnew = zeros(ndom+1, nth);
	int idom = 0;
	for (int n=0; n<nzones; n++) {
		Rnew.setrow(idom, Rzone.row(n));
		double weight_norm = 0;
		for (int i=0; i<ndom_zone[n]; i++) weight_norm += 1./weight[idom+i];
		double pif = pzone(n);
		for (int i=1; i<ndom_zone[n]; i++) {
			idom++;
			pif += (pzone(n+1) - pzone(n)) / weight_norm / weight[idom-1];
			Rnew.setrow(idom, find_boundary(pif));
		}
		idom++;
	}
	Rnew.setrow(ndom, Rzone.row(nzones));
	return Rnew;
}

matrix star2d::find_boundary(double pif) {

    matrix zi(1, nth);
    for(int j=0;j<nth;j++) {
    	double zj, dzj, yj, dyj;
        matrix TT;
        zj = 0;
        int k = 0;
        while (log(PRES(k,j)) > pif) {
        	k++;
        	if (k == nr) break;
        }

        zj = z(k-1);
        int fin = 0, nit = 0;
        while(fin<3) {
            yj = log(map.gl.eval(PRES.col(j),zj,TT))(0);
            dyj = (TT,(D,log(PRES).col(j)))(0);
            dzj = -(yj - pif)/dyj;
            if ( fabs(dzj) < 1e-10 ) fin++;
            dzj = dzj < -zj/2 ? -zj/2 : dzj;
            dzj = dzj > (1-zj)/2 ? (1-zj)/2 : dzj;
            zj += dzj;
            nit++;
            if(nit>100) {
                fprintf(stderr,"Error: (star2d) No convergence in find_boundaries\n");
                exit(1);
            }
        }
        zi.setcol(j, zj*ones(1,1));
    }
    return map.zeta_to_r(zi);
}

// check_convec detects the appearance of a convective core but is not
// used to move the core boundary
int star2d::check_convec(double &p_cc, matrix &Rcc) {
    if(!core_convec) return 0; // core_covec: input param to disable CC

    if(conv) {
        int j=0;
        for(int n=0;n<conv;n++) j+=map.gl.npts[n]; // number of grid in CC
        if(z(j)<min_core_size) {
            if(config.verbose) printf("Size(convective core) < min_core_size. Removing...\n");
            return 0;
        } else return conv;
    }
    // else
    int i=0;
    while(z(i)<1.05*min_core_size) i++;
// "i" is the first grid point where Schwarzschild is tested (to avoid
// too small cores)
    matrix schw;
// Schwarzschild criterion
// schw = -(grad P).(grad s)
// if schw > 0 then stable
// if schw < 0 then unstable = convection zone

    schw=-(map.gzz*(D,p)+map.gzt*(p,Dt))*((D,log(T))-eos.del_ad*(D,log(p)))
        -(map.gzt*(D,p)+map.gtt*(p,Dt))*((log(T),Dt)-eos.del_ad*(log(p),Dt));
    schw.setrow(0,zeros(1,nth));
    schw=schw/r/r;
    schw.setrow(0,zeros(1,nth));
    schw.setrow(0,-(D.row(0),schw)/D(0,0));
    if(map.leg.eval_00(schw.row(i),0)(0)>=0) return 0;
// if Sch > 0 no CC (or CC too small) and return

    if(config.verbose) printf("Found convective core\n");
    if(ndomains==1) {
        fprintf(stderr,"Warning: At least 2 domains are needed to deal with core convection\n");
        return 0;
    }

    while(schw(i,-1)<0) i++; // look for change of sign of schw
    double zi=z(i);
    matrix dschw;
    dschw=(D,schw);
    Rcc=zeros(1,nth);
    matrix pcore(Rcc);
    for(int j=nth-1;j>=0;j--) {
        int fin=0;
        double schwi,dschwi,dzi;
        matrix TT;
        while(fin<3) {
            schwi=map.gl.eval(schw.col(j),zi,TT)(0);
            dschwi=(TT,dschw.col(j))(0);
            dzi=-schwi/dschwi;
            if(fabs(dzi)<1e-9) fin++;
            zi+=dzi;
        }
        Rcc(j)=map.gl.eval(r.col(j),zi,TT)(0); // R(theta) of the CC
        pcore(j)=(TT,PRES.col(j))(0); // PRES on R_CC(theta)
    }

    p_cc = log(map.leg.eval_00(pcore,0)(0)); // p_cc=log(PRES) at theta=0
    return 1;

}
