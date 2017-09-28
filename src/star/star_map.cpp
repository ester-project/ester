#include "ester-config.h"
#include "star.h"

extern "C" {
#include <stdlib.h>
}

void star2d::new_check_map() {
        DEBUG_FUNCNAME;
        matrix pif,R_inter,p_inter;
        remapper *red;
        matrix R(ndomains+1,nth);
	R=zeros(ndomains+1,nth);
//printf("new_check_map config.input_file=%d\n",config.input_file);
//for (int n=0;n<ndomains;n++) printf("domain_type[%d] = %d\n",n,domain_type[n]);
	if (glit == 1 && config.input_file==0) {
// At first iteration set the number of zone_type to 1
                zone_type = std::vector<int>(1);
	}

	if (global_err < 0.01) { // look for new convective regions and eventually remap the star
	   // Find the zone boundaries and the associated pressures
	   // and output zone_type as global var.
	   int nzones_av=zone_type.size();
	   find_zones(R_inter, p_inter);
	   int nzones_ap=zone_type.size();
	   if (nzones_av == nzones_ap) {
		printf("Number of zones unchanged: do nothing\n");
		return;
	   }
	   // Redistribute the domains and output izif (index of zone interface)
	   // as a global var.
	   pif=New_distribute_domains(ndomains,p_inter);
	   // fill the matrix R with the radii of the first interface
	   // to the lat interface, so 1==> -2, for all theta 0==>-1
	   R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
	   red=new remapper(map);
	   red->set_R(R);
	// Install the new mapping and do interpolation for this mapping
	   map=red->get_map(); // generate r,r_z,...
	   interp(red);
FILE *fic=fopen("new_R.txt", "a");
for (int k=0;k<=ndomains;k++) fprintf(fic,"k= %d R= %e \n",k,R(k,0));
fclose(fic);
	   delete red;
	} else return; // else do nothing
}

matrix star2d::New_distribute_domains(int ndom,matrix p_inter) {
// ndom ==> input
// called by check_map to redistribute domain

    DEBUG_FUNCNAME;
    double p_s;
    matrix pif(ndom,1); // pressure at domain interfaces
    int nzones=zone_type.size(); // zone_type is initialized by find_zones
FILE *fic=fopen("New_pif.txt", "a");
printf("Start NEW distribute domains \n");
printf("nzones =  %d \n",nzones);
printf("p_inter nrows %d \n",p_inter.nrows());
printf("p_inter ncols %d \n",p_inter.ncols());
printf("p_inter(nzones-1) %e \n",p_inter(nzones-1));

    p_s=map.leg.eval_00(PRES.row(-1),0)(0);
fprintf(fic,"p_s= %e\n",p_s);

	if (nzones == 1) { // A single domain (purely radiative)
           for(int n=0;n<ndom;n++) pif(n)=exp((n+1)*log(p_s)/ndom);
        } else { // convection is now taken into account

// We attribute one domain to each zone and then distribute other domains
// so as to minimize the "PRES" jump
// 1/ Check there is enough domains
	if ( ndom < nzones ) printf("Stop not enough domains!!");
	if ( ndom == nzones ) {
         printf("Likely not enough domains, but I try");
        } else { // This is the general case
		matrix dlog(nzones,1);
		int ndom_left=ndom-nzones;
		int iz,izmax,k,ki,ndom_check;
		double dlogmax;
		std::vector <int> ndz,lifi;
		for (int k=0;k<nzones;k++) {ndz.push_back(1);izif[k]=1;}
		for (iz=0; iz<nzones; iz++) {
                 if (iz==0) { dlog(iz)=fabs(log(p_inter(0)));
			} else {
			    dlog(iz)=fabs(log(p_inter(iz)/p_inter(iz-1))) ;
			}
		}
//for (int k=0;k<nzones;k++) fprintf(fic,"%d %e \n",k,dlog(k));
// insert an additional domain where the PRES drop is max.
// until no domain is left!
		while (ndom_left >=1) {
		  ndom_left=ndom_left-1;
                  dlogmax=max(dlog);
		  izmax=99;
                  for (iz=0; iz<nzones; iz++) {
                    if (fabs(dlog(iz)-dlogmax)<1e-15) izmax=iz;
                  }
		  if (izmax==99) printf("izmax = %d, ndom_left= %d\n",izmax,ndom_left);
                  dlog(izmax)=dlog(izmax)*ndz[izmax]/(ndz[izmax]+1);
		  ndz[izmax]=ndz[izmax]+1;
        	}             
for (int k=0;k<nzones;k++) fprintf(fic,"zone %d ==> %d domains \n",k,ndz[k]);
// check total number of domains
		ndom_check=0;
		for (iz=0; iz<nzones;iz++) ndom_check=ndom_check+ndz[iz];
		printf("ndom_check = %d\n",ndom_check);
		if (ndom_check != ndom) {
    			std::cerr << "ERROR: ndoms do not match in new_distrib...\n";
    			std::terminate();
		}
// Now each zone has the optimal number of domains.

// Let's compute the pif with a law of the form exp(a*k+b)
           // First deal with the first zone
		for (k=0;k<ndz[0];k++) pif(k)=exp((k+1)*log(p_inter(0))/ndz[0]);
		ki=ndz[0]-1;
                izif[0]=ki;
		printf("izif (%d) = %d\n",0,ki);
           // First zone done
		for (iz=1; iz<nzones; iz++) {
			for (k=ki;k<ki+ndz[iz];k++) {
				double a=log(p_inter(iz)/p_inter(iz-1))/ndz[iz];
				pif(k)=exp(a*(k-ki)+log(p_inter(iz-1)));
			}
 			ki=ki+ndz[iz];
			izif[iz]=ki;
			printf("izif (%d) = %d\n",iz,ki);
		}
           // All zones done
// last interface is the surface and not account by the foregoing algo
		pif(ndom-1)=p_inter(nzones-1);
// Now we can tag each domain with its type
		// First zone is special
		for(int n=0;n<=izif[0];n++) domain_type[n] = zone_type[0];
		// Other zones
		for (iz=1; iz<nzones; iz++) {
		  for(int n=izif[iz-1]+1;n<=izif[iz];n++) domain_type[n]=zone_type[iz];
		}
	for (int n=0;n<ndomains;n++) printf("domain_type[%d] = %d\n",n,domain_type[n]);


// Finally we need the index of the interfaces between the zones once the
// domains have been distributed. These are stored in izif.

for (int k=0;k<nzones;k++) fprintf(fic,"k= %d p_inter %e \n",k,p_inter(k));
for (int k=0;k<ndomains;k++) fprintf(fic,"k= %d pif %e \n",k,pif(k));
fclose(fic);
	}
}

    return pif;
}


// take a model and modify the resolution (nth, ndomains, number of
// points in each domain) 
// object == remapper defined here
// used to move from one grid to another
// be careful of discontinuities that might fall in the middle of a
// domain.

void star2d::remap(int ndom_in,int *npts_in,int nth_in,int nex_in) {
    remapper red(map); // declaration object of class remapper 

	printf("    Enter remap in star_map\n");
    red.set_ndomains(ndom_in);
    red.set_npts(npts_in);
    red.set_nt(nth_in);
    red.set_nex(nex_in);

    if(ndom_in!=ndomains) 
        remap_domains(ndom_in,red); // the new R_i are now known

    map=red.get_map(); // update the mapping
    interp(&red); // interpolate the variable on the new update

	printf("    Leave remap in star_map\n");
}

// Some domains have boundaries imposed by the physics and these
// boundaries cannot be moved (ex. CC & RZ interface).
// domain_type = integer for CC RZ CZ (see star.h)
bool star2d::remap_domains(int ndom, remapper &red) {
    //	Count zones
    int nzones=1;
    std::vector<int> index;	
// Here look for interface between zones of different type.
if (config.verbose == 19) printf("++++++ start of remap_domains\n");
    for(int n=1,type=domain_type[0];n<ndomains;n++) {
        if(domain_type[n]!=type) {
            index.push_back(n-1);
            nzones++;
            type=domain_type[n];
        }
    }
// index == index of interface sepa. zone of diff type from center to
// surface
    index.push_back(ndomains-1);

// zif = zeta of interfaces of all zones before remapping 
    matrix zif(nzones,nth);
    for(int n=0,j0=-1,izone=0;n<ndomains;n++) {
        j0+=map.gl.npts[n];
        if(n==index[izone]) {
            zif.setrow(izone,z(j0)*ones(1,nth));
            izone++;
        }
    }

    std::vector<int> index_new;

// test if interpolation necessary
    index_new=distribute_domains(ndom,zif,true);

// if nothing has changed return
    if(index_new==index&&ndom==ndomains) {
     if (config.verbose == 19) printf("quit remap_domain with nothing changed\n"); return false;
    }

// compute new indices of interfaces between zones and recompute the zeta of the new
// domains (zif)
printf("++++++ remap_domains and redistribute domains\n");
    index_new=distribute_domains(ndom,zif);

// Now that the zif are known we can recompute the radii R of the new domains inside
// the zones
    red.set_R(zeros(1,nth).concatenate(map.zeta_to_r(zif)));
    for(int n=0;n<nzones;n++) 
        red.set_fixed(index_new[n]+1,index[n]+1); // do not update the interface between zones
// Update domain_type (each domain is tagged with the zone type)
    std::vector<int> domain_type_new(ndom,0);
    for(int n=0,izone=0;n<ndom;n++) {
        domain_type_new[n]=domain_type[index[izone]];
        if(n==index_new[izone]) izone++;
    }

    domain_type=domain_type_new;
    conv=0;
    int n=0;
    while(domain_type[n++]==CORE) conv++; // update of the conv variable

printf("++++++ remap_domains end with redistributed domains\n");
    return true;

}

// if check_only 'True', return the indices of the zones interfaces
// usually used with check_only "false" to 
std::vector<int> star2d::distribute_domains(int ndom,matrix &zif,bool check_only) const {
    matrix dlogT;
    int nzones=zif.nrows();

//printf("Start distribute domains with check_only, nzones = %d check = %d\n",nzones,check_only);
    if(nzones>ndom) {
        fprintf(stderr,"Error: At least %d domains are needed for this model\n",nzones);
        exit(1);
    }

    // Calculate Delta(log(T)) in each zone at theta=0
    dlogT.zero(nzones,1);
    for(int n=0;n<nzones;n++) {
        dlogT(n)=-log(map.gl.eval(PRES.col(-1),zif(n,-1))(0));
        if(n) dlogT(n)+=log(map.gl.eval(PRES.col(-1),zif(n-1,-1))(0));
    }

// Distribute the domains (a limited number) into the different zones so as to
// decrease optimally the LogT jump between two following domains interfaces.

    int ndomi[nzones]; // Number of domains in each zone
    for(int n=0;n<nzones;n++) ndomi[n]=1;
    // Distribute domains based on dlogT
    for(int n=nzones;n<ndom;n++) {
        double dTi=0;
        int k=0;
        for(int i=0;i<nzones;i++) {
            if(dlogT(i)/ndomi[i]>dTi) {
                k=i;
                dTi=dlogT(i)/ndomi[i];
            }
        }
        ndomi[k]++;
    }

// Calculate the new indices of zone interfaces
    std::vector<int> index;
    index.resize(nzones);

    for(int n=0,k=1,izone=0;n<ndom-1;n++,k++) {
        if(k==ndomi[izone]) {
            index[izone]=n;
            k=0;
            izone++;
        }
    }
    index[nzones-1]=ndom-1;
    if(check_only) return index;
// END OF PART check_only = .true.

// Find new boundaries (ie the new zetai= zif_new) on the old zeta-grid
    matrix zif_new(ndom,nth);
    matrix logTi(ndom-1,nth);
    matrix logT0,logT1(1,nth);
    logT0=zeros(1,nth);
    for(int j=0;j<nth;j++) logT1(j)=log(map.gl.eval(PRES.col(j),zif(0,j))(0));
    for(int n=0,k=1,izone=0;n<ndom-1;n++,k++) {
        if(k==ndomi[izone]) {
            k=0;
            izone++;
            for(int j=0;j<nth;j++) { // evaluate PRES on the interfaces bounding the domain
                logT0(j)=log(map.gl.eval(PRES.col(j),zif(izone-1,j))(0));
                logT1(j)=log(map.gl.eval(PRES.col(j),zif(izone,j))(0));
            }
            logTi.setrow(n,logT0);
        } else {
// For the interfaces of domains inside the same zone attribute a value of the
// log(temperature) as a linear function of the domain rank
            logTi.setrow(n,logT0+(logT1-logT0)*((double) k)/ndomi[izone]);		
        }
    }

// A priori the "temperature" (in fact the variable PRES)  on the interfaces of zones is not constant; hence the
// interfaces of domains inside a zone have not a constant temperature. But we like to
// have isothermal interfaces of domain inside a zone for (eg) numerical stability so
// the option ifdef....

#ifdef T_CONSTANT_DOMAINS
    logTi=map.leg.eval_00(logTi,0)*ones(1,nth);
#endif
    zif_new.setblock(0,ndom-2,0,-1,find_boundaries(logTi));

    for(int izone=0;izone<nzones;izone++)
        zif_new.setrow(index[izone],zif.row(izone));

    zif=zif_new;
    return index;

}

// give a function logTi(theta) and find the associated zeta_i(theta_k) of the surface
// PRES(zeta,theta_k)=logTi(theta_k)
// the zeta_k(theta_k) are estimated with a Newton solver
matrix star2d::find_boundaries(const matrix &logTi) const {
    matrix zi(logTi.nrows(),nth);

    for(int j=0;j<nth;j++) {
        matrix zj,dzj,yj,dyj,TT;
        zj=zeros(logTi.nrows(),1);
        int l=0;
        for(int k=0;k<nr;k++) {
            if(PRES(k,j)<logTi(l,j)) {
                zj(l)=z(k-1);
                l++;
                if(l==logTi.nrows()) break;
            }
        }
        if(l<logTi.nrows())
            for(int k=l;k<logTi.nrows();k++) zj(k)=z(-1);
        int fin=0,nit=0;
        while(fin<3) {
            yj=log(map.gl.eval(PRES.col(j),zj,TT));
            dyj=(TT,(D,log(PRES).col(j)));
            dzj=-(yj-logTi.col(j))/dyj;
            if(max(abs(dzj))<1e-10) fin++;
            dzj=max(dzj,-zj/2);
            dzj=min(dzj,(1-zj)/2);
            zj+=dzj;
            nit++;
            if(nit>100) {
                fprintf(stderr,"Error: (star2d) No convergence in find_boundaries\n");
                exit(1);
            }
        }
        zi.setcol(j,zj);
    }
    return zi;
}

matrix star2d::distribute_domains(int ndom,int &conv_new,double p_cc) const {
// conv_new ==> output
// p_cc ==> input
// ndom ==> input
// called by check_map to redistribute domain when conv has changed (CC has appeared or
// disappeared)

    DEBUG_FUNCNAME;
        int j,ndc;
    double p_s;
    matrix pif(ndom,1); // pressure at domain interfaces

printf("Start distribute domains with conv\n");
    p_s=map.leg.eval_00(PRES.row(-1),0)(0);
    if(p_cc==0) { // pcc=0 is default value at start  the core may exist or not (depend on conv
        j=0;
        for(int n=0;n<conv;n++) j+=map.gl.npts[n];
        p_cc=map.leg.eval_00(PRES.row(j),0)(0); // p_cc=central pressure if no CC
    }

/*
// Here star redistribute domains as in "distribute_domains(....)"
    conv_new=conv;
    double drad=(log(p_cc)-log(p_s));
    double dconv=(0.-log(p_cc));
    if(!conv) dconv=0;
    //conv_new=conv==0?0:1;
    conv_new=1;
    if (conv == 0) conv_new=0;
    for(int n=1+conv_new;n<ndom;n++) {
        if(dconv>drad) {
            conv_new++;
            dconv=(0.-log(p_cc))/conv_new;
        } else
            drad=(log(p_cc)-log(p_s))/(n+1-conv_new);
    }

    for(int n=0;n<conv_new;n++) {
        pif(n)=exp(-(n+1)*dconv);
    }
    for(int n=0;n<ndom-conv_new;n++) {
        pif(n+conv_new)=exp(log(p_cc)-(n+1)*drad);
    }
FILE *fic=fopen("pif.txt", "a");
fprintf(fic,"p_cc %e p_s= %e\n",p_cc,p_s);
for (int k=0;k<ndom;k++) fprintf(fic,"%d %e \n",k,pif(k));
*/
	if (conv == 0) {
           for(int n=0;n<ndom;n++) pif(n)=exp((n+1)*log(p_s)/ndom);
        } else {
           ndc=round(ndom*log(p_cc)/log(p_s));
           if (ndc == 0) {ndc=1; printf("ndc = %d\n",ndc);}
           for(int n=0;n<ndc;n++) pif(n)=exp((n+1)*log(p_cc)/ndc);
           double a=log(p_s/p_cc)/(ndom-ndc);
	   double b=log(p_cc)-a*(ndc-1);
 	   for(int n=ndc;n<ndom;n++) pif(n)=exp(a*n+b);
	}
conv_new=ndc;
printf("conv_new = %d ",ndc);
printf("of round = %e\n",ndom*log(p_cc)/log(p_s));

    return pif;
}

matrix star2d::new_distribute_domains(int ndom,matrix p_inter) { //,std::vector<int> zone_type) 
// ndom ==> input
// called by check_map to redistribute domain when conv has changed (CC has appeared or
// disappeared)

    DEBUG_FUNCNAME;
    double p_s;
    matrix pif(ndom,1); // pressure at domain interfaces
    int nzones=zone_type.size();
FILE *fic=fopen("pif.txt", "a");
printf("Start new distribute domains \n");
printf("nzones =  %d \n",nzones);
printf("p_inter nrows %d \n",p_inter.nrows());
printf("p_inter ncols %d \n",p_inter.ncols());
printf("p_inter(nzones-1) %e \n",p_inter(nzones-1));

    p_s=map.leg.eval_00(PRES.row(-1),0)(0);
fprintf(fic,"p_s= %e\n",p_s);

	if (conv == 0) { // Purely radiative case
           for(int n=0;n<ndom;n++) pif(n)=exp((n+1)*log(p_s)/ndom);
        } else { // convection is now taken into account

// We attribute one domain to each zone and then distribute other domains
// so as to minimize the "PRES" jump
// 1/ Check there is enough domains
	if ( ndom < nzones ) printf("Stop not enough domains!!");
	if ( ndom == nzones ) {
         printf("Likely not enough domains, but I try");
        } else { // This is the general case
		matrix dlog(nzones,1);
		int ndom_left=ndom-nzones;
		int iz,izmax,k,ki,ndom_check;
		double dlogmax;
		std::vector <int> ndz,lifi;
               	//izif.clear(); // izif = index of zones interfaces
		for (int k=0;k<nzones;k++) {ndz.push_back(1);izif[k]=1;}
	//for (int k=0;k<nzones;k++) {ndz.push_back(1);izif.push_back(1);}
	//for (int k=0;k<nzones;k++) {ndz.push_back(1);lifi.push_back(1);}
		for (iz=0; iz<nzones; iz++) {
                 if (iz==0) { dlog(iz)=fabs(log(p_inter(0)));
			} else {
			    dlog(iz)=fabs(log(p_inter(iz)/p_inter(iz-1))) ;
			}
		}
//for (int k=0;k<nzones;k++) fprintf(fic,"%d %e \n",k,dlog(k));
// insert an additional domain where the PRES drop is max.
// until no domain is left!
	while (ndom_left >=1) {
		ndom_left=ndom_left-1;
                dlogmax=max(dlog);
		izmax=99;
                for (iz=0; iz<nzones; iz++) {
                    if (fabs(dlog(iz)-dlogmax)<1e-15) izmax=iz;
                }
		if (izmax==99) printf("izmax = %d, ndom_left= %d\n",izmax,ndom_left);
                dlog(izmax)=dlog(izmax)*ndz[izmax]/(ndz[izmax]+1);
		ndz[izmax]=ndz[izmax]+1;
        }             
for (int k=0;k<nzones;k++) fprintf(fic,"zone %d ==> %d domains \n",k,ndz[k]);
// check total number of domains
		ndom_check=0;
		for (iz=0; iz<nzones;iz++) ndom_check=ndom_check+ndz[iz];
		printf("ndom_check = %d\n",ndom_check);
		if (ndom_check != ndom) {
    			std::cerr << "ERROR: ndoms do not match in new_distrib...\n";
    			std::terminate();
		}
// Now each zone has the optimal number of domains. Let's compute the pif
// with a law of the form exp(a*k+b)
		for (k=0;k<ndz[0];k++) pif(k)=exp((k+1)*log(p_inter(0))/ndz[0]);
		ki=ndz[0]-1;
                izif[0]=ki;
			printf("izif (%d) = %d\n",0,ki);
		for (iz=1; iz<nzones; iz++) {
			for (k=ki;k<ki+ndz[iz];k++) {
				double a=log(p_inter(iz)/p_inter(iz-1))/ndz[iz];
				pif(k)=exp(a*(k-ki)+log(p_inter(iz-1)));
			}
 			ki=ki+ndz[iz];
			izif[iz]=ki;
			printf("izif (%d) = %d\n",iz,ki);
		}
// last interface is the surface and not account by the foregoing algo
		pif(ndom-1)=p_inter(nzones-1);
// Now we can tag each domain with its type
		// First zone is special
		for(int n=0;n<=izif[0];n++) domain_type[n] = zone_type[0];
		// Other zones
		for (iz=1; iz<nzones; iz++) {
		  for(int n=izif[iz-1]+1;n<=izif[iz];n++) domain_type[n]=zone_type[iz];
		}
	for (int n=0;n<ndomains;n++) printf("domain_type[%d] = %d\n",n,domain_type[n]);


// Finally we need the index of the interfaces between the zones once the
// domains have been distributed. These are stored in izif.

for (int k=0;k<nzones;k++) fprintf(fic,"k= %d p_inter %e \n",k,p_inter(k));
for (int k=0;k<ndomains;k++) fprintf(fic,"k= %d pif %e \n",k,pif(k));
fclose(fic);
	}
}

    return pif;
}

// First version of find_boundaries but still used in check_map
// to do could be replaced by the new find_boundaries...?
// Find iso "PRES" from the value of pif(idom)
// pif = input param

matrix star2d::find_boundaries_old(matrix pif) const {
    matrix R(pif.nrows(),nth);

    for(int j=0;j<nth;j++) {
        matrix zj,dzj,pj,dpj,TT;
        zj=zeros(pif.nrows(),1);
        int l=0;
        for(int k=0;k<nr;k++) {
            if(PRES(k,j)<pif(l)) {
                zj(l)=z(k-1);
                l++;
                if(l==pif.nrows()) break;
            }
        }
        if(l<pif.nrows())
            for(int k=l;k<pif.nrows();k++) zj(k)=z(-1);
        int fin=0,nit=0;
        while(fin<3) {
            pj=log(map.gl.eval(PRES.col(j),zj,TT));
            dpj=(TT,(D,log(PRES).col(j)));
            dzj=-(pj-log(pif))/dpj;
            if(max(abs(dzj))<1e-10) fin++;
            dzj=max(dzj,-zj/2);
            dzj=min(dzj,(1-zj)/2);
            zj+=dzj;
            nit++;
            if(nit>100) {
                fprintf(stderr,"Error: (star2d) No convergence in find_boundaries\n");
                exit(1);
            }
        }

        R.setcol(j,map.gl.eval(r.col(j),zj));
    }
    return R;
}

void star2d::check_map() {
	DEBUG_FUNCNAME;
	double pcc;
	matrix Rcc,pif, R_inter,p_inter;
//	std::vector <int> zone_type;
	int conv_new;
	remapper *red;
		matrix R(ndomains+1,nth);
//	printf("Start of check_map conv=%d core_convec=%d\n",conv,core_convec);
FILE *fic=fopen("R.txt", "a");
	if(check_CC(pcc,Rcc)!=conv) { // does the following if CC appears or disappears
		red=new remapper(map);
		if(conv) { // CC has disappeared !
			conv=0;
			for(int n=0;n<ndomains;n++) {
				if(n<conv) domain_type[n]=CORE; // unnecessary
				else domain_type[n]=RADIATIVE; // always true
			}
			pif=distribute_domains(ndomains,conv_new);
			R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
			red->set_R(R);
			domain_type.resize(ndomains);
			izif.resize(ndomains);
		} else { // There is a CC that has been discovered by check_conv
		     printf("CORE appeared pcc= %e\n",pcc);
			conv=1;
			pif=distribute_domains(ndomains,conv_new,pcc);
			conv=conv_new; // conv_new may be higher than 1 if big core!
			for(int n=0;n<ndomains;n++) {
				if(n<conv) domain_type[n]=CORE;
				else domain_type[n]=RADIATIVE;
			}
			R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0))); 
			R.setrow(conv,Rcc);
			red->set_R(R);
		}
	} else {
// Check if domain boundaries still verify the rule of not more than
// factor "10" of PRES between two domain boundaries
// called after one Newton iteration
		red=new remapper(map);
		//printf("call remap_domains\n");
		if(!remap_domains(ndomains,*red)) {delete red;return;} //WARNING : prevents from
									//searching new CZ
// If remap_domain false (no remapping needed) then returns
// If remapping true do interp below
// note that remap_domain call distribute_domain(with check option)
	}

	// New stuf to program....
	//find_zones(R_inter, zone_type, p_inter);
	//if (core_convec == 0) return;
	find_zones(R_inter, p_inter);
	pif=new_distribute_domains(ndomains,p_inter); //,zone_type);
	R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
	red->set_R(R);


	if(config.verbose) {printf("Remapping...");fflush(stdout);}
	// Install the new mapping and do interpolation for this mapping
	map=red->get_map();
	interp(red);
for (int k=0;k<=ndomains;k++) fprintf(fic,"k= %d R= %e \n",k,R(k,0));
fclose(fic);
	delete red;
	if(config.verbose) printf("Done\n");

}

// check_CC detects the appearance of a convective core but is not used
// to move the core boundary
// Recall: conv = nb of domains in CC

int star2d::check_CC(double &p_cc,matrix &Rcc) {
    DEBUG_FUNCNAME;
	
if (config.verbose == 19) printf("++++ Start of CHECK_CONVEC, core_convec = %d\n",core_convec);
    if(!core_convec) return 0; // core_covec: input param to disable CC

    if(conv) {
        int j=0;
        for(int n=0;n<conv;n++) j+=map.gl.npts[n]; // number of grid in CC
        if(z(j)<min_core_size) {
        //    if(config.verbose)
          printf("Size(convective core) < min_core_size. Removing...\n");
            return 0;
        } else return conv;
    }
    // else if no CC at calling check_CC
    int i=0;
    while(z(i)<1.05*min_core_size) i++;
// "i" is the first grid point where Schwarzschild is tested (to avoid too small cores)
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
    FILE *fic = fopen("old_schwi.txt", "a");
    for (int k=0;k<nr;k++) fprintf(fic,"i= %d schwi= %e \n",k,schw(k,-1));
    fclose(fic);

    if(map.leg.eval_00(schw.row(i),0)(0)>=0) return 0; 
// if Sch > 0 no CC (or CC too small) and return

    //if(config.verbose)
      printf("Found convective core with check_CC\n");
    if(ndomains==1) {
        fprintf(stderr,"Warning: At least 2 domains are needed to deal with core convection\n");
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
            dzi=-schwi/dschwi;  // Newton method to find the place of the CC-RZ boundary
            if(fabs(dzi)<1e-9) fin++;
            zi+=dzi;
        }
        Rcc(j)=map.gl.eval(r.col(j),zi,TT)(0); // R(theta) of the CC
        pcore(j)=(TT,PRES.col(j))(0);  // PRES on R_CC(theta)
    }

    p_cc=map.leg.eval_00(pcore,0)(0); // p_cc PRES at theta=0
    return 1;

}

matrix star2d::solve_temp_rad() {
    symbolic S;
    S.set_map(map);

    sym t = S.regvar("T");
    sym xi = S.regvar("xi");
    sym eq = div(xi*grad(t))/xi;
    S.set_value("xi", this->opa.xi);

    solver op;
    op.init(map.ndomains, 1, "full");
    op.regvar("T");
    op.set_nr(this->map.npts);

    eq.add(&op, "T", "T");

    matrix rhs = -this->Lambda * this->rho * (this->nuc.eps / this->opa.xi);

    op.bc_bot2_add_l(0, "T", "T", ones(1, this->nth), this->map.D.block(0).row(0));
    rhs.setrow(0, zeros(1, this->nth));

    int ir = this->map.npts[0];
    for (int n=1; n<ndomains; n++) {
        op.bc_top1_add_d(n-1, "T", "T", ones(1, this->nth));
        op.bc_top2_add_d(n-1, "T", "T", -ones(1, this->nth));
        rhs.setrow(ir-1, zeros(1, this->nth));

        op.bc_bot1_add_l(n, "T", "T", ones(1, this->nth), this->map.D.block(n-1).row(-1));
        op.bc_bot2_add_l(n, "T", "T", -ones(1, this->nth), this->map.D.block(n).row(0));
        rhs.setrow(ir, zeros(1, this->nth));
        ir += this->map.npts[n];
    }

    op.bc_top1_add_d(ndomains-1, "T", "T", ones(1, this->nth));
    rhs.setrow(-1, this->T.row(-1));

    op.set_rhs("T", rhs);

    op.solve();

    return op.get_var("T");
}

//int star2d::find_zones(matrix& r_inter, std::vector<int>& zone_type, matrix& p_inter)
int star2d::find_zones(matrix& r_inter, matrix& p_inter) {
    int n = 1;
    matrix schw, dschw;

    matrix T_schw = solve_temp_rad();

#if 0
    FILE *temp_file = fopen("temp.txt", "w");
    for (int i=0; i<nr; i++) {
        for (int j=0; j<nth; j++) {
            fprintf(temp_file, "%e %d %e %e\n", r(i, j), j, T_schw(i, j), T(i, j));
        }
    }
    fclose(temp_file);
#endif
	printf("Start of find zone\n");


// Paco's version
//    schw=-(map.gzz*(D,p)+map.gzt*(p,Dt))*((D,log(T))-eos.del_ad*(D,log(p)))
//        -(map.gzt*(D,p)+map.gtt*(p,Dt))*((log(T),Dt)-eos.del_ad*(log(p),Dt));

// Bertrand version
//    schw = -(map.gzz*(D, p)+map.gzt*(p, Dt))*((D, log(T_schw))-eos.del_ad*(D, log(p))) -
//        (map.gzt*(D, p)+map.gtt*(p, Dt))*((log(T_schw), Dt)-eos.del_ad*(log(p), Dt));
    //schw = -((D, log(T_schw))-eos.del_ad*(D, log(p)));
    schw = log(T_schw);

    FILE *fic = fopen("schwi.txt", "a");
    for (int k=0;k<nr;k++) fprintf(fic,"i= %d schw= %e T=%e, p=%e, na=%e, gzz=%e\n",k,schw(k,-1),T_schw(k,-1),p(k,-1),eos.del_ad(k,-1),map.gzz(k,-1));
    fclose(fic);

    schw.setrow(0, zeros(1, nth));
    schw = schw/r/r;
    schw.setrow(0, zeros(1, nth));
    schw.setrow(0, -(D.row(0), schw)/D(0, 0));

    dschw = (D, schw);


    for (int i=1; i<nr; i++) {
        if (schw(i-1, -1) * schw(i, -1) < 0) {
            n++;
        }
    }

    r_inter = zeros(n, nth);
    p_inter = zeros(n, nth);
    n = 0;
    double last_zi = 0.0;
    for (int i=1; i<nr; i++) {
        if (schw(i-1, -1) * schw(i, -1) < 0) {
            matrix TT;
            double dzi, schwi, dschwi;
            for (int j=0; j<nth; j++) {
                int end = 0, done = 0;
                double zi = z(i-1);
                while (end < 3 && done++ < 100) {
                    schwi = map.gl.eval(schw.col(j), zi, TT)(0);
                    dschwi = (TT, dschw.col(j))(0);
                    dzi = -schwi/dschwi;
                    if (fabs(schwi) < 1e-9) {
                        end = 3;
                    }
                    else {
                        if (fabs(dzi) < 1e-9) end++;
                        zi += dzi;
                    }
                }
                r_inter(n, j) = map.gl.eval(r.col(j), zi)(0);
                p_inter(n, j) = map.gl.eval(PRES.col(j), zi)(0);
                last_zi = zi;
            }
            n++;
        }
    }
    r_inter.setrow(-1, z(-1) * ones(1, nth));
    p_inter.setrow(-1, PRES(-1) * ones(1, nth));

    std::cout << "CONV: " << CONVECTIVE << ", ";
    std::cout << "RAD: " << RADIATIVE << "\n";

    last_zi = 0.0;
    zone_type = std::vector<int>(n+1);
    for (int i=0; i<n+1; i++) {
        double mid_zone = (last_zi+r_inter(i, 0))/2.0;
        double schwi = map.gl.eval(schw, mid_zone)(0);
        if (schwi < 0 && last_zi <=1e-5) {
            zone_type[i] = CORE;
        } else if (schwi < 0 && last_zi >1e-5) {
            zone_type[i] = CONVECTIVE;
        } else {
            zone_type[i] = RADIATIVE;
        }
        printf("ZONE: [%e, %e]: convection (schw: %e): %d\n",
                last_zi, r_inter(i, 0),
                schwi,
                zone_type[i]);
        last_zi = r_inter(i, 0);
    }
// Check if there are contiguous zones and suppress them if true
	int nsz=0;
	redo:
	printf("In find_zones n+1=%d\n",n+1);
	int flag=0;
	for (int iz=0;iz<n;iz++){
		if (zone_type[iz]==zone_type[iz+1]) {
		printf("iz =%d flag=%d\n",iz,flag);
		  for (int k=iz;k<n;k++) {
		    r_inter(k,0)=r_inter(k+1,0);
		    p_inter(k,0)=p_inter(k+1,0);
		    zone_type[k]=zone_type[k+1];
		  }
		  nsz++; // increase the number of suppressed zones
		  flag=1;
		  goto exit_loop;
		}
         }
	exit_loop:
for (int i=0; i<n+1;i++) printf("i=%d, zone_type=%d\n",i,zone_type[i]);
	if (flag) {
	   n=n-1; // reduce the number of zones
	   goto redo; // and check again
	}
	printf("In find_zones nsz=%d \n",nsz);
	printf("In find_zones n+1=%d after reduction\n",n+1);
        //std::vector<int> zone_type_bis(n+1);	
	//for (int i=0;i<n+1;i++) zone_type_bis[i]=zone_type[i];
	zone_type.resize(n+1);
	//for (int i=0;i<n+1;i++) zone_type[i]=zone_type_bis[i];

    return n+1; // n+1 is the number of zones
}
