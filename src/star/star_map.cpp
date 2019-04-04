// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "star.h"
#include "matplotlib.h"

extern "C" {
#include <stdlib.h>
}

void star2d::new_check_map() {
        matrix pif,R_inter,p_inter;
	double seuil;
        remapper *red;
        matrix R(ndomains+1,nth);
	R=zeros(ndomains+1,nth);
	if (glit == 1 && config.input_file==0) {
// At first iteration set the number of zone_type to 1
                zone_type = std::vector<int>(1);
	}

	if (n_essai == 0) seuil=1e-2;
	if (n_essai > 0) seuil=1e-2;
//printf("global errors %e %e\n",prev_global_err,global_err);

// Check the "PRESSURE" drop from bot CZ to surface
// At the beginning izif=0, thus n=0 and j0=number of pts of 1st domain
// hence p_drop > p_dm, always.
	//printf("izif size %d\n",izif.size());
	//for (int iz=0;iz<ndomains-1;iz++) printf(" %d ",izif[iz]);
	//printf("\n");
	//printf("nzones %d\n",nzones);
	int n=izif[nzones-2],j0,k; // n is the index of the domain just below the interface
        for(k=0,j0=0;k<n+1;k++) j0+=map.gl.npts[k]; // j0 is the radial index of the interface
	
	double p_drop=PRES(-1,0)/PRES(j0,0);
	int ndom=ndomains-n;
	double p_dm=exp(ndom*log(PRES(-1,0))/ndomains); //expected drop
	//printf("j0 = %d ndom= %d\n",j0,ndom);
	//printf("pdrop = %e\n",p_drop);
	//printf("pdm   = %e\n",p_dm);

	if (p_dm > p_drop) { // PRES drop is too important in CZ, redist domain
           if (details) printf("check_map: PRES drop high : REDISTRIBUTE DOMAINS\n");
    	   p_inter = zeros(nzones, nth);
	   for (int iz=0;iz<nzones-1;iz++) { // we scan the zone interfaces
                n=izif[iz]; // n is the index of the domain just below the interface
                for(k=0,j0=0;k<n+1;k++) j0+=map.gl.npts[k];
	  	p_inter(iz,0)=PRES(j0,0);
	   }
	   p_inter(nzones-1,0)=PRES(-1,0);
	   pif=New_distribute_domains(ndomains,p_inter);
	   R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
	   red=new remapper(map);
	   red->set_R(R);
	   map=red->get_map(); // generate r,r_z,...
	   interp(red);
           delete red;
	   return;
	}


	if (global_err < seuil) { // look for new convective regions and eventually remap the star
	   // Find the zone boundaries and the associated pressures
	   // and output zone_type as global var.
 	   printf("seuil %e\n",seuil);
	   int nzones_av=zone_type.size();

	   nzones=find_zones(R_inter, p_inter); // defines R,p_inter

	   int nzones_ap=zone_type.size();
	   if (nzones_av == nzones_ap) {
		n_essai=n_essai+1;
		  if (details) printf("n_essai=%d\n",n_essai);
		  //if (global_err < prev_global_err || global_err<2e-3) {
		  if (global_err < 2e-3 || global_err<2e-3) {
		  if (details) printf("Number of zones unchanged: do nothing\n");
		  return;
		}
	   }

// The global error is low enough, find_zone was called
// and the number of zones has changed so we redistribute the domains
// Also output izif (index of zone interface) as a global var.

           if (details) printf("check_map: REDISTRIBUTE DOMAINS\n");
	   n_essai=0;
	   pif=New_distribute_domains(ndomains,p_inter);
	   // fill the matrix R with the radii of the first interface
	   // to the lat interface, so 1==> -2, for all theta 0==>-1
	   R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
	   red=new remapper(map);
	   red->set_R(R);
	// Install the new mapping and do interpolation for this mapping
	   map=red->get_map(); // generate r,r_z,...
	   interp(red);
	   printf("DISTRIBUTION OF DOMAINS IN ZONES nzones=%d\n",nzones);
	   for (int iz=0;iz<nzones;iz++) {
		if (iz == 0) ndom=izif[iz]+1;
		if (iz != 0) ndom=izif[iz]-izif[iz-1];
		printf(" iz=%d ndom=%d, ",iz,ndom);
	   }
	   printf("\n");

if (details) {
FILE *fic=fopen("new_R.txt", "a");
fprintf(fic,"it= %d\n",glit);
for (int k=0;k<=ndomains;k++) fprintf(fic,"k= %d R= %e \n",k,R(k,0));
fclose(fic);
}
	   delete red;

	} else return; // else do nothing
}
//-----------------------------------------------------------------------------

matrix star2d::New_distribute_domains(int ndom,matrix p_inter) {
// ndom ==> input
// called by new_check_map to redistribute domain

    double p_s;
    matrix pif(ndom,1); // pressure at domain interfaces
    int nzo=zone_type.size(); // zone_type is initialized by find_zones
FILE *fic=fopen("New_pif.txt", "a");
if (details) {
printf("Start NEW distribute domains \n");
printf("nzs =  %d \n",nzo);
printf("p_inter nrows %d \n",p_inter.nrows());
printf("p_inter ncols %d \n",p_inter.ncols());
for (int k=0; k<nth;k++) printf("p_inter(0,%d)= %e, p_inter(1,%d)= %e \n",k,p_inter(0,k),k,p_inter(nzo-1,k));
}

    p_s=map.leg.eval_00(PRES.row(-1),0)(0);
fprintf(fic,"p_s= %e\n",p_s);

	if (nzo == 1) { // A single domain (purely radiative)
           for(int n=0;n<ndom;n++) pif(n)=exp((n+1)*log(p_s)/ndom);
        } else { // convection is now taken into account

// We attribute one domain to each zone and then distribute other domains
// so as to minimize the "PRES" jump
// 1/ Check there is enough domains
	if ( ndom < nzo ) printf("Stop not enough domains!!");
	if ( ndom == nzo ) {
         printf("Likely not enough domains, but I try");
        } else { // This is the general case
		matrix dlog(nzo,1);
		int ndom_left=ndom-nzo;
		int iz,izmax,k,ki,ndom_check;
		double dlogmax;
		std::vector <int> ndz,lifi;
		for (int k=0;k<nzo;k++) {ndz.push_back(1);izif[k]=1;}
		for (iz=0; iz<nzo; iz++) {
                 if (iz==0) { dlog(iz)=fabs(log(p_inter(0)));
			} else {
			    dlog(iz)=fabs(log(p_inter(iz)/p_inter(iz-1))) ;
			}
		}
	if (details) printf("dlogp zone 1 = %e, dlogp zone 2 = %e\n",dlog(0),dlog(1));
// insert an additional domain where the PRES drop is max.
// until no domain is left!
		while (ndom_left >=1) {
		  ndom_left=ndom_left-1;
                  dlogmax=max(dlog);
		  izmax=99;
                  for (iz=0; iz<nzo; iz++) {
                    if (fabs(dlog(iz)-dlogmax)<1e-15) izmax=iz;
                  }
		  if (izmax==99) printf("izmax = %d, ndom_left= %d\n",izmax,ndom_left);
                  dlog(izmax)=dlog(izmax)*ndz[izmax]/(ndz[izmax]+1);
		  ndz[izmax]=ndz[izmax]+1;
        	}             
for (int k=0;k<nzo;k++) fprintf(fic,"zone %d ==> %d domains \n",k,ndz[k]);
// check total number of domains
		ndom_check=0;
		for (iz=0; iz<nzo;iz++) ndom_check=ndom_check+ndz[iz];
		if (details) printf("ndom_check = %d\n",ndom_check);
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
		if (details) printf("izif (%d) = %d\n",0,ki);
           // First zone done
	   // Deal with the following zones
		for (iz=1; iz<nzo; iz++) {
			for (k=ki;k<ki+ndz[iz];k++) {
				double a=log(p_inter(iz)/p_inter(iz-1))/ndz[iz];
				pif(k)=exp(a*(k-ki)+log(p_inter(iz-1)));
			}
 			ki=ki+ndz[iz];
			izif[iz]=ki;
			if (details) printf("izif (%d) = %d\n",iz,ki);
		}
		
           // All zones done
// last interface is the surface and not account by the foregoing algo
		pif(ndom-1)=p_inter(nzo-1);
// Now we can tag each domain with its type
		// First zone is special
		for(int n=0;n<=izif[0];n++) domain_type[n] = zone_type[0];
		// Other zones
		for (iz=1; iz<nzo; iz++) {
		  for(int n=izif[iz-1]+1;n<=izif[iz];n++) domain_type[n]=zone_type[iz];
		}
// Compute the number of domains in the core
	nd_core=0;
	for (int n=0;n<ndomains;n++) {
		if (details) printf("domain_type[%d] = %d\n",n,domain_type[n]);
		if (domain_type[n] == CORE) nd_core++;
	}
		if (details) printf("nd_core=%d\n",nd_core);


// Finally we need the index of the interfaces between the zones once the
// domains have been distributed. These are stored in izif.

if (details) for (int k=0;k<nzo;k++) fprintf(fic,"k= %d p_inter %e \n",k,p_inter(k));
if (details) for (int k=0;k<ndomains;k++) fprintf(fic,"k= %d pif %e \n",k,pif(k));
fclose(fic);
	}
}

    return pif;
}

//-----------------------------------------------------------------------------

// take a model and modify the resolution (nth, ndomains, number of
// points in each domain) 
// object == remapper defined here
// used to move from one grid to another
// be careful of discontinuities that might fall in the middle of a
// domain.

void star2d::remap(int ndom_in,int *npts_in,int nth_in,int nex_in) {
details=1;
    remapper red(map); // declaration object of class remapper 

    if (details) printf("    Enter remap in star_map 1\n");
	if (details) printf("ndom_old=%d, ndom_new=%d\n",ndomains,ndom_in);
    red.set_ndomains(ndom_in);
    red.set_npts(npts_in);
    red.set_nt(nth_in);
    red.set_nex(nex_in);
    if (details) printf("    sets done in remap in star_map 1\n");

    if(ndom_in!=ndomains) 
        remap_domains(ndom_in,red); // the new R_i are now known

    map=red.get_map(); // update the mapping
    if (details) for (int n=0;n<ndomains;n++)
            printf("%d, domain_type = %d, r=%e\n",n,domain_type[n],map.gl.xif[n]);
    if (details) printf("mapping updated in remap in star_map 1\n");
    interp(&red); // interpolate the variable on the new update

	if (details) printf("    Leave remap in star_map 1\n");
}

//-----------------------------------------------------------------------------
// Some domains have boundaries imposed by the physics and these
// boundaries cannot be moved (ex. CC & RZ interface).
// domain_type = integer for CC RZ CZ (see star.h)
bool star2d::remap_domains(int ndom, remapper &red) {
//	Count zones
details=1;
    int nzo=1;
    std::vector<int> index;	
// Here look for interface between zones of different type.
if (details) printf("++++++ start of remap_domains\n");
    for(int n=1,type=domain_type[0];n<ndomains;n++) {
        if(domain_type[n]!=type) {
            index.push_back(n-1);
            nzo++;
            type=domain_type[n];
        }
    }
// index == index of interface sepa. zone of diff type from center to
// surface
    index.push_back(ndomains-1);

// zif = zeta of interfaces of all zones before remapping 
    matrix zif(nzo,nth);
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
     if (details) printf("quit remap_domain with nothing changed\n");
     return false;
    }

// compute new indices of interfaces between zones and recompute the zeta of the new
// domains (zif)
if (details) printf("++++++ remap_domains and redistribute domains\n");
    index_new=distribute_domains(ndom,zif);

// Now that the zif are known we can recompute the radii R of the new domains inside
// the zones
    red.set_R(zeros(1,nth).concatenate(map.zeta_to_r(zif)));
    for(int n=0;n<nzo;n++) 
        red.set_fixed(index_new[n]+1,index[n]+1); // do not update the interface between zones
// Update domain_type (each domain is tagged with the zone type)
    std::vector<int> domain_type_new(ndom,0);
    for(int n=0,izone=0;n<ndom;n++) {
        domain_type_new[n]=domain_type[index[izone]];
        if(n==index_new[izone]) izone++;
    }
    domain_type=domain_type_new;

// Update index of zone interface (izif)
    izif.resize(ndom);
	for(int n=0,izone=0;n<ndom-1;n++) {
	   if (domain_type[n] != domain_type[n+1]) {
		izif[izone]=n;
		izone++;
	   }
	}
    conv=0;
    int n=0;
    while(domain_type[n++]==CORE) conv++; // update of the conv variable

if (details) printf("++++++ remap_domains end with redistributed domains\n");
    return true;
}
//-----------------------------------------------------------------------------

// if check_only 'True', return the indices of the zones interfaces
// usually used with check_only "false" to 
std::vector<int> star2d::distribute_domains(int ndom,matrix &zif,bool check_only) const {
    matrix dlogT;
    int nzo=zif.nrows();

    if(nzo>ndom) {
        fprintf(stderr,"Error: At least %d domains are needed for this model\n",nzo);
        exit(1);
    }

    // Calculate Delta(log(T)) in each zone at theta=0
    dlogT.zero(nzo,1);
    for(int n=0;n<nzo;n++) {
        dlogT(n)=-log(map.gl.eval(PRES.col(-1),zif(n,-1))(0));
        if(n) dlogT(n)+=log(map.gl.eval(PRES.col(-1),zif(n-1,-1))(0));
    }

// Distribute the domains (a limited number) into the different zones so as to
// decrease optimally the LogT jump between two following domains interfaces.

    int ndomi[nzo]; // Number of domains in each zone
    for(int n=0;n<nzo;n++) ndomi[n]=1;
    // Distribute domains based on dlogT
    for(int n=nzo;n<ndom;n++) {
        double dTi=0;
        int k=0;
        for(int i=0;i<nzo;i++) {
            if(dlogT(i)/ndomi[i]>dTi) {
                k=i;
                dTi=dlogT(i)/ndomi[i];
            }
        }
        ndomi[k]++;
    }

// Calculate the new indices of zone interfaces
    std::vector<int> index;
    index.resize(nzo);

    for(int n=0,k=1,izone=0;n<ndom-1;n++,k++) {
        if(k==ndomi[izone]) {
            index[izone]=n;
            k=0;
            izone++;
        }
    }
    index[nzo-1]=ndom-1;
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

    for(int izone=0;izone<nzo;izone++) 
        zif_new.setrow(index[izone],zif.row(izone));

    zif=zif_new;
    return index;

}
//-----------------------------------------------------------------------------

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
       //fprintf(stderr,"Error: (star2d) No convergence in find_boundaries\n");
       //exit(1);
                plt::clf();
                plt::plot(r, log(T), "$T$");
                for (int i=0; i<logTi.nrows(); i++) {
                    plt::axvline(zj(i));
                    plt::axhline(logTi(i));
                    LOGE("ri%d zj=%e, dzj=%e\n", i, zj(i), dzj(i));
                }
                LOGE("No convergence in find_boundaries");
                plt::show(true);
                ester_err("No convergence in find_boundaries\n");
            }
        }
        zi.setcol(j,zj);
    }
    return zi;
}

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------


matrix star2d::solve_temp_rad() {
    symbolic S;
    S.set_map(map);

    sym t = S.regvar("T");
    sym xi = S.regvar("xi");
    sym eq = div(xi*grad(t))/xi;
    S.set_value("xi", opa.xi);

    solver op;
    op.init(map.ndomains, 1, "full");
    op.regvar("T");
    op.set_nr(map.npts);

    eq.add(&op, "T", "T");

    matrix rhs = -Lambda * rho * (nuc.eps / opa.xi);

    op.bc_bot2_add_l(0, "T", "T", ones(1, nth), map.D.block(0).row(0));
    rhs.setrow(0, zeros(1, nth));
    //op.bc_bot2_add_d(0, "T", "T", ones(1, nth));
    //rhs.setrow(0, ones(1, nth));

    int ir = map.npts[0];
    for (int n=1; n<ndomains; n++) {
        op.bc_top1_add_d(n-1, "T", "T", ones(1, nth));
        op.bc_top2_add_d(n-1, "T", "T", -ones(1, nth));
        rhs.setrow(ir-1, zeros(1, nth));

        op.bc_bot1_add_l(n, "T", "T", ones(1, nth), map.D.block(n-1).row(-1));
        op.bc_bot2_add_l(n, "T", "T", -ones(1, nth), map.D.block(n).row(0));
        rhs.setrow(ir, zeros(1, nth));
        ir += map.npts[n];
    }

    op.bc_top1_add_d(ndomains-1, "T", "T", ones(1, nth));
    rhs.setrow(-1, T.row(-1));

    op.set_rhs("T", rhs);

    op.solve();

    return op.get_var("T");
}
//-----------------------------------------------------------------------------

//int star2d::find_zones(matrix& r_inter, std::vector<int>& zone_type, matrix& p_inter)
int star2d::find_zones(matrix& r_inter, matrix& p_inter) {
    double last_zi, zi;
    matrix schw, dschw,bv2,schwBP,schwFEL;

    matrix T_rad = solve_temp_rad();

	if (details) {
            printf("Start of find zone, min_core= %e\n",min_core_size);
	    FILE *coco=fopen("les_z.txt", "w");
	    for (int i=0;i<nr;i++) fprintf(coco,"%d  %e\n",i,z(i));
	    fclose(coco);
	}

    int itest_sch=0;
    while(r(itest_sch,-1)<1.05*min_core_size) itest_sch++;
// "itest_sch" is the first grid point where Schwarzschild is tested
// (to avoid too small cores)

/** Schwarzschild criterion
    schw = -(grad P).(grad s)
    if schw > 0 then stable
    if schw < 0 then unstable = convection zone
**/

// Paco's version
      schwFEL=-(map.gzz*(D,p)+map.gzt*(p,Dt))*((D,log(T))-eos.del_ad*(D,log(p)))
        -(map.gzt*(D,p)+map.gtt*(p,Dt))*((log(T),Dt)-eos.del_ad*(log(p),Dt));
      schwFEL.setrow(0, zeros(1, nth));

// We evaluate the Schwarzschild criterion from the radiative solution T_rad
      schwBP = -(map.gzz*(D, p)+map.gzt*(p, Dt))*((D, T_rad)/T_rad-eos.del_ad*(D, log(p))) -
          (map.gzt*(D, p)+map.gtt*(p, Dt))*((T_rad, Dt)/T_rad-eos.del_ad*(log(p), Dt));
      schwBP.setrow(0, zeros(1, nth));


    //for (int k=0;k<nr;k++) fprintf(fic,"i= %d schw= %e T=%e, p=%e, na=%e, gzz=%e\n",k,schw(k,-1),T_rad(k,-1),p(k,-1),eos.del_ad(k,-1),map.gzz(k,-1));

    //schw = schw/r/r;
    //schw.setrow(0, zeros(1, nth));
    //schw.setrow(0, -(D.row(0), schw)/D(0, 0));
	//schwarz=schw;
	////printf("size of schwarz %d,%d\n",schwarz.nrows(),schwarz.ncols());

    //for (int k=0;k<nr;k++) fprintf(fic,"i= %d schw= %e T_rad=%e, T=%e\n",k,schw(k),T_rad(k),T(k));

// Compute the square of the BV frequency
	bv2=(D,p)/p/eos.G1-(D,rho)/rho;
	bv2(0,-1)=0.;
/*
matrix ss=entropy();
FILE *entrop;
entrop=fopen("sm_entropie.txt","a");
fprintf(entrop,"entropy it = %d\n",glit);
for (int k=0;k<nr;k++) fprintf(entrop,"%d %e %e\n",k,map.r(k),ss(k));
fclose(entrop);
int k_rad=0;
    for (int k=0;k<nr;k++) {
	if (ss(k,-1)-ss(nr-1,-1) < 0.0 ) {
           k_rad=k; // store the last point where entropy is below surface entropy
        }
    }
	   printf("k_rad = %d -------------------------\n",k_rad);
*/

// filtrage chebyshev
	matrix sp_bv2,bv2f;
	sp_bv2=(map.gl.P,bv2);
	int j0=0, n_filter=10;
	for (int idom=0; idom< ndomains;idom++) {
	   int ndom=map.gl.npts[idom];
           int j1=j0+ndom-1;
	   for (int n=j0+n_filter; n<=j1; n++) sp_bv2(n,-1) = 0.;
	   j0+=ndom;
	}
	bv2f = (map.gl.P1, sp_bv2); 

// filtrage brutal
    schw=bv2f;
    for (int k=0;k<nr;k++) {
	if (fabs(bv2f(k,-1)) < 1.0e-3 ) schw(k,-1)=-1e-3;
    }
// pas de filtrage
	//schw=schwBP;
	schw=schwFEL;



    FILE *fic = fopen("schwi.txt", "a");
	fprintf(fic,"it= %d\n",glit);
    for (int k=0;k<nr;k++) fprintf(fic,"%d r=%e schw=%e bv2=%e %e\n",k,r(k,-1),schw(k,-1),bv2(k,-1),bv2f(k,-1));
    fclose(fic);

    //dschw = (D, schw);

    int nc=999,nl;
    int n=1;
    for (int i=itest_sch; i<nr; i++) {
        if (schw(i-1, -1) * schw(i, -1) < 0 ) {
            n++;
	    nc=i; // memo of the last sign change
        }
    }
// We need to remove the outest (surface) convection zone (generally very thin)
// since we want models that terminate with a radiative zone
	nl=nc;
	if (details) printf(" av. n= %d, nl=%d\n",n,nl);

	if (schw(nr-2,-1) < 0) {
		printf("Last domain Conv  schw<0\n");
		//n=n-1;
		//nl=nc-1;
	}

	if (details) printf(" ap. n= %d, nl=%d\n",n,nl);


// We introduce les_zi that are the zeta of the interfaces
// We require that les_zi[0]=0. (center)
// If there are n zones there are n+1 interfaces which encompass centre and
// surface.

    std::vector<double> les_zi(n+2);
    for (int i=0; i< n+2;i++) les_zi[i]=0.;

// We first detect the number of true interfaces using the last latitude
// point (-1), which is the closest one to the pole.
    n = 0;
    for (int i=itest_sch; i<nr; i++) {  
        if (schw(i-1, -1) * schw(i, -1) < 0 ) {
                n++;
		zi=z(i-1)-schw(i-1,-1)*(z(i)-z(i-1))/(schw(i,-1)-schw(i-1,-1));
		les_zi[n]=zi;
        }
    }
les_zi[n+1]=1.;
if (details) printf("nb dinterface detecte= %d\n",n);
if (details) for (int i=0; i<n+2;i++) printf("i= %d  z= %e\n",i,les_zi[i]);

// Now we select the thick zones by removing the thin ones (less than
// 0.2%)
int is=0;
for (int i=0; i< n+1;i++) {
if (les_zi[i+1]-les_zi[i] < 0.001) {
   les_zi[i+1]=0.;
   les_zi[i]=0.;
   is=is+2;
   i++;
   }
}

if (details) printf("number of interf suppressed %d \n",is);
if (details) for (int i=0; i< n+2;i++) printf("i= %d  z= %e\n",i,les_zi[i]);

// Now we count the number of true interfaces left ==> n_interf
// and recollect their positions ==> the_zi
int n_interf=0;
std::vector<double> the_zi(n);

for (int i=0; i< n+2;i++) {
if (les_zi[i] !=0. && les_zi[i] !=1.) {
	the_zi[n_interf]=les_zi[i];
	   n_interf++;
          }
}
if (details) printf("Number of interfaces calcule %d \n",n_interf);
matrix bb=ones(2,1);
bb(0,0)=1; bb(1,0)=n_interf;
n_interf=min(bb);
if (details) printf("number of interfaces calcule %d \n",n_interf);

// Now we set the values of r_inter and p_inter
// which include the true surface
    r_inter = zeros(n_interf+1, nth);
    p_inter = zeros(n_interf+1, nth);
    for (int j=0; j<nth;j++) {
      for (int i=0; i< n_interf;i++) {
	   zi=the_zi[i];
           r_inter(i, j) = map.gl.eval(r.col(j), zi)(0);
           p_inter(i, j) = map.gl.eval(PRES.col(j), zi)(0);
      }
    }
// We add the surface as the last r_inter
    r_inter.setrow(-1, z(-1) * ones(1, nth));
    p_inter.setrow(-1, PRES(-1) * ones(1, nth));
//if (details) printf("p_inter i= %d, %e\n",n_interf,p_inter(-1,0));

	if (details) { 
           std::cout << "CONV: " << CONVECTIVE << ", ";
           std::cout << "CORE: " << CORE << ", ";
           std::cout << "RAD: " << RADIATIVE << "\n";
        }

    int j=-1;
    zone_type = std::vector<int>(n_interf+1);
    double mid_zone = r_inter(0, j)/2.0; // middle of first zone
    double schwi = map.gl.eval(schw, mid_zone)(0);
    if (schwi < 0) {
            zone_type[0] = CORE;
    } else {
            zone_type[0] = RADIATIVE;
    }
    printf("ZONE: [%e, %e], type: %d\n",
                0., r_inter(0, j), zone_type[0]);

    last_zi = r_inter(0, j);
    for (int i=1; i<n_interf+1; i++) {
        if ( zone_type[i-1] == CORE || zone_type[i-1] == CONVECTIVE) {
            zone_type[i] = RADIATIVE;
        } else if (zone_type[i-1] == RADIATIVE) {
            zone_type[i] = CONVECTIVE;
        }
        printf("ZONE: [%e, %e], type: %d\n",
                last_zi, r_inter(i, j), zone_type[i]);
        last_zi = r_inter(i, j);
    }

//exit(0);
if (details) printf("End of find_zones: number of zones =%d \n",n_interf+1);

    return n_interf+1; // return the number of zones
}
//-----------------------------------------------------------------------------
