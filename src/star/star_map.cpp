#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "star.h"
#include "matplotlib.h"

extern "C" {
#include <stdlib.h>
}
// take a model and modify the resolution (nth, ndomains, number of
// points in each domain)
// object == remapper defined here
// used to move from one grid to another
// be careful of discontinuities that might fall in the middle of a
// domain.

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

// Some domains have boundaries imposed by the physics and these
// boundaries cannot be moved (ex. CC & RZ interface).
// domain_type = integer for CC RZ CZ (see star.h)
bool star2d::remap_domains(int ndom, remapper &red) {
    //	Count zones
    int nzones=1;
    std::vector<int> index;
    // Here look for interface between zones of different type.
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
    if(index_new==index&&ndom==ndomains) return false;

    // compute new indices of interfaces between zones and recompute the
    // zeta of the new domains (zif)
    index_new=distribute_domains(ndom,zif);

    red.set_R(zeros(1,nth).concatenate(map.zeta_to_r(zif)));
    for(int n=0;n<nzones;n++)
        // do not update the interface between zones:
        red.set_fixed(index_new[n]+1,index[n]+1);
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

    return true;

}

// if check_only 'True', return the indices of the zones interfaces
// usually used with check_only "false"
std::vector<int> star2d::distribute_domains(int ndom,matrix &zif,bool check_only) const {
    matrix dlogT;
    int nzones=zif.nrows();

    if(nzones>ndom) {
        ester_err("Error: At least %d domains are needed for this model\n",nzones);
    }

    // Calculate Delta(log(T)) in each zone at theta=0
    dlogT.zero(nzones,1);
    for(int n=0;n<nzones;n++) {
        dlogT(n)=-log(map.gl.eval(PRES.col(-1),zif(n,-1))(0));
        if(n) dlogT(n)+=log(map.gl.eval(PRES.col(-1),zif(n-1,-1))(0));
    }

    // Distribute the domains (a limited number) into the different zones
    // so as to decreases optimally the LogT jump between two following
    // domains interfaces.


    int ndomi[nzones]; // Number of domains in each zone
    for(int n=0;n<nzones;n++) ndomi[n]=1;
    // Distribute domains based on dlogT
    for(int n=nzones;n<ndom;n++) {
        double dTi=0;
        int k=0;
        for(int i=0;i<nzones;i++)
            if(dlogT(i)/ndomi[i]>dTi) {
                k=i;
                dTi=dlogT(i)/ndomi[i];
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
            // evaluate PRES on the interfaces bounding the domain:
            for(int j=0;j<nth;j++) {
                logT0(j)=log(map.gl.eval(PRES.col(j),zif(izone-1,j))(0));
                logT1(j)=log(map.gl.eval(PRES.col(j),zif(izone,j))(0));
            }
            logTi.setrow(n,logT0);
        } else {
            // For the interfaces of domains inside the same zone attribute a value
            // of the log(temperature) as a linear function of the domain rank
            logTi.setrow(n,logT0+(logT1-logT0)*((double) k)/ndomi[izone]);
        }
    }

    // A priori the "temperature" (in fact the variable PRES)  on the
    // interfaces of zones is not constant; hence the interfaces of domains
    // inside a zone have not a constant temperature.  But we like to have
    // isothermal interfaces of domain inside a zone for (eg) numerical stability
    // so the option ifdef....

#ifdef T_CONSTANT_DOMAINS
    logTi=map.leg.eval_00(logTi,0)*ones(1,nth);
#endif
    zif_new.setblock(0,ndom-2,0,-1,find_boundaries(logTi));

    for(int izone=0;izone<nzones;izone++)
        zif_new.setrow(index[izone],zif.row(izone));

    zif=zif_new;
    return index;

}

// give a function logTi(theta) and find the associated zeta_i(theta_k)
// of the surface
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

matrix star2d::distribute_domains(int ndom,int &conv_new,double p_cc) const {
    // conv_new ==> output
    // p_cc ==> input
    // ndom ==> input
    // called by check_map to redistribute domain when conv has changed (CC
    // has appeared or disappeared)

    int j;
    double p_s;

    conv_new=conv;
    p_s=map.leg.eval_00(PRES.row(-1),0)(0);

    // pcc=0 is default value at start  the core may exist or not (depend on
    // conv)
    if(p_cc==0) {
        j=0;
        for(int n=0;n<conv;n++) j+=map.gl.npts[n];
        p_cc=map.leg.eval_00(PRES.row(j),0)(0);
    }

    // Here star redistribute domains as in "distribute_domains(other arg...)"
    double drad=(log(p_cc)-log(p_s));
    double dconv=(0.-log(p_cc));
    if(!conv) dconv=0;
    conv_new=conv==0?0:1;
    for(int n=1+conv_new;n<ndom;n++) {
        if(dconv>drad) {
            conv_new++;
            dconv=(0.-log(p_cc))/conv_new;
        } else
            drad=(log(p_cc)-log(p_s))/(n+1-conv_new);
    }

    matrix pif(ndom,1); // pressure at domain interfaces
    for(int n=0;n<conv_new;n++) {
        pif(n)=exp(-(n+1)*dconv);
    }
    for(int n=0;n<ndom-conv_new;n++) {
        pif(n+conv_new)=exp(log(p_cc)-(n+1)*drad);
    }

    return pif;
}

// First version of find_boundaries but still used in check_map
// to do: could be replaced by the new find_boundaries...???
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
                ester_err("No convergence in find_boundaries\n");
            }
        }

        R.setcol(j,map.gl.eval(r.col(j),zj));
    }
    return R;
}

void star2d::check_map() {
    double pcc;
    matrix Rcc,pif;
    int conv_new;
    remapper *red;

    if(check_convec(pcc,Rcc)!=conv) {
        // does the following if a CC appears or disappears
        matrix R(ndomains+1,nth);
        red=new remapper(map);
        if(conv) { // CC has disappeared !
            conv=0;
            for(int n=0;n<ndomains;n++) {
                if(n<conv) domain_type[n]=CORE; //unnecessary
                else domain_type[n]=RADIATIVE; // always true
            }
            pif=distribute_domains(ndomains,conv_new);
            R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
            red->set_R(R);
            domain_type.resize(ndomains);
        } else { // There is a CC that has been discovered by check_conv
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
        if(!remap_domains(ndomains,*red)) {delete red;return;}
    }
    if(config.verbose) {printf("Remapping...");fflush(stdout);}
    // Install the new mapping and do interpolation for this mapping
    map=red->get_map();
    interp(red);
    delete red;
    if(config.verbose) printf("Done\n");

}

// check_convec detects the appearance of a convective core but is not
// used to move the core boundary

int star2d::check_convec(double &p_cc,matrix &Rcc) {
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

    p_cc=map.leg.eval_00(pcore,0)(0); // p_cc=PRES at theta=0
    return 1;

}


