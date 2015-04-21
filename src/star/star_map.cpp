#include "ester-config.h"
#include "star.h"
#include "utils.h"

extern "C" {
#include <stdlib.h>
}

void star2d::remap(int ndom_in,int *npts_in,int nth_in,int nex_in) {
    DEBUG_FUNCNAME;
    remapper red(map);

    red.set_ndomains(ndom_in);
    red.set_npts(npts_in);
    red.set_nt(nth_in);
    red.set_nex(nex_in);

    if(ndom_in!=ndomains) 
        remap_domains(ndom_in,red);

    map=red.get_map();
    interp(&red);
}

bool star2d::remap_domains(int ndom, remapper &red) {
    DEBUG_FUNCNAME;
    //	Count zones
    int nzones=1;
    std::vector<int> index;	
    for(int n=1,type=domain_type[0];n<ndomains;n++) {
        if(domain_type[n]!=type) {
            index.push_back(n-1);
            nzones++;
            type=domain_type[n];
        }
    }
    index.push_back(ndomains-1);

    matrix zif(nzones,nth);
    for(int n=0,j0=-1,izone=0;n<ndomains;n++) {
        j0+=map.gl.npts[n];
        if(n==index[izone]) {
            zif.setrow(izone,z(j0)*ones(1,nth));
            izone++;
        }
    }

    std::vector<int> index_new;

    index_new=distribute_domains(ndom,zif,true);

    if(index_new==index&&ndom==ndomains) return false;

    index_new=distribute_domains(ndom,zif);

    red.set_R(zeros(1,nth).concatenate(map.zeta_to_r(zif)));
    for(int n=0;n<nzones;n++) 
        red.set_fixed(index_new[n]+1,index[n]+1);
    // Update domain_type
    std::vector<int> domain_type_new(ndom,0);
    for(int n=0,izone=0;n<ndom;n++) {
        domain_type_new[n]=domain_type[index[izone]];
        if(n==index_new[izone]) izone++;
    }

    domain_type=domain_type_new;
    conv=0;
    int n=0;
    while(domain_type[n++]==CORE) conv++;

    return true;

}

std::vector<int> star2d::distribute_domains(int ndom,matrix &zif,bool check_only) const {
    DEBUG_FUNCNAME;
    matrix dlogT;
    int nzones=zif.nrows();

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

    // Find new boundaries
    matrix zif_new(ndom,nth);
    matrix logTi(ndom-1,nth);
    matrix logT0,logT1(1,nth);
    logT0=zeros(1,nth);
    for(int j=0;j<nth;j++) logT1(j)=log(map.gl.eval(PRES.col(j),zif(0,j))(0));
    for(int n=0,k=1,izone=0;n<ndom-1;n++,k++) {
        if(k==ndomi[izone]) {
            k=0;
            izone++;
            for(int j=0;j<nth;j++) {
                logT0(j)=log(map.gl.eval(PRES.col(j),zif(izone-1,j))(0));
                logT1(j)=log(map.gl.eval(PRES.col(j),zif(izone,j))(0));
            }
            logTi.setrow(n,logT0);
        } else {
            logTi.setrow(n,logT0+(logT1-logT0)*((double) k)/ndomi[izone]);		
        }
    }

#ifdef T_CONSTANT_DOMAINS
    logTi=map.leg.eval_00(logTi,0)*ones(1,nth);
#endif
    zif_new.setblock(0,ndom-2,0,-1,find_boundaries(logTi));

    for(int izone=0;izone<nzones;izone++) 
        zif_new.setrow(index[izone],zif.row(izone));

    zif=zif_new;
    return index;

}

matrix star2d::find_boundaries(const matrix &logTi) const {
    DEBUG_FUNCNAME;
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
                ester_err("(star2d) No convergence in find_boundaries");
                exit(1);
            }
        }
        zi.setcol(j,zj);
    }
    return zi;
}

matrix star2d::distribute_domains(int ndom,int &conv_new,double p_cc) const {
    DEBUG_FUNCNAME;
        int j;
    double p_s;

    conv_new=conv;
    p_s=map.leg.eval_00(PRES.row(-1),0)(0);
    if(p_cc==0) {
        j=0;
        for(int n=0;n<conv;n++) j+=map.gl.npts[n];
        p_cc=map.leg.eval_00(PRES.row(j),0)(0);
    }

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

    matrix pif(ndom,1);
    for(int n=0;n<conv_new;n++) {
        pif(n)=exp(-(n+1)*dconv);
    }
    for(int n=0;n<ndom-conv_new;n++) {
        pif(n+conv_new)=exp(log(p_cc)-(n+1)*drad);
    }

    return pif;
}

matrix star2d::find_boundaries_old(matrix pif) const {
    DEBUG_FUNCNAME;
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
                ester_err("Error: (star2d) No convergence in find_boundaries\n");
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
    matrix Rcc,pif;
    int conv_new;
    remapper *red;

    if(check_convec(pcc,Rcc)!=conv) {
        matrix R(ndomains+1,nth);
        red=new remapper(map);
        if(conv) {
            conv=0;
            for(int n=0;n<ndomains;n++) {
                if(n<conv) domain_type[n]=CORE;
                else domain_type[n]=RADIATIVE;
            }
            pif=distribute_domains(ndomains,conv_new);
            R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0)));
            red->set_R(R);
            domain_type.resize(ndomains);
        } else {
            conv=1;
            pif=distribute_domains(ndomains,conv_new,pcc);
            conv=conv_new;
            for(int n=0;n<ndomains;n++) {
                if(n<conv) domain_type[n]=CORE;
                else domain_type[n]=RADIATIVE;
            }
            R.setblock(1,-2,0,-1,find_boundaries_old(pif.block(0,-2,0,0))); 
            R.setrow(conv,Rcc);
            red->set_R(R);
        }
    } else {
        red=new remapper(map);
        if(!remap_domains(ndomains,*red)) {
            delete red;
            return;
        }
    }
    if(config.verbose) {printf("Remapping...");fflush(stdout);}
    map=red->get_map();
    interp(red);
    delete red;
    if(config.verbose) printf("Done\n");

}

int star2d::check_convec(double &p_cc,matrix &Rcc) {
    DEBUG_FUNCNAME;
    if(!core_convec) return 0;

    if(conv) {
        int j=0;
        for(int n=0;n<conv;n++) j+=map.gl.npts[n];
        if(z(j)<min_core_size) {
            if(config.verbose)
                ester_warn("Size(convective core) < min_core_size. Removing...");
            return 0;
        } else return conv;
    }
    // else
    int i=0;
    while(z(i)<1.05*min_core_size) i++;
    matrix schw;
    schw=-(map.gzz*(D,p)+map.gzt*(p,Dt))*((D,log(T))-eos.del_ad*(D,log(p)))
        -(map.gzt*(D,p)+map.gtt*(p,Dt))*((log(T),Dt)-eos.del_ad*(log(p),Dt));
    schw.setrow(0,zeros(1,nth));
    schw=schw/r/r;
    schw.setrow(0,zeros(1,nth));
    schw.setrow(0,-(D.row(0),schw)/D(0,0));
    if(map.leg.eval_00(schw.row(i),0)(0)>=0) return 0;

    if(config.verbose) printf("Found convective core\n");
    if(ndomains==1) {
        fprintf(stderr,"Warning: At least 2 domains are needed to deal with core convection\n");
    }

    while(schw(i,-1)<0) i++;
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
        Rcc(j)=map.gl.eval(r.col(j),zi,TT)(0);
        pcore(j)=(TT,PRES.col(j))(0);
    }

    p_cc=map.leg.eval_00(pcore,0)(0);
    return 1;

}


