#include"star.h"
#include<stdlib.h>

void star2d::remap(int ndom_in,int *npts_in,int nth_in,int nex_in) {

	mapping_redist red(map);
	
	red.set_ndomains(ndom_in);
	red.set_npts(npts_in);
	red.set_nth(nth_in);
	red.set_nex(nex_in);
	
	if(ndom_in!=ndomains) {
		matrix pif;
		int conv_new;
		pif=check_boundaries(ndom_in,conv_new);
		matrix R(ndom_in,nth);
		if(ndom_in>1) R.setblock(0,-2,0,-1,find_boundaries(pif.block(0,-2,0,0))); //faster
		//R=find_boundaries(pif);
		red.set_R(R);
		if(conv) red.set_fixed(conv_new-1,conv-1);
		conv=conv_new;
	}
	
	map=red.get_map();
	interp(&red);

}

matrix star2d::check_boundaries(int ndom,int &conv_new,double p_cc) const {

	int j;
	double prad0,prad1,pconv0,pconv1,dprad,dpconv;
	double p_s;
	
	conv_new=conv;
	p_s=map.leg.eval_00(p.row(-1),0)(0);
	j=0;
	if(p_cc==0) {
		for(int n=0;n<conv;n++) j+=map.gl.npts[n];
		p_cc=map.leg.eval_00(p.row(j),0)(0);
	}
	
	if(conv) {
		int fin=0;
		while(!fin) {
			dprad=(log(p_cc)-log(p_s))/(ndom-conv_new);
			dpconv=-log(p_cc)/conv_new;
			fin=1;
			if(dpconv>1.5*dprad&&conv_new<ndom-1) {
				conv_new++;fin=0;
			} else if(dpconv<0.6*dprad&&conv_new>1) {
				conv_new--;fin=0;
			}
		}
	}
	
	dprad=(log(p_cc)-log(p_s))/(ndom-conv_new);
	dpconv=-log(p_cc)/conv_new;
	matrix pif(ndom,1);
	for(int n=0;n<conv_new;n++) {
		pif(n)=exp(-(n+1)*dpconv);
	}
	for(int n=conv_new;n<ndom;n++) {
		pif(n)=exp(log(p_cc)-(n-conv_new+1)*dprad);
	}
	
	return pif;
}

matrix star2d::find_boundaries(matrix pif) const {

	matrix R(pif.nrows(),nth);

	for(int j=0;j<nth;j++) {
		matrix zj,dzj,pj,dpj,T;
		zj=0.5*ones(pif.nrows(),1);
		int fin=0,nit=0;
		while(fin<3) {
			pj=log(map.gl.eval(p.col(j),zj,T));
			dpj=(T,(D,log(p).col(j)));
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

	double pcc;
	matrix Rcc,pif;
	int conv_new;
	mapping_redist *red;

	if(check_convec(pcc,Rcc)!=conv) {
		matrix R(ndomains,nth);
		red=new mapping_redist(map);
		if(conv) {
			conv=0;
			pif=check_boundaries(ndomains,conv_new);
			R.setblock(0,-2,0,-1,find_boundaries(pif.block(0,-2,0,0))); //faster
			//R=find_boundaries(pif);
			red->set_R(R);
		} else {
			conv=1;
			pif=check_boundaries(ndomains,conv_new,pcc);
			conv=conv_new;
			R.setblock(0,-2,0,-1,find_boundaries(pif.block(0,-2,0,0))); //faster
			//R=find_boundaries(pif);
			R.setrow(conv-1,Rcc);
			red->set_R(R);
		}
	} else {
		pif=check_boundaries(ndomains,conv_new);
		if(conv_new==conv) return;
		
		red=new mapping_redist(map);
		matrix R(ndomains,nth);
		R.setblock(0,-2,0,-1,find_boundaries(pif.block(0,-2,0,0))); //faster
		//R=find_boundaries(pif);
		red->set_R(R);
		red->set_fixed(conv_new-1,conv-1);
		conv=conv_new;
	}
	if(config.verbose) {printf("Remapping...");fflush(stdout);}
	map=red->get_map();
	interp(red);
	if(config.verbose) printf("Done\n");

}

int star2d::check_convec(double &p_cc,matrix &Rcc) {

	if(!core_convec) return 0;

	if(conv) {
		int j=0;
		for(int n=0;n<conv;n++) j+=map.gl.npts[n];
		if(z(j)<0.9*min_core_size) {
			if(config.verbose) printf("Size(convective core) < min_core_size. Removing...\n");
			return 0;
		} else return conv;
	}
	// else
	int i=0;
	while(z(i)<min_core_size) i++;
	matrix n2;
	n2=N2();
	if(map.leg.eval_00(n2.row(i),0)(0)>=0) return 0;
	
	if(config.verbose) printf("Found convective core\n");
	if(ndomains==1) {
		fprintf(stderr,"Warning: At least 2 domains are needed to deal with core convection\n");
	}
	
	while(n2(i,-1)<0) i++;
	double zi=z(i);
	matrix dn2;
	dn2=(D,n2);
	Rcc=zeros(1,nth);
	matrix pcore(Rcc);
	for(int j=nth-1;j>=0;j--) {
		int fin=0;
		double n2i,dn2i,dzi;
		matrix T;
		while(fin<3) {
			n2i=map.gl.eval(n2.col(j),zi,T)(0);
			dn2i=(T,dn2.col(j))(0);
			dzi=-n2i/dn2i;
			if(fabs(dzi<1e-9)) fin++;
			zi+=dzi;
		}
		Rcc(j)=map.gl.eval(r.col(j),zi,T)(0);
		pcore(j)=(T,p.col(j))(0);
	}
	
	p_cc=map.leg.eval_00(pcore,0)(0);
	return 1;

}


