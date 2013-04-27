#include"physics.h"

static bool init=false;

extern "C" {
	void nuc_cesam_init_();
	void nuc_cesam_init_abon_(double *X,double *Z,double *comp);
	void nuc_cesam_eps_(double *t,double *ro,double *comp,
			double *epsilon,double *et,double *ero,double *ex);
}

int nuc_cesam(const matrix &X,double Z,const matrix &T,const matrix &rho,
		nuc_struct &nuc) {
		
	if(!init) {
		nuc_cesam_init_();
		double XX=0.7,ZZ=0.02;
		matrix ab(10);
		nuc_cesam_init_abon_(&XX,&ZZ,ab.data());
		init=true;
	}
	nuc.eps=zeros(T.nrows(),T.ncols());
	nuc.pp=zeros(T.nrows(),T.ncols());
	nuc.cno=zeros(T.nrows(),T.ncols());
	nuc.dlneps_lnrho=zeros(T.nrows(),T.ncols());
	nuc.dlneps_lnT=zeros(T.nrows(),T.ncols());
	matrix epsilon(4),ex(10),ab(10);
	double et,ero,t,ro;
	for(int j=0;j<T.ncols();j++) {
		for(int i=0;i<T.nrows();i++) {
			double_map comp;
			comp=init_comp(X(i,j),Z);
			ab(0)=comp["H"];
			ab(1)=comp["He3"];
			ab(2)=comp["He4"];
			ab(3)=comp["C12"];
			ab(4)=comp["C13"];
			ab(5)=comp["N14"];
			ab(6)=comp["N15"];
			ab(7)=comp["O16"];
			ab(8)=comp["O17"];
			ab(9)=1.-comp["H"]-comp["He3"]-comp["He4"]-comp["C12"]-comp["C13"]-comp["N14"]
				-comp["N15"]-comp["O16"]-comp["O17"];
			t=T(i,j);ro=rho(i,j);
			nuc_cesam_eps_(&t,&ro,ab.data(),epsilon.data(),&et,&ero,ex.data());
			nuc.eps(i,j)=epsilon(0);
			nuc.pp(i,j)=epsilon(1);
			nuc.cno(i,j)=epsilon(2);
			if(epsilon(0)==0) {
				nuc.dlneps_lnrho(i,j)=0;
				nuc.dlneps_lnT(i,j)=0;
			} else {
				nuc.dlneps_lnrho(i,j)=ero/epsilon(0)*rho(i,j);
				nuc.dlneps_lnT(i,j)=et/epsilon(0)*T(i,j);
			}
		}
	}
	return 0;
}
