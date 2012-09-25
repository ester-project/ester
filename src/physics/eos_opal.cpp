#include"physics.h"
#include"constants.h"

static double Z_table=-99;

extern"C" {
	void zfs_interp_eos5_(double *z);
	void eos5_xtrin_(double *x,double *ztab,double *t6,double *p,double *r);
	extern struct{
		double esact,eos[10];
	} eeos_;
	extern struct{
		int itime;
	} lreadco_;
}

int eos_opal_init(double Z) {
    
    double X;
    FILE *fp;
    char filename[512];
    
    //fprintf(stderr,"Initializing OPAL EOS table Z=%f\n",Z);
    sprintf(filename,"%s/tables/opal/eos_tables/EOS5_data",ESTER_ROOT);
    if(fp=fopen(filename,"rt")) {
    	fscanf(fp," X=%lf Z= %lf",&X,&Z_table);
    	fclose(fp);
    }
    if(Z_table==Z) return 0;
    zfs_interp_eos5_(&Z);
    Z_table=Z;
    lreadco_.itime=0;
    
    return 0;
}

int eos_opal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos) {
    
    matrix t6,p_mb;
    int i,N,error=0;
    
    t6=T*1e-6;
    p_mb=p*1e-12;
    
    if(Z!=Z_table) error=eos_opal_init(Z);
    
    if(error) {
    	printf("Can't initialize OPAL EOS table\n");
    	return error;
    }
    
    rho.dim(T.nrows(),T.ncols());
    eos.G1.dim(T.nrows(),T.ncols());
    eos.del_ad.dim(T.nrows(),T.ncols());
    eos.G3_1.dim(T.nrows(),T.ncols());
    eos.d.dim(T.nrows(),T.ncols());
    eos.cp.dim(T.nrows(),T.ncols());
    eos.cv.dim(T.nrows(),T.ncols());
    eos.chi_rho.dim(T.nrows(),T.ncols());
    eos.chi_T.dim(T.nrows(),T.ncols());
        
    eos.prad=A_RAD/3*pow(T,4);
    
    N=T.nrows()*T.ncols();

    for(i=0;i<N;i++) {
    	eos5_xtrin_(X.data()+i,&Z,&t6(i),&p_mb(i),&rho(i));
    	if(rho(i)==-9e99) {
   			printf("Values outside OPAL eos table\n");
    	}
		eos.G1(i)=*(eeos_.eos+7);
		eos.del_ad(i)=1/(*(eeos_.eos+8));
    	eos.G3_1(i)=*(eeos_.eos+7)/(*(eeos_.eos+8));
    	eos.d(i)=(*(eeos_.eos+6))/(*(eeos_.eos+5));
    	eos.cp(i)=1e6*(*(eeos_.eos+7))*(*(eeos_.eos+4))/(*(eeos_.eos+5));
    	eos.cv(i)=1e6*(*(eeos_.eos+4));
    	eos.chi_rho(i)=*(eeos_.eos+5);
    	eos.chi_T(i)=*(eeos_.eos+6);
    }

    return 0;
	
}

