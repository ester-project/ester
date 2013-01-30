#include"star.h"
#include"version.h"
#include<string.h>
#include<stdlib.h>

star2d::star2d() :r(map.r),z(map.z),D(map.D),th(map.th),Dt(map.Dt),Dt2(map.Dt2)
		,zex(map.ex.z),Dex(map.ex.D),rex(map.ex.r),
		nr(map.gl.N),nth(map.leg.npts),nex(map.ex.gl.N),ndomains(map.gl.ndomains) {

	config.newton_dmax=0.5;
	config.verbose=0;
	*version='\0';

}

star2d::~star2d() {
}

star2d::star2d(const star2d &A) :r(map.r),z(map.z),D(map.D),th(map.th),Dt(map.Dt),Dt2(map.Dt2)
	,zex(map.ex.z),Dex(map.ex.D),rex(map.ex.r),
	nr(map.gl.N),nth(map.leg.npts),nex(map.ex.gl.N),ndomains(map.gl.ndomains){

	copy(A);

}

star2d &star2d::operator=(const star2d &A) {

	copy(A);
	
	return *this;

}

void star2d::copy(const star2d &A) {

	phi=A.phi;
	p=A.p;
	T=A.T;
	rho=A.rho;
	Xr=A.Xr;
	
	opa=A.opa;
	eos=A.eos;
	nuc=A.nuc;
	strcpy(atm_name,A.atm_name);
	config=A.config;
	units=A.units;
	
	R=A.R;M=A.M;
	Tc=A.Tc;pc=A.pc;rhoc=A.rhoc;
	X=A.X;Y=A.Y;Z=A.Z;
	surff=A.surff;
	conv=A.conv;
	Xc=A.Xc;
	
	map=A.map;
	
	ps=A.ps;Ts=A.Ts;
	m=A.m;pi_c=A.pi_c;Lambda=A.Lambda;
	
	phiex=A.phiex;
	w=A.w;
	vr=A.vr;vt=A.vt;G=A.G;
	
	Omega=A.Omega;Omega_bk=A.Omega_bk;
	Ekman=A.Ekman;
	
	core_convec=A.core_convec;
	min_core_size=A.min_core_size;
	
}

void star2d::write(const char *output_file,char mode) const {

	OUTFILE fp;

	char cur_version[32];	
	strncpy(cur_version,ESTER_VERSION,31);cur_version[31]='\0';
	
	fp.open(output_file,mode);
	write_tag(&fp,mode);
	if(mode=='b') {
		fp.write("ndomains",&ndomains);
		fp.write("npts",map.gl.npts,ndomains);
		fp.write("xif",map.gl.xif,ndomains+1);
		fp.write("nth",&map.leg.npts,1);
		fp.write("nex",map.ex.gl.npts,1);
		fp.write("M",&M);
		fp.write("R",&R);
		fp.write("X",&X);
		fp.write("Z",&Z);
		fp.write("Xc",&Xc);
		fp.write("conv",&conv);
		fp.write("surff",&surff);
		fp.write("Tc",&Tc);
		fp.write("pc",&pc);
		fp.write("opa.name",opa.name,strlen(opa.name)+1);
		fp.write("eos.name",eos.name,strlen(eos.name)+1);
		fp.write("nuc.name",nuc.name,strlen(nuc.name)+1);
		fp.write("atm_name",atm_name,strlen(atm_name)+1);
		fp.write("Omega",&Omega);
		fp.write("Omega_bk",&Omega_bk);
		fp.write("Ekman",&Ekman);
		fp.write("core_convec",&core_convec);
		fp.write("min_core_size",&min_core_size);
		fp.write("version",cur_version,strlen(cur_version)+1);
	} else {
		fp.write_fmt("ndomains","%d",&ndomains);
		fp.write_fmt("npts","%d",map.gl.npts,ndomains);
		fp.write_fmt("xif","%.16e",map.gl.xif,ndomains+1);
		fp.write_fmt("nth","%d",&map.leg.npts,1);
		fp.write_fmt("nex","%d",map.ex.gl.npts,1);
		fp.write_fmt("M","%.16e",&M);
		fp.write_fmt("R","%.16e",&R);
		fp.write_fmt("X","%.16e",&X);
		fp.write_fmt("Z","%.16e",&Z);
		fp.write_fmt("Xc","%.16e",&Xc);
		fp.write_fmt("conv","%d",&conv);
		fp.write_fmt("surff","%.16e",&surff);
		fp.write_fmt("Tc","%.16e",&Tc);
		fp.write_fmt("pc","%.16e",&pc);
		fp.write_fmt("opa.name","%s",&opa.name);
		fp.write_fmt("eos.name","%s",&eos.name);
		fp.write_fmt("nuc.name","%s",&nuc.name);
		fp.write_fmt("atm_name","%s",&atm_name);
		fp.write_fmt("Omega","%.16e",&Omega);
		fp.write_fmt("Omega_bk","%.16e",&Omega_bk);
		fp.write_fmt("Ekman","%.16e",&Ekman);
		fp.write_fmt("core_convec","%d",&core_convec);
		fp.write_fmt("min_core_size","%.16e",&min_core_size);
		fp.write_fmt("version","%s",&cur_version);
	}
	
	fp.write("phi",&phi);
	fp.write("p",&p);
	fp.write("T",&T);
	fp.write("phiex",&phiex);
	fp.write("map.R",&map.R);
	fp.write("w",&w);
	fp.write("G",&G);
	
	fp.close();

}

int star2d::read(const char *input_file){

	char tag[32],mode;
	int ndom;
	INFILE fp;
	
	if(fp.open(input_file,'b')) mode='b';
	else if(fp.open(input_file,'t')) mode='t';
	else return read_old(input_file);
	
	if(mode=='t') fp.read_fmt("tag","%s",tag);
	else {
		tag[0]='\0';
		if(fp.len("tag")<=16) fp.read("tag",tag);
	}
	tag[16]='\0';
	if(!check_tag(tag)) {
		fp.close();
		return 0;
	}
	
	if(mode=='b') {
		fp.read("ndomains",&ndom);
		map.gl.set_ndomains(ndom);
		fp.read("npts",map.gl.npts);
		fp.read("xif",map.gl.xif);
		if(!fp.read("nth",&map.leg.npts)) map.leg.npts=1;
		fp.read("nex",map.ex.gl.npts);
		fp.read("M",&M);
		fp.read("R",&R);
		fp.read("X",&X);
		fp.read("Z",&Z);
		fp.read("Xc",&Xc);
		fp.read("conv",&conv);
		fp.read("surff",&surff);
		fp.read("Tc",&Tc);
		fp.read("pc",&pc);
		fp.read("opa.name",opa.name);
		fp.read("eos.name",eos.name);
		fp.read("nuc.name",nuc.name);
		fp.read("atm_name",atm_name);
		if(!fp.read("Omega",&Omega)) Omega=0;
		if(!fp.read("Omega_bk",&Omega_bk)) Omega_bk=0;
		if(!fp.read("Ekman",&Ekman)) Ekman=0;
		if(!fp.read("core_convec",&core_convec)) core_convec=1;
		if(!fp.read("min_core_size",&min_core_size)) min_core_size=0.01;
		if(!fp.read("version",version)) strcpy(version,"0");
	} else {
		fp.read_fmt("ndomains","%d",&ndom);
		map.gl.set_ndomains(ndom);
		fp.read_fmt("npts","%d",map.gl.npts);
		fp.read_fmt("xif","%le",map.gl.xif);
		if(!fp.read_fmt("nth","%d",&map.leg.npts)) map.leg.npts=1;
		fp.read_fmt("nex","%d",map.ex.gl.npts);
		fp.read_fmt("M","%le",&M);
		fp.read_fmt("R","%le",&R);
		fp.read_fmt("X","%le",&X);
		fp.read_fmt("Z","%le",&Z);
		fp.read_fmt("Xc","%le",&Xc);
		fp.read_fmt("conv","%d",&conv);
		fp.read_fmt("surff","%le",&surff);
		fp.read_fmt("Tc","%le",&Tc);
		fp.read_fmt("pc","%le",&pc);
		fp.read_fmt("opa.name","%s",opa.name);
		fp.read_fmt("eos.name","%s",eos.name);
		fp.read_fmt("nuc.name","%s",nuc.name);
		fp.read_fmt("atm_name","%s",atm_name);
		if(!fp.read_fmt("Omega","%le",&Omega)) Omega=0;
		if(!fp.read_fmt("Omega_bk","%le",&Omega_bk)) Omega_bk=0;
		if(!fp.read_fmt("Ekman","%le",&Ekman)) Ekman=0;
		if(!fp.read_fmt("core_convec","%d",&core_convec)) core_convec=1;
		if(!fp.read_fmt("min_core_size","%le",&min_core_size)) min_core_size=0.01;
		if(!fp.read_fmt("version","%s",version)) strcpy(version,"0");
	}
		
	map.init();
	
	fp.read("phi",&phi);
	fp.read("p",&p);
	fp.read("T",&T);
	if(!fp.read("phiex",&phiex)) phiex=zeros(nex,nth);
	fp.read("map.R",&map.R);
	if(!fp.read("w",&w)) w=zeros(nr,nth);
	if(!fp.read("G",&G)) G=zeros(nr,nth);
	
	map.remap();
	fp.close();
	fill();
	return 1;
	
}

void star2d::write_tag(OUTFILE *fp,char mode) const {

	char tag[7]="star2d";
	
	if(mode=='b') fp->write("tag",tag,7);
	else fp->write_fmt("tag","%s",&tag);
		
}

int star2d::check_tag(const char *tag) const {

	if(strcmp(tag,"star2d")) return 0;
	return 1;

}

int star2d::read_old(const char *input_file){

	FILE *fp,*fp2;
	char tag[7],mode,*c;
	int ndom,i;
	
	if(!(fp=fopen(input_file,"rb"))) {
		return 0;
	}
	fread(tag,1,7,fp);
	tag[6]='\0';
	if(strcmp(tag,"star2d")) {
		return 0;
	}
	fread(&mode,1,1,fp);
	fclose(fp);

	if(mode=='b')
		fp=fopen(input_file,"rb");
	else
		fp=fopen(input_file,"rt");
	
	fread(tag,1,7,fp);
	fread(&mode,1,1,fp);
	
	if(mode=='b') {
		fread(&ndom,sizeof(int),1,fp);
		map.gl.set_ndomains(ndom);
		fread(map.gl.npts,sizeof(int),ndom,fp);
		fread(map.gl.xif,sizeof(double),ndom+1,fp);
		map.gl.init();
		fread(map.ex.gl.npts,sizeof(int),1,fp);
		map.ex.gl.init();
		fread(&map.leg.npts,sizeof(int),1,fp);
		map.leg.init();
		fread(&M,sizeof(double),1,fp);
		fread(&R,sizeof(double),1,fp);
		fread(&X,sizeof(double),1,fp);
		fread(&Z,sizeof(double),1,fp);
		fread(&Xc,sizeof(double),1,fp);
		fread(&conv,sizeof(int),1,fp);
		fread(&surff,sizeof(double),1,fp);
		fread(&Omega,sizeof(double),1,fp);
		fread(&Omega_bk,sizeof(double),1,fp);
		fread(&Tc,sizeof(double),1,fp);
		fread(&pc,sizeof(double),1,fp);
		c=opa.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=eos.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=nuc.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=atm_name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
	} else {
		fscanf(fp,"\n%d ",&ndom);
		map.gl.set_ndomains(ndom);
		for(i=0;i<ndom;i++) fscanf(fp,"%d ",(map.gl.npts+i));
		for(i=0;i<ndom+1;i++) fscanf(fp,"%le ",(map.gl.xif+i));
		map.gl.init();
		fscanf(fp,"%d %d",map.ex.gl.npts,&map.leg.npts);
		map.ex.gl.init();
		map.leg.init();
		fscanf(fp,"\n%le %le %le %le\n",&M,&R,&X,&Z);
		fscanf(fp,"%le %d %le\n",&Xc,&conv,&surff);
		fscanf(fp,"%le %le\n",&Omega,&Omega_bk);
		fscanf(fp,"%le %le\n",&Tc,&pc);
		fscanf(fp,"%s\n",opa.name);
		fscanf(fp,"%s\n",eos.name);
		fscanf(fp,"%s\n",nuc.name);
		fscanf(fp,"%s\n",atm_name);
	}
	map.init();
	map.R.read(ndomains,nth,fp,mode);
	map.remap();
	phi.read(nr,nth,fp,mode);
	phiex.read(nex,nth,fp,mode);
	p.read(nr,nth,fp,mode);
	T.read(nr,nth,fp,mode);
	w.read(nr,nth,fp,mode);
	G.read(nr,nth,fp,mode);
	matrix psi;
	psi.read(nr,nth,fp,mode);
	if(mode=='b') {
		fread(&Ekman,sizeof(double),1,fp);
	} else {
		fscanf(fp,"%le\n",&Ekman);
	}
	core_convec=1;
	min_core_size=0.01;
	strcpy(version,"0");
	fclose(fp);
	fill();
	return 1;
		
}

int star2d::init(const char *input_file,const char *param_file,int argc,char *argv[]) {

	cmdline_parser cmd;
	file_parser fp;
	char *arg,*val,default_params[256];
	mapping map0;
	int i,j,k,change_grid=0,nt=-1,next=-1;
	star1d in1d;
	diff_leg leg_new;
	matrix Tr,m0;

	sprintf(default_params,"%s/config/2d_default.par",ESTER_ROOT);

	if(*input_file) {
		if(!read(input_file)) {
			if(in1d.read(input_file)) {
				if(*param_file) {
					if(fp.open(param_file)) { 
						while(k=fp.get(arg,val)) {
							if(!strcmp(arg,"nth")&&val) nt=atoi(val);		
							if(!strcmp(arg,"nex")&&val) next=atoi(val);
						}
						fp.close();
					}
				}
				cmd.open(argc,argv);
				while(cmd.get(arg,val)) {
					if(!strcmp(arg,"nth")&&val) nt=atoi(val);		
					if(!strcmp(arg,"nex")&&val) next=atoi(val);
				}
				cmd.close();
				init1d(in1d,nt,next);
			} else {
				fprintf(stderr,"Error reading input file: %s\n",input_file);
				return 0;
			}
		}
		map0=map;
	} else {
		if(!fp.open(default_params)) { 
			fprintf(stderr,"Can't open default parameters file %s\n",default_params);
			return 0;
		}
		else {
			while(k=fp.get(arg,val)) {			
				if(i=check_arg(arg,val,&change_grid)) {
					fprintf(stderr,"Sintax error in parameters file %s, line %d\n",default_params,k);
					if(i==2) {
						fprintf(stderr,"Error: Argument to '%s' missing\n",arg);
						return 0;
					}
					if(i==1) {
						fprintf(stderr,"Unknown parameter %s\n",arg);
						return 0;
					}
				}
			}
			fp.close();
		}
		change_grid=0;
	}
	
	if(*param_file) {
		if(!fp.open(param_file)) { 
			fprintf(stderr,"Can't open parameters file %s\n",param_file);
			return 0;
		}
		else {
			while(k=fp.get(arg,val)) {			
				if(i=check_arg(arg,val,&change_grid)) {
					fprintf(stderr,"Sintax error in parameters file %s, line %d\n",param_file,k);
					if(i==2) {
						fprintf(stderr,"Error: Argument to '%s' missing\n",arg);
						return 0;
					}
					if(i==1) {
						fprintf(stderr,"Unknown parameter %s\n",arg);
						return 0;
					}
				}
			}
			fp.close();
		}
	}

	cmd.open(argc,argv);
	while(int err_code=cmd.get(arg,val)) {
		if(err_code==-1) exit(1);
		err_code=check_arg(arg,val,&change_grid);
		if(err_code==2) {
			fprintf(stderr,"Error: Argument to '%s' missing\n",arg);
			return 0;
		}
		if(err_code==1) {
			fprintf(stderr,"Unknown parameter '%s'\n",arg);
			return 0;
		}
		cmd.ack(arg,val);
	}
	cmd.close();
	
	if((change_grid&1)&&!(change_grid&2)) {
		fprintf(stderr,"Must specify number of points per domain (npts)\n");
		exit(1);
	}
	
	if (*input_file) {
		if(change_grid) {
			mapping map_new;
			map_new=map;
			map=map0;
			remap(map_new.ndomains,map_new.gl.npts,map_new.nth,map_new.nex);
		}
	} else {
		for(i=0;i<=ndomains;i++) map.gl.xif[i]=i*1./ndomains;
		map.init();
		T=1-0.5*r*r;
		p=T;
		phi=-T;
		phiex=zeros(nex,nth);
		w=zeros(nr,nth);
		G=zeros(nr,nth);
		conv=0;
	}
	fill();
	return 1;
	
}

void star2d::init1d(const star1d &A,int npts_th,int npts_ex) {

	matrix thd;
	char *arg,*val,default_params[256];
	int k;
	file_parser fp;

	sprintf(default_params,"%s/config/2d_default.par",ESTER_ROOT);

	Omega_bk=0;	
	Omega=0;
	
	if(fp.open(default_params)) { 
		while(k=fp.get(arg,val)) {
			if(!strcmp(arg,"nth")&&val&&npts_th==-1) npts_th=atoi(val);		
			if(!strcmp(arg,"nex")&&val&&npts_ex==-1) npts_ex=atoi(val);
			if(!strcmp(arg,"Omega_bk")&&val) Omega_bk=atof(val);
			if(!strcmp(arg,"Ekman")&&val) Ekman=atof(val);
		}
		fp.close();
	}

	if(npts_th==-1) npts_th=8;
	if(npts_ex==-1) npts_ex=8;

	map.gl=A.map.gl;
	map.leg.npts=npts_th;
	map.ex.gl.set_npts(npts_ex);
	map.init();

	thd=ones(1,npts_th);
	
	phi=A.phi*thd;
	p=A.p*thd;
	T=A.T*thd;
	w=zeros(T.nrows(),T.ncols());
	G=zeros(T.nrows(),T.ncols());
	
	surff=A.surff;
	conv=A.conv;
	strcpy(opa.name,A.opa.name);
	strcpy(eos.name,A.eos.name);
	strcpy(nuc.name,A.nuc.name);
	strcpy(atm_name,A.atm_name);
	core_convec=A.core_convec;
	min_core_size=A.min_core_size;
	Tc=A.Tc;pc=A.pc;
	R=A.R;M=A.M;
	X=A.X;Z=A.Z;
	Xc=A.Xc;
	
	phiex=A.phi(A.nr-1)/map.ex.r;	

	fill();

}

void star2d::interp(mapping_redist *red) {

	p=red->interp(p);
	phi=red->interp(phi);
	T=red->interp(T);
	w=red->interp(w);
	G=red->interp(G,11);
	phiex=red->interp_ex(phiex);
	fill();

}

int star2d::check_arg(char *arg,char *val,int *change_grid) {

	int err=0,i;
	char *tok;

	if(!strcmp(arg,"ndomains")) {
		if(val==NULL) return 2;
		map.gl.set_ndomains(atoi(val));
		*change_grid=*change_grid|1;
		if(*change_grid&2) {
			fprintf(stderr,"ndomains must be specified before npts\n");
			exit(1);
		}
	}
	else if(!strcmp(arg,"npts")) {
		if(val==NULL) return 2;
		tok=strtok(val,",");
		i=0;
		while(tok!=NULL) {
			*(map.gl.npts+i)=atoi(tok);
			tok=strtok(NULL,",");
			i++;
		}	
		if(i==1) {
			for(i=1;i<map.gl.ndomains;i++) {
				*(map.gl.npts+i)=*map.gl.npts;
			}
		}
		*change_grid=*change_grid|2;
	}
	else if(!strcmp(arg,"xif")) {
		fprintf(stderr,"Warning: Parameter xif is now automatically handled by the code and cannot be modified by the user\n"); 
		return 1;
	}
	else if(!strcmp(arg,"nth")) {
		if(val==NULL) return 2;
		map.leg.npts=atoi(val);
		*change_grid=*change_grid|8;
	}
	else if(!strcmp(arg,"nex")) {
		if(val==NULL) return 2;
		map.ex.gl.set_npts(atoi(val));
		*change_grid=*change_grid|4;
	}
	else if(!strcmp(arg,"M")) {
		if(val==NULL) return 2;
		M=atof(val)*M_SUN;
	}
	else if(!strcmp(arg,"X")) {
		if(val==NULL) return 2;
		X=atof(val);
	}
	else if(!strcmp(arg,"Z")) {
		if(val==NULL) return 2;
		Z=atof(val);
	}
	else if(!strcmp(arg,"Xc")) {
		if(val==NULL) return 2;
		Xc=atof(val);
	}
	else if(!strcmp(arg,"conv")) {
		fprintf(stderr,"Param. conv is no longer modifiable. Disable core convection with core_convec 0.\n");
		return 1;
	}
	else if(!strcmp(arg,"surff")) {
		if(val==NULL) return 2;
		surff=atof(val);
	}
	else if(!strcmp(arg,"Omega_bk")) {
		if(val==NULL) return 2;
		Omega_bk=atof(val);
	}
	else if(!strcmp(arg,"Tc")) {
		if(val==NULL) return 2;
		Tc=atof(val);
	}
	else if(!strcmp(arg,"pc")) {
		if(val==NULL) return 2;
		pc=atof(val);
	}
	else if(!strcmp(arg,"opa")) {
		if(val==NULL) return 2;
		strcpy(opa.name,val);
	}
	else if(!strcmp(arg,"eos")) {
		if(val==NULL) return 2;
		strcpy(eos.name,val);
	}
	else if(!strcmp(arg,"nuc")) {
		if(val==NULL) return 2;
		strcpy(nuc.name,val);
	}
	else if(!strcmp(arg,"atm")) {
		if(val==NULL) return 2;
		strcpy(atm_name,val);
	}
	else if(!strcmp(arg,"Ekman")) {
		if(val==NULL) return 2;
		Ekman=atof(val);
	}
	else if(!strcmp(arg,"core_convec")) {
		if(val==NULL) return 2;
		core_convec=atoi(val);
	}
	else if(!strcmp(arg,"min_core_size")) {
		if(val==NULL) return 2;
		min_core_size=atof(val);
	}
	else err=1;

	return err;

}

void star2d::dump_info() {

	printf("\n2d ESTER model file  (Version %s)\n\n",version);
	
	printf("General parameters:\n\n");
	printf("\tMass = %.5f Msun (%e g)\n",M/M_SUN,M);
	printf("\tRadius (p) = %.5f Rsun (%e cm)\n",R/R_SUN,R);
	double re=map.leg.eval_00(r.row(-1),PI/2)(0);
	printf("\tRadius (e) = %.5f Rsun (%e cm)  (flat.=%.3f)\n",R/R_SUN*re,R*re,1.-1./re);
	printf("\tLuminosity = %.4f Lsun (%e erg/s)\n",luminosity()/L_SUN,luminosity());
	printf("\tTeff (p) = %.2f\n",map.leg.eval_00(Teff(),0)(0));
	printf("\tTeff (e) = %.2f\n",map.leg.eval_00(Teff(),PI/2)(0));
	printf("\tlog(geff) (p) = %.4f\n",log10(map.leg.eval_00(gsup(),0)(0)));
	printf("\tlog(geff) (e) = %.4f\n",log10(map.leg.eval_00(gsup(),PI/2)(0)));
	printf("\tX=%.4f   Y=%.4f   Z=%.4f\n",X,Y,Z);
	printf("\n");
	
	printf("Rotation:\n\n");
	double we=map.leg.eval_00(w.row(-1),PI/2)(0);
	double ve=we*re*units.Omega*units.r/1e5;
	printf("\tv_eq = %.3f km/s\n",ve);
	printf("\tOmega (e) = %e rad/s (%.2f%%)\n",we*units.Omega,we/Omegac*100);
	double wp=map.leg.eval_00(w.row(-1),0)(0);
	printf("\tOmega (p) = %e rad/s\n",wp*units.Omega);
	printf("\tOmega (c) = %e rad/s\n",w(0,0)*units.Omega);
	printf("\tPeriod (e) = %.5f days\n",2*PI/we/units.Omega/3600./24.);
	printf("\tPeriod (p) = %.5f days\n",2*PI/wp/units.Omega/3600./24.);
	printf("\tPeriod (c) = %.5f days\n",2*PI/w(0,0)/units.Omega/3600./24.);
	printf("\n");
	
	if(conv==0) printf("No convective core\n\n");
	else {
		printf("Convective core:\n\n");
		int jcc=0;
		for(int n=0;n<conv;n++) jcc+=map.gl.npts[n];
		jcc--;
		double mcc=2*PI*(map.gl.I.block(0,0,0,jcc),
			(rho*r*r*map.rz).block(0,jcc,0,-1),map.leg.I_00)(0)*units.rho*units.r*units.r*units.r;
		printf("\tMass_core = %.5f Msun (%e g)\n",mcc/M_SUN,mcc);
		double rcc_p=map.leg.eval_00(r.row(jcc),0)(0)*units.r;
		printf("\tRadius_core (p) = %.5f Rsun (%e cm)\n",rcc_p/R_SUN,rcc_p);
		double rcc_e=map.leg.eval_00(r.row(jcc),PI/2)(0)*units.r;
		printf("\tRadius_core (e) = %.5f Rsun (%e cm)  (flat.=%.3f)\n",rcc_e/R_SUN,rcc_e,1.-rcc_p/rcc_e);
		printf("\tX_core/X_env = %.4f\n",Xc);
		printf("\n");
	}
	printf("Central values:\n\n");
	printf("\tTemperature = %e K\n",Tc);
	printf("\tDensity = %e g/cm3\n",rhoc);
	printf("\tPressure = %e dyn/cm2\n",pc);
	printf("\n");
	
	printf("Grid parameters:\n\n");
	printf("\t # of domains = %d\n",ndomains);
	printf("\t # of domains in convective core = %d\n",conv);
	printf("\t nr = %d    (",nr);
	for(int n=0;n<ndomains;n++) {
		printf("%d",map.gl.npts[n]);
		if(n<ndomains-1) printf(",");
	}
	printf(")\n");
	printf("\t nth = %d\n",nth);
	printf("\t nex = %d\n",nex);
	printf("\n");
	
	printf("Additional parameters:\n\n");
	printf("\tOpacity = %s\n",opa.name);
	printf("\tEquation of state = %s\n",eos.name);
	printf("\tNuclear reactions = %s\n",nuc.name);
	printf("\tAtmosphere = %s\n",atm_name);
	printf("\tsurff = %e\n",surff);
	printf("\tcore_convec = %d\n",core_convec);
	printf("\tmin_core_size = %e\n",min_core_size);
	printf("\n");
	
	printf("Tests:\n\n");
	printf("\tVirial test = %e\n",virial());
	printf("\tEnergy test = %e\n",energy_test());
	printf("\n");
	
}




