#include"star.h"
#include<string.h>
#include<stdlib.h>


star2d::star2d():r(map.r),z(map.z),D(map.D),th(map.th),Dt(map.Dt),Dt2(map.Dt2)
		,zex(map.ex.z),Dex(map.ex.D),rex(map.ex.r),
		nr(map.gl.N),nth(map.leg.npts),nex(map.ex.gl.N),ndomains(map.gl.ndomains) {

	config.newton_dmax=0.5;
	config.verbose=0;
}

star2d::~star2d() {
}

star2d::star2d(const star2d &A):r(map.r),z(map.z),D(map.D),th(map.th),Dt(map.Dt),Dt2(map.Dt2),zex(map.ex.z),Dex(map.ex.D),rex(map.ex.r),
	map(A.map),phi(A.phi),p(A.p),T(A.T),phiex(A.phiex),
	vr(A.vr),vt(A.vt),G(A.G),w(A.w),psi(A.psi),
	nr(map.gl.N),nth(map.leg.npts),nex(map.ex.gl.N),ndomains(map.gl.ndomains) {
		
	Omega=A.Omega;Omega_bk=A.Omega_bk;
	surff=A.surff;
	conv=A.conv;
	
	strcpy(opa.name,A.opa.name);
	strcpy(eos.name,A.eos.name);
	strcpy(nuc.name,A.nuc.name);
	strcpy(atm_name,A.atm_name);
	Tc=A.Tc;pc=A.pc;
	R=A.R;M=A.M;
	X=A.X;Z=A.Z;
	Xc=A.Xc;
	
	Ekman=A.Ekman;
	
	config=A.config;
	
	fill();
	
}

star2d &star2d::operator=(const star2d &A) {

	map=A.map;
	phi=A.phi;
	p=A.p;
	T=A.T;
	phiex=A.phiex;
	w=A.w;
	
	vr=A.vr;vt=A.vt;G=A.G;psi=A.psi;

	Omega=A.Omega;Omega_bk=A.Omega_bk;
	surff=A.surff;
	conv=A.conv;
	strcpy(opa.name,A.opa.name);
	strcpy(eos.name,A.eos.name);
	strcpy(nuc.name,A.nuc.name);
	strcpy(atm_name,A.atm_name);
	Tc=A.Tc;pc=A.pc;
	R=A.R;M=A.M;
	X=A.X;Z=A.Z;
	Xc=A.Xc;
	
	Ekman=A.Ekman;
	
	config=A.config;
	
	fill();
	
	return *this;

}



void star2d::write(const char *output_file,char mode) const{

	FILE *fp;
	const char *tag="star2d/";
	int ndom,i;
	
	if(mode=='b')
		fp=fopen(output_file,"wb");
	else
		fp=fopen(output_file,"wt");
	
	fwrite(tag,1,7,fp);
	fwrite(&mode,1,1,fp);
	
	ndom=ndomains;
	
	if(mode=='b') {
		fwrite(&ndom,sizeof(int),1,fp);
		fwrite(map.gl.npts,sizeof(int),ndom,fp);
		fwrite(map.gl.xif,sizeof(double),ndom+1,fp);
		fwrite(map.ex.gl.npts,sizeof(int),1,fp);
		fwrite(&map.leg.npts,sizeof(int),1,fp);
		fwrite(&M,sizeof(double),1,fp);
		fwrite(&R,sizeof(double),1,fp);
		fwrite(&X,sizeof(double),1,fp);
		fwrite(&Z,sizeof(double),1,fp);
		fwrite(&Xc,sizeof(double),1,fp);
		fwrite(&conv,sizeof(int),1,fp);
		fwrite(&surff,sizeof(double),1,fp);
		fwrite(&Omega,sizeof(double),1,fp);
		fwrite(&Omega_bk,sizeof(double),1,fp);
		fwrite(&Tc,sizeof(double),1,fp);
		fwrite(&pc,sizeof(double),1,fp);
		fwrite(opa.name,sizeof(char),strlen(opa.name)+1,fp);
		fwrite(eos.name,sizeof(char),strlen(eos.name)+1,fp);
		fwrite(nuc.name,sizeof(char),strlen(nuc.name)+1,fp);
		fwrite(atm_name,sizeof(char),strlen(atm_name)+1,fp);
	} else {
		fprintf(fp,"\n%d ",ndom);
		for(i=0;i<ndom;i++) fprintf(fp,"%d ",*(map.gl.npts+i));
		for(i=0;i<ndom+1;i++) fprintf(fp,"%16.16e ",*(map.gl.xif+i));
		fprintf(fp,"%d %d",*map.ex.gl.npts,map.leg.npts);
		fprintf(fp,"\n%16.16e %16.16e %16.16e %16.16e\n",M,R,X,Z);
		fprintf(fp,"%16.16e %d %16.16e\n",Xc,conv,surff);
		fprintf(fp,"%16.16e %16.16e\n",Omega,Omega_bk);
		fprintf(fp,"%16.16e %16.16e\n",Tc,pc);
		fprintf(fp,"%s\n",opa.name);
		fprintf(fp,"%s\n",eos.name);
		fprintf(fp,"%s\n",nuc.name);
		fprintf(fp,"%s\n",atm_name);	
	}
	
	map.R.write(fp,mode);
	phi.write(fp,mode);
	phiex.write(fp,mode);
	p.write(fp,mode);
	T.write(fp,mode);
	w.write(fp,mode);
	G.write(fp,mode);
	psi.write(fp,mode);
	if(mode=='b') {
		fwrite(&Ekman,sizeof(double),1,fp);
	} else {
		fprintf(fp,"\n%16.16e\n",Ekman);
	}
	fclose(fp);
		
}

int star2d::read(const char *input_file){

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
	map.R.read(ndomains,nth,fp,mode);
	map.remap();
	phi.read(nr,nth,fp,mode);
	phiex.read(nex,nth,fp,mode);
	p.read(nr,nth,fp,mode);
	T.read(nr,nth,fp,mode);
	w.read(nr,nth,fp,mode);
	G.read(nr,nth,fp,mode);
	psi.read(nr,nth,fp,mode);
	if(mode=='b') {
		fread(&Ekman,sizeof(double),1,fp);
	} else {
		fscanf(fp,"%le\n",&Ekman);
	}
	
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
				init(in1d,nt,next);
			} else {
				printf("Error reading input file: %s\n",input_file);
				return 0;
			}
		}
		map0=map;
	} else {
		if(!fp.open(default_params)) { 
			printf("Can't open default parameters file %s\n",default_params);
			return 0;
		}
		else {
			while(k=fp.get(arg,val)) {			
				if(i=check_arg(arg,val,&change_grid)) {
					printf("Sintax error in parameters file %s, line %d\n",default_params,k);
					if(i==2) {
						printf("Error: Argument to '%s' missing\n",arg);
						return 0;
					}
					if(i==1) {
						printf("Unknown parameter %s\n",arg);
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
			printf("Can't open parameters file %s\n",param_file);
			return 0;
		}
		else {
			while(k=fp.get(arg,val)) {			
				if(i=check_arg(arg,val,&change_grid)) {
					printf("Sintax error in parameters file %s, line %d\n",param_file,k);
					if(i==2) {
						printf("Error: Argument to '%s' missing\n",arg);
						return 0;
					}
					if(i==1) {
						printf("Unknown parameter %s\n",arg);
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
	
	if (*input_file) {
		if(change_grid) {
			map.gl.init();
			map.ex.gl.init();
			map.leg.init();
			leg_new=map.leg;
			map.leg=map0.leg;
			map.remap(leg_new);
			interp(map0);
		}
	} else {
		map.init();
		T=1-0.5*r*r;
		p=T;
		phi=-T;
		phiex=zeros(nex,nth);
		w=zeros(nr,nth);
		G=zeros(nr,nth);
		psi=zeros(nr,nth);
		Omega2=0;
	}
	fill();
	return 1;
	
}

void star2d::init(const star1d &A,int npts_th,int npts_ex) {

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

	map.gl=A.gl;
	map.leg.npts=npts_th;
	map.ex.gl.set_npts(npts_ex);
	map.init();

	thd=ones(1,npts_th);
	
	phi=A.phi*thd;
	p=A.p*thd;
	T=A.T*thd;
	w=zeros(T.nrows(),T.ncols());
	G=zeros(T.nrows(),T.ncols());
	psi=zeros(T.nrows(),T.ncols());
	
	surff=A.surff;
	conv=A.conv;
	strcpy(opa.name,A.opa.name);
	strcpy(eos.name,A.eos.name);
	strcpy(nuc.name,A.nuc.name);
	strcpy(atm_name,A.atm_name);
	Tc=A.Tc;pc=A.pc;
	R=A.R;M=A.M;
	X=A.X;Z=A.Z;
	Xc=A.Xc;
	
	phiex=A.phi(A.nr-1)/map.ex.r;	

	fill();

}


void star2d::interp(mapping map_old) {

	matrix Tr,Tex,Tt_00,Tt_01,Tt_10,Tt_11;
	
	map.interps(map_old,Tr,Tex,Tt_00,Tt_01,Tt_10,Tt_11);
	
	p=(Tr,p,Tt_00);
	T=(Tr,T,Tt_00);
	phi=(Tr,phi,Tt_00);
	phiex=(Tex,phiex,Tt_00);
	w=(Tr,w,Tt_00);
	G=(Tr,G,Tt_11);
	psi=(Tr,psi,Tt_00);
}

int star2d::check_arg(char *arg,char *val,int *change_grid) {

	int err=0,i;
	char *tok;

	if(!strcmp(arg,"ndomains")) {
		if(val==NULL) return 2;
		map.gl.set_ndomains(atoi(val));
		*change_grid=*change_grid|1;
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
		if(val==NULL) return 2;
		tok=strtok(val,",");
		i=0;
		while(tok!=NULL) {
			*(map.gl.xif+i)=atof(tok);
			tok=strtok(NULL,",");
			i++;
		}
		if(i==1) {
			double gamma=*map.gl.xif;
			*map.gl.xif=0;
			for(i=1;i<map.gl.ndomains;i++) 
				*(map.gl.xif+i)=1.-pow(1-(double) i/map.gl.ndomains,gamma);
			*(map.gl.xif+map.gl.ndomains)=1;
		}		
		*change_grid=*change_grid|1;
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
		if(val==NULL) return 2;
		conv=atoi(val);
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
	else err=1;

	return err;

}


