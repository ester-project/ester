#include"star.h"
#include<string.h>
#include<stdlib.h>

star1d::star1d():r(gl.x),D(gl.D),nr(gl.N),ndomains(gl.ndomains) {}

star1d::~star1d() {
}

star1d::star1d(const star1d &A) : star(A)
	,r(gl.x),D(gl.D),nr(gl.N),ndomains(gl.ndomains) {

	copy(A);
	
}

star1d &star1d::operator=(const star1d &A) {

	star::operator=(A);
	copy(A);
	
	return *this;

}

void star1d::copy(const star1d &A) {

	gl=A.gl;

	surff=A.surff;
	conv=A.conv;

	strcpy(atm_name,A.atm_name);
	Xc=A.Xc;

	fill();

}

void star1d::write(const char *output_file,char mode) const {

	OUTFILE fp;
	char tag[7]="star1d";
	
	fp.open(output_file,mode);
	if(mode=='b') {
		fp.write("tag",tag,7);
		fp.write("ndomains",&ndomains);
		fp.write("npts",gl.npts,ndomains);
		fp.write("xif",gl.xif,ndomains+1);
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
	} else {
		fp.write_fmt("tag","%s",&tag);
		fp.write_fmt("ndomains","%d",&ndomains);
		fp.write_fmt("npts","%d",gl.npts,ndomains);
		fp.write_fmt("xif","%.16e",gl.xif,ndomains+1);
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
	}
	
	fp.write("phi",&phi);
	fp.write("p",&p);
	fp.write("T",&T);
	
	fp.close();

}

int star1d::read(const char *input_file){

	char tag[1024],mode;
	int ndom;
	INFILE fp;
	
	if(fp.open(input_file,'b')) mode='b';
	else if(fp.open(input_file,'t')) mode='t';
	else return read_old(input_file);
	
	if(mode=='t') fp.read_fmt("tag","%s",tag);
	else {
		if(fp.len("tag")>16) tag[0]='\0';
		else fp.read("tag",tag);
	}
	tag[16]='\0';
	if(strcmp(tag,"star1d")) {
		fp.close();
		return 0;
	}
	
	if(mode=='b') {
		fp.read("ndomains",&ndom);
		gl.set_ndomains(ndom);
		fp.read("npts",gl.npts);
		fp.read("xif",gl.xif);
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
	} else {
		fp.read_fmt("ndomains","%d",&ndom);
		gl.set_ndomains(ndom);
		fp.read_fmt("npts","%d",gl.npts);
		fp.read_fmt("xif","%le",gl.xif);
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
	}
		
	gl.init();
	
	fp.read("phi",&phi);
	fp.read("p",&p);
	fp.read("T",&T);
	
	fp.close();
	fill();
	return 1;
	
}

int star1d::read_old(const char *input_file){

	FILE *fp;
	char tag[7],mode,*c;
	int ndom,i;
	
	if(!(fp=fopen(input_file,"rb"))) {
		return 0;
	}
	fread(tag,1,7,fp);
	tag[6]='\0';
	if(strcmp(tag,"star1d")) {
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
		gl.set_ndomains(ndom);
		fread(gl.npts,sizeof(int),ndom,fp);
		fread(gl.xif,sizeof(double),ndom+1,fp);
		fread(&M,sizeof(double),1,fp);
		fread(&R,sizeof(double),1,fp);
		fread(&X,sizeof(double),1,fp);
		fread(&Z,sizeof(double),1,fp);
		fread(&Xc,sizeof(double),1,fp);
		fread(&conv,sizeof(int),1,fp);
		fread(&surff,sizeof(double),1,fp);
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
		gl.set_ndomains(ndom);
		for(i=0;i<ndom;i++) fscanf(fp,"%d ",(gl.npts+i));
		for(i=0;i<ndom+1;i++) fscanf(fp,"%le ",(gl.xif+i));
		fscanf(fp,"\n%le %le %le %le\n",&M,&R,&X,&Z);
		fscanf(fp,"%le %d %le\n",&Xc,&conv,&surff);		
		fscanf(fp,"%le %le\n",&Tc,&pc);
		fscanf(fp,"%s\n",opa.name);
		fscanf(fp,"%s\n",eos.name);
		fscanf(fp,"%s\n",nuc.name);
		fscanf(fp,"%s\n",atm_name);
	}
	gl.init();
	
	phi.read(nr,1,fp,mode);
	p.read(nr,1,fp,mode);
	T.read(nr,1,fp,mode);
	
	fclose(fp);
	fill();
	return 1;
		
}

int star1d::init(const char *input_file,const char *param_file,int argc,char *argv[]) {

	cmdline_parser cmd;
	file_parser fp;
	char *arg,*val,default_params[256];
	diff_gl gl0;
	int i,k,change_grid=0;
	matrix Tr;

	sprintf(default_params,"%s/config/1d_default.par",ESTER_ROOT);

	if(*input_file) {
		if(!read(input_file)) {
			printf("Error reading input file: %s\n",input_file);
			return 0;
		}
		gl0=gl;
	} else {
		if(!fp.open(default_params)) { 
			printf("Can't open default parameters file %s\n",default_params);
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
			gl.init();
			T=gl0.eval(T,gl.x,Tr);
			p=(Tr,p);
			phi=(Tr,phi);
		}
	} else {
		gl.init();
		T=1-0.5*r*r;
		p=T;
		phi=-T;
	}
	
	fill();
	
	return 1;
	
}

int star1d::check_arg(char *arg,char *val,int *change_grid) {

	int err=0,i;
	char *tok;

	if(!strcmp(arg,"ndomains")) {
		if(val==NULL) return 2;
		gl.set_ndomains(atoi(val));
		*change_grid=1;
	}
	else if(!strcmp(arg,"npts")) {
		if(val==NULL) return 2;
		tok=strtok(val,",");
		i=0;
		while(tok!=NULL) {
			*(gl.npts+i)=atoi(tok);
			tok=strtok(NULL,",");
			i++;
		}
		if(i==1) {
			for(i=1;i<gl.ndomains;i++) {
				*(gl.npts+i)=*gl.npts;
			}
		}	
		*change_grid=1;
	}
	else if(!strcmp(arg,"xif")) {
		if(val==NULL) return 2;
		tok=strtok(val,",");
		i=0;
		while(tok!=NULL) {
			*(gl.xif+i)=atof(tok);
			tok=strtok(NULL,",");
			i++;
		}
		if(i==1) {
			double gamma=*gl.xif;
			*gl.xif=0;
			for(i=1;i<gl.ndomains;i++) 
				*(gl.xif+i)=1.-pow(1-(double) i/gl.ndomains,gamma);
			*(gl.xif+gl.ndomains)=1;
		}		
		*change_grid=1;
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
	else err=1;

	return err;

}


