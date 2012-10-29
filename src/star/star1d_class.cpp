#include"star.h"
#include<string.h>
#include<stdlib.h>

star1d::star1d():r(gl.x),D(gl.D),nr(gl.N),ndomains(gl.ndomains) {
	config.newton_dmax=0.5;
	config.verbose=0;
}

star1d::~star1d() {
}

star1d::star1d(const star1d &A):r(gl.x),D(gl.D),gl(A.gl),phi(A.phi),
		p(A.p),T(A.T),nr(gl.N),ndomains(gl.ndomains) {
		
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
	config=A.config;
	fill();
	
}

star1d &star1d::operator=(const star1d &A) {

	gl=A.gl;
	phi=A.phi;
	p=A.p;
	T=A.T;

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
	config=A.config;
	fill();
	
	return *this;

}

void star1d::write(const char *output_file,char mode) const{

	FILE *fp;
	const char *tag="star1d/";
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
		fwrite(gl.npts,sizeof(int),ndom,fp);
		fwrite(gl.xif,sizeof(double),ndom+1,fp);
		fwrite(&M,sizeof(double),1,fp);
		fwrite(&R,sizeof(double),1,fp);
		fwrite(&X,sizeof(double),1,fp);
		fwrite(&Z,sizeof(double),1,fp);
		fwrite(&Xc,sizeof(double),1,fp);
		fwrite(&conv,sizeof(int),1,fp);
		fwrite(&surff,sizeof(double),1,fp);
		fwrite(&Tc,sizeof(double),1,fp);
		fwrite(&pc,sizeof(double),1,fp);
		fwrite(opa.name,sizeof(char),strlen(opa.name)+1,fp);
		fwrite(eos.name,sizeof(char),strlen(eos.name)+1,fp);
		fwrite(nuc.name,sizeof(char),strlen(nuc.name)+1,fp);
		fwrite(atm_name,sizeof(char),strlen(atm_name)+1,fp);
	} else {
		fprintf(fp,"\n%d ",ndom);
		for(i=0;i<ndom;i++) fprintf(fp,"%d ",*(gl.npts+i));
		for(i=0;i<ndom+1;i++) fprintf(fp,"%le ",*(gl.xif+i));
		fprintf(fp,"\n%le %le %le %le\n",M,R,X,Z);
		fprintf(fp,"%le %d %le\n",Xc,conv,surff);
		fprintf(fp,"%le %le\n",Tc,pc);
		fprintf(fp,"%s\n",opa.name);
		fprintf(fp,"%s\n",eos.name);
		fprintf(fp,"%s\n",nuc.name);
		fprintf(fp,"%s\n",atm_name);
	}
	phi.write(fp,mode);
	p.write(fp,mode);
	T.write(fp,mode);
	
	fclose(fp);
		
}

int star1d::read(const char *input_file){

	FILE *fp,*fp2;
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


