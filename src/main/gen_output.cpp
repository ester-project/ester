#include "ester-config.h"
#include "utils.h"
#include "star.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void output2d(const star2d &);
void write(const star2d &A,char *var,char *fmt);
void readconf(char *line);
void matrix_fmt(char *fmt,matrix a);

int dim=1,pole=0,equator=0,transpose=0;

int main(int argc,char *argv[]) {

	char *input_file;
	star1d m1d;
	star2d m2d;

	input_file=argv[1];
	if(!m1d.read(input_file, 1)) output2d(m1d);
	else if (!m2d.read(input_file)) output2d(m2d);
	else {
		ester_err("Error reading input file: %s",input_file);
		return 1;
	}
	
	return 0;
}

void readconf(char *line) {

	char *var,*val;

	if(strncmp(line,"\\conf",5)) return;
	line+=5;
	while(*line!='{') {
		if(*line=='\0') return;
		line++;
	}
	line++;
	var=strtok(line,"=");
	while(var) {
		val=strtok(NULL,",}");
		if(!strcmp(var,"dim")) dim=atoi(val);
		if(!strcmp(var,"pole")) pole=atoi(val);
		if(!strcmp(var,"equator")) equator=atoi(val);
		if(!strcmp(var,"transpose")) transpose=atoi(val);
		var=strtok(NULL,"=");
	}	
}

void matrix_fmt(char *fmt,matrix a) {

	int i,j;
	double *p;

	if(transpose) a=a.transpose();
	p=a.data();
	for(i=0;i<a.ncols();i++) {
		for(j=0;j<a.nrows();j++) {
			fprintf(stdout,fmt,*(p++));
			if(j<a.nrows()-1) fprintf(stdout," ");
		}
		fprintf(stdout,"\n");
	}
}


void output2d(const star2d &A) {

	FILE *fp;
	char line[1024],*p,*p2,*var,*fmt;
	char *saveptr1,*saveptr2,*saveptr3;
	matrix rho0;
	
	fp=stdin;
	while(fgets(line,1024,fp)) {
		if(*line=='\n') fprintf(stdout,"\n");
		if(*line=='\\') {
			readconf(line);
			continue;
		}
		if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';
		if(*line=='$') {
			p=strtok_r(line,"$",&saveptr1);
		} else {
			p=strtok_r(line,"$",&saveptr1);
			if(p) fprintf(stdout,"%s",p);
			p=strtok_r(NULL,"$",&saveptr1);
		}
		while(p) {
			p++;
			var=strtok_r(p,"}",&saveptr2);
			var=strtok_r(var,",",&saveptr3);
			fmt=strtok_r(NULL,",",&saveptr3);
			p2=strtok_r(NULL,"}",&saveptr2);
			write(A,var,fmt);
			if(p2) fprintf(stdout,"%s",p2);
			p=strtok_r(NULL,"$",&saveptr1);
		}
	}
}

void write(const star2d &A,char *var,char *fmt) {

	int i;
	double d;
	matrix m,m2,T,T_odd;

	if(pole||equator) {
		m2=zeros(1,A.nth+pole+equator);
		m2.setblock(0,0,equator,A.nth+equator-1,A.th);
		if(equator) m2(0)=PI/2;
		m=A.map.leg.eval_00(A.th,m2,T);
		m=A.map.leg.eval_11(A.th,m2,T_odd);
	} else {
		T=eye(A.nth);
		T_odd=T;
	}
	
	if(!strcmp(var,"nr")) {
		if(fmt) fprintf(stdout,fmt,A.nr);
		else {
			i=A.nr;
			fwrite(&i,sizeof(int),1,stdout);
		}
	} else if(!strcmp(var,"ndomains")) {
		if(fmt) fprintf(stdout,fmt,A.ndomains);
		else {
			i=A.ndomains;
			fwrite(&i,sizeof(int),1,stdout);
		}
	} else if(!strcmp(var,"npts")) {
		if(fmt) {
			for(i=0;i<A.ndomains;i++) {	
				fprintf(stdout,fmt,*(A.map.gl.npts+i));
				if(i<A.ndomains-1) fprintf(stdout,",");
			}
		} else {
			fwrite(A.map.gl.npts,sizeof(int),A.ndomains,stdout);
		}
	} else if(!strcmp(var,"xif")) {
		if(fmt) {
			for(i=0;i<A.ndomains+1;i++) {	
				fprintf(stdout,fmt,*(A.map.gl.xif+i));
				if(i<A.ndomains) fprintf(stdout,",");
			}
		} else {
			fwrite(A.map.gl.xif,sizeof(double),A.ndomains+1,stdout);
		}
	} else if(!strcmp(var,"nth")) {
		if(fmt) fprintf(stdout,fmt,A.nth);
		else {
			i=A.nth;
			fwrite(&i,sizeof(int),1,stdout);
		}
	} else if(!strcmp(var,"nex")) {
		if(fmt) fprintf(stdout,fmt,A.nex);
		else {
			i=A.nex;
			fwrite(&i,sizeof(int),1,stdout);
		}
	} else if(!strcmp(var,"ps")) {
		m=A.ps;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"m")) {
		if(fmt) fprintf(stdout,fmt,A.m);
		else {
			fwrite(&A.m,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"surff")) {
		if(fmt) fprintf(stdout,fmt,A.surff);
		else {
			fwrite(&A.surff,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"conv")) {
		if(fmt) fprintf(stdout,fmt,A.conv);
		else {
			fwrite(&A.conv,sizeof(int),1,stdout);
		}
	} else if(!strcmp(var,"Omega")) {
		if(dim) d=A.Omega*A.units.Omega;
		else d=A.Omega;
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Omega_bk")) {
		if(fmt) fprintf(stdout,fmt,A.Omega_bk);
		else {
			fwrite(&A.Omega_bk,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Omegac")) {
		if(dim) d=A.Omegac*A.units.Omega;
		else d=A.Omegac;
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Xc")) {
		if(fmt) fprintf(stdout,fmt,A.Xc);
		else {
			fwrite(&A.Xc,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"rhoc")) {
		if(fmt) fprintf(stdout,fmt,A.rhoc);
		else {
			fwrite(&A.rhoc,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Tc")) {
		if(fmt) fprintf(stdout,fmt,A.Tc);
		else {
			fwrite(&A.Tc,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"pc")) {
		if(fmt) fprintf(stdout,fmt,A.pc);
		else {
			fwrite(&A.pc,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"R")) {
		if(fmt) fprintf(stdout,fmt,A.R);
		else {
			fwrite(&A.R,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Rp")) {
		if(fmt) fprintf(stdout,fmt,A.R);
		else {
			fwrite(&A.R,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Re")) {
		d=A.map.leg.eval_00(A.r.row(A.nr-1),PI/2)(0)*A.units.r;
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"M")) {
		if(fmt) fprintf(stdout,fmt,A.M);
		else {
			fwrite(&A.M,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Lz")) {
		d=A.Lz();
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Mcore")) {
		d=A.Mcore();
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Lzcore")) {
		d=A.Lzcore();
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"X")) {
		if(fmt) fprintf(stdout,fmt,A.X0);
		else {
			fwrite(&A.X0,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Y")) {
		if(fmt) fprintf(stdout,fmt,A.Y0);
		else {
			fwrite(&A.Y0,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Z")) {
		if(fmt) fprintf(stdout,fmt,A.Z0);
		else {
			fwrite(&A.Z0,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"R/R_SUN")) {
		if(fmt) fprintf(stdout,fmt,A.R/R_SUN);
		else {
			d=A.R/R_SUN;
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Rp/R_SUN")) {
		if(fmt) fprintf(stdout,fmt,A.R/R_SUN);
		else {
			d=A.R/R_SUN;
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Re/R_SUN")) {
		d=A.map.leg.eval_00(A.r.row(A.nr-1),PI/2)(0)*A.units.r/R_SUN;
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"M/M_SUN")) {
		if(fmt) fprintf(stdout,fmt,A.M/M_SUN);
		else {
			d=A.M/M_SUN;
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"opa")) {
		if(fmt) fprintf(stdout,fmt,A.opa.name);
		else {
			fwrite(A.opa.name,sizeof(char),strlen(A.opa.name)+1,stdout);
		}
	} else if(!strcmp(var,"eos")) {
		if(fmt) fprintf(stdout,fmt,A.eos.name);
		else {
			fwrite(A.eos.name,sizeof(char),strlen(A.eos.name)+1,stdout);
		}
	} else if(!strcmp(var,"r")) {
		if(dim) m=A.r*A.units.r;
		else m=A.r;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"z")) {
		if(dim) m=A.z*A.units.r;
		else m=A.z;
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"D")) {
		if(dim) m=A.D/A.units.r;
		else m=A.D;
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"I")) {
		if(dim) m=A.map.gl.I*A.units.r;
		else m=A.map.gl.I;
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"rex")) {
		if(dim) m=A.map.ex.r*A.units.r;
		else m=A.map.ex.r;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"phiex")) {
		if(dim) m=A.phiex*A.units.phi;
		else m=A.phiex;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"Dex")) {
		if(dim) m=A.map.ex.D/A.units.r;
		else m=A.map.ex.D;
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"th")) {
		m=zeros(1,A.nth+pole+equator);
		m.setblock(0,0,equator,A.nth-1+pole,A.th);
		if(equator) m(0)=PI/2;
		if(pole) m(m.ncols()-1)=0;
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"Dt")) {
		m2=(A.Dt,T_odd);
		m=zeros(A.nth+pole+equator,A.nth+pole+equator);
		m.setblock(equator,A.nth+equator-1,0,A.nth+pole+equator-1,m2);
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"Dtodd")) {
		m2=(A.map.leg.D_11,T);
		m=zeros(A.nth+pole+equator,A.nth+pole+equator);
		m.setblock(equator,A.nth+equator-1,0,A.nth+pole+equator-1,m2);
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"Dt2")) {
		m2=(A.Dt2,T);
		m=zeros(A.nth+pole+equator,A.nth+pole+equator);
		m.setblock(equator,A.nth+equator-1,0,A.nth+pole+equator-1,m2);
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"It")) {
		m=zeros(A.nth+pole+equator,1);
		m.setblock(equator,A.nth+equator-1,0,0,A.map.leg.I_00);
		if(fmt) matrix_fmt(fmt,m);
		else m.write(stdout,'b');
	} else if(!strcmp(var,"Ts")) {
		m=A.Ts;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"map.R")) {
		if(dim) m=A.map.R.block(1,-1,0,-1)*A.units.r;
		else m=A.map.R.block(1,-1,0,-1);
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"rho")) {
		if(dim) m=A.rho*A.units.rho;
		else m=A.rho;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"phi")) {
		if(dim) m=A.phi*A.units.phi;
		else m=A.phi;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"p")) {
		if(dim) m=A.p*A.units.p;
		else m=A.p;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"Xr")) {
		m=A.comp.X();
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"Yr")) {
		m=A.comp.Y();
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"Zr")) {
		m=A.comp.Z();
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strncmp(var,"X_",2)) {
		std::string elem(var+2);
		if(A.comp.count(elem)) m=A.comp[elem];
		else m=zeros(A.nr,A.nth);
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"w")) {
		if(dim) m=A.w*A.units.Omega;
		else m=A.w;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"T")) {
		if(dim) m=A.T*A.units.T;
		else m=A.T;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"vr")) {
		if(dim) m=A.vr*A.units.v;
		else m=A.vr;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"vt")) {
		if(dim) m=A.vt*A.units.v;
		else m=A.vt;
		if(fmt) matrix_fmt(fmt,(m,T_odd));
		else (m,T_odd).write(stdout,'b');
	} else if(!strcmp(var,"N2")) {
		m=A.N2();
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"opa.k")) {
		m=A.opa.k;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"opa.xi")) {
		m=A.opa.xi;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"opa.dlnxi_lnT")) {
		m=A.opa.dlnxi_lnT;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"opa.dlnxi_lnrho")) {
		m=A.opa.dlnxi_lnrho;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"nuc.eps")) {
		m=A.nuc.eps;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"nuc.pp")) {
		m=A.nuc.pp;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"nuc.cno")) {
		m=A.nuc.cno;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.G1")) {
		m=A.eos.G1;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.cp")) {
		m=A.eos.cp;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.del_ad")) {
		m=A.eos.del_ad;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.G3_1")) {
		m=A.eos.G3_1;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.cv")) {
		m=A.eos.cv;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.d")) {
		m=A.eos.d;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.prad")) {
		m=A.eos.prad;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.chi_T")) {
		m=A.eos.chi_T;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.chi_rho")) {
		m=A.eos.chi_rho;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"eos.s")) {
		m=A.eos.s;
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"L")) {
		d=A.luminosity();
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"L/L_SUN")) {
		d=A.luminosity()/L_SUN;
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Teff")) {
		m=A.Teff();
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"Teff_eq")) {
		m=A.Teff();
		d=A.map.leg.eval_00(m,PI/2)(0);
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"Teff_pol")) {
		m=A.Teff();
		d=A.map.leg.eval_00(m,0)(0);
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"gsup")) {
		m=A.gsup();
		if(fmt) matrix_fmt(fmt,(m,T));
		else (m,T).write(stdout,'b');
	} else if(!strcmp(var,"gsup_eq")) {
		m=A.gsup();
		d=A.map.leg.eval_00(m,PI/2)(0);
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"gsup_pol")) {
		m=A.gsup();
		d=A.map.leg.eval_00(m,0)(0);
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"eps")) {
		d=1-1./A.map.leg.eval_00(A.r,PI/2)(0);
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"eps_c")) {
		d=0;
		if(A.conv) d=1-A.map.eta(A.conv)/A.map.leg.eval_00(A.map.R.row(A.conv),PI/2)(0);
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"virial")) {
		d=A.virial();
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else if(!strcmp(var,"energy_test")) {
		d=A.energy_test();
		if(fmt) fprintf(stdout,fmt,d);
		else {
			fwrite(&d,sizeof(double),1,stdout);
		}
	} else {
		ester_err("Unknown variable %s", var);
		exit(1);
	}
}

