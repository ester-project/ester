#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "star.h"
#include "matplotlib.h"

#include <string.h>
#include <stdlib.h>
#include <iomanip>

star1d::star1d() {
}

star1d::~star1d() {
}

star1d::star1d(const star1d &A) : star2d(A) {
}

star1d &star1d::operator=(const star1d &A) {
	star2d::operator=(A);
	
	return *this;

}

void star1d::write_tag(OUTFILE *fp) const {
	char tag[7]="star1d";
	
	fp->write("tag",tag,7);
}

std::string star1d::get_tag() const {
    return std::string("star1d");
}

double star1d::roundToPrecision(double value, int decimalPlaces){
return star2d::roundToPrecision(value,decimalPlaces);
}

bool star1d::check_tag(const char *tag) const {
	if(strcmp(tag,"star1d")) return false;
	return true;
}
/*
int star1d::read_old(const char *input_file) {
	FILE *fp;
	char tag[7],mode,*c;
	int ndom,i;
	
	if((fp=fopen(input_file,"rb")) == NULL) {
		return 1;
	}
	if (fread(tag,1,7,fp) < 7)
        ester_warn("could not read tag");
	tag[6]='\0';
	if(strcmp(tag,"star1d")) {
		return 1;
	}
	if (fread(&mode,1,1,fp) < 1)
        ester_warn("could not read mode");
	fclose(fp);

	if(mode=='b')
		fp=fopen(input_file,"rb");
	else
		fp=fopen(input_file,"rt");
	
	if (fread(tag,1,7,fp) < 7)
        ester_warn("could not read tag");
	if (fread(&mode,1,1,fp) < 1)
        ester_warn("could not read mode");
	
	if(mode=='b') {
		if (fread(&ndom,sizeof(int),1,fp) < 1)
            ester_warn("could not read ndom");
		map.gl.set_ndomains(ndom);
		if (fread(map.gl.npts,sizeof(int),ndom,fp) < (size_t) ndom)
            ester_warn("could not read npts");
		map.leg.npts=1;
		map.init();
		if (fread(map.R.data(),sizeof(double),ndom+1,fp) < (size_t) ndom+1)
            ester_warn("could not read map.R");
		if (fread(&M,sizeof(double),1,fp) < 1)
            ester_warn("could not read M");
		if (fread(&R,sizeof(double),1,fp) < 1)
            ester_warn("could not read R");
		if (fread(&X0,sizeof(double),1,fp) < 1)
            ester_warn("could not read X0");
		if (fread(&Z0,sizeof(double),1,fp) < 1)
            ester_warn("could not read Z0");
		if (fread(&Xc,sizeof(double),1,fp) < 1)
            ester_warn("could not read Xc");
		if (fread(&conv,sizeof(int),1,fp) < 1)
            ester_warn("could not read conv");
		if (fread(&surff,sizeof(double),1,fp) < 1)
            ester_warn("could not read surff");
		if (fread(&Tc,sizeof(double),1,fp) < 1)
            ester_warn("could not read Tc");
		if (fread(&pc,sizeof(double),1,fp) < 1)
            ester_warn("could not read pc");
		c=opa.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=eos.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=nuc.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=atm.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
	} else {
		if (fscanf(fp,"\n%d ",&ndom) < 1)
            ester_warn("could not read ndom");
		map.gl.set_ndomains(ndom);
		for(i=0;i<ndom;i++)  {
            if (fscanf(fp,"%d ",(map.gl.npts+i)) < 1)
                ester_warn("could not read %dth map.gl.npts", i);
        }
		map.leg.npts=1;
		map.init();
		for(i=0;i<ndom+1;i++) {
            if (fscanf(fp,"%le ",&map.R(i)) < 1)
                ester_warn("could not %dth map.R", i);
        }
		if (fscanf(fp,"\n%le %le %le %le\n",&M,&R,&X0,&Z0) < 4)
            ester_warn("could not read M, R X0 or Z0");
		if (fscanf(fp,"%le %d %le\n",&Xc,&conv,&surff) < 3)
            ester_warn("could not read Xc, conv or surff");
		if (fscanf(fp,"%le %le\n",&Tc,&pc) < 2)
            ester_warn("could not read Tc, pc");
		if (fscanf(fp,"%s\n",opa.name) < 1)
            ester_warn("could not read opa.name");
		if (fscanf(fp,"%s\n",eos.name) < 1)
            ester_warn("could not read eos.name");
		if (fscanf(fp,"%s\n",nuc.name) < 1)
            ester_warn("could not read nuc.name");
		if (fscanf(fp,"%s\n",atm.name) < 1)
            ester_warn("could not read atm.name");
	}
	
	map.remap();
	phi.read(nr,1,fp,mode);
	p.read(nr,1,fp,mode);
	T.read(nr,1,fp,mode);
	
	core_convec=1;
	env_convec=0;
	min_core_size=0.03;
	version.major=0;
	version.minor=0;
	version.rev=0;
	version.svn=1;
    version.name = std::string("0.0 svn 1");
	domain_type.resize(ndomains);
	for(int n=0;n<ndomains;n++) {
		if(n<conv) domain_type[n]=CORE;
		else domain_type[n]=RADIATIVE;
	}
	
	fclose(fp);
	fill();
	return 0;
}*/

int star1d::read(const char *input_file, int dim) {
    return star2d::read(input_file, 1);
}

int star1d::init(const char *input_file,const char *param_file,int argc,char *argv[]) {
	cmdline_parser cmd;
	file_parser fp;
	char *arg,*val,default_params[256];
	mapping map0;
	int i,k,change_grid=0;
	matrix Tr;

	sprintf(default_params,"%s/ester/1d_default.par", ESTER_DATADIR);

	if(*input_file) {
		if(read(input_file, 1)) {
			printf("Error reading input file: %s\n",input_file);
			return 1;
		}
		map0=map;
	} else {
		if(!fp.open(default_params)) { 
			printf("Can't open default parameters file %s\n",default_params);
			return 1;
		}
		else {
			while((k=fp.get(arg,val))) {
             	if (strcmp(arg, "Z") == 0){
					    				
					initial_composition(0.7,0.02); // just so Zmix is defined in the global map 
					
					std::string strVal = std::to_string(roundToPrecision(global_abundance_map.Zmix, 9)); // converting to char* for check_arg()
					val = strdup(strVal.c_str());             		
             	}

				if((i=check_arg(arg,val,&change_grid))) {
					printf("Syntax error in parameters file %s, line %d\n",param_file,k);
					if(i==2) {
						printf("Error: Argument to '%s' missing\n",arg);
						exit(EXIT_FAILURE);
					}
					if(i==1) {
						printf("Unknown parameter %s\n",arg);
						exit(EXIT_FAILURE);
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
			return 1;
		}
		else {
			while((k=fp.get(arg,val))) {
				if((i=check_arg(arg,val,&change_grid))) {
					printf("Sintax error in parameters file %s, line %d\n",param_file,k);
					if(i==2) {
						printf("Error: Argument to '%s' missing\n",arg);
						exit(EXIT_FAILURE);
					}
					if(i==1) {
						printf("Unknown parameter %s\n",arg);
						exit(EXIT_FAILURE);
					}
				}
			}
			fp.close();
		}
	}

	cmd.open(argc,argv);
	while(int err_code=cmd.get(arg,val)) {
				
             	if (strcmp(arg, "mixture") == 0 && age == 0){
				
				if (strcmp(mixture.name,val)==0) {
    				printf("[INFO] Mixture '%s' matches input argument.\n", val);
					
				} else {
				printf("mixture in cmd argument is different from default/input file\n");
				user_specified_mix = true;
				}

             	}
			 
			 	if (strcmp(arg, "Z") == 0 && age == 0){

					if (fabs(atof(val) - Z0) > (1e-12 + 1e-12 * Z0)) {
						// if the user has kept the argument Z input same as before
						// this avoids re-doing comp unnecessarily
						user_specified_Z = true;
					}
					
				}

			 	if (strcmp(arg, "X") == 0 && age == 0){

					if (fabs(atof(val) - X0) > (1e-12 + 1e-12 *X0)) {
						// if the user has kept the argument X input same as before
						// this avoids re-doing comp unnecessarily
						user_specified_X = true;
					}
				}

		if(err_code==-1) exit(1);
		err_code=check_arg(arg,val,&change_grid);
		if(err_code==2) {
			fprintf(stderr,"Error: Argument to '%s' missing\n",arg);
            exit(EXIT_FAILURE);
		}
		if(err_code==1) {
			fprintf(stderr,"Unknown parameter '%s'\n",arg);
            exit(EXIT_FAILURE);
		}
		cmd.ack(arg,val);
	}
	cmd.close();

	if (user_specified_mix && age == 0) { // i.e. mixture has changed
		global_abundance_map.mixture_name = mixture.name; // update the global mixture name
	}
	if ((user_specified_X + user_specified_Z + user_specified_mix) > 0 && age == 0 && loaded_from_file==1) {

		printf("X,Z, or mixture has changed to ZAMS file, updating comp...\n");

		(void)parse_composition_data();  // updates global_abundance_map.comp_data
    	init_comp();               // rebuild stratified + cno core composition

	} else if ((user_specified_X + user_specified_Z + user_specified_mix) > 0 && age != 0 && loaded_from_file==1) {

		printf("mixture, X, or Z were given in command line that is different from *evolved* input model\n");
		printf("this option is not allowed, using input model's parameters\n");
	}
	
	if((change_grid&1)&&!(change_grid&2)) {
		fprintf(stderr,"Must specify number of points per domain (npts)\n");
		exit(1);
	}
	if (*input_file) {
		if(change_grid) {
			mapping map_new;
			map_new=map;
			map=map0;
			remap(map_new.ndomains,map_new.gl.npts,map_new.nt,map_new.nex);
		}
	} else {
		map.leg.npts=1;
		map.init();
		T=1-0.5*r*r;
		p=T;
		phi=-T;
		vr=0*T;
		vt=0*T;
		w=0*T;
		conv=0;
		domain_type.resize(ndomains);
		for(int n=0;n<ndomains;n++) domain_type[n]=RADIATIVE;
		domain_weight = init_domain_weight(domain_type);
		phiex=zeros(map.nex,map.nt);
	}
	
	//init_comp();
	//fill();

	if (!loaded_from_file) {
		(void)parse_composition_data();  // updates/creates global_abundance_map.comp_data
		if (!user_specified_Z && user_specified_mix){
			printf("mixture is given but Z is not on cmd line!\n");
		Z0 = roundToPrecision(global_abundance_map.Zmix, 9); 
		}
		init_comp();
	}

	fill();

	return 0;
}

int star1d::check_arg(char *arg,char *val,int *change_grid) {
	if(!strcmp(arg,"nth")) {
		return 1;
	} else if(!strcmp(arg,"nex")) {
		return 1;
	} else if(!strcmp(arg,"Omega_bk")) {
		return 1;
	} else if(!strcmp(arg,"Ekman")) {
		return 1;
	}
	return star2d::check_arg(arg,val,change_grid);

}

void star1d::dump_info() {
    printf("ESTER 1d model file");
    printf(" (Version %s)\n", version.name.c_str());
	// printf("\n1d ESTER model file  (Version %d.%d rev %d",version.major,version.minor,version.rev);
	// if(version.svn) printf(".svn");
	// printf(")\n\n");
	
	printf("General parameters:\n\n");
	printf("\tMass = %.5f Msun (%e g)\n",M/M_SUN,M);
	printf("\tRadius = %.5f Rsun (%e cm)\n",R/R_SUN,R);
	printf("\tLuminosity = %.4f Lsun (%e erg/s)\n",luminosity()/L_SUN,luminosity());
	printf("\tTeff = %.2f\n",Teff()(0));
	printf("\tlog(geff) = %.4f\n",log10(gsup())(0));
	printf("\tX0=%.4f   Y0=%.4f   Z0=%.4f\n",X0,Y0,Z0);
	printf("\n");
	
	if(conv==0) printf("No convective core\n\n");
	else {
		printf("Convective core:\n\n");
		double mcc=Mcore();
		printf("\tMass_core = %.5f Msun (%e g)\n",mcc/M_SUN,mcc);
		double rcc=Rcore()(0);
		printf("\tRadius_core (p) = %.5f Rsun (%e cm)\n",rcc/R_SUN,rcc);
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
	printf("\n");
	
	printf("Additional parameters:\n\n");
	printf("\tOpacity = %s\n",opa.name);
    printf("\tAbundance Mixture = %s\n", mixture.name);    
	printf("\tEquation of state = %s\n",eos.name);
	printf("\tNuclear reactions = %s\n",nuc.name);
	printf("\tAtmosphere = %s\n",atm.name);
	printf("\tsurff = %e\n",surff);
	printf("\tcore_convec = %d\n",core_convec);
	printf("\tenv_convec = %d\n",core_convec);
	printf("\tmin_core_size = %e\n",min_core_size);
	printf("\n");
	
	printf("Tests:\n\n");
	printf("\tVirial test = %e\n",virial());
	printf("\tEnergy test = %e\n",energy_test());
	printf("\n");
	
}

matrix star1d::spectrum(const matrix& var) {
    matrix spec;
    spec = (map.gl.P, var, map.leg.P_00);
    int j = 0;
    for (int i=0; i<ndomains; i++) {
        spec.setblock(j, j+map.gl.npts[i]-1, 0, -1,
                spec.block(j, j+map.gl.npts[i]-1, 0, -1)/max(spec.block(j, j+map.gl.npts[i]-1, 0, -1))
                );
        j += map.gl.npts[i];
    }
    return abs(spec);
}

void star1d::plot(const matrix_map& error) {

    plt::clf();

    plt::subplot(231);
    // plt::title(std::string("iter: ") + std::to_string(nit));
    plt::plot(r, rho, "$\\rho$");
    plt::plot(r, T, "$T$");
    plt::plot(r, p, "$p$");
    plt::legend();
    for (int i=0; i<ndomains; i++) {
        plt::axvline(map.gl.xif[i]);
    }
    plt::axvline(1.0);

    plt::subplot(232);
    // plt::title(std::string("iter: ") + std::to_string(nit));
    plt::plot(r, phi, "$\\Phi$");
    plt::legend();

    plt::subplot(233, true);
    std::ostringstream str_stream;

    str_stream.clear();
    str_stream.str("");
    str_stream << Tc;
    plt::text(0.0, .3, std::string("$T_c$:   ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pc;
    plt::text(0.0, .2, std::string("$p_c$:   ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << rhoc;
    plt::text(0.0, 0.1, std::string("$\\rho_c$:  ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pi_c;
    plt::text(0.0, 0.0, std::string("$\\pi_c$: ") + str_stream.str());

    if (error["Phi"].ncols()*error["Phi"].nrows() > 0 && error["Phi"](0) > .0) {
        plt::subplot(223);
        plt::title("Error");
        plt::semilogy(error["Phi"], "$\\Phi$");
        plt::semilogy(error["log_p"], "$\\ln p$");
        plt::semilogy(error["log_T"], "$\\ln T$");
        // plt::semilogy(error["log_pc"], "error $log_{p_c}$");
        plt::semilogy(error["log_Tc"], "$\\ln T_c$");
        // plt::semilogy(error["Ri"], "error $R_i$");
        plt::legend("lower left");
    }

    plt::subplot(224);
    plt::title("Spectrum");
    plt::semilogy(spectrum(rho), "$\\rho$");
    int n = 0;
    plt::axvline(n);
    for (int i=0; i<ndomains; i++) {
        n += map.gl.npts[i];
        plt::axvline(n);
    }
    plt::legend();

    plt::draw();
    plt::pause();
}
