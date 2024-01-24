#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

#include "debug.h"
#include "utils.h"
#include "parser.h"
#include "read_config.h"
#include "matplotlib.h"

configuration::configuration() {
    // These are default values aimed to be replaced by config file values then command line ones
	verbose=1;
	strcpy(plot_device,"/NULL");
	plot_interval=10;
	strcpy(output_file,"star.out");
	*input_file=0;
	*param_file=0;
	minit=1;
	maxit=200;
	tol=1e-8;
	newton_dmax=0.5;
    noplot = false;
}


void configuration::read_config(int argc, char *argv[]) {
    // To ensure command line argument's precendence (complying with documentation)
    // over the configuration file argument these methods must be called IN THIS ORDER

    read_config_file();            // 1st --> default values
    read_command_line(argc, argv); // 2nd --> can erase default config file values 
}


void configuration::read_config_file() {
	char *arg,*val;
	int err_code, line;
	char file[256];
	file_parser fp;

	// Write in file the path to the default config file
	sprintf(file, "%s/ester/star.cfg", ESTER_DATADIR);

	// Try opening the file, get 1 on error
	if(fp.open(file)){
		printf("Can't open configuration file %s\n", file);
		ester_err(strerror(errno));
	}

	// iterate over each line/arg of the config file
	while((line = fp.get(arg, val))) {
		if((err_code = parse_arg(arg, val))) {
			printf("Syntax error in configuration file %s, line %d\n", file, line);
			if(err_code == 2)
				missing_argument(arg);
			if(err_code == 1)
				unknown_parameter(arg);
		}
	}
	fp.close();
}


void configuration::read_command_line(int argc, char *argv[]) {
    char *arg,*val;
    int err_code;
    cmdline_parser cmd;

    cmd.open(argc, argv);
    while(err_code = cmd.get(arg, val)) {
        if(err_code == -1)
            ester_err("Invalid argument %s", arg);

        err_code = parse_arg(arg,val);
        // We ignore Error Code 1 (unknown parameter) on purpose: it might be interpreted by star1d::init later
        if(err_code == 2)
            missing_argument(arg);
        if(err_code == 0)
            cmd.ack(arg,val);
    }
    cmd.close();

    if (noplot == false)
        plt::init();
}


int configuration::parse_arg(const char *arg,const char *val) {
    /*
    Return codes:
        0: No Error
        1: Unknown paramter
        2: Missing value
    */
	int err=0;

	if(!strcmp(arg,"v0"))
		verbose=0;
	else if(!strcmp(arg,"v1"))
		verbose=1;
	else if(!strcmp(arg,"v2"))
		verbose=2;
	else if(!strcmp(arg,"v3"))
		verbose=3;
	else if(!strcmp(arg,"v4"))
		verbose=4;
	else if(!strcmp(arg,"verbose")) {
		if(val==NULL) return 2;
		verbose=atoi(val);
		verbose=verbose>4?4:verbose;
		verbose=verbose<0?0:verbose;
	}
	else if(!strcmp(arg,"o")||!strcmp(arg,"output_file")) {
		if(val==NULL) return 2;
		strcpy(output_file,val);
	}
	else if(!strcmp(arg,"i")||!strcmp(arg,"input_file")) {
		if(val==NULL) return 2;
		strcpy(input_file,val);
	}
	else if(!strcmp(arg,"p")||!strcmp(arg,"param_file")) {
		if(val==NULL) return 2;
		strcpy(param_file,val);
	}
	else if(!strcmp(arg,"plot_interval")) {
		if(val==NULL) return 2;
		plot_interval=atof(val);
	}
	else if(!strcmp(arg,"plot_device")) {
		if(val==NULL) return 2;
		strcpy(plot_device,val);
	}
	else if(!strcmp(arg,"noplot")) {
		strcpy(plot_device,"/NULL");
        noplot = true;
        plt::init(true);
	}
	else if(!strcmp(arg,"maxit")) {
		if(val==NULL) return 2;
		maxit=atoi(val);
	}
	else if(!strcmp(arg,"minit")) {
		if(val==NULL) return 2;
		minit=atoi(val);
	}
	else if(!strcmp(arg,"tol")) {
		if(val==NULL) return 2;
		tol=atof(val);
	}
	else if(!strcmp(arg,"newton_dmax")) {
		if(val==NULL) return 2;
		newton_dmax=atof(val);
	}
	// else if(!strcmp(arg,"fpe")) {
    //     this->sigfpe = true;
	// }
	else err=1;

	return err;

}


void configuration::missing_argument(const char *arg) {
	ester_err("Argument to '%s' missing", arg);
}


void configuration::unknown_parameter(const char *arg) {
	ester_err("Unknown parameter '%s'", arg);
}
