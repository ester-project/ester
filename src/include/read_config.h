#include <stdlib.h>
#include <string.h>
#include <stdio.h>

class configuration {
public:
	int minit,maxit;
	double tol,newton_dmax;
	int verbose;
	char plot_device[64];
	double plot_interval;
	char input_file[256];
	char param_file[256];
	char output_file[256];

	configuration(int argc, char *argv[]);
	~configuration(){};
	void read_config_file();
	void read_command_line(int argc, char *argv[]);
	void missing_argument(const char *arg);
	int check_arg(const char *arg,const char *val);
    bool noplot;
};

