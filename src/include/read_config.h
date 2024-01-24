#include <stdlib.h>
#include <string.h>
#include <stdio.h>

class configuration {
public:
    int minit,maxit;
    double tol,newton_dmax;
    int verbose;
    bool noplot;
    char plot_device[64];
    double plot_interval;
    char input_file[256];
    char param_file[256];
    char output_file[256];

    configuration();
    ~configuration(){};

    void read_config(int argc, char *argv[]); // wrapper around the following 2 methods
    void read_config_file();
    void read_command_line(int argc, char *argv[]);

    int parse_arg(const char *arg,const char *val);

    void missing_argument(const char *arg);
    void unknown_parameter(const char *arg);
};

