#include "ester-config.h"
#include "utils.h"
#include "parser.h"
#include "read_config.h"

configuration::configuration(int argc, char *argv[]) {
    int i, k;
    char *arg, *val;
    char file[256];
    cmdline_parser cmd;
    file_parser fp;

    verbose=1;
    strcpy(plot_device, "/NULL");
    plot_interval=10;
    strcpy(output_file, "star.out");
    *input_file=0;
    *param_file=0;
    output_mode='b';
    minit=1;
    maxit=200;
    tol=1e-8;
    newton_dmax=0.5;

    sprintf(file, "%s/ester/star.cfg", ESTER_DATADIR);
    if (!fp.open(file)) 
        printf("Can't open configuration file %s\n", file);
    else {
        while ((k=fp.get(arg, val))) {
            if ((i=check_arg(arg, val))) {
                printf("Syntax error in configuration file %s, line %d\n",
                        file, k);
                if (i==2) missing_argument(arg);
                if (i==1) {
                    printf("Unknown parameter %s\n", arg);
                    exit(1);
                }
            }
        }
        fp.close();
    }

    cmd.open(argc, argv);
    while (int err_code=cmd.get(arg, val)) {
        if (err_code == -1) exit(1);
        err_code = check_arg(arg, val);
        if (err_code == 2) missing_argument(arg);
        if (err_code == 0) cmd.ack(arg, val);
    }
    cmd.close();

}

int configuration::check_arg(const char *arg, const char *val) {

    int err=0;

    if (!strcmp(arg, "v0"))
        verbose=0;
    else if (!strcmp(arg, "v1"))
        verbose=1;
    else if (!strcmp(arg, "v2"))
        verbose=2;
    else if (!strcmp(arg, "v3"))
        verbose=3;
    else if (!strcmp(arg, "v4"))
        verbose=4;
    else if (!strcmp(arg, "verbose")) {
        if (val==NULL) return 2;
        verbose=atoi(val);
        verbose=verbose>4?4:verbose;
        verbose=verbose<0?0:verbose;
    } 
    else if (!strcmp(arg, "o")||!strcmp(arg, "output_file")) {
        if (val==NULL) return 2;
        strcpy(output_file, val);
    }
    else if (!strcmp(arg, "i")||!strcmp(arg, "input_file")) {
        if (val==NULL) return 2;
        strcpy(input_file, val);
    }  
    else if (!strcmp(arg, "p")||!strcmp(arg, "param_file")) {
        if (val==NULL) return 2;
        strcpy(param_file, val);
    }  
    else if (!strcmp(arg, "ascii")) 
        output_mode='t';
    else if (!strcmp(arg, "binary"))
        output_mode='b';
    else if (!strcmp(arg, "output_mode")) {
        if (val==NULL) return 2;
        if (val[0]!='b'&&val[0]!='t') 
            printf("Ignoring unknown output_mode %s\n", val);
        else output_mode=val[0];
    }
    else if (!strcmp(arg, "plot_interval")) {
        if (val==NULL) return 2;
        plot_interval=atof(val);
    }
    else if (!strcmp(arg, "plot_device")) {
        if (val==NULL) return 2;
        strcpy(plot_device, val);
    }
    else if (!strcmp(arg, "noplot")) {
        strcpy(plot_device, "/NULL");
    }
    else if (!strcmp(arg, "maxit")) {
        if (val==NULL) return 2;
        maxit=atoi(val);
    }
    else if (!strcmp(arg, "minit")) {
        if (val==NULL) return 2;
        minit=atoi(val);
    }
    else if (!strcmp(arg, "tol")) {
        if (val==NULL) return 2;
        tol=atof(val);
    }
    else if (!strcmp(arg, "newton_dmax")) {
        if (val==NULL) return 2;
        newton_dmax=atof(val);
    }
    else err=1;

    return err;

}

void configuration::missing_argument(const char *arg) {
    ester_err("Error: Argument to '%s' missing", arg);
    exit(1);
}
