#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "star.h"
#include "matplotlib.h"

#include <string.h>
#include <stdlib.h>

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

bool star1d::check_tag(const char *tag) const {
    if(strcmp(tag, "star1d")) return false;
    return true;
}

int star1d::read(const char *input_file, int dim) {
    return star2d::read(input_file, 1);
}

int star1d::init(const char *input_file, const char *param_file, int argc, char *argv[]) {
    cmdline_parser cmd;
    file_parser fp;
    char *arg, *val, default_params[256];
    mapping map0;
    int i, k, change_grid=0;
    matrix Tr;

    sprintf(default_params, "%s/ester/1d_default.par", ESTER_DATADIR);

    if(*input_file) {
        if(read(input_file, 1)) {
            printf("Error reading input file: %s\n", input_file);
            return 1;
        }
        map0=map;
    } else {
        if(!fp.open(default_params)) {
            printf("Can't open default parameters file %s\n", default_params);
            return 1;
        }
        else {
            while((k=fp.get(arg, val))) {
                if((i=check_arg(arg, val, &change_grid))) {
                    printf("Syntax error in parameters file %s, line %d\n", param_file, k);
                    if(i==2) {
                        printf("Error: Argument to '%s' missing\n", arg);
                        exit(EXIT_FAILURE);
                    }
                    if(i==1) {
                        printf("Unknown parameter %s\n", arg);
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
            printf("Can't open parameters file %s\n", param_file);
            return 1;
        }
        else {
            while((k=fp.get(arg, val))) {
                if((i=check_arg(arg, val, &change_grid))) {
                    printf("Syntax error in parameters file %s, line %d\n", param_file, k);
                    if(i==2) {
                        printf("Error: Argument to '%s' missing\n", arg);
                        exit(EXIT_FAILURE);
                    }
                    if(i==1) {
                        printf("Unknown parameter %s\n", arg);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            fp.close();
        }
    }

    cmd.open(argc, argv);
    while(int err_code=cmd.get(arg, val)) {
        if(err_code==-1) exit(1);
        err_code=check_arg(arg, val, &change_grid);
        if(err_code==2) {
            fprintf(stderr, "Error: Argument to '%s' missing\n", arg);
            exit(EXIT_FAILURE);
        }
        if(err_code==1) {
            fprintf(stderr, "Unknown parameter '%s'\n", arg);
            exit(EXIT_FAILURE);
        }
        cmd.ack(arg, val);
    }
    cmd.close();

    if((change_grid&1)&&!(change_grid&2)) {
        fprintf(stderr, "Must specify number of points per domain (npts)\n");
        exit(1);
    }
    if (*input_file) {
        if(change_grid) {
            mapping map_new;
            map_new=map;
            map=map0;
            remap(map_new.ndomains, map_new.gl.npts, map_new.nt, map_new.nex);
        }
    } else {
        map.leg.npts=1;
        map.init();
        T=1-0.5*r*r;
        p=T;
        phi=-T;
        if (config.init_poly) {
            double n = 1.5;
            n = 3.;
            matrix h = solve_poly1d(n, 1e-10, 40, 1e-1);

            mapping map;
            map.set_ndomains(1);
            map.set_npts(40);
            map.gl.set_xif(0., 1.);
            map.set_nt(1);
            map.init();

            h = map.gl.eval(h, this->map.r);

            p = pow(h, n+1);
            rho = pow(h, n);
            T = p/rho;
        }
        G=0*T;
        w=0*T;
        conv=0;
        domain_type.resize(ndomains);
        for(int n=0;n<ndomains;n++) domain_type[n]=RADIATIVE;
        phiex=zeros(map.nex, map.nt);
    }

    init_comp();
    fill();

    phi = solve_phi();

    return 0;
}

int star1d::check_arg(char *arg, char *val, int *change_grid) {
    if(!strcmp(arg, "nth")) {
        return 1;
    } else if(!strcmp(arg, "nex")) {
        return 1;
    } else if(!strcmp(arg, "Omega_bk")) {
        return 1;
    } else if(!strcmp(arg, "Ekman")) {
        return 1;
    }
    return star2d::check_arg(arg, val, change_grid);

}

void star1d::dump_info() {
    printf("ESTER 1d model file");
    printf(" (Version %s)\n", version.name.c_str());

    printf("General parameters:\n\n");
    printf("\tMass = %.5f Msun (%e g)\n", M/M_SUN, M);
    printf("\tRadius = %.5f Rsun (%e cm)\n", R/R_SUN, R);
    printf("\tLuminosity = %.4f Lsun (%e erg/s)\n", luminosity()/L_SUN, luminosity());
    printf("\tTeff = %.2f\n", Teff()(0));
    printf("\tlog(geff) = %.4f\n", log10(gsup())(0));
    printf("\tX0=%.4f   Y0=%.4f   Z0=%.4f\n", X0, Y0, Z0);
    printf("\n");

    if(conv==0) printf("No convective core\n\n");
    else {
        printf("Convective core:\n\n");
        double mcc=Mcore();
        printf("\tMass_core = %.5f Msun (%e g)\n", mcc/M_SUN, mcc);
        double rcc=Rcore()(0);
        printf("\tRadius_core (p) = %.5f Rsun (%e cm)\n", rcc/R_SUN, rcc);
        printf("\tX_core/X_env = %.4f\n", Xc);
        printf("\n");
    }
    printf("Central values:\n\n");
    printf("\tTemperature = %e K\n", Tc);
    printf("\tDensity = %e g/cm3\n", rhoc);
    printf("\tPressure = %e dyn/cm2\n", pc);
    printf("\n");

    printf("Grid parameters:\n\n");
    printf("\t # of domains = %d\n", ndomains);
    printf("\t # of domains in convective core = %d\n", conv);
    printf("\t nr = %d    (", nr);
    for(int n=0;n<ndomains;n++) {
        printf("%d", map.gl.npts[n]);
        if(n<ndomains-1) printf(", ");
    }
    printf(")\n");
    printf("\n");

    printf("Additional parameters:\n\n");
    printf("\tOpacity = %s\n", opa.name);
    printf("\tEquation of state = %s\n", eos.name);
    printf("\tNuclear reactions = %s\n", nuc.name);
    printf("\tAtmosphere = %s\n", atm.name);
    printf("\tsurff = %e\n", surff);
    printf("\tcore_convec = %d\n", core_convec);
    printf("\tenv_convec = %d\n", core_convec);
    printf("\tmin_core_size = %e\n", min_core_size);
    printf("\n");

    printf("Tests:\n\n");
    printf("\tVirial test = %e\n", test_virial);
    printf("\tEnergy test = %e\n", test_energy);
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
    plt::plot(r, rho, "rho");
    plt::plot(r, T, "T");
    plt::plot(r, p, "p");
    plt::legend();
    for (int i=0; i<ndomains; i++) {
        plt::axvline(map.gl.xif[i]);
    }
    plt::axvline(1.0);

    plt::subplot(232);
    // plt::title(std::string("iter: ") + std::to_string(nit));
    plt::plot(r, phi, "Phi");
    plt::legend();

    plt::subplot(233, true);
    std::ostringstream str_stream;

    str_stream.clear();
    str_stream.str("");
    str_stream << Tc;
    plt::text(0.0, .3, std::string("T_c:   ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pc;
    plt::text(0.0, .2, std::string("p_c:   ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << rhoc;
    plt::text(0.0, 0.1, std::string("rho_c:  ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pi_c;
    plt::text(0.0, 0.0, std::string("pi_c: ") + str_stream.str());

    if (error["Phi"].ncols()*error["Phi"].nrows() > 2 && error["Phi"](0) > .0) {
        plt::subplot(223);
        plt::semilogy(error["Phi"], "Phi");
        plt::semilogy(error["log_p"], "ln p");
        plt::semilogy(error["log_T"], "ln T");
        // plt::semilogy(error["log_pc"], "error $log_{p_c}$");
        plt::semilogy(error["log_Tc"], "ln T_c");
        // plt::semilogy(error["Ri"], "error $R_i$");
        plt::legend("lower left");
        plt::title("Error");
    }

    plt::subplot(224);
    plt::semilogy(spectrum(rho), "rho");
    plt::title("Spectrum");
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
