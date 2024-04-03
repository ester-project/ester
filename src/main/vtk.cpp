#include "star.h"

#include <cmath>
#include <stdio.h>

void help() {
    fprintf(stderr, "Usage: ester vtk <input-model> -o <output-file>\n");
}

int main(int argc, char *argv[]) {

    // model resolution radius, theta phi
    // odd number in theta to have value on the equator
    int n[3] = {32, 129, 256};

    char c;

    FILE *f = nullptr;
    char *input_model = nullptr;

    while ((c = getopt(argc, argv, "o:")) != -1) {
        switch (c) {
            case 'o':
                f = fopen(optarg, "w");
                break;
            default:
                help();
                return 1;
        }
    }

    for (auto index=optind; index<argc; index++) {
        if (input_model != nullptr) {
            help();
            return 1;
        }
        input_model = argv[index];
    }

    if (f == nullptr || input_model == nullptr) {
        help();
        return 1;
    }




    star2d A;
    if(A.read(input_model) == 0) {

        matrix Tr, Tt; // interpolation matrices
        matrix zetas, thetas; // matrices grid points coordinates
        zetas = ones(n[0], 1);
        thetas = ones(1, n[1]);

        fprintf(f, "# vtk DataFile Version 3.1\n");
        fprintf(f, "ESTER model %s\n", argv[1]);
        fprintf(f, "ASCII\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(f, "POINTS %d FLOAT\n", n[0]*n[1]*n[2]);

        for (auto i=0; i<n[0]; i++) {
            double zeta = i/(double) (n[0]-1);
            zetas(i) = zeta;
            for (auto j=0; j<n[1]; j++) {
                double theta = M_PI * j/(double) (n[1]-1);
                double r = zeta*A.map.leg.eval_00(A.map.r.row(-1), theta)(0);
                if (i == 0)
                    thetas(j) = theta;
                for (auto k=0; k<n[2]; k++) {

                    double phi = 2*M_PI * k/(double) (n[2]-1);

                    double x = r*sin(theta)*cos(phi);
                    double y = r*sin(theta)*sin(phi);
                    double z = r*cos(theta);

                    fprintf(f, "%f %f %f\n", x, y, z);
                }
            }
        }

        A.map.gl.eval(A.T, zetas, Tr);
        A.map.leg.eval_00(A.T, thetas, Tt);

        fprintf(f, "\nCELLS %d %d\n",
                (n[0]-1)*(n[1]-1)*(n[2]-1),
                9*(n[0]-1)*(n[1]-1)*(n[2]-1));

        for (auto i=0; i<n[0]-1; i++) {
            for (auto j=0; j<n[1]-1; j++) {
                for (auto k=0; k<n[2]-1; k++) {
                    fprintf(f, "8 %d %d %d %d %d %d %d %d\n",
                            (i)*(n[1]*n[2]) + (j)*n[2] + k,
                            (i+1)*(n[1]*n[2]) + (j)*n[2] + k,
                            (i+1)*(n[1]*n[2]) + (j+1)*n[2] + k,
                            (i)*(n[1]*n[2]) + (j+1)*n[2] + k,
                            (i)*(n[1]*n[2]) + (j)*n[2] + k+1,
                            (i+1)*(n[1]*n[2]) + (j)*n[2] + k+1,
                            (i+1)*(n[1]*n[2]) + (j+1)*n[2] + k+1,
                            (i)*(n[1]*n[2]) + (j+1)*n[2] + k+1
                           );
                }
            }
        }

        fprintf(f, "\nCELL_TYPES %d\n", (n[0]-1)*(n[1]-1)*(n[2]-1));
        for (auto i=0; i<n[0]-1; i++) {
            for (auto j=0; j<n[1]-1; j++) {
                for (auto k=0; k<n[2]-1; k++) {
                    fprintf(f, "12\n");
                }
            }
        }


        std::map<std::string, matrix> exportedFields;
        exportedFields["T"] = A.T*A.Tc;
        exportedFields["rho"] = A.rho*A.rhoc;
        exportedFields["p"] = A.p*A.pc;
        exportedFields["eps"] = A.nuc.eps;
        exportedFields["w"] = A.w;
        exportedFields["kappa"] = A.opa.k;
        exportedFields["Phi"] = A.phi;
        exportedFields["G"] = A.G;

        fprintf(f, "\nPOINT_DATA %d\n", n[0]*n[1]*n[2]);

        for (auto field: exportedFields) {
            fprintf(f, "SCALARS %s FLOAT\n", field.first.c_str());
            fprintf(f, "LOOKUP_TABLE DEFAULT\n");

            // interpolate fields on the new Cartesian grid
            matrix T = (Tr, field.second, Tt);

            for (auto i=0; i<n[0]; i++) {
                for (auto j=0; j<n[1]; j++) {
                    for (auto k=0; k<n[2]; k++) {
                        fprintf(f, "%f\n", T(i, j));
                    }
                }
            }
            fprintf(f, "\n");
        }

        fclose(f);
    }
    else {
        ester_err("Couldn't open model %s (is is a 2D model?)", input_model);
        return 1;
    }

    return 0;
}



