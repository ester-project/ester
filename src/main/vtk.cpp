#include "star.h"

#include <cmath>
#include <stdio.h>

void help() {
    fprintf(stderr, "Usage: ester vtk <input-model> -o <output-file>\n");
}

int main(int argc, char *argv[]) {

    if (argc !=  4) {
        help();
        return 1;
    }
    if (argv[2][0] != '-' && argv[2][1] != 'o') {
        help();
        return 1;
    }

    star2d A;
    if(A.read(argv[1]) == 0) {

        int n[3] = {16, 32, 32}; // Add equator and pole


        FILE * f = fopen(argv[3], "w");
        fprintf(f, "# vtk DataFile Version 3.1\n");
        fprintf(f, "ESTER model %s\n", argv[1]);
        fprintf(f, "ASCII\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(f, "POINTS %d FLOAT\n", n[0]*n[1]*n[2]);

        for (auto i=0; i<n[0]; i++) {
            double r = i/(double) (n[0]-1);
            for (auto j=0; j<n[1]; j++) {
                double theta = M_PI * j/(double) (n[1]-1);
                for (auto k=0; k<n[2]; k++) {
                    double phi = 2*M_PI * k/(double) (n[2]-1);

                    double x = r*sin(theta)*cos(phi);
                    double y = r*sin(theta)*sin(phi);
                    double z = r*cos(theta);
                    fprintf(f, "%f %f %f\n", x, y, z);
                }
            }
        }

        fprintf(f, "\nCELLS %d %d\n", (n[0]-1)*(n[1]-1)*(n[2]-1), 9*(n[0]-1)*(n[1]-1)*(n[2]-1));
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
        exportedFields["temperature"] = A.T;
        exportedFields["density"] = A.rho;
        exportedFields["eps"] = A.nuc.eps;
        exportedFields["w"] = A.w;

        fprintf(f, "\nPOINT_DATA %d\n", n[0]*n[1]*n[2]);

        for (auto field: exportedFields) {
            fprintf(f, "SCALARS %s FLOAT\n", field.first.c_str());
            fprintf(f, "LOOKUP_TABLE DEFAULT\n");

            for (auto i=0; i<n[0]; i++) {
                double r = i/(double) (n[0]-1);
                matrix T = A.map.gl.eval(field.second, r);
                for (auto j=0; j<n[1]; j++) {
                    double theta = M_PI * j/(double) (n[1]-1);
                    double t = A.map.leg.eval_00(T, theta)(0);
                    for (auto k=0; k<n[2]; k++) {
                        fprintf(f, "%f\n", t);
                    }
                }
            }
            fprintf(f, "\n");
        }

        fclose(f);
    }
    else {
        fprintf(stderr, "Couldn't open model %s\n", argv[1]);
        return 1;
    }

    return 0;
}



