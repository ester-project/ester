#include "ester-config.h"
#include <cmath>
#include <string.h>
#include "matrix.h"
#include "constants.h"
#include "physics.h"

extern"C" {
    int init_cesam_opa_();
}

int opa_cesam(const matrix& X, double Z, const matrix& T, const matrix& rho,
		opa_struct& opa) {

    int error = 0;
    static bool init = false;

    if (!init) {
        init_cesam_opa_();
    }



	return error;
		
}
		
