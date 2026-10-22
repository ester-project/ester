#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <string.h>
#include "matrix.h"
#include "constants.h"
#include "physics.h"
#include <iostream>
//#ifdef __cplusplus
extern "C" {
//#endif
#include "c_kap_interface_lib.h"
#include "c_kap_interface_external.h"
#include "c_kap_interface_def.h"
#include "c_utils_interface_lib.h"
#include "c_chem_interface_lib.h"
#include "c_chem_interface_def.h"
//#ifdef __cplusplus
}
//#endif

int opa_mesa_wrapper(const composition_map &comp,const matrix& X, double Z, const matrix& T, const matrix& rho,const eos_struct &eos,
		opa_struct& opa){

int error = 0;


int ierr = 0;

// MESA KAP must be initialized once and the same handle reused for kap_get.
static bool mesa_kap_initialized = false;
static int handle = -1;

char *mesa_dir = std::getenv("MESA_DIR");
if (mesa_dir == nullptr) {
    std::fprintf(
        stderr,
        "opa_mesa_wrapper: MESA_DIR is not set\n"
    );
    return 1;
}

char *libmesac_dir = std::getenv("LIBMESAC_ROOT");
if (libmesac_dir == nullptr) {
    std::fprintf(
        stderr,
        "opa_mesa_wrapper: LIBMESAC_ROOT is not set\n"
    );
    return 1;
}

// ------------------------------------------------------------
// Allocate/setup KAP handle
// ------------------------------------------------------------

char path_cache[1024];
snprintf(path_cache, sizeof(path_cache),
         "%s/data/kap_data/cache/", libmesac_dir);

char isotopes_file[] = "isotopes.data";

char path_init[1024];
snprintf(path_init, sizeof(path_init),
         "%s/install/share/kap_test_data/inlist_sample",
         libmesac_dir);

if (!mesa_kap_initialized) {

    mesa_kap_full_init(
        mesa_dir,
        isotopes_file,
        false,
        path_cache,
        &ierr);

    if (ierr != 0) {
        fprintf(stderr,
                "MESA KAP: mesa_kap_full_init failed, ierr=%d\n",
                ierr);
        return 1;
    }

    handle = kap_setup(path_init, &ierr);

    if (ierr != 0 || handle < 1) {
        fprintf(stderr,
                "ERROR: kap_setup failed, handle=%d ierr=%d\n",
                handle, ierr);
        return 1;
    }

    c_Kap_General_Info kap_info;

    kap_ptr(handle, &kap_info, &ierr);
    if (ierr != 0) {
        fprintf(stderr,
                "MESA KAP: failed to retrieve Kap_General_Info, ierr=%d\n",
                ierr);
        return 1;
    }

    kap_info.Zbase = opa.Zbase;

    kap_set_ptr(handle, &kap_info, &ierr);
    if (ierr != 0) {
        fprintf(stderr,
                "MESA KAP: failed to set Zbase=%g, ierr=%d\n",
                opa.Zbase, ierr);
        return 1;
    }

    mesa_kap_initialized = true;

}


const int maxpts = 2000; 

std::vector<int> chem_id_array;
std::vector<std::string> name_order;

const int num_chem_isos = get_num_chem_isos();

std::vector<int> net_iso(num_chem_isos, 0);

for (const auto &kv : comp) {

    const std::string &name = kv.first;

    // Skip aggregate and bookkeeping entries.
    if (name == "H" ||
        name == "X" ||
        name == "Y" ||
        name == "Z" ||
        name == "Xsol" ||
        name == "Ysol" ||
        name == "Zsol" ||
        name == "Ex")
        continue;

    const std::string mesa_name =
        mesa_nuclide_name(name);

    const int id =
        get_nuclide_index_from_name(mesa_name.c_str());

    if (id <= 0) {

        if (std::fabs(kv.second.data()[0]) > 1.0e-14) {
            std::fprintf(
                stderr,
                "opa_mesa_wrapper: unmapped isotope %s -> %s\n",
                name.c_str(),
                mesa_name.c_str()
            );
            return 1;
        }

        continue;
    }

    const int offset = id - 1;

    if (offset < 0 || offset >= num_chem_isos) {
        std::fprintf(
            stderr,
            "opa_mesa_wrapper: invalid chem_id=%d\n",
            id
        );
        return 1;
    }

    if (net_iso[offset] != 0) {
        std::fprintf(
            stderr,
            "opa_mesa_wrapper: duplicate chem_id=%d\n",
            id
        );
        return 1;
    }

    name_order.push_back(name);
    chem_id_array.push_back(id);

    net_iso[offset] =
        static_cast<int>(name_order.size());
}

const int NSpec_selected =
    static_cast<int>(chem_id_array.size());

if (NSpec_selected == 0) {
    std::fprintf(
        stderr,
        "opa_mesa_wrapper: no isotopes mapped\n"
    );
    return 1;
}

int num_kap_fracs = get_num_kap_fracs();
std::vector<double> kap_fracs(num_kap_fracs);
std::vector<double> dlnkap_dxa(NSpec_selected);
double kap, dlnkap_dlnRho, dlnkap_dlnT;


// Prepare opa matrices (flattened indexing like other opa_* functions)
const int Ncells = T.nrows() * T.ncols(); 
opa.k.dim(T.nrows(), T.ncols());
opa.xi.dim(T.nrows(), T.ncols());
opa.dlnxi_lnrho.dim(T.nrows(), T.ncols());
opa.dlnxi_lnT.dim(T.nrows(), T.ncols());


// Per-point loop: build xa, call kap_get, fill opa
std::vector<double> xa(NSpec_selected);

for (int ipt = 0; ipt < Ncells; ++ipt) { // flattened grid-cell index

    double lnfree_e = 0.0;
    double d_lnfree_e_dlnRho = 0.0;
    double d_lnfree_e_dlnT = 0.0;
    double eta = 0.0;
    double d_eta_dlnRho = 0.0;
    double d_eta_dlnT = 0.0;

    // fill xa for this point
    for (int j = 0; j < NSpec_selected; ++j) {
        const std::string &ename = name_order[j];
        xa[j] = comp.at(ename).data()[ipt];
    }

    double sum_xa = 0.0;

    for (double x : xa) {
        sum_xa += x;
    }

    if (!std::isfinite(sum_xa) ||
        std::fabs(sum_xa - 1.0) > 1.0e-8) {

        std::fprintf(
            stderr,
            "opa_mesa_wrapper: invalid composition sum "
            "at cell %d: %.16e\n",
            ipt,
            sum_xa
        );

        return 1;
    }

	double logT_val   = log10(T.data()[ipt]); 
	double logRho_val = log10(rho.data()[ipt]);

    if (!strcmp(eos.name, "mesa")) {

        lnfree_e =
            eos.lnfree_e.data()[ipt];

        d_lnfree_e_dlnRho =
            eos.d_lnfree_e_dlnRho.data()[ipt];

        d_lnfree_e_dlnT =
            eos.d_lnfree_e_dlnT.data()[ipt];

        eta =
            eos.eta.data()[ipt];

        d_eta_dlnRho =
            eos.d_eta_dlnRho.data()[ipt];

        d_eta_dlnT =
            eos.d_eta_dlnT.data()[ipt];
    }


	kap_get(handle, NSpec_selected,
			chem_id_array.data(),
			net_iso.data(),
			xa.data(),
			NSpec_selected,
			maxpts,
			logRho_val,
			logT_val,
			lnfree_e,
			d_lnfree_e_dlnRho,
			d_lnfree_e_dlnT,
			eta,
			d_eta_dlnRho,
			d_eta_dlnT,
			kap_fracs.data(),
			&kap,
			&dlnkap_dlnRho,
			&dlnkap_dlnT,
			dlnkap_dxa.data(),
			&ierr);

    if (ierr) {
        std::cerr << "kap_get returned ierr=" << ierr << " at ipt=" << ipt << std::endl;
        return 1;
    }

    opa.k.data()[ipt] = kap; 

    double Tval = T.data()[ipt];
    double rhoval = rho.data()[ipt];

    opa.xi.data()[ipt] = 16.0 * SIG_SB * pow(Tval, 3) / (3.0 * kap * rhoval);
    opa.dlnxi_lnrho.data()[ipt] = -1.0 - dlnkap_dlnRho;
    opa.dlnxi_lnT.data()[ipt] = 3.0 - dlnkap_dlnT;

}

return error;

}
