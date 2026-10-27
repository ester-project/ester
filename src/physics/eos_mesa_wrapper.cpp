#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "physics.h"
#include "constants.h"

extern "C" {
#include "c_eos_interface_lib.h"
#include "c_chem_interface_lib.h"
#include "c_chem_interface_def.h"
}

namespace {

struct MesaEosContext {
    bool initialized = false;
    bool mapping_ready = false;

    int handle = -1;
    int nres = 0;

    std::vector<std::string> name_order;
    std::vector<int> chem_id;
    std::vector<int> net_iso;

    std::vector<double> res;
    std::vector<double> d_dlnd;
    std::vector<double> d_dlnT;
};

MesaEosContext mesa_eos;


/* Initialise MESA EOS once and keep one EOS handle for the process lifetime. */
int ensure_mesa_eos_ready()
{
    if (mesa_eos.initialized)
        return 0;

    const char *mesa_dir_env = std::getenv("MESA_DIR");

    if (mesa_dir_env == nullptr) {
        std::fprintf(stderr, "eos_mesa_wrapper: MESA_DIR is not set\n");
        return 1;
    }

    char *mesa_dir = const_cast<char *>(mesa_dir_env);
    char isotopes[] = "isotopes.data";

    int ierr = 0;

    mesa_eos_full_init(
        mesa_dir,
        isotopes,
        true,
        &ierr
    );

    if (ierr != 0) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: mesa_eos_full_init failed, ierr=%d\n",
            ierr
        );
        return 1;
    }

    mesa_eos.handle = alloc_eos_handle(&ierr);

    if (ierr != 0 || mesa_eos.handle < 0) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: alloc_eos_handle failed, handle=%d ierr=%d\n",
            mesa_eos.handle,
            ierr
        );
        return 1;
    }

    mesa_eos.nres = get_num_eos_basic_results();

    if (mesa_eos.nres < 26) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: unexpected num_eos_basic_results=%d\n",
            mesa_eos.nres
        );
        return 1;
    }

    mesa_eos.res.assign(mesa_eos.nres, 0.0);
    mesa_eos.d_dlnd.assign(mesa_eos.nres, 0.0);
    mesa_eos.d_dlnT.assign(mesa_eos.nres, 0.0);

    mesa_eos.initialized = true;

    return 0;
}


// build the ESTER isotope -> MESA chem_id mapping once

int build_mesa_composition_mapping(const composition_map &comp)
{
    if (mesa_eos.mapping_ready)
        return 0;

    if (!mesa_eos.initialized) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: composition mapping requested before MESA init\n"
        );
        return 1;
    }

   /* c_chem_ids cid;
    get_some_isos(&cid);

    const std::map<std::string, int> key_to_id = {
        {"H1",   cid.ih1},
        {"H2",   cid.ih2},

        {"He3",  cid.ihe3},
        {"He4",  cid.ihe4},

        {"Li6",  cid.ili6},
        {"Li7",  cid.ili7},

        {"Be7",  cid.ibe7},
        {"Be9",  cid.ibe9},

        {"B8",   cid.ib8},
        {"B10",  cid.ib10},
        {"B11",  cid.ib11},

        {"C12",  cid.ic12},
        {"C13",  cid.ic13},

        {"N14",  cid.in14},
        {"N15",  cid.in15},

        {"O16",  cid.io16},
        {"O17",  cid.io17},
        {"O18",  cid.io18},

        {"F19",  cid.if19},

        {"Ne20", cid.ine20},
        {"Ne21", cid.ine21},
        {"Ne22", cid.ine22},

        {"Na22", cid.ina22},
        {"Na23", cid.ina23},

        {"Mg24", cid.img24},
        {"Mg25", cid.img25},
        {"Mg26", cid.img26},

        {"Al26", cid.ial26},
        {"Al27", cid.ial27},

        {"Si28", cid.isi28},
        {"Si29", cid.isi29},
        {"Si30", cid.isi30},

        {"P31",  cid.ip31},
        {"S32",  cid.is32},

        {"Ar36", cid.iar36},
        {"Ar38", cid.iar38},
        {"Ar40", cid.iar40},

        {"K39",  cid.ik39},
        {"K40",  cid.ik40},
        {"K41",  cid.ik41},

        {"Ca40", cid.ica40},
        {"Ca42", cid.ica42},
        {"Ca43", cid.ica43},
        {"Ca44", cid.ica44},
        {"Ca46", cid.ica46},
        {"Ca48", cid.ica48},

        {"Sc45", cid.isc45},

        {"Ti46", cid.iti46},
        {"Ti47", cid.iti47},
        {"Ti48", cid.iti48},
        {"Ti49", cid.iti49},
        {"Ti50", cid.iti50},

        {"V50",  cid.iv50},
        {"V51",  cid.iv51},

        {"Cr50", cid.icr50},
        {"Cr52", cid.icr52},
        {"Cr53", cid.icr53},
        {"Cr54", cid.icr54},

        {"Mn55", cid.imn55},

        {"Fe54", cid.ife54},
        {"Fe56", cid.ife56},
        {"Fe57", cid.ife57},
        {"Fe58", cid.ife58},

        {"Co59", cid.ico59},

        {"Ni58", cid.ini58},
        {"Ni60", cid.ini60},
        {"Ni61", cid.ini61},
        {"Ni62", cid.ini62},
        {"Ni64", cid.ini64},

        {"Cu63", cid.icu63},
        {"Cu65", cid.icu65},

        {"Zn64", cid.izn64},
        {"Zn66", cid.izn66},
        {"Zn67", cid.izn67},
        {"Zn68", cid.izn68},
        {"Zn70", cid.izn70}
    };*/

    mesa_eos.name_order.clear();
    mesa_eos.chem_id.clear();

    int unmapped_nonzero = 0;
    double unmapped_first_cell_mass = 0.0;

    for (const auto &kv : comp) {

        const std::string &name = kv.first;

        /* ESTER aggregate/bookkeeping entries: not MESA species. */
        if (name == "H" ||
            name == "X" ||
            name == "Y" ||
            name == "Z" ||
            name == "Xsol" ||
            name == "Ysol" ||
            name == "Zsol" ||
            name == "Ex")
            continue;

        /* Convert ESTER isotope name to MESA naming convention. */
        const std::string mesa_name = mesa_nuclide_name(name);

        /* Obtain the isotope ID directly from MESA chemistry. */
        const int id = get_nuclide_index_from_name(mesa_name.c_str());

        if (id > 0) {

            mesa_eos.name_order.push_back(name);
            mesa_eos.chem_id.push_back(id);

        } else {

            const double x0 = kv.second.data()[0];

            if (std::fabs(x0) > 1.0e-14) {

                ++unmapped_nonzero;
                unmapped_first_cell_mass += x0;

                if (unmapped_nonzero <= 20) {

                    std::printf(
                        "eos_mesa_wrapper: unmapped ESTER isotope "
                        "%s -> MESA %s, chem_id=%d, X=%.16e\n",
                        name.c_str(),
                        mesa_name.c_str(),
                        id,
                        x0
                    );
                }
            }
        }
    }

    if (mesa_eos.chem_id.empty()) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: no ESTER isotopes mapped to MESA\n"
        );
        return 1;
    }

    const int num_chem_isos = get_num_chem_isos();

    mesa_eos.net_iso.assign(num_chem_isos, 0);

    for (int i = 0;
         i < static_cast<int>(mesa_eos.chem_id.size());
         ++i) {

        const int offset = mesa_eos.chem_id[i] - 1;

        if (offset < 0 || offset >= num_chem_isos) {
            std::fprintf(
                stderr,
                "eos_mesa_wrapper: chem_id out of range: %d\n",
                mesa_eos.chem_id[i]
            );
            return 1;
        }

        if (mesa_eos.net_iso[offset] != 0) {
            std::fprintf(
                stderr,
                "eos_mesa_wrapper: duplicate MESA chem_id=%d "
                "for ESTER isotope %s\n",
                mesa_eos.chem_id[i],
                mesa_eos.name_order[i].c_str()
            );
            return 1;
        }

        mesa_eos.net_iso[offset] = i + 1;
    }

    mesa_eos.mapping_ready = true;

    if (unmapped_nonzero > 0) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: mapped %zu ESTER isotopes to MESA; "
            "unmapped non-zero first-cell isotopes=%d, mass=%.6e\n",
            mesa_eos.name_order.size(),
            unmapped_nonzero,
            unmapped_first_cell_mass
        );
    }

    return 0;
}

} // namespace


int eos_mesa_wrapper(
    const composition_map &comp,
    const matrix &X,
    double Z,
    const matrix &T,
    const matrix &p,
    matrix &rho,
    eos_struct &eos)
{
    /*
     * X and Z remain in the signature so this wrapper stays close to ESTER's
     * existing EOS interface.  MESA receives the isotope-resolved composition
     * from comp.
     */
    (void)X;
    (void)Z;

    if (ensure_mesa_eos_ready() != 0)
        return 1;

    if (build_mesa_composition_mapping(comp) != 0)
        return 1;

    const int nrows = T.nrows();
    const int ncols = T.ncols();
    const int ncells = nrows * ncols;

    if (p.nrows() != nrows || p.ncols() != ncols) {
        std::fprintf(
            stderr,
            "eos_mesa_wrapper: T and p dimensions differ: "
            "T=(%d,%d) p=(%d,%d)\n",
            nrows, ncols, p.nrows(), p.ncols()
        );
        return 1;
    }

    const int species =
        static_cast<int>(mesa_eos.name_order.size());

    rho.dim(nrows, ncols);

    eos.G1.dim(nrows, ncols);
    eos.cp.dim(nrows, ncols);
    eos.del_ad.dim(nrows, ncols);
    eos.G3_1.dim(nrows, ncols);
    eos.cv.dim(nrows, ncols);
    eos.d.dim(nrows, ncols);
    eos.prad.dim(nrows, ncols);
    eos.chi_rho.dim(nrows, ncols);
    eos.chi_T.dim(nrows, ncols);
    eos.s.dim(nrows, ncols);

    eos.lnfree_e.dim(nrows, ncols);
    eos.d_lnfree_e_dlnRho.dim(nrows, ncols);
    eos.d_lnfree_e_dlnT.dim(nrows, ncols);

    eos.eta.dim(nrows, ncols);
    eos.d_eta_dlnRho.dim(nrows, ncols);
    eos.d_eta_dlnT.dim(nrows, ncols);

    std::vector<double> xa(species, 0.0);

    for (int ipt = 0; ipt < ncells; ++ipt) {

        double sum_xa = 0.0;

        for (int j = 0; j < species; ++j) {

            const std::string &name = mesa_eos.name_order[j];

            const auto it = comp.find(name);

            if (it == comp.end()) {
                std::fprintf(
                    stderr,
                    "eos_mesa_wrapper: mapped isotope '%s' "
                    "is no longer present in comp\n",
                    name.c_str()
                );
                return 1;
            }

            const double x = it->second.data()[ipt];

            if (!std::isfinite(x) || x < 0.0) {
                std::fprintf(
                    stderr,
                    "eos_mesa_wrapper: invalid abundance "
                    "ipt=%d %s=%.16e\n",
                    ipt,
                    name.c_str(),
                    x
                );
                return 1;
            }

            xa[j] = x;
            sum_xa += x;
        }

        /*
          Do not silently hide missing composition. Roundoff-level differences are normalised away, a real missing
          mass fraction stops the test so the mapping can be completed
         */
        if (std::fabs(sum_xa - 1.0) > 1.0e-8) {
            std::fprintf(
                stderr,
                "eos_mesa_wrapper: mapped composition does not sum to 1 "
                "at ipt=%d: sum=%.16e missing=%.16e\n",
                ipt,
                sum_xa,
                1.0 - sum_xa
            );
            return 1;
        }

        if (sum_xa != 1.0) {
            for (double &x : xa)
                x /= sum_xa;
        }

        const double Tval = T.data()[ipt];
        const double Ptot = p.data()[ipt];

        if (!std::isfinite(Tval) || !std::isfinite(Ptot) ||
            Tval <= 0.0 || Ptot <= 0.0) {
            std::fprintf(
                stderr,
                "eos_mesa_wrapper: invalid ESTER state "
                "ipt=%d T=%.16e Ptot=%.16e\n",
                ipt,
                Tval,
                Ptot
            );
            return 1;
        }

        const double Prad =
            A_RAD * std::pow(Tval, 4.0) / 3.0;

        const double Pgas =
            Ptot - Prad;

        if (!std::isfinite(Pgas) || Pgas <= 0.0) {
            std::fprintf(
                stderr,
                "eos_mesa_wrapper: invalid gas pressure "
                "ipt=%d T=%.16e Ptot=%.16e Prad=%.16e Pgas=%.16e\n",
                ipt,
                Tval,
                Ptot,
                Prad,
                Pgas
            );
            return 1;
        }

        std::fill(
            mesa_eos.res.begin(),
            mesa_eos.res.end(),
            0.0
        );
        std::fill(
            mesa_eos.d_dlnd.begin(),
            mesa_eos.d_dlnd.end(),
            0.0
        );
        std::fill(
            mesa_eos.d_dlnT.begin(),
            mesa_eos.d_dlnT.end(),
            0.0
        );

        double rho_val = 0.0;
        double log10rho = 0.0;
        double dlnRho_dlnPgas_const_T = 0.0;
        double dlnRho_dlnT_const_Pgas = 0.0;

        int ierr = 0;

        eosPT_get(
            mesa_eos.handle,
            species,
            mesa_eos.chem_id.data(),
            mesa_eos.net_iso.data(),
            xa.data(),
            Pgas,
            std::log10(Pgas),
            Tval,
            std::log10(Tval),
            &rho_val,
            &log10rho,
            &dlnRho_dlnPgas_const_T,
            &dlnRho_dlnT_const_Pgas,
            mesa_eos.res.data(),
            mesa_eos.d_dlnd.data(),
            mesa_eos.d_dlnT.data(),
            &ierr
        );

        if (ierr != 0 || !std::isfinite(rho_val) || rho_val <= 0.0) {
            std::fprintf(
                stderr,
                "eos_mesa_wrapper: eosPT_get failed "
                "ipt=%d ierr=%d T=%.16e Pgas=%.16e rho=%.16e\n",
                ipt,
                ierr,
                Tval,
                Pgas,
                rho_val
            );
            return 1;
        }

        rho.data()[ipt] = rho_val;

        /*
         * MESA result indices are C indices:
         *   2  lnS
         *   6  grad_ad
         *   7  chiRho
         *   8  chiT
         *   9  Cp
         *   10 Cv
         *   14 Gamma1
         *   15 Gamma3
         */
        eos.s.data()[ipt]       = std::exp(mesa_eos.res[2]);
        eos.del_ad.data()[ipt]  = mesa_eos.res[6];
        eos.chi_rho.data()[ipt] = mesa_eos.res[7];
        eos.chi_T.data()[ipt]   = mesa_eos.res[8];
        eos.cp.data()[ipt]      = mesa_eos.res[9];
        eos.cv.data()[ipt]      = mesa_eos.res[10];
        eos.G1.data()[ipt]      = mesa_eos.res[14];
        eos.G3_1.data()[ipt]    = mesa_eos.res[15] - 1.0;

        eos.d.data()[ipt] =
            mesa_eos.res[8] / mesa_eos.res[7];

        eos.prad.data()[ipt] =
            Prad;

        // MESA EOS electron state for MESA KAP.
        // MESA Fortran indices are one-based; C arrays are zero-based.

        eos.lnfree_e.data()[ipt] =
            mesa_eos.res[4];

        eos.d_lnfree_e_dlnRho.data()[ipt] =
            mesa_eos.d_dlnd[4];

        eos.d_lnfree_e_dlnT.data()[ipt] =
            mesa_eos.d_dlnT[4];

        eos.eta.data()[ipt] =
            mesa_eos.res[5];

        eos.d_eta_dlnRho.data()[ipt] =
            mesa_eos.d_dlnd[5];

        eos.d_eta_dlnT.data()[ipt] =
            mesa_eos.d_dlnT[5];

    }

    return 0;
}