#include "ester-config.h"

        subroutine init_cesam_opa()
            use mod_opa
            use mod_donnees, only: langue, f_opa, nchim, nom_opa, &
                x0, y0, z0, ihe4, NOM_CHEMIN

            implicit none

            langue = 'PSE'
            NOM_CHEMIN = ESTER_DATADIR//'/tables/cesam/'
            nom_opa = 'opa_yveline.bin'
            nchim = 3

            print*, "Opacity table: " &
                & // trim(NOM_CHEMIN) // trim(nom_opa)

        end subroutine init_cesam_opa
