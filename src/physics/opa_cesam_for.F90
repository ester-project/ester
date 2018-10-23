#include "ester-config.h"

        subroutine init_cesam_opa()
            use mod_opa
            use mod_donnees, only: langue, f_opa, nchim, nom_opa, &
                x0, y0, z0, ihe4, NOM_CHEMIN

            implicit none

            langue = 'PSE'
            NOM_CHEMIN = &
&ESTER_DATADIR//&
'/ester/tables/cesam/'

            nom_opa = 'opa_yveline'
            f_opa(1) = 'opa_yveline.data'
            nchim = 3
            z0 = 0.02d0
            ihe4 = 2

        end subroutine init_cesam_opa

        subroutine opa_cesam(xchim, t, ro, kap, dkapt, dkapro, dkapx)

            use mod_opa
            use mod_kind
            use mod_donnees, only: nchim
            implicit none

            real(kind=dp), intent(in), dimension(nchim) :: xchim
            real(kind=dp), intent(in) :: t, ro
            real(kind=dp), intent(out) :: kap, dkapt, dkapro, dkapx

            call opa(xchim, t, ro, kap, dkapt, dkapro, dkapx)

        end subroutine opa_cesam
