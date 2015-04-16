#include "ester-config.h"

    subroutine nuc_cesam_init

        USE mod_donnees, ONLY : lit_nl, nchim, nom_fich2, w_rot
        USE mod_kind
        USE mod_nuc, ONLY : nuc

        implicit none

        double precision, dimension(0, 0) :: jac
        double precision, dimension(0) :: comp, dcomp, ex
        double precision :: et, ero, hhe, be7e, b8e, n13e, o15e, f17e
        double precision, dimension(5) :: epsilon

        integer :: fait

        ! lecture du fichier de donnees
        nom_fich2 = ESTER_DATADIR//'/ester/tables/cesam_nuc/test'
        call lit_nl(w_rot)

        fait = 0
        call nuc(0d0, 0d0,  comp, dcomp, jac, .FALSE., fait, &
            epsilon, et, ero, ex, hhe, be7e, b8e, n13e, o15e, f17e)

    end

    subroutine nuc_cesam_init_abon(X, Z, comp)

        USE mod_donnees, ONLY : x0, y0, z0, nchim, nucleo
        USE mod_kind
        USE mod_nuc, ONLY : nuc

        implicit none

        double precision, dimension(0, 0) :: jac
        double precision, dimension(0) :: dcomp, ex
        double precision, dimension(nchim) :: comp
        double precision :: et, ero, hhe, be7e, b8e, n13e, o15e, f17e, X, Z
        double precision, dimension(5) :: epsilon

        integer :: i, fait

        x0 = X; z0 = Z; y0 = 1.-X-Z

        fait = 1

        call nuc(0d0, 0d0, comp, dcomp, jac, .FALSE., fait, &
            epsilon, et, ero, ex, hhe, be7e, b8e, n13e, o15e, f17e)

        comp = comp * nucleo

    end


    subroutine nuc_cesam_eps(t, ro, comp, epsilon, et, ero, ex)

        USE mod_donnees, ONLY : nucleo, nchim
        USE mod_kind
        USE mod_nuc, ONLY : nuc

        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(0, 0) :: jac
        DOUBLE PRECISION, DIMENSION(0) :: dcomp
        DOUBLE PRECISION, DIMENSION(nchim) :: comp, ex
        DOUBLE PRECISION :: t, ro, et, ero, hhe, be7e, b8e, n13e, o15e, f17e
        DOUBLE PRECISION, DIMENSION(5) :: epsilon
        INTEGER :: fait

        comp = comp/nucleo
        fait = 3
        CALL nuc(t, ro, comp, dcomp, jac, .TRUE., fait, &
            epsilon, et, ero, ex, hhe, be7e, b8e, n13e, o15e, f17e)
        comp = comp * nucleo
        ! print *, "epsilon:", epsilon
    end


    subroutine nuc_cesam_dcomp(t, ro, comp, dcomp, jac, pout)

        USE mod_donnees, ONLY : nucleo, nchim
        USE mod_kind
        USE mod_nuc, ONLY : nuc

        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(nchim, nchim) :: jac
        DOUBLE PRECISION, DIMENSION(nchim) :: dcomp
        DOUBLE PRECISION, DIMENSION(nchim) :: comp, ex
        DOUBLE PRECISION :: t, ro, et, ero, hhe, be7e
        DOUBLE PRECISION :: b8e, n13e, o15e, f17e
        DOUBLE PRECISION, DIMENSION(5) :: epsilon
        INTEGER :: fait, i, j
        INTEGER pout

        comp = comp/nucleo
        fait = 2
        CALL nuc(t, ro, comp, dcomp, jac, .TRUE., fait, &
            epsilon, et, ero, ex, hhe, be7e, b8e, n13e, o15e, f17e)
        comp = comp * nucleo
        dcomp = dcomp * nucleo
        do i = 1, nchim
            do j = 1, nchim
                jac(i, j) = nucleo(i) / nucleo(j) * jac(i, j)
            end do
        end do

    end

