#include "ester-config.h"

      subroutine init_evol(order)
          use mod_donnees, only: ordre, nom_nuc, langue
          implicit none
          integer, intent(in) :: order

          ordre = order
          nom_nuc = 'ppcno9'
          langue = 'english'
      end subroutine

      subroutine update_comp(t0, t1, rho0, rho1, comp, r, &
              nconv, dt)
          use mod_evol, only: rk_imps
          use mod_donnees, only: nchim, nucleo
          use mod_kind
          implicit none

          integer, intent(in) :: nconv
          real(kind=dp), intent(in) :: dt
          real(kind=dp), intent(in), dimension(nconv) :: t0, t1
          real(kind=dp), intent(in), dimension(nconv) :: rho0, rho1
          real(kind=dp), intent(in), dimension(nconv) :: r
          real(kind=dp), intent(inout), dimension(nchim, nconv) :: comp

          logical ret
          integer i
          real(kind=dp), dimension(nchim) :: comp1, esti
          real(kind=dp), dimension(nconv) :: dm

          comp1 = 0d0
          dm(:) = r(:)**2*rho0(:) !/ntheta

          do i=1, nconv
              comp(:,i) = comp(:,i) / nucleo
          enddo

          call rk_imps(t0(:), rho0(:), comp, &
                       t1(:), rho1(:), comp1,  &
                       dt, esti, ret, nconv, dm, .FALSE.)

          do i=1, nconv
              ! comp(:,i) = comp(:,i) * nucleo
              comp(:,i) = comp1(:) * nucleo / sum(dm)
          enddo

      end subroutine update_comp
