#include "ester-config.h"

      subroutine update_comp(T, rho, comp, r, nr, nchim)
          use  mod_evol
          implicit none

          integer, intent(in) :: nr, nchim
          double precision, intent(in), dimension(nr) :: T, rho, r
          double precision, intent(inout), dimension(nr, nchim) :: comp

          integer e

          ! print *, "** nr:", nr

          ! print *, "** r:", r(1), r(nr)
          ! print *, "** T:", T(1), T(nr)
          ! print *, "** rho:", rho(1), rho(nr)

          do e = 1, nchim
              ! print *, "** ", e, ": ", comp(1, e),  " ", comp(nr, e)
              comp(:,e) = comp(:,e) - comp(:,e)*T*0.01
          end do

      end subroutine update_comp

