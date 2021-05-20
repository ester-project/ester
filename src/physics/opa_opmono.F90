#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

          subroutine opa_opmono(xchim, t, ro, kap, dkapt, dkapro, dkapx)
            use opmono_modules
            use mod_kind
            implicit none

            real(kind=dp), intent(in), dimension(3) :: xchim
            real(kind=dp), intent(in) :: t, ro
            real(kind=dp), intent(out) :: kap, dkapt, dkapro, dkapx

            real(kind=dp), save :: X_prev = -1
            logical :: loaded_op_master = .false.
            logical :: use_OPmono = .false.

            integer , pointer, save :: izz(:),ite(:),jne(:)
            real(kind=dp), pointer, save :: sig(:,:,:)
            real(kind=dp), pointer, save :: epatom(:,:),amamu(:),eumesh(:,:,:)
            real(kind=dp), pointer, save :: lkap_ross_pcg(:), logT_pcg(:), logRho_pcg(:)

            integer :: ierr, i
            real(dp) :: fk(17), A(17), eps(3:17), logT, logRho, lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, Zsun, Z, X
            ierr = 0
            Zsun = 0.0134
            eps = (/8.43, 7.83, 8.69, 7.93, 6.24, 7.60, 6.45, 7.51, 7.12, 6.40, 6.34, 5.64, 5.43, 7.50, 6.22 /)
            A   = (/1.008, 4.0026, 12.011, 14.007, 15.999, 20.180, 22.990, 24.305, 26.982, 28.085, 32.06, 39.948, 40.078, &
                    51.9961, 54.938043, 55.845, 58.693/)           
            ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
            X = xchim(1)
            Z = xchim(3)
            logT   = log10(t)
            logRho = log10(ro)
                        

            do i=3,17
              fk(i) = (10**(eps(i) - 12))*(Z/Zsun) !/A(i)
            enddo
            fk(1) = X/A(1)
            fk(2) = (1-X-Z)/A(2)
            fk = fk /sum(fk)
            if (use_OPmono) then
              if (.not. loaded_op_master) then
                call load_op_master(izz,ite,jne,epatom,amamu,sig,eumesh,ierr)
                loaded_op_master = .true.
              endif  
              if (X .ne. X_prev) then 
                write(*,*) 'compute_kappa_grid'
                call compute_kappa_grid(fk, lkap_ross_pcg, logT_pcg, logRho_pcg, ierr, ite,jne,epatom,amamu,sig)
                X_prev = X 
              endif
            else
              if (.not. loaded_op_master) then           
                call load_OP_table_for_mixture(ite, jne, logT_pcg, logRho_pcg, lkap_ross_pcg, ierr)
                loaded_op_master = .true.
              endif
            endif
            call compute_kappa_fast(fk, logT, logRho, &
              lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ierr,&
              lkap_ross_pcg, logT_pcg, logRho_pcg, ite, jne)

            kap    = lkap_ross_cell
            dkapt  = dlnkap_rad_dlnT !* (t/kap)
            dkapro = dlnkap_rad_dlnRho !* (ro/kap)
            dkapx  = 0.0


        end subroutine opa_opmono          
