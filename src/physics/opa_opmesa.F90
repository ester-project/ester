#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

          subroutine opa_opmesa(xchim, t, ro, kap, dkapt, dkapro, dkapx)
            use mesa_opacity_module
            use mod_kind
            implicit none

            real(kind=dp), intent(in), dimension(3) :: xchim
            real(kind=dp), intent(in) :: t, ro
            real(kind=dp), intent(out) :: kap, dkapt, dkapro, dkapx

            logical :: loaded_op_master = .false.
            logical :: loaded_op_mix = .false.

            integer , pointer, save :: izz(:),ite1(:),jne1(:),ite2(:),jne2(:)
            real(kind=dp), pointer, save :: sig(:,:,:)
            real(kind=dp), pointer, save :: epatom(:,:),amamu(:),eumesh(:,:,:)
            real(kind=dp), pointer, save :: lkap_ross_pcg1(:), logT_pcg1(:), logRho_pcg1(:)
            real(kind=dp), pointer, save :: lkap_ross_pcg2(:), logT_pcg2(:), logRho_pcg2(:)

            real(dp), pointer :: logT_lowT1(:), logR_lowT1(:), logkap_lowT1(:,:)
            real(dp), pointer :: logT_lowT2(:), logR_lowT2(:), logkap_lowT2(:,:)
            character(len=256) :: filename1, filename2
            character(len=5) :: Xstr1_om, Xstr2_om, Zstr
            character(len=:), allocatable :: Xstr1, Xstr2
            
            
            integer :: ierr, i, ix
            real(dp) :: fk(17), A(17), eps(3:17), logT, logRho, logkap_highT, dlnkap_rad_dlnT_highT, dlnkap_rad_dlnRho_highT
            real(dp) :: logkap_highT1, dlnkap_rad_dlnT_highT1, dlnkap_rad_dlnRho_highT1
            real(dp) :: logkap_highT2, dlnkap_rad_dlnT_highT2, dlnkap_rad_dlnRho_highT2
            real(dp) :: Zsun, Z, X, X1, X2
            real(dp) ::  logkap_lowT, dlnkap_rad_dlnT_lowT, dlnkap_rad_dlnR_lowT
            real(dp) :: dlnkap_rad_dlnT, dlnkap_rad_dlnRho
            real(dp) ::  logR, lkap_ross_cell_lowT1, dlnkap_rad_dlnT_lowT1, dlnkap_rad_dlnR_lowT1
            real(dp) ::  alpha, lkap_ross_cell_lowT2, dlnkap_rad_dlnT_lowT2, dlnkap_rad_dlnR_lowT2
            real(dp) :: log_kap, kap_lowT, kap_lowT1, kap_lowT2, kap_highT, kap_highT1, kap_highT2
            real(dp) :: log_kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT
            real(dp) :: log_kap_rad, dlnkap_dlnRho, dlnkap_dlnT
            real(dp) :: logk, dlogK_dlogT, dlogK_dlogRho
            real(dp) :: Zi(17), zbar, Xt(21)
            logical :: use_mono
            ierr = 0

            
            ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
            
            !fk = (/0.9208, 0.0783, 10**-3.65, 10**-4.26, 10**-3.38, 10**-4.2, 10**-5.87, 10**-4.51, 10**-5.67, 10**-4.53, 10**-4.9,&
            ! 10**-5.86, 10**-5.73, 10**-6.4, 10**-6.65, 10**-4.59, 10**-5.18/)
            X             = xchim(1)
            Z             = xchim(3)
            logT          = log10(t)
            logRho        = log10(ro)
            use_mono = .false.
            
            if (use_mono .and. .not. loaded_op_master) then
             call load_op_master(izz,ite1,jne1,epatom,amamu,sig,eumesh,ierr)
             loaded_op_master = .true.
            end if

            
            !write(*,*) 'fk', fk
            
            Zsun = 0.0134 ! A09
            !Zsun = 0.0179 ! GN93
            eps = (/8.43, 7.83, 8.69, 7.93, 6.24, 7.60, 6.45, 7.51, 7.12, 6.40, 6.34, 5.64, 5.43, 7.50, 6.22 /) !A09
            !eps  = (/8.55, 7.97, 8.87, 8.08, 6.33, 7.58, 6.47, 7.55, 7.20, 6.52, 6.36, 5.67, 5.39, 7.50, 6.25 /) !GN93
            A   = (/1.008, 4.0026, 12.011, 14.007, 15.999, 20.180, 22.990, 24.305, 26.982, 28.085, 32.06, 39.948, 40.078, &
                    51.9961, 54.938043, 55.845, 58.693/)
            Zi   = (/1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, &
                    24, 25, 26, 28/)
            Xt   = (/0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,&
                     0.85, 0.90, 0.95, 1.00/)        
            do i=3,17
              fk(i) = 10**(eps(i) - 12)*(Z/Zsun) !/A(i)
            enddo
            
            
            do ix=1,1,1!1,11,1
              !X = (ix-1)*0.1
              fk(1) = X/A(1)
              fk(2) = (1-X-Z)/A(2)
              fk = fk /sum(fk)
              zbar = sum(fk*Zi)
              write(Zstr,'(f5.3)')Z
              do i=21,1,-1
                if(Xt(i) < X) exit 
              end do 
              write(Xstr1_om,'(f5.3)')Xt(i)
              write(Xstr2_om,'(f5.3)')Xt(i+1)

              filename1 = '/users/p0107/mombarg/ester/tables/op_mono/OP_lkap_A09_X' // Xstr1_om // 'Z' // Zstr // '.txt'
              filename2 = '/users/p0107/mombarg/ester/tables/op_mono/OP_lkap_A09_X' // Xstr2_om // 'Z' // Zstr // '.txt'
              if (use_mono) then
            
                call compute_kappa_grid(fk, lkap_ross_pcg1, logT_pcg1, logRho_pcg1, ierr, ite1,jne1,epatom,amamu,sig)
            

                write(*,*) X, Z, filename1
                if (.not. loaded_op_mix) then
                 call create_OP_table_for_mixture(filename1, ite1, jne1, logT_pcg1, logRho_pcg1, lkap_ross_pcg1)
                 loaded_op_mix = .true.
                endif
              else

                call load_OP_table_for_mixture(filename1, ite1, jne1, logT_pcg1, logRho_pcg1, lkap_ross_pcg1, ierr)
                call load_OP_table_for_mixture(filename2, ite2, jne2, logT_pcg2, logRho_pcg2, lkap_ross_pcg2, ierr)

                ! Tables available for X = 0.0, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98
                if (X .le. 0.1) then
                  Xstr1 = '0.0'
                  Xstr2 = '0.1'
                else if (X > 0.1 .and. X .le. 0.2) then
                  Xstr1 = '0.1'
                  Xstr2 = '0.2'
                else if (X > 0.2 .and. X .le. 0.35) then
                  Xstr1 = '0.2'
                  Xstr2 = '0.35'
                else if (X > 0.35 .and. X .le. 0.5) then
                  Xstr1 = '0.35'
                  Xstr2 = '0.5'  
                else if (X > 0.5 .and. X .le. 0.7) then
                  Xstr1 = '0.5'
                  Xstr2 = '0.7'
                else if (X > 0.7 .and. X .le. 0.8) then
                  Xstr1 = '0.7'
                  Xstr2 = '0.8'
                else if (X > 0.8 .and. X .le. 0.90) then
                  Xstr1 = '0.8'
                  Xstr2 = '0.9'
                else if (X > 0.90 .and. X .le. 0.95) then
                  Xstr1 = '0.9'
                  Xstr2 = '0.95'
                else if (X > 0.95 .and. X .le. 0.98) then
                  Xstr1 = '0.95'
                  Xstr2 = '0.98'
                else
                  write(*,*) 'X outside table ranges!'
                end if
            
                filename1 = '/users/p0107/mombarg/ester/tables/lowT_F05/lowT_fa05_a09p_z0.02_x' // Xstr1 // '.data'
                call load_table_lowT(filename1, logT_lowT1, logR_lowT1, logkap_lowT1,ierr)
                filename2 = '/users/p0107/mombarg/ester/tables/lowT_F05/lowT_fa05_a09p_z0.02_x'  // Xstr2 // '.data'
                call load_table_lowT(filename2, logT_lowT2, logR_lowT2, logkap_lowT2,ierr)
            
              endif
            enddo


            if (logT .ge. 3.75) then
              read(Xstr1_om,*)X1
              read(Xstr2_om,*)X2
              alpha = (X - X1)/(X2 - X1)

              call compute_kappa_fast(fk, logT, logRho, &
                logkap_highT1, dlnkap_rad_dlnT_highT1, dlnkap_rad_dlnRho_highT1, ierr,&
                lkap_ross_pcg1, logT_pcg1, logRho_pcg1, ite1, jne1)
                call compute_kappa_fast(fk, logT, logRho, &
                logkap_highT2, dlnkap_rad_dlnT_highT2, dlnkap_rad_dlnRho_highT2, ierr,&
                lkap_ross_pcg2, logT_pcg2, logRho_pcg2, ite2, jne2)

                kap_highT1 = 10**logkap_highT1
                kap_highT2 = 10**logkap_highT2
            
                kap_highT = (alpha)*kap_highT1 + (1-alpha)*kap_highT2
            
                dlnkap_rad_dlnT_highT   = (alpha*kap_highT1*dlnkap_rad_dlnT_highT1   + &
                                          (1-alpha)*kap_highT2*dlnkap_rad_dlnT_highT2)/kap_highT
                dlnkap_rad_dlnRho_highT = (alpha*kap_highT1*dlnkap_rad_dlnRho_highT1 + &
                                          (1-alpha)*kap_highT2*dlnkap_rad_dlnRho_highT2)/kap_highT
                logkap_highT = log10(kap_highT)


              if (ierr == 1) then
                logkap_highT = 1d99
                dlnkap_rad_dlnT_highT = 0d0
                dlnkap_rad_dlnRho_highT = 0d0
              endif  
            endif
            
            logR = logRho - 3*logT + 18d0
 

            if (logT .le. 4.5) then
                read(Xstr1,*)X1
                read(Xstr2,*)X2

                alpha = (X - X1)/(X2 - X1)

                call compute_kappa_fast_lowT(&
                    logT, logR, &
                    lkap_ross_cell_lowT1, dlnkap_rad_dlnT_lowT1, dlnkap_rad_dlnR_lowT1, ierr,&
                    logkap_lowT1, logT_lowT1, logR_lowT1)
            

                call compute_kappa_fast_lowT(&
                    logT, logR, &
                    lkap_ross_cell_lowT2, dlnkap_rad_dlnT_lowT2, dlnkap_rad_dlnR_lowT2, ierr,&
                    logkap_lowT2, logT_lowT2, logR_lowT2)
            
                kap_lowT1 = 10**lkap_ross_cell_lowT1
                kap_lowT2 = 10**lkap_ross_cell_lowT2
            
                kap_lowT = (alpha)*kap_lowT1 + (1-alpha)*kap_lowT2
            
                dlnkap_rad_dlnT_lowT = (alpha*kap_lowT1*dlnkap_rad_dlnT_lowT1 + &
                                       (1-alpha)*kap_lowT2*dlnkap_rad_dlnT_lowT2)/kap_lowT
                dlnkap_rad_dlnR_lowT = (alpha*kap_lowT1*dlnkap_rad_dlnR_lowT1 + &
                                       (1-alpha)*kap_lowT2*dlnkap_rad_dlnR_lowT2)/kap_lowT
                logkap_lowT = log10(kap_lowT)
            
            
            endif
           
           if (logT .ge. 3.75 .and. logT .le. 4.5  .and. logR .le. 1d0) then
             call do_blend_in_T(logT, logkap_lowT, logkap_highT,&
                       dlnkap_rad_dlnT_highT, dlnkap_rad_dlnRho_highT, dlnkap_rad_dlnT_lowT, dlnkap_rad_dlnR_lowT, &
                          log_kap_rad, dlnkap_rad_dlnT, dlnkap_rad_dlnRho)
           else if (logT < 3.75  .and. logR .le. 1d0) then
             log_kap_rad       = logkap_lowT
             dlnkap_rad_dlnT   =  dlnkap_rad_dlnT_lowT
             dlnkap_rad_dlnRho = dlnkap_rad_dlnR_lowT
           else if (logT > 4.5  .and. logR .le. 1d0) then
               log_kap_rad       = logkap_highT
               dlnkap_rad_dlnT   =  dlnkap_rad_dlnT_highT
               dlnkap_rad_dlnRho = dlnkap_rad_dlnRho_highT
           else
             log_kap_rad       = 1d99
             dlnkap_rad_dlnT   = 0d0
             dlnkap_rad_dlnRho = 0d0
         endif
          
          call init_potekhin(ierr)
          !write(*,*) 'ester', zbar, logRho, logT 
          call do_electron_conduction_potekhin( &
                zbar, logRho, logT, log_kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
          
          call combine_rad_with_conduction( &
                logRho, logT, &
                log_kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
                log_kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, &
                log_kap, dlnkap_dlnRho, dlnkap_dlnT, ierr) 
          !write(*,*) 'log_kap_rad/log_kap_cond', log_kap_rad/log_kap_ec
          !log_kap = log_kap_rad
          !dlnkap_dlnRho = dlnkap_rad_dlnRho  
          !dlnkap_dlnT = dlnkap_rad_dlnT    
          
          
          !write(*,*) 'logT', logT, 'logRho', logRho, 'log_kap_rad_highT', logkap_highT,&
          !'log_kap_rad_lowT', logkap_lowT, 'log_kap_T_blend', log_kap_rad, 'log_kap_cond', log_kap_ec,&
          !'log_kap_total', log_kap      
            
!            write(*,*) '**** high logT ****'
!            write(*,*) 'logT             ', logT
!            write(*,*) 'logRho           ', logRho
!            write(*,*) 'lkapp_Ross_rad   ', logkap_highT
!            write(*,*) 'dlnkap_rad_dlnT  ', dlnkap_rad_dlnT_highT
!            write(*,*) 'dlnkap_rad_dlnRho', dlnkap_rad_dlnRho_highT
!            write(*,'(A)')
!            
!            write(*,*) '**** low  logT ****'
!            write(*,*) 'logT             ', logT
!            write(*,*) 'logR             ', logR
!            write(*,*) 'lkapp_Ross_rad   ', logkap_lowT
!            write(*,*) 'dlnkap_rad_dlnT  ', dlnkap_rad_dlnT_lowT
!            write(*,*) 'dlnkap_rad_dlnR  ', dlnkap_rad_dlnR_lowT
!            write(*,'(A)')
!            
!            write(*,*) '**** T blend  ****'
!            write(*,*) 'log_kap_rad_blend', log_kap_rad
!            write(*,*) 'dlnkap_rad_dlnT  ', dlnkap_rad_dlnT
!            write(*,*) 'dlnkap_rad_dlnRho', dlnkap_rad_dlnRho
!            write(*,'(A)')
!            
!            write(*,*) '**** conduction  ****'
!            write(*,*) 'logK, dT, dRho  ', logK, dlogK_dlogT, dlogK_dlogRho
!            write(*,*) 'log_kap_ec      ', log_kap_ec
!            write(*,*) 'dlnkap_ec_dlnT  ', dlnkap_ec_dlnT
!            write(*,*) 'dlnkap_ec_dlnRho', dlnkap_ec_dlnRho
!            write(*,'(A)')
!            
!            write(*,*) '**** total combined  ****'
!            write(*,*) 'log_kap_rad_blend', log_kap
!            write(*,*) 'dlnkap_rad_dlnT  ', dlnkap_dlnT
!            write(*,*) 'dlnkap_rad_dlnRho', dlnkap_dlnRho

            kap    = log_kap
            dkapt  = dlnkap_dlnT !* (t/kap)
            dkapro = dlnkap_dlnRho !* (ro/kap)
            dkapx  = 0.0

            if (ISNAN(kap)) write(*,*) 'NaN in kap', logT, logRho


        end subroutine opa_opmesa          
