module opmono_modules            
contains
        subroutine load_op_master(iz,ite,jne,epatom,amamu,sig,eumesh,ierr)
            use mod_kind
            implicit none
  
  
            integer, intent(inout) :: ierr
            integer , pointer, intent(out) :: iz(:),ite(:),jne(:)
            real(dp), pointer, intent(out) :: sig(:,:,:)
            real(dp), pointer, intent(out):: epatom(:,:),amamu(:),eumesh(:,:,:)
          
            integer :: n, m, ke, ik
            CHARACTER(LEN=72) :: FMT !84
            integer :: nel, nptot, np
            parameter(nel = 17, nptot = 10000, np=1648) !number of elements and number of u-mesh points.
            real(dp), allocatable :: amamu_f(:,:)
            integer, allocatable  :: iz_f(:,:)
          
          
            allocate(iz_f(nel,np),iz(nel),ite(np),jne(np),stat=ierr)
            allocate(sig(nel,np,nptot),stat=ierr)
            allocate(epatom(nel,np),amamu_f(nel,np),amamu(nel),eumesh(nel,np,nptot),stat=ierr)
          
            FMT = '(i2,1x,i3,1x,i3,1x,F14.10,1x,F14.10,10000(1x,E12.6E3),10000(1x,E13.6E3))'
          
            write(*,*) 'Opening file...'
            open(1, file = 'OP_mono_master_grid_MESA_emesh.txt', &
            !open(1, file = '/gpfs/work/p0107/mombarg/OPmono_table/OP_mono_master_grid_MESA_emesh.txt', &
            form = 'formatted', action ='read')
            write(*,*) 'Loading OP mono data...'
            do ke =1, nel
              do n =1,np
              read(1,FMT)iz_f(ke,n),ite(n),jne(n),epatom(ke,n),amamu_f(ke,n),&
              (sig(ke,n,m), m=1,nptot),(eumesh(ke,n,m), m=1,nptot)
              enddo
            enddo
            close(1)
          
            do ke=1,nel
              amamu(ke) = amamu_f(ke,1)
              iz(ke)    = iz_f(ke,1)
            enddo
          
            write(*,*) 'OP mono data loaded.'
            ierr = 0
            end subroutine load_op_master

            subroutine load_OP_table_for_mixture(ite, jne, logT_pcg, logRho_pcg, lkap_ross_pcg, ierr)
              use mod_kind
              implicit none
              real :: X, Z
              integer, pointer, intent(out) :: ite(:),jne(:)
              real(dp), pointer, intent(out) :: logT_pcg(:), logRho_pcg(:), lkap_ross_pcg(:)
              integer, intent(inout) :: ierr
          
              integer :: n
              CHARACTER(LEN=38) :: FMT
              integer :: np
              parameter(np=1648)
          
          
              allocate(ite(np), jne(np), logT_pcg(np), logRho_pcg(np), lkap_ross_pcg(np),stat=ierr)
          
          
              FMT = '(i3,1x,i3,1x,F4.2,1x,F19.15,1x,F19.15)'
          
              write(*,*) 'Opening file...'
              !open(1, file = '/gpfs/work/p0107/mombarg/OPmono_table/OP_lkap_A09_X700Z200.txt', &
              !form = 'formatted', action ='read')              
              open(1, file = 'OP_lkap_A09_X700Z200.txt', &
              form = 'formatted', action ='read')
              write(*,*) 'Loading OP quick data...'
              do n =1,np
                read(1,FMT) ite(n), jne(n), logT_pcg(n), logRho_pcg(n), lkap_ross_pcg(n)
              enddo
              close(1)
          
            end subroutine load_OP_table_for_mixture

          subroutine compute_kappa_grid(fk, &
              lkap_ross_pcg, logT_pcg, logRho_pcg, ierr,&
              ite,jne,epatom,amamu,sig)
              ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
              
            use mod_kind
            implicit none
  
         
            integer, intent(inout) :: ierr
            real(dp), intent(in) :: fk(:)
            integer :: nel, nptot
            parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
            real(dp), pointer, intent(out) :: lkap_ross_pcg(:)!, fk_pcg(:)
            real(dp), pointer, intent(out) :: logT_pcg(:), logRho_pcg(:)
          
            integer , pointer :: ite(:),jne(:)
            real(dp), pointer :: sig(:,:,:)!, ak_f(:,:)
            real(dp), pointer :: epatom(:,:),amamu(:)
            integer :: imin, imax
          
            integer :: n, ke, nz, id, m, ik, i
          
            real(dp):: epa_mix_cell(1648), amu_mix_cell ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
            !integer ::  eid(nel)
            real(dp) :: log10_bohr_radius_sqr = -16.55280
            real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
            real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
            real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
            real(dp) :: mH
          
            !!!! For interpolator.
            integer ::  delta_min_idx
          
            real(dp) :: sig_Ross(1648) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
            real(dp), parameter :: pi = 3.141592653589793239
            real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light
          
            real, allocatable :: sig_mix_cell(:,:),sig_int(:)
            real(dp) :: log_amu_mix_cell
          
          
            ierr = 0
          
            !!! Compute an estimated temperature range.
            imin = 1 !-1
            imax = 1648 !-1
          
            !!! Compute the number of electrons per atom for the local mixture.
            epa_mix_cell = 0d0
            do i=imin,imax
                epa_mix_cell(i) = dot_product(fk,epatom(1:nel,i))
            enddo
          
            amu_mix_cell = dot_product(fk,amamu)
          
            mH = 1.00784 * 1.660538782d-24
            log_amu_mix_cell = log10(amu_mix_cell * mH)
            !logRho = 0.25*jne + log_amu_mix_cell - log10(epa_mix_cell)
          
            !!! Compute the monochromatic cross-section for the local mixture.
            allocate(sig_mix_cell(1648,nptot),sig_int(nptot-1), stat=ierr)
            sig_mix_cell = 0d0
          !$OMP PARALLEL DO PRIVATE(i,m) SCHEDULE(guided)
          
                do i=1,1648
                    do m=1,nptot
                    sig_mix_cell(i,m) = dot_product(fk,sig(:,i,m))
                    enddo
                enddo
          !$OMP END PARALLEL DO
          
            !!! Compute the Rossland mean cross-section by integrating over variable v (mesh equally spaced in v).
            sig_Ross = 0
          !$OMP PARALLEL DO PRIVATE(i,m,sig_int) SCHEDULE(guided)
              do i=1,1648
                    sig_int(1:nptot-1) = (1/sig_mix_cell(i,1:nptot-1) + 1/sig_mix_cell(i,2:nptot))/2. * dv !inv_sig_mix_cell(ii,jj,1:nptot-1) + inv_sig_mix_cell(ii,jj,2:nptot)
                    do m=1, nptot-1
                      sig_Ross(i) = sig_Ross(i) + sig_int(m)
                    enddo
              enddo
          !$OMP END PARALLEL DO
          
            deallocate(sig_mix_cell,sig_int)
          
            allocate(logT_pcg(1648),stat=ierr)
            logT_pcg   = 0.025*ite
            allocate(logRho_pcg(1648),stat=ierr)
            logRho_pcg = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) - logNa
            allocate(lkap_ross_pcg(1648),stat=ierr)
            lkap_ross_pcg =  log10_bohr_radius_sqr - log_amu_mix_cell - log10(sig_Ross)
          
            end subroutine compute_kappa_grid
          
          
            !!! Routine to compute Rossland opacity from a precomputed grid of the entire OP mono data.
            subroutine compute_kappa_fast(&
              fk, logT_cntr, logRho_cntr, &
              lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ierr,&
              lkap_ross_pcg, logT, logRho, ite, jne)
              ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
            use cubic_interpolator, only: interpolator
            use mod_kind
            implicit none
  
         
            integer, intent(inout) :: ierr
            real(dp), intent(in) :: fk(:), logT_cntr, logRho_cntr
            integer, pointer, intent(in) :: ite(:), jne(:)
            real(dp), pointer, intent(in) :: lkap_ross_pcg(:), logT(:), logRho(:)
            integer :: nel, nptot
            parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
            real(dp), intent(out) :: lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho
          
            integer :: n, ke, nz, id, m, ik, i
          
            !real(dp) :: fk(nel), fk_norm_fac !Local fractional abudance per element and normalization factor.
            real(dp):: epa_mix_cell(1648), amu_mix_cell ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
            real(dp) :: delta(1648)
            real(dp) :: log10_bohr_radius_sqr = -16.55280
            real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
            real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
            real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
            real(dp) :: mH
          
            !!!! For interpolator.
            integer ::  delta_min_idx
            real(dp) :: lkap_Ross(4,4),sig_Ross(4,4) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
            real(dp) :: sf, flux !'scale factor' for g_rad, local flux.
            real(dp), parameter :: pi = 3.141592653589793239
            real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light
          
            integer :: ii, jj, ite_min, jne_min, ii_min, jj_min, ite_step, jne_step
            integer :: ite_i, jne_i, dite, djne, i_grid(4,4)
            real(dp) :: logT_min, logRho_min, logT_grid(4,4), logRho_grid(4,4)
            integer ::  offset1, offset2, tries, missing_point(4,4)
            real(dp) :: log_amu_mix_cell
            integer :: imin, imax
            logical :: retry, do_difficult_point
          
            type(interpolator) :: rossl_interpolator
          
            ierr = 0
          
            imin = 1 !-1
            imax = 1648 !-1
          
          
            imin = 1 !-1
            imax = 1648 !-1
            ite_step = 2
            jne_step = 2
          
            !!! Select nearest points in logT,logRho for interpolation.
            !!! First, find the nearest OP data point, and check which of the four possible positionings this minimum has wrt to (T,Rho)_cell.
            !!! Acquire the remaining 15 points of the 4x4 grid, where (T,Rho)_cell is located in the inner square.
          
            delta = sqrt((logRho - logRho_cntr)*(logRho - logRho_cntr)/0.25 +&
             (logT-logT_cntr)*(logT-logT_cntr)/0.025)
          
            delta_min_idx = MINLOC(delta, DIM=1)
            !delta_min(1)     = MINVAL(delta)(1)
            ite_min   = ite(delta_min_idx)
            jne_min   = jne(delta_min_idx)
            logT_min   = logT(delta_min_idx)
            logRho_min = logRho(delta_min_idx)
            if (logT_min <= logT_cntr .and. logRho_min <= logRho_cntr) then
              ii_min = 2
              jj_min = 2
            else if (logT_min < logT_cntr .and. logRho_min > logRho_cntr) then
              ii_min = 2!3
              jj_min = 3!2
            else if (logT_min > logT_cntr .and. logRho_min < logRho_cntr) then
              ii_min = 3
              jj_min = 2
            else if (logT_min > logT_cntr .and. logRho_min > logRho_cntr) then
              ii_min = 3
              jj_min = 3
            endif
          
            offset1 = 0
            offset2 = 0
            missing_point = 1
            tries = 0
            retry = .true.
            !do while (point_not_found .ne. 0) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
            do while (retry) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
            missing_point = 1
            retry = .false.
            do jj=1,4
                  do ii=1,4
                    dite = (ii - ii_min + offset1)*ite_step !(ii - ii_min)*ite_step + offset2
                    djne = (jj - jj_min + offset2)*jne_step !(jj - jj_min)*jne_step + offset1
                   do i =imin,imax
                     ite_i = ite(i)
                     jne_i = jne(i)
          
                      if (ite_i .eq. ite_min +  dite .and. jne_i .eq. jne_min + djne) THEN
                        logT_grid(ii,jj) = logT(i)
                        logRho_grid(ii,jj) = logRho(i)
                        i_grid(ii,jj) = i
                        missing_point(ii,jj) = 0
  
                      endif
                    enddo
                  enddo
            enddo
          
            if (SUM(missing_point) > 0) then
              retry = .true.
              if (SUM(missing_point(2:4,1:3)) == 0) then
                offset1 = offset1 + 1
                offset2 = offset2 - 1
              else if (SUM(missing_point(1:3,1:3)) == 0) then
                offset1 = offset1 - 1
                offset2 = offset2 - 1
              else if (SUM(missing_point(2:4,2:4)) == 0) then
                offset1 = offset1 + 1
                offset2 = offset2 + 1
              else if (SUM(missing_point(1:3,2:4)) == 0) then
                offset1 = offset1 - 1
                offset2 = offset2 + 1
              else
                if (ii_min == 3 .and. jj_min == 2) then
                  offset1 = offset1 + 2
                  offset2 = offset2 - 2
                else if (ii_min == 2 .and. jj_min == 3) then
                  offset1 = offset1 - 2
                  offset2 = offset2 + 2
                else if (ii_min == 2 .and. jj_min == 2) then
                  offset1 = offset1 - 2
                  offset2 = offset2 - 2
                else if (ii_min == 3 .and. jj_min == 3) then
                  offset1 = offset1 + 2
                  offset2 = offset2 - 1
                endif
          
              endif
            endif
          
            if (tries > 2) THEN ! To prevent loop from getting stuck.
              do_difficult_point = .false.
              if (ite_min == 226 .and. jne_min == 90 .and. ii_min == 2 .and. jj_min == 3) do_difficult_point = .true.
              if (ite_min == 226 .and. jne_min == 88 .and. ii_min == 2 .and. jj_min == 2) do_difficult_point = .true.
              if (ite_min == 230 .and. jne_min == 88 .and. ii_min == 3 .and. jj_min == 2) do_difficult_point = .true.
              if (ite_min == 230 .and. jne_min == 90 .and. ii_min == 3 .and. jj_min == 3) do_difficult_point = .true.
          
              if (do_difficult_point) then
                 do jj=1,4
                   do ii =1,4
                     dite = (ii - ii_min)*2*ite_step
                     djne = (jj - jj_min - 1)*jne_step
                     do i=imin,imax
                       ite_i = ite(i)
                       jne_i = jne(i)
                       if (ite_i == ite_min + dite .and. jne_i == jne_min + djne) then
                         logT_grid(ii,jj) = logT(i)
                         logRho_grid(ii,jj) = logRho(i)
                         i_grid(ii,jj) = i
                       endif
                     enddo
                    enddo
                   enddo
                   retry = .false.
               else
              write(*,*) 'Cannot find points for interpolation compute_kappa_fast', ite_min, jne_min, logT_min, logRho_min, &
              logT_cntr, logRho_cntr, missing_point, ii_min, jj_min, offset1, offset2,&
              imin, imax
              ierr = 1
              return
              endif
            endif
            tries = tries + 1
            enddo
          
            do jj=1,4
              do ii=1,4
                ik = i_grid(ii,jj)
                lkap_Ross(ii,jj) = lkap_ross_pcg(ik)
              enddo
            enddo
          
            call rossl_interpolator% initialize()
            do jj = 1, 4
               do ii = 1, 4
                  call rossl_interpolator% add_point(logT_grid(ii,jj), logRho_grid(ii,jj) &
                  , lkap_Ross(ii,jj))
               enddo
            enddo
          
            lkap_ross_cell  = rossl_interpolator% evaluate(logT_cntr,logRho_cntr)
            dlnkap_rad_dlnT = rossl_interpolator% evaluate_deriv(logT_cntr, &
             logRho_cntr, .true., .false.)
            dlnkap_rad_dlnRho = rossl_interpolator% evaluate_deriv(logT_cntr, &
             logRho_cntr, .false., .true.)
          
          
            end subroutine compute_kappa_fast
end module opmono_modules            
