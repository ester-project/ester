
module mesa_opacity_module
  use mod_kind
  logical :: initialized = .false.
  logical :: OP_initialized = .false.
  logical :: lowT_initialized = .false.
  integer, parameter :: num_logTs=29, num_logRhos=71, num_logzs=15
  real(dp) :: logTs(num_logTs), logRhos(num_logRhos), logzs(num_logzs)
  real(dp), target :: f_ary(4*num_logRhos*num_logTs*num_logzs) ! for bicubic splines
  real(dp), pointer :: f(:,:,:,:)
  integer :: ilinx(num_logzs), iliny(num_logzs)

  integer , pointer, save :: ite(:),jne(:)
  real(dp), pointer, save :: lkap_ross_pcg(:,:), logT_pcg(:), logRho_pcg(:)
  real(dp), pointer, save :: lkap_pcg_lowT(:,:,:), logT_pcg_lowT(:),logR_pcg_lowT(:)

  contains
  
    subroutine load_table_lowT(filename, logT, logR, logkap,ierr)
      use mod_kind
      implicit none 

    character(len=256), intent(in) :: filename
  
    integer, intent(inout) :: ierr
    real(dp), pointer, intent(out) :: logT(:), logR(:), logkap(:,:)
  
    integer :: n, m, ke, ik
    CHARACTER(LEN=72) :: FMT !84
    integer :: nr, nt
    parameter(nr=46) !number of logR points.
    parameter(nt=105) !number of logT points A09
    !parameter(nt=75) !number of logT points GN93

    allocate(logT(nt),logR(nr),logkap(nt,nr),stat=ierr)
  
  
    FMT = '(3x,F5.3,46(1x,F7.3))'
  
    open(1, file = filename, &
    form = 'formatted', action ='read')
  
    do n= 1,5
      read(1,*)
    end do
  
    read(1,*) (logR(m), m=1,nr)
    read(1,*)
  
    do n =1,nt
      read(1,FMT) logT(n), (logkap(n,m), m=1,nr)
    enddo
    close(1)
  
   
    if (ierr /= 0) write(*,*) 'ERROR in load_table_lowT'
    ierr = 0
  end subroutine load_table_lowT

    subroutine init_highT_tables(Z, ierr)
      use mod_kind
      implicit none  
  
      real(dp), intent(in) :: Z
      integer, intent(inout) :: ierr

      real(dp) :: Xt(21)
      character(len=5) :: Xstr_om, Zstr

      character(len=256) :: filename

      integer :: i,n
      CHARACTER(LEN=38) :: FMT
      integer :: np
      parameter(np=1648)

      if (OP_initialized) return
    
      allocate(ite(np), jne(np), logT_pcg(np), logRho_pcg(np), lkap_ross_pcg(21,np),stat=ierr)
  
  
      FMT = '(i3,1x,i3,1x,F4.2,1x,F19.15,1x,F19.15)'

      Xt   = (/0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,&
      0.85, 0.90, 0.95, 1.00/)        


      write(Zstr,'(f5.3)')Z
      do i=1,21,1
         write(Xstr_om,'(f5.3)')Xt(i)
         filename = '/users/p0107/mombarg/ester/tables/op_mono/OP_lkap_A09_X' // Xstr_om // 'Z' // Zstr // '.txt'
         open(1, file = filename, &
         form = 'formatted', action ='read')
         do n =1,np
           read(1,FMT) ite(n), jne(n), logT_pcg(n), logRho_pcg(n), lkap_ross_pcg(i,n)
         enddo
         close(1)
      end do 

      write(*,*) 'OP tables initialized'
      OP_initialized = .true.

    end subroutine init_highT_tables  
 
    subroutine init_lowT_tables(Z, ierr)
      use mod_kind
      implicit none  

      real(dp), intent(in) :: Z
      integer, intent(inout) :: ierr

      character(len=256) :: filename
      character(len=:), allocatable :: Xstr
      
      integer :: n, m, ke, ik, i
      CHARACTER(LEN=72) :: FMT !84
      integer :: nr, nt
      parameter(nr=46) !number of logR points.
      parameter(nt=105) !number of logT points A09
      !parameter(nt=75) !number of logT points GN93

      if (lowT_initialized) return

      allocate(logT_pcg_lowT(nt),logR_pcg_lowT(nr),lkap_pcg_lowT(10,nt,nr),stat=ierr)
    
    
      FMT = '(3x,F5.3,46(1x,F7.3))' 
      
      do i=1,10,1
         select case (i)
         case(1)
            Xstr = '0.0'
         case(2)
            Xstr = '0.1'           
         case(3)
            Xstr = '0.2'
         case(4)
            Xstr = '0.35'         
         case(5)
            Xstr = '0.5'         
         case(6)
            Xstr = '0.7'         
         case(7)
            Xstr = '0.8'         
         case(8)
            Xstr = '0.9'         
         case(9)
            Xstr = '0.95'         
         case(10)
            Xstr = '0.98'         
         end select   
         filename = '/users/p0107/mombarg/ester_test/ester/tables/lowT_F05/lowT_fa05_a09p_z0.02_x' // Xstr // '.data'
         open(1, file = filename, &
         form = 'formatted', action ='read')
      
         do n= 1,5
         read(1,*)
         end do
      
         read(1,*) (logR_pcg_lowT(m), m=1,nr)
         read(1,*)
      
         do n =1,nt
         read(1,FMT) logT_pcg_lowT(n), (lkap_pcg_lowT(i,n,m), m=1,nr)
         enddo
         close(1)
      enddo
    
     
      if (ierr /= 0) write(*,*) 'ERROR in load_table_lowT'
      ierr = 0

      write(*,*) 'low-T tables initialized'
      lowT_initialized = .true.

    end subroutine init_lowT_tables  


  
    !!! Routine to compute Rossland opacity from a precomputed grid of the entire OP mono data.
    subroutine compute_kappa_fast(&
      fk, logT_cntr, logRho_cntr, &
      lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ix, ierr)!,&
      !lkap_ross_pcg, logT, logRho, ite, jne)
      ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
      
    use cubic_interpolator, only: interpolator
    use mod_kind
    implicit none 

    integer, intent(inout) :: ierr
    real(dp), intent(in) :: fk(:), logT_cntr, logRho_cntr
    integer, intent(in) :: ix

    integer :: nel, nptot
    parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
    real(dp), intent(out) :: lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho
  
    integer :: n, ke, nz, id, m, ik, i
  
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
  
    delta = sqrt((logRho_pcg - logRho_cntr)*(logRho_pcg - logRho_cntr)/0.25 +&
     (logT_pcg-logT_cntr)*(logT_pcg-logT_cntr)/0.025)
  
    delta_min_idx = MINLOC(delta, DIM=1)
    !delta_min(1)     = MINVAL(delta)(1)
    ite_min   = ite(delta_min_idx)
    jne_min   = jne(delta_min_idx)
    logT_min   = logT_pcg(delta_min_idx)
    logRho_min = logRho_pcg(delta_min_idx)
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
                logT_grid(ii,jj) = logT_pcg(i)
                logRho_grid(ii,jj) = logRho_pcg(i)
                i_grid(ii,jj) = i
                !write(*,*) ite_i, jne_i, ii, jj, i_grid(ii,jj), ',',logT_grid(ii,jj),&
                !',', logRho_grid(ii,jj)
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
        !write(*,*) 'warning interpolation quality in cell kappa_fast', k,ite_min,jne_min, &
        !logT_cntr, logRho_cntr, ii_min, jj_min
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
                 logT_grid(ii,jj) = logT_pcg(i)
                 logRho_grid(ii,jj) = logRho_pcg(i)
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
      !stop
      endif
    endif
    tries = tries + 1
    enddo
  
    do jj=1,4
      do ii=1,4
        ik = i_grid(ii,jj)
        lkap_Ross(ii,jj) = lkap_ross_pcg(ix,ik)
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
  
    !!! Routine to compute Rossland opacity from a precomputed grid of the entire OP mono data.
    subroutine compute_kappa_fast_lowT(&
      logT_cntr, logR_cntr, &
      lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnR, Xstr, ierr)
      ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
  
    use cubic_interpolator, only: interpolator
    use mod_kind
    implicit none

    integer, intent(inout) :: ierr
    real(dp), intent(in) :: logT_cntr, logR_cntr
    character(len=:), allocatable, intent(in) :: Xstr
    real(dp), intent(out) :: lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnR
  
    real(dp) :: lkap_Ross(4,4)
  
  
    integer :: ii,jj, ii_sel, jj_sel, ii_min, jj_min, ix
    !real(dp) :: logT_min   = 3.0 ! A09: 2.5
    !real(dp) :: delta_logT = 0.05
    real(dp) :: logR_min   = -8.0
    real(dp) :: delta_logR = 0.2
  
    type(interpolator) :: rossl_interpolator
  
    ierr = 0

    select case (Xstr)
    case('0.0')
      ix = 1
    case('0.1')
      ix = 2
    case('0.2')
      ix = 3
    case('0.35')
      ix = 4
    case('0.5')
      ix = 5
    case('0.7')
      ix = 6
    case('0.8')
      ix = 7
    case('0.9')
      ix = 8
    case('0.95')
      ix = 9
    case('0.98')
      ix = 10
   end select 

  
    ii_sel = 1
    do while(logT_pcg_lowT(ii_sel) .le. logT_cntr)
      ii_sel = ii_sel + 1
      !write(*,*) logT(ii_sel)
    end do
  

    ii_min = ii_sel - 2
    ii_min = MIN(ii_min,102) ! Can't center the 4x4 grid at the edge.

    jj_sel = NINT((logR_cntr - logR_min) / delta_logR) + 1


    jj_min = jj_sel - 2
    jj_min = MIN(jj_min, 42)

    do jj=1,4
      do ii=1,4
        lkap_Ross(ii,jj) = lkap_pcg_lowT(ix,ii+ii_min-1,jj+jj_min-1)
      enddo
    enddo
  
    !write(*,*) 'ii_min', ii_min, 'jj_min', jj_min
    !write(*,*) 'logT', logT(ii_min:ii_min+3)
    !write(*,*) 'logR', logR(jj_min:jj_min+3)
    !write(*,*) 'lkap_Ross', lkap_Ross(:,:)
  
    call rossl_interpolator% initialize()
    do jj = 1, 4
       do ii = 1, 4
          call rossl_interpolator% add_point(logT_pcg_lowT(ii+ii_min-1), logR_pcg_lowT(jj+jj_min-1) &
          , lkap_Ross(ii,jj))
       enddo
    enddo
  
    lkap_ross_cell  = rossl_interpolator% evaluate(logT_cntr,logR_cntr)
    dlnkap_rad_dlnT = rossl_interpolator% evaluate_deriv(logT_cntr, &
     logR_cntr, .true., .false.)
    dlnkap_rad_dlnR = rossl_interpolator% evaluate_deriv(logT_cntr, &
     logR_cntr, .false., .true.)
  
  
   end subroutine compute_kappa_fast_lowT
  
    subroutine do_blend_in_T(logT, logkap_lowT, logkap_highT, dlnkap_dlnT_highT, &
             dlnkap_dlnRho_highT, dlnkap_dlnT_lowT, dlnkap_dlnRho_lowT, &
             logkap, dlnkap_dlnT, dlnkap_dlnRho)
    use mod_kind
    implicit none             
  
    real(dp), intent(in) :: logT, logkap_lowT, logkap_highT
    real(dp), intent(in) :: dlnkap_dlnT_highT, dlnkap_dlnRho_highT, dlnkap_dlnT_lowT, dlnkap_dlnRho_lowT
    real(dp), intent(out) :: logkap, dlnkap_dlnT,dlnkap_dlnRho
  
    real(dp) :: F, S, d_S_dlnT, d_1mS_dlnT
    real(dp), parameter :: pi = 3.141592653589793239
    real(dp), parameter :: logT_l = 3.75
    real(dp), parameter :: logT_u = 4.50
  
    if (logT > logT_u) then
      logkap         = logkap_highT
      dlnkap_dlnT    = dlnkap_dlnT_highT
      dlnkap_dlnRho  = dlnkap_dlnRho_highT
    else if (logT < logT_l) then
      logkap = logkap_lowT
      dlnkap_dlnT    = dlnkap_dlnT_lowT
      dlnkap_dlnRho  = dlnkap_dlnRho_lowT
    else
    F = (logT - logT_l)/(logT_u - logT_l)
    S = 0.5*(1 - cos(F*pi))
    logkap = S*logkap_highT + (1-S)*logkap_lowT
    dlnkap_dlnRho = S*dlnkap_dlnRho_highT + (1-S)*dlnkap_dlnRho_lowT
    d_S_dlnT = pi/2d0*sin(F*pi)/(logT_u - logT_l) !removed minus sign
    d_1mS_dlnT = -d_S_dlnT
    dlnkap_dlnT = &
         S*dlnkap_dlnT_highT + (1-S)*dlnkap_dlnT_lowT + &
         d_S_dlnT*logkap_highT*log(10d0) + d_1mS_dlnT*logkap_lowT*log(10d0)
    endif
  
  end subroutine do_blend_in_T
  
   subroutine compute_kappa_conduction(logT, logRho, logK, dlogK_dlogT, dlogK_dlogRho, logkap, dlogkap_dlogT, dlogkap_dlogRho )
    use mod_kind
    implicit none 
     real(dp), intent(in) :: logT, logRho, logK, dlogK_dlogT, dlogK_dlogRho
     real(dp), intent(out) :: logkap, dlogkap_dlogT, dlogkap_dlogRho
     real(dp) :: boltz_sigma = 5.670374419184426e-05
  
  
     logkap = 3d0*logT - logRho - logK + log10(16d0 * boltz_sigma / 3d0)
  
     dlogkap_dlogRho = -1d0 - dlogK_dlogRho
     dlogkap_dlogT = 3d0 - dlogK_dlogT
  
   end subroutine compute_kappa_conduction
  
  
   subroutine combine_rad_with_conduction( &
         logRho, logT, &
         log_kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
         log_kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, &
         log_kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use mod_kind
         implicit none 

      real(dp), intent(in) :: logRho, logT
      real(dp), intent(in) :: log_kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, log_kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT
      real(dp), intent(out) :: log_kap, dlnkap_dlnRho, dlnkap_dlnT
      integer, intent(out) :: ierr ! 0 means AOK.
  
  
      real(dp) :: kap_rad, kap_ec, kap
      ierr = 0
  
      kap_rad = 10**log_kap_rad
      kap_ec = 10**log_kap_ec
      kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)
      log_kap = log10(kap)
  
      dlnkap_dlnRho = (kap/kap_rad) * dlnkap_rad_dlnRho + (kap/kap_ec) * dlnkap_ec_dlnRho
  
      dlnkap_dlnT = (kap/kap_rad) * dlnkap_rad_dlnT + (kap/kap_ec) * dlnkap_ec_dlnT
  
  
   end subroutine combine_rad_with_conduction

   subroutine init_potekhin(ierr)
    ! use kap_def, only: kap_dir
    !use interp_2d_lib_db, only: interp_mkbicub_db
    use mod_kind

    integer, intent(out) :: ierr

    character (len=256) :: filename
    integer :: read_err, iz, it, ir, shift
    integer :: ibcxmin                   ! bc flag for x=xmin
    real(dp) :: bcxmin(num_logTs)               ! bc data vs. y at x=xmin
    integer :: ibcxmax                   ! bc flag for x=xmax
    real(dp) :: bcxmax(num_logTs)               ! bc data vs. y at x=xmax
    integer :: ibcymin                   ! bc flag for y=ymin
    real(dp) :: bcymin(num_logRhos)               ! bc data vs. x at y=ymin
    integer :: ibcymax                   ! bc flag for y=ymax
    real(dp) :: bcymax(num_logRhos)               ! bc data vs. x at y=ymax
    real(dp) :: Z
    real(dp), pointer :: f1(:)

    !include 'formats'

    ierr = 0
    if (initialized) return

    shift = 4*num_logRhos*num_logTs
    f(1:4,1:num_logRhos,1:num_logTs,1:num_logzs)=>f_ary(1:shift*num_logzs)

    filename = '/gpfs/work/p0107/mombarg/mesa-r22.05.1/data/kap_data/condtabl.data'
    open(1,file=trim(filename),status='OLD',iostat=ierr)
    if (ierr /= 0) then
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,*) 'NOTICE: missing ' // trim(filename)
       write(*,*) 'Please remove the directory mesa/data/kap_data,'
       write(*,*) 'and rerun the mesa ./install script.'
       write(*,'(A)')
       !call mesa_error(__FILE__,__LINE__)
    end if
    !print*,'Reading thermal conductivity data...'
    read_err = 0
    read(1,'(A)',iostat=read_err) ! skip the first line
    if (read_err /= 0) ierr = read_err
    do iz = 1, num_logzs
       read (1,*,iostat=read_err)z, (logTs(it),it=1,num_logTs)
       if (read_err /= 0) ierr = read_err
       if (z .eq. 1d0) then
          logzs(iz) = 0d0
       else
          logzs(iz) = log10(z)
       endif
       do ir = 1, num_logRhos
          read(1,*,iostat=read_err) logRhos(ir), (f(1,ir,it,iz),it=1,num_logTs)
          if (read_err /= 0) ierr = read_err
       end do
    end do
    close(1)
    if (ierr /= 0) then
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,*) 'NOTICE: error trying to read ' // trim(filename)
       write(*,*) 'Please remove the directory mesa/data/kap_data,'
       write(*,*) 'and rerun the mesa ./install script.'
       write(*,'(A)')
       !call mesa_error(__FILE__,__LINE__)
    end if
    ! boundary condition is slope=0 at edges of tables
    ! this ensures continuous derivatives when we clip
    ibcxmin = 3; bcxmin(1:num_logTs) = 0d0
    ibcxmax = 3; bcxmax(1:num_logTs) = 0d0
    ibcymin = 3; bcymin(1:num_logRhos) = 0d0
    ibcymax = 3; bcymax(1:num_logRhos) = 0d0
    do iz = 1, num_logzs
       f1(1:shift) => f_ary(1+(iz-1)*shift:iz*shift)
       call interp_mkbicub_db( &
          logRhos, num_logRhos, logTs, num_logTs, f1, num_logRhos, &
          ibcxmin, bcxmin, ibcxmax, bcxmax, &
          ibcymin, bcymin, ibcymax, bcymax, &
          ilinx(iz), iliny(iz), ierr)
       if (ierr /= 0) then
          write(*,'(A)')
          write(*,'(A)')
          write(*,'(A)')
          write(*,'(A)')
          write(*,'(A)')
          write(*,*) 'NOTICE: error in ' // trim(filename)
          write(*,*) 'Please report the problem.'
          write(*,'(A)')
          !call mesa_error(__FILE__,__LINE__)
       end if
    end do
    initialized = .true.
    write(*,*) 'Thermal conductivity tables initialized'
 end subroutine init_potekhin

 subroutine interp_mkbicub_db(x,nx,y,ny,f1,nf2,&
         ibcxmin,bcxmin,ibcxmax,bcxmax,&
         ibcymin,bcymin,ibcymax,bcymax,&
         ilinx,iliny,ier)
         use mod_kind

    use bicub_db, only: do_mkbicub_db
    integer, intent(in) :: nx                        ! length of x vector
    integer, intent(in) :: ny                        ! length of y vector
    real(dp), intent(in) :: x(:) ! (nx)            ! x vector, strict ascending
    real(dp), intent(in) :: y(:) ! (ny)            ! y vector, strict ascending
    integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
    real(dp), intent(inout), pointer :: f1(:) ! =(4,nf2,ny)   ! data & spline coefficients
    integer, intent(in) :: ibcxmin                   ! bc flag for x=xmin
    real(dp), intent(in) :: bcxmin(:) ! (ny)       ! bc data vs. y at x=xmin
    integer, intent(in) :: ibcxmax                   ! bc flag for x=xmax
    real(dp), intent(in) :: bcxmax(:) ! (ny)       ! bc data vs. y at x=xmax
    integer, intent(in) :: ibcymin                   ! bc flag for y=ymin
    real(dp), intent(in) :: bcymin(:) ! (nx)       ! bc data vs. x at y=ymin
    integer, intent(in) :: ibcymax                   ! bc flag for y=ymax
    real(dp), intent(in) :: bcymax(:) ! (nx)       ! bc data vs. x at y=ymax
    integer, intent(out) :: ilinx                    ! =1: x grid is "nearly" equally spaced
    integer, intent(out) :: iliny                    ! =1: y grid is "nearly" equally spaced
    integer, intent(out) :: ier                      ! =0 on exit if there is no error.
    call do_mkbicub_db(x,nx,y,ny,f1,nf2,&
         ibcxmin,bcxmin,ibcxmax,bcxmax,&
         ibcymin,bcymin,ibcymax,bcymax,&
         ilinx,iliny,ier)
    end subroutine interp_mkbicub_db


     subroutine do_electron_conduction_potekhin( &
           zbar, logRho_in, logT_in, logkap, dlogkap_dlogRho, dlogkap_dlogT, ierr)

           use mod_kind

        real(dp) :: boltz_sigma = 5.670374419184426e-05
        real(dp), intent(in) :: zbar, logRho_in, logT_in
        real(dp), intent(out) :: logkap, dlogkap_dlogRho, dlogkap_dlogT
        integer, intent(out) :: ierr

        integer :: iz, iz1, iz2, shift
        real(dp) :: zlog, logRho, logT
        real(dp) :: alfa, beta, &
           logK1, kap1, dlogK1_dlogRho, dlogK1_dlogT, &
           logK2, kap2, dlogK2_dlogRho, dlogK2_dlogT, &
           logK, dlogK_dlogRho, dlogK_dlogT

        logical :: clipped_logRho, clipped_logT

        !include 'formats'

        ierr = 0
        shift = 4*num_logRhos*num_logTs

        if (logRho_in .lt. logRhos(1)) then
           logRho = logRhos(1)
           clipped_logRho = .true.
        else if (logRho_in .gt. logRhos(num_logRhos)) then
           logRho = logRhos(num_logRhos)
           clipped_logRho = .true.
        else
           logRho = logRho_in
           clipped_logRho = .false.
        end if

        if (logT_in .lt. logTs(1)) then
           logT = logTs(1)
           clipped_logT = .true.
        else if (logT_in .gt. logTs(num_logTs)) then
           logT = logTs(num_logTs)
           clipped_logT = .true.
        else
           logT = logT_in
           clipped_logT = .false.
        end if

        zlog = max(logzs(1),min(logzs(num_logzs),log10(max(1d-30,zbar))))
        !write(*,*) 'zlog', zlog, logzs(1), logzs(num_logzs)

        if (zlog <= logzs(1)) then ! use 1st
           call get1(1, logK, dlogK_dlogRho, dlogK_dlogT, ierr, shift, logRho, logT)
           logkap = logK
           return
        end if

        if (zlog >= logzs(num_logzs)) then ! use last
           call get1(num_logzs, logK, dlogK_dlogRho, dlogK_dlogT, ierr, shift, logRho, logT)
           logkap = logK
           return
        end if

        iz1 = -1
        do iz = 2, num_logzs
           if (zlog >= logzs(iz-1) .and. zlog <= logzs(iz)) then
              iz1 = iz-1; iz2 = iz; exit
           end if
        end do
        if (iz1 < 0) then
           write(*,*) 'num_logzs', num_logzs
           do iz = 1, num_logzs
              write(*,*) 'logzs(iz)', iz, logzs(iz)
           end do
           write(*,*) 'zlog', zlog
           write(*,*) 'confusion in do_electron_conduction'
           !call mesa_error(__FILE__,__LINE__)
        end if

        call get1(iz1, logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr, shift, logRho, logT)
        if (ierr /= 0) then
           write(*,*) 'interp failed for iz1 in do_electron_conduction', iz1, logRho, logT
           !call mesa_error(__FILE__,__LINE__)
        end if

        call get1(iz2, logK2, dlogK2_dlogRho, dlogK2_dlogT, ierr, shift, logRho, logT)
        if (ierr /= 0) then
           write(*,*) 'interp failed for iz2 in do_electron_conduction', iz2, logRho, logT
           !call mesa_error(__FILE__,__LINE__)
        end if

        ! linear interpolation in zlog
        alfa = (zlog - logzs(iz1)) / (logzs(iz2) - logzs(iz1))
        beta = 1d0-alfa
        logK = alfa*logK2 + beta*logK1
        if (clipped_logRho) then
           dlogK_dlogRho = 0
        else
           dlogK_dlogRho = alfa*dlogK2_dlogRho + beta*dlogK1_dlogRho
        end if
        if (clipped_logT) then
           dlogK_dlogT = 0
        else
           dlogK_dlogT = alfa*dlogK2_dlogT + beta*dlogK1_dlogT
        end if

        logkap = 3d0*logT_in - logRho_in - logK + log10(16d0 * boltz_sigma / 3d0)
        !write(*,*) 'ec', logT_in, logRho_in, logK

        !kap = exp10(logkap)
        dlogkap_dlogRho = -1d0 - dlogK_dlogRho
        dlogkap_dlogT = 3d0 - dlogK_dlogT

     end subroutine do_electron_conduction_potekhin


             subroutine get1(iz, logK, dlogK_dlogRho, dlogK_dlogT, ierr, shift, logRho, logT)
              use mod_kind

                !use kap_eval_support, only: Do_Kap_Interpolations
                integer, intent(in) :: iz, shift
                real(dp), intent(in) :: logRho, logT
                real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
                integer, intent(out) :: ierr
                logical, parameter :: dbg = .false.
                real(dp) :: fval, df_dx, df_dy, &
                   logRho0, logRho1, logT0, logT1
                integer :: i_logRho, j_logT, k

                real(dp), pointer :: f1(:)


                !include 'formats'
                ierr = 0
                f1(1:shift) => f_ary(1+(iz-1)*shift:iz*shift)

                if (logRho < logRhos(2)) then
                   i_logRho = 1
                else
                   i_logRho = num_logRhos-1
                   do k=2,num_logRhos-1
                      if (logRho >= logRhos(k) .and. logRho < logRhos(k+1)) then
                         i_logRho = k; exit
                      end if
                   end do
                end if
                logRho0 = logRhos(i_logRho)
                logRho1 = logRhos(i_logRho+1)

                if (logT < logTs(2)) then
                   j_logT = 1
                else
                   j_logT = num_logTs-1
                   do k=2,num_logTs-1
                      if (logT >= logTs(k) .and. logT < logTs(k+1)) then
                         j_logT = k; exit
                      end if
                   end do
                end if
                logT0 = logTs(j_logT)
                logT1 = logTs(j_logT+1)

                call Do_Kap_Interpolations( &
                   f1, num_logRhos, num_logTs, i_logRho, j_logT,&
                   logRho0, logRho, logRho1, logT0, logT, logT1, &
                   logK, dlogK_dlogRho, dlogK_dlogT)
                if (ierr /= 0) return

             end subroutine get1

             subroutine Do_Kap_Interpolations( &
                   fin1, nx, ny, i, j, x0, xget, x1, y0, yget, y1, fval, df_dx, df_dy)
                ! derived from routines in the PSPLINE package written by Doug McCune

                real(dp), dimension(:), pointer :: fin1 ! the spline data array, dimensions (4, nx, ny)
                integer, intent(in) :: nx, ny, i, j           ! target cell in the spline data
                real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
                real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
                real(dp), intent(out) :: fval, df_dx, df_dy

                real(dp), parameter :: z36th = 1d0 / 36d0
                real(dp), parameter :: one_sixth = 1d0 / 6d0

                real(dp), pointer :: fin(:,:,:)

                real(dp) :: xp, xpi, xp2, xpi2, cx, cxi, hx2, cxd, cxdi, hx, hxi
                real(dp) :: yp, ypi, yp2, ypi2, cy, cyi, hy2, cyd, cydi, hy, hyi

                logical, parameter :: dbg = .false.

                !include 'formats'

                fin(1:4,1:nx,1:ny) => fin1(1:4*nx*ny)

                hx=x1-x0
                hxi=1d0/hx
                hx2=hx*hx

                xp=(xget-x0)*hxi
                xpi=1d0-xp
                xp2=xp*xp
                xpi2=xpi*xpi

                cx=xp*(xp2-1d0)
                cxi=xpi*(xpi2-1d0)
                cxd=3d0*xp2-1d0
                cxdi=-3d0*xpi2+1d0

                hy=y1-y0
                hyi=1d0/hy
                hy2=hy*hy

                yp=(yget-y0)*hyi
                ypi=1d0-yp
                yp2=yp*yp
                ypi2=ypi*ypi

                cy=yp*(yp2-1d0)
                cyi=ypi*(ypi2-1d0)
                cyd=3d0*yp2-1d0
                cydi=-3d0*ypi2+1d0

                ! bicubic spline interpolation
                fval = &
                   xpi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1)) &
                   +xp*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)) &
                   +one_sixth*hx2*( &
                      cxi*(ypi*fin(2,i,j) +yp*fin(2,i,j+1))+ &
                      cx*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1))) &
                   +one_sixth*hy2*( &
                      xpi*(cyi*fin(3,i,j) +cy*fin(3,i,j+1))+ &
                      xp*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1))) &
                   +z36th*hx2*hy2*( &
                      cxi*(cyi*fin(4,i,j) +cy*fin(4,i,j+1))+ &
                      cx*(cyi*fin(4,i+1,j)+cy*fin(4,i+1,j+1)))

                ! derivatives of bicubic splines
                df_dx = &
                   hxi*( &
                      -(ypi*fin(1,i,j)  +yp*fin(1,i,j+1)) &
                      +(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1))) &
                   +one_sixth*hx*( &
                      cxdi*(ypi*fin(2,i,j) +yp*fin(2,i,j+1))+ &
                      cxd*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1))) &
                   +one_sixth*hxi*hy2*( &
                      -(cyi*fin(3,i,j)  +cy*fin(3,i,j+1)) &
                      +(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1))) &
                   +z36th*hx*hy2*( &
                      cxdi*(cyi*fin(4,i,j) +cy*fin(4,i,j+1))+ &
                      cxd*(cyi*fin(4,i+1,j)+cy*fin(4,i+1,j+1)))

                df_dy = &
                   hyi*( &
                      xpi*(-fin(1,i,j) +fin(1,i,j+1))+ &
                      xp*(-fin(1,i+1,j)+fin(1,i+1,j+1))) &
                   +one_sixth*hx2*hyi*( &
                      cxi*(-fin(2,i,j) +fin(2,i,j+1))+ &
                      cx*(-fin(2,i+1,j)+fin(2,i+1,j+1))) &
                   +one_sixth*hy*( &
                      xpi*(cydi*fin(3,i,j) +cyd*fin(3,i,j+1))+ &
                      xp*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1))) &
                   +z36th*hx2*hy*( &
                      cxi*(cydi*fin(4,i,j) +cyd*fin(4,i,j+1))+ &
                      cx*(cydi*fin(4,i+1,j)+cyd*fin(4,i+1,j+1)))

             end subroutine Do_Kap_Interpolations

!!!! =================================  Currently not used. ==================================================

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
   !real(dp) :: logT_pcg(1648), logRho_pcg(1648)
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
   open(1, file = '/gpfs/work/p0107/mombarg/OPmono_table/OP_mono_master_grid_MESA_emesh.txt', &
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
  

   subroutine create_OP_table_for_mixture(filename, ite, jne, logT, logRho, logkap)
      use mod_kind
      implicit none  

      character(len=256), intent(in) :: filename
      real(dp), pointer, intent(in) :: logT(:), logRho(:), logkap(:)
      integer, pointer, intent(in) :: ite(:), jne(:)
      integer :: np, i
      CHARACTER(LEN=38) :: FMT
   
      np = 1648 !SIZEOF(logT)
      FMT = '(i3,1x,i3,1x,F4.2,1x,F19.15,1x,F19.15)'
   
      !open(1, file = 'OP_lkap_A09_X710Z140.txt', status = 'new')
      open(1, file = filename, status = 'new')
      do i=1,np
         write(1,FMT) ite(i), jne(i), logT(i), logRho(i), logkap(i)
      enddo
      close(1)
   
      end subroutine create_OP_table_for_mixture
  end  module mesa_opacity_module
  