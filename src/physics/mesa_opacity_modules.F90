module mesa_opacity_module
  use mod_kind

  logical :: initialized = .false.
  integer, parameter :: num_logTs=29, num_logRhos=71, num_logzs=15
  real(dp) :: logTs(num_logTs), logRhos(num_logRhos), logzs(num_logzs)
  real(dp), target :: f_ary(4*num_logRhos*num_logTs*num_logzs) ! for bicubic splines
  real(dp), pointer :: f(:,:,:,:)
  integer :: ilinx(num_logzs), iliny(num_logzs)

  contains
    
    subroutine load_op_master_binary(iz, ite, jne, epatom, amamu, sig, eumesh, ierr)
    use mod_kind
    implicit none

    integer, intent(inout) :: ierr
    integer, pointer, intent(out) :: iz(:), ite(:), jne(:)
    real(dp), pointer, intent(out) :: sig(:,:,:)
    real(dp), pointer, intent(out):: epatom(:,:), amamu(:), eumesh(:,:,:)

    integer, allocatable :: iz_f(:,:)
    real(dp), allocatable :: amamu_f(:,:)
    integer :: nel, nptot, np
    integer :: i, j    ! Declare loop variables
    character(len=256) :: fixed_path, abs_path,full_path
    
    parameter(nel = 17, nptot = 10000, np = 1648)
    
    !allocate(iz(nel), ite(np), jne(np), stat=ierr)
    allocate(iz(nel), iz_f(nel, np), ite(np), jne(np), stat=ierr)
    !if (ierr /= 0) return
    allocate(sig(nel, np, nptot), stat=ierr)
    !if (ierr /= 0) return
    allocate(epatom(nel, np), amamu_f(nel, np), amamu(nel), eumesh(nel, np, nptot), stat=ierr)
    !allocate(epatom(nel, np), amamu(nel), eumesh(nel, np, nptot), stat=ierr)
    !if (ierr /= 0) return

!    open(unit=1, file='/home/mgent/Documents/Ester/tables/op_mono/OP_mono_master_grid_MESA_emesh.bin', status='old', action='read', form='unformatted')
    !open(unit=1, file='OP_mono_master_grid_MESA_emesh.bin', status='old', action='read', form='unformatted')

    !!! obtain absolute path to where grids are saved
    fixed_path = 'tables/op_mono'
    call absolute_path(fixed_path,abs_path)
    full_path = TRIM(abs_path) // '/OP_mono_master_grid_MESA_emesh.bin'
    !write(*,*) 'full_path ',full_path
    write(*,*) 'Loading OP mono data...'
    !open(1, file = '../tables/op_mono/OP_mono_master_grid_MESA_emesh.bin', &
    !open(1, file = '/home/mgent/Documents/Ester/tables/op_mono/OP_mono_master_grid_MESA_emesh.bin', &
    open(1, file = full_path, &       
    form = 'unformatted', action ='read',status='old')
    read(1) iz_f, ite, jne, epatom, amamu_f, sig, eumesh
    !read(1) iz, ite, jne, epatom, amamu, sig, eumesh
    close(1)

    amamu = amamu_f(:, 1)
    deallocate(amamu_f)
    iz = iz_f(:, 1)
    deallocate(iz_f)
    
    write(*,*) 'OP mono data loaded.'

    ierr = 0
    
    end subroutine load_op_master_binary

    subroutine load_op_master(iz,ite,jne,epatom,amamu,sig,eumesh,ierr)
      !!! Routine to load in the large table with all OP monochromatic data.
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
    character(len=256) :: fixed_path, abs_path,full_path

    allocate(iz_f(nel,np),iz(nel),ite(np),jne(np),stat=ierr)
    allocate(sig(nel,np,nptot),stat=ierr)
    allocate(epatom(nel,np),amamu_f(nel,np),amamu(nel),eumesh(nel,np,nptot),stat=ierr)

    FMT = '(i2,1x,i3,1x,i3,1x,F14.10,1x,F14.10,10000(1x,E12.6E3),10000(1x,E13.6E3))'


    !!! obtain absolute path to where grids are saved
    fixed_path = 'tables/op_mono'
    call absolute_path(fixed_path,abs_path)
    full_path = TRIM(abs_path) // '/OP_mono_master_grid_MESA_emesh.txt'
    !write(*,*) 'full_path ',full_path
    write(*,*) 'Opening file...'
    !open(1, file = '/gpfs/work/p0107/mombarg/OPmono_table/OP_mono_master_grid_MESA_emesh.txt', &
    !open(1, file = '/home/mgent/Documents/Ester/tables/op_mono/OP_mono_master_grid_MESA_emesh.txt', &
    open(1, file = full_path, &
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
    
    ! iz_f, ite, jne, epatom, amamu_f, sig, eumesh
    end subroutine load_op_master


    subroutine create_OP_table_for_mixture(filename, ite, jne, logT, logRho, logkap)
      !write(*,*) 'grid being created'
      !!! Subroutine to save the precomputed Rossland opacities for a specific mixture.
      use mod_kind
      implicit none

      character(len=256), intent(in) :: filename
      real(dp), pointer, intent(in) :: logT(:), logRho(:), logkap(:)
      integer, pointer, intent(in) :: ite(:), jne(:)
      integer :: np, i
      CHARACTER(LEN=38) :: FMT

      np = 1648 !SIZEOF(logT)
      FMT = '(i3,1x,i3,1x,F4.2,1x,F19.15,1x,F19.15)'

      write(*,*) 'grid being created' ! new addition - Mg 

      open(1, file = filename, status = 'new')
      do i=1,np
         write(1,FMT) ite(i), jne(i), logT(i), logRho(i), logkap(i)
      enddo
      close(1)
      write(*,*) 'grid created' ! new addition - Mg 
    end subroutine create_OP_table_for_mixture

    subroutine load_OP_table_for_mixture(filename, ite, jne, logT_pcg, logRho_pcg, lkap_ross_pcg,ierr)
      !!! Subroutine to load in a precomputed table created by the create_OP_table_for_mixture routine.
      use mod_kind
      implicit none

      integer, intent(inout) :: ierr
      character(len=256), intent(in) :: filename
      integer, pointer, intent(out) :: ite(:),jne(:)
      real(dp), pointer, intent(out) :: logT_pcg(:), logRho_pcg(:), lkap_ross_pcg(:)

      real :: X, Z

      integer :: n
      CHARACTER(LEN=38) :: FMT
      integer :: np
      parameter(np=1648)

      allocate(ite(np), jne(np), logT_pcg(np), logRho_pcg(np), lkap_ross_pcg(np),stat=ierr)


      FMT = '(i3,1x,i3,1x,F4.2,1x,F19.15,1x,F19.15)'


      open(1, file = filename, &
      form = 'formatted', action ='read')
      do n =1,np
        read(1,FMT) ite(n), jne(n), logT_pcg(n), logRho_pcg(n), lkap_ross_pcg(n)
      enddo
      close(1)

      if (ierr /= 0) write(*,*) 'ERROR in load_OP_table_for_mixture'

    end subroutine load_OP_table_for_mixture

  subroutine compute_kappa_grid(fk, &
      lkap_ross_pcg, logT_pcg, logRho_pcg, ierr,&
      ite,jne,epatom,amamu,sig)
      !!! Routine to compute the Rossland mean opacity from the monochromatic ones for a given mixture defined by the fractional
      !!! abundances fk, where sum(fk) = 1. The table is spaced in log(T) and log(rho).
      ! OP mono data for: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
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


    subroutine compute_kappa_fast(&
      fk, logT_cntr, logRho_cntr, &
      lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ierr,&
      lkap_ross_pcg, logT, logRho, ite, jne)
      !!! Routine to interpolate Rossland opacity from a precomputed grid of the entire OP mono data.
      ! OP mono data for: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.

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
    !real(dp) :: logT(1648), logRho(1648)!, ite_grid(4,4), jne_grid(4,4)
    logical :: retry, do_difficult_point

    type(interpolator) :: rossl_interpolator

    ierr = 0

    imin = 1 !-1
    imax = 1648 !-1

    !!! Compute the number of electrons per atom for the local mixture.
  !  epa_mix_cell = 0d0
  !    do i=imin,imax
  !      epa_mix_cell(i) = dot_product(fk,epatom(:,i))
  !    enddo
  !
  !  !!! Compute the mean atomic mass for the local mixture.
  !  amu_mix_cell = dot_product(fk,amamu)


  !  logT   = 0.025*ite
  !  logRho = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) - logNa

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
                !write(*,*) ite_i, jne_i, ii, jj, i_grid(ii,jj), ',',logT_grid(ii,jj),&
                !',', logRho_grid(ii,jj)
                missing_point(ii,jj) = 0
                !ite_grid(ii,jj) = ite_i
                !jne_grid(ii,jj) = jne_i
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
      !stop
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

    subroutine absolute_path(relative_path, full_path)
        character(len=*), intent(in) :: relative_path
        character(len=256),intent(out) :: full_path
        character(len=256) :: temp_path
        logical :: found
        integer :: status,pos, len_e
    	character(len=256) :: cwd,cmd
    	
    	    !write(*,*) 'relative path ',relative_path
        
	    cmd = 'pwd > current_directory.txt'
	    ! Execute the command
	    call execute_command_line(cmd, wait=.true., exitstat=status)
	    if (status /= 0) then
		print *, "Failed to execute command."
		stop
	    endif
	    ! Open the temporary file to read the directory
	    open(unit=10, file='current_directory.txt', status='old', action='read')
	    read(10, '(A)') cwd
	    close(10)
	    ! Remove any trailing newline characters
	    cwd = adjustl(trim(cwd))
	    ! Print the current working directory
	    !print *, 'Current working directory is: ', cwd
	    ! Delete the temporary file
	    call execute_command_line('rm current_directory.txt')
	    ! Fixed path relative to the CWD

        ! Initialize the full_path
        full_path = ''
        

    ! Length of "Ester" including leading '/'
    len_e = len_trim('/Ester/')

    ! Check for "ester", "Ester", or "ESTER" in cwd
    pos = max(index(cwd, '/ester/'), max(index(cwd, '/Ester/'), index(cwd, '/ESTER/')))
    
    if (pos > 0) then
        ! Correct truncation: ensure the length includes "Ester" fully
        temp_path = cwd(1:pos + len_e - 2) ! very janky
        !temp_path = adjustl(trim(temp_path))
        
        full_path = trim(temp_path) // '/' // trim(relative_path)
    else
        print *, "Directory 'Ester' (case-insensitive) not found in the current working directory."
        stop
    endif



        !! Check for "ester", "Ester", or "ESTER" in cwd --> variations have been made by several users
        !found = .false.
        
        !!write(*,*) 'cwd: ', cwd

        !if (index(cwd, '/ester/') > 0) then
        !    found = .true.
        !    temp_path = trim(cwd) // '/' // relative_path
        !else if (index(cwd, '/Ester/') > 0) then
        !    found = .true.
        !    temp_path = trim(cwd) // '/' // relative_path
        !else if (index(cwd, '/ESTER/') > 0) then
        !    found = .true.
        !    temp_path = trim(cwd) // '/' // relative_path
        !endif

        !if (found) then
        !    full_path = trim(temp_path)
        !else
        !    print *, "Directory 'Ester' (case-insensitive) not found in the current working directory."
        !    stop
        !endif
    end subroutine absolute_path


  end  module mesa_opacity_module
