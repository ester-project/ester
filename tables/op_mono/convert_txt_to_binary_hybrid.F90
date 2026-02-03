module mod_kind
    implicit none
    integer, parameter :: sp = selected_real_kind(6, 37)   ! Single precision
    integer, parameter :: dp = selected_real_kind(15, 307) ! Double precision
end module mod_kind

program convert_to_binary_v2

    !!! converts text file Mono opacity grid to binary so loading is faster when creating new grids 

    use mod_kind
    implicit none
    
    !integer :: n, m, ke, ik, unit_in, unit_out, nel, np, nptot
    !parameter(nel = 17, nptot = 10000, np = 1648)
    !integer, allocatable :: iz_f(:,:)
    !integer, allocatable :: ite(:), jne(:)
    !real(dp), allocatable :: epatom(:,:), amamu_f(:,:)
    !real(dp), allocatable :: sig(:,:,:), eumesh(:,:,:)
    !real(dp), allocatable :: amamu(:)
    !integer, allocatable  :: iz(:)
    !integer(8) :: size_collect
    !real(dp) :: size_total_text, percentage
    !CHARACTER(LEN=72) :: FMT
    
    !allocate(iz_f(nel, np))
    !allocate(ite(np), jne(np))
    !allocate(epatom(nel, np), amamu_f(nel, np))
    !allocate(sig(nel, np, nptot))
    !allocate(eumesh(nel, np, nptot))

    integer , allocatable :: iz(:),ite(:),jne(:)
    real(dp), allocatable :: sig(:,:,:)
    real(dp), allocatable :: epatom(:,:),amamu(:),eumesh(:,:,:)

    integer :: n, m, ke, ik, unit_in, unit_out
    CHARACTER(LEN=72) :: FMT !84
    integer :: nel, nptot, np
    parameter(nel = 17, nptot = 10000, np=1648) !number of elements and number of u-mesh points.
    real(dp), allocatable :: amamu_f(:,:)
    integer, allocatable  :: iz_f(:,:)
    integer(8) :: size_collect
    real(dp) :: size_total_text, percentage
    character(len=256) :: filename_inp
    
    size_total_text = 7565468656_dp ! bytes

    allocate(iz_f(nel,np),iz(nel),ite(np),jne(np))
    allocate(sig(nel,np,nptot))
    allocate(epatom(nel,np),amamu_f(nel,np),amamu(nel),eumesh(nel,np,nptot))
    
    write(*,*) 'arrays allocated dynamically, opening text file'

    FMT = '(i2,1x,i3,1x,i3,1x,F14.10,1x,F14.10,10000(1x,E12.6E3),10000(1x,E13.6E3))'

    write(*,*) 'Opening file...'
    
    ! Open the input and output files with proper unit numbers
    unit_in = 10
    unit_out = 20

    filename_inp = ''
    open(unit=unit_in, file=trim(filename_inp)//'OP_mono_master_grid_MESA_emesh.txt',&
     status='old', action='read', form='formatted')
    write(*,*) 'opening binary file to write in'    
    open(unit=unit_out, file='OP_mono_master_grid_MESA_emesh.bin', status='new', action='write', form='unformatted')
    
    write(*,*) 'Loading OP mono data...'
    
    size_collect = 0
    do ke =1, nel
      do n =1,np
       read(unit_in,FMT)iz_f(ke,n),ite(n),jne(n),epatom(ke,n),amamu_f(ke,n),&
       (sig(ke,n,m), m=1,nptot),(eumesh(ke,n,m), m=1,nptot)
      
       size_collect = size_collect + len(trim(FMT))
       ! Print the progress every 100 lines
       if (mod(n, 1000) == 0) then
                percentage = (real(size_collect) / size_total_text) * 100.0_dp
                write(*, '(I0, " bytes collected, ", F5.2, "% completed")') size_collect, percentage
       end if
      enddo
    enddo
    close(unit_in)
    
    write(*,*) 'doing amamu and iz loop'

    !do ke=1,nel
    !  amamu(ke) = amamu_f(ke,1)
    !  iz(ke)    = iz_f(ke,1)
    !enddo
    
    !deallocate(iz_f, amamu_f)

    write(*,*) 'OP mono data loaded. Writing to binary'
    !ierr = 0
    
    write(unit_out) iz_f, ite, jne, epatom, amamu_f, sig, eumesh
    !write(unit_out) iz, ite, jne, epatom, amamu, sig, eumesh
        
    close(unit_out)
    
    write(*,*) 'written to binary'
    
    !deallocate(ite, jne, epatom, sig, eumesh)
    deallocate(iz_f, amamu_f,ite, jne, epatom, sig, eumesh)
    
end program convert_to_binary_v2



