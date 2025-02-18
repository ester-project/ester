#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

module path_module
  use iso_c_binding
  implicit none
  type, bind(c) :: PathData
    character(kind=c_char), dimension(256) :: full_path_cstr  ! C-compatible string as an array
    integer(c_int) :: length                                   ! C-compatible integer
  end type PathData
end module path_module

      	  subroutine CheckFile(filename,file_exists)
	    !!! checks the existence of a file, I think this is overcomplicating it, I really just need to use the INQUIRE line, at least for initial testing 
	    implicit none
	    CHARACTER(LEN=*), INTENT(IN) :: filename
	    LOGICAL, INTENT(OUT) :: file_exists

	    INQUIRE(file=filename, exist=file_exists)
	    !if (file_exists) then
	    !  write(*,*) filename, ' exists '
	    !else
	    !  write(*,*) filename, ' does not exist! '
	    !end if
    	  end subroutine CheckFile
    	  
    	  !contains

          !subroutine opa_opmesa(xchim, t, ro, kap, dkapt, dkapro, dkapx,abund,a_weights,abund_name,abund_name_len,full_path_cstr,full_path_cstr_len)
          
          !! original statement
          !subroutine opa_opmesa(xchim, t, ro, kap, dkapt, dkapro, dkapx,abund,a_weights,abund_name,abund_name_len)
          !! new statement including path 
	  !subroutine opa_opmesa(xchim, t, ro, kap, dkapt, dkapro, dkapx,abund,a_weights,abund_name,abund_name_len,full_path_cstr,full_path_cstr_len)        
          
          !subroutine opa_opmesa(xchim, t, ro, kap, dkapt, dkapro, dkapx,abund,a_weights,abund_name,abund_name_len,full_path_cstr,full_path_cstr_length)
          
          subroutine opa_opmesa(xchim, t, ro, kap, dkapt, dkapro, dkapx, abund, a_weights, abund_name, abund_name_len, path_data)
            !!! Routine to call subroutines for computing Rossland mean opacity tables and interpolating in X. Input and output consistent with ESTER.
            use mesa_opacity_module
            use mod_kind
            use, intrinsic :: iso_c_binding
            use path_module 
            
            implicit none

            
            integer, intent(in), value :: abund_name_len ! feels kind of dumb to pass the length of the name, but I'm trying to make the allocation flexible on the C++ end. 
            character(len=abund_name_len), intent(in) :: abund_name
            
            !integer :: full_path_cstr_len = 27
            
            !integer, intent(in), value :: full_path_cstr_length ! feels kind of dumb to pass the length of the name, but I'm trying to make the allocation flexible on the C++ end. 
            !character(len=256), intent(in) :: full_path_cstr
            
            !character(len=*) :: full_path_cstr  ! Deferred-length string
	    !integer :: full_path_cstr_length   ! Length passed from C++
            
            type(PathData), intent(inout) :: path_data  ! INTENT(INOUT) for modification
	    character(len=256) :: received_path
	    	
            !character(len=256) :: valid_path


            
            !integer, intent(in),value :: full_path_cstr_len
            !character(len=full_path_cstr_len)           
	    !character(len=256), intent(in) :: full_path_cstr            
	    !character(len=256) :: path_up_to_relative
	    !integer :: pos,pos_i
    
            real(kind=dp), intent(in), dimension(6) :: xchim
            real(kind=dp), intent(in), dimension(17) :: abund
            real(kind=dp), intent(in), dimension(17) :: a_weights

                       
            real(kind=dp), intent(in) :: t, ro
            real(kind=dp), intent(out) :: kap, dkapt, dkapro, dkapx

            logical :: loaded_op_master = .false.
            logical :: loaded_op_mix = .false.
            logical :: file_exists_1  
            logical :: file_exists_2  
            
            integer , pointer, save :: izz(:),ite1(:),jne1(:),ite2(:),jne2(:)
            real(kind=dp), pointer, save :: sig(:,:,:)
            real(kind=dp), pointer, save :: epatom(:,:),amamu(:),eumesh(:,:,:)
            real(kind=dp), pointer, save :: lkap_ross_pcg1(:), logT_pcg1(:), logRho_pcg1(:)
            real(kind=dp), pointer, save :: lkap_ross_pcg2(:), logT_pcg2(:), logRho_pcg2(:)

            real(dp), pointer :: logT_lowT1(:), logR_lowT1(:), logkap_lowT1(:,:)
            real(dp), pointer :: logT_lowT2(:), logR_lowT2(:), logkap_lowT2(:,:)
            !character(len=512) :: filename1, filename2
            character(len=1024) :: filename1, filename2,filename1_test,filename2_test
            !character(len=256), allocatable :: filename1_test,filename1_test_part_1,filename1_test_part_2
            
            INTEGER :: start_time, end_time, count_rate
	    REAL :: elapsed_time
            character(len=5) :: Xstr1_om, Xstr2_om, Zstr
            character(len=:), allocatable :: Xstr1, Xstr2

            integer :: ierr, i, ix
            real(dp) :: fk(17), A(17), eps(3:17), logT, logRho, logkap_highT, dlnkap_rad_dlnT_highT, dlnkap_rad_dlnRho_highT
            real(dp) :: fk_test(17)
            real(dp) :: logkap_highT1, dlnkap_rad_dlnT_highT1, dlnkap_rad_dlnRho_highT1
            real(dp) :: logkap_highT2, dlnkap_rad_dlnT_highT2, dlnkap_rad_dlnRho_highT2
            real(dp) :: Zsun, Z, X, X1, X2
            real(dp) :: dlnkap_rad_dlnT, dlnkap_rad_dlnRho
            real(dp) ::  logR, lkap_ross_cell_lowT1, dlnkap_rad_dlnT_lowT1, dlnkap_rad_dlnR_lowT1
            real(dp) ::  alpha, lkap_ross_cell_lowT2, dlnkap_rad_dlnT_lowT2, dlnkap_rad_dlnR_lowT2
            real(dp) :: log_kap, kap_lowT, kap_lowT1, kap_lowT2, kap_highT, kap_highT1, kap_highT2
            real(dp) :: log_kap_rad, dlnkap_dlnRho, dlnkap_dlnT
            real(dp) :: logk, dlogK_dlogT, dlogK_dlogRho
            real(dp) :: Zi(17), zbar, Xt(21)
            logical :: use_mono
            character(len=256) :: path_prefix,suffix, ext
	    !character(len=256) :: cmd, dir, fixed_path, full_path,cwd
	    !integer :: status
	    character(len=256) :: fixed_path, abs_path
	    character(len=20) :: abund_path_str
	    
	    ierr = 0
	    
	! Allocate the full_path_cstr string based on the passed length
    	!allocate(character(len=path_data%length) :: path_data%full_path_cstr)




	! Convert the C-style string into a Fortran string
	received_path = transfer(path_data%full_path_cstr, received_path)

	!print *, 'Received Full Path: ', trim(adjustl(received_path))
	!print *, 'Received Path Length (len_trim): ', len_trim(received_path)
	!print *, 'Received Path Length: ', path_data%length
	
	! Use path_data%length to extract only the valid part of the path
	received_path = adjustl(trim(received_path(1:path_data%length)))



	! Print the trimmed path and its actual length
	!print *, 'Valid Full Path: ', trim(adjustl(received_path))
	!print *, 'Valid Path Length: ', len_trim(received_path)
	    
	! Initialize path_up_to_relative
	!path_up_to_relative = ' '

	! Find the position of the relative path in the full_path_cstr
	!pos = index(full_path_cstr_test, fixed_path)

	!if (pos > 0) then
	! Extract the portion up to and including the relative path
	!path_up_to_relative = full_path_cstr_test(1:pos + len(fixed_path) - 1)

	! Print the cleaned path
	!print *, 'F90 Cleaned path is: ', path_up_to_relative
	!else
	!print *, 'Relative path not found in full_path_cstr.'
	!end if

            ! Trim the string to remove trailing spaces
            !full_path_trimmed = adjustl(trim(full_path_cstr))
	    
            !ierr = 0

            ! OP mono data for: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.

            X             = xchim(1)
            Z             = xchim(3)
            logT          = log10(t)
            logRho        = log10(ro)
            use_mono = .false. ! I think should be user defined outside of the script i.e. parser argument, but default is false i.e. we assume the grids exist 
            !use_mono = .true. 

	    Zsun = xchim(6)
            
            Zi   = (/1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, &
                    24, 25, 26, 28/)
            Xt   = (/0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,&
                     0.85, 0.90, 0.95, 1.00/) 
                     
            do i =3,17

            	fk(i) = 10**(abund(i) - 12)*(Z/Zsun)

            enddo
                     

            fk(1) = X/a_weights(1) ! fractional abundance for H
            fk(2) = (1-X-Z)/a_weights(2) ! fractional abundance for He, but calculated based on X and Z, not Y
            fk = fk /sum(fk) ! fractional abundance of the all elements 
            
            zbar = sum(fk*Zi) ! fractional abundances multiplied by atomic number? hmm what is this used for? --> Literally only exists here, must be something joey uses
            write(Zstr,'(f5.3)')Z
            do i=21,1,-1
                if(Xt(i) < X) exit
            end do
            write(Xstr1_om,'(f5.3)')Xt(i)
            write(Xstr2_om,'(f5.3)')Xt(i+1)            
			
            !!! obtain absolute path to where grids are saved
	    fixed_path = 'tables/op_mono'
	    !call absolute_path(fixed_path,abs_path) ! was blocked
	    !print *, 'F90 Abs path is: ', abs_path  
	    
    	    !print *, 'C++ Abs path is: ',path_data%full_path_cstr // "/" // fixed_path
    	    !print *, 'C++ Abs path_new is: ',TRIM(path_data%full_path_cstr(1:path_data%length)) // "/" // TRIM(fixed_path) 
    	    
    	    !abs_path = path_data%full_path_cstr // "/" // fixed_path
    	    
	    !abs_path = '/home/mgent/Documents/Ester/tables/op_mono'
	    !abs_path = '/home/mgent/Documents/Ester/'
	    !print *, 'F90 Abs path is: ', abs_path    

	    path_prefix = '/OP_lkap_'
	    !path_prefix = '/home/mgent/Documents/Ester/tables/op_mono/OP_lkap_' ! need to make this flexible to the user... 
	    suffix = '_X'
	    ext = '.txt'

            write(abund_path_str,'(A)') abund_name

	    !filename1 = TRIM(abs_path) // TRIM(path_prefix) // TRIM(abund_path_str) //&
	    ! TRIM(suffix) // TRIM(Xstr1_om) // TRIM(Zstr) // TRIM(ext)
	    !filename2 = TRIM(abs_path) // TRIM(path_prefix) // TRIM(abund_path_str) //&
	    ! TRIM(suffix) // TRIM(Xstr2_om) // TRIM(Zstr) // TRIM(ext)
	     
	     
	    ! Use received_path as part of the filename
	    filename1_test = trim(received_path) // '/' // trim(fixed_path) //&
	    trim(path_prefix) // trim(abund_path_str) // trim(suffix) //&
	    trim(Xstr1_om) // trim(Zstr) // trim(ext)

	    filename2_test = trim(received_path) // '/' // trim(fixed_path) //&
	    trim(path_prefix) // trim(abund_path_str) // trim(suffix) //&
	    trim(Xstr2_om) // trim(Zstr) // trim(ext)
	    
	    !print *, ''
	    !print *, 'Filename1: ', trim(filename1)
	    !print *, 'Generated Filename1: ', trim(filename1_test)
	    !print *, 'Filename2: ', trim(filename2)
	    !print *, 'Generated Filename2: ', trim(filename2_test)
	    
	    filename1 = filename1_test
	    filename2 = filename2_test

	    ! Print final filenames
	    !write(*,*) 'Filename1:', TRIM(filename1)
	    !write(*,*) 'Filename2:', TRIM(filename2)

    
              !! BOTH FILES NEED TO EXIST
            call CheckFile(filename1,file_exists_1)
            call CheckFile(filename2,file_exists_2)
            
            
              
            if (.not. file_exists_1 .or. .not. file_exists_2) then ! if either file does not exist, they need to be created, go to use_mono section
              	use_mono = .true.
              	loaded_op_mix = .false.
              	loaded_op_master = .false. ! not sure if these are necessary
            else
              	!write(*,*) 'both files exist! move to loading them'
              	use_mono = .false. ! this should be the default, thus unnecessary.   
            end if
              
              !!! time 
            CALL system_clock(count_rate=count_rate)

            if (use_mono) then
              
                ! if this was activated, its because either grid 1 or grid 2 doesn't exist 
                
                 ! check grid 2 --> Doing it first incase grid 1 has to be the last one loaded (idky)
                 
                if (.not. file_exists_2) then
                        !write(*,*) X, Z, Xstr2_om, Zstr, filename2
                        
                        ! Get the start time
	      		CALL system_clock(start_time)
	     	 	!call load_op_master(izz,ite2,jne2,epatom,amamu,sig,eumesh,ierr) 
	     	 	call load_op_master_binary(izz,ite2,jne2,epatom,amamu,sig,eumesh,ierr) 
	     	 	!call load_op_master_h5(izz,ite2,jne2,epatom,amamu,sig,eumesh,ierr) 
	     	        ! Get the end time
	                CALL system_clock(end_time)
	                ! Calculate the elapsed time in seconds
	                elapsed_time = REAL(end_time - start_time) / REAL(count_rate)
	                WRITE(*,*) elapsed_time, ' --- seconds'
	                
                	call compute_kappa_grid(fk, lkap_ross_pcg2, logT_pcg2, logRho_pcg2, ierr, ite2,jne2,epatom,amamu,sig)
      
                	if (.not. loaded_op_mix) then
		 		call create_OP_table_for_mixture(filename2, ite2, jne2, logT_pcg2, logRho_pcg2, lkap_ross_pcg2)
		 		
                 	end if
                 end if      
                
                ! check grid 1 
                
                if (.not. file_exists_1) then
                	
                	if (.not. file_exists_2) then
                		! this is if file2 was created just before i.e. load grid was done
                		deallocate(izz)
                		deallocate(epatom)
                		deallocate(amamu)
                		deallocate(sig)
                		deallocate(eumesh)
                	endif
                	
                	!! by blocking^ and adding deallocated lines to end of grid 2 block, when the grid2 was created, the code stopped after grid1, this is unpredictable!
                	
                        !write(*,*) X, Z, Xstr1_om, Zstr, filename1
                	!call load_op_master(izz,ite1,jne1,epatom,amamu,sig,eumesh,ierr) ! load master grid, is it applicable for all Xs? My tests showed otherwise
                	call load_op_master_binary(izz,ite1,jne1,epatom,amamu,sig,eumesh,ierr) ! load master grid, is it applicable for all Xs? My tests showed otherwise
                	!call load_op_master_h5(izz,ite1,jne1,epatom,amamu,sig,eumesh,ierr)                  	
	                call compute_kappa_grid(fk, lkap_ross_pcg1, logT_pcg1, logRho_pcg1, ierr, ite1,jne1,epatom,amamu,sig) ! create grid
                	!call compute_kappa_grid(fk, lkap_ross_pcg2, logT_pcg2, logRho_pcg2, ierr, ite2,jne2,epatom,amamu,sig) ! testing
                	loaded_op_master = .true.                	
                	if (.not. loaded_op_mix) then
                 		call create_OP_table_for_mixture(filename1, ite1, jne1, logT_pcg1, logRho_pcg1, lkap_ross_pcg1) ! table should have been created from grid
		 		!call create_OP_table_for_mixture(filename1, ite2, jne2, logT_pcg2, logRho_pcg2, lkap_ross_pcg2) ! testing             		
                 		loaded_op_mix = .true. 
                 		
                 		deallocate(izz) ! incase the grid needs to be created again 
                		deallocate(epatom)
                		deallocate(amamu)
                		deallocate(sig)
                		deallocate(eumesh)
                 		
                 	end if
                 end if
                 
                 use_mono = .false. ! now it should be set to false right? 
                 
                 
                 
        
            else
              
              	!write(*,*) 'loading two tables'
              	
              	!if (X > 0.7) then
              	!   write(*,*) filename1
              	!   write(*,*) filename2
              	!   write(*,*) ite1
              	!   write(*,*) ite2
              	!   write(*,*) jne1
              	!   write(*,*) jne2
              	!end if                 	   

                call load_OP_table_for_mixture(filename1, ite1, jne1, logT_pcg1, logRho_pcg1, lkap_ross_pcg1, ierr)
                call load_OP_table_for_mixture(filename2, ite2, jne2, logT_pcg2, logRho_pcg2, lkap_ross_pcg2, ierr)

                ! Tables available for X = 0.0, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98 in MESA. ! not here they are not 
                if (X .le. 0.1) then
                  Xstr1 = '0.0' ! will Xstr1_om not change naturlly? 
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

              endif

            !!! Interpolate in X between two precomputed tables. If Z changes, it should also be interpolated.
            
            !write(*,*) 'got to interpolation stage'
            if (logT .ge. 3.5) then
              !write(*,*) 'logT > 3.5'
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
                
                ! end of kappa table 2 info here (I think?) - MG
             
                logkap_highT = logkap_highT1
                dlnkap_rad_dlnT_highT = dlnkap_rad_dlnT_highT1
                dlnkap_rad_dlnRho_highT = dlnkap_rad_dlnRho_highT1

		! end of kappa table 1 info here (I think) - MG
		
		

              if (ierr == 1) then
                logkap_highT = 1d99
                dlnkap_rad_dlnT_highT = 0d0
                dlnkap_rad_dlnRho_highT = 0d0
              endif
            endif

            logR = logRho - 3*logT + 18d0
	    !write(*,*) 'logT is < 3.5' 


            if (logT .ge. 3.5 .and. logT .le. 8.0  .and. logR .le. 1d0) then
               !write(*,*) 'logT > 3.5 < 8 and logR < 1'

               log_kap_rad       = logkap_highT
               dlnkap_rad_dlnT   =  dlnkap_rad_dlnT_highT
               dlnkap_rad_dlnRho = dlnkap_rad_dlnRho_highT
            else
             write(*,*) 'Cell outside OP mono table!'
             log_kap_rad       = 1d99
             dlnkap_rad_dlnT   = 0d0
             dlnkap_rad_dlnRho = 0d0
            endif
            
            !write(*,*) 'defining kap'

            kap    = log_kap_rad
            dkapt  = dlnkap_rad_dlnT !* (t/kap)
            dkapro = dlnkap_rad_dlnRho !* (ro/kap)
            dkapx  = 0.0

            if (ISNAN(kap)) write(*,*) 'NaN in kap', logT, logRho
	    !endif ! added, was it missing? Nope it stops working some how! 
	    
	    !write(*,*) 'end of subroutine'

        end subroutine opa_opmesa
 
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


     
    
