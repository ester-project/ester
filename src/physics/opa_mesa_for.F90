#ifndef WITH_CMAKE
#include "ester-config.h"
#endif


!       integer :: species
!       integer, pointer, dimension(:) :: net_iso, chem_id
!       integer ::  kap_handle

        subroutine init_mesa_opa()
         use var_mesa_ester!, only:xa,kap_handle,
         use chem_def, only: chem_isos
         use chem_lib, only: chem_init
            use kap_lib
            use kap_lib
         use math_lib
         
            implicit none
         integer :: ierr
         ierr=0
         
        if (.not.init) then 
            call math_init() 
            call kap_init(.true., 'kap_cache', ierr)

        !inlist_fname : fichier avec les parametre d'entrée pour un RUN de MESA.
!                            print*,"init_mesa_opa kap_handle= ",kap_handle
            kap_handle = alloc_kap_handle_using_inlist( 'inlist', ierr)
!                            print*,"init_mesa_opa kap_handle= ",kap_handle
        
            species = 1000 !
!         net_iso = !NOT NEEDED in KAP get ==> no value should work.
        
        ! CHECK THE CONTENT OF 'isotopes.data'
            call chem_init('isotopes.data', ierr)
!                        print*,"toto c"
            allocate(chem_id(size(chem_isos% chem_id)))

            chem_id = chem_isos% chem_id
!                        print*,"toto d"

            if (.not.allocated(xa)) then
             allocate(xa(species))
            endif
            init=.true.
        endif
!                print*,"toto e"

        end subroutine init_mesa_opa

        subroutine opa_mesa(xchim, t, ro, kap, dlnkap_dlnT, dlnkap_dlnd, dkapx)
            use var_mesa_ester
            use mod_opa
            use mod_kind
!             use mod_donnees, only: nchim
            use eos_def, only: i_lnfree_e, i_eta
            use kap_lib
            use kap_def, only : num_kap_fracs  
            use math_lib
         use chem_def, only: chem_isos
            
            implicit none

            real(kind=8), intent(in), dimension(species) :: xchim
            real(kind=8), intent(in)  :: t, ro
            real(kind=8), intent(out) :: kap, dlnkap_dlnT, dlnkap_dlnd, dkapx
         integer :: ierr
         integer :: i
         character(len=4) :: ename
!             real(kind=8), dimension(nchim) :: xchim
            real(kind=8)   :: logRho, logT
            real(kind=8)   :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            real(kind=8)   :: eta, d_eta_dlnRho, d_eta_dlnT
            real(kind=8)   :: opacity,  dlnkap_dxa(species)
            real(8) :: kap_fracs(num_kap_fracs) ! pas d'utilité pour ESTER
!                        print*,"toto AX"

             
         do i=1,1000
         
         ename=chem_isos%name(chem_isos%chem_id(i))
!          print*,"ename=",ename
         if(ename == 'h1')   xa(i)=xchim(2)
         if(ename == 'he3')   xa(i)=xchim(3)
         if(ename == 'he4')   xa(i)=xchim(4)
         if(ename == 'c12')   xa(i)=xchim(5)
         if(ename == 'c13')   xa(i)=xchim(6)
         if(ename == 'n14')   xa(i)=xchim(7)
         if(ename == 'n15')   xa(i)=xchim(8)
         if(ename == 'o16')   xa(i)=xchim(9)
         if(ename == 'o17')   xa(i)=xchim(10)
         
         
         enddo
         
!          
!          xa(1)=xchim(1) !Mass fraction of H in the Sun.
!          xa(2)=xchim(2) !Mass fraction of He in the Sun.
!          xa(3)=xchim(3) !Mass fraction of H in the Sun.
!          xa(4)=xchim(4) !Mass fraction of He in the Sun.
!          xa(5)=xchim(5) !Mass fraction of He in the Sun.
!          xa(6)=xchim(6) !Mass fraction of H in the Sun.
!          xa(7)=xchim(7) !Mass fraction of He in the Sun.
!          xa(8)=xchim(8) !Mass fraction of H in the Sun.
!          xa(9)=xchim(9) !Mass fraction of He in the Sun.
!          xa(10)=xchim(10) !Mass fraction of He in the Sun.
!          xa(11)=xchim(11) !Mass fraction of H in the Sun.
!          xa(12)=xchim(12) !Mass fraction of He in the Sun.
!          xa(13)=xchim(13) !Mass fraction of H in the Sun.
!          xa(14)=xchim(14) !Mass fraction of He in the Sun.
!          xa(15)=xchim(15) !Mass fraction of He in the Sun.
!          xa(16)=xchim(16) !Mass fraction of H in the Sun.
!          xa(17)=xchim(17) !Mass fraction of He in the Sun.
!          xa(18)=xchim(18) !Mass fraction of H in the Sun.
!          xa(19)=xchim(19) !Mass fraction of He in the Sun.
!          xa(20)=xchim(20) !Mass fraction of He in the Sun.
!          xa(6)=xchim(6) !Mass fraction of He in the Sun.

!                        print*,"toto BX"
!                        print*,"xa(1)=",xa(1)
!                        print*,"xa(2)=",xa(2)
!                        print*,"xa(3)=",xa(3)
!                        print*,"xa(4)=",xa(4)
!                        print*,"xa(5)=",xa(5)
!                        print*,"xa(6)=",xa(6)


           logRho=log10(ro)
           logT=log10(t)
!                        print*,"logRho=",logRho
!                        print*,"logT=",logT
!                        print*,"toto CX"

!           lnfree_e = res(i_lnfree_e)
!           d_lnfree_e_dlnRho = d_dlnd(i_lnfree_e)
!           d_lnfree_e_dlnT = d_dlnT(i_lnfree_e)

!           eta = res(i_eta)
!           d_eta_dlnRho = d_dlnd(i_eta)
 !          d_eta_dlnT = d_dlnT(i_eta)
!                        print*,"toto DX"

           call kap_get( &
                kap_handle, species, chem_id, net_iso, xa, &
                logRho, logT, &
                lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                eta, d_eta_dlnRho, d_eta_dlnT, &
                kap_fracs, opacity, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)
!                        print*,"xa(1)=",xa(1)
!                        print*,"xa(2)=",xa(2)
!                        print*,"xa(3)=",xa(3)
!                        print*,"xa(4)=",xa(4)
!                        print*,"xa(5)=",xa(5)
!                        print*,"opacity=",opacity
!                        print*,"toto EX"

           kap    = opacity
!            dkapt  = dlnkap_dlnT*(kap/t)
!            dkapro = dlnkap_dlnd*(kap/ro)
           dkapx  = 0.0_8
           
!            print*,"t=",t
!            print*,"opacity=",opacity
!            print*,"kap=",kap
!            print*,"dkapro=",dkapro
!            print*,"dlnkap_dlnd=",dlnkap_dlnd
!            print*,"dkapt=",dkapt
!            print*,"dlnkap_dlnT=",dlnkap_dlnT

!            print*,"kap=",kap
!            print*,"dkapt=",dkapt
!            print*,"dkapro=",dkapro
!            print*,"dkapx=",dkapx
!            
!                        print*,"toto FX"


        end subroutine opa_mesa
