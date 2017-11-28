
      PROGRAM ascii2bin_opa

c     traduit une table d'opacite du binaire en ASCII
c     ou l'inverse (Yveline- octobre 1996)
c     adaptation F95 P.morel, juin 2001

c---------------------------------------------------------

      USE mod_kind
      USE mod_numerique, ONLY : pause

      IMPLICIT NONE

      REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: kap
      REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: T6, logR
      REAL (kind=dp) :: X, Z

      INTEGER :: nt, nr, nx, nz, ix, iz, it, ir

      LOGICAL :: ok

      CHARACTER (len=1) :: answer
      CHARACTER (len=40) :: nomasc, nombin

c-------------------------------------------------------------------

      PRINT*,'traduction ASCII ---> binaire (o/n?)'
      READ*, answer

      IF(answer == 'o')THEN
       PRINT*,'ASCII ---> binaire'
       PRINT*,'enter the name of the existing ASCII opacity table: opa_yveline.data'
       READ*, nomasc
       PRINT*,'enter the name of the binary table: opa_yveline.bin'
       READ*, nombin
       INQUIRE(file=trim(nombin),exist=ok)
       IF(ok)THEN
        PRINT*,'le fichier ',trim(nombin),' existe, a ecraser? o/n'
        READ*, answer
        IF(answer /= 'o') THEN
         PRINT*,'arret' ; STOP
        ENDIF
       ENDIF
       OPEN(unit=8,form='formatted',status='old',action='read',file=TRIM(nomasc))
       OPEN(unit=9,form='unformatted',status='unknown',file=TRIM(nombin))
       READ(8,1)nz,nx,nt,nr ; WRITE(9)nz,nx,nt,nr
c       WRITE(*,1)nz,nx,nt,nr ; CALL pause('en 1')
1      FORMAT(4i3)
       ALLOCATE(T6(nt),logR(nr),kap(nt,nr))
       READ(8,2)T6 ; WRITE(9)T6
2      FORMAT(3(1x,es23.15))      
       READ(8,2)logR ; WRITE(9)logR
       DO iz=1,nz
        READ(8,2)Z ; WRITE(9)Z
        DO ix=1,nx
         READ(8,2)X ; WRITE(9)X
         READ(8,2)((kap(it,ir),ir=1,nr),it=1,nt)
         WRITE(9)((kap(it,ir),ir=1,nr),it=1,nt)
        ENDDO
       ENDDO
 
      ELSE
       PRINT*,'binaire ---> ASCII' 
       PRINT*,'name of the existing binary table: '
       READ*,nombin
       PRINT*,'name of the ASCII opacity table: '
       READ*,nomasc
       INQUIRE(file=trim(nomasc),exist=ok)
       IF(ok)THEN
        PRINT*,'le fichier ',trim(nomasc),' existe, a ecraser? o/n'
        READ*, answer
        IF(answer /= 'o') THEN
         PRINT*,'arret' ; STOP
        ENDIF
       ENDIF
       OPEN(unit=8,form='formatted',status='unknown',file=nomasc)
       OPEN(unit=9,form='unformatted',status='old',action='read',file=nombin)
       READ(9)nz,nx,nt,nr ; WRITE(8,*)nz,nx,nt,nr
       ALLOCATE(T6(nt),logR(nr),kap(nt,nr))
       READ(9)T6 ; WRITE(8,*)T6
       READ(9)logR ; WRITE(8,*)logR
       DO iz=1,nz
        READ(9)Z ; WRITE(8,*)Z
        DO ix=1,nx
         READ(9)X ; WRITE(8,*) X
         READ(9)((kap(it,ir),ir=1,nr),it=1,nt)
         WRITE(8,*)((kap(it,ir),ir=1,nr),it=1,nt)
        ENDDO
       ENDDO
      ENDIF

      CLOSE(unit=8) ; CLOSE(unit=9)

      CALL pause('OK c''est fait')

      STOP

      END PROGRAM ascii2bin_opa
