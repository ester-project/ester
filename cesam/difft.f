         	
c****************************************************************
	 
	SUBROUTINE difft(melange,nu,t,ro,drox,kap,dkapx,deff,d,dd)

c	subroutine private du module mod_evol
	 	
c	subroutine générique de calcul des coefficients de
c	diffusion turbulente

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : nom_difft
	USE mod_kind	 

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: deff, dkapx, kap, nu, t
	LOGICAL, INTENT(in) :: melange			
	REAL (kind=dp), INTENT(inout) :: drox, ro
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d
	 
c-------------------------------------------------------------------------
	 
	SELECT CASE(nom_difft)
	CASE('difft_nu')
	 CALL difft_nu(melange,t,ro,drox,kap,dkapx,deff,d,dd)
	CASE('difft_ventura')
	 CALL difft_ventura(melange,nu,t,ro,drox,kap,dkapx,deff,d,dd)
	CASE DEFAULT
	 PRINT*,'routine de diffusion turbulente inconnue: ',nom_difft
	 PRINT*,'routines connues: difft_nu' ; PRINT*,'arrêt' ; STOP
	END SELECT
	
	RETURN
	 
	END SUBROUTINE difft
