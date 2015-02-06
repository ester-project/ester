          	
c********************************************************************

	SUBROUTINE diffm(p,t,r,l,m,ro,drox,kap,dkapx,w,
	1 gradrad,dgradradx,xi,d,dd,v,dv)

c	subroutine private du module mod_evol
	
c	subroutine générique de calcul des coefficients de
c	diffusion microscopique

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : nom_diffm
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xi
	REAL (kind=dp), INTENT(in) :: p, t, m, l, r, ro, drox, kap, dkapx,
	1 gradrad, dgradradx, w
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d, dv
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: v

c---------------------------------------------------------------------

	SELECT CASE(nom_diffm)
	CASE('diffm_br')
	 CALL diffm_br(t,r,l,m,ro,drox,kap,dkapx,gradrad,dgradradx,w,
	1 xi,d,dd,v,dv)
	CASE('diffm_mp')
	 CALL diffm_mp(p,t,r,m,ro,drox,w,gradrad,dgradradx,xi,d,dd,v,dv)
	CASE('diffm_0')
	 dd=0.d0 ; dv=0.d0 ; d=0.d0 ; v=0.d0
	CASE DEFAULT
	 PRINT*,'routine de diffusion microscopique inconnue: ',nom_diffm
	 PRINT*,'routines connues: diffm_br, diffm_mp, diffm_0' 
	 PRINT*,'arrêt' ; STOP
	END SELECT
 
	RETURN

	END SUBROUTINE diffm
