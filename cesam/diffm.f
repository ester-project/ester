          	
c********************************************************************

	SUBROUTINE diffm(p,t,r,l,m,ro,drox,kap,dkapx,w,
	1 gradrad,dgradradx,xi,d,dd,v,dv)

c	subroutine private du module mod_evol
	
c	subroutine générique de calcul des coefficients de
c	diffusion microscopique

c entrées
c     p,t,r,lum,ltot,m : pression, temperature, lum. locale et totale, masse
c     ro,drox,kap,dkapx : densite, dérivée/ X, opacite, dérivée / X 
c     gradrad,dgradradx : gradient rad. et dérivées / X
c     xi : comp. chim. par mole
c     les dérivées / X, issues de thermo, sont par gramme 

c sorties
c     d, dd : coefficients d_ij de d x_j / d m et dérivées / x_k
c     v, dv : coefficients v_i de x_i et dérivées / x_k

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

2000 	FORMAT(8es10.3)
c	WRITE(*,2000)r,w

	SELECT CASE(nom_diffm)
	CASE('diffm_br')
	 CALL diffm_br(t,r,l,m,ro,drox,kap,dkapx,w,gradrad,dgradradx,
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
