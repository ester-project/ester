
c****************************************************************

	SUBROUTINE dnunl(r,c,n,nu0,dnu02,dnu13,a)
	
croutine private du module mod_cesam	
c calcul de nu0 et delta nu  n l, formule 100, 101, 102 p.389 Schatzman et
c Praderie,
c Formule 102: delta nu 02 ~ 6A/20, delta nu 13 ~ 10A/20, J.Provost

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c entrées:
c	r(n): rayons du centre à l'extérieur (Rsol)
c	c(n): vitesse du son (cgs)
c	n: nombre de points
c	rtot: rayon (peut différer de r(n))

c sorties
c	nu0: fondamental
c	dnu02, dnu13: écarts approximatifs de fréquence delta nu 02 et 13
c	a: formule 100, ~ - somme 1/r dc/dr = 2 somme dc /dr**2

c----------------------------------------------------------------

	USE mod_donnees, ONLY : pi
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, intgauss, no_croiss
	
	IMPLICIT NONE

	INTEGER, INTENT(in) :: n
	REAL (kind=dp), INTENT(in), DIMENSION(n) :: c, r
	REAL (kind=dp), INTENT(out) :: nu0, dnu02, dnu13, a


	INTEGER, PARAMETER :: m=4, mg=2	!ordre interp. m/2 pour int. Gauss	
	REAL (kind=dp), DIMENSION(1,n) :: ci  
	REAL (kind=dp), DIMENSION(n) :: r2
	REAL (kind=dp), DIMENSION(n+m) :: r2t
	REAL (kind=dp), DIMENSION(mg) :: rg, wg
	REAL (kind=dp), DIMENSION(1) :: cr, dcdr	

	INTEGER :: i, ig, knotr, l
		
c----------------------------------------------------------------

2000	FORMAT(8es10.3)

c les abscisses doivent être strictement croissantes, sinon pas de calcul
	DO i=1,n-1
	 IF(r(i) >= r(i+1))THEN
	  nu0=0.d0 ; dnu02=0.d0 ; dnu13=0.d0 ; a=0.d0 ; RETURN
	 ENDIF
	ENDDO
		
c tabulation de c(r**2)
	r2=r**2
	ci=RESHAPE(c,SHAPE(ci))		!ie. : ci(1,1:n)=c(1:n)
	
	CALL bsp1dn(1,ci,r2,r2t,n,m,knotr,.FALSE.,r2(1),l,cr,dcdr)
        IF(no_croiss)THEN
         PRINT*,'Arrêt 1 dans dnunl' ; STOP
        ENDIF

c recherche de la limite et calcul des intégrales
c calcul approximatif: on ne tient pas compte de l'atmosphère
	a=0.d0 ; nu0=0.d0
	DO i=1,n-1
	 CALL intgauss(r(i),r(i+1),rg,wg,mg) !int.Gauss
	  DO ig=1,mg
	   CALL bsp1dn(1,ci,r2,r2t,n,m,knotr,.TRUE.,rg(ig)**2,l,cr,dcdr)
	   nu0=nu0+wg(ig)/cr(1) ; a=a+wg(ig)*dcdr(1)
	  ENDDO
	ENDDO
	a=2.d0*a ; nu0=2.d0*nu0 ; nu0=1.d6/nu0
	
	CALL bsp1dn(1,ci,r2,r2t,n,m,knotr,.TRUE.,r(n)**2,l,cr,dcdr)
	a=(cr(1)/r(n)-a)/(2.*pi)**2*1.d6	!formule 100, p.389, modif
	dnu02=6.d0*a/20.d0 ; dnu13=10.d0*a/20.d0

	RETURN

	END SUBROUTINE dnunl
