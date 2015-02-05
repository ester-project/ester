
c******************************************************************

      SUBROUTINE iben(t,dcomp,jac,deriv,fait,
     1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c	routine private du module mod_nuc

c	fausses reactions thermonucléaires pour modèle initial de
c	pre-main sequence: epsilon=cT d'après Iben

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées :
c	t : température
c	fait=1 : initialisation de la composition chimique
c	    =2 : derivée première
c	    =3 : produit c*t

c sorties
c	dcomp : dérivée temporelle (unité de temps : 10**6 ans)

c----------------------------------------------------------------

	USE mod_kind
	USE mod_variables, ONLY : c_iben
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: t 
	INTEGER, INTENT(in) :: fait	  	  
	LOGICAL, INTENT(in) :: deriv  
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: jac
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dcomp, epsilon, ex
	REAL (kind=dp), INTENT(out) :: et, ero, hhe, be7e, b8e, n13e,
	1 o15e, f17e
	
c-------------------------------------------------------------------------

2000	FORMAT(8es10.3)
	
	SELECT CASE(fait)
	CASE(1)
	 WRITE(*,1) ; WRITE(2,1)
1	 FORMAT(/,'méthode de Iben pour recherche du modèle initial PMS',/)
	CASE(2)
	 dcomp=0.d0
	 IF(deriv)jac=0.d0	  
	CASE(3)
	 epsilon=0.d0 ; epsilon(1)=c_iben*t
c	 PRINT*,'Iben : epsilon(1),c_iben,t'
c	 WRITE(*,2000)epsilon(1),c_iben,t
	 IF(deriv)then
	  ero=0.d0 ; et=c_iben ; ex=0.d0
	 ENDIF
	CASE(4)  !neutrino
	 hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0 ; n13e=0.d0 ; o15e=0.d0
	 f17e=0.d0
	END SELECT
	
	RETURN

	END SUBROUTINE iben
