		
c**************************************************************

	SUBROUTINE matinv(a,n,inversible)

	
c	inversion de matrice par elimination de Gauss, pivot partiel
c	utilisation de gauss_band

c entree:
c	n: ordre

c entree/sortie:
c	a: matrice/matrice inverse

c sortie:
c	inversible=.false.: matrice non inversible

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.

c----------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: n	
	REAL (kind=dp), INTENT(inout), DIMENSION(n,n) :: a	
	LOGICAL, INTENT(out) :: inversible
			
	REAL (kind=dp), DIMENSION(n,n) :: am1
	INTEGER, DIMENSION(n) :: indpc	
	INTEGER :: i

c---------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'matinv',n

c	identite dans am1

	am1=0.d0 ; indpc=1
	DO i=1,n
	 am1(i,i)=1.d0
	ENDDO
	
	CALL gauss_band(a,am1,indpc,n,n,n,n,inversible)
		
	IF(.NOT.inversible)THEN
	 PRINT*,'dans matinv matrice non inversible' ; STOP
	ENDIF
	
c	transposition car avec gauss_band on a resolu A t(A^-1) = I
c	la matrice second membre etant stockee (colonnes,lignes)

	a=transpose(am1)
		
	RETURN
	
	END SUBROUTINE matinv
