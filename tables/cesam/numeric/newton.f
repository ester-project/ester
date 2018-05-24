
c********************************************************************

	SUBROUTINE newton(l,k,f,x,z,p,der)

c	routine public du module mod_numerique

c	polynome de newton

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
c	F95: 29 08 00

c entree
c	l : premier indice
c	k : degre du polynome
c	f : table des f
c	x : abscisses
c	z : point d'interpolation
c	der : ordre de la derivee

c sortie
c	p : valeur des derivees

c-----------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(0:) :: f, x
	REAL (kind=dp), INTENT(in) :: z	
	INTEGER, INTENT(in) :: l, k, der
	REAL (kind=dp), INTENT(out), DIMENSION(0:) :: p	
	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: a, ah, xh
	INTEGER :: j
	
c-------------------------------------------------------------------

	ALLOCATE(a(0:k),xh(0:k),ah(0:k))
	
c	PRINT*,l,k,der,z

	CALL difdiv(l,k,f,x,a)

	DO j=0,k
	 xh(j)=x(l+j) ; ah(j)=a(j)
	ENDDO
c	PRINT*,(f(j+l),j=0,k) ; PRINT*,(xh(j),j=0,k) ; PRINT*,(ah(j),j=0,k)
	CALL horner(k,der,z,ah,xh,p)
	
	DEALLOCATE(a, ah, xh)

	RETURN

	END SUBROUTINE newton
