
c******************************************************************

	SUBROUTINE difdiv(l,k,f,x,a)

c	differences divisees

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	version: 03 06 91
c	F95 29 08 00

c entree
c	l : premier indice
c	k : degre du polynome
c	f : fonction
c	x : abscisses

c sortie
c	a : differences divisees

c---------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(0:*) :: f,x
	INTEGER, INTENT(in) :: l,k
	REAL (kind=dp), INTENT(out), DIMENSION(0:*) :: a
		
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: d
	
	INTEGER :: i,j
	
c-------------------------------------------------------------------------

c	PRINT*,(f(l+j),j=0,k) ; PRINT*,(x(l+j),j=0,k)

	ALLOCATE(d(0:k))

	DO j=0,k
	 d(j)=f(l+j)
	 DO i=j-1,0,-1
	  d(i)=(d(i+1)-d(i))/(x(l+j)-x(l+i))
	 ENDDO
	 a(j)=d(0)
	ENDDO	
c	PRINT*,'les a' ; PRINT*,(a(j),j=0,k)

	DEALLOCATE(d)

	RETURN

	END SUBROUTINE difdiv
