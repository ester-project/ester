
c***************************************************

	FUNCTION neville(x,xi,fi,n)

c	routine public du module mod_numerique

c	algorithme de Neville

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	inspire de Stoer & Burlirsh 2.1.2 p. 42
c	18 12 00: version F95

c entrees:
c	x: point d'interpolation
c	xi: abscisses
c	fi: fonction
c	n: degré du polynome

c----------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(0:) :: xi, fi
	REAL (kind=dp), INTENT(in) :: x
	INTEGER, INTENT(in) :: n
	
	REAL (kind=dp) :: neville
	
	REAL (kind=dp), DIMENSION(0:n) :: t
	INTEGER :: i,j

c----------------------------------------------------------

2000	format(8es10.3)
c	WRITE(*,2000)xi(0:n) ; WRITE(*,2000)fi(0:n)

	DO i=0,n
	 t(i)=fi(i)
	 DO j=i-1,0,-1
c	  WRITE(*,2000)xi(i),xi(j),float(i),float(j)
	  t(j)=t(j+1)+(t(j+1)-t(j))*(x-xi(i))/(xi(i)-xi(j))
	 ENDDO	!j
	ENDDO	!i
	neville=t(0)

	RETURN

	END FUNCTION neville
