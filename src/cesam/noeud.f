
c**************************************************************

	SUBROUTINE noeud(x,y,mi,n,knot)
	
c	subroutine public du module mod_numerique	

c	formation de la suite croissante des noeuds y(knot) à partir
c	de la suite strictement croissante des points de raccordement
c	x(n)
c	mi(n) : multiplicité des noeuds définie par :
c		mi(1):=m, mi(n):=m
c	au point x(i) 1<i<n, la dérivée d'ordre j=m-(mi(i)+1) est continue
c	mi(i)=m -> j=-1 (discontinuité)	* mi(1)=1 -> j=m-2=2 (pour m=4)
c	     =m-1    0   		*       2      m-3=1
c	      m-2    1	(dérivée 1ière)	*       3      m-4=0
c	      m-3    2  (dérivée 2de)   *       4      m-5=-1

c entrees:
c	x : suite strictement croissante des points de raccordement
c	mi : vecteur de multiplicité des noeuds
c	n : nombre de points

c sorties:
c	y : vecteur nodal
c	knot : nb. de points du vecteur nodal

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: x
	INTEGER, INTENT(in), DIMENSION(:) :: mi
	INTEGER, INTENT(in) :: n
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: y
	INTEGER, INTENT(out) ::	knot
	
	INTEGER :: i, k
	
c-------------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(3es22.15)

c	vérification de la stricte croissance de la suite des x

	DO i=1,n-1
	 IF(x(i) >= x(i+1))THEN
	  PRINT*,'dans noeud la suite des abscisses n''est pas',
	1 ' strictement croissante en i=',i
	  PRINT*,'nombre de points: ',n ; PRINT*,'abscisse x=',x(i)
	  WRITE(*,2001)x(1:i-1) ; WRITE(*,2001)x(i:i+1) ; WRITE(*,2001)x(i+2:n) 
	  no_croiss=.TRUE. ; RETURN
	 ENDIF
	ENDDO

	knot=0

	DO k=1,mi(1)
	 knot=knot+1 ; y(knot)=x(1)
	ENDDO

	DO i=2,n-1
	 DO k=1,mi(i)
	  knot=knot+1 ; y(knot)=x(i)
	 ENDDO
	ENDDO

	DO k=1,mi(n)
	 knot=knot+1 ; y(knot)=x(n)
	ENDDO

c	PRINT*,'noeud',n,knot ; PRINT*,x(1:n) ; PRINT*,mi(1:n)
c	WRITE(*,2000)y(1:knot) ; PRINT*,y(1:knot) ; pause

	RETURN

	END SUBROUTINE noeud
