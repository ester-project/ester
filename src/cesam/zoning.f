
c***********************************************************************

	SUBROUTINE zoning(n,m,p,newn,newm)
	
c	routine public du module mod_numerique

c	détermine les abscisses newm de façon a répartir la fonction p
c	a peu près uniformement

c	on utilise une interpolation linéaire inverse

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	CESAM95 

c entree
c	n : nombre de points
c	m : abscisses
c	p : fonction de répartition
c	newn : nombre de points de la nouvelle répartition

c sortie
c	newm : abscisses assurant une répartition uniforme de p

c----------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: m, p
	INTEGER, INTENT(in) :: n, newn
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: newm

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: y
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: xt, x
	REAL (kind=dp), DIMENSION(1) :: fx, dfxdx
	REAL (kind=dp) :: xi, pas
	INTEGER :: i, knot, l
	
c------------------------------------------------------------

2000	FORMAT(8es10.3)

	ALLOCATE(x(n),xt(n+2),y(1,n)) ; y(1,1:n)=m(1:n) ; x(1:n)=p(1:n)

c	PRINT*,'zoning.f: x,y' ; WRITE(*,2000)x(1:n)
c	WRITE(*,2000)y(1,1:n) ; PRINT*,'newn',newn
	pas=(x(n)-x(1))/(newn-1)
	newm(1)=m(1)
	newm(newn)=m(n)
	l=2
	CALL bsp1dn(1,y,x,xt,n,2,knot,.FALSE.,x(1),l,fx,dfxdx)
        IF(no_croiss)THEN
         PRINT*,'Arrêt dans zoning' ; STOP
        ENDIF
	DO i=2,newn-1
	 xi=x(1)+pas*(i-1)
	 CALL bsp1dn(1,y,x,xt,n,2,knot,.TRUE.,xi,l,fx,dfxdx)		 
	 newm(i)=fx(1)
c	 WRITE(*,2000)xi,newm(i)
	ENDDO
	
	DEALLOCATE(x,xt,y)

	RETURN

	END SUBROUTINE zoning
