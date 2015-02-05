
c*******************************************************************

	SUBROUTINE newspl(n,x,x1t,kno1,m1,x2t,kno2,m2,s1,s2)

c	subroutine public du module mod_splines95

c	pour n fonctions
c	transforme la spline s1 d'ordre m1 sur le réseau x1t \ kno1,
c	en la spline s2 d'ordre m2 sur le réseau x2t \ kno2
c	on n'utilise pas la base duale

c	le cas ou la spline s2 a des discontinuités est envisagé:
c	la diccontinuité ne peut être localisée avec la précision
c	machine, on se déplace légèrement à gauche et à droite

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k (version f95 de newspl1.f)

c entrées
c	n: nombre de fonctions
c	x, x1t, x2t, kno1, kno2: abcisses et vecteurs nodaux
c	m1, m2: ordres des splines
c	s1: première spline

c sorties:
c	s2: nouvelle spline

c--------------------------------------------------------------------------

	USE mod_kind
		
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION (:) :: x, x2t
	INTEGER, INTENT(in) :: kno2, m1, m2, n
	REAL (kind=dp), INTENT(inout), DIMENSION (:,:) :: s1
	REAL (kind=dp), INTENT(inout), DIMENSION (:) :: x1t
	INTEGER, INTENT(inout) :: kno1
	REAL (kind=dp), INTENT(out), DIMENSION (:,:) :: s2
	
	REAL (kind=dp), DIMENSION (kno2-m2,m2) :: a
	REAL (kind=dp), DIMENSION (MAX(m1,m2)+1) :: q
	REAL (kind=dp), DIMENSION (n) :: dfxdx, fx	
	
	REAL (kind=dp), PARAMETER :: dx=1.d-3 , un_dx=1.d0-dx
	REAL (kind=dp) :: p
			
	INTEGER, DIMENSION(kno2-m2) :: indpc
		
	INTEGER :: i, l=1, n1=0, n2
	LOGICAL :: inversible

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)
	
c dimension de la spline 2
	n2=kno2-m2
		
c calcul des tau selon de-Boor p.123

	a=0.d0 ; i=0
	B1: DO
	 i=i+1 ; IF(i > n2)EXIT B1
	 p=SUM(x2t(i+1:i+m2-1))/REAL(m2-1,dp)	!abscisses de de Boor
	 
c détection d'une discontinuité pour s2	 
	 IF(i > 1 .AND. i < n2 .AND. p == x2t(i+1)) THEN
	 	  
c on se place à gauche de la discontinuité 	  
	  p=x2t(i+1)*un_dx+x2t(i)*dx	  	  
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,p,l,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. en 1 dans newspl'	  	  
	  s2(:,i)=fx
	  CALL linf(p,x2t,kno2,l) ; CALL bval0(p,x2t,m2,l,q)
	  a(i,1:m2)=q(1:m2) ; indpc(i)=l-m2+1
c	  PRINT*,'à gauche, i=',i
c	  PRINT*,x2t(i),x2t(i+1),p
	  
c on se place à droite de la discontinuité
	  i=i+1
	  p=x2t(i)*un_dx+x2t(i+m2)*dx	  	  
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,p,l,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. en 2 dans newspl'	  	  
	  s2(:,i)=fx
	  CALL linf(p,x2t,kno2,l) ; CALL bval0(p,x2t,m2,l,q)
	  a(i,1:m2)=q(1:m2) ; indpc(i)=l-m2+1
c	  PRINT*,'à droite, i=',i
c	  PRINT*,x2t(i),x2t(i+m2),p
	  	  
c sans discontinuité	  	  
	 ELSE	  
c	  PRINT*,i,p,x1t(1),x2t(1)
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,
	1   MAX(x1t(1),MIN(p,x1t(kno1))),l,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. en 3 dans newspl'
	  s2(:,i)=fx
	  CALL linf(p,x2t,kno2,l) ; CALL bval0(p,x2t,m2,l,q)
	  a(i,1:m2)=q(1:m2) ; indpc(i)=l-m2+1	  
	 ENDIF 
	ENDDO	B1

c coefficients de s2

	CALL gauss_band(a,s2,indpc,n2,n2,m2,n,inversible)
	IF(.NOT.inversible)PRINT*,'matrice non inversible dans newspl'
		
	RETURN

	END SUBROUTINE newspl
