	
c***********************************************************************

	SUBROUTINE bvald(x,xt,m,l,r,d)
	
c	subroutine public du module mod_numerique	

c	d(r,j):= derivee d'ordre r,  0 <= r < m-1, de
c	la l-m+j-ieme B-splines, 1 <= j <= m, d'ordre m non nulle
c	au point x, xt(l) <=  x <  xt(l+1).
c	les derivees de 0 a r sont calculees
c	on doit avoir d(0:r,m) dans le programme appelant
c	n'a d'interet que pour r >= 2
c	pour r=0 / r=1, on utilise bval0 / bval1 qui sont plus efficaces

c	version F95 de bvald

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	optimisation Th. Corbard

c entrees:
c	x: point de calcul
c	xt: vecteur nodal
c	m: ordre des splines
c	r: ordre de derivation
c	l:  xt(l) <=  x <  xt(l+1)
c sortie:
c	d: tableau des derivees d(r,m)

c-----------------------------------------------------------------------
	
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xt
	REAL (kind=dp), INTENT(in) :: x
	INTEGER, INTENT(in) :: l, m, r
	REAL (kind=dp), INTENT(out), DIMENSION(0:,:) :: d
	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: c, n
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: dd, q
	
	REAL (kind=dp) :: a1, a2, denom
			
	INTEGER :: j, i, k, maxi, mini, mj, p
	
c	DEALLOCATE(c,n,dd,q)

c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

	ALLOCATE(q(m+1))

	IF(r == 0)THEN
	 CALL bval0(x,xt,m,l,q) ; d(0,1:m)=q(1:m) ; DEALLOCATE(q)
	ELSEIF(r == 1)THEN
	 ALLOCATE(dd(m))
	 CALL bval1(x,xt,m,l,q,dd) ; d(0,1:m)=q(1:m) ; d(1,1:m)=dd(1:m)
	 DEALLOCATE(dd,q) 	
	ELSE
	
c	 n(i,j) matrice  triangulaire inferieure
c	 des valeurs en x des B-splines d'ordre i=m, m-1, ..., 1 non nulles
c	 d'indices l-m+p, ..., l, p=m, m-1, ..., 1
c	 adapte de Schumaker p. 192 theoreme 5-9 en prenant pour spline
c	 chacune des B-splines normalissees non identiquement nulles en x
c	 avec pour coefficients c_ij=delta_ij dans la somme de Schumaker
c	 N_j^m(x)= sum_i=1,n c_ij N_j^m(x)

c	 au fur et a mesure on forme les B-splines normalisees N
c	 a partir des Q en mulitipliant par y(i+m)-y(i)  formule 4.29
c	 les q etant les Q de cette formule

	 ALLOCATE(n(m,m)) ; n=0.d0 ; q=0.d0 ; q(m)=1.d0/(xt(l+1)-xt(l))
	 n(m,m)=1.d0	!pour la splines d'ordre m=1
	 DO j=2,m	!pour des splines d'ordre 2 de a m
	  DO i=m-j+1,m	!pour les j splines non nulles
	   denom=xt(i+l-m+j)-xt(i+l-m) ; a1=(x-xt(i+l-m))/denom
	   a2=1.d0-a1 ; q(i)=a1*q(i)+a2*q(i+1)
	   n(m-j+1,i)=q(i)*denom  !N(m-ordre+1,spline l-k+1,..., l)
	  ENDDO
	 ENDDO

c	 On forme la matrice
c	 c(j,i) pour chacune des B-splines non nulle sur [xt(l),xt(l+1)[
c	 adapte de Schumaker 5-9 en prenant c(0,0)=1
c	 j=0,r est l'indice de "derivation", i=1,r est l'indice de spline

c	 il faut dimensionner la seconde coordonnee de c de -1 a m
c	 avec c(*,0)=0 qui est la valeur du coefficient=0 de
c	 la B-spline precedent la premiere non nulle ie. celle d'indice l-m 

	 ALLOCATE(c(0:r,-1:m))	!c(j,i)=c_i^j

	 d=0.d0 ; d(0,:)=n(1,:)	!derivee d'ordre 0
	 DO p=1,m
	  k=l-m+p ; c=0.d0 ; c(0,0)=1.d0	!coefficient de la spline N_k^m
	  DO j=1,r
	   mj=m-j
	   DO i=0,min(j,m-p)	!optimisation Th. Corbard
	    denom=xt(k+i+m-j)-xt(k+i)
	    IF(denom > 0.d0)c(j,i)=mj/denom*(c(j-1,i)-c(j-1,i-1))
	   ENDDO	!i
	  ENDDO	!j
	  DO j=1,min(m-2,r)	!pour les derivees opt. Th. Corbard
	   mini=max(j+1,p)	!optimisation Th. Corbard
	   maxi=min(j+p,m)	!optimisation Th. Corbard
	   d(j,p)=sum(c(j,mini-p:maxi-p)*n(j+1,mini:maxi))	!Th. Corbard
	  ENDDO		!j
	  IF(r == m-1)d(r,p)=c(r,m-p)	!optimisation Th. Corbard
	 ENDDO		!k
	 DEALLOCATE(c,q,n)
	ENDIF

	RETURN

	END SUBROUTINE bvald
