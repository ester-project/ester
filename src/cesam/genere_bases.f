
c*******************************************************

	SUBROUTINE genere_bases(id,m,mdis,n,nd, x,knot,xdt,xct)

c	routine PUBLIC du module mod_numerique
	
c	Formation de la base avec discontinuités
c	et de la base continue associée

c	la fonction d'interpolation est dérivable jusqu'à l'ordre m-1

c	1 <= mdis(id) <= m
c	mdis(id)=m : discontinuité
c	mdis(id)=m-1 : fonction continue non dérivable
c	mdis(id)=m-2 : fonction continue, 1 fois dérivable
c	mdis(id)=m-3 : fonction continue, 2 fois dérivable
c	.	.	.	.	.	.	.
c	mdis(id)=m-(m-1)=1 : fonction continue, m-1 fois dérivable

c	nd=0 : fonction continue, m-1 fois dérivable 

c entrées
c	id(0:nd+1) : indices des discontinuités
c	m : ordre des splines
c	mdis(nd) : ordre de la discontinuité
c	n : nombre d'abscisses	
c	nd : nombre de discontinuités
c	x(n) : abscisses

c sorties:
c	knot : nombre de point des vecteurs nodaux
c	xct : vecteur nodal de la base continue associée
c	xdt : vecteur nodal de la base discontinue

c	Auteur: P. Morel, Département Cassiopée, O.C.A.

c	CESAM2k

c------------------------------------------------------------------

	use mod_kind
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: x
	INTEGER, INTENT(in), DIMENSION(0:) :: id
	INTEGER, INTENT(in), DIMENSION(:) :: mdis	
	INTEGER, INTENT(in) :: m, n, nd
	
	REAL (kind=dp), INTENT(out), ALLOCATABLE, OPTIONAL, DIMENSION(:) :: xct
	REAL (kind=dp), INTENT(out), ALLOCATABLE, DIMENSION(:) :: xdt
	
	INTEGER, INTENT(out) :: knot
	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: xc
	REAL (kind=dp) :: pas
	INTEGER, ALLOCATABLE, DIMENSION(:) :: mult
	INTEGER :: dimc, i, j, k	
	
c------------------------------------------------------------------

2000	FORMAT(8es10.3)

	ALLOCATE(mult(n))
	mult=1 ; mult(1)=m ; mult(n)=m
	IF(nd == 0)THEN
	
c nd=0, pas de discontinuité, xct, xdt sont identiques	
	 knot=SUM(mult)
	 ALLOCATE(xdt(knot))
	 CALL noeud(x,xdt,mult,n,knot)
	 IF(PRESENT(xct))THEN
	  ALLOCATE(xct(knot)) ; xct=xdt
	 ENDIF
	
	ELSE
c formation de la base discontinue
	 mult(id(1:nd))=mdis(1:nd)
	 knot=SUM(mult)
	 ALLOCATE(xdt(knot))
	 CALL noeud(x,xdt,mult,n,knot)

c formation de la base continue
	 IF(PRESENT(xct))THEN
	  dimc=knot-2*(m-1)
	  ALLOCATE(xc(dimc),xct(knot))
	  xc(1:id(1)-1)=x(1:id(1)-1) ; k=id(1)-1	
	  B1: DO j=1,nd
	   pas=(x(id(j)+1)-x(id(j)-1))/REAL(mult(id(j))+1,dp)
	   DO i=1,mult(id(j))
	    k=k+1 ; xc(k)=x(id(j)-1)+pas*REAL(i,dp)
	   ENDDO 	 
	   IF(j == nd)EXIT B1
	   DO i=id(j)+1,id(j+1)-1
	    k=k+1 ; xc(k)=x(i)
	   ENDDO
	  ENDDO B1
	  DO i=id(nd)+1,n
	   k=k+1 ; xc(k)=x(i)
	  ENDDO
	  DEALLOCATE(mult) ; ALLOCATE(mult(dimc))	
	  mult=1 ; mult(1)=m ; mult(dimc)=m
	  CALL noeud(xc,xct,mult,dimc,knot)
	  DEALLOCATE(xc)
	 ENDIF
	ENDIF
	DEALLOCATE(mult)

	RETURN
	
	END SUBROUTINE genere_bases
