
c**********************************************************************

	SUBROUTINE noeu_dis(eps,id,knot,m,nd,nx,x,xt)	
	
c	subroutine private du module mod_numerique

c	détermine la séquence de noeuds de raccord xt(knot)
c	pour une interpolation "optimale"
c	en tenant compte de discontinuités
c	identique à noein s'il n'y a pas de discontinuité i.e. nd=0

c	par B-splines d'ordre m sur la suite strictement
c	croissante x de n points de données, cf. de Boor p.219 formule (10)
c	aux limites le polynome d'interpolation s'appuie sur m points
c	de table

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.,
c	version: 06 06 03
c	CESAM2k

c entrees :
c	eps : écart à droite, parameter définit avant l'appel à bsp_dis
c	x : abscisses strictement croissantes
c	nx : nombre de points
c	id(0:nd+1) : indice des discon. avec, en sortie id(0)=1, id(nd+1)=n
c	nd : nombre de discontimuites
c	m : ordre des splines

c sorties :
c	xt : points de table (nodaux)
c	knot : nombre de points de table, knot=nx+nd+m

c-------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: x
	REAL (kind=dp), INTENT(in) :: eps
	INTEGER, INTENT(in) :: m, nd, nx
	INTEGER, INTENT(inout), DIMENSION(0:) :: id	
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: xt
		
	INTEGER, INTENT(out) :: knot	
	
	REAL (kind=dp) :: mm1

	INTEGER :: i, j, ij
	
c-----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'NOEU_DIS nx,nd,m,id/x',nx,nd,m,id(0:nd+1)
c	WRITE(*,2000)x(1:nx) ; WRITE(*,2000)eps ; PAUSE

c	cas nd <= 0, on construit le vecteur nodal de noein

	IF(nd <= 0)THEN
	 PRINT*,'dans noeu_dis le nombre de discontinuités nd =',nd
	 PRINT*,'nd <= 0, on ne tient pas compte de discontinuité'
	 CALL noein(x,xt,nx,m,knot)

c	 cas nd > 0, on construit le vecteur nodal avec discontinuités
	 
	ELSE
	
c	 vérification de la stricte croissance de la suite des x

	 DO i=1,nx-1
c	  PRINT*,i,nx ; WRITE(*,2000)x(i) ; WRITE(*,2000)x(i+1)
	  IF(x(i) >= x(i+1))THEN
	   PRINT*,'dans noeu_dis la suite des abscisses n''est pas',
	1  ' strictement croissante, en i=',i
	   PRINT*,'nombre de points: ',nx ; PRINT*,'abscisses: '
	   PRINT*,x(1:i); PRINT* ; PRINT*,x(i+1:nx)
	   no_croiss=.TRUE. ; RETURN	   
	  ENDIF
	 ENDDO

c	 pour l'interpolation spline il faut n >= m

	 IF(nx < m)THEN
	  WRITE(*,11)nx,m ; STOP
11	  FORMAT('dans noeu_dis nx=',i3,' < ',i3,'=m')
	 ENDIF
	 	
	 DO i=2,nd
	  IF(id(i)-id(i-1) < m)THEN
	   PRINT*,'dans noeu_dis on impose m+1 points entre 2 discon.'
	   PRINT*,'nombre de discontinuités:',nd
	   PRINT*,'indices des discontinuités:',id(1:nd) ; STOP
	  ENDIF
	 ENDDO
	 
c	 initialisations 

	 mm1=m-1 ; id(0)=1 ; id(nd+1)=nx ; knot=0
	 
c	 petit écart normalisé à gauche, pour cohérence
c	 eps est un parameter global définit dans bsp_dis	 
c	 m points de table en x(1)

	 DO i=1,m
	  knot=knot+1 ; xt(knot)=x(1)-eps*(x(2)-x(1))
	 ENDDO

	 DO ij=1,nd+1			!entre chaque discontinuité
	  DO i=id(ij-1),id(ij)-m	!entre 2 discontinuités
	   knot=knot+1 ; xt(knot)=0.d0
	   DO j=i+1,i+m-1		!moyenne de Schomberg
	    xt(knot)=xt(knot)+x(j)
	   ENDDO
	   xt(knot)=xt(knot)/mm1
	  ENDDO		!i
c	  PRINT*,ij,id(ij),id(ij)+1 ; PAUSE'ij'

c	  à chaque discontinuité petit écart normalisé à droite	  
	  
	  IF(ij /= nd+1)THEN
	   DO j=1,m			!à chaque discontinuité
	    knot=knot+1 ; xt(knot)=x(id(ij))+eps*(x(id(ij)+1)-x(id(ij)))
	   ENDDO	!j
	  ELSE
	   DO j=1,m	!à l'extémité
	    knot=knot+1 ; xt(knot)=x(nx)+eps*(x(nx)-x(nx-1))
	   ENDDO	  
	  ENDIF
	   
	 ENDDO		!ij
	ENDIF

c	la dernière couche a été traitée comme discontinuité en id(nd+1)

c	PRINT*,knot,nx+nd+m,nx,m,nd ; WRITE(*,2000)xt
	
	RETURN

	END SUBROUTINE noeu_dis
