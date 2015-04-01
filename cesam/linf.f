
c********************************************************************

	SUBROUTINE linf(x,y,n,l)

c subroutine public du module mod_numerique   

c localisation dans une table
c inspire de l'algorithme HUNT p. 91 de Numerical Recipes

c y(n) : suite croissante (non strictement) des points nodaux
c recherche de l'indice l tel que :
c	    y(l) .le. x < y(l+1)

c s'il y a débordement à gauche ou à droite on prend:
c x < y(1) <= y(l) < y(l+1) ou y(l) < y(l+1) <= y(n) < x

c on commence la recherche au voisinage de l

c Auteur: P. Morel, Département J.D. Cassini, O.C.A., CESAM2k

c entrées:
c	x: abscisse a localiser
c	y: table
c	n: nombre de points

c entrée/sortie
c	l: indice

c----------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
    
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: y
	REAL (kind=dp), INTENT(in) :: x
	INTEGER, INTENT(in) :: n    
	INTEGER, INTENT(inout) :: l

	INTEGER :: inc, h, m

c----------------------------------------------------------

	IF(x <= y(1))THEN   	!test aux limites
	 no_croiss= x < y(1)
	 IF(no_croiss)PRINT*,'linf : x =',x,' < y(1) =',y(1)
	 l=1
	 DO WHILE(y(l+1) <= y(1))
	  l=l+1     !à cause des points mutiples aux extrèmités
	 ENDDO
	ELSEIF(x >= y(n))THEN
	 no_croiss= x > y(n)
	 IF(no_croiss)PRINT*,'linf : x =',x,' > y(n) =',y(n)
	 l=n-1 
	 DO WHILE(y(l) >= y(n))
	  l=l-1
	 ENDDO  
	ELSE
	 IF(l >= 1 .AND. l < n)THEN
	  inc=1
	  IF(x >= y(l))THEN
	   h=l+inc
	   DO WHILE (x >= y(h))
	    l=h ; inc=inc+1 ; h=MIN(n,l+inc)
	   ENDDO
	  ELSE
	   h=l ; l=l-inc
	   DO WHILE(x < y(l))
	    h=l ; inc=inc+1 ; l=MAX(1,l-inc)
	   ENDDO
	  ENDIF 
	 ELSE
	  l=1 ; h=n
	 ENDIF
 
	 DO WHILE(h-l > 1)
	  m=(l+h)/2
	  IF(x >= y(m))THEN
	   l=m
	  ELSE
	   h=m
	  ENDIF
	 ENDDO
	ENDIF

	RETURN

	END SUBROUTINE linf
