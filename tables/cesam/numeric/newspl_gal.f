	
c******************************************************************

	SUBROUTINE newspl_gal(n,x,x1t,knot1,m1,s1,x2t,knot2,m2,s2)

c routine public du module mod_numerique
	
c On utilise le formalisme des éléments finis pour effectuer un changement de
c base. Les bases de spline sont quelquonques
c en particulier, elles peuvent comporter des discontinuités de tout ordre
c Il y a lissage.

c entrées
c	n: nombre de fonctions
c	x, x1t, x2t, knot1, knot2: abcisses et vecteurs nodaux
c	m1, m2: ordres des splines
c	s1: première spline

c sorties:
c	s2: nouvelle spline

c Auteur: P.Morel, Département Cassiopée, O.C.A.

c----------------------------------------------------------------

	USE mod_kind
		
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION (:) :: x, x2t
	INTEGER, INTENT(in) :: knot2, m1, m2, n
	REAL (kind=dp), INTENT(inout), DIMENSION (:,:) :: s1
	REAL (kind=dp), INTENT(inout), DIMENSION (:) :: x1t
	INTEGER, INTENT(inout) :: knot1
	REAL (kind=dp), INTENT(out), DIMENSION (:,:) :: s2
			
	REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a
	REAL(kind=dp), DIMENSION(n) :: df, f	
	REAL(kind=dp), DIMENSION(m2+1) :: q
	REAL(kind=dp), DIMENSION(m1) :: qg, wg
	REAL(kind=dp) :: nu		
	
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc
	INTEGER, SAVE :: bloc, nx=0
	INTEGER :: col, colonne, dim, i, ig, j, k, l=1, l0=1, ligne, spi

	LOGICAL :: inversible
	
c----------------------------------------------------------------

2000	FORMAT(8es10.3)

c initialisations, bloc: longueur du bloc des coeff. non idt. nuls
	bloc=2*m2-1

c dim : dimension de l'espace des splines
c nl = rang : nombre de lignes du système linéaire
	dim=knot2-m2

c allocations et initialisations		
	ALLOCATE(a(dim,bloc),indpc(dim))
	a=0.d0 ; s2=0.d0
	 
c indices de première colonne pour les produits scalaires
c ils seront changés dans gauss_band
	DO i=1,dim      	!pour chaque spline de la base
	 indpc(i)=MAX(1,i-m2+1)
	ENDDO
	 
c---------------début de la construction du sytème linéaire------------	
	    
c les produits scalaires, pour chaque intervalle de m à dim
	ligne=0
	B1: DO k=m2,dim

c on évite les points multiples
	 IF(x2t(k) >= x2t(k+1))CYCLE B1	 
	  	 
c spi : indice-1 de la première spline non id.nulle  
c comme x2t(l) <= qg(ig) < x2t(l+1), l est le même pour toutes
c les abscisses de Gauss dans l'intervalle [x2t(k),x2t(k+1)]
	 CALL linf(x2t(k),x2t,knot2,l) ; spi=l-m2	  

c poids, wg et abscisses, qg pour intégration Gauss d'ordre 2m1-1
	 CALL intgauss(x2t(k),x2t(k+1),qg,wg,m1)

c pour chaque abscisse de GAUSS
	 DO ig=1,m1

c variables, nx=0 ne sert pas car .TRUE.
	  nu=MAX(x1t(1),MIN(qg(ig),x1t(knot1)))	   
	  CALL bsp1dn(n,s1,x,x1t,nx,m1,knot1,.TRUE.,nu,l0,f,df)
	  	   		   
c les B-splines en qg(ig), dérivées 0 à 1
	  CALL bval0(qg(ig),x2t,m2,l,q)   
   	   
c contribution au système linéaire
	  DO i=1,m2		!pour chaque spline d'indice i
	   ligne=spi+i
	   s2(:,ligne)=s2(:,ligne)+wg(ig)*q(i)*f(:)
		
c la matrice compressée est le jacobien 'diagonal' ie. sans les
c éléments 'non diagonaux' identiquement nuls
	   DO j=1,m2		!pour chaque spline j
	    colonne=spi+j 
	    col=colonne-indpc(ligne)+1   !colonne matrice compressée	     
	    a(ligne,col)=a(ligne,col)+wg(ig)*q(i)*q(j)
	   ENDDO	!j
	  ENDDO		!i
	 ENDDO 		!ig
	ENDDO B1
		 
c----------------fin de la construction du système linéaire-------------- 

c	 DO i=1,nl
c	  PRINT*,i,indpc(i)
c	  WRITE(*,2000)a(i,:),b(1,i)
c	 ENDDO
c	 PAUSE'a'
	 	
c résolution du système linéaire
	CALL gauss_band(a,s2,indpc,dim,dim,bloc,n,inversible)
	IF(.NOT.inversible)THEN
	 PRINT*,'ARRRET, matrice singulière dans newspl_gal' ; STOP
	ENDIF	 
c	WRITE(*,2000)b(1,:) ; PAUSE'solution'
	
c sortie	
	DEALLOCATE(a,indpc)
	
	RETURN
	
	END SUBROUTINE newspl_gal
