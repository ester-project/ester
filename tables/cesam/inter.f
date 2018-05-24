
c********************************************************************

	SUBROUTINE inter(m23_ou_r2,bp,q,qt,n,knot,x,f,dfdx,r2,m23)

c	routine public du module mod_variables
	
c	en (x=m**23 ou x=r2) ou bien (x=m ou x=r)
c	on fait une interpolation inverse x-->q-->f
c	valable pour variables lagrangiennes alors m23=m^2/3, r2=r^2
c	ou eulériennes alors m23=m, r2=r

c	la localisation de x en indice q0 est faite :
c		soit dans m23 si m23_ou_r2 = 'm23',
c		soit dans r2 si m23_ou_r2 = 'r2 ' ATTENTION r2 + blanc

c	remplacement du Newton par une dichotomie 16 03 01 (P. Morel)
c	adaptation de l'algorithme de Brent Num. Receip. p.354

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	ATTENTION : il faut initialiser et utiliser les bons bp, m23 et r2 !!!

c	en_masse = .TRUE. variables lagrangiennes m23=m**23, r2=r**2
c	en_masse = .FALSE. variables eulériennes m23=m, r2=r

c entrées
c	m23_ou_r2 = 'm23' : en m23, = 'r2 ' : en r2
c	ou bien en m et r
c	ne, bp,  : nb de fonctions (6)
c	bp,q,qt,n,m,knot : spline a interpoler
c	x : point d'interpolation
c	r2, m23 : abscisses 

c sorties
c	f, dfdx : valeurs et dérivées/r2 ou m23
c	ATTENTION f(6)=q_int : point d'interpolation

c------------------------------------------------------------------------

	USE mod_donnees, ONLY: ne, ord_qs
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, linf

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: x
	INTEGER, INTENT(in) :: n		
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: m23, r2
	CHARACTER (len=3), INTENT(in) :: m23_ou_r2		
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: bp
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: q, qt
	INTEGER, INTENT(inout) :: knot
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: f, dfdx
		
	REAL (kind=dp) :: corr, qi, qs, q_int	
	INTEGER :: i, ntour, inc, l=1	
	LOGICAL :: iterate
	
c----------------------------------------------------------

2000	FORMAT(8es10.3)

	SELECT CASE(m23_ou_r2)
	CASE('m23')		!interpolation inverse en m^2/3 ou m
	 inc=5	!indice de l'inconnue pour l'interpolation inverse
	 IF(x <= m23(1))THEN
	  q_int=1.d0 ; iterate=.FALSE.
	 ELSEIF(x >= m23(n))THEN
	  q_int=n ; iterate=.FALSE.	  
	 ELSE
	  CALL linf(x,m23,n,l)		!entre q_int et q_int+1
	  q_int=l	  
	  IF(x == m23(l))THEN
	   iterate=.FALSE.
	  ELSEIF(x == m23(l+1))THEN
	   q_int=l+1 ; iterate=.FALSE.
	  ELSE
	   iterate=.TRUE.
	  ENDIF	
	 ENDIF

	CASE('r2 ')			!interpolation inverse en r^2 ou r
	 inc=3
	 IF(x <= r2(1))THEN
	  q_int=1.d0 ; iterate=.FALSE.	  
	 ELSEIF(x >= r2(n))THEN
	  q_int=n ; iterate=.FALSE.	  
	 ELSE
	  CALL linf(x,r2,n,l)		!entre q_int et q_int+1
	  q_int=l	  
	  IF(x == r2(l))THEN
	   iterate=.FALSE.
	  ELSEIF(x== r2(l+1))THEN
	   q_int=l+1 ; iterate=.FALSE.
	  ELSE   
	   iterate=.TRUE.
	  ENDIF
	 ENDIF
	 
	CASE DEFAULT
	 WRITE(*,1)m23_ou_r2 ; WRITE(2,1)m23_ou_r2 ; CALL sortie
1	 FORMAT('ARRET appel a inter avec m23_ou_r2=',a3)
	END SELECT
	
c détermination de l'indice par algorithme de Brent pour interpolation inverse
	IF(iterate)THEN		!dichotomie
	 l=q_int ; corr=1 ; ntour=0
	 qi=q_int	!q_int inferieur
	 qs=q_int+1	!q_int superieur
	 DO WHILE(ntour < 30 .AND. corr > 1.d-6)
	  q_int=(qs+qi)/2.d0
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qs,knot,.TRUE.,q_int,l,f,dfdx)
	  IF(f(inc) <= x)THEN
	   qi=q_int
	  ELSE
	   qs=q_int
	  ENDIF
	  corr=qs-qi	  	  	 	 
	  ntour=ntour+1
c	  WRITE(*,2000)corr,qs,qi,x,f(inc)
	 ENDDO
	 IF(ntour >= 30)THEN
	  WRITE(*,"('pas de conv. dicho. dans inter en',a3)")m23_ou_r2
	  WRITE(*,"('q_int=',es10.3,', corr=',es10.3)")q_int,corr,x,f(inc)
	 ENDIF
	ELSE 
	 CALL bsp1dn(ne,bp,q,qt,n,ord_qs,knot,.TRUE.,q_int,l,f,dfdx)
	ENDIF	!iterate
c	PRINT*,ntour,iterate ; WRITE(*,2000)q_int,corr,x,f(inc),x-f(inc)
c	PAUSE

c dérivées par rapport à la variable d'interpolation inverse
c inc=3 R^2 ou R, inc=5 M^2/3 ou M	
	DO i=1,ne
	 IF(i /= inc)dfdx(i)=dfdx(i)/dfdx(inc)
	ENDDO
	
c astuce pour passer l'indice f(6) est fonction de répartition
	f(6)=q_int
	
	RETURN
	
	END SUBROUTINE inter
