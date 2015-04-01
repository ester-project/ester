
c********************************************************************

	SUBROUTINE inter_atm(m_ou_r,ne_atm,bp_atm,x_atm,xt_atm,n_atm,
	1 ord_atm,knot_atm,m_atm,r_atm,x,f,dfdx)
	
c	routine public du module mod_exploit
	
c	en x=m ou x=r
c	on fait une interpolation inverse x-->q-->f
c	valable pour variables lagrangiennes alors m23=m^2/3, r2=r^2
c	ou eulériennes alors m23=m, r2=r

c	la localisation de x en indice q0 est faite :
c		soit en m si m_ou_r = 'm',
c		soit dans r si m_ou_r = 'r'

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées
c	m23_ou_r2 = 'm' : en m, = 'r' : en r
c	ne_atm, bp_atm, x_atm, xt_atm, n_atm, ord_atm, knot_atm: spline
c	d'interpolation de l'atmosphère
c	m_atm, r_atm : masses et rayons de l'atmosphère
c	x: point d'interpolation


c sorties
c	f, dfdx : valeurs et dérivées/r ou m
c	ATTENTION f(ne_atm+1)=x_int : point d'interpolation

c------------------------------------------------------------------------

	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, linf

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: m_atm, r_atm
	REAL (kind=dp), INTENT(in) :: x
	INTEGER, INTENT(in) :: ne_atm, n_atm, ord_atm
	CHARACTER (len=1), INTENT(in) :: m_ou_r	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: bp_atm	
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: x_atm, xt_atm
	INTEGER, INTENT(inout) :: knot_atm			
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: f, dfdx
		
	REAL (kind=dp) :: corr, qi, qs, q_int	
	INTEGER :: i, ntour, inc, l=1	
	LOGICAL :: iterate
	
c----------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(m_ou_r == 'm')THEN	!on donne la masse
	 inc=5				!5: indice de la masse
	 IF(x <= m_atm(1))THEN
	  q_int=1.d0 ; iterate=.FALSE.
	 ELSEIF(x >= m_atm(n_atm))THEN
	  q_int=n_atm ; iterate=.FALSE.	  
	 ELSE
	  CALL linf(x,m_atm,n_atm,l)	!entre q_int et q_int+1
	  q_int=l	  
	  IF(x == m_atm(l))THEN
	   iterate=.FALSE.
	  ELSEIF(x == m_atm(l+1))THEN
	   q_int=l+1 ; iterate=.FALSE.
	  ELSE
	   iterate=.TRUE.
	  ENDIF	
	 ENDIF 
	ELSEIF(m_ou_r == 'r')THEN	!on donne le rayon
	 inc=3				!3: indice du rayon
	 IF(x <= r_atm(1))THEN
	  q_int=1.d0 ; iterate=.FALSE.	  
	 ELSEIF(x >= r_atm(n_atm))THEN
	  q_int=n_atm ; iterate=.FALSE.	  
	 ELSE
	  CALL linf(x,r_atm,n_atm,l)		!entre q_int et q_int+1
	  q_int=l	  
	  IF(x == r_atm(l))THEN
	   iterate=.FALSE.
	  ELSEIF(x == r_atm(l+1))THEN
	   q_int=l+1 ; iterate=.FALSE.
	  ELSE   
	   iterate=.TRUE.
	  ENDIF
	 ENDIF	 
	ELSE
	 WRITE(*,1)m_ou_r ; STOP
1	 FORMAT('ARRET appel a inter_atm avec m23_ou_r2=',a)
	ENDIF
	
c	détermination de l'indice par algorithme de Brent pour
c	interpolation inverse
		
	IF(iterate)THEN		!dichotomie
	 l=q_int ; ntour=0
	 qi=q_int	!q_int inferieur
	 qs=q_int+1	!q_int superieur
	 
	 B6: DO
	  q_int=(qs+qi)/2.d0
	  CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
	1 .TRUE.,q_int,l,f,dfdx)
	  IF(f(inc) <= x)THEN
	   qi=q_int
	  ELSE
	   qs=q_int
	  ENDIF
	  corr=qs-qi ; ntour=ntour+1
	  IF(ntour < 30 .AND. corr > 1.d-6)THEN
	   CYCLE B6
	  ELSE
	   EXIT B6
	  ENDIF
	 ENDDO B6 
	 
	 IF(ntour >= 30)THEN
	  WRITE(*,"('pas de conv. dicho. dans inter_atm en',a)")m_ou_r
	  WRITE(*,"('q_int=',es10.3,', corr=',es10.3)")q_int,corr,x,f(inc)
	 ENDIF
	ELSE 
	 CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
	1 .TRUE.,q_int,l,f,dfdx)
	ENDIF	!iterate

c	IF(iterate)THEN
c	 PRINT*,ntour ; WRITE(*,2000)q_int,corr,x,f(inc),x-f(inc)
c	ELSE
c	 WRITE(*,2000)q_int,x,f(inc),x-f(inc)
c	ENDIF
c	PAUSE'inter_atm'
	
	DO i=1,ne_atm	!inc=3 ou 5 indice de l'inc. pour l'interp. inverse
	 IF(i /= inc)dfdx(i)=dfdx(i)/dfdx(inc)
	ENDDO
	
c	astuce pour passer l'indice
	
	f(ne_atm+1)=q_int
	
	RETURN
	
	END SUBROUTINE inter_atm
