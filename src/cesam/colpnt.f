
c*****************************************************************

	 FUNCTION colpnt(i,l,m,x)

c fonction PRIVATE du module mod_numerique
	
c calcul de l'abscisses du i-ieme point de collocation
c dans l'intervalle ]x(l),x(l+1)[
c pour la résolution d'une équation différentielle d'ordre r 
c avec m=ordre des splines il y a m points de collocation
c ces points sont les abscisses des formules d'intégration de Gauss
c d'ordre m sur l'ouvert ]x(l) , xt(l)[
c d'après de Boor p.293

c Auteur: P.Morel, Département J.D. Cassini, O.C.A, CESAM2k

c---------------------------------------------------------------------

c entrées:
c	 i: indice du point de collocation

c----------------------------------------------------------------------

	 USE mod_kind

	 IMPLICIT NONE
	 
	 REAL (kind=dp), INTENT(in), DIMENSION(:) :: x	 
	 INTEGER, INTENT(in) :: i, l, m
	 
	 REAL (kind=dp) :: colpnt	 
	
	 REAL (kind=dp), SAVE, DIMENSION(8) :: rho
	 REAL (kind=dp) :: xm, xp
	 
	 INTEGER, SAVE :: m0=-100

c----------------------------------------------------------------------

c	 initialisation des rho si m /= m0

	 IF(m /= m0)THEN	
	  SELECT CASE(m) 
	  CASE(1)
	   rho(1)=0.d0
	  CASE(2)
	   rho(2)=.57735 02691 89626 d0
	   rho(1)=-rho(2)
	  CASE(3)
	   rho(3)=.77459 66692 41483 d0
	   rho(1)=-rho(3)
	   rho(2)=0.d0
	  CASE(4)
	   rho(3)=.33998 10435 84856 d0
	   rho(2)=-rho(3)
	   rho(4)=.86113 63115 94053 d0
	   rho(1)=-rho(4)
	  CASE(5)
	   rho(4)=.53846 93101 05683 d0
	   rho(2)=-rho(4)
	   rho(5)=.90617 98459 38664 d0
	   rho(1)=-rho(5)
	   rho(3)=0.d0
	  CASE(6)
	   rho(4)=.23861 91860 83197 d0
	   rho(3)=-rho(4)
	   rho(5)=.66120 93864 66265 d0
	   rho(2)=-rho(5)
	   rho(6)=.93246 95142 03152 d0
	   rho(1)=-rho(6)
	  CASE(7)
	   rho(5)=.40584 51513 77397 d0
	   rho(3)=-rho(5)
	   rho(6)=.74153 11855 99394 d0
	   rho(2)=-rho(6)
	   rho(7)=.94910 79123 42759 d0
	   rho(1)=-rho(7)
	   rho(4)=0.d0
	  CASE(8)
	   rho(5)=.18343 46424 95650 d0
	   rho(4)=-rho(5)
	   rho(6)=.52553 24099 16329 d0
	   rho(3)=-rho(6)
	   rho(7)=.79666 64774 13627 d0
	   rho(2)=-rho(7)
	   rho(8)=.96028 98564 97536 d0
	   rho(1)=-rho(8)
	  CASE default
	   PRINT*,'dans colpnt, m < 0 ou  m > 8'
	   PRINT*,'entrer les poids de gauss pour m=',m
	   stop
	  END SELECT
	  m0=m
	 ENDIF 

	 xp=(x(l+1)+x(l))/2.d0 ; xm=(x(l+1)-x(l))/2.d0
	 colpnt=xp+rho(i)*xm

	 RETURN

	 END FUNCTION colpnt
