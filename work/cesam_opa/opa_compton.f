
c*********************************************************************

	SUBROUTINE opa_compton(xchim,t,ro,kap,dkapt,dkapro,dkapx)

c routine public du module mod_opa	
c Opacités free-free Compton utilisées pour T > 7d7K ~ 7Kev

c Utilisation de la tabulation de Cox & Guili §16.a

c entrées :
c        xchim(1)=X : comp. chim. par gramme
c        t : température K
c        ro : densité cgs
  
c sorties :
c        kappa : opacité gr / cm2)
c        dkapdt : kappa / d t
c        dkapdr : kappa / d densité              
c        dkapdx : kappa / d xchim(1)

c	 Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	 CESAM2k
	
c--------------------------------------------------------------------

	 USE mod_donnees, ONLY : amu, clight, echarg, eve, kbol, langue, me,
	1 nchim, nucleo, pi, zi
	 USE mod_kind
	 
	 USE mod_numerique, ONLY : bsp1dn

	 IMPLICIT NONE
   
         REAL (kind=dp), INTENT(in), DIMENSION(nchim) :: xchim
	 REAL (kind=dp), INTENT(in) :: t, ro
	 REAL (kind=dp), INTENT(out) ::	kap, dkapt, dkapro, dkapx

	 INTEGER, PARAMETER :: mkt=4, nkt=11
	 REAL (kind=dp), DIMENSION(1,nkt) :: g	 
	 REAL (kind=dp), SAVE, DIMENSION(nkt+mkt) :: at	  	 
	 REAL (kind=dp), SAVE, DIMENSION(nkt) :: a, g1	 
         REAL (kind=dp), DIMENSION(nchim) :: xh
	 REAL (kind=dp), DIMENSION(1) :: corr, dcorr
	 
	 REAL (kind=dp), SAVE :: cte1, cte2, cte3
	 REAL (kind=dp) :: kt	 
	 
	 INTEGER, SAVE :: knot, l=1
	 	 
	 LOGICAL, SAVE :: init=.TRUE.	 
	
c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 cte1=8.d0*pi/3.d0*((echarg/clight)**2/me)**2/amu
	 cte2=cte1*zi(1)**2/nucleo(1)
	 cte3=kbol/eve*1.d-3

c températures en Kev (1Kev ~ 10 T6)
	 a=(/ 1.d0, 2.d0, 4.d0, 6.d0, 9.d0, 14.d0, 20.d0, 30.d0, 50.d0,
	1 80.d0, 125.d0 /)
	 g1=(/ 1.d0, 0.95055d0,0.9044d0, 0.8626d0, 0.8087d0, 0.7279d0, 0.6525d0,
	1 0.5590d0, 0.4408d0, 0.3411d0, 0.2579d0 /)
	 g=reshape(g1,shape(g))

	 CALL bsp1dn(1,g,a,at,nkt,mkt,knot,.FALSE.,a(1),l,corr,dcorr)
	 SELECT CASE(langue)	  
	 CASE('english')
	  WRITE(*,1001)t ; WRITE(2,1001)t
1001	  FORMAT('Use of Compton opacity, at least once, T=',es10.3) 	 
	 CASE DEFAULT	 
	  WRITE(*,1)t ; WRITE(2,1)t
1	  FORMAT('Utilisation de l''opacité Compton au moins 1 fois, T=',es10.3)
	 END SELECT
	ENDIF

c < Z / A >
	xh=xchim*zi/nucleo
	
c diffusion Thomson
	kap=cte1*SUM(xh)
	
c dérivée / X 
	dkapx=cte2

c correction Compton
	kt=MIN(cte3*t,a(nkt))		!Kt en Kev	
	CALL bsp1dn(1,g,a,at,nkt,mkt,knot,.TRUE.,kt,l,corr,dcorr)

c dérivée T
	dkapt=kap*dcorr(1)*cte3

c diffusion Compton	
	kap=kap*corr(1)	

c dérivées ro	
	dkapro=0.d0

	RETURN
	
	END SUBROUTINE opa_compton
