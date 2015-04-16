
c*******************************************************************

	SUBROUTINE eq_ini_rota4(fait,nu,as,bs)	
	
c routine PRIVATE du module mod_evol

c formation des équations à résoudre pour l'initialisation
c Omega= w_rot=cte, dOmega=0, 
c Les inconnues sont Omega(fictif), U, Tau, Phi et Pi vérifient leurs équations
c Psi, Lambda, Upsilon sont nuls
c le système est linéaire, as : coefficients des variables et dérivées

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c entrées :
c	fait=0 : point courant, fait=1 : extérieur
c	nu: point de calcul en m**2/3
c	dt : pas temporel

c sorties :	
c	as, ad : coefficients de la spline/dérivée
c	bs : second membre

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : g, lsol, msol, m_ch, nrot, nvth, ord_rot, pi,
	1 rsol, w_rot
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : knotc, knotr, mc, mct, mrot, mrott,
	1 n_ch, n_rot, rota, vth
				
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: nu
	INTEGER, INTENT(in) :: fait
	
	REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: as
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: bs
	REAL (kind=dp), DIMENSION(nvth) :: dfvth, fvth
	REAL (kind=dp), SAVE, DIMENSION(11) :: cte_0
	REAL (kind=dp), DIMENSION(11) :: fr
	REAL (kind=dp), SAVE, DIMENSION(2) :: cts_0
	REAL (kind=dp), DIMENSION(2) :: cts
	
	REAL (kind=dp), SAVE ::	cte1_0, cte2_0	
	REAL (kind=dp) :: cp, delta, dgrad, dlldlm, dlngrav, dlnro,
	1 grad, gradad, grav, hp, lum, nu12, nu32, Omega, p, r, ro, t

	INTEGER, SAVE :: i, l=1
		
	LOGICAL, SAVE :: init=.TRUE.
		
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c initialisations
	IF(init)THEN
	 init=.FALSE.
	 
c pour la gravité
	 cte1_0=g*msol/rsol**2		!cte3_0
	 cte2_0=4.d0*pi*rsol**3/msol	!cte5_0
	 
	 	 
c pour les fr (F)
	 cte_0(1)=8.d0*pi*rsol**3/3.d0/msol 	!cte_0(25)	 
	 cte_0(2)=1.d0				!cte_0(27)
	 cte_0(3)=4.d0*pi*g*rsol**3/3.d0	!cte_0(29)
	 cte_0(4)=-rsol/3.d0/pi/g		!cte_0(1)
	 cte_0(5)=2.d0*rsol/3.d0		!cte_0(3)
	 cte_0(6)=8.d0*rsol**2/3.d0/g/msol	!cte_0(5)
	 cte_0(7)=-1.d0/pi/g/rsol		!cte_0(6)
	 cte_0(8)=2.d0/rsol			!cte_0(9)
	 cte_0(9)=-16.d0*pi*rsol**2/3.d0/msol	!cte_0(10)
	 cte_0(10)=msol/lsol			!cte_0(16)
	 cte_0(11)=8.d0*pi*rsol**2/15.d0/msol	!cte_0(31)
	 
	 	 
c pour les F* 	 
	 cts_0(1)=9.d0*msol/4.d0/pi/rsol**3	!cts_0(10)
	 cts_0(2)=4.d0*pi*g*rsol		!cts_0(11)
	 
	ENDIF

c~~~~~~~~~~~~~~~~fin des initialisations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c les coefficients au temps t+dt	
	
c nu12 etc...
	nu12=SQRT(nu) ; nu32=nu12**3
	
c utilisation de la tabulation vth	
c fvth(1)=lnP, fvth(2)=lnT, fvth(3)=r**2, fvth(4)=l**2/3, fvth(5)=ln ro, 
c fvth(6)=cp, fvth(7)=delta, fvth(8)=gamma1, fvth(9)=ln µ, fvth(10)=ln kap
	CALL bsp1dn(nvth,vth,mc,mct,n_ch,m_ch,knotc,.TRUE.,
	1 MAX(mc(1),MIN(nu,mc(n_ch))),l,fvth,dfvth)
	IF(no_croiss)PRINT*,'Pb. at 2 in eq_ini_rota4'
	
c quantités directement obtenues de la tabulation
	cp=fvth(6)
	delta=fvth(7)
	dlnro=dfvth(5)
	lum=SQRT(ABS(fvth(4)))**3
	r=SQRT(ABS(fvth(3)))
	
c quantités dérivées de la tabulation
	grad=dfvth(2)/dfvth(1)
	dlldlm=dfvth(4)*nu/fvth(4)	!d ln L / d ln m
	p=EXP(fvth(1))
	ro=EXP(fvth(5))
	t=EXP(fvth(2))
	
c hp
	hp=r**2*p/ro/nu32/cte1_0

c gravité	
	grav=cte1_0*nu32/r**2
	dlngrav=cte2_0*r**3*ro/nu32-2.d0
	
c gradient adiabatique
	gradad=p*delta/ro/cp/t
	
c dgrad
	 dgrad=ABS(gradad-grad)
				
c les cts*
	cts(1)=nu12/r**3/ro		!F*10
	cts(2)=r*ro/grav*dlnro		!F*11
	cts=cts*cts_0

c les coefficients fr(1)=F25, fr(2)=F27, fr(3)=F29, fr(4)=F1, fr(5)=F3,
c  fr(6)=F5, fr(7)=F6, fr(8)=F9, fr(9)=F10, fr(10)=F16, fr(11)=F31

	fr(1)=r**3*ro/nu12				!F25
	fr(2)=cts(1)+cts(2)				!F27
	fr(3)=r**3*ro/grav*dlnro			!F29
	fr(4)=r/ro/nu/grav*(dlngrav-2.d0)		!F1
	fr(5)=r/nu/grav*(1.d0-dlldlm)*(dlngrav-2.d0)	!F3
	fr(6)=r**2/grav/nu32				!F5
	fr(7)=dlngrav/grav/ro/nu/r			!F6
	fr(8)=(1.d0-dlldlm)/grav/nu/r*dlngrav		!F9
	fr(9)=r**2*ro/grav/nu32*(1.d0-dlldlm)		!F10
	fr(10)=cp*t/lum/hp*nu12*dgrad			!F16	
	fr(11)=r**4*ro					!F31	
	fr=fr*cte_0

c Omega est fixé
	Omega=ABS(w_rot)
		
c mises à 0
	as=0.d0 ; bs=0.d0	

	SELECT CASE(fait)
	CASE(0)

c Omega	(pour faciliter l'algorithme)
	 bs(1)=Omega ; as(1,1,0)=1.d0
	 
c (F6 Omega^2 + F9)Phi + (F5 Omega^2 + F10)dPhi + F16 U =
c -(F1 Omega^2 + F3)Omega^2
	 bs(2)=-(fr(4)*Omega**2+fr(5))*Omega**2
	 as(2,2,0)=fr(10)
	 as(2,7,0)=(fr(7)*Omega**2+fr(8))
	 as(2,7,1)=(fr(6)*Omega**2+fr(9))
	 
c Psi, Lambda, Upsilon =0.d0		
	 bs(3:4)=0.d0
	 DO i=3,4
	  as(i,i,0)=1.d0
	 ENDDO

c Tau - F31 Omega U = 0
	bs(5)=0.d0 ; as(5,2,0)=-fr(11)*Omega ; as(5,5,0)=1.d0

c Upsilon
	bs(6)=0.d0 ; as(6,6,0)=1.d0		
	
c dPi - F27 Phi = F29 Omega^2
	 bs(7)=fr(3)*Omega**2 ; as(7,7,0)=-fr(2) ; as(7,8,1)=1.d0
	 
c Pi - Phi - F25 dPhi = 0
	 bs(8)=0.d0 ; as(8,7,0)=-1.d0 ; as(8,7,1)=-fr(1) ; as(8,8,0)=1.d0

c condition limite au centre
	CASE(1)
	
c Omega	(pour faciliter l'algorithme)
	 bs(1)=Omega ; as(1,1,0)=1.d0	 

c U, Psi =0
	 bs(2:3)=0.d0 ; as(2,2,0)=1.d0 ; as(3,3,0)=1.d0
	 
c Phi = 0	
	 bs(4)=0.d0 ; as(4,7,0)=1.d0
	 	
c condition limite en surface
	CASE(2)
	
c Lambda, Tau, Upsilon =0.d0
	 bs(1:3)=0 ; as(1,4,0)=1.d0 ; as(2,5,0)=1.d0 ; as(3,6,0)=1.d0

c Phi + F25/3 dPhi=0		
	 bs(4)=0.d0 ; as(4,7,0)=1.d0 ; as(4,7,1)=fr(1)/3.d0	 
	 
	END SELECT

	RETURN

	END SUBROUTINE eq_ini_rota4
