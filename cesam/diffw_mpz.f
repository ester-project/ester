
c*************************************************************

	SUBROUTINE diffw_mpz(dlnro,grand_k,nu,n2mu,n2t,r,ro,dfrot,frot,
	1 deff,dh,dv)

c	routine private du module mod_evol

c	calcul des coefficients de diffusion turbulente pour la rotation
c	hve(1:2) = Dh, Dv, Deff et dérivées dhve(2,2,0:1) / w et U
c	suivant Matis, Palacio & Zahn A&A 2004

c entrées :
c	nu : m^2/3
c	frot, dfrot : variables propres à la rotation et dérivées / nu

c sorties :
c	hve(1:2) : Dh, Dv
c	dhve(2,2,0:1) : dérivées / w, U : dw, dU

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : d_conv, msol, nom_rot, nrot, pi, rsol 
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: dfrot, frot 	
	REAL (kind=dp), INTENT(in) :: dlnro, grand_k, nu, n2mu, n2t, r, ro
	REAL (kind=dp), INTENT(out) :: deff, dh, dv
	
	REAL (kind=dp), PARAMETER :: betta= 1.5d-6, ric=1.d0/6.d0	
	REAL (kind=dp), SAVE :: cte1_0, cte3_0, cte5_0, cte6_0, cts1_0, cts2_0
	REAL (kind=dp) :: cte1, cte2, cte3, cte5, cts1, cts2, gksdh, nu12

	LOGICAL, SAVE :: init=.TRUE.
			
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 cte1_0=8.d0*pi*rsol**6*betta/9.d0/msol
	 cte3_0=4.d0*pi*rsol**6*betta/3.d0/msol	 
	 cte5_0=64.d0*pi**2*ric*rsol**6/9.d0/msol**2
	 cte6_0=rsol**2/30.d0	 
	 cts1_0=8.d0*pi*rsol**6*betta/9.d0/msol
	 cts2_0=rsol**3*betta/3.d0
	ENDIF
	
	nu12=SQRT(ABS(nu))
	cts1=cts1_0*r**6*ro/nu12*dlnro
	cts2=cts2_0*r**3
		
	cte1=cte1_0*r**6*ro/nu12	
	cte2=cts1-cts2
	cte3=cte3_0*r**6*ro/nu12
	
c	Dh
	dh=cte1*frot(1)*dfrot(2)+cte2*frot(1)*frot(2)-cte3*dfrot(1)*frot(2)
	dh=MAX(SQRT(ABS(dh)),1.d-10)	! pour éviter dh <= 0

c	Dv  	  	
	gksdh=1.d0+grand_k/dh
	cte5=cte5_0*r**6*ro**2/nu  
	dv=cte5*dh*gksdh*dfrot(1)**2/(n2t+n2mu*gksdh)
	dv=MAX(ABS(dv),1.d-10)	! pour éviter dv <= 0

c	Deff
	deff=cte6_0*r**2*frot(2)**2/dh+dv
	
c	dh=1.d3 ; dv=1.d3 ; deff=1.d3 

	RETURN
	
	END SUBROUTINE diffw_mpz
