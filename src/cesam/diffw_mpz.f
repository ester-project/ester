
c*************************************************************

	SUBROUTINE diffw_mpz(dlnro,grand_k,nu,nu_cin,n2mu,n2t,r,ro,y,
	1 deff,dh,dv)

c routine private du module mod_evol

c calcul des coefficients de diffusion turbulente pour la rotation
c Dh, Dv, Deff, suivant Matis, Palacio & Zahn A&A 2004
c en optional dérivées de Dh et Dv / Omega, U, dOmega, dU

c entrées :
c	dlnro : d ln ro / dnu
c	grand_k : K
c	nu : m^2/3
c	nu_cin : viscosité cinématique
c	n2mu : Vaîssala mu
c	n2t : Vaîssala température
c 	r : rayon

c	y(:,0:1) : variables propres à la rotation et dérivées / nu

c sorties :
c	Deff, Dh, Dv : notations évidentes

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY :  langue, msol, nom_diffw, pi, rsol 
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(2,0:1) :: y	
	REAL (kind=dp), INTENT(in) :: dlnro, grand_k, nu, nu_cin, n2mu, n2t,
	1 r, ro
	REAL (kind=dp), INTENT(out) :: deff, dh, dv
	
	REAL (kind=dp), PARAMETER :: betta=1.5d-6, ric=1.d0/6.d0
	REAL (kind=dp), SAVE :: cte1_0, cte3_0, cte5_0, cte6_0, cts1_0, cts2_0
	REAL (kind=dp) :: cte1, cte2, cte3, cte5, den, dh2, gksdh, nu12
	
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
	 
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001)nom_diffw ; WRITE(2,1001)nom_diffw
1001	  FORMAT('Diffusion of angular momentum according to mpz',a,/,
	1 'Matis et al. A&A ')
	 CASE DEFAULT
	  WRITE(*,1)nom_diffw ; WRITE(2,1)nom_diffw
1	  FORMAT('Coefficients de diffusion du moment cinétique de mpz',a,/,
	1 'Matis et al. A&A 425, 243, 2004')
	 END SELECT
	 
	ENDIF
	
c initialisations
	nu12=SQRT(ABS(nu))		
	cte1=cte1_0*r**6*ro/nu12	
	cte2=(cts1_0*r**3*ro/nu12*dlnro-cts2_0)*r**3
	cte3=cte3_0*r**6*ro/nu12

c Dh
	dh2=cte1*y(1,0)*y(2,1)+cte2*y(1,0)*y(2,0)-cte3*y(1,1)*y(2,0)
	IF(SQRT(ABS(dh2)) <= nu_cin)THEN	
c	 PRINT*,'nu,nu_cin,dh2'; WRITE(*,2000)nu,nu_cin,dh2	
	 dh=nu_cin ; dv=nu_cin ; deff=nu_cin ; RETURN
	ELSE
	 IF(dh2 < 0.d0)THEN
	  dh2=-dh2
	 ENDIF
	 dh=SQRT(dh2)

c Dv  	  	
	 gksdh=1.d0+grand_k/dh
	 cte5=cte5_0*r**6*ro**2/nu	 
	 den=n2t+n2mu*gksdh 
	 dv=cte5*dh*gksdh*y(1,1)**2/den	 
c	 PRINT*,'nu,nu_cin,dh,dv'; WRITE(*,2000)nu,nu_cin,dh,dv	 
	 IF(dv > nu_cin)THEN
	  dv=dv+nu_cin
	 ELSE
	  dv=nu_cin
	 ENDIF

c Deff
	 deff=cte6_0*(r*y(2,0))**2/dh+dv

	ENDIF		!dh <= nu_cin

	RETURN
	
	END SUBROUTINE diffw_mpz
