
c******************************************************************

	SUBROUTINE diffw_p03(dlnro,grand_k,nu,nu_cin,n2mu,n2t,r,ro,y,
	1 deff,dh,dv)

c routine private du module mod_evol

c calcul des coefficients de diffusion turbulente pour la rotation
c Dh, Dv, Deff, suivant Palacios & al A&A 2003

c entrées :
c	dlnro : d ln ro / dnu
c	grand_k : K
c	nu : m^2/3
c	nu_cin : viscosité cinématique
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

	REAL (kind=dp), PARAMETER :: ch=1.d0, ric=2.d0/5.d0
	REAL (kind=dp), SAVE :: cte1_0, cte2_0, cte3_0, cte5_0, cte6_0,
	1 cts1_0, cts2_0
	REAL (kind=dp) :: cte1, cte2, cte3, cte5, cts1, cts2, den, gksdh, nu12

	LOGICAL, SAVE :: init=.TRUE.

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 cts1_0=8.d0*pi*rsol**4/9.d0/msol
	 cts2_0=rsol/3.d0

	 cte1_0=1.d0/ch
	 cte2_0=8.d0*pi*rsol**4/9.d0/msol
	 cte3_0=4.d0*pi*rsol**4/3.d0/msol
	 cte5_0=64.d0*pi**2*ric*rsol**6/9.d0/msol**2
	 cte6_0=rsol**2/30.d0
	 
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001)nom_diffw ; WRITE(2,1001)nom_diffw
1001	  FORMAT('Diffusion of angular momentum according to p03',a,/,
	1 'Palacios et al. 399, 603, 2003')
	 CASE DEFAULT
	  WRITE(*,1)nom_diffw ; WRITE(2,1)nom_diffw
1	  FORMAT('Coefficients de diffusion du moment cinétique de p03',a,/,
	1 'Palacios et al. 399, 603, 2003')
	 END SELECT
	 
	ENDIF

c initialisations
	nu12=SQRT(ABS(nu))
	cts1=cts1_0*r**4/nu12
	cts2=cts2_0*r

	cte1=(cts1-cts2)*r**4*ro*dlnro/nu12*cte1_0
	cte2=cte2_0*r**4*ro/nu12*cte1_0
	cte3=cte3_0*r**4*ro/nu12*cte1_0

c Dh
	dh=cte1*y(2,0)+cte2*y(2,1)-cte3*y(2,0)*y(1,1)/y(1,0)
	dh=ABS(dh)

	IF(dh <= nu_cin)THEN
	 dh=nu_cin ; dv=nu_cin ; deff=nu_cin ; RETURN
	ELSE

c Dv
	 gksdh=1.d0+grand_k/dh
	 cte5=cte5_0*r**6*ro**2/nu
	 den=n2t+n2mu*gksdh
	 dv=cte5*dh*gksdh*y(1,1)**2/den
	 IF(dv > nu_cin)THEN
	  dv=dv+nu_cin
	 ELSE
	  dv=nu_cin
	 ENDIF

c Deff
	 deff=cte6_0*(r*y(2,0))**2/dh+dv

	ENDIF		!dh <= nu_cin

	RETURN

	END SUBROUTINE diffw_p03
