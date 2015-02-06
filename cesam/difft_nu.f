
c---------------------------------------------------------------

	SUBROUTINE difft_nu(melange,t,ro,drox,kap,dkapx,deff,d,dd)

c	routine private du module mod_evol

c	calcul du coefficient de diffusion turbulente, d_turb + nu_rad
c	et du Moment Angulaire, sauf dans les ZC mélangées

c	Dimensions et initialisations dans le programme appelant
c	d(nchim,nchim), dd(nchim,nchim,nchim), v(nchim),
c	dv(nchim,nchim)

c	convention de notation :
c	équation de diffusion dXi/dt=dFi/dm + nuclear, i=1,nchim
c	Fi=4 pi r**2 ro (4 pi r**2 ro D.dX/dm - Vi Xi)

c	d=D=(di1, di2,... din) avec Dij coefficient de d Xj / dm
c	dans le produit scalaire D.dX/dm=sum d_ij d Xj / dm

c	pour ligne d'indice i
c	v(i) coefficient de x_i,
c	dv(i,k)=dv(nchim*(k-1)+i)=dérivée v_i / x_k
c	seule la première colonne de dv
c	est non nulle (pas de dérivées / Xi, i .ne. 1)
c	d(i,j)=coefficient d_ij de d x_j / dm
c	dd(i,j,k)= dérivée de d_ij / x_k

c	Auteur: P.Morel, OCA
c	CESAM2k

c entrées
c	melange=.TRUE.: on est dans une ZC
c	p, t, r, l, m, ro: données au point de calcul
c	xi: composition chimique, par mole
c	kap: opacité 
c	gradad, gradrad: gradients
c	terminaisons x : dérivées/ X1 (ie H) par gramme
c	mstar: masse avec perte de masse
c	m_zc, r_zc, lim : masses, rayons, nombre de limites ZR/ZC
c	age, gamma1, cp, delta: notations evidentes

c sorties
c	d0, dd : coefficients d_ij de d x_j / d m et dérivées / x_k
c	v0, dv : coefficients v_i de x_i et dérivées / x_k

c-------------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, d_conv, d_turb, Krot,
	1 langue, nchim, re_nu
	USE mod_kind

	IMPLICIT NONE
      
	REAL (kind=dp), INTENT(in) :: deff, dkapx, drox, kap, ro, t
	LOGICAL, INTENT(in) :: melange	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d
	REAL (kind=dp), SAVE :: cte2
	REAL (kind=dp) :: dnu_radx, nu_rad

	INTEGER :: i

	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 WRITE(2,*)
	 cte2=re_nu*aradia/clight*4.d0/15.d0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1010)d_conv,d_turb,re_nu
	  WRITE(2,1010)d_conv,d_turb,re_nu
1010	  FORMAT('Turbulent diffusion : Dturb + Dradiative',/,
	1  'In CZ, Dconv=',es10.3,'. In RZ, Dturb=',es10.3,
	2  ', Re_nu=',es10.3)	 
	 CASE DEFAULT	 
	  WRITE(*,10)d_conv,d_turb,re_nu ; WRITE(2,10)d_conv,d_turb,re_nu
10	  FORMAT('Diffusion turbulente: Dturb + Dradiative',/,
	1  'Dans ZC, Dconv=',es10.3,'. Dans ZR, Dturb=',es10.3,
	2  ', Re_nu=',es10.3,' + Deff')
	 END SELECT
	ENDIF
	
	IF(melange)THEN
	 DO i=1,nchim
	  d(i,i)=d_conv
	 ENDDO
	ELSE
	
c	 coefficient de diffusivité radiative
	 nu_rad=cte2*t**4/kap/ro**2
	 dnu_radx=-nu_rad*(dkapx/kap+2.d0*drox/ro)
	 
c	 contributions des diverses diffusivités turbulentes
	 DO i=1,nchim
	  d(i,i)=d(i,i)+d_turb+nu_rad+deff ; dd(i,i,1)=dd(i,i,1)+dnu_radx
	 ENDDO
	ENDIF
	
	RETURN

	END SUBROUTINE difft_nu
