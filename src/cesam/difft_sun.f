
c---------------------------------------------------------------

	SUBROUTINE difft_sun(melange,deff,gradad,gradrad,m,d)

c routine private du module mod_evol

c coefficient de diffusion turbulente, sous la ZC externe solaire
c suivant M. Gabriel 1997, A&A 327, 771
c d_turb + nu_gab + Deff

c Dimensions et initialisations dans le programme appelant
c d(nchim,nchim), dd(nchim,nchim,nchim), v(nchim),dv(nchim,nchim)

c convention de notation :
c équation de diffusion dXi/dt=dFi/dm + nuclear, i=1,nchim
c Fi=4 pi r**2 ro (4 pi r**2 ro D.dX/dm - Vi Xi)

c d=D=(di1, di2,... din) avec Dij coefficient de d Xj / dm
c dans le produit scalaire D.dX/dm=sum d_ij d Xj / dm

c entrées
c	melange=.TRUE.: on est dans une ZC
c	xchim
c	deff : diffusivité turbulente due à la rotation
c	gradad, gradrad : gradients
c	m : abscisse lagrangienne

c sorties
c	d: coefficients d_ij de d x_j / d m

c Auteur: P.Morel, OCA
c CESAM2k

c-------------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, d_conv, d_turb, langue, nchim
	USE mod_kind
	USE mod_variables, ONLY : age, lim, mstar, m_zc

	IMPLICIT NONE
      
	REAL (kind=dp), INTENT(in) :: deff, gradad, gradrad, m
	LOGICAL, INTENT(in) :: melange	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d
	REAL (kind=dp), SAVE :: log25
	REAL (kind=dp) :: la, lf, ld, lq, nu_gab, q, qstar
	INTEGER :: i

	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 WRITE(2,*)
	 log25=LOG10(2.5d0)
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1010)d_conv,d_turb
	  WRITE(2,1010)d_conv,d_turb
1010	  FORMAT('Solar turbulent diffusion of M.Gabriel. In CZ, Dconv=',es10.3,
	1 ', in RZ, Dturb=',es10.3,' + Nu_gab + Deff with rotation')	 
	 CASE DEFAULT	 
	  WRITE(*,10)d_conv,d_turb ; WRITE(2,10)d_conv,d_turb
10	  FORMAT('Diffusion turbulente solaire de M.Gabriel. Dans ZC, Dconv=',
	1 es10.3,', dans ZR, Dturb=',es10.3,' + Nu_gab + Deff avec rotation')
	 END SELECT
	ENDIF

c dans une zone de mélange	
	IF(melange)THEN
	 DO i=1,nchim
	  d(i,i)=d_conv
	 ENDDO
	ELSE
	
c coefficient de diffusivité de M.Gabriel
	q=m/mstar
	qstar=m_zc(lim)/mstar-0.01d0
	la=log10(age)-3.d0
	lq=(q-qstar)*10.d0**(1.16d0+0.48d0*la)
	ld=3.42d0-0.95d0*la-log25	
	lf=LOG10(5.d0*MIN(0.2d0,MAX(gradad-gradrad,0.15d0)))
	nu_gab=10.d0**(lf+ld+lq)
	 
c contributions des diverses diffusivités turbulentes
	 DO i=1,nchim
	  d(i,i)=d(i,i)+d_turb+nu_gab+deff
	 ENDDO
	ENDIF
	
	RETURN

	END SUBROUTINE difft_sun
