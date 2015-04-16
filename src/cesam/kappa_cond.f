
c****************************************************************

	SUBROUTINE kappa_cond(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)

c	routine private du module mod_opacite

c	ajout de l'opacite conductive

c	adaptation de la partie concernée du code de Genève procurée par
c	Y. Lebreton et selon ses instructions

c	expression of Iben(1975  AP.j. 196,525 Appendix A) and numerical
c	evaluations of the derivatives with respect to density and
c	temperature assuming dlro=0.001 and dlt=0.0001. Derivatives
c	are calculated assuming vmye=2./(1.+x)=constante=nb. nucleons par
c	electron libre

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées:
c	xh(1)=X : comp. chim. en fraction de masse
c	t : température K
c	ro : densité cgs

c entrées/sorties :
c	kapa : opacité gr / cm2)
c	dkapdt : kappa / d t
c	dkapdr : kappa / d densité              
c       dkapdx : kappa / d xchim(1)

c--------------------------------------------------------------------------

	USE mod_donnees, ONLY : ihe4, langue, nchim, nom_elem, pi, z0
	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xh
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(inout) :: kappa, dkapdt, dkapdr, dkapdx

	REAL (kind=dp), PARAMETER, DIMENSION(8) ::
	1 zr=(/ 1.d0, 2.d0, 6.d0, 7.d0, 8.d0, 10.d0, 10.d0, 12.d0/),
	2 a=(/1.d0, 4.d0, 12.d0, 14.d0, 16.d0, 20.d0, 22.d0, 24.d0/)	
	REAL (kind=dp), DIMENSION(8) :: xx
	REAL (kind=dp), PARAMETER :: dx=0.00001d0, unpdx=1.d0+dx

	REAL (kind=dp) :: vmye, vmye_log, za, za_log, zb, xz2, xhh,
	1 rol6, tl6, dellg, del, eta0, eta2, t1, ro1,
	2 a1, a2, b1, b2, rnedne, flg, blamr2, alfa, zc, f, thxlg,
	3 thx, penktl, ef, efm, gam, glg, vkcc, thy, thc, thylg, thclg,
	4 vcond, vkchr, kappa_c, dkapdr_c, dkapdt_c, dkapdx_c,
	5 store, store0, dstore

	INTEGER, SAVE :: j12=0, j14=0, j16=0, j20=0, j22=0, j24=0
	INTEGER ::  j, i, k
	
	LOGICAL, SAVE :: init=.TRUE.
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(nchim <= 3)RETURN !pas d'opacités conductives

	IF(init)THEN			!indices des elements

c	 WRITE(*,2000)xh(1:nchim),t,ro ; WRITE(*,*)nom_elem(1:nchim)

	 init=.FALSE.
	
	 DO j=3,nchim
	  IF(nom_elem(j) == 'C12 ')THEN
	   j12=j
	  ELSEIF(nom_elem(j) == 'N14 ')THEN
	   j14=j
	  ELSEIF(nom_elem(j) == 'O16 ')THEN
	   j16=j
	  ELSEIF(nom_elem(j) == 'Ne20')THEN
	   j20=j
	  ELSEIF(nom_elem(j) == 'Ne22')THEN
	   j22=j
	  ELSEIF(nom_elem(j) == 'Mg24')THEN
	   j24=j
	  ENDIF
	 ENDDO
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1001) ; WRITE(2,1001)	 
1001	  FORMAT(/,'Conductives Opacities (Iben Apj 196,515,1975)')	 
	 CASE DEFAULT	  	 
	  WRITE(*,1) ; WRITE(2,1)	 
1	  FORMAT(/,'Opacités conductives (Iben Apj 196,515,1975)')
	 END SELECT
	ENDIF

c	test whether electron conduction is likely to be important
c	calculation of conductive opacity according to the analytic

	IF(ro <= 1.d-4 .OR. t <= 2.d4)RETURN

c	expression of Iben(1975  AP.j. 196,525 Appendix A) and numerical
c	evaluations of the derivatives with respect to density and
c	temperature and chemical composition X.

	ro1=ro ; t1=t ; xhh=xh(1)
	DO k=1,4
	 IF(k == 2)THEN
	  store0=t1
	  store=store0*unpdx
	  dstore=store-store0
	  t1=store
	 ELSEIF(k == 3)THEN
	  store0=ro1
	  store=store0*unpdx
	  dstore=store-store0
	  ro1=store
	 ELSEIF(k == 4)THEN
	  store0=xh(1)
	  store=max(1.d-5,store0*unpdx)
	  dstore=store-store0
	  xhh=store
	 ENDIF

c 	 valeurs de T6 et ro6 en Log10

	 tl6=log10(t1)-6.d0
	 rol6=log10(ro1)-6.d0

c	 abondances en fraction de masse de H, He4, C12,
c	 N14, O16, Ne20, Ne22 et Mg24

	 xx=0.d0
	 xx(1)=xhh
	 xx(2)=1.d0-(xhh+SUM(xh(ihe4+1:nchim)))
	 IF(j12 >= 3)xx(3)=xh(j12)
	 IF(j14 >= 3)xx(4)=xh(j14)
	 IF(j16 >= 3)xx(5)=xh(j16)
	 IF(j20 >= 3)xx(6)=xh(j20)
	 IF(j22 >= 3)xx(7)=xh(j22)
	 IF(j24 >= 3)xx(8)=xh(j24)

c	 l'opacité conductive

	 vmye=2.d0/(1.d0+xx(1)) ; vmye_log=log10(vmye) ; za=0.d0 ; zb=0.d0
	 DO i=1,8
          xz2=xx(i)*zr(i)**2 ; za=za+xz2/(a(i)**(1.d0/3.d0))
	  zb=zb+xz2/a(i)
	 ENDDO
	 za_log=log10(za)

	 IF(rol6 <= 0.3d0) THEN		!Hubbart pour ro<2.d6
	  dellg=rol6+6.d0-1.5d0*tl6-vmye_log ; del=10.d0**(dellg)  !A1
	  eta0=10.d0**(-0.52255d0+2.d0*dellg/3.d0) ; eta2=eta0**2  !A6

c	  logique de A2 a A5

	  a1=-3.29243d0+log10(del*(1.d0+0.02804d0*del))		!A3
	  b1=-4.80946d0+log10(del*del*(1.d0+9.376d0/eta2))	!A4
	  IF(dellg <= 0.645d0) THEN
	   flg=-3.2862d0+log10(del*(1.d0+0.024417d0*del))	!A2	
	  ELSEIF(dellg <= 2.d0)THEN
	   flg=a1					!A3
	  ELSEIF(dellg <= 2.5d0) THEN
	   flg=2.d0*a1*(2.5d0-dellg)+2.d0*b1*(dellg-2.d0)	!A5
	  ELSE
	   flg=b1					!A4
	  ENDIF

c	  logique de A8 a A10

	  a2=log10(1.d0+0.021876d0*del)			!A8
	  b2=log10(0.4d0*eta0*(1.d0+4.1124d0/eta2))	!A9
	  IF(dellg <= 1.5d0)THEN	  
	   penktl=a2
	  ELSEIF(dellg <= 2.d0) THEN
	   penktl=2.d0*a2*(2.d0-dellg)+2.d0*b2*(dellg-1.5d0)	!A10
	  ELSE
	   flg=b1					!A9
	   penktl=2.d0*a2*(2.d0-dellg)+2.d0*b2*(dellg-1.5d0)	!à vérifier
          ENDIF

c	  logique de A11, A12

	  IF(del <= 40.d0)THEN
	   rnedne=1.0d0-0.01d0*del*(2.8966d0-0.034838d0*del)	!A11
	  ELSE
	   rnedne=(1.5d0/eta0)*(1.d0-0.8225d0/eta2)		!A12
	  ENDIF

c	  logique de A13 a A21

	  blamr2=9.24735d-3*10.d0**(dellg-0.5d0*tl6-penktl)	!A7
	1 *(vmye*zb+rnedne)                    !lamba ou lambda_barre 
	  alfa = log10(blamr2)                               ! corrigé 
	  IF(alfa <= -3.d0)THEN		!pour H
	   thxlg=1.048d0-0.124d0*alfa				!A13
	  ELSEIF(alfa <= -1.d0)THEN
	   thxlg=0.13d0-alfa*(0.745d0+0.105d0*alfa)		!A14
	  ELSE
	   thxlg=0.185d0-0.585d0*alfa			!A15 corrigée
	  ENDIF
	  IF(alfa <= -3.d0)THEN		!pour He4
	   thylg=0.937d0-0.111d0*alfa				!A16
	   ELSEIF(alfa <= 0.d0)THEN
	    thylg=0.24d0-alfa*(0.55d0+0.0689d0*alfa)		!A17
	   ELSE
	    thylg=0.24d0-0.6d0*alfa				!A18
	  ENDIF
	  IF(alfa <= -2.5d0)THEN				!pour C
	   thclg=1.27d0-0.1d0*alfa				!A19
	  ELSEIF(alfa <= .5)THEN
	   thclg=0.727d0-alfa*(0.511d0+0.0778d0*alfa)		!A20
	  ELSE
	   thclg=0.843d0-0.785d0*alfa				!A21
	  ENDIF
	  thx=10.d0**(thxlg) ; thy=10.d0**(thylg) ; thc=10.d0**(thclg)

c	  zc is an effective abundance of C,N,O,Ne weighted by the
c	  squares of the charges and normalised with respect to C since
c	  electron conduction has been calculated by Hubbard+Lampe only
c	  for pure H,pure He and pure C

	  zc=(3.d0*xx(3)+3.5d0*xx(4)+4.d0*xx(5)+5.d0*xx(6)+4.54d0*xx(7)+
	1 6.d0*xx(8))/3.d0 !A23
	  vkchr=(xx(1)*thx+xx(2)*thy+zc*thc)*10.d0**(-tl6-flg)	   !A22

	  IF(rol6 <= 0.d0)THEN	 !Hubbart en dessous de ro=1.d6
	   vcond=vkchr
	  ELSE	!Hubbart + Canuto pour 1.d6< ro < 2.d6
	   ef=sqrt(1.d0+10.d0**(2.d0*(rol6-vmye_log)/3.d0)) - 1.!A24
	   gam=22.76d0*10.d0**(rol6/3.d0-tl6)*za		!A25
	   efm=min(1.d0,0.5d0+log10(ef))			!A30
	   glg=(0.873d0-0.298d0*za_log+(0.333d0-0.168d0*za_log)*efm)* !A29
	1 (1.d0-(1.d0+gam)**(-0.85d0))
	   vkcc=6.753d-8*10.d0**(2.d0*tl6-glg)*zb/(ef**2.5d0*(1.d0+ef))!A28
	   f=0.5d0*(1.d0-cos(pi*rol6/0.3d0))			!A32
	   vcond=10.d0**((1.d0-f)*log10(vkchr)+f*log10(vkcc))	!A31
	  ENDIF
	 ELSE	!Canuto pour ro >2.d6
	  ef=sqrt(1.d0+10.d0**(2.d0*(rol6-vmye_log)/3.d0)) - 1.d0    !A24
	  gam=22.76d0*10.d0**(rol6/3.d0-tl6)*za			!A25
	  efm=min(1.d0,0.5d0+log10(ef))				!A30
	  glg=(0.873d0-0.298d0*za_log+(0.333d0-0.168d0*za_log)*efm)*   !A29
	1	(1.d0-(1.d0+gam)**(-0.85d0))
	  vkcc=6.753d-8*10.d0**(2.d0*tl6-glg)*zb/(ef**2.5d0*(1.d0+ef)) !A28
	  vcond=vkcc
	 ENDIF
	 IF(k == 1)THEN
	  kappa_c=vcond
	 ELSEIF(k == 2)THEN
	  t1=t ; dkapdt_c=(vcond-kappa_c)/dstore
	 ELSEIF(k == 3)THEN
	  ro1=ro ; dkapdr_c=(vcond-kappa_c)/dstore
	 ELSEIF(k == 4)THEN
	  xhh=xh(1) ; dkapdx_c=(vcond-kappa_c)/dstore
	 ENDIF
	ENDDO		!k

c	WRITE(*,2000)kappa,kappa_c,t,ro,xh(1)

	dkapdr=(kappa**2*dkapdr_c+kappa_c**2*dkapdr)/(kappa+kappa_c)**2
	dkapdt=(kappa**2*dkapdt_c+kappa_c**2*dkapdt)/(kappa+kappa_c)**2
	dkapdx=(kappa**2*dkapdx_c+kappa_c**2*dkapdx)/(kappa+kappa_c)**2
	kappa=kappa*kappa_c/(kappa+kappa_c)
	
	RETURN

	END SUBROUTINE kappa_cond
