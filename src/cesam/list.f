		
c********************************************************************	

	SUBROUTINE list(alfa,anub8,anube7,anun13,anuo15,anupep,anupp,beta,
	1 compg,cp,delta,dcapdr,dcapdt,depsdr,depsdt,d2p,d2ro,
	2 chaine,convec,ecritout,epsilon,gamma,gamma_atm,gradad,grada_atm,
	3 gradconv,gradc_atm,gradrad,gradr_atm,hp,i_cno,i_gr,i_pp,i_3a,
	4 kap,l,m,mu,mue,m_atm,p,pt,pt_atm,p_atm,r,ro,ro_atm,r_atm,t,tau,
	5 teff,tx_ioni,t_atm,u,vaissala,w,z,degene)

c	routine private du module mod_cesam
C
c	LISTING des résultats
C
c	Auteur : P.Morel, Département J.D. Cassini, O.C.A.
C	CESAM2k

c	07 03 97 : écriture de "Z"
c	18 08 97 : correction de l'écriture de la gravité (N.Audard)
c	05 11 97 : new taux de capture neutrinos solaires (G.Berthomieu)
c	20 06 99 : modif du taux de capture de B8 pour le Cl selon
c	Bahcall & Wolstein Phys.Rev. C 33, 2121, 1996 (G. Berthomieu)
c	23 03 00 : impression de [Fe/H]
c	05 04 00 : mise en en tête des déplétions, Z/X, Fe/H etc...
c	25 04 00 : mise en en tête des abondances/masse
c	17 07 00 : écriture des abondances dans Z / masse
C
c	n_qs: nombre de couches, n_qs: raccord à l'enveloppe, 1: le centre
c	m(.) : masse/msol
c	pt(.) : Ptot
c	p(.) : Pgas
c	ro(.) : densité
c	t(.) : température
c	r(.) : rayon/rsoleil
c	l(.) : luminosité/Lsoleil
c	compg(nchim,.) : composition chimique par gramme
c	teff : température effective
c	mtot : masse totale initiale
c	alpha : l/Hp
c	nchim : nombre d'éléments chimiques
c	ecritout : .true. / .false. pour listing complet / partiel
c	pt_atm : Ptot dans l'atmosphère
c	p_atm : Pgas dans l'atmosphère
c	t_atm : température
c	m_atm : masse
c	tau : épaisseur optique
c	r_atm : rayon
c	n_atm : nombre de couches
c	epsilon : énergie thermonucléaire et gravitationnelle
c	u : énergie interne spécifique
c	degene: dégénérescence
c	mstar: masse au temps t+dt, avec perte de masse
c	w(.): rotation
c	degene(.): dégénérescence
c	z(.) : le Z

c----------------------------------------------------------------

	USE mod_donnees, ONLY : alpha, amu, fesh_sol, g, gmsol, ife56,
	1 ihe4, langue, msol, mtot, nchim, nom_elem, nucleo,
	2 n_atm, pi, rsol, x0, y0, z0, zsx_sol
	USE mod_kind	
	USE mod_variables, ONLY : age, mstar, n_qs, rstar

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: epsilon, compg
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: alfa, anub8, anube7,
	1 anun13, anuo15, anupep, anupp, beta, delta, cp, dcapdr, dcapdt,
	2 depsdr, depsdt, gamma, gamma_atm, gradad, gradconv, grada_atm,
	3 gradc_atm, gradrad, gradr_atm, hp, kap, l, m, mu, mue, m_atm,
	4 p, pt, pt_atm, p_atm, r, ro, ro_atm, r_atm, t, tau, t_atm, u,
	5 vaissala, w, z, degene		
	REAL (kind=dp), INTENT(in) :: d2p, d2ro, teff	
	INTEGER, INTENT(in) :: i_cno, i_gr, i_pp, i_3a	

	LOGICAL, INTENT(in), DIMENSION(:) :: convec
	LOGICAL, INTENT(in) :: ecritout

        CHARACTER (len=*), INTENT(in), DIMENSION(:,:) :: tx_ioni		
	CHARACTER (len=*), INTENT(in) :: chaine

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: r_son, v_son
	REAL (kind=dp), DIMENSION(nchim) :: x	
	REAL (kind=dp), SAVE :: cte2, cte3, cte4, fes, logg	
	REAL (kind=dp) :: a, anufp, anufq, anupps, anupeps, anub8s, anube7s,
	1 anun13s, anuo15s, anufs, cv, dm, dnu02, dnu13, fesh,
	2 nu0, nz, nm, yc, zsx
	INTEGER, PARAMETER :: ch_max=5		
	INTEGER, SAVE :: j3=0, j4=0, j12=0, j13=0, j14=0, j15=0, j16=0,
	1 j17=0, j18=0
c	j2=0, j7=0, j7b=0, 	
	
	INTEGER :: i, iln, j, n_tot
	LOGICAL, SAVE :: abon, init=.TRUE.	
	LOGICAL :: convp

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE. ; cte2=beta(1)		!instruction fictive
	 cte2=4.*pi*(1.5d13)**2 !4 pi (R orb. terre)**2, pour neutrinos
	 cte4=0.107d-6/cte2	!pour neutrino du B8 (G. Berthomieu)
	 fes=-LOG10(zsx_sol)	![Fe/H]_sol

c identification des indices de certains éléments pour
c l'écritures de rapports isotopiques
	 DO j=2,nchim
c	  IF(nom_elem(j) == ' H2 ')THEN
c	   j2=j	
c	  ELSEIF(nom_elem(j) == 'He3 ')THEN
	   j3=j	   
	  IF(nom_elem(j) == 'He3 ')THEN
	   j3=j   
	  ELSEIF(nom_elem(j) == 'He4 ')THEN
	   j4=j
c	  ELSEIF(nom_elem(j) == 'Li7 ')THEN
c	   j7=j	  
c	  ELSEIF(nom_elem(j) == 'Be7 ')THEN
c	   j7b=j	  
	  ELSEIF(nom_elem(j) == 'C12 ')THEN
	   j12=j
	  ELSEIF(nom_elem(j) == 'C13 ')THEN
	   j13=j
	  ELSEIF(nom_elem(j) == 'N14 ')THEN
	   j14=j
	  ELSEIF(nom_elem(j) == 'N15 ')THEN
	   j15=j
	  ELSEIF(nom_elem(j) == 'O16 ')THEN
	   j16=j
	  ELSEIF(nom_elem(j) == 'O17 ')THEN
	   j17=j
	  ELSEIF(nom_elem(j) == 'O18 ')THEN
	   j18=j	   
	  ENDIF
	 ENDDO
	 abon=j3*j4*j12*j13*j14*j15*j16*j17*j18 > 0
	 logg=LOG10(gmsol/rsol/rsol)	 
	ENDIF	!initialisation

	cte3=mstar*msol/2.d0	!pour neutrinos

c ecritout=.TRUE. on forme le listing complet,
c ecritout=.FALSE. on n'écrit que l'en-tête	
	IF(ecritout)THEN	!listing ou entête, suivant le cas

c nombre total de lignes avec atmosphère
	 n_tot=n_qs
	 IF(n_atm > 0)n_tot=n_tot+n_atm-1
	
c calcul de v_son et r_son		
	 ALLOCATE(r_son(n_tot),v_son(n_tot))
	 DO i=1,n_qs
	  v_son(i)=sqrt(gamma(i)*p(i)/ro(i)) ; r_son(i)=r(i)*rsol
	 ENDDO
	 IF(n_atm > 0)THEN
	  DO i=2,n_atm
	   v_son(n_qs+i-1)=sqrt(gamma_atm(i)*p_atm(i)/ro_atm(i))
	   r_son(n_qs+i-1)=r_atm(i)*rsol
	  ENDDO
	 ENDIF	
c Yc	
	 IF(nchim > 1)THEN
	  yc=compg(ihe4,1)
	 ELSE
	  yc=1.d0-compg(1,1)-z0
	 ENDIF	
	 WRITE(2,13)age*1.d-3,LOG10(teff),LOG10(l(n_qs)),rstar,
	1 logg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),
	2 compg(1,1),yc,teff,l(n_qs),
	3 d2p,d2ro,i_pp,i_cno,i_3a,i_gr,chaine
13	 FORMAT(/,'age(milliard années)=',es10.3,', Log(temp. eff.)=',
	1 es11.4,', Log(luminosité/Lsol)=',es11.4,', Rstar/Rsol=',
	2 es11.4,/,'Log10g=',es10.3,', Pc=',es10.3,
	3 ', Tc=',es10.3,', roc=',es10.3,', Xc=',es10.3,', Yc=',es10.3,
	4 ', Teff=',es10.3,', L/Lsol=',es10.3,/,	
	5 '(Rstar**2/p d2p/dr2)c=',es10.3,
	6 ', (Rstar**2/ro d2ro/dr2)c=',es10.3,/,
	7 'ePP/eNUC=',i3,'%, eCNO/eNUC=',i3,'%, e3a+C+O/eNUC=',i3,
	8 '%, eGRAV/eNUC=',i4,'%, ',a30,/)	
	WRITE(2,26)mstar/mtot-1.d0,mstar,mtot
	 
c abondances en surface
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(2,1034)(nom_elem(i),i=1,MIN(nchim,12))
1034	  FORMAT('Abundances, volume/mass/dex(12+Log10(Ni/Nh)) at surface :',
	1 /,12(3x,a4,3x))
	  WRITE(2,31)(compg(i,n_qs)*ro(n_qs)/amu/nucleo(i),
	1 i=1,MIN(nchim,12))
	  WRITE(2,31)(compg(i,n_qs),i=1,MIN(nchim,12))
	  WRITE(2,31)(12.d0-LOG10(compg(1,n_qs)/nucleo(1))
	1 +LOG10(compg(i,n_qs)/nucleo(i)),i=1,MIN(nchim,12))
 	 CASE DEFAULT
	  WRITE(2,34)(nom_elem(i),i=1,MIN(nchim,12))
34	  FORMAT('Abondances, volume/masse/dex(12+Log10(Ni/Nh)) en surface :',
	1 /,12(3x,a4,3x))
	  WRITE(2,31)(compg(i,n_qs)*ro(n_qs)/amu/nucleo(i),
	1 i=1,MIN(nchim,12))
31	  FORMAT(1x,12es10.3)
	  WRITE(2,31)(compg(i,n_qs),i=1,MIN(nchim,12))
	  WRITE(2,31)(12.d0-LOG10(compg(1,n_qs)/nucleo(1))
	1 +LOG10(compg(i,n_qs)/nucleo(i)),i=1,MIN(nchim,12))
	 END SELECT		 	
	 IF(nchim > 12)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(2,37)(nom_elem(i),i=13,nchim)	 	
	   WRITE(2,31)(compg(i,n_qs)*ro(n_qs)/amu/nucleo(i),i=13,nchim)
	   WRITE(2,31)(compg(i,n_qs),i=13,nchim)
	   WRITE(2,31)(12.d0-LOG10(compg(1,n_qs)/nucleo(1))
	1  +LOG10(compg(i,n_qs)/nucleo(i)),i=13,nchim)
	  CASE DEFAULT	 
	   WRITE(2,37)(nom_elem(i),i=13,nchim)
37	   FORMAT(1x,12(3x,a4,3x))	 	
	   WRITE(2,31)(compg(i,n_qs)*ro(n_qs)/amu/nucleo(i),i=13,nchim)
	   WRITE(2,31)(compg(i,n_qs),i=13,nchim)
	   WRITE(2,31)(12.d0-LOG10(compg(1,n_qs)/nucleo(1))
	1  +LOG10(compg(i,n_qs)/nucleo(i)),i=13,nchim)
	  END SELECT  		
	 ENDIF
	 
c abondances au centre
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(2,1038)(nom_elem(i),i=1,MIN(nchim,12))
1038	  FORMAT(/,'Abundances / mass at center',/,12(3x,a4,3x))
	  WRITE(2,31)(compg(i,1),i=1,MIN(nchim,12))	 
	 CASE DEFAULT
	  WRITE(2,38)(nom_elem(i),i=1,MIN(nchim,12))
38	  FORMAT(/,'Abondances / masse au centre',/,12(3x,a4,3x))
	  WRITE(2,31)(compg(i,1),i=1,MIN(nchim,12))
	 END SELECT	  	 
	 IF(nchim > 12)THEN	  
	  WRITE(2,37)(nom_elem(i),i=13,nchim)
	  WRITE(2,31)(compg(i,1),i=13,nchim)	 
	 ENDIF
	 
c si on suit Fe56	 
	 IF(ife56 > 0)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(2,1050)LOG10(compg(ife56,n_qs)/nucleo(ife56))-
	1  LOG10(compg(1,n_qs)/nucleo(1))-fesh_sol
1050	   FORMAT(/,'with the abundance of Fe56 : [Fe/H]=',es10.3)	 
	  CASE DEFAULT
	   WRITE(2,50)LOG10(compg(ife56,n_qs)/nucleo(ife56))-
	1  LOG10(compg(1,n_qs)/nucleo(1))-fesh_sol
50	   FORMAT(/,'avec l''abondance de Fe56 : [Fe/H]=',es10.3)	 
	  END SELECT	 
	 ENDIF	
	 	 
c rapports isotopiques en surface
	 IF(abon)THEN	 
	  WRITE(2,35)compg(j3,n_qs)/compg(j4,n_qs)*nucleo(j4)/nucleo(j3),
	1 compg(j13,n_qs)/compg(j12,n_qs)*nucleo(j12)/nucleo(j13),
	2 compg(j15,n_qs)/compg(j14,n_qs)*nucleo(j14)/nucleo(j15),	
	3 compg(j17,n_qs)/compg(j16,n_qs)*nucleo(j16)/nucleo(j17),
	4 compg(j18,n_qs)/compg(j16,n_qs)*nucleo(j16)/nucleo(j18)	
35	  FORMAT(/,'Rapports isotopiques en nombre à la surface :',/,
	1 'He3/He4=',es10.3,', C13/C12=',es10.3,
	2 ', N15/N14=',es10.3,', O17/O16=',es10.3,', O18/O16=',es10.3)
	 ENDIF

c le Z/X
	 IF(ihe4 <= 1)THEN
	  zsx=z0
	 ELSE
	  zsx=1.d0-compg(ihe4,n_qs)-compg(ihe4-1,n_qs)-compg(1,n_qs)
	 ENDIF
	 zsx=zsx/compg(1,n_qs) ; fesh=fes+LOG10(zsx)	![Fe/H]
		
	 IF(w(n_qs) > 0.d0)THEN
	  WRITE(2,32)z(n_qs),zsx,fesh,2.d0*pi/24.d0/3600.d0/w(n_qs),
	1 rstar*w(n_qs)*rsol*1.d-5
32	  FORMAT(/,'Eléments lourds en surface, Z :',es10.3,' Z / X :',
	1 es10.3,', [Z/X]=',es10.3,//,'Période de rotation=',es10.3,' jours,',
	2 ' vitesse de rotation à la surface=',es10.3,' km/s')
	 ELSE
	  WRITE(2,33)z(n_qs),zsx,fesh
33	  FORMAT(/,'Eléments lourds en surface, Z :',es10.3,', Z / X :',
	2 es10.3,', [Z/X]=',es10.3)
	 ENDIF
	
c abondances relatives en nombre dans Z
	 IF(nchim > 1)THEN
	  nz=0.d0 ; nm=0.d0
	  DO i=ihe4+1,nchim	!somme des abondances en nombre dans Z
	   x(i)=compg(i,n_qs)/nucleo(i) ; nz=nz+x(i) ; nm=nm+compg(i,n_qs)
	  ENDDO	
	  WRITE(2,40)(nom_elem(i),i=ihe4+1,nchim)
40	  FORMAT(/,'Abondances relatives (en nombre/masse) dans Z des éléments en surface :',
	1 //,1x,12(3x,a4,3x))
	  WRITE(2,31)(x(i)/nz,i=ihe4+1,nchim)
	  WRITE(2,31)(compg(i,n_qs)/nm,i=ihe4+1,nchim)	  
	 ENDIF
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(2,1042)
1042	  FORMAT(/,10x,'------------- the model -------------')	 
	 CASE DEFAULT	 
	  WRITE(2,42)
42	  FORMAT(/,10x,'------------- le modèle -------------')
	 END SELECT
	 
c écriture réduite
	ELSE
	 WRITE(2,14)age*1.d-3,LOG10(teff),LOG10(l(n_qs)),rstar,
	1 logg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2 1.d0-compg(1,1)-z0,teff,l(n_qs),i_pp,i_cno,i_3a,i_gr,chaine
14	 FORMAT(/,'age(milliard années)=',es10.3,', Log(temp. eff.)=',
	1 es11.4,', Log(luminosité/lsol)=',es11.4,', Rstar/Rsol=',
	2 es10.3,/,'Log10g=',es10.3,', Pc=',es10.3,
	3 ', Tc=',es10.3,', roc=',es10.3,', Xc=',es10.3,', Yc=',es10.3,
	4 ', Teff=',es10.3,', L/Lsol=',es10.3,/,
	5 ' ePP/eNUC=',i3,'%, eCNO/eNUC=',i3,'%, e3a+C+O/eNUC=',i3,
	6 '%, eGRAV/eNUC=',i4,'%, ',a30,/)
	 WRITE(2,26)mstar/mtot-1.d0,mstar,mtot
26	 FORMAT('Var. rel. de masse=',es10.3,', mstar=',es10.3,
	1 ', m_ini=',es010.3,/)

c abondances en surface
	 WRITE(2,34)(nom_elem(i),i=1,MIN(nchim,12))
	 WRITE(2,31)(compg(i,n_qs)*ro(n_qs)/amu/nucleo(i),
	1 i=1,MIN(nchim,12))
	 WRITE(2,31)(compg(i,n_qs),i=1,MIN(nchim,12))	  
	 WRITE(2,31)(12.d0-LOG10(compg(1,n_qs)/nucleo(1))
	1 +LOG10(compg(i,n_qs)/nucleo(i)),i=1,MIN(nchim,12))	 	
	 IF(nchim > 12)THEN
	  WRITE(2,37)(nom_elem(i),i=13,nchim) 	
	  WRITE(2,31)(compg(i,n_qs)*ro(n_qs)/amu/nucleo(i),i=13,nchim)
	  WRITE(2,31)(compg(i,n_qs),i=13,nchim)
	  WRITE(2,31)(12.d0-LOG10(compg(1,n_qs)/nucleo(1))
	1 +LOG10(compg(i,n_qs)/nucleo(i)),i=13,nchim)
	 ENDIF
	 
c abondances au centre
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(2,1038)(nom_elem(i),i=1,MIN(nchim,12))
	  WRITE(2,31)(compg(i,1),i=1,MIN(nchim,12))	 
	 CASE DEFAULT
	  WRITE(2,38)(nom_elem(i),i=1,MIN(nchim,12))
	  WRITE(2,31)(compg(i,1),i=1,MIN(nchim,12))
	 END SELECT	  	 
	 IF(nchim > 12)THEN	  
	  WRITE(2,37)(nom_elem(i),i=13,nchim)
	  WRITE(2,31)(compg(i,1),i=13,nchim)	 
	 ENDIF
	 	
c rapports isotopiques en surface	 	 	
	 IF(abon)THEN	 
	  WRITE(2,35)compg(j3,n_qs)/compg(j4,n_qs)*nucleo(j4)/nucleo(j3),
	1 compg(j13,n_qs)/compg(j12,n_qs)*nucleo(j12)/nucleo(j13),
	2 compg(j15,n_qs)/compg(j14,n_qs)*nucleo(j14)/nucleo(j15),
	3 compg(j17,n_qs)/compg(j16,n_qs)*nucleo(j16)/nucleo(j17),
	4 compg(j18,n_qs)/compg(j16,n_qs)*nucleo(j16)/nucleo(j18)
	 ENDIF

c le Z/X
	 IF(ihe4 <= 1)THEN
	  zsx=z0
	 ELSE
	  zsx=1.d0-compg(ihe4,n_qs)-compg(ihe4-1,n_qs)-compg(1,n_qs)
	 ENDIF
	 zsx=zsx/compg(1,n_qs) ; fesh=fes+LOG10(zsx)	![Fe/H]
		
	 IF(w(n_qs) > 0.d0)THEN
	  WRITE(2,32)z(n_qs),zsx,fesh,2.d0*pi/24.d0/3600.d0/w(n_qs),
	1 rstar*w(n_qs)*rsol*1.d-5
	 ELSE
	  WRITE(2,33)z(n_qs),zsx,fesh
	 ENDIF
	
c abondances relatives en nombre dans Z
	 IF(nchim > 1)THEN
	  nz=0.d0 ; nm=0.d0
	  DO i=ihe4+1,nchim	!somme des abondances en nombre dans Z
	   x(i)=compg(i,n_qs)/nucleo(i) ; nz=nz+x(i) ; nm=nm+compg(i,n_qs)
	  ENDDO
	
	  WRITE(2,40)(nom_elem(i),i=ihe4+1,nchim)
	  WRITE(2,31)(x(i)/nz,i=ihe4+1,nchim)
	  WRITE(2,31)(compg(i,n_qs)/nm,i=ihe4+1,nchim)	  
	 ENDIF
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(2,1041)
1041	  FORMAT(//,'******------- next time-step------******',/)	 
	 CASE DEFAULT	 
	  WRITE(2,41)
41	  FORMAT(//,'******------- pas temporel suivant------******',/)
	 END SELECT	
	 RETURN
	ENDIF		!ecritout

c suite du listing dans le cas où on "ecritout"
	convp=.false.	!pour repérer les zones convectives
	
c pour intégrer le flux de neutrinos		
	anupps=0. ; anupeps=0. ; anub8s=0. ; anube7s=0.
	anun13s=0. ; anuo15s=0. ;anufs=0.	!modif GB
	
c iln: décompte du nombre de couches écrites
c il y a écriture de l'entête toutes les ch_max couches
	iln=ch_max+1	!pour écrire l'entête au début
	DO i=n_qs,1,-1	!pour chaque couche
	 IF(i /= 1)THEN	!intégration de la production de neutrinos
	  dm=(m(i)-m(i-1))*mstar	!formule des trapèzes
	  anupps=anupps+(anupp(i-1)+anupp(i))*dm
	  anupeps=anupeps+(anupep(i-1)+anupep(i))*dm
	  anub8s=anub8s+(anub8(i-1)+anub8(i))*dm
	  anube7s=anube7s+(anube7(i-1)+anube7(i))*dm
	  anun13s=anun13s+(anun13(i-1)+anun13(i))*dm
          anuo15s=anuo15s+(anuo15(i-1)+anuo15(i))*dm
	  anufq=anupp(i)*1.102d-4*(1.+0.02d0*t(i)*1.d-6)*ro(i)*	!modif GB
	1 (1.d0+compg(1,i))/2.d0/dsqrt(t(i)*1.d-6)
	  anufp=anupp(i-1)*1.102d-4*(1.d0+0.02d0*t(i-1)*1.d-6)*ro(i-1)*
	1 (1.d0+compg(1,i-1))/2.d0/dsqrt(t(i-1)*1.d-6)
	  anufs=anufs+(anufp+anufq)*dm
	 ENDIF

c écriture des entêtes toutes les ch_max couches
	 IF(iln >= ch_max)THEN
	  WRITE(2,8)		!en tête
8	  FORMAT(/,' couche',t10,'m/Mstar',t21,'r/Rstar',t30,'pression',
	1 t42,'temp.',t50,'densité',t60,'L/Lsol',t69,'en. int.',t80,
	2 'gradrad',t91,'kappa',t99,'epsilon',t112,'TdS',t121,'eps3AL',/,
	3 t9,'1-m/Mstar',t20,'1-r/Rstar',t31,'mu_mol',t42,'mu_e',
	4 t52,'cv',t60,'dégené',t70,'gamma1',t80,'grad_cvt',
	5 t89,'(dK/dT)ro',t99,'desp/dT',t111,'epsPP',t123,'Z',/,
	6 t10,'alpha',t21,'delta',t32,'cp',t40,'ech.ht.p',t49,'Pgaz/Ptot',
	7 t60,'v. son',t69,'vaissala',t81,'gradad',t89,
	8 '(dK/dro)T',t99,'deps/dro',t111,'epsCNO',t121,'Omega')	  
	  WRITE(2,9)(nom_elem(j),j=1,MIN(nchim,12))
9	  FORMAT(7x,12(3x,a4,3x))
	  IF(nchim > 12)WRITE(2,9)(nom_elem(j),j=13,nchim)
	  iln=0
	 ENDIF
	 IF(convec(i) .NEQV. convp)THEN
	  IF(convec(i))THEN
	   WRITE(2,*)
	   WRITE(2,*)'---------------début de zone convective------------'
	  ELSE
	   WRITE(2,*)
	   WRITE(2,*)'---------------fin de zone convective--------------'
	  ENDIF
	 ENDIF
	 convp=convec(i)
	 cv=cp(i)-p(i)*delta(i)**2/ro(i)/t(i)/alfa(i)
	 
c iln augmente de 1 à chaque couche écrite
	 iln=iln+1	 
	 WRITE(2,12)i,m(i)/mstar,r(i)/rstar,pt(i),t(i),ro(i),
	1 l(i),u(i),gradrad(i),kap(i),epsilon(1,i),epsilon(5,i),
	2 epsilon(4,i),1.d0-m(i)/mstar,1.d0-r(i)/rstar,mu(i),mue(i),cv,
	3 degene(i),gamma(i),gradconv(i),dcapdt(i),depsdt(i),epsilon(2,i),
	4 z(i),alfa(i),delta(i),cp(i),hp(i),p(i)/pt(i),v_son(i),
	5 vaissala(i),gradad(i),dcapdr(i),depsdr(i),epsilon(3,i),w(i)
	 WRITE(2,16)(compg(j,i),j=1,MIN(nchim,12))
	 WRITE(2,160)(tx_ioni(j,i),j=1,MIN(nchim,12))
160	 FORMAT(7x,12a10)	 
	 IF(nchim > 12)THEN
	  WRITE(2,16)(compg(j,i),j=13,nchim)
	  WRITE(2,160)(tx_ioni(j,i),j=13,nchim)
	 ENDIF
12	 FORMAT(/,1x,i6,12es10.3,/,(7x,12es10.3))
16	 FORMAT(7x,12es10.3)
	ENDDO	!i

	WRITE(2,*)' ' ; WRITE(2,*)'------------fin du modele-------------'

	IF(n_atm > 1)THEN
	 WRITE(2,23)
23	 FORMAT(//,t40,'ATMOSPHERE',//,t6,'ep.opt.',3x,'pression',1x,
	1 'température',2x,'R/R* - 1',3x,'M/M* - 1',5x,'ro',5x,
	2 'Pgaz/Ptot',2x,'grad.rad.',3x,'grad.conv',2x,'grad.ad',
	3 3x,'convec',/)
	 DO i=1,n_atm
	  WRITE(2,24)tau(i),p_atm(i),t_atm(i),r_atm(i)/rstar-1.d0,
	1 m_atm(i)/mstar-1.d0,ro_atm(i),p_atm(i)/pt_atm(i),gradr_atm(i),
	2 gradc_atm(i),grada_atm(i),(gradr_atm(i) > grada_atm(i))
24	  FORMAT(1x,10es11.4,3X,l1)
	  IF(mod(i,5) == 0)WRITE(2,*)' '
	 ENDDO
	ENDIF

c calcul des flux de neutinos à une distance de 1UA
	anupps=anupps*cte3
	anupeps=anupeps*cte3
	anub8s=anub8s*cte3
	anube7s=anube7s*cte3
	anun13s=anun13s*cte3
	anuo15s=anuo15s*cte3
	anufs=anufs*cte3	!modif GB

c !!!!!!!attention décalage entre nom du flux et physique
c !!!!!!! exemple anun13s correspond à la réaction O15 et anuo15s au fluor
c!!!!!!!! anupeps correspond au neutrino du Be7......
c
c	WRITE(2,15)anupps,anupeps,anub8s,anube7s,anun13s,anuo15s,
c	1 anupps/cte2,anupeps/cte2,anub8s/cte2,anube7s/cte2,
c	2 anun13s/cte2,anuo15s/cte2
	WRITE(2,15)anupps,anufs,anupeps,anub8s,anube7s,anun13s, !modif GB
	1 anuo15s,anupps/cte2,anufs/cte2,anupeps/cte2,anub8s/cte2,
	2 anube7s/cte2,anun13s/cte2,anuo15s/cte2,anupps*11.72d-10/cte2,
	3 anufs*204.d-10/cte2,anupeps*71.7d-10/cte2,anub8s*2.4d-6/cte2,
	4 anube7s*60.4d-10/cte2,anun13s*113.7d-10/cte2,
	5 anuo15s*113.9d-10/cte2,anupps*0.d0/cte2,anufs*16.d-10/cte2,
	6 anupeps*2.38d-10/cte2,anub8s*1.14d-6/cte2,	!modif 06/1999
	7 anube7s*1.66d-10/cte2,anun13s*6.61d-10/cte2,
	8 anuo15s*6.67d-10/cte2
15	 FORMAT(/,t30,'N E U T R I N O S',//,			!modif GB
	1 t24,'anupp',t34,'anupep',t44,'anube7',t54,'anub8',t64,'anun13',
	2 t74,'anuo15',t84,'anuf17',/,1x,'production',9x,7es10.3,/
	3 1x,'flux sur terre',5x,7es10.3,/1x,'pour le gallium',4x,
	4 7es10.3,/1x,'pour le chlore',5x,7es10.3)
	
c flux en SNU pour le chlore et le gallium, sections efficaces selon
c	G.Berthomieu taux de capture pour Ga d'après BahCALL 1997
	WRITE(2,17)(anupps*11.72d-10+anupeps*71.7d-10+anub8s*2.4d-6+ !GB 
	1 anube7s*60.4d-10+anun13s*113.7d-10+anufs*204.d-10
	2 +anuo15s*113.9d-10)/cte2,
	3 (anupps*0.d0+anupeps*2.38d-10+anub8s*1.14d-6+	!modif 06/1999
	4 anube7s*1.66d-10+anun13s*6.61d-10 +anufs*16.d-10
	5 +anuo15s*6.67d-10)/cte2,anub8s*cte4
17	FORMAT(/,'Pour le Gallium:',es10.3,' SNU,',' pour le Chlore: ',
	1 es10.3,' SNU,', ' pour B8/modèle:',es10.3,/)
	
c estimation de Nu0 et dnu02 et dnu13
	CALL dnunl(r_son,v_son,n_tot,nu0,dnu02,dnu13,a)	
	WRITE(2,18)nu0,dnu02,dnu13,a
18	FORMAT('Valeurs approx. de Nu0, delta Nu02, delta Nu13 et A',/,
	1 'Nu0=',es10.3,', dnu02=',es10.3,', dnu13=',es10.3,', A=',es10.3)
	
	WRITE(2,19)
19	FORMAT(/,'---------------------------------------------------',/)

	DEALLOCATE(r_son,v_son)
	
	RETURN
	
	END SUBROUTINE list

