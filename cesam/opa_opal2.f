
c***********************************************************************

	SUBROUTINE opa_opal2(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx,cno)

c	routine private du module mod_opa

c Calcul de l'opacité basée sur les tables d'opacité de type 2 de Livermore
c	Adaptation à CESAM2k de la routine d'interpolation xcotrin21 de
c	A.I. Boothroyd

c	S'il le néon est présent, appel à opal_x_cno_fu de façon à tenir
c	compte dans l'opacité des écarts d'abondances de C, N, O passés
c	dans le néon

c	Sans néon, appel à opac en tenant compte des variations de C et O
c	il n'y a pas interpolation en Z, mais correction pour C et O par
c	rapport au Z de la table OPAL utilisée

c	Auteur : P.Morel, Département Cassiopée, O.C.A.
c	Conseils : A.I. Boothroyd
c	CESAM2k	

c entrées :
c	xchim: comp. chim. en H + C & O
c	t6 : température millions de K
c	r  : densité cgs/(t6**3)
c	cno=.TRUE. : cas général avec la routine opal_x_cno_fu
c sorties :
c	kapa : opacité gr / cm2
c	dkapdt : kappa / d T
c	dkapdr : kappa / d densite              
c       dkapdx : kappa / d X

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY: abon_m, ab_min, ife56, langue, nchim,
	1 nom_abon, nom_elem, nom_xheavy
	USE mod_kind
	USE mod_variables, ONLY : sortie
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xh
	REAL (kind=dp), INTENT(in) :: t, ro
	LOGICAL, INTENT(in) :: cno	
	REAL (kind=dp), INTENT(out) :: dkapdt, dkapdr, dkapdx, kappa
	
	INTEGER, PARAMETER :: iulow=23, max_iso=5, nmet=19
	REAL (kind=sp), DIMENSION(nmet) :: xmet=0., xisz_absent=0.
	REAL (kind=sp), PARAMETER :: mfe=55.847, mo=15.9994, dx=0.001,
	1 unpdx=1.+dx
	REAL (kind=sp), SAVE :: c0, dslr_max, dslr_min, dslt_max, dslt_min,
	1 osfe_sol, o0, slr_max, slr_min, slt_max, slt_min, z_table
	REAL (kind=sp) :: ex_c, ex_o, fedge, fcnone_opal,
	1 fcn_opal, fcon_opal, ftredge, fu_opal, fzedge, dopactd, dopacr,
	2 dopact, ofebrack, opact, slr, slr_opal, slt,
	3 slt_opal, xc_opal, x, xheavy, xo_opal, x_opal, y, z, zmax, zmin,
	4 z_absent, z_opal
	INTEGER, SAVE, DIMENSION(nmet,max_iso) :: iso
	INTEGER, SAVE, DIMENSION(nmet) :: niso=0	
	INTEGER, SAVE, DIMENSION(max_iso) :: isoh=0, isohe=0	
	INTEGER, SAVE :: i_xheavy, nisoh=0, nisohe=0, nzmax
	INTEGER :: i, khighz

	LOGICAL, SAVE, DIMENSION(nmet) :: absent
	LOGICAL, SAVE :: init=.TRUE., pass=.TRUE.
	LOGICAL :: ok, sort, sort_r, sort_t

	CHARACTER (len=80) :: nom_chemin = "/data1/sdeheuve/local/src/cesam2k_v1.1.8_ESTA/SUN_STAR_DATA/"
	CHARACTER (len=100), SAVE, DIMENSION(16) :: chaine	
	CHARACTER (len=100) :: chain
	     
c-----------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 
c directory où se trouvent les fichiers de données OPAL	 

	 CALL set_opal_dir(nom_chemin)	 
	 
c identification des indices de C, N, O, 
c	 nom_xheavy est l'élément de complément du réseau nucléaire
c	 i_xheavy est son indice parmi les xmet
c	 il est supposé qu'un élément a au plus max_iso isotopes

	 DO i=1,nchim
	  IF(nom_elem(i)(1:2) == 'He')THEN
	   nisohe=nisohe+1 ; isohe(nisohe)=i	  	 
	  ELSEIF(nom_elem(i)(1:2) == ' H')THEN
	   nisoh=nisoh+1 ; isoh(nisoh)=i
	  ELSEIF(nom_elem(i)(1:2) == 'C1')THEN
	   niso(1)=niso(1)+1 ; iso(1,niso(1))=i
	  ELSEIF(nom_elem(i)(1:2) == 'N1')THEN
	   niso(2)=niso(2)+1 ; iso(2,niso(2))=i
	  ELSEIF(nom_elem(i)(1:2) == 'O1')THEN
	   niso(3)=niso(3)+1 ; iso(3,niso(3))=i
	  ELSEIF(nom_elem(i)(1:2) == 'Ne')THEN
	   niso(4)=niso(4)+1 ; iso(4,niso(4))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=4
	  ELSEIF(nom_elem(i)(1:2) == 'Na')THEN
	   niso(5)=niso(5)+1 ; iso(5,niso(5))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=5	   
	  ELSEIF(nom_elem(i)(1:2) == 'Mg')THEN
	   niso(6)=niso(6)+1 ; iso(6,niso(6))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=6	   
	  ELSEIF(nom_elem(i)(1:2) == 'Al')THEN
	   niso(7)=niso(7)+1 ; iso(7,niso(7))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=7	   
	  ELSEIF(nom_elem(i)(1:2) == 'Si')THEN
	   niso(8)=niso(8)+1 ; iso(8,niso(8))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=8	   
	  ELSEIF(nom_elem(i)(1:1) == 'P')THEN
	   niso(9)=niso(9)+1 ; iso(9,niso(9))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=9	   
	  ELSEIF(nom_elem(i)(1:1) == 'S')THEN
	   niso(10)=niso(10)+1 ; iso(10,niso(10))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=10	   
	  ELSEIF(nom_elem(i)(1:2) == 'Cl')THEN
	   niso(11)=niso(11)+1 ; iso(11,niso(11))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=11	   	   
	  ELSEIF(nom_elem(i)(1:2) == 'Ar')THEN
	   niso(12)=niso(12)+1 ; iso(12,niso(12))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=12	   
	  ELSEIF(nom_elem(i)(1:1) == 'K')THEN
	   niso(13)=niso(13)+1 ; iso(13,niso(13))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=13	   
	  ELSEIF(nom_elem(i)(1:2) == 'Ca')THEN
	   niso(14)=niso(14)+1 ; iso(14,niso(14))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=14	   
	  ELSEIF(nom_elem(i)(1:2) == 'Ti')THEN
	   niso(15)=niso(15)+1 ; iso(15,niso(15))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=15	   
	  ELSEIF(nom_elem(i)(1:2) == 'Cr')THEN
	   niso(16)=niso(16)+1 ; iso(16,niso(16))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=16	   
	  ELSEIF(nom_elem(i)(1:2) == 'Mn')THEN
	   niso(17)=niso(17)+1 ; iso(17,niso(17))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=17	   
	  ELSEIF(nom_elem(i)(1:2) == 'Fe')THEN
	   niso(18)=niso(18)+1 ; iso(18,niso(18))=i
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=18	   
	  ELSEIF(nom_elem(i)(1:2) == 'Ni')THEN
	   niso(19)=niso(19)+1 ; iso(19,niso(19))=i 
	   IF(nom_elem(i) == nom_xheavy)i_xheavy=19	   
	  ENDIF	  	 	  	
	 ENDDO
	 
c	 PRINT*,nom_xheavy ; PRINT*,i_xheavy ; PAUSE'nom_xheavy'
	 
c les éléments absents du réseau nucléaire, xheavy élément moyen est absent

	 DO i=1,nmet
	  absent(i)=niso(i) <= 0
	 ENDDO
	 absent(i_xheavy)=.TRUE.
c	 PRINT*,absent
	 
c le Z des absents, xmet commence à C=6, abon_m	à H  
	 
	  z_absent=0.d0
	  DO i=1,nmet
	   IF(absent(i))z_absent=z_absent+abon_m(i+5) 
	  ENDDO
c	  WRITE(*,2000)z_absent
	  
c les rapports xisz_absent=Xi/z_absent de la mixture
	 
	 DO i=1,nmet
	  IF(absent(i))xisz_absent(i)=abon_m(i+5)/z_absent
	 ENDDO
c	 WRITE(*,2000)xisz_absent ; PAUSE'absent'	 
	 
c détermination de ofebrack=[O/Fe]=LOG10(O/Fe_mixture)-LOG10(O/Fe_solaire)
c	en nombre, mo, mfe : masses atomiques de O et Fe
	 osfe_sol=8.87-7.5	!LOG10 O/Fe_solaire en nombre de GN93
	 IF(ife56 > 0)THEN 
	  ofebrack=LOG10(SUM(xh(iso(3,1):iso(3,niso(i))))*mo/
	1 SUM(xh(iso(18,1):iso(18,niso(i))))*mfe)-osfe_sol
	 ELSE
	  ofebrack=0.
	 ENDIF	
	 
	 IF(cno)THEN
	  SELECT CASE(langue)	 
	  CASE('english')	
	   WRITE(*,1009) ; WRITE(2,1009)
1009	   FORMAT(/,'Use of OPAL_X_CNO_FU')	  
	  CASE DEFAULT	  
	   WRITE(*,9) ; WRITE(2,9)
9	   FORMAT(/,'Utilisation de OPAL_X_CNO_FU')
	  END SELECT
	 ELSE	 
	  SELECT CASE(langue)	 
	  CASE('english')	
	   WRITE(*,1008) ; WRITE(2,1008)
1008	   FORMAT(/,'Use of the subroutine OPAL')	  
	  CASE DEFAULT	  
	   WRITE(*,8) ; WRITE(2,8)
8	   FORMAT(/,'Utilisation de la routine OPAL')
	  END SELECT
	 ENDIF	 

c rapports C/Z, O/Z pour les excès de C et O
	 IF(nom_abon == 'enhan_cha')THEN
	  khighz=13. ; c0=9.194E-02 ; o0=6.426E-01
	  chain=TRIM(nom_chemin)//'C95hz_withCNO.bin'	  	  	  	  
	 ELSEIF(nom_abon == 'meteorites')THEN
	  khighz=11 ; c0=1.615E-01 ; o0=5.067E-01
	  chain=TRIM(nom_chemin)//'GN93hz_withCNO.bin'
	 ELSEIF(nom_abon == 'enhan_w')THEN
	  khighz=14 ; c0=7.646E-02 ; o0=6.728E-01
	  chain=TRIM(nom_chemin)//'W95hz_withCNO.bin'
	 ELSEIF(nom_abon == 'enhan_al')THEN
	  khighz=12 ; c0=1.027E-01 ; o0=5.702E-01
	  chain=TRIM(nom_chemin)//'Alrd96_withCNO.bin'
	 ELSEIF(nom_abon == 'solaire_gs')THEN	
	  khighz=-11 ; c0=1.733E-01 ; o0=4.822E-01
	  chain=TRIM(nom_chemin)//'GS98hz_withCNO.bin'
	 ELSE	!(solaire_gn par défaut Grevesse-Noels 93)
	  khighz=11 ; c0=1.733E-01 ; o0=4.822E-01
	  chain=TRIM(nom_chemin)//'GN93hz_withCNO.bin'
	 ENDIF
	  
c lecture de la table en binaire si elle existe, sinon création
	  
	 INQUIRE(file=TRIM(chain),exist=ok)
	 IF(ok)THEN
	  SELECT CASE(langue)	  
	  CASE('english')	
	   WRITE(*,1005)TRIM(chain) ; WRITE(2,1005)TRIM(chain)
1005	   FORMAT('with the opacity binary data file : ',a)	  
	  CASE DEFAULT	  
	   WRITE(*,5)TRIM(chain) ; WRITE(2,5)TRIM(chain)
5	   FORMAT('avec le fichier binaire :',
	1  /,a,/)
	  END SELECT
	  CALL read_opal_dump(23,TRIM(chain))
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1006)TRIM(chain) ; WRITE(2,1006)TRIM(chain)
1006	   FORMAT('CESAM constructs the binary opacity file : ',a,
	1  /,'that will take a very long long time')	  
	  CASE DEFAULT	  
	   WRITE(*,6)TRIM(chain) ; WRITE(2,6)TRIM(chain)
6	   FORMAT('Pour les opacités CESAM crée le fichier binaire : ',
	1  a,/,'ce qui sera très très long')
	  END SELECT
	  CALL ask_z_limits(nzmax,zmin,zmax )	 
	  CALL readzexco(nzmax,0.0,0.05,0.1,khighz,iulow,ofebrack)
	  CALL dump_opal_opac(23,chain)
	  SELECT CASE(langue)	  
	  CASE('english')	
	   WRITE(*,1007)TRIM(chain) ; WRITE(2,1007)TRIM(chain)
1007	   FORMAT('End of the building of the binary file : ',a)	  
	  CASE DEFAULT	  
	   WRITE(*,7)TRIM(chain) ; WRITE(2,7)TRIM(chain)
7	   FORMAT('Fin de la création du fichier binaire : ',a)
	  END SELECT
	 ENDIF

c les limites	 
	 CALL ask_logt6_limits(slt_min,dslt_min,slt_max,dslt_max)
	 CALL ask_logr_limits(slr_min,dslr_min,slr_max,dslr_max)
c	 WRITE(*,2000)slt_min,dslt_min,slt_max,dslt_max	 	 
c	 WRITE(*,2000)slr_min,dslr_min,slr_max,dslr_max

	 slt_min=slt_min-dslt_min ; slt_max=slt_max+dslt_max
	 slr_min=slr_min-dslr_min ; slr_max=slr_max+dslr_max
	 
c	 WRITE(*,2000)slt_min,slt_max,slr_min,slr_max	 	 
c	 PAUSE'limites'	
	ENDIF
	
c calcul de l'opacité	
	
	x=SUM(xh(isoh(1):isoh(nisoh))) 		!fraction totale de H
	y=SUM(xh(isohe(1):isohe(nisohe))) 	!fraction totale de He
	z=1.d0-x-y
	xmet=0.
	DO i=1,nmet
	 IF(niso(i) > 0)xmet(i)=SUM(xh(iso(i,1):iso(i,niso(i))))
	ENDDO
c	WRITE(*,2000)x,y,z ; WRITE(*,2000)xmet ; PAUSE'xmet'

c abondance de xheavy : xmet(i_xheavy)

	xheavy=xmet(i_xheavy)
	DO i=1,nmet
	 IF(absent(i))xmet(i)=xisz_absent(i)*xheavy
	ENDDO

c au premier passage on écrit xh, xmet, c0, O0, z_table 
	
	IF(init)THEN
	 init=.FALSE.
	 DO i=1,nchim
	  CHAINE(i)=TRIM(ADJUSTL(nom_elem(i)))//'='
	  IF(MOD(i,5) /= 1)CHAINE(i)=', '//CHAINE(i)
	 ENDDO	 	 	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1014)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))
	  WRITE(2,1014)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))	  
1014	  FORMAT('at the first call',/,5(a,es10.3)) 
	  IF(nchim > 5)WRITE(*,15)(TRIM(chaine(i)),xh(i),
	1 i=6,MIN(10,nchim))
	  IF(nchim > 10)WRITE(*,15)(TRIM(chaine(i)),xh(i),
	1 i=11,MIN(15,nchim))	 	
	  WRITE(*,1016)xmet,c0,O0,x,y,z
	  WRITE(2,1016)xmet,c0,O0,x,y,z	  
1016	  FORMAT('C,N,O,Ne,Na,Mg,Al,Si',/,
	1 8es10.3,/,'P,S,Cl,Ar,K,Ca,Ti,Cr',/,8es10.3,/,'Mn,Fe,Ni',/,
	2 3es10.3,/,'c0=',es10.3,', O0=',es10.3,', X=',es10.3,
	3 ', Y=',es10.3', Z=',es10.3)
	 CASE DEFAULT
	  WRITE(*,14)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))
	  WRITE(2,14)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))	  
14	  FORMAT('au premier appel',/,5(a,es10.3)) 
	  IF(nchim > 5)THEN
	   WRITE(*,15)(TRIM(chaine(i)),xh(i),i=6,MIN(10,nchim))
	   WRITE(2,15)(TRIM(chaine(i)),xh(i),i=6,MIN(10,nchim))
	  ENDIF 
	  IF(nchim > 10)THEN
	   WRITE(*,15)(TRIM(chaine(i)),xh(i),i=11,MIN(15,nchim))
	   WRITE(2,15)(TRIM(chaine(i)),xh(i),i=11,MIN(15,nchim))
	  ENDIF	  
15	  FORMAT(5(a,es10.3))	
	  WRITE(*,16)xmet,c0,O0,x,y,z
	  WRITE(2,16)xmet,c0,O0,x,y,z	  
16	  FORMAT('C,N,O,Ne,Na,Mg,Al,Si',/,8es10.3,/,
	1 'P,S,Cl,Ar,K,Ca,Ti,Cr',/,8es10.3,/,'Mn,Fe,Ni',/,3es10.3,/,
	2 'c0=',es10.3,', O0=',es10.3,', X=',es10.3,
	3 ', Y=',es10.3', Z=',es10.3)
	 END SELECT		
	ENDIF

	slt=LOG10(t)-6. ; slr=LOG10(ro)-3.*slt
	
	IF(pass)THEN	
	 sort_t=slt < slt_min .OR. slt > slt_max
	 sort_r=slr < slr_min .OR. slr > slr_max
	 sort=sort_t .OR. sort_r		
	 IF(sort)THEN
	  pass=.FALSE.	 !on ne passera qu'une fois
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1010)slt,slt_min,slt_max,slr,slr_min,slr_max
	   WRITE(2,1010)slt,slt_min,slt_max,slr,slr_min,slr_max
1010	   FORMAT(/,'Out of range in opal2, at least one time : ',/,
	1  'or slt = ',es10.3,' out of [',es10.3,',',es10.3,']',/,
	2  'and/or slr = ',es10.3,' out of [',es10.3,',',es10.3,']',/,
	3  'one uses the values at the limits',/)	  
	  CASE DEFAULT	  
	   WRITE(*,10)slt,slt_min,slt_max,slr,slr_min,slr_max
	   WRITE(2,10)slt,slt_min,slt_max,slr,slr_min,slr_max
10	   FORMAT(/,'Au moins une sortie de table dans opal2 : ',/,
	1  'ou slt = ',es10.3,' en dehors de [',es10.3,',',es10.3,']',/,
	2  'et/ou slr = ',es10.3,' en dehors de [',es10.3,',',es10.3,']',/,
	3  'on utilise les valeurs limites aux bord des tables',/)
	  END SELECT
	 ENDIF
	ENDIF

c pour éviter de sortir	
	slt=MIN(slt_max,MAX(slt,slt_min))
	slr=MIN(slr_max,MAX(slr,slr_min))

c calcul de l'opacité	
	IF(cno)THEN	
	 CALL opal_x_cno_fu(x,slt,slr,xmet,nmet,1.)	!interpolation
	 CALL ask_last_opac(opact,dopact,dopacr,dopactd,
	1 fedge,ftredge,fzedge) 	!opacité, dérivées, tests du calcul

c vérifications pour tests
	 IF(.FALSE.)THEN !IF(.TRUE.)THEN	
	  CALL ask_last_xcnou(z_opal,x_opal,xc_opal,xo_opal,slt_opal,
	1 slr_opal,fcn_opal,fcon_opal,fcnone_opal,fu_opal)
	  WRITE(*,1)z,nmet ;  WRITE(*,1)z,nmet
	  WRITE(*,2)xmet ; WRITE(2,2)xmet
	  WRITE(*,4) z_opal, x_opal, xc_opal, xo_opal, fcn_opal, fcon_opal,
	1 fcnone_opal, fu_opal
	  WRITE(2,4) z_opal, x_opal, xc_opal, xo_opal, fcn_opal, fcon_opal,
	1 fcnone_opal, fu_opal
	  WRITE(*,11)xh(1:MIN(8,nchim)) ; WRITE(*,2000)xh(9:MIN(16,nchim))		
	  WRITE(2,11)xh(1:MIN(8,nchim)) ; WRITE(*,2000)xh(9:MIN(16,nchim))
	  WRITE(*,12)ftredge,fzedge,slt,slr,t,ro
	  WRITE(2,12)ftredge,fzedge,slt,slr,t,ro
	  WRITE(*,13)opact,dopact,dopacr,dopactd
	  WRITE(2,13)opact,dopact,dopacr,dopactd	  	
	 ENDIF
	 
c sortie de table 
	 IF(fedge == 0.d0)THEN
	  CALL ask_last_xcnou(z_opal,x_opal,xc_opal,xo_opal,slt_opal,
	1 slr_opal,fcn_opal,fcon_opal,fcnone_opal,fu_opal)
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1011) ; WRITE(2,1011)	
1011	   FORMAT('STOP, out of opacity table limits')	  
	  CASE DEFAULT	
	   WRITE(*,11) ; WRITE(2,11)	
11	   FORMAT('ARRET, sortie de table d''opacité')
	  END SELECT	
	  WRITE(*,1)z,nmet ; WRITE(2,1)z,nmet
1	  FORMAT('Z=',es10.3,', nmet=',i4)
	  WRITE(*,2)xmet ; WRITE(2,2)xmet
2	  FORMAT('C,N,O,Ne,Na,Mg,Al,Si',/,8es10.3,/,
	1 'P,S,Cl,Ar,K,Ca,Ti,Cr',/,8es10.3,/,'Mn,Fe,Ni',/,3es10.3)
	  WRITE(*,4)z_opal, x_opal, xc_opal, xo_opal, fcn_opal, fcon_opal,
	1 fcnone_opal, fu_opal
	  WRITE(2,4)z_opal, x_opal, xc_opal, xo_opal, fcn_opal, fcon_opal,
	1 fcnone_opal, fu_opal	
4	  FORMAT('Z_opal=',es10.3,', X_opal=',es10.3,', exC,O=',2es10.3 ,/,
	1 'fCN,fCON,fCNONe,fu=',4es10.3)
	  WRITE(*,15)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))
	  WRITE(2,15)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))	  
	  IF(nchim > 5)THEN
	   WRITE(*,15)(TRIM(chaine(i)),xh(i),i=6,MIN(10,nchim))
	   WRITE(2,15)(TRIM(chaine(i)),xh(i),i=6,MIN(10,nchim))
	  ENDIF
	  IF(nchim > 10)THEN
	   WRITE(*,15)(TRIM(chaine(i)),xh(i),i=11,MIN(15,nchim))
	   WRITE(2,15)(TRIM(chaine(i)),xh(i),i=11,MIN(15,nchim))
	  ENDIF
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1012)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro
	   WRITE(2,1012)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro	
1012	   FORMAT('ftredge=',es10.3,', fzedge=',es10.3,/,'slt=', es10.3,
	1  ' out of [',es10.3,',',es10.3,']',/,'and/or slr=',es10.3,
	2  ' out of [',es10.3,',',es10.3,']',/,'T=',es10.3,', ro=',es10.3)
	  CASE DEFAULT		  
	   WRITE(*,12)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro
	   WRITE(2,12)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro	
12	   FORMAT('ftredge=',es10.3,', fzedge=',es10.3,/,'slt=', es10.3,
	1  ' en dehors de [',es10.3,',',es10.3,']',/,'et/ou slr=',es10.3,
	2  ' en dehors de [',es10.3,',',es10.3,']',/,'T=',es10.3,', ro=',
	3  es10.3)
	  END SELECT
	  WRITE(*,13)opact,dopact,dopacr,dopactd
	  WRITE(2,13)opact,dopact,dopacr,dopactd	  	  
13	  FORMAT('opact=',es10.3,', dopact=',es10.3,', dopacr=',es10.3,
	1 ', dopactd=',es10.3)
	  CALL sortie !STOP	  
	 ENDIF
	  
	 kappa=10.d0**opact ; dkapdt=kappa*dopact/t
	 dkapdr=kappa*dopacr/ro ; dkapdx=0.d0

c d kappa / dx	 
	 IF(.FALSE.)THEN 
c	 IF(.TRUE.)THEN		
	  IF(x > ab_min(1))THEN
	   x=x*unpdx
	   CALL opal_x_cno_fu(x,slt,slr,xmet,nmet,1.)
	   CALL ask_last_opac(opact,dopact,dopacr,dopactd,
	1  fedge,ftredge,fzedge) 
	   dkapdx=(10.d0**opact-kappa)/x/dx
	  ENDIF
	 ENDIF
	 
c utilisation simplifiée, appel à opal  	
	ELSE
	
	 z_table=MAX(0.0,MIN(0.1,z))	
	 ex_c=xmet(1)-c0*z_table ; ex_o=xmet(3)-o0*z_table
	 CALL opal(z_table,x,ex_c,ex_o,slt,slr)
	 CALL ask_last_opac(opact,dopact,dopacr,dopactd,
	1 fedge,ftredge,fzedge)		
	 kappa=10.d0**opact ; dkapdt=kappa*dopact/t
	 dkapdr=kappa*dopacr/ro ; dkapdx=0.d0

c sortie de table 
	 IF(fedge == 0.d0)THEN
	  CALL ask_last_xcnou(z_opal,x_opal,xc_opal,xo_opal,slt_opal,
	1 slr_opal,fcn_opal,fcon_opal,fcnone_opal,fu_opal)
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1011) ; WRITE(2,1011)	
	  CASE DEFAULT	
	   WRITE(*,11) ; WRITE(2,11)	
	  END SELECT	
	  WRITE(*,1)z,nmet ; WRITE(2,1)z,nmet
	  WRITE(*,2)xmet ; WRITE(2,2)xmet
	  WRITE(*,4)z_opal, x_opal, xc_opal, xo_opal, fcn_opal, fcon_opal,
	1 fcnone_opal, fu_opal
	  WRITE(2,4)z_opal, x_opal, xc_opal, xo_opal, fcn_opal, fcon_opal,
	1 fcnone_opal, fu_opal	
	  WRITE(*,15)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))
	  WRITE(2,15)(TRIM(chaine(i)),xh(i),i=1,MIN(5,nchim))	  
	  IF(nchim > 5)THEN
	   WRITE(*,15)(TRIM(chaine(i)),xh(i),i=6,MIN(10,nchim))
	   WRITE(2,15)(TRIM(chaine(i)),xh(i),i=6,MIN(10,nchim))
	  ENDIF
	  IF(nchim > 10)THEN
	   WRITE(*,15)(TRIM(chaine(i)),xh(i),i=11,MIN(15,nchim))
	   WRITE(2,15)(TRIM(chaine(i)),xh(i),i=11,MIN(15,nchim))
	  ENDIF
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1012)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro
	   WRITE(2,1012)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro	
	  CASE DEFAULT		  
	   WRITE(*,12)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro
	   WRITE(2,12)ftredge,fzedge,slt,slt_min,slt_max,slr,slr_min,
	1  slr_max,t,ro	
	  END SELECT
	  WRITE(*,13)opact,dopact,dopacr,dopactd
	  WRITE(2,13)opact,dopact,dopacr,dopactd	  	  
	  CALL sortie	!STOP	  
	 ENDIF
	 
c d kappa / dx	 
	 IF(.FALSE.)THEN
c	 IF(.TRUE.)THEN
	  IF(x > ab_min(1))THEN
	   x=x*unpdx
	   CALL opal(z_table,x, ex_c, ex_o,slt,slr) 
	   CALL ask_last_opac(opact,dopact,dopacr,dopactd,
	1  fedge,ftredge,fzedge) 
	   dkapdx=(10.d0**opact-kappa)/x/dx
	  ENDIF
	 ENDIF
	ENDIF
	
	RETURN
	
	END SUBROUTINE opa_opal2
