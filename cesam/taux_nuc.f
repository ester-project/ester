
c******************************************************************

	SUBROUTINE taux_nuc(t,total,rt,zit,izzt,qt,nom_react,nucleot,
	1 at,nom_elemt)

c	routine private du module mod_nuc
	
c	calcul des taux des r�actions thermonucl�aires

c	Auteur: P.Morel, D�partement J.D. Cassini, O.C.A.
c	CESAM2k

c r�f�rences :

c	Caughlan et Al. Atomic Data and Nuclear Data Tables, 40, 283, 1988
c	Compilation Adelberger
c	Compilation NACRE
c	Landr� et al. 1990 A&A 240, 85
c	Nomoto et al. astro-ph/9706024

c	B. Pichon DARC Observatoire de Paris-Meudon : calcul detaille des
c	�nergies lib�r�es par les r�actions thermonucl�aires dans les
c	chaines PP et CNO, 1989

c entr�es d�finies dans tabul_nuc:
c	t : temp�rature en K
c	total=1, 2, 3 : compilations de Caughlan & Fowler,
c	 d'Adelberger, de NACRE
c sorties d�finies dans tabul_nuc :
c	rt : ln(taux des r�actions)
c	zit : charges des noyaux
c	izzt :charges des noyaux Z1, Z2 des noyaux de la r�action
c	qt : �nergies
c	nom_react : noms des r�actions
c	nucleot : masses des isotopes
c	at : masses r�duites
c	nom_elemt : noms des �l�ments

c	dans reac_nuc on calcule toutes les r�actions
c	lors de la tabulation dans tabul_nuc on ne retiendra que
c	celles qui sont n�cessaires
c	le fait qu'il n'y a pas d'ordre logique dans la suite
c	des r�actions/�l�ments, est une cons�quence
c	des additions successives

c pour ajouter 1 isotope: 
c	1) dans le module mod_nuc augmenter le PARAMETER niso_tot d'une
c	 unite, Ex: niso_tot=28
c	2) dans taux_nuc:
c	 - ajouter l'exc�s de masse de l'isotope et son
c	symbole dans les d�clarations, Ex: Fe56=-60.6054d0
c	 - ajouter la masse de l'isotope, (�ventuellement, ce qui est
c	 conseill� l'ajouter dans les constantes physiques ctes_85,
c	 ctes_95 et dans mod_donnees),
c	 Ex: nucleot(27)=afe56 (afe56=55.847d0 et afe56, af18...)
c	 - ajouter la charge de l'isotope, Ex: zit(28)=11
c	-ajouter le nom de l'isotope, Ex: 'B11 ' (ces chaines
c	de 4 caracteres sont � cadrer � gauche sauf pour les
c	noms avec 2 caract�res qui sont centr�s) 

c Pour ajouter une r�action,
c	1)dans taux_nuc:
c	 - compl�ter la liste de commentaires ci-dessous,
c	  Ex: r�action 46 : B11(p,g)12C
c	 - initialiser nom_reac, izz(i,1)..
c	  Ex:	 nom_react(2)='H2(p,g)He3' ; nuc=H2+p-He3
c	  qt(2)=nuc ; izzt(2,1)=1 ; izzt(2,2)=1
c	- entrer la r�action rt(i) en log
c	  penser a inclure dans rt(i) les 1! ou 2! ou 3! du
c	  d�nominateur, Ex: r�action 8 : C12(p,g)N13(e+ nu)C13
c	  z0=6, z1=1........
c	 on recalcule au besoin les coefficients � partir
c	 des S(0), S'(0),
c	 S"(0) en utilisant Lang Astrophysical formulae eq.4-49

c	 le fait que les calculs soient effectu�s en ln est historique
c	 pour d'une part, �viter les erreurs de troncature et surtout pour
c	 �viter les d�passements de capacit� avec la limitation � 10.d38

c	2) dans tabul_nuc:
c	 entrer les parametres du calcul des taux,
c	 nombre de r�actions et correspondance des indices
c	 de la routine de calcul des r�actions, Ex: 
c	 CASE('ppcno12')
c	 WRITE(2,6) ; WRITE(*,6)....

c	3)dans le module mod_nuc: 
c	 - introduire par un INCLUDE le nom de la routine de calcul des
c	 taux, Ex: INCLUDE 'ppcno9.f'

c	nom_elem : symboles des �l�ments chimiques utilis�s
c	H1, H2, He3, He4, Li7, Be7, C12, C13, N13, N14, N15, O16, O17, O18,
c	F19, Ne20, Ne21, Ne22, Ne23, Na23, Mg23, Mg24, Mg25, Mg26, n, Be9,
c	Li6, Fe56, B11

c	r�actions disponibles:

c	r�action 1 : H1(p,e+ nu)H2		PP
c	r�action 2 : H2(p,g)He3
c	r�action 3 : He3(He3,2p)He4
c	r�action 4 : He4(He3,g)Be7
c	r�action 5 : Li7(p,a)He4
c	r�action 6 : Be7(e-,nu g)Li7
c	r�action 7 : Be7(p,g)B8(e+ nu)Be8(a)He4

c	r�action 8 : C12(p,g)N13(e+ nu)C13	CNO
c	r�action 9 : C13(p,g)N14
c	r�action 10 : N14(p,g)O15(e+ nu)N15
c	r�action 11 : N15(p,g)O16
c	r�action 12 : N15(p,a)C12
c	r�action 13 : O16(a,g)F17(e+ nu)O17
c	r�action 14 : O17(p,a)N14

c	r�action 15 : He4(2a,g)C12		3 alpha
c	r�action 16 : C12(a,g)O16
c	r�action 17 : O16(a,g)Ne20

c	r�action 18 : C12(C12,g)Mg24		carbone
c	r�action 19 : C12(C12,n)Mg23
c	r�action 20 : C12(C12,p)Na23
c	r�action 21 : C12(C12,a)Ne20
c	r�action 22 : C12(C12,2a)O16
c	r�action 23 : N13(p,g)O14(e+nu)N14
c	r�action 24 : O17(p,g)F18(e+nu)O18
c	r�action 25 : O18(p,g)F19
c	r�action 26 : F19(p,a)O16
c	r�action 27 : F19(p,g)Ne20
c	r�action 28 : O18(p,a)N15
c	r�action 29 : Ne20(a,g)Mg24

c	r�action 30 : C13(a,n)O16		oxyg�ne
c	r�action 31 : O17(a,n)Ne20
c	r�action 32 : N14(a,g)F18(e+nu)O18
c	r�action 33 : O18(a,g)Ne22
c	r�action 34 : Ne22(a,n)Mg25
c	r�action 35 : Ne22(a,g)Mg26
c	r�action 36 : C12(a,n)O15(e+nu)N15
c	r�action 37 : Ne21(a,n)Mg24
c	r�action 38 : He4(an,g)Be9

c	r�action 39 : Be9(p,d)2He4		divers
c	r�action 40 : Be9(a,n)C12
c	r�action 41 : Li6(p,He3)He4
c	r�action 42 : Li6(p,g)Be7
c	r�action 43 : Be9(p,a)Li6
c	r�action 44 : B11(p,a)2He4
c	r�action 45 : B11(p,g)C12

c---------------------------------------------------------------------

	USE mod_donnees, only : an, ah, ah2, ahe3, ahe4, ali6, ali7, 
	1 abe7, abe9, ab11, ac12, ac13, an13, an14, an15, ao16, ao17, ao18,
	2 af18, af19, ane20, ane21, ane22, ana23, amg23, amg24, amg25,
	3 amg26, afe56, langue, nom_nuc_cpl
	USE mod_kind	
	USE mod_numerique, ONLY : polyder

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: t	
	INTEGER, INTENT(in) :: total		
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: at, qt, rt	
	REAL (kind=dp), INTENT(out), DIMENSION(-1:) :: nucleot, zit
	INTEGER, INTENT(out), DIMENSION(:,:) :: izzt
	CHARACTER (len=*), INTENT(out), DIMENSION(-1:) :: nom_elemt
	CHARACTER (len=*), INTENT(out), DIMENSION(:) :: nom_react	


	REAL (kind=dp) :: ai(0:8), po(0:8), z1z2a, c1, c2, s0, sq, ss

	REAL (kind=dp) :: t9, lnt9a, t923, t953, lnt9, lnt912, t912,
	1 lnt932, lnt913, t913, lnt923, lnt934, lnt953, v0, v1,
	2 v2, v3, v4, v5,
	3 he4abe8, be8agc12, t9a13, t9a, t932,
	4 lnfpt9a, lngt9, lnft9a, t6, p, a, n, e, H1, H2, He3, He4,
	5 Li7, Be7, B11, C12, C13, N13, N14, N15, O16,
	6 O17, O18, F18, F19, Ne20, Ne21, Ne22, Na23, Mg23, Mg24,
	7 Mg25, Mg26, Be9, Li6, dw, enu, nuc, qd, neu
c	  Be8, B8, Fe56, 	
	
	REAL (kind=dp), SAVE :: ln38, zero21c, zero21l	
	
	INTEGER :: i
	
	LOGICAL, SAVE :: init=.TRUE.

c----------------------------------------------------------------------
	
	dw(qd)=1.d0/(qd/0.51099906d0+1.d0)
	neu(qd)=(1.d0+dw(qd))*(1.d0-dw(qd)*(0.25d0+dw(qd)/9.d0))/2.d0
	
2000	FORMAT(8es10.3)
2001	FORMAT(8es11.4)

	IF(init)THEN
	 init=.FALSE.
	 
	 zero21c=0.1d0		!0 to 1 pour Caughlan
	 zero21l=0.5d0		!0 to 1 pour Landre	 

	 ln38=LOG(1.d-100)	!LOG(0)
	
c	 exc�s de masse Clayton, p 289
	
	 p=7.28899d0
	 a=2.42475d0
	 n=8.07144d0
	 e=0.5110034d0
	 
	 H1=p
	 H2=13.13591d0
	 He3=14.93134d0
	 He4=a
	 Li6=14.08840d0
	 Li7=14.90730d0
	 Be7=15.76890d0
c	 Be8=4.94420d0
c	 B8=22.92310d0
	 B11=8.66768
	 C12=0.d0
	 C13=3.12460d0
	 N13=5.34520d0
	 N14=2.86373d0
	 N15=0.10040d0
	 O16=-4.73655d0
	 O17=-0.80770d0
	 O18=-0.78243d0
	 F18=0.87240d0
	 F19=-1.48600d0	
	 Ne20=-7.04150d0
	 Ne21=-5.7299d0
	 Ne22=-8.02490d0
	 Na23=-9.52830d0	
	 Mg23=-5.47240d0
	 Mg24=-13.93330d0
	 Mg25=-13.19070d0
	 Mg26=-16.21420d0
	 Be9=11.35050d0
c	 Fe56=-60.6054d0
	 
c	 initialisation des nucleot : masses atomiques

	 nucleot(-1)=0.d0	!Ex fictif
	 nucleot(0)=an
	 nucleot(1)=ah
	 nucleot(2)=ah2
	 nucleot(3)=ahe3
	 nucleot(4)=ahe4
	 nucleot(5)=ali7
	 nucleot(6)=abe7
	 nucleot(7)=ac12
	 nucleot(8)=ac13
	 nucleot(9)=an13	 
	 nucleot(10)=an14
	 nucleot(11)=an15
	 nucleot(12)=ao16
	 nucleot(13)=ao17
         nucleot(14)=ao18
         nucleot(15)=af19
         nucleot(16)=ane20
         nucleot(17)=ane22
         nucleot(18)=amg24         
         nucleot(19)=amg25
         nucleot(20)=amg26
         nucleot(21)=af18
         nucleot(22)=ane21
         nucleot(23)=ana23
         nucleot(24)=amg23
         nucleot(25)=abe9
         nucleot(26)=ali6
         nucleot(27)=afe56
	 nucleot(28)=ab11     
                     
c	 initialisation des charges

	 zit(-1)=0	!Ex fictif
	 zit(0)=0	!n
	 zit(1)=1	!H1
	 zit(2)=1	!H2
	 zit(3)=2	!He3
	 zit(4)=2	!He4
	 zit(5)=3	!Li7
	 zit(6)=4	!Be7
	 zit(7)=6	!C12
	 zit(8)=6	!C13
	 zit(9)=7	!N13
	 zit(10)=7	!N14
	 zit(11)=7	!N15
	 zit(12)=8	!O16
	 zit(13)=8	!O17
         zit(14)=8	!O18
         zit(15)=9	!F19
         zit(16)=10	!Ne20
         zit(17)=10	!Ne22
         zit(18)=12	!Mg24
         zit(19)=12	!Mg25
         zit(20)=12	!Mg26
         zit(21)=9	!F18
         zit(22)=10	!Ne21
         zit(23)=11	!Na23
         zit(24)=12	!Mg23
         zit(25)=4	!Be9
         zit(26)=3	!Li6
         zit(27)=26	!Fe56
	 zit(28)=5	!B11        
        	 
c	 initialisation des noms des �l�ments

	 nom_elemt(-1)='Ex '
	 nom_elemt(0)=' n  '	 
	 nom_elemt(1)=' H1 '
	 nom_elemt(2)=' H2 '
	 nom_elemt(3)='He3 '
	 nom_elemt(4)='He4 '
	 nom_elemt(5)='Li7 '
	 nom_elemt(6)='Be7 '
	 nom_elemt(7)='C12 '
	 nom_elemt(8)='C13 '
	 nom_elemt(9)='N13 '	 
	 nom_elemt(10)='N14 '
	 nom_elemt(11)='N15 '
	 nom_elemt(12)='O16 '
	 nom_elemt(13)='O17 '
	 nom_elemt(14)='O18 '	 
	 nom_elemt(15)='F19 '	 
	 nom_elemt(16)='Ne20'
	 nom_elemt(17)='Ne22'
	 nom_elemt(18)='Mg24'
	 nom_elemt(19)='Mg25'	 
	 nom_elemt(20)='Mg26'
	 nom_elemt(21)='F18 '	 
	 nom_elemt(22)='Ne21'	 
	 nom_elemt(23)='Na23'
	 nom_elemt(24)='Mg23'
	 nom_elemt(25)='Be9 '
	 nom_elemt(26)='Li6 '
	 nom_elemt(27)='Fe56'
	 nom_elemt(28)='B11 '	 	 	 	 
	 	 
c initialisation des nom_reac, izz(i,1)=Z1, izz(i,2)=Z2, m et qt

c r�action 1 : H1(p,e+ nu)H2	 
	 nom_react(1)='H1(p,e+ nu)H2'	!pichon q=1.172
	 nuc=2.d0*H1-H2
	 enu=0.263d0			!cours Pichon p.287
	 qt(1)=nuc-enu
         izzt(1,1)=1
         izzt(1,2)=1	 
c        at(1)=1.d0*1.d0/(1.d0+1.d0)	!A : masse atomique r�duite
         at(1)=nucleot(1)*nucleot(1)/(nucleot(1)+nucleot(1))	!m r�duite

c r�action 2 : H2(p,g)He3         
	 nom_react(2)='H2(p,g)He3' ; nuc=H2+H1-He3 ; qt(2)=nuc
         izzt(2,1)=1 ; izzt(2,2)=1
c        at(2)=2.d0*1.d0/(2.d0+1.d0)	!A : masse atomique r�duite
         at(2)=nucleot(2)*nucleot(1)/(nucleot(2)+nucleot(1))	!m r�duite         

c r�action 3 : He3(He3,2p)He4         
	 nom_react(3)='He3(He3,2p)He4'
	 nuc=2.d0*He3-(2.d0*H1+He4)
	 qt(3)=nuc
         izzt(3,1)=2
         izzt(3,2)=2
c        at(3)=3.d0*3.d0/(3.d0+3.d0)	!A : masse atomique r�duite
         at(3)=nucleot(3)*nucleot(3)/(nucleot(3)+nucleot(3))	!m r�duite

c r�action 4 : He4(He3,g)Be7         
	 nom_react(4)='He3(a,g)Be7'
	 nuc=He3+a-Be7
	 qt(4)=nuc
         izzt(4,1)=2
         izzt(4,2)=2         
c        at(4)=3.d0*4.d0/(3.d0+4.d0)	!A : masse atomique r�duite
         at(4)=nucleot(3)*nucleot(4)/(nucleot(3)+nucleot(4))	!m r�duite

c r�action 5 : Li7(p,a)He4                 	 
	 nom_react(5)='Li7(p,a)He4'
	 nuc=Li7+H1-(a+He4)
	 qt(5)=nuc
         izzt(5,1)=3
         izzt(5,2)=1
c        at(5)=7.d0*1.d0/(7.d0+1.d0)	!A : masse atomique r�duite
         at(5)=nucleot(5)*nucleot(1)/(nucleot(5)+nucleot(1))	!m r�duite

c r�action 6 : Be7(e-,nu g)Li7                  	 
	 nom_react(6)='Be7(e-,nu g)Li7'
	 qt(6)=0.0497d0	!cf. B. Pichon le nu emporte pratiquement tout: Be7-Li7
         izzt(6,1)=0	!factice
         izzt(6,2)=0	!avec z1=z2(6)=0 l'effet d'ecran sera 1
c        at(6)=7.d0*1.d0/(7.d0+1.d0)	!A : masse atomique r�duite fictif
         at(6)=nucleot(6)*nucleot(1)/(nucleot(6)+nucleot(1))	!m r�duite         

c r�action 7 : Be7(p,g)B8(e+ nu)Be8(a)He4                  	 
	 nom_react(7)='Be7(p,a)He4'	!Be7(p,g)B8(e+ nu)Be8(a)He4
c	 qt(7)=0.137+11.261+0.092	!cf. B. Pichon
	 nuc=Be7+H1-(a+He4)
	 enu=6.71d0		!Bahcall, cf cours Pichon p.287, 305
	 qt(7)=nuc-enu
         izzt(7,1)=4
         izzt(7,2)=1
c        at(7)=7.d0*1.d0/(7.d0+1.d0)	!A : masse atomique r�duite
         at(7)=nucleot(6)*nucleot(1)/(nucleot(6)+nucleot(1))	!m r�duite

c r�action 8 : C12(p,g)N13(e+ nu)C13                           	 
	 nom_react(8)='C12(p,g)C13'	!C12(p,g)N13(e+ nu)C13
c	 qt(8)=1.9435d0+1.508d0	!cf. B. Pichon
	 nuc=C12+H1-C13
	 enu=0.713d0		!cours Pichon p. 287
	 qt(8)=nuc-enu
         izzt(8,1)=6
         izzt(8,2)=1
c        at(8)=12.d0*1.d0/(12.d0+1.d0)	!A : masse atomique r�duite
         at(8)=nucleot(7)*nucleot(1)/(nucleot(7)+nucleot(1))	!m r�duite

c r�action 9 : C13(p,g)N14                           	 
	 nom_react(9)='C13(p,g)N14'
c	 qt(9)=7.55063d0		!cf. B. Pichon
	 nuc=C13+H1-N14
	 qt(9)=nuc
         izzt(9,1)=6
         izzt(9,2)=1
c        at(9)=13.d0*1.d0/(13.d0+1.d0)	!A : masse atomique r�duite
         at(9)=nucleot(8)*nucleot(1)/(nucleot(8)+nucleot(1))	!m r�duite

c r�action 10 : N14(p,g)O15(e+ nu)N15         	
	 nom_react(10)='N14(p,g)N15'	!N14(p,g)O15(e+ nu)N15
c	 qt(10)=7.2971d0+1.757d0	!cf. B. Pichon
	 nuc=N14+H1-N15
	 enu=0.997d0		!cours Pichon p.287
	 qt(10)=nuc-enu
         izzt(10,1)=7
         izzt(10,2)=1
c        at(10)=14.d0*1.d0/(14.d0+1.d0)	!A : masse atomique r�duite
         at(10)=nucleot(10)*nucleot(1)/(nucleot(10)+nucleot(1))	!m r�duite

c r�action 11 : N15(p,g)O16         	
	 nom_react(11)='N15(p,g)O16'
	 nuc=N15+H1-O16
	 qt(11)=nuc
         izzt(11,1)=7
         izzt(11,2)=1
         at(11)=15.d0*1.d0/(15.d0+1.d0)	!A : masse atomique r�duite
         at(11)=nucleot(11)*nucleot(1)/(nucleot(11)+nucleot(1))	!m r�duite

c r�action 12 : N15(p,a)C12                  	 	 
	 nom_react(12)='N15(p,a)C12'
c	 nuc=4.96561d0		!cf. B. Pichon
	 nuc=N15+H1-(a+C12)
	 qt(12)=nuc
         izzt(12,1)=7
         izzt(12,2)=1
c        at(12)=15.d0*1.d0/(15.d0+1.d0)	!A : masse atomique r�duite
         at(12)=nucleot(11)*nucleot(1)/(nucleot(11)+nucleot(1))	!m r�duite

c r�action 13 : O16(a,g)F17(e+ nu)O17                           	 	 
	 nom_react(13)='O16(p,g)O17'	!O16(p,g)F17(e+ nu)O17
c	 qt(13)=0.60035d0+1.761d0	!cf. B. Pichon
	 nuc=O16+H1-O17
	 enu=1.001d0		!cours Pichon p.278
	 qt(13)=nuc-enu
         izzt(13,1)=8
         izzt(13,2)=1
c        at(13)=16.d0*1.d0/(16.d0+1.d0)	!A : masse atomique r�duite
         at(13)=nucleot(12)*nucleot(1)/(nucleot(12)+nucleot(1))	!m r�duite

c r�action 14 : O17(p,a)N14         	
	 nom_react(14)='O17(p,a)N14'
c	 qt(14)=1.1907d0		!cf. B. Pichon
	 nuc=O17+H1-(a+N14)
	 qt(14)=nuc
         izzt(14,1)=8
         izzt(14,2)=1
c        at(14)=17.d0*1.d0/(17.d0+1.d0)	!A : masse atomique r�duite
         at(14)=nucleot(13)*nucleot(1)/(nucleot(13)+nucleot(1))	!m r�duite

c r�action 15 : He4(2a,g)C12         	
	 nom_react(15)='He4(2a,g)C12'	!3 alpha q=7.275
	 nuc=3.*a-C12
	 qt(15)=nuc
         izzt(15,1)=2
         izzt(15,2)=2
c        at(15)=4.d0*4.d0/(4.d0+4.d0)	!A : masse atomique r�duite ?????
         at(15)=nucleot(4)*nucleot(4)/(nucleot(4)+nucleot(4))	!m r�duite

c r�action 16 : C12(a,g)O16         	
	 nom_react(16)='C12(a,g)O16'		!q=7.162
	 nuc=C12+a-O16
	 qt(16)=nuc
         izzt(16,1)=6
         izzt(16,2)=2
c        at(16)=12.d0*4.d0/(12.d0+4.d0)	!A : masse atomique r�duite
         at(16)=nucleot(7)*nucleot(4)/(nucleot(7)+nucleot(4))	!m r�duite

c r�action 17 : O16(a,g)Ne20         	
	 nom_react(17)='O16(a,g)Ne20'	!q=4.734
	 nuc=O16+a-Ne20
	 qt(17)=nuc
         izzt(17,1)=8
         izzt(17,2)=2
c        at(17)=16.d0*4.d0/(16.d0+4.d0)	!A : masse atomique r�duite
         at(17)=nucleot(12)*nucleot(4)/(nucleot(12)+nucleot(4))	!m r�duite

c r�action 18 : C12(C12,g)Mg24         	
	 nom_react(18)='C12(C12,g)Mg24'
	 nuc=2.d0*C12-Mg24
	 qt(18)=nuc
         izzt(18,1)=6
         izzt(18,2)=6
c        at(18)=12.d0*12.d0/(12.d0+12.d0)	!A : masse atomique r�duite
         at(18)=nucleot(7)*nucleot(7)/(nucleot(7)+nucleot(7))	!m r�duite

c r�action 19 : C12(C12,n)Mg23         	
	 nom_react(19)='C12(C12,n)Mg23'
	 nuc=C12-(n+Mg23)
	 qt(19)=nuc
         izzt(19,1)=6
         izzt(19,2)=6
c        at(19)=12.d0*12.d0/(12.d0+12.d0)	!A : masse atomique r�duite
         at(19)=nucleot(7)*nucleot(7)/(nucleot(7)+nucleot(7))	!m r�duite

c r�action 20 : C12(C12,p)Na23         	
	 nom_react(20)='C12(C12,p)Na23'
	 nuc=2.d0*C12-(H1+Na23)
	 qt(20)=nuc
         izzt(20,1)=6
         izzt(20,2)=6
c        at(20)=12.d0*12.d0/(12.d0+12.d0)	!A : masse atomique r�duite
         at(20)=nucleot(7)*nucleot(7)/(nucleot(7)+nucleot(7))	!m r�duite

c r�action 21 : C12(C12,a)Ne20         	
	 nom_react(21)='C12(C12,a)Ne20'
	 nuc=2.d0*C12-(a+Ne20)
	 qt(21)=nuc
         izzt(21,1)=6
         izzt(21,2)=6
c        at(21)=12.d0*12.d0/(12.d0+12.d0)	!A : masse atomique r�duite
         at(21)=nucleot(7)*nucleot(7)/(nucleot(7)+nucleot(7))	!m r�duite

c r�action 22 : C12(C12,2a)O16         	
	 nom_react(22)='C12(C12,2a)O16'
	 nuc=2.d0*C12-(2.d0*a+O16)
	 qt(22)=nuc
         izzt(22,1)=6
         izzt(22,2)=6
         at(22)=12.d0*12.d0/(12.d0+12.d0)	!A : masse atomique r�duite

c r�action 23 : N13(p,g)O14(e+nu)N14	
	 nom_react(23)='N13(p,g)N14'	!N13(p,g)O14(e+nu)N14 
	 nuc=N13+H1-N14
	 enu=1.037d0		!cours Pichon p.278
	 qt(23)=nuc-enu
         izzt(23,1)=7
         izzt(23,2)=1
c        at(23)=13.d0*1.d0/(13.d0+1.d0)	!A : masse atomique r�duite
         at(23)=nucleot(9)*nucleot(1)/(nucleot(9)+nucleot(1))	!m r�duite

c r�action 24 : O17(p,g)F18(e+nu)O18         	
	 nom_react(24)='O17(p,g)O18'	!O17(p,g)F18(e+nu)O18
	 nuc=O17+H1-O18
	 qd=F18-O18-2.d0*e
	 enu=qd*neu(qd)	!neutrino, approximation formule Chiu cours Pichon
	 qt(24)=nuc-enu
         izzt(24,1)=8
         izzt(24,2)=1
c        at(24)=17.d0*1.d0/(17.d0+1.d0)	!A : masse atomique r�duite
         at(24)=nucleot(13)*nucleot(1)/(nucleot(13)+nucleot(1))	!m r�duite

c r�action 25 : O18(p,g)F19         	
	 nom_react(25)='O18(p,g)F19'
	 nuc=O18+H1-F19
	 qt(25)=nuc
         izzt(25,1)=8
         izzt(25,2)=1
c        at(25)=18.d0*1.d0/(18.d0+1.d0)	!A : masse atomique r�duite
         at(25)=nucleot(14)*nucleot(1)/(nucleot(14)+nucleot(1))	!m r�duite

c r�action 26 : F19(p,a)O16         	
	 nom_react(26)='F19(p,a)O16'
	 nuc=F19+H1-(a+O16)
	 qt(26)=nuc
         izzt(26,1)=9
         izzt(26,2)=1
c        at(26)=19.d0*1.d0/(19.d0+1.d0)	!A : masse atomique r�duite
         at(26)=nucleot(15)*nucleot(1)/(nucleot(15)+nucleot(1))	!m r�duite

c r�action 27 : F19(p,g)Ne20         	
	 nom_react(27)='F19(p,g)Ne20'
	 nuc=F19+H1-Ne20
	 qt(27)=nuc
         izzt(27,1)=9
         izzt(27,2)=1
c        at(27)=19.d0*1.d0/(19.d0+1.d0)	!A : masse atomique r�duite
         at(27)=nucleot(15)*nucleot(1)/(nucleot(15)+nucleot(1))	!m r�duite

c r�action 28 : O18(p,a)N15         	
	 nom_react(28)='O18(p,a)N15'
	 nuc=O18+H1-(a+N15)
	 qt(28)=nuc
         izzt(28,1)=8
         izzt(28,2)=1
c        at(28)=18.d0*1.d0/(18.d0+1.d0)	!A : masse atomique r�duite
         at(28)=nucleot(14)*nucleot(1)/(nucleot(14)+nucleot(1))	!m r�duite

c r�action 29 : Ne20(a,g)Mg24         	
	 nom_react(29)='Ne20(a,g)Mg24'
	 nuc=Ne20+a-Mg24
	 qt(29)=nuc
         izzt(29,1)=16
         izzt(29,2)=4
c        at(29)=20.d0*4.d0/(20.d0+4.d0)	!A : masse atomique r�duite
         at(29)=nucleot(16)*nucleot(4)/(nucleot(16)+nucleot(4))	!m r�duite

c r�action 30 : C13(a,n)O16         	
	 nom_react(30)='C13(a,n)O16'
	 nuc=C13+a-(n+O16)
	 qt(30)=nuc
         izzt(30,1)=10
         izzt(30,2)=2
c        at(30)=13.d0*4.d0/(13.d0+4.d0)	!A : masse atomique r�duite
         at(30)=nucleot(8)*nucleot(4)/(nucleot(8)+nucleot(4))	!m r�duite

c r�action 31 : O17(a,n)Ne20         	
	 nom_react(31)='O17(a,n)Ne20'
	 nuc=O17+a-(n+Ne20)
	 qt(31)=nuc
         izzt(31,1)=8
         izzt(31,2)=2
c        at(31)=17.d0*4.d0/(17.d0+4.d0)	!A : masse atomique r�duite
         at(31)=nucleot(13)*nucleot(4)/(nucleot(13)+nucleot(4))	!m r�duite

c r�action 32 : N14(a,g)F18(e+nu)O18         	
	 nom_react(32)='N14(a,g)O18'	!N14(a,g)F18(e+nu)O18
	 nuc=N14+a-O18
	 qd=F18-O18-2.d0*e
	 enu=qd*neu(qd)	!neutrino, approximation formule Chiu cours Pichon
	 qt(32)=nuc-enu
         izzt(32,1)=7
         izzt(32,2)=2
c        at(32)=14.d0*4.d0/(14.d0+4.d0)	!A : masse atomique r�duite
         at(32)=nucleot(10)*nucleot(4)/(nucleot(10)+nucleot(4))	!m r�duite

c r�action 33 : O18(a,g)Ne22         	
	 nom_react(33)='O18(a,g)Ne22'
	 nuc=O18+a-Ne22
	 qt(33)=nuc
         izzt(33,1)=8
         izzt(33,2)=2
c        at(33)=18.d0*4.d0/(18.d0+4.d0)	!A : masse atomique r�duite
         at(33)=nucleot(14)*nucleot(4)/(nucleot(14)+nucleot(4))	!m r�duite

c r�action 34 : Ne22(a,n)Mg25         	
	 nom_react(34)='Ne22(a,n)Mg25'
	 nuc=Ne22+a-(n+Mg25)
	 qt(34)=nuc
         izzt(34,1)=10
         izzt(34,2)=2
c        at(34)=22.d0*4.d0/(22.d0+4.d0)	!A : masse atomique r�duite
         at(34)=nucleot(17)*nucleot(4)/(nucleot(17)+nucleot(4))	!m r�duite

c r�action 35 : Ne22(a,g)Mg26         	
	 nom_react(35)='Ne22(a,g)Mg26'
	 nuc=Ne22+a-Mg26
	 qt(35)=nuc
         izzt(35,1)=10
         izzt(35,2)=2
c        at(35)=22.d0*4.d0/(22.d0+4.d0)	!A : masse atomique r�duite
         at(35)=nucleot(17)*nucleot(4)/(nucleot(17)+nucleot(4))	!m r�duite
         
c r�action 36 : C12(a,n)O15(e+nu)N15	 	
	 nom_react(36)='C12(a,n)N15'	!C12(a,n)O15(e+nu)N15
	 nuc=C12+a-(n+N15)
	 enu=0.997d0
	 qt(36)=nuc-enu
         izzt(36,1)=6
         izzt(36,2)=2
c        at(36)=12.d0*4.d0/(12.d0+4.d0)	!A : masse atomique r�duite
         at(36)=nucleot(7)*nucleot(4)/(nucleot(7)+nucleot(4))	!m r�duite
        
c r�action 37 : Ne21(a,n)Mg24	 	
	 nom_react(37)='Ne21(a,n)Mg24'
	 nuc=Ne21+a-(n+Mg24)
	 qt(37)=nuc
         izzt(37,1)=10
         izzt(37,2)=2
         at(37)=21.d0*4.d0/(21.d0+4.d0)	!A : masse atomique r�duite
         at(37)=nucleot(22)*nucleot(4)/(nucleot(22)+nucleot(4))	!m r�duite

c r�action 38 : He4(an,g)Be9         	
	 nom_react(38)='He4(an,g)Be9'
	 nuc=He4+a+n-Be9
	 qt(38)=nuc
         izzt(38,1)=2
         izzt(38,2)=2
c        at(38)=4.d0*4.d0/(4.d0+4.d0)	!A : masse atomique r�duite???
         at(38)=nucleot(4)*nucleot(4)/(nucleot(4)+nucleot(4))	!m r�duite

c r�action 39 : Be9(p,d)2He4         	
	 nom_react(39)='Be9(p,d)2He4'
	 nuc=Be9+H1-(H2+2.d0*He4)
	 qt(39)=nuc
         izzt(39,1)=4
         izzt(39,2)=1
c        at(39)=9.d0*1.d0/(9.d0+1.d0)	!A : masse atomique r�duite
         at(39)=nucleot(25)*nucleot(1)/(nucleot(25)+nucleot(1))	!m r�duite
         	
c r�action 40 : Be9(a,n)C12
	 nom_react(40)='Be9(a,n)C12'
	 nuc=Be9+a-(n+C12)
	 qt(40)=nuc
         izzt(40,1)=4
         izzt(40,2)=2
c        at(40)=9.d0*4.d0/(9.d0+4.d0)	!A : masse atomique r�duite
         at(40)=nucleot(25)*nucleot(4)/(nucleot(25)+nucleot(4))	!m r�duite
         
c r�action 41 : Li6(p,He3)He4
	 nom_react(41)='Li6(p,He3)He4'
	 nuc=Li6+H1-(He3+He4)
	 qt(41)=nuc
         izzt(41,1)=3
         izzt(41,2)=1
c        at(41)=6.d0*1.d0/(6.d0+1.d0)	!A : masse atomique r�duite
         at(41)=nucleot(26)*nucleot(1)/(nucleot(26)+nucleot(1))	!m r�duite

c r�action 42 : Li6(p,g)Be7
	 nom_react(42)='Li6(p,g)Be7'
	 nuc=Li6+H1-Be7 ; qt(42)=nuc ; izzt(42,1)=3 ; izzt(42,2)=1
c        at(42)=6.d0*1.d0/(6.d0+1.d0)	!A : masse atomique r�duite
         at(42)=nucleot(26)*nucleot(1)/(nucleot(26)+nucleot(1))	!m r�duite

c r�action 43 : Be9(p,a)Li6
	 nom_react(43)='Be9(p,a)Li6'
	 nuc=Be9+H1-(He4+Li6) ; qt(43)=nuc ; izzt(43,1)=4 ; izzt(43,2)=1
c        at(43)=9.d0*1.d0/(9.d0+1.d0)	!A : masse atomique r�duite
         at(43)=nucleot(25)*nucleot(1)/(nucleot(25)+nucleot(1))	!m r�duite

c r�action 44 : B11(p,a)2He4
	 nom_react(44)='B11(p,a)2He4'
	 nuc=B11+H1-3.d0*He4 ; qt(44)=nuc ; izzt(44,1)=5 ; izzt(44,2)=1
c        at(44)=11.d0*1.d0/(11.d0+1.d0)	!A : masse atomique r�duite
         at(44)=nucleot(28)*nucleot(1)/(nucleot(28)+nucleot(1))	!m r�duite

c r�action 45 : B11(p,g)C12
	 nom_react(45)='B11(p,g)C12'
	 nuc=B11+H1-C12 ; qt(45)=nuc ; izzt(45,1)=5 ; izzt(45,2)=1
c        at(45)=11.d0*1.d0/(11.d0+1.d0)	!A : masse atomique r�duite
         at(45)=nucleot(28)*nucleot(1)/(nucleot(28)+nucleot(1))	!m r�duite

c	 do i=1,nreac_tot
c	  WRITE(*,10)i,nom_react(i),qt(i),izzt(i,1:2)
10	  FORMAT('r�action ',i2,1x,a20,es13.6,2es10.3,2i3)
c	 ENDDO

	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1005) ; WRITE(2,1005)
1005	  FORMAT('thermonuclear reactions of the compilation :')	 	 	
	 CASE DEFAULT
	  WRITE(*,5) ; WRITE(2,5)
5	  FORMAT('r�actions thermonucl�aires de la compilation :')	 
	 END SELECT	 
	 SELECT CASE (total)
	 CASE(0)
	  WRITE(*,1) ; WRITE(2,1)
1	  FORMAT('Caughland & Fowler')
	 CASE(1)
	  WRITE(*,2) ; WRITE(2,2)
2	  FORMAT('Aldeberger et al. 1997')
	 CASE(2)
	  WRITE(*,3) ; WRITE(2,3)
3	  FORMAT('NACRE')
	 CASE DEFAULT
	  PRINT*,'total ne peut prendre que les valeurs 0, 1, 2 et non:',
	1 total
	  PRINT*,'ARRET dans reac_c' ; STOP 
	 END SELECT
c	 pause	 
	 
	ENDIF
	
c	les facteurs t9

c	WRITE(*,2000)t

	t6=t*1.d-6 ; t9=t*1.d-9 ; lnt9=LOG(t9) ; lnt912=lnt9/2.d0
	lnt932=lnt912*3.d0 ; lnt913=lnt9/3.d0 ; lnt923=lnt913*2.d0
	lnt934=lnt932/2.d0 ; lnt953=lnt913*5.d0
	
c	WRITE(*,2000)t9,lnt9,lnt912,lnt932,lnt913,lnt923
c	pause	
	
	t913=EXP(lnt913) ; t912=EXP(lnt912); t923=EXP(lnt923)
	t932=EXP(lnt932) ; t953=EXP(lnt953)
	
c	WRITE(*,2000)t913,t912,t923,t932
c	pause		
	
c=======================================================================

c	cycle PP

c========================================================================

c r�action 1 : H1(p,e+ nu)H2	z0=1, z1=1

	IF(t9 > 3.d0)THEN
	 rt(1)=ln38
	ELSE
	
	 SELECT CASE (total)	
	 CASE(0)
	  ai(0)=1.d0		!coefficients du polynome en t913
	  ai(1)=0.123d0
	  ai(2)=1.09d0
	  ai(3)=0.938d0
	  CALL polyder(ai,3,0,t913,po)	!algorithme de Horner
	  rt(1)=LOG(po(0))+LOG(4.01d-15)-lnt923-3.38d0/t913-LOG(2.d0)  !/2!

	 CASE(1)       !formulation en S(0), S'(0), S"(0) Brun et al. 1998 
	  z1z2a=((izzt(1,1)*izzt(1,2))**2*at(1))**(1./3.)    !z1**2 z2**2 A
	  s0=4.d-25		!S(0)
	  sq=4.48d-24		!S'(0)
	  c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(1))
	  c2=4.2475d0*z1z2a
	  ai(0)=1.d0
	  ai(1)=9.81d-2/z1z2a
	  ai(2)=0.122d0*sq/s0*z1z2a
	  ai(3)=8.377d-2*sq/s0
	  CALL polyder(ai,3,0,t913,po)	!algorithme de Horner	
	  rt(1)=LOG(c1)-lnt923-c2/t913+LOG(po(0))-LOG(2.d0)	!/2!
c	  WRITE(*,2001)(ai(i),i=0,3),c1,c2

	 CASE(2)	!formule de NACRE
	  ai(0)=1.d0		!coefficients du polynome en t9
	  ai(1)=3.82d0
	  ai(2)=1.51d0
	  ai(3)=0.144d0
	  ai(4)=-1.14d-2	  
	  CALL polyder(ai,4,0,t9,po)	!algorithme de Horner
	  rt(1)=LOG(po(0))+LOG(4.08d-15)-lnt923-3.381d0/t913-LOG(2.d0) !/2!
	 END SELECT
	ENDIF
c	WRITE(*,2000)t,rt(1)	 
c	pause'H1(p,e+ nu)H2'

c r�action 2 : H2(p,g)He3	z0=1, z1=1

	SELECT CASE (total)
	CASE(0,1)
	 ai(0)=1.d0		!coefficients du polynome en t913
	 ai(1)=0.112d0
	 ai(2)=3.38d0
	 ai(3)=2.65d0
	 CALL polyder(ai,3,0,t913,po)	!algorithme de Horner
	 rt(2)=LOG(po(0))+LOG(2.24d+03)-lnt923-3.720d0/t913

	CASE(2)		!formule de NACRE
	 IF(t9 <= 0.11d0)THEN
	  ai(0)=1.d0		!coefficients du polynome en t9
	  ai(1)=14.3d0
	  ai(2)=-90.5d0
	  ai(3)=395.d0
	  CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	  rt(2)=LOG(po(0))+LOG(1.81d3)-lnt923-3.721d0/t913
	 ELSE
	  ai(0)=1.d0		!coefficients du polynome en t913
	  ai(1)=3.96d0
	  ai(2)=0.116d0
	  CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	  rt(2)=LOG(po(0))+LOG(2.58d3)-lnt923-3.721d0/t913
	 ENDIF
	END SELECT
c	WRITE(*,2000)t,rt(2)
c	pause'H2(p,g)He3'
	
c r�action 3 : He3(He3,2p)He4	z0=2, z1=2

	SELECT CASE (total)
	CASE(0)	
	 ai(0)=1.d0		!coefficients du polynome en t913
	 ai(1)=0.034d0
	 ai(2)=-0.522d0
	 ai(3)=-0.124d0
	 ai(4)=0.353d0
	 ai(5)=0.213d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 rt(3)=LOG(po(0))+LOG(6.04d10)-lnt923-12.276d0/t913-LOG(2.d0)  !/2!
c	 PRINT*,'r3'	
c	 WRITE(*,2000)rt(3)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Brun et al. 1998	 
	 z1z2a=((izzt(3,1)*izzt(3,2))**2*at(3))**(1./3.)     !z1**2 z2**2 A
	 s0=5.4d0	!S(0)
	 sq=-4.1d0	!S'(0)
	 ss=4.6d0	!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(3))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(3)=LOG(c1)-lnt923-c2/t913+LOG(po(0))-LOG(2.d0)	!/2!
c	 WRITE(*,2000)rt(3)	 
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'He3(He3,2p)He4'

	CASE(2)		!NACRE formule
	 ai(0)=1.d0		!coefficients du polynome en t9
	 ai(1)=-0.135d0
	 ai(2)=2.54d-2
	 ai(3)=-1.29d-3
	 CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	 rt(3)=LOG(po(0))+LOG(5.59d10)-lnt923-12.277d0/t913-LOG(2.d0) !/2!
c	 PRINT*,'r3'	
c	 WRITE(*,2000)rt(3)
	 
	END SELECT
c	WRITE(*,2000)t,rt(3)	 
c	pause'He3(He3,2p)He4'
	
c r�action 4 : He3(a,g)Be7	z0=2, z1=2

	SELECT CASE (total)
	CASE(0)
	 lnt9a=lnt9-LOG(1.d0+4.95d-02*t9)
	 rt(4)=LOG(5.61d6)+5.d0/6.d0*lnt9a-lnt932-12.826d0/EXP(lnt9a/3.d0)
c	 PRINT*,'r4'	
c	 WRITE(*,2000)rt(4)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Brun et al. 1998	 
	 z1z2a=((izzt(4,1)*izzt(4,2))**2*at(4))**(1.d0/3.d0) !z1**2 z2**2 A
	 s0=5.3d-4		!S(0)
	 sq=-3.d-4		!S'(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(4))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 CALL polyder(ai,3,0,t913,po)	!algorithme de Horner	
	 rt(4)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(4)	 
c	 WRITE(*,2001)(ai(i),i=0,3),c1,c2
c	 pause'He3(a,g)Be7'

	CASE(2)		!NACRE formule	 
	 ai(0)=1.d0		!coefficients du polynome en t9
	 ai(1)=-0.307d0
	 ai(2)=8.81d-2
	 ai(3)=-1.06d-2
	 ai(4)=4.46d-4	 
	 CALL polyder(ai,4,0,t9,po)	!algorithme de Horner
	 rt(4)=LOG(po(0))+LOG(5.46d6)-lnt923-12.827d0/t913
	END SELECT
c	WRITE(*,2000)t,rt(4)	 
c	pause'He3(a,g)Be7'

c	r�action 5 : Li7(p,a)He4	z0=3, z1=1

	SELECT CASE (total)
	CASE(0,1)
	 lnt9a=lnt9-LOG(1.d0+.759d0*t9)
	 v0=LOG(1.096d9)-lnt923-8.472d0/t913
	 v1=LOG(4.830d8)+5.d0/6.d0*lnt9a-lnt932-8.472d0/EXP(lnt9a/3.d0)
	 v2=LOG(1.06d10)-lnt932-30.442d0/t9
	 
c r�action Li7(p,g)Be8(a)He4 (B. Pichon)

	 ai(0)=1.d0		!coefficients du polynome en t9
	 ai(1)=0.049d0
	 ai(2)=2.498d0
	 ai(3)=0.86d0
	 ai(4)=3.518d0
	 ai(5)=3.080d0	 	 
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v3=LOG(po(0))+LOG(1.56d5)-lnt923-8.472/t913-(t9/1.696d0)**2
	 v4=LOG(1.55d6)-lnt932-4.478d0/t9

	 rt(5)=EXP(v0)-EXP(v1)+EXP(v2)+EXP(v3)+EXP(v4)
	 IF(rt(5) > 0.d0)THEN
	  rt(5)=LOG(rt(5))
	 ELSE
	  rt(5)=ln38
	 ENDIF
c	 PRINT*,'r5'	
c	 WRITE(*,2000)rt(5)

	CASE(2)			!NACRE formule	 
	 ai(0)=1.d0		!coefficients du polynome en t9
	 ai(1)=1.05d0
	 ai(2)=-0.653d0
	 ai(3)=0.185d0
	 ai(4)=-2.12d-2
	 ai(5)=9.3d-4	 	 
	 CALL polyder(ai,5,0,t9,po)	!algorithme de Horner
	 v1=LOG(po(0))+LOG(7.2d8)-lnt923-8.473d0/t913-(t9/6.5d0)**2
	 v2=LOG(9.85d6)+0.576d0*lnt9-10.415d0/t9
	 
c	 on ajoute Li7(p,g)Be8(a)He4 (B. Pichon)

	 ai(0)=1.d0		!coefficients du polynome en t9
	 ai(1)=-1.47d0
	 ai(2)=4.43d0
	 CALL polyder(ai,2,0,t9,po)	!algorithme de Horner	 
	 v3=LOG(po(0))+LOG(1.75d7)-lnt923-8.473d0/t913-(t9/0.8d0)**2
	 v4=LOG(1.6d6)-lnt923-4.441d0/t9
	 v5=LOG(4.32d4)+0.309d0*lnt9-2.811d0/t9
	 
	 rt(5)=EXP(v1)+EXP(v2)+EXP(v3)+EXP(v4)+EXP(v5)
	 IF(rt(5) > 0.d0)THEN
	  rt(5)=LOG(rt(5))
	 ELSE
	  rt(5)=ln38
	 ENDIF	 
	END SELECT
c	WRITE(*,2000)t,rt(5)	 
c	pause'Li7(p,a)He4'

c r�action 6 : Be7(e-,nu g)Li7     pas de saturation t9 >= .001

	SELECT CASE (total)
	CASE(0)
	 ai(0)=1.d0		!coefficients du polynome en t913
	 ai(1)=-0.537d0
	 ai(2)=3.86d0
	 CALL polyder(ai,2,0,t913,po)	!algorithme de Horner
	 po(0)=LOG(po(0)+.0027/t9*EXP(2.515d-3/t9))
	 rt(6)=po(0)+LOG(1.34d-10)-lnt912
	 IF(t9 > 3.d0)rt(6)=ln38
c	 PRINT*,'r6'	
c	 WRITE(*,2000)rt(6)

	CASE(1,2)		!pour le soleil Adelberger et al. 1998
	 rt(6)=LOG(5.6d-9/sqrt(t6)*(1.d0+0.004d0*(t6-16.d0)))	!et NACRE
c	 PRINT*,'r6'	
c	 WRITE(*,2000)rt(6)
	END SELECT
c	WRITE(*,2000)t,rt(6)	 
c	pause'Be7(e-,nu g)Li7'
		
c r�action 7 : Be7(p,g)B8(e+ nu)Be8(a)He4	z0=4, z1=1

	SELECT CASE (total)
	CASE(0)
	 v0=EXP(LOG(3.11d5)-lnt923-10.262d0/t913)
	 v1=EXP(LOG(2.53d3)-lnt932-7.306d0/t9)
	 rt(7)=v0+v1		!Be7(H,g)B8
	 IF(rt(7) > 0.d0)THEN
	  rt(7)=LOG(rt(7))
	 ELSE
	  rt(7)=ln38
	 ENDIF
c	 PRINT*,'r7'	
c	 WRITE(*,2000)rt(7)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Brun et al. 1998	 
	 z1z2a=((izzt(7,1)*izzt(7,2))**2*at(7))**(1./3.)     !z1**2 z2**2 A
	 s0=1.9d-5		!S(0)
	 sq=-1.35d-5		!S'(0)
	 ss=7.33d-5		!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(7))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(7)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(7)	 
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'Be7(p,g)B8(e+ nu)Be8(a)He4'

	CASE(2)		!NACRE formules	 
	 ai(0)=1.d0		!coefficients du polynome en t9
	 ai(1)=-5.11d-2
	 ai(2)=4.68d-2
	 ai(3)=-6.6d-3
	 ai(4)=3.12d-4
	 CALL polyder(ai,4,0,t9,po)	!algorithme de Horner
	 v1=LOG(po(0))+LOG(2.61d5)-lnt923-10.264d0/t913
	 v2=LOG(2.05d3)-lnt932-7.345d0/t9
	 rt(7)=LOG(EXP(v1)+EXP(v2))
	END SELECT
c	WRITE(*,2000)t,rt(7)	 
c	pause'Be7(p,g)B8(e+ nu)Be8(a)He4'

c======================================================================

c	cycle CNO

c======================================================================

c r�action 8 : C12(p,g)N13(e+ nu)C13	z0=6, z1=1

	SELECT CASE (total)
	CASE(0)
	 ai(0)=1.d0
	 ai(1)=0.03d0
	 ai(2)=1.19d0
	 ai(3)=0.254d0
	 ai(4)=2.06d0
	 ai(5)=1.12d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(2.04d7)-lnt923-13.690/t913-(t9/1.5)**2)
	 v1=EXP(LOG(1.08d5)-lnt932-4.925/t9)
	 v2=EXP(LOG(2.15d5)-lnt932-18.179/t9)
	 rt(8)=v0+v1+v2
	 IF(rt(8) > 0.d0)THEN
	  rt(8)=LOG(rt(8))
	 ELSE
	  rt(8)=ln38
	 ENDIF	
c	 PRINT*,'r8'
c	 WRITE(*,2000)rt(8)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Adelberger et al. 1998	 
	 z1z2a=((izzt(8,1)*izzt(8,2))**2*at(8))**(1./3.)     !z1**2 z2**2 A
	 s0=1.34*1.d-3		!S(0)
	 sq=2.6d-3		!S'(0)
	 ss=8.3d-5*1.d3		!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(8))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(8)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(8)	 
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'C12(p,g)N13(e+ nu)C13'

	CASE(2)		!NACRE formule	 
	 ai(0)=1.d0
	 ai(1)=9.89d0
	 ai(2)=-59.8d0
	 ai(3)=266.d0
	 CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(2.d7)-lnt923-13.692d0/t913-(t9/0.46d0)**2)
	 v1=EXP(LOG(1.d5)-lnt932-4.913d0/t9)
	 v2=EXP(LOG(4.24d5)-lnt932-21.62d0/t9)
	 rt(8)=v0+v1+v2	 
	 IF(rt(8) > 0.d0)THEN
	  rt(8)=LOG(rt(8))
	 ELSE
	  rt(8)=ln38
	 ENDIF	

	END SELECT
c	WRITE(*,2000)t,rt(8)	 
c	pause'C12(p,g)N13(e+ nu)C13'

c r�action 9 : C13(p,g)N14	z0=6, z1=1

	SELECT CASE (total)
	CASE(0)
	 ai(0)=1.d0
	 ai(1)=0.03d0
	 ai(2)=0.958d0
	 ai(3)=0.204d0
	 ai(4)=1.39d0
	 ai(5)=0.753d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(8.01d7)-lnt923-13.717/t913-(t9/2.)**2)
	 v1=EXP(LOG(1.21d6)-6./5.*lnt9-5.701/t9)
	 rt(9)=v0+v1
	 IF(rt(9) > 0.d0)THEN
	  rt(9)=LOG(rt(9))
	 ELSE
	  rt(9)=ln38
	 ENDIF
c	 PRINT*,'r9'
c	 WRITE(*,2000)rt(9)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Adelberger et al. 1998
	 z1z2a=((izzt(9,1)*izzt(9,2))**2*at(9))**(1./3.)     !z1**2 z2**2 A
	 s0=7.6*1.d-3		!S(0)
	 sq=-7.8d-3		!S'(0)
	 ss=7.3d-4*1.d3		!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(9))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(9)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(9)	 
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'C13(p,g)N14'

	CASE(2)		!NACRE formule
	 ai(0)=1.d0
	 ai(1)=3.56d0
	 CALL polyder(ai,1,0,t9,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(9.57d7)-lnt923-13.72d0/t913-t9**2)
	 v1=EXP(LOG(1.5d6)-lnt932-5.93d0/t9)
	 v3=EXP(LOG(6.83d5)-0.864d0*lnt9-12.057d0/t9)
	 rt(9)=(v0+v1+v3)*(1.d0-EXP(LOG(2.07d0)-37.938/t9))	!N_A(sv)_tt
	 IF(rt(9) > 0.d0)THEN
	  rt(9)=LOG(rt(9))
	 ELSE
	  rt(9)=ln38
	 ENDIF
	END SELECT
c	WRITE(*,2000)t,rt(9)	 
c	pause'C13(p,g)N14'

c r�action 10 : N14(p,g)O15(e+ nu)N15	z0=7, z1=1

	SELECT CASE (total)
	CASE(0)
	 ai(0)=1.d0
	 ai(1)=0.027d0
	 ai(2)=-0.778d0
	 ai(3)=-0.149d0
	 ai(4)=0.261d0
	 ai(5)=0.127d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(4.9d7)-lnt923-15.228/t913-(t9/3.294)**2)
	 v1=EXP(LOG(2.37d3)-lnt932-3.011/t9)
	 v2=EXP(LOG(2.19d4)-12.53/t9)
	 rt(10)=v0+v1+v2
	 IF(rt(10) > 0.d0)THEN
	  rt(10)=LOG(rt(10))
	 ELSE
	  rt(10)=ln38
	 ENDIF
c	 PRINT*,'r10'
c	 WRITE(*,2000)rt(10)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Brun et al. 1998	 
	 z1z2a=((izzt(10,1)*izzt(10,2))**2*at(10))**(1./3.)   !z1**2z2**2 A
	 s0=3.5d-3		!S(0)
	 sq=-0.0128d0		!S'(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(10))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 CALL polyder(ai,3,0,t913,po)	!algorithme de Horner	
	 rt(10)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(10)	 
c	 WRITE(*,2001)(ai(i),i=0,3),c1,c2
c	 pause'N14(p,g)O15(e+ nu)N15'

	CASE(2)		!NACRE formulation
	 IF(nom_nuc_cpl == 'NACRE_LUNA') THEN ! taux de LUNA pour N14
c	 Imbriani et al., Eur. Phys. J. A. 25, 455-466 (2005), eq. (15), pag. 465
	  ai(0)=0.782d0
	  ai(1)=-1.5d0
	  ai(2)=17.97d0
	  ai(3)=-3.32d0
	  CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(po(0))+LOG(3.12d7)-lnt923-15.193d0/t913-(t9/0.486d0)**2)
	  v1=EXP(LOG(2.11d3)-lnt932-2.998d0/t9)
	  v2=EXP(LOG(8.42d2)+0.0682d0*lnt9-4.891d0/t9)
	 ELSE ! taux de NACRE
	  ai(0)=1.d0
	  ai(1)=-2.d0
	  ai(2)=3.41d0
	  ai(3)=-2.43d0
	  CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(po(0))+LOG(4.83d7)-lnt923-15.231d0/t913-(t9/0.8d0)**2)
	  v1=EXP(LOG(2.36d3)-lnt932-3.01d0/t9)
	  v2=EXP(LOG(6.72d3)+0.38d0*lnt9-9.53d0/t9)
	 ENDIF ! NACRE ou LUNA

	 rt(10)=v0+v1+v2
c********************************************************************
C	 rt(10)=.65d0*rt(10)	! 70% of NACRE
c********************************************************************
	 IF(rt(10) > 0.d0)THEN
	  rt(10)=LOG(rt(10))
	 ELSE
	  rt(10)=ln38
	 ENDIF
	END SELECT
c	WRITE(*,2000)t,rt(10)	 
c	pause'N14(p,g)O15(e+ nu)N15'
	
c r�action 11 : N15(p,g)O16	z0=7, z1=1

	SELECT CASE (total)
	CASE(0)
	 ai(0)=1.d0
	 ai(1)=0.027d0
	 ai(2)=0.219d0
	 ai(3)=0.042d0
	 ai(4)=6.83d0
	 ai(5)=3.32d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(9.78d8)-lnt923-15.251/t913-(t9/0.450)**2)
	 v1=EXP(LOG(1.11d4)-lnt932-3.328/t9)
	 v2=EXP(LOG(1.49d4)-lnt932-4.665/t9)
	 v3=EXP(LOG(3.80d6)-lnt932-11.048/t9)
	 rt(11)=v0+v1+v2+v3
	 IF(rt(11) > 0.)THEN
	  rt(11)=LOG(rt(11))
	 ELSE
	  rt(11)=ln38
	 ENDIF
c	 PRINT*,'r11'
c	 WRITE(*,2000)rt(11)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Adelberger et al. 1998
	 z1z2a=((izzt(11,1)*izzt(11,2))**2*at(11))**(1./3.)   !z1**2z2**2 A
	 s0=64.*1.d-3		!S(0)
	 sq=2.1d-2		!S'(0)
	 ss=4.1d-3*1.d3		!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(11))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(11)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(11)	 
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'N15(p,g)O16'

	CASE(2)		!NACRE formulation
	 IF(t9 <= 3.5d0)THEN
	  ai(0)=1.d0
	  ai(1)=6.15d0
	  ai(2)=16.4d0
	  CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(po(0))+LOG(1.08d9)-lnt923-15.254d0/t913
	1 -(t9/0.34d0)**2)
	  v1=EXP(LOG(9.23d3)-lnt932-3.597d0/t9)
	  v2=EXP(LOG(3.27d6)-lnt932-11.024/t9)
	  rt(11)=v0+v1+v2
	 ELSE
	  rt(11)=EXP(LOG(3.54d4)+0.095d0*lnt9-2.306d0/t9)
	 ENDIF
	 IF(rt(11) > 0.d0)THEN
	  rt(11)=LOG(rt(11))
	 ELSE
	  rt(11)=ln38
	 ENDIF
c	 PRINT*,'r11'
c	 WRITE(*,2000)rt(11)
	END SELECT
c	WRITE(*,2000)t,rt(11)	 
c	pause'N15(p,g)O16'

c r�action 12 : N15(p,a)C12	z0=7, z1=1

	SELECT CASE (total)
	CASE(0)
	 ai(0)=1.d0
	 ai(1)=0.027d0
	 ai(2)=2.62d0
	 ai(3)=0.501d0
	 ai(4)=5.36d0
	 ai(5)=2.60d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(MAX(ln38,LOG(po(0))+LOG(1.08d12)-lnt923-15.251d0/t913
	1	-(t9/0.522d0)**2))
	 v1=EXP(MAX(ln38,LOG(1.19d8)-lnt932-3.676d0/t9))
	 v2=EXP(MAX(ln38,LOG(5.41d8)-lnt912-8.926d0/t9))
	 v3=EXP(MAX(ln38,LOG(4.72d8)-lnt932-7.721d0/t9))
	 v4=EXP(MAX(ln38,LOG(2.2d9)-lnt932-11.418d0/t9))
	 rt(12)=v0+v1+v2+zero21c*(v3+v4) 
	 IF(rt(12) > 0.d0)THEN
	  rt(12)=LOG(rt(12))
	 ELSE
	  rt(12)=ln38
	 ENDIF
c	 PRINT*,'r12'
c	 WRITE(*,2000)rt(12)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Adelberger et al. 1998
	 z1z2a=((izzt(12,1)*izzt(12,2))**2*at(12))**(1./3.)  !z1**2z2**2 A
	 s0=6.75d4*1.d-3	!S(0)
	 sq=310.d0		!S'(0)
	 ss=12.d0*1.d3		!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(12))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(12)=LOG(c1)-lnt923-c2/t913+LOG(po(0))
c	 WRITE(*,2000)rt(12)	 
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'N15(p,a)C12'

	CASE(2)		!NACRE formulation
	 IF(t9 <= 2.5d0)THEN
	  ai(0)=1.d0
	  ai(1)=4.95d0
	  ai(2)=143.d0
	  CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(po(0))+LOG(1.12d12)-lnt923-15.253/t913-(t9/0.28d0)**2)
	  v1=EXP(LOG(1.01d8)-lnt932-3.643d0/t9)
	  v2=EXP(LOG(1.19d9)-lnt932-7.406d0/t9)
	  rt(12)=v0+v1+v2
	 ELSE
	  rt(12)=EXP(LOG(4.17d7)+0.917d0*lnt9-3.292d0/t9)
	 ENDIF
	 IF(rt(12) > 0.d0)THEN
	  rt(12)=LOG(rt(12))
	 ELSE
	  rt(12)=ln38
	 ENDIF
c	 PRINT*,'r12'
c	 WRITE(*,2000)rt(12)
	END SELECT
c	WRITE(*,2000)t,rt(12)	 
c	pause'N15(p,a)C12'

c r�action 13 : O16(p,g)F17(e+ nu)O17	z0=8, z1=1

	SELECT CASE (total)
	CASE(0)
	 rt(13)=LOG(1.5d8)-LOG(t923*(1.d0+2.13d0*(1.d0-EXP(-0.728d0*
	1 t923))))-16.692d0/t913
c	 PRINT*,'r13'
c	 WRITE(*,2000)rt(13)

	CASE(1)	!formulation en S(0), S'(0), S"(0) Adelberger et al. 1998
	 z1z2a=((izzt(13,1)*izzt(13,2))**2*at(13))**(1./3.)  !z1**2z2**2 A
	 s0=9.4*1.d-3		!S(0)
	 sq=-2.4d-2		!S'(0)
	 ss=5.7d-5*1.d3		!S"(0)
	 c1=7.8324d9*sqrt(z1z2a)*s0/sqrt(at(13))
	 c2=4.2475d0*z1z2a
	 ai(0)=1.d0
	 ai(1)=9.81d-2/z1z2a
	 ai(2)=0.122d0*sq/s0*z1z2a
	 ai(3)=8.377d-2*sq/s0
	 ai(4)=7.442d-3*ss/s0*z1z2a**2
	 ai(5)=1.299d-2*ss/s0*z1z2a
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner	
	 rt(13)=LOG(c1)-lnt923-c2/t913+LOG(po(0))-LOG(2.d0)	!/2!
c	 WRITE(*,2000)rt(13)
c	 WRITE(*,2001)(ai(i),i=0,5),c1,c2
c	 pause'O16(p,g)F17(e+ nu)O17'

	CASE(2)		!NACRE formulation
	 rt(13)=LOG(7.37d7)-16.696d0/t913-0.82d0*lnt9
	1 +LOG(1.d0+EXP(LOG(202.d0)-70.348d0/t9-0.161d0*t9))!N_A(sv)_tt
	END SELECT
c	WRITE(*,2000)t,rt(13)	 
c	pause'O16(p,g)F17(e+ nu)O17'

c r�action 14 : O17(p,a)N14	z0=8, z1=1

	SELECT CASE (total)
	CASE(0,1)
	 ai(0)=1.d0	!modif Landre et al. A&A 1990, 240, 85
	 ai(1)=0.025d0
	 ai(2)=5.39d0
	 ai(3)=0.940d0
	 ai(4)=13.5d0
	 ai(5)=5.98d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(po(0))+LOG(1.53d7)-lnt923-16.712/t913-(t9/0.565)**2)
	 v1=EXP(LOG(2.92d6)+lnt9-4.247/t9)
c	 v2=EXP(LOG(4.81d10)+lnt9-16.712/t913-(t9/0.04)**2) !Cauhglan
c	 v3=EXP(LOG(5.05d-5)-lnt932-.723/t9)
c	 v4=EXP(LOG(1.31d1)-lnt932-1.961/t9)
c	 v5=0.d0

	 v5=EXP(LOG(1.78d5)-lnt923-16.67d0/t913-		!Landre
	1	(2.d0*LOG(EXP(LOG(0.479)+lnt923)+3.12d-3)))
	 v2=EXP(LOG(2.8d11)+lnt9-16.67/t913-(t9/0.04)**2) 	!Landre
	 v3=EXP(LOG(2.94d-3)-lnt932-.767/t9)
	 v4=EXP(LOG(9.8d1)-lnt932-2.077/t9)
	 rt(14)=v0+v1+v5+(v3+v2+v4)*zero21l
	 IF(rt(14) > 0.d0)THEN
	  rt(14)=LOG(rt(14))
	 ELSE
	  rt(14)=ln38
	 ENDIF
	 
	CASE(2)	!NACRE formulation
	 IF(t9 <= 6.d0)THEN
	  ai(0)=1.d0
	  ai(1)=-80.31d0
	  ai(2)=2211.d0
	  CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(po(0))+LOG(9.2d8)-lnt923-16.715d0/t913-(t9/0.06d0)**2)
	  v1=EXP(LOG(9.13d-4)-lnt932-0.7667d0/t9)
	  v2=EXP(LOG(9.68d0)-lnt932-2.083d0/t9)
	  v3=EXP(LOG(8.13d6)-lnt932-5.685d0/t9)
	  v4=EXP(LOG(1.85d6)+1.591d0*lnt9-4.848d0/t9)
	  rt(14)=v0+v1+v2+v3+v4
	 ELSE
	  rt(14)=EXP(LOG(8.73d6)+0.95d0*lnt9-7.508d0/t9)
	 ENDIF	 
	 rt(14)=rt(14)*(1.d0+EXP(LOG(1.033d0)-10.034/t9-0.165*t9))
	 IF(rt(14) > 0.d0)THEN
	  rt(14)=LOG(rt(14))
	 ELSE
	  rt(14)=ln38
	 ENDIF
	END SELECT
c	WRITE(*,2000)t,rt(14)	 
c	pause'O17(p,a)N14'
	 
c====================================================================

c	cycle 3 alpha

c=====================================================================

c r�action 15 : He4(2a,g)C12	z0=2, z1=2, z2=2

	IF(t9 <= 0.08d0)THEN
	 ai(0)=1.d0
	 ai(1)=0.031d0
	 ai(2)=8.009d0
	 ai(3)=1.732d0
	 ai(4)=49.883d0
	 ai(5)=27.426d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(7.4d5)-lnt932-1.0663/t9)
	 v1=EXP(LOG(po(0))+LOG(4.164d9)-lnt923-13.49d0/t913
	1 -(t9/0.098d0)**2)
	 he4abe8=LOG(v0+v1)		!ecran?
c	 PRINT*,'he4abe8'
c	 WRITE(*,2000)he4abe8

	 ai(0)=1.d0
	 ai(1)=0.018d0
	 ai(2)=5.249d0
	 ai(3)=0.65d0
	 ai(4)=19.176d0
	 ai(5)=6.034d0
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner

	 v0=EXP(LOG(1.3d2)-lnt932-3.3364d0/t9)
	 v1=EXP(LOG(po(0))+LOG(2.51d7)-lnt923-23.57d0/t913
	1 -(t9/0.0235d0)**2)
	 be8agc12=LOG(v0+v1)
c	 PRINT*,'be8agc12'
c	 WRITE(*,2000)be8agc12
	 v0=LOG(2.9d-16)+he4abe8+be8agc12+
	1 LOG(0.01d0+0.2d0*(1.d0+4.d0*EXP(-(0.025d0/t9)**3.263d0)))-
	2 LOG(1.d0+4.d0*EXP(-(t9/0.025d0)**9.227d0))
	 v1=LOG(1.35d-7)-lnt932-24.811d0/t9
	 v0=MAX(v0,ln38) ; v1=MAX(v1,ln38)
	 rt(15)=EXP(v0)+zero21c*EXP(v1)

	ELSE
	 v0=EXP(LOG(2.79d-8)-lnt9*3.d0-4.4027d0/t9)
	 v1=EXP(LOG(1.35d-7)-lnt932-24.811d0/t9)
	 rt(15)=v0+zero21c*v1
	ENDIF
	rt(15)=rt(15)/6.d0		!1/3!
	IF(rt(15) > 0.d0)THEN
	 rt(15)=LOG(rt(15))
	ELSE
	 rt(15)=ln38
	ENDIF
c	PRINT*,'r15'
c	WRITE(*,2000)rt(15)

c r�action 16 : C12(a,g)O16	z0=6, z1=2
c	multiplication par 1.7 de la r�action 16 C12(a,g)O16
c	(Weaver & Woosley 1993, Phys. Rept., 227, 65)

	v0=EXP(LOG(1.04d8)-2.*lnt9-2.d0*LOG(1.d0+0.0489d0/t923)
	1	-32.120d0/t913-(t9/3.496d0)**2)
	v1=EXP(LOG(1.76d8)-2.d0*lnt9-2.d0*LOG(1.d0+0.2654d0/t923)
	1 -32.120d0/t913)
	v2=EXP(LOG(1.25d3)-lnt932-27.499d0/t9)
	v3=EXP(LOG(1.43d-2)+5.d0*lnt9-15.541d0/t9)
	rt(16)=v0+v1+v2+v3
	IF(rt(16) > 0.d0)THEN
	 rt(16)=LOG(rt(16)*1.7)	!Weaver & Woosley 1993, Phys.Rept., 227, 65
	ELSE
	 rt(16)=ln38
	ENDIF
c	PRINT*,'r16'
c	WRITE(*,2000)rt(16)

c	taux de Caughlan et Fowler 1985 qui est meilleur
c	Nomoto et al. astro-ph/9706024

	v0=EXP(LOG(2.96d8)-2.d0*lnt9-2.d0*LOG(1.d0+0.0489d0/t923)
	1	-32.120d0/t913-(t9/3.496d0)**2)
	v1=EXP(LOG(3.14d8)-2.d0*lnt9-2.d0*LOG(1.d0+0.2654d0/t923)
	1 -32.120d0/t913)
	v2=EXP(LOG(1.25d3)-lnt932-27.499d0/t9)
	v3=EXP(LOG(1.43d-2)+5.d0*lnt9-15.541d0/t9)

	rt(16)=v0+v1+v2+v3
	IF(rt(16) > 0.d0)THEN
	 rt(16)=LOG(rt(16))
	ELSE
	 rt(16)=ln38
	ENDIF
c	PRINT*,'r16'
c	WRITE(*,2000)rt(16)

c r�action 17 : O16(a,g)Ne20	z0=8, z1=2

	v0=EXP(LOG(9.37d9)-lnt923-39.757d0/t913-(t9/1.586d0)**2)
	v1=EXP(LOG(6.21d1)-lnt932-10.297d0/t9)
	v2=EXP(LOG(5.38d2)-lnt932-12.226d0/t9)
	v3=EXP(LOG(1.3d1)+2.d0*t9-20.093d0/t9)
	rt(17)=v0+v1+v2+v3
	IF(rt(17) > 0.d0)THEN
	 rt(17)=LOG(rt(17))
	ELSE
	 rt(17)=ln38
	ENDIF
c	PRINT*,'r17'
c	WRITE(*,2000)rt(17)

c=======================================================================

c	cycle carbone, le 1/2! est dans rt(18)

c========================================================================

c r�action 18 : C12(C12,g)Mg24	z0=6, z1=6

	lnt9a=lnt9-LOG(1.d0+0.0396d0*t9)
	t9a13=EXP(lnt9a/3.d0)

	rt(18)=LOG(4.27d26)+5.d0/6.d0*lnt9a-lnt932-
	1	84.165d0/t9a13-2.12d-3*t9**3-LOG(2.d0)		!1/2!
	
	rt(18)=EXP(rt(18))
c	PRINT*,'r18'
c	WRITE(*,2000)rt(18)

c r�action 19 : C12(C12,n)Mg23	z0=6, z1=6 "N-yield"

	IF(t9 < 1.75d0)THEN
	 rt(19)=ln38
	ELSEIF(t9 < 3.3d0)THEN
	 rt(19)=rt(18)*0.05d0
	ELSEIF(t9 < 6.d0)THEN
	 rt(19)=rt(18)*0.07d0	
	ELSE
	 PRINT*,'pas de taux pour la r�action 19 pour t=',t9
	 rt(19)=1.d-100
	ENDIF
c	PRINT*,'r19'
c	WRITE(*,2000)rt(19)	
	
c r�action 20 : C12(C12,p)Na23	z0=6, z1=6 "P-yield"

	IF(t9 < 1.75)THEN
	 rt(20)=rt(18)*0.44d0
	ELSEIF(t9 < 3.3d0)THEN
	 rt(20)=rt(18)*0.45d0
	ELSEIF(t9 < 6.d0)THEN
	 rt(20)=rt(18)*0.40d0	
	ELSE
	 PRINT*,'pas de taux pour la r�action 20 pour t=',t9
	 rt(20)=1.d-100
	ENDIF
c	PRINT*,'r20'
c	WRITE(*,2000)rt(20)	

c r�action 21 : C12(C12,a)Ne20	z0=6, z1=6	"A-yield"

	IF(t9 < 1.75d0)THEN
	 rt(21)=rt(18)*0.56d0
	ELSEIF(t9 < 3.3)THEN
	 rt(21)=rt(18)*0.5d0
	ELSEIF(t9 < 6.d0)THEN
	 rt(21)=rt(18)*0.53d0	
	ELSE
	 PRINT*,'pas de taux pour la r�action 20 pour t=',t9
	 rt(21)=1.d-100
	ENDIF
c	PRINT*,'r21'
c	WRITE(*,2000)rt(21)	
	
	do i=18,21
	 rt(i)=MAX(rt(i),1.d-100)
	 rt(i)=LOG(rt(i))
	ENDDO
c	PRINT*,'apres log'

c r�action 22 : C12(C12,2a)O16	z0=6, z1=6

c!!!!!!!!!!!!!!!A CREER!!!!!!!!!!!!!!!!!!

	rt(22)=ln38	
c	PRINT*,'r22'
c	WRITE(*,2000)rt(22)

c=======================================================================

c	autres r�actions

c========================================================================

c r�action 23 : N13(p,g)O14(e+nu)N14 z0=7, z1=1

	ai(0)=1.d0
	ai(1)=0.658d0
	ai(2)=0.379d0
	CALL polyder(ai,2,0,t912,po)	!algorithme de Horner
	rt(23)=LOG(po(0))+LOG(6.74d7)-68.762/t9
c	PRINT*,'r23'
c	WRITE(*,2000)rt(23)

c r�action 24: O17(p,g)F18(e+nu)O18 (Landre V. et al 1990) z0=8, z1=8
c	modif Landre et al. A&A 1990, 240, 85

	ai(0)=1.d0
	ai(1)=0.025d0
	ai(2)=-0.051d0
	ai(3)=-8.82d-3
	CALL polyder(ai,3,0,t913,po)	!algorithme de Horner
	lnt9a=(lnt9-LOG(1.+2.69*t9))/3.
	t9a=EXP(lnt9a)
	v0=EXP(LOG(7.79d7)+2.5*lnt9a-lnt932-16.712/t9a)
	v1=EXP(LOG(1.51d8)-lnt923-16.712/t913+LOG(po(0)))
	v2=EXP(LOG(1.56d0)-lnt9-6.272/t9)
	v3=EXP(LOG(3.16d-5)-t932-0.767/t9)
	v4=EXP(LOG(98.d0)-t932-2.077/t9)
	rt(24)=v0+v1+v2+(v3+v4)*zero21l
	IF(rt(24) > 0.d0)THEN
	 rt(24)=LOG(rt(24))
	ELSE
	 rt(24)=ln38
	ENDIF		
c	PRINT*,'r24'
c	WRITE(*,2000)rt(24)

c r�action 25 : O18(p,g)F19	z0=8, z1=1

	ai(0)=1.d0
	ai(1)=0.025d0
	ai(2)=2.26d0
	ai(3)=0.394d0
	ai(4)=30.56d0
	ai(5)=13.55d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(3.45d8)-lnt923-16.729d0/t913-(t9/0.139d0)**2
	1 +LOG(po(0)))
	v1=EXP(LOG(1.25d-15)-lnt932-0.231d0/t9)
	v2=EXP(LOG(1.64d2)-lnt932-1.67d0/t9)
	v3=EXP(LOG(1.28d4)+lnt912-5.098d0/t9)
	rt(25)=v0+v1+v2+v3
	IF(rt(25) > 0.d0)THEN
	 rt(25)=LOG(rt(25))
	ELSE
	 rt(25)=ln38
	ENDIF		
c	PRINT*,'r25'
c	WRITE(*,2000)rt(25)

c r�action 26 : F19(p,a)O16	z0=9, z1=1

	ai(0)=1.d0
	ai(1)=0.023d0
	ai(2)=1.96d0
	ai(3)=0.316d0
	ai(4)=2.86d0
	ai(5)=1.17d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(3.55d11)-lnt923-18.113d0/t913-(t9/0.845d0)**2
	1 +LOG(po(0)))
	v1=EXP(LOG(3.67d6)-lnt932-3.752d0/t9)
	v2=EXP(LOG(3.07d8)-6.019d0/t9)
	v3=LOG(1.d0+EXP(LOG(4.d0)-2.090d0/t9)+EXP(LOG(7.d0)-16.44d0/t9))	
	rt(26)=v0+v1+v2
	IF(rt(26) > 0.d0)THEN
	 rt(26)=LOG(rt(26))-v3
	ELSE
	 rt(26)=ln38
	ENDIF		
c	PRINT*,'r26'
c	WRITE(*,2000)rt(26)
	
c r�action 27 : F19(p,g)Ne20	z0=9, z1=1

	ai(0)=1.d0
	ai(1)=0.023d0
	ai(2)=2.06d0
	ai(3)=0.332d0
	ai(4)=3.16d0
	ai(5)=1.3d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(6.04d7)-lnt923-18.113d0/t913-(t9/0.416d0)**2
	1 +LOG(po(0)))	
	v1=EXP(LOG(6.32d2)-lnt932-3.752d0/t9)
	v2=EXP(LOG(7.56d4)-2.d0/7.d0*lnt9-5.722d0/t9)
	rt(27)=v0+v1+v2	!!!!v3 precedent
	IF(rt(27) > 0.d0)THEN
	 rt(27)=LOG(rt(27))-v3
	ELSE
	 rt(27)=ln38
	ENDIF		
c	PRINT*,'r27'
c	WRITE(*,2000)rt(27)
	
c r�action 28 : O18(p,a)N15	z0=8, z1=1

	ai(0)=1.d0
	ai(1)=0.025d0
	ai(2)=1.88d0
	ai(3)=0.327d0
	ai(4)=4.66d0
	ai(5)=2.06d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(3.63d11)-lnt923-16.729d0/t913-(t9/1.361d0)**2
	1 +LOG(po(0)))
	v1=EXP(LOG(9.9d-14)-lnt932-0.231d0/t9)
	v2=EXP(LOG(2.66d4)-lnt932-1.67d0/t9)
	v3=EXP(LOG(2.41d9)-lnt932-7.638d0/t9)
	v4=EXP(LOG(1.46d9)-lnt9-8.31d0/t9)
	rt(28)=v0+v1+v2+v3+v4
	IF(rt(28) > 0.d0)THEN
	 rt(28)=LOG(rt(28))
	ELSE
	 rt(28)=ln38
	ENDIF		
c	PRINT*,'r28'
c	WRITE(*,2000)rt(28)
	
c r�action 29 : Ne20(a,g)Mg24	z0=10, z1=2

	ai(0)=1.d0
	ai(1)=0.009d0
	ai(2)=0.882d0
	ai(3)=0.055d0
	ai(4)=0.749d0
	ai(5)=0.119d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(MAX(ln38,LOG(4.11d11)-lnt923-46.766d0/t913-(t9/2.219d0)**2+
	1	LOG(po(0))))
	v1=EXP(MAX(ln38,LOG(5.27d3)-lnt932-15.869d0/t9))
	v2=EXP(MAX(ln38,LOG(6.51d3)-lnt912-16.223d0/t9))
	v3=EXP(MAX(ln38,LOG(4.21d1)-lnt932-9.115d0/t9))
	v4=EXP(MAX(ln38,LOG(3.2d1)-lnt923-9.383d0/t9))
	v5=LOG(MAX(ln38,1.d0+5.d0*EXP(-18.96d0/t9)))
	rt(29)=v0+v1+v2+zero21c*(v3+v4)
	IF(rt(29) > 0.d0)THEN
	 rt(29)=LOG(rt(29))-v5
	ELSE
	 rt(29)=ln38
	ENDIF		
c	PRINT*,'r29'
c	WRITE(*,2000)rt(29)

c r�action 30 : C13(a,n)O16	z0=6, z1=2

	ai(0)=1.d0
	ai(1)=0.013d0
	ai(2)=2.04d0
	ai(3)=0.184d0
	CALL polyder(ai,3,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(6.77d15)-lnt923-32.329d0/t913-(t9/1.284d0)**2
	1 +LOG(po(0)))
	v1=EXP(LOG(3.82d5)-lnt932-9.373d0/t9)	
	v2=EXP(LOG(1.41d6)-lnt932-11.873d0/t9)
	v3=EXP(LOG(2.d9)-lnt932-20.409d0/t9)
	v4=EXP(LOG(2.92d9)-lnt932-29.283d0/t9)
	rt(30)=v0+v1+v2+v3+v4	
	IF(rt(30) > 0.d0)THEN
	 rt(30)=LOG(rt(30))
	ELSE
	 rt(30)=ln38
	ENDIF		
c	PRINT*,'r30'
c	WRITE(*,2000)rt(30)
	
c r�action 31 : O17(a,n)Ne20	z0=8, z1=2

	lnt9a=LOG(0.0232d0)+lnt913*5.d0-LOG(1.d0+.0268d0*t9)*2.d0/3.d0
	lnt9a=lnt9-LOG(1.d0+.0268d0*t9+EXP(lnt9a))
	lngt9=-10.106d0/t9-LOG(3.d0)
	lngt9=LOG(1.d0+EXP(lngt9))
	rt(31)=LOG(1.03d18)-lngt9+lnt9a*5.d0/6.d0-lnt932-39.914d0
	1	*EXP(-lnt9a/3.d0)
c	PRINT*,'r31'
c	WRITE(*,2000)rt(31)
	
c r�action 32 : N14(a,g)F18(e+nu)O18	z0=7, z1=2

	ai(0)=1.d0
	ai(1)=0.012d0
	ai(2)=1.45d0
	ai(3)=0.117d0
	ai(4)=1.97d0
	ai(5)=0.406d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(7.78d9)-lnt923-36.031d0/t913-(t9/0.881d0)**2
	1 +LOG(po(0)))	
	v1=EXP(LOG(2.36d-10)-lnt932-2.798d0/t9)
	v2=EXP(LOG(2.03d0)-lnt932-5.054d0/t9)	
	v3=EXP(LOG(1.15d4)-lnt932-12.31d0/t9)	
	rt(32)=v0+v1+v2+v3			
	IF(rt(32) > 0.d0)THEN
	 rt(32)=LOG(rt(32))
	ELSE
	 rt(32)=ln38
	ENDIF		
c	PRINT*,'r32'
c	WRITE(*,2000)rt(32)

c r�action 33 : O18(a,g)Ne22	z0=8, z1=2

	ai(0)=1.d0
	ai(1)=0.01d0
	ai(2)=0.988d0
	ai(3)=0.072d0
	ai(4)=3.17d0
	ai(5)=0.586d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(1.82d12)-lnt923-40.057d0/t913-(t9/0.343d0)**2
	1 +LOG(po(0)))	
	v1=EXP(LOG(7.54d0)-lnt932-6.228d0/t9)
	v2=EXP(LOG(34.8d0)-lnt932-7.301d0/t9)
	v3=EXP(LOG(6.23d3)+lnt9-16.987d0/t9)	
	v4=EXP(LOG(1.d-11)-lnt932-1.994d0/t9)
	rt(33)=v0+v1+v2+v3+zero21c*v4		
	IF(rt(33) > 0.d0)THEN
	 rt(33)=LOG(rt(33))
	ELSE
	 rt(33)=ln38
	ENDIF		
c	PRINT*,'r33'
c	WRITE(*,2000)rt(33)
		
c r�action 34 : Ne22(a,n)Mg25	z0=10, z1=2

	lnt9a=lnt9-LOG(1.d0+.0548d0*t9)
	lngt9=LOG(1.d0+5.d0*EXP(-14.791d0/t9))
	lnft9a=-EXP((LOG(0.197d0)-lnt9a)*4.82d0)
	v0=-lngt9+lnt9a*5.d0/6.d0-lnt932-47.004d0/EXP(lnt9a/3.d0)
	v1=EXP(v0+LOG(4.16d19)+lnft9a)	
	v2=EXP(LOG(1.44d-4)-lngt9-5.577d0/t9)
	rt(34)=v1+v2
	IF(rt(34) > 0.d0)THEN
	 rt(34)=LOG(rt(34))
	ELSE
	 rt(34)=ln38
	ENDIF		
c	PRINT*,'r34'
c	WRITE(*,2000)rt(34)
	
c r�action 35 : Ne22(a,g)Mg26	z0=10, z1=2

	lnfpt9a=-EXP((lnt9a-LOG(0.249d0))*2.31d0)
	v1=EXP(v0+LOG(4.16d19)+lnfpt9a)
	v2=EXP(v0+LOG(2.08d16)+lnft9a)	
	rt(35)=v1+v2
	IF(rt(35) > 0.d0)THEN
	 rt(35)=LOG(rt(35))
	ELSE
	 rt(35)=ln38
	ENDIF		
c	PRINT*,'r35'
c	WRITE(*,2000)rt(35)
	
c r�action 36 : C12(a,n)O15(e+nu)N15 z0=6, z1=2

	ai(0)=1.d0
	ai(1)=0.188d0
	ai(2)=0.015d0
	CALL polyder(ai,5,0,t912,po)	!algorithme de Horner
	rt(36)=2.48d7*po(0)
	IF(rt(36) > 0.d0)THEN
	 rt(36)=LOG(rt(36))
	ELSE
	 rt(36)=ln38
	ENDIF		
c	PRINT*,'r36'
c	WRITE(*,2000)rt(36)

c r�action 37 : Ne21(a,n)Mg24	z0=10, z1=2

	t9a=(t9/(1.d0+0.0537d0*t9))**(1./3.)
	rt(37)=4.94d19*t9a**2.5d0*EXP(-46.89d0/t9a)+2.66d7/t932*
	1 EXP(-22.049d0/t9)
	IF(rt(37) > 0.d0)THEN
	 rt(37)=LOG(rt(37)/(1.d0+1.5d0*EXP(-4.068d0/t9)+2.d0*
	1 EXP(-20.258d0/t9)))
	ELSE
	 rt(37)=ln38
	ENDIF
c	PRINT*,'r37'
c	WRITE(*,2000)rt(37)
	
c r�action 38 : He4(an,g)Be9; la division par 6=3! est incluse

	rt(38)=LOG(2.59d-6)-(2.d0*lnt9+LOG(1.d0+0.344d0*t9))-1.062d0/t9
	1 -LOG(6.d0)
c	PRINT*,'r38'
c	WRITE(*,2000)rt(38)
	
c r�action 39 : Be9(p,d)2He4

	SELECT CASE (total)	
	CASE(0,1)	 

	 ai(0)=1.d0
	 ai(1)=0.04d0
	 ai(2)=1.09d0
	 ai(3)=0.307d0
	 ai(4)=3.21d0
	 ai(5)=2.3d0	
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(2.11d11)-lnt923-10.359d0/t913-(t9/0.52d0)**2
	1 +LOG(po(0)))	
	 v1=EXP(LOG(5.79d8)-lnt9-3.046d0/t9)
	 v2=EXP(LOG(8.5d8)-lnt934-5.8d0/t9)
	 rt(39)=v0+v1+v2		
	 IF(rt(39) > 0.d0)THEN
	  rt(39)=LOG(rt(39))
	 ELSE
	  rt(39)=ln38
	 ENDIF
	 
	CASE(2)			!NACRE formulation	 
	 ai(0)=1.d0
	 ai(1)=-0.427d0
	 ai(2)=34.055d0
	 CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	 v0=EXP(LOG(2.18d11)-lnt923-10.361d0/t913-(t9/0.42d0)**2
	1 +LOG(po(0)))	
	 v1=EXP(LOG(6.24d8)-lnt932-3.446d0/t9)
	 v2=EXP(LOG(3.53d8)-0.205d0*lnt9-3.889d0/t9)
	 rt(39)=v0+v1+v2		
	 IF(rt(39) > 0.d0)THEN
	  rt(39)=LOG(rt(39))
	 ELSE
	  rt(39)=ln38
	 ENDIF
	 
	END SELECT	 
c	PRINT*,'r39'
c	WRITE(*,2000)rt(39)

c r�action 40 : Be9(a,n)C12
	SELECT CASE (total)
	CASE(0,1)	
	 ai(0)=1.d0
	 ai(1)=0.017d0
	 ai(2)=8.57d0
	 ai(3)=1.05d0
	 ai(4)=74.51d0
	 ai(5)=23.15d0	
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(4.62d13)-lnt923-23.87d0/t913-(t9/0.049d0)**2
	1 +LOG(po(0)))
	 v1=EXP(LOG(7.34d-5)-lnt932-1.184d0/t9)
	 v2=EXP(LOG(2.27d-1)-lnt932-1.834d0/t9)
	 v3=EXP(LOG(1.26d5)-lnt932-4.179d0/t9)	
	 v4=EXP(LOG(2.4d8)-12.732d0/t9)	
	 rt(40)=v0+v1+v2+v3+v4		
	 IF(rt(40) > 0.d0)THEN
	  rt(40)=LOG(rt(40))
	 ELSE
	  rt(40)=ln38
	 ENDIF

	CASE(2)		!NACRE formulation	
	 ai(0)=1.d0
	 ai(1)=27.3d0
	 ai(2)=1632.d0
	 CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	 v0=EXP(LOG(5.d13)-lnt923-23.872d0/t913-(t9/0.154d0)**2
	1 +LOG(po(0)))
	 v1=EXP(LOG(0.7d0)-lnt932-1.832d0/t9)
	 v2=EXP(LOG(1.77d5)-lnt932-4.385d0/t9)
	 v3=EXP(LOG(4.12d7)-0.65d0*lnt9-10.06d0/t9)	
	 rt(40)=v0+v1+v2+v3		
	 IF(rt(40) > 0.d0)THEN
	  rt(40)=LOG(rt(40))
	 ELSE
	  rt(40)=ln38
	 ENDIF
	END SELECT		
c	PRINT*,'r40'
c	WRITE(*,2000)rt(40)

c r�action 41 : Li6(p,He3)He4
	ai(0)=1.d0
	ai(1)=0.050d0
	ai(2)=-0.061d0
	ai(3)=-0.021d0
	ai(4)=0.006d0
	ai(5)=0.005d0	
	CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	v0=EXP(LOG(3.73d10)-lnt923-8.413d0/t913-(t9/5.5d0)**2+LOG(po(0)))
	v1=EXP(LOG(1.33d10)-lnt932-17.763d0/t9)
	v2=EXP(LOG(1.29d9)-lnt9-21.82d0/t9)
	rt(41)=v0+v1+v2		
	IF(rt(41) > 0.d0)THEN
	 rt(41)=LOG(rt(41))
	ELSE
	 rt(41)=ln38
	ENDIF		
c	PRINT*,'r41' ; WRITE(*,2000)rt(41)

c r�action 42 : Li6(p,g)Be7
	SELECT CASE (total)
	CASE(0,1)
	 t9a=t9/(1.d0-9.69d-2*t9+2.84d-2*t953/
	1 (1.d0-9.69d-2*t9)**(2.d0/3.d0))	
	 rt(42)=6.69d5*t9a**(5.d0/6.d0)/t932*
	1 EXP(-8.413d0/t9a**(1.d0/3.d0))
	CASE(2)		!NACRE formulation
	 ai(0)=1.d0
	 ai(1)=-0.252d0
	 ai(2)=5.19d-2
	 ai(3)=-2.92d-3
	 CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	 rt(42)=EXP(LOG(1.25d6)-lnt923-8.415d0/t913+LOG(po(0)))	
	END SELECT
	IF(rt(42) > 0.d0)THEN
	 rt(42)=LOG(rt(42))
	ELSE
	 rt(42)=ln38
	ENDIF	 		
c	PRINT*,'r42' ; WRITE(*,2000)rt(42)

c r�action 43 : Be9(p,a)Li6
	SELECT CASE (total)
	CASE(0,1)	
	 ai(0)=1.d0
	 ai(1)=0.040d0
	 ai(2)=1.09d0
	 ai(3)=3.07d0
	 ai(4)=3.21d0
	 ai(5)=2.30d0	
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(2.11d11)-lnt923-10.359d0/t913-(t9/0.520d0)**2+
	1 LOG(po(0)))
	 v1=EXP(LOG(4.51d8)-lnt9-3.046d0/t9)
	 v2=EXP(LOG(6.70d8)-lnt934-5.160d0/t9)
	 rt(43)=v0+v1+v2		
	CASE(2)		!NACRE formulation	
	 ai(0)=1.d0
	 ai(1)=-0.189d0
	 ai(2)=35.2d0
	 CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	 v0=EXP(LOG(2.11d11)-lnt923-10.361d0/t913-(t9/0.4d0)**2+
	1 LOG(po(0)))
	 v1=EXP(LOG(5.24d8)-lnt932-3.446d0/t9)
	 v2=EXP(LOG(4.65d8)-lnt9*0.293d0-4.396d0/t9)
	 rt(43)=v0+v1+v2
	END SELECT
	IF(rt(43) > 0.d0)THEN
	 rt(43)=LOG(rt(43))
	ELSE
	 rt(43)=ln38
	ENDIF	 		
c	PRINT*,'r43' ; WRITE(*,2000)rt(43)

c	r�action 44 : B11(p,a)2He4
	SELECT CASE (total)
	CASE(0,1)	
	 ai(0)=1.d0
	 ai(1)=0.034d0
	 ai(2)=0.140d0
	 ai(3)=0.034d0
	 ai(4)=0.190d0
	 ai(5)=0.116d0	
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(2.20d12)-lnt923-12.095d0/t913-(t9/1.644d0)**2+
	1 LOG(po(0)))
	 v1=EXP(LOG(4.03d6)-lnt932-1.734d0/t9)
	 v2=EXP(LOG(6.73d9)-lnt932-6.262d0/t9)
	 v3=EXP(LOG(3.88d9)-lnt9-14.154d0/t9)	 
	 rt(44)=v0+v1+v2+v3		
	CASE(2)		!NACRE formulation
	 IF(t9 <= 2.d0)THEN	
	  ai(0)=1.d0
	  ai(1)=1.62d0
	  ai(2)=-1.31d0
	  ai(3)=2.60d0
	  CALL polyder(ai,3,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(2.68d12)-lnt923-12.097d0/t913+LOG(po(0)))
	  v1=EXP(LOG(2.12d6)-lnt932-1.724d0/t9)
	  rt(44)=v0+v1
	 ELSE
	  ai(0)=-1.d0
	  ai(1)=0.883d0
	  ai(2)=0.012d0
	  CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	  v0=EXP(LOG(5.84d11)-lnt923-12.097d0/t913+LOG(po(0)))
	  v1=(t923-1.47d0)**2+0.187
	  rt(44)=v0/v1	  	  
	 ENDIF
	END SELECT
	IF(rt(44) > 0.d0)THEN
	 rt(44)=LOG(rt(44))
	ELSE
	 rt(44)=ln38
	ENDIF	 		
c	PRINT*,'r44' ; WRITE(*,2000)rt(44)

c	r�action 45 : B11(p,g)C12
	SELECT CASE (total)
	CASE(0,1)	
	 ai(0)=1.d0
	 ai(1)=0.035d0
	 ai(2)=3.0d0
	 ai(3)=0.723d0
	 ai(4)=9.91d0
	 ai(5)=6.07d0	
	 CALL polyder(ai,5,0,t913,po)	!algorithme de Horner
	 v0=EXP(LOG(4.62d7)-lnt923-12.095d0/t913-(t9/0.239d0)**2+
	1 LOG(po(0)))
	 v1=EXP(LOG(7.89d3)-lnt932-1.733d0/t9)
	 v2=EXP(LOG(9.68d4)-lnt9/5.d0-5.617d0/t9)
	 rt(45)=v0+v1+v2		
	CASE(2)		!NACRE formulation	
	 ai(0)=1.d0
	 ai(1)=0.353d0
	 ai(2)=0.842d0
	 CALL polyder(ai,2,0,t9,po)	!algorithme de Horner
	 v0=EXP(LOG(4.58d7)-lnt923-12.097d0/t913-(t9/0.6d0)**2+
	1 LOG(po(0)))
	 v1=EXP(LOG(2.80d4)+lnt9*0.104d0-3.892d0/t9)
	 rt(45)=v0+v1
	END SELECT
	IF(rt(45) > 0.d0)THEN
	 rt(45)=LOG(rt(45))
	ELSE
	 rt(45)=ln38
	ENDIF	 		
c	PRINT*,'r45' ; WRITE(*,2000)rt(45)

	RETURN

	END SUBROUTINE taux_nuc
