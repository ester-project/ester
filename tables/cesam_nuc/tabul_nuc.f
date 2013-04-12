
!*************************************************************************

	SUBROUTINE tabul_nuc
	
! routine private du module mod_nuc

! tabulation des réactions thermonucléaires (dérivé du programme tab_reac)
! le nom de la routine de réactions nucléaires
! (par Ex. pour ppcno9) détermine les réactions, les éléments à
! utiliser l'intervalle de températures pour la tabulation

! dans un premier appel à taux_nuc on initialise les énergies,
! les masses réduites etc... de toutes les nreac_tot réac. nucléaires
! considérées, puis on ne garde que celles à utiliser

! des appels à taux_reac_nuc calculent les taux, on tabule en f(ln T)

! les variables avec terminaison `t' telles que nucleot, zit etc..
! sont dimensionnées à nreac_tot >= nreac car toutes les réactions
! ne sont pas toujours toutes utilisées, idem pour les éléments chimiques

! voir dans taux_nuc la procédure pour ajouter des éléments, des
! isotopes, des réactions

! Auteur: P.Morel, Département J.D. Cassini, O.C.A.
! CESAM2k

!--------------------------------------------------------------------

	USE mod_donnees, ONLY : amu, eve, ihe4, langue, nchim,
	1 nom_elem, nom_nuc, nom_nuc_cpl, nucleo, t_inf, t_stop, zi
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn
	
	IMPLICIT NONE
	
	REAL (kind=dp), DIMENSION(nreac_tot) :: rt, qt, at
	REAL (kind=dp), DIMENSION(-1:niso_tot) :: nucleot, zit	
	REAL (kind=dp) :: ti, pas

	INTEGER, DIMENSION(nreac_tot,2) :: izzt
	INTEGER, DIMENSION(nreac_tot) :: nb
	INTEGER, DIMENSION(-1:niso_tot) :: ind
	INTEGER :: i, j, nchim_reac, total
	
	LOGICAL :: lisse
					
	CHARACTER (len=4), DIMENSION(-1:niso_tot) :: nom_elemt
	CHARACTER (len=20), DIMENSION(nreac_tot) ::  nom_react

!-----------------------------------------------------------------	 

2000	FORMAT(8es10.3)
2002	FORMAT(20i4)

! nchim_reac est le nombre d'isotopes utilisés dans les réactions
! nchim=nchim_reac + Ex, nchim est défini dans la routine de
! réactions nucléaires, Ex: ppcno9
! pour éviter de reclasser les isotopes et de mettre Ex comme dernier
! élément on a aussi nchim_reac=nchim, auquel cas l'indice de Ex est -1
! exemple : ppcno10BeBFe
! PRINT*,nchim ; PAUSE'entrée tabul'

	SELECT CASE (langue)
	CASE('english')
	 WRITE(2,1001)nom_nuc ; WRITE(*,1001)nom_nuc
1001	 FORMAT('Thermonuclear reactions tabulated for : ',a20)	
	CASE DEFAULT
	 WRITE(2,1)nom_nuc ; WRITE(*,1)nom_nuc
1	 FORMAT('Réactions nucléaires tabulées pour: ',a20)
	END SELECT 
	 
	SELECT CASE(nom_nuc)
	CASE ('pp1')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1002) ; WRITE(*,1002) 	
1002	  FORMAT('No tabulation for the thermonuclear reaction PP1',/)
	 CASE DEFAULT		
	  WRITE(2,2) ; WRITE(*,2) 	
2	  FORMAT('Pas de tabulation de la réaction PP1 simplifié',/)
	 END SELECT
	 
	CASE ('ppcno9','ppcno9Fe')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1003) ; WRITE(*,1003) 
1003	  FORMAT('PP+CNO, 9 species, H2, Li7, Be7 at equilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,3) ; WRITE(*,3) 
3	  FORMAT('PP+CNO, 9 éléments, H2, Li7, Be7 à l''équilibre',/)
	 END SELECT
	 	  
! nombre de réactions et correspondance des indices
	 nreac=14	!nombre de réactions	 
	 nchim_reac=9
	 ihe4=3	!indice de He4

! indices, dans taux_nuc, des isotopes utilisés dans les réactions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
 	 
! ordre des réactions: on garde l'ordre de reac_c 
 	 nb(1:nreac)=(/ (i, i=1,nreac) /)
 	 
! températures max/min de tabulation
	 t_inf=1.d6 ; t_sup=80.d6 ; lisse=.FALSE.
	 
	CASE ('pp3')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1021) ; WRITE(*,1021) 
1021	  FORMAT('PP, 3 species, H2, Li7, Be7 at equilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,21) ; WRITE(*,21) 
21	  FORMAT('PP, 3 éléments, H2, Li7, Be7 à l''équilibre',/)
	 END SELECT
	 	  
! nombre de réactions et correspondance des indices
	 nreac=7	!nombre de réactions	 
	 nchim_reac=3
	 ihe4=3	!indice de He4

! indices, dans taux_nuc, des isotopes utilisés dans les réactions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 
! ordre des réactions: on garde l'ordre de reac_c 
 	 nb(1:nreac)=(/ (i, i=1,nreac) /)
 	 
! températures max/min de tabulation
	 t_inf=1.0d6 ; t_sup=40.d6 ; lisse=.FALSE.
	  
	CASE('ppcno10','ppcno10Fe','ppcno10K')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1004) ; WRITE(*,1004) 	
1004	  FORMAT('PP+CNO, 10 species, H2, Be7 at equilibrium',/)
	 CASE DEFAULT		
	  WRITE(2,4) ; WRITE(*,4) 	
4	  FORMAT('PP+CNO, 10 éléments, H2, Be7 à l''équilibre',/)
	 END SELECT
	 	 
! nombre de réactions et correspondance des indices 
	 nreac=14	!nombre de réactions
	 nchim_reac=10	 
	 ihe4=3		!indice de He4
	 	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=3	!He3 
 	 ind(3)=4	!He4
	 ind(4)=5	!Li7 
	 ind(5)=7	!C12
	 ind(6)=8	!C13	 
	 ind(7)=10	!N14	 
 	 ind(8)=11	!N15
 	 ind(9)=12	!O16
 	 ind(10)=13	!O17

! ordre des réactions: on garde l'ordre de reac_c
	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.
	 
	CASE('ppcno10BeBFe')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1016) ; WRITE(*,1016) 	
1016	  FORMAT('PP + CNO + Be + B + Fe with H2 and Be7 at equilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,16) ; WRITE(*,16) 	
16	  FORMAT('PP + CNO + Be + B + Fe avec H2, Be7 à l''équilibre',/)
	 END SELECT
	 
! nombre de réactions et correspondance des indices 
	 nreac=20	!nombre de réactions
	 nchim_reac=15	 
	 ihe4=3		!indice de He4
	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=3	!He3 
 	 ind(3)=4	!He4
	 ind(4)=5	!Li7
	 ind(5)=7	!C12
	 ind(6)=8	!C13	 
	 ind(7)=10	!N14	 
 	 ind(8)=11	!N15
 	 ind(9)=12	!O16
 	 ind(10)=13  	!O17
 	 ind(11)=27	!Fe56
 	 ind(12)=-1	!Ex
 	 ind(13)=26	!Li6	 	 	 
 	 ind(14)=25	!Be9	 
 	 ind(15)=28	!B11
	  
! ordre des réactions
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=39	!Be9(p,d)2He4
	 nb(16)=41	!Li6(p,He3)He4
	 nb(17)=42	!Li6(p,g)Be7
	 nb(18)=43	!Be9(p,a)Li6
 	 nb(19)=44	!B11(p,a)2He4
	 nb(20)=45	!B11(p,g)C12
  	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.

	CASE('ppcno11')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1005) ; WRITE(*,1005) 	
1005	  FORMAT('PP+CNO, 11 species, Be7 at equilibrium',/)	
	 CASE DEFAULT	
	  WRITE(2,5) ; WRITE(*,5) 	
5	  FORMAT('PP+CNO, 11 éléments, Be7 à l''équilibre',/)	
	 END SELECT
	 
! nombre de réactions et correspondance des indices 
	 nreac=14	!nombre de réactions
	 nchim_reac=11	 
	 ihe4=4		!indice de He4
	 	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7 
	 ind(6)=7	!C12
	 ind(7)=8	!C13	 
	 ind(8)=10	!N14	 
 	 ind(9)=11	!N15
 	 ind(10)=12	!O16
 	 ind(11)=13	!O17
 
! ordre des réactions: on garde l'ordre de reac_c
	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.	 
		  
	CASE('ppcno12')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1006) ; WRITE(*,1006) 	
1006	  FORMAT('PP+CNO, 12 species',/)	
	 CASE DEFAULT	
	  WRITE(2,6) ; WRITE(*,6) 	
6	  FORMAT('PP+CNO, 12 éléments',/)
	 END SELECT
	 
! nombre de réactions et correspondance des indices 
	 nreac=14	!nombre de réactions
	 nchim_reac=12	 
	 ihe4=4		!indice de He4
	 	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 
! ordre des réactions: on garde l'ordre de reac_c
	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.	 
	
	CASE('ppcno12Be')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1007) ; WRITE(*,1007) 	
1007	  FORMAT('PP+CNO+Be9, 12 species+Be9',/)
	 CASE DEFAULT	
	  WRITE(2,7) ; WRITE(*,7) 	
7	  FORMAT('PP+CNO+Be9, 12 éléments+Be9',/)
	 END SELECT
	 
! nombre de réactions et correspondance des indices 
	 nreac=16	!nombre de réactions
	 nchim_reac=13	 
	 ihe4=4		!indice de He4
 	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 	 ind(13)=25	!Be9
 
! ordre des réactions: on garde l'ordre de reac_c
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=39
	 nb(16)=40
 	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.
	 
	CASE('ppcno12Li')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1008) ; WRITE(*,1008) 	
1008	  FORMAT('PP+CNO+Li6, 12 species',/)	
	 CASE DEFAULT	
	  WRITE(2,8) ; WRITE(*,8) 	
8	  FORMAT('PP+CNO+Li6, 12 éléments',/)
	 END SELECT
	 
! nombre de réactions et correspondance des indices
	 nreac=15	!nombre de réactions
	 nchim_reac=13	 
	 ihe4=4		!indice de He4
	 	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 	 ind(13)=26	!Li6
 
! ordre des réactions: on garde l'ordre de reac_c
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=41	!réaction 41 Li6(p,He3)He4
 	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.	 

	CASE('ppcno12BeBFe')
	 SELECT CASE (langue)
	 CASE('english')
	 WRITE(2,1013) ; WRITE(*,1013) 	
1013	 FORMAT('PP + CNO + Be + B + Fe, none species at equilibrium',/)	
	 CASE DEFAULT	
	  WRITE(2,13) ; WRITE(*,13) 	
13	  FORMAT('PP + CNO + Be + B + Fe, aucun élément à l''équilibre',/)
	 END SELECT
	 
! nombre de réactions et correspondance des indices 
	 nreac=20	!nombre de réactions
	 nchim_reac=17	 
	 ihe4=4		!indice de He4
	 	 
! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 	 ind(13)=25	!Be9
 	 ind(14)=-1	!Ex	 
 	 ind(15)=28	!B11
 	 ind(16)=27	!Fe56	 
 	 ind(17)=26	!Li6
	 
! ordre des réactions
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=39	!Be9(p,d)2He4
	 nb(16)=41	!Li6(p,He3)He4
	 nb(17)=42	!Li6(p,g)Be7
	 nb(18)=43	!Be9(p,a)Li6
 	 nb(19)=44	!B11(p,a)2He4
	 nb(20)=45	!B11(p,g)C12
  	 
! températures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; lisse=.FALSE.
	 
	CASE('ppcno3a9')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1009) ; WRITE(*,1009) 	
1009	  FORMAT('PP+CNO+3alpha, 9 species, H2, Li7, Be7 at quilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,9) ; WRITE(*,9) 	
9	  FORMAT('PP+CNO+3al, 9 éléments, H2, Li7, Be7 à l''équilibre',/)
	 END SELECT
	 
! nombre de réactions et correspondance des indices
	 nreac=17	!nombre de réactions
	 nchim_reac=9	 
	 ihe4=3		!indice de He4
	 i3al=15	!indice de la réaction 3 alpha

! indices, dans taux_nuc, des isotopes utilisés dans les réactions
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
 	 
! ordre des réactions: on garde l'ordre de reac_c 
 	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
! températures max/min de tabulation
	 t_inf=1.0d6 ; t_sup=500.d6 ; lisse=.TRUE.	 
	  	 
	CASE('ppcno3a12Ne') 
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1022) ; WRITE(*,1022) 	
1022	  FORMAT('PP+CNO+3al+C, 12 species, H2, Li7, Be7 at equilibrium')
	 CASE DEFAULT	
	  WRITE(2,22) ; WRITE(*,22) 	
22	  FORMAT('PP+CNO+3al+C, 12 éléments, H2, Li7, Be7 en équilibre')
	 END SELECT
	 
! nombre de réactions et correspondance des indices
	 nreac=25	!nombre de réactions
	 nchim_reac=12	 
	 ihe4=3		!indice de He4
	 i3al=15	!indice de la réaction 3 alpha

! indices, dans taux_nuc, des isotopes utilisés dans les réactions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
 	 ind(10)=16		!Ne20
	 ind(11)=-1		!Ex
	 ind(12)=14		!O18
	 
! ordre des réactions: on garde l'ordre de reac_c
 	 nb(1:17)=(/ (i, i=1,17) /)
	 nb(18)=32 ; nb(19)=33 ; nb(20)=24 ; nb(21)=28 ; nb(22)=29
	 nb(23)=18 ; nb(24)=20 ; nb(25)=21
 	 
! températures max/min de tabulation
	 t_inf=1.0d6 ; t_sup=0.8d9 ; lisse=.TRUE.
	 	 
	CASE('ppcno3aco') 
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1027) ; WRITE(*,1027) 	
1027	  FORMAT('PP+CNO+3al+C+0, 17 species, H2, Li7, Be7 at equilibrium')
	 CASE DEFAULT	
	  WRITE(2,27) ; WRITE(*,27) 	
27	  FORMAT('PP+CNO+3al+C+O, 17 éléments, H2, Li7, Be7 en équilibre')
	 END SELECT

! nombre de réactions et correspondance des indices
	 nreac=33	!nombre de réactions
	 nchim_reac=17	 
	 ihe4=3		!indice de He4
	 i3al=15	!indice de la réaction 3alpha

! indices, dans taux_nuc, des isotopes utilisés dans les réactions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
	 ind(10)=-1		!Ex	 
 	 ind(11)=16		!Ne20
	 ind(12)=23		!Na23
	 ind(13)=18		!Mg24
	 ind(14)=30		!Al27	 	 
	 ind(15)=29		!Si28
	 ind(16)=32		!P31	 
	 ind(17)=31		!S32
	 	 
! correspondance de l'ordre des réactions ppcno3acos <---> taux_nuc
 	 nb(1:17)=(/ (i, i=1,17) /)
	 nb(18)=29 ; nb(19)=18 ; nb(20)=20 ; nb(21)=21 ; nb(22)=56
	 nb(23)=57 ; nb(24)=48 ; nb(25)=50 ; nb(26)=51 ; nb(27)=52
	 nb(28)=54 ; nb(29)=55 ; nb(30)=58 ; nb(31)=22 ; nb(32)=47
 	 nb(33)=64 
	 
! températures max/min de tabulation
	 t_inf=1.0d6 ; t_sup=3.0d9 ; lisse=.TRUE.
	 	 	 	  
	CASE DEFAULT
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1017)nom_nuc ; WRITE(*,1017)nom_nuc
1017	  FORMAT('STOP, unknown routine of thermonuclear reactions : ',a,/,
	1 'known routines :')
	 CASE DEFAULT
	  WRITE(2,17)nom_nuc ; WRITE(*,17)nom_nuc
17	  FORMAT('ARRET, routine de réactions nucléaires inconnue : ',a,/,
	1 'routines connues :')
	 END SELECT
	 WRITE(2,18) ; WRITE(*,18) ; STOP	
18	 FORMAT('pp1, pp3, ppcno9, ppcno10, ppcno10Fe, ppcno10K, ppcno11',/,
	1 'ppcno12, ppcno12Be, ppcno12Li, ppcno3a9, ppcno3a12Ne, ppcno3aco')
	END SELECT
	
! 200 points de tabulation par 100d6 K, sauf pour PP1
	IF(TRIM(nom_nuc) /= 'pp1')THEN	
	 n_temp=NINT(t_sup*1.d-6)*2
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1024)n_temp,t_inf,t_sup ; WRITE(*,1024)n_temp,t_inf,t_sup
1024	  FORMAT('Number of tabulated  points : ',i5,' from',es10.3,'K to',
	1 es10.3,'K')
	 CASE DEFAULT
	  WRITE(2,24)n_temp,t_inf,t_sup ; WRITE(*,24)n_temp,t_inf,t_sup
24	  FORMAT('Nombre de points de tabulation : ',i5,' de',es10.3,'K à',
	1 es10.3,'K')
	 END SELECT
	ENDIF
	
! on limite t_stop à 0.95 t_sup pour éviter des dépassements de tabulation
! t_stop est lu dans mon_modele.don
	t_stop=MIN(0.95d0*t_sup,t_stop)

! sélection du type de compilation à utiliser	
	SELECT CASE(nom_nuc_cpl)
	CASE('Cau-Fow')
	 total=0
	CASE('Adelb')
	 total=1	
	CASE('NACRE')
	 total=2	
	CASE DEFAULT
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1019)nom_nuc_cpl ; WRITE(*,1019)nom_nuc_cpl	
1019	  FORMAT('STOP, unknown compilation of thermonuclear reactions : ',
	1 a,/,'known compilations :')
	 CASE DEFAULT	
	  WRITE(2,19)nom_nuc_cpl ; WRITE(*,19)nom_nuc_cpl	
19	  FORMAT('ARRET, compilation de réactions nucléaires inconnue : ',
	1 a,/,'compilations connues :')
	END SELECT	 
	WRITE(2,20) ; WRITE(*,20) ; STOP
20	FORMAT('Cau-Fow, Adelb, NACRE')	
	END SELECT

! initialisations (ti=20.d6) des énergies: q0, masses réduites: ar,
! masses des noyaux: nucleo, charges des noyaux: zi
! charges des noyaux Z1, Z2 des noyaux de la réaction : izz
! noms des réactions: nom_react 	
! appel pour remplir le common/reac_nuc/
! identification des q, nom_react, z1, z2
	ALLOCATE(q0(nreac),ar(nreac),izz(nreac,2),
	1 nom_elem(nchim),nucleo(nchim),zi(nchim))
	
	ti=20.d6 ; CALL taux_nuc(ti,total,rt,zit,izzt,qt,nom_react,nucleot,
	1 at,nom_elemt)	
	DO i=1,nreac
	 q0(i)=qt(nb(i))*eve*1.d6/amu	!énergie des réac. en erg/réaction
	 ar(i)=at(nb(i)) ; izz(i,:)=izzt(nb(i),:)
	ENDDO
	
	DO i=1,nchim_reac	 
	 nom_elem(i)=nom_elemt(ind(i)) ; nucleo(i)=nucleot(ind(i))
	 zi(i)=zit(ind(i))
	ENDDO

! tabulation des points en ln(température)	
	ALLOCATE(temp(n_temp))
	pas=(t_sup-t_inf)/REAL(n_temp-1,dp)
	DO i=1,n_temp
	 temp(i)=t_inf+pas*REAL(i-1,dp)	 
!	 WRITE(*,2000)temp(i),EXP(t(i))
	ENDDO

! tabulation des réactions
	ALLOCATE(taux_reac(nreac,n_temp))
	DO i=1,n_temp
	 CALL taux_nuc(temp(i),total,rt,zit,izzt,qt,nom_react,
	1 nucleot,at,nom_elemt)
	 DO j=1,nreac
	  taux_reac(j,i)=rt(nb(j))
	 ENDDO	!j
!	 WRITE(*,2000)temp(i),EXP(temp(i)),rt(1) ; WRITE(*,2000)rt(2:9)
!	 WRITE(*,2000)rt(10:17)
	ENDDO	!i
	
	ALLOCATE(ttemp(n_temp+m_temp))
	temp=LOG(temp)
	CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .FALSE.,temp(1),j,rt,qt,lisse)		!rt, qt: Vt
	
! écritures diverses	 
	SELECT CASE (langue)
	CASE('english')
	 WRITE(2,1011) ; WRITE(*,1011)
1011	 FORMAT(/,'isotopes employed :')
	CASE DEFAULT
	 WRITE(2,11) ; WRITE(*,11)
11	 FORMAT(/,'isotopes utilisés dans les réactions :')
	END SELECT	
	DO i=1,nchim_reac
	 WRITE(*,15)i,nom_elem(i) ; WRITE(2,15)i,nom_elem(i)
15	 FORMAT(i3,' : ',a)
	ENDDO

	SELECT CASE (langue)
	CASE('english')	
	 WRITE(2,1012) ; WRITE(*,1012)
1012	 FORMAT(/,'thermonuclear reactions at work :')
	CASE DEFAULT	  
	 WRITE(2,12) ; WRITE(*,12)
12	 FORMAT(/,'réactions nucléaires utilisées :')
	END SELECT	
	DO i=1,nreac
	 WRITE(*,15)i,nom_react(nb(i)) ; WRITE(2,15)i,nom_react(nb(i))
	ENDDO	
		
	RETURN
	 
	END SUBROUTINE tabul_nuc
