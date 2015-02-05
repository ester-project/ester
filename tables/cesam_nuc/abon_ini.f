
!******************************************************************

	SUBROUTINE abon_ini

! routine private du module mod_nuc
	
! initialisation des abondances, des rapports isotopiques,
!	des charges des éléments. Les abondances initiales des isotopes
!	utilisés dans le réseau nucléaire seront déduites de ces valeurs
!	initiales suivant la chaine de réactions utilisées et les isotopes
!	seulement diffusés.
!	source : ftp:://www-phys.llnl.gov/pub/opal

! pour ajouter un élément:
!	dans le module mod_nuc augmenter le PARAMETER nelem_ini d'une
!	unité, Ex: nelem_ini=28

!	Auteur: P.Morel, Département Cassiopée, O.C.A.,CESAM2k
!	Conseils de F.Thévenin, Département Cassiopée, O.C.A. 
!	contributions d'Yveline Lebreton

!-----------------------------------------------------------------------

	USE mod_donnees, ONLY : abon_m, garde_xish, langue, li_ini, ln_Tli,
	1 modif_chim, nom_abon, nom_fich2, x0, y0, z0, zsx_sol
	USE mod_kind
	USE mod_variables, ONLY : sortie

	IMPLICIT NONE
	
	REAL(kind=dp), DIMENSION(28) :: abon_relam, add, nisnh
	REAL(kind=dp) :: add_Li, add_Be, add_B, add_C,
	1 add_N, add_O, add_F, add_Ne, add_Na, add_Mg, add_Al, add_Si,
	2 add_P, add_S, add_Cl, add_Ar, add_K, add_Ca, add_Sc, add_Ti,
	3 add_V, add_Cr, add_Mn, add_Fe, add_Co, add_Ni, add_Z, fesh,
	4 sniai, zsx

	INTEGER :: i
	
	LOGICAL :: ok
	
	CHARACTER (len=256) :: chain

	
	NAMELIST/nl_mixture/ab
	NAMELIST/nl_rap_iso/be7sbe9, be7sz, c13sc12, h2sh1,
	1 he3she4, he3she4z, li6sli7, mg25smg24, mg26smg25, ne22sne20,
	2 n15sn14, o17so16, o18so16
	NAMELIST/nl_modif_mix/ add_Li, add_Be, add_B, add_C,
	1 add_N, add_O, add_F, add_Ne, add_Na, add_Mg, add_Al, add_Si,
	2 add_P, add_S, add_Cl, add_Ar, add_K, add_Ca, add_Sc, add_Ti,
	3 add_V,add_Cr, add_Mn, add_Fe, add_Co, add_Ni, add_Z
      
!-------------------------------------------------------------------

2000	FORMAT(8es10.3)

!	ln_Tli : température au dessus de laquelle il n'y a pas de lithium
!	inital pour le modèle de ZAMS Li initial

	ln_Tli=LOG(3.d6) ; li_ini=1.d-15

! charges et masses (OPAL) des éléments
	c=(/ (i, i=1,nelem_ini) /)
	elem(1)=' H'  ; m(1)=1.0079d0
	elem(2)='He'  ; m(2)=4.0026d0
	elem(3)='Li'  ; m(3)=6.941d0    
	elem(4)='Be'  ; m(4)=9.0122d0
	elem(5)=' B'  ; m(5)=10.811d0       
	elem(6)=' C'  ; m(6)=12.011d0
	elem(7)=' N'  ; m(7)=14.0067d0  
	elem(8)=' O'  ; m(8)=15.9994d0  
	elem(9)=' F'  ; m(9)=18.9984d0  
	elem(10)='Ne' ; m(10)=20.179d0  
	elem(11)='Na' ; m(11)=22.98977d0
	elem(12)='Mg' ; m(12)=24.305d0  
	elem(13)='Al' ; m(13)=26.98154d0 
	elem(14)='Si' ; m(14)=28.0855d0  
	elem(15)=' P' ; m(15)=30.97376d0
	elem(16)=' S' ; m(16)=32.06d0   
	elem(17)='Cl' ; m(17)=35.453d0  
	elem(18)='Ar' ; m(18)=39.948d0
	elem(19)=' K' ; m(19)=39.0983d0          
	elem(20)='Ca' ; m(20)=40.08d0   
	elem(21)='Sc' ; m(21)=44.956d0  
	elem(22)='Ti' ; m(22)=47.90d0
	elem(23)=' V' ; m(23)=50.9414d0
	elem(24)='Cr' ; m(24)=51.996d0
	elem(25)='Mn' ; m(25)=54.938d0
	elem(26)='Fe' ; m(26)=55.847d0
	elem(27)='Co' ; m(27)=58.9332d0
	elem(28)='Ni' ; m(28)=58.7d0

! abondances initiales
	SELECT CASE(nom_abon)

	CASE('meteorites_ag')	!Z/X=2.672E-02, [Fe/H]= 7.509-12, Z= 1.890E-02
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1001) ; WRITE(2,1001)
1001	  FORMAT('Meteoritic mixture of Anders & Grevesse 89')	
	 CASE DEFAULT
	  WRITE(*,1) ; WRITE(2,1)
1	  FORMAT('Mixture des météorites de Anders & Grevesse 89')
	 END SELECT	!la normalisation à H=1.d12 effectuée plus loin
	 ab(1)=2.79d10 ; ab(2)=2.72d9   ; ab(3)=57.1d0  ; ab(4)=0.73d0   ; ab(5)=21.2d0
	 ab(6)=1.01d7  ; ab(7)=3.13d6   ; ab(8)=2.38d7  ; ab(9)=843.d0   ; ab(10)=3.44d6
	 ab(11)=5.74d4 ; ab(12)=1.074d6 ; ab(13)=8.49d4 ; ab(14)=1.d6    ; ab(15)=1.04d4
	 ab(16)=5.15d5 ; ab(17)=5240.d0 ; ab(18)=1.01d5 ; ab(19)=3370.d0 ; ab(20)=6.11d4
	 ab(21)=34.2d0 ; ab(22)=2400.d0 ; ab(23)=293.d0 ; ab(24)=1.35d4  ; ab(25)=9550.d0
	 ab(26)=9.d5   ; ab(27)=2250.d0 ; ab(28)=4.93d4	 

	CASE('meteorites_gs')	!Z/X=2.293E-02, [Fe/H]=7.5-12, Z= 1.685E-02	
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1026) ; WRITE(2,1026)
1026	  FORMAT('Meteoritic mixture of Grevesse & Sauval 98')	
	 CASE DEFAULT
	  WRITE(*,26) ; WRITE(2,26)
26	  FORMAT('Mixture des météorites de Grevesse & Sauval 98')
	 END SELECT
	 ab(1)=12.d0   ; ab(2)=10.93d0 ; ab(3)=3.31d0  ; ab(4)=1.42d0  ; ab(5)=2.79d0
	 ab(6)=8.52d0  ; ab(7)=7.920d0 ; ab(8)=8.83d0  ; ab(9)=4.48d0  ; ab(10)=8.08d0
	 ab(11)=6.32d0 ; ab(12)=7.58d0 ; ab(13)=6.49d0 ; ab(14)=7.56d0 ; ab(15)=5.56d0
	 ab(16)=7.20d0 ; ab(17)=5.28d0 ; ab(18)=6.40d0 ; ab(19)=5.13d0 ; ab(20)=6.35d0
	 ab(21)=3.10d0 ; ab(22)=4.94d0 ; ab(23)=4.02d0 ; ab(24)=5.69d0 ; ab(25)=5.53d0
	 ab(26)=7.50d0 ; ab(27)=4.91d0 ; ab(28)=6.25d0
	 ab=10.d0**ab	!car données en LOG10	 
	 	 
	CASE('meteorites_ags05')	!Z/X=46.72, [Fe/H]=11.20-12, Z=9.790E-01
	 SELECT CASE(langue)		!X= 2.095E-02, Y= 9.125E-09
	 CASE('english')	
	  WRITE(*,1030) ; WRITE(2,1030)	!table d'Yveline Lebreton
1030	  FORMAT('Meteoritic mixture of Asplund Grevesse & Sauval 05')	
	 CASE DEFAULT			
	  WRITE(*,30) ; WRITE(2,30)
30	  FORMAT('Mixture des météorites de Asplund Grevesse & Sauval 05')
	 END SELECT
	 ab(1)=8.25d0 ; ab(2)=1.29d0 ; ab(3)=3.25d0 ; ab(4)=1.38d0 ; ab(5)=2.75d0
	 ab(6)=7.40d0 ; ab(7)=6.25d0 ; ab(8)=8.39d0 ; ab(9)=4.43d0 ; ab(10)=-1.06d0
	 ab(11)=6.27d0 ; ab(12)=7.53d0 ; ab(13)=6.43d0 ; ab(14)=7.51d0 ; ab(15)=5.40d0
	 ab(16)=7.16d0 ; ab(17)=5.23d0 ; ab(18)=-0.45d0 ; ab(19)=5.06d0 ; ab(20)=6.29d0
	 ab(21)=3.04d0 ; ab(22)=4.89d0 ; ab(23)=3.97d0 ; ab(24)=5.63d0 ; ab(25)=5.47d0
	 ab(26)=7.45d0 ; ab(27)=4.86d0 ; ab(28)=6.19d0
	 ab=ab-ab(1)+12.d0	!normalisation à H=1.d12
	 ab=10.d0**ab		!car données en LOG10	 	 
	 
	CASE('solaire_gn')	!Z/X=2.440E-02, [Fe/H]= 7.5-12, Z= 1.727E-02	
	 SELECT CASE(langue)
	 CASE('english')	
c	  WRITE(*,1002) ; WRITE(2,1002)
1002	  FORMAT('Solar mixture of Grevesse & Noels 93')	
	 CASE DEFAULT
c	  WRITE(*,2) ; WRITE(2,2)
2	  FORMAT('Mixture solaire de Grevesse & Noels 93')
	 END SELECT
	 ab(1)=12.00d0 ; ab(2)=10.99d0 ; ab(3)=1.16d0  ; ab(4)=1.15d0  ; ab(5)=2.60d0
	 ab(6)=8.55d0  ; ab(7)=7.97d0  ; ab(8)=8.87d0  ; ab(9)=4.56d0  ; ab(10)=8.08d0
	 ab(11)=6.33d0 ; ab(12)=7.58d0 ; ab(13)=6.47d0 ; ab(14)=7.55d0 ; ab(15)=5.45d0
	 ab(16)=7.21d0 ; ab(17)=5.5d0  ; ab(18)=6.52d0 ; ab(19)=5.12d0 ; ab(20)=6.36d0
	 ab(21)=3.17d0 ; ab(22)=5.02d0 ; ab(23)=4.d0   ; ab(24)=5.67d0 ; ab(25)=5.39d0
	 ab(26)=7.5d0  ; ab(27)=4.92d0 ; ab(28)=6.25d0
	 ab=10.d0**ab	!car données en LOG10

	CASE('solaire_gs')	!Z/X=2.307E-02, [Fe/H]=7.5-12, Z= 1.695E-02	
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1017) ; WRITE(2,1017)
1017	  FORMAT('Solar mixture of Grevesse & Sauval 98')	
	 CASE DEFAULT			
	  WRITE(*,17) ; WRITE(2,17)
17	  FORMAT('Mixture solaire de Grevesse & Sauval 98')
	 END SELECT
	 ab(1)=12.d0   ; ab(2)=10.93d0 ; ab(3)=1.10d0  ; ab(4)=1.40d0  ; ab(5)=2.55d0
	 ab(6)=8.52d0  ; ab(7)=7.920d0 ; ab(8)=8.83d0  ; ab(9)=4.56d0  ; ab(10)=8.08d0
	 ab(11)=6.33d0 ; ab(12)=7.58d0 ; ab(13)=6.47d0 ; ab(14)=7.55d0 ; ab(15)=5.45d0
	 ab(16)=7.33d0 ; ab(17)=5.50d0 ; ab(18)=6.40d0 ; ab(19)=5.12d0 ; ab(20)=6.36d0
	 ab(21)=3.17d0 ; ab(22)=5.02d0 ; ab(23)=4.00d0 ; ab(24)=5.67d0 ; ab(25)=5.39d0
	 ab(26)=7.50d0 ; ab(27)=4.92d0 ; ab(28)=6.25d0
	 ab=10.d0**ab	!car données en LOG10
	 	 
	CASE('solaire_ags03')	!Z/X=0.0171, [Fe/H]=7.45-12, Z= 1.262E-02	
	 SELECT CASE(langue)	!Asplund 2003, ASP Conf. Series vol 304.
	 CASE('english')	!données communiquées par Yveline
	  WRITE(*,1029) ; WRITE(2,1029)
1029	  FORMAT('Solar mixture of Asplund Grevesse & Sauval 03')	
	 CASE DEFAULT	
	  WRITE(*,29) ; WRITE(2,29)
29	  FORMAT('Mixture solaire de Asplund Grevesse & Sauval 03')
	 END SELECT
	 ab(1)=12.00d0 ; ab(2)=10.93d0 ; ab(3)=1.10d0  ; ab(4)=1.40d0  ; ab(5)=2.55d0
	 ab(6)=8.41d0  ; ab(7)=7.80d0  ; ab(8)=8.66d0  ; ab(9)=4.56d0  ; ab(10)=7.84d0
	 ab(11)=6.33d0 ; ab(12)=7.58d0 ; ab(13)=6.47d0 ; ab(14)=7.51d0 ; ab(15)=5.45d0
	 ab(16)=7.33d0 ; ab(17)=5.5d0  ; ab(18)=6.18d0 ; ab(19)=5.12d0 ; ab(20)=6.36d0
	 ab(21)=3.17d0 ; ab(22)=5.02d0 ; ab(23)=4.d0   ; ab(24)=5.67d0 ; ab(25)=5.39d0
	 ab(26)=7.45d0  ; ab(27)=4.92d0 ; ab(28)=6.25d0
	 ab=10.d0**ab	!car données en LOG10
	 
	CASE('solaire_ags05')	!Z/X=0.0166, [Fe/H]=7.45-12, Z= 1.220E-02	
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1031) ; WRITE(2,1031)	!table de ASP conf. 336, 25, 2005
1031	  FORMAT('Solar mixture of Asplund Grevesse & Sauval 05')	
	 CASE DEFAULT			
	  WRITE(*,31) ; WRITE(2,31)
31	  FORMAT('Mixture solaire de Asplund Grevesse & Sauval 05')
	 END SELECT
	 ab(1)=12.00d0 ; ab(2)=10.93d0 ; ab(3)=1.05d0 ; ab(4)=1.38d0 ; ab(5)=2.70d0
	 ab(6)=8.39d0  ; ab(7)=7.78d0  ; ab(8)=8.66d0  ; ab(9)=4.56d0  ; ab(10)=7.84d0
	 ab(11)=6.17d0 ; ab(12)=7.53d0 ; ab(13)=6.37d0 ; ab(14)=7.51d0 ; ab(15)=5.36d0
	 ab(16)=7.14d0 ; ab(17)=5.50d0 ; ab(18)=6.18d0 ; ab(19)=5.08d0 ; ab(20)=6.31d0
	 ab(21)=3.05d0 ; ab(22)=4.90d0 ; ab(23)=4.00d0 ; ab(24)=5.64d0 ; ab(25)=5.39d0
	 ab(26)=7.45d0  ; ab(27)=4.92d0 ; ab(28)=6.23d0
	 ab=10.d0**ab	!car données en LOG10

	CASE('enhan_cha')	!Z/X=4.6E-02, [Fe/H]=7.5-12, Z= 3.207E-02	
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1003) ; WRITE(2,1003)
1003	  FORMAT('Alpha enhanced mixture of Chaboyer 95')	
	 CASE DEFAULT			
	  WRITE(*,3) ; WRITE(2,3)
3	  FORMAT('Mixture alpha enhanced de Chaboyer 95')
	 END SELECT
	 ab(1)=12.d0   ; ab(2)=10.99d0 ; ab(3)=1.16d0  ; ab(4)=1.15d0  ; ab(5)=2.60d0
	 ab(6)=8.55d0  ; ab(7)=7.970d0 ; ab(8)=9.27d0  ; ab(9)=4.56d0  ; ab(10)=8.08d0
	 ab(11)=6.33d0 ; ab(12)=7.98d0 ; ab(13)=6.47d0 ; ab(14)=7.95d0 ; ab(15)=5.45d0
	 ab(16)=7.61d0 ; ab(17)=5.50d0 ; ab(18)=6.52d0 ; ab(19)=5.52d0 ; ab(20)=6.76d0
	 ab(21)=3.17d0 ; ab(22)=5.42d0 ; ab(23)=4.00d0 ; ab(24)=5.67d0 ; ab(25)=4.99d0
	 ab(26)=7.50d0 ; ab(27)=4.92d0 ; ab(28)=6.25d0
	 ab=10.d0**ab	!car données en LOG10

	CASE('enhan_w')		!Z/X=5.531E-02, [Fe/H]=7.5-12, Z= 3.832E-02
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1004) ; WRITE(2,1004)
1004	  FORMAT('Alpha enhanced mixture of Weiss 95')	
	 CASE DEFAULT			
	  WRITE(*,4) ; WRITE(2,4)
4	  FORMAT('Mixture alpha enhanced de Weiss 95')
	 END SELECT
	 ab(1)=12.0d0  ; ab(2)=10.99d0 ; ab(3)=1.16d0  ; ab(4)=1.15d0  ; ab(5)=2.60d0
	 ab(6)=8.55d0  ; ab(7)=7.970d0 ; ab(8)=9.37d0  ; ab(9)=4.56d0  ; ab(10)=8.37d0
	 ab(11)=6.33d0 ; ab(12)=7.98d0 ; ab(13)=6.47d0 ; ab(14)=7.85d0 ; ab(15)=5.45d0
	 ab(16)=7.54d0 ; ab(17)=5.50d0 ; ab(18)=6.52d0 ; ab(19)=5.12d0 ; ab(20)=6.86d0
	 ab(21)=3.17d0 ; ab(22)=5.65d0 ; ab(23)=4.00d0 ; ab(24)=5.67d0 ; ab(25)=5.39d0
	 ab(26)=7.50d0 ; ab(27)=4.92d0 ; ab(28)=6.27d0
	 ab=10.d0**ab	!car données en LOG10	 

	CASE('enhan_al')	!Z/X=2.064E-03, [Fe/H]=6.2-12, Z= 1.484E-03	
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1005) ; WRITE(2,1005)
1005	  FORMAT('Alpha enhanced mixture of Allard 95')	
	 CASE DEFAULT
	  WRITE(*,5) ; WRITE(2,5)
5	  FORMAT('Mixture alpha enhanced de Allard 95')
	 END SELECT
	 ab(1)=12.0d0  ; ab(2)=10.99d0 ; ab(3)=-0.14d0 ; ab(4)=-0.15d0 ; ab(5)=1.30d0
	 ab(6)=7.25d0  ; ab(7)=6.670d0 ; ab(8)=7.87d0  ; ab(9)=3.26d0  ; ab(10)=7.08d0
	 ab(11)=5.33d0 ; ab(12)=6.58d0 ; ab(13)=4.87d0 ; ab(14)=6.55d0 ; ab(15)=4.45d0
	 ab(16)=6.21d0 ; ab(17)=4.50d0 ; ab(18)=5.52d0 ; ab(19)=3.82d0 ; ab(20)=5.36d0
	 ab(21)=1.87d0 ; ab(22)=4.02d0 ; ab(23)=2.70d0 ; ab(24)=4.37d0 ; ab(25)=3.94d0
	 ab(26)=6.20d0 ; ab(27)=3.92d0 ; ab(28)=4.95d0
	 ab=10.d0**ab	!car données en LOG10	 

	CASE('mixture')
	 chain=TRIM(nom_fich2)//'.mix'	 
	 INQUIRE(file=TRIM(chain),exist=ok)	 
	 IF(.NOT.ok)THEN
	  SELECT CASE(langue)	 
	  CASE('english')	
	   WRITE(*,1008)TRIM(chain) ; WRITE(2,1008)TRIM(chain)
1008	   FORMAT('The file of initial abundances : ',a,/,
	1  'is unknown, CESAM looks for the file mixture') 	 
	  CASE DEFAULT
	   WRITE(*,8)TRIM(chain) ; WRITE(2,8)TRIM(chain)
8	   FORMAT('le fichier des abondances initiales : ',a,/,
	1 'non trouvé, recherche du fichier mixture')
	  END SELECT		 
	  chain='mixture' ; INQUIRE(file=TRIM(chain),exist=ok)
	  IF(.NOT.ok)THEN
	   SELECT CASE(langue)	  
	   CASE('english')	
	    WRITE(*,1009)TRIM(chain) ; WRITE(2,1009)TRIM(chain)
1009	    FORMAT('STOP, the file : ',a,/,'is unknown')	  
	   CASE DEFAULT	  
	    WRITE(*,9)TRIM(chain) ; WRITE(2,9)TRIM(chain)
9	    FORMAT('ARRET, le fichier : ',a,/,'non trouvé')
	   END SELECT
	   STOP 	   
	  ENDIF
	 ENDIF
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1007)TRIM(chain) ; WRITE(2,1007)TRIM(chain)	 
1007	  FORMAT('Initial mixture in DeX read on the file : ',a)
	 CASE DEFAULT	 
	  WRITE(*,7)TRIM(chain) ; WRITE(2,7)TRIM(chain)
7	  FORMAT('Mixture initiale en DeX du fichier : ',a)
	 END SELECT
	 READ(3,nl_mixture) ; WRITE(*,nl_mixture) ; WRITE(2,nl_mixture)
	 CLOSE(unit=3)
	 IF(ab(1) /= 12.d0)THEN
	  SELECT CASE(langue)	  
	  CASE('english')	
	   WRITE(*,1013) ; WRITE(2,1013)
1013	   FORMAT('STOP, the initial mixture is not normalized to H=12') 
	  CASE DEFAULT	  
	   WRITE(*,13) ; WRITE(2,13)
13	   FORMAT('ARRET, la mixture initiale doit être normalisée à H=12')
	 END SELECT
	 STOP 	   
	 ENDIF	 
	 ab=10.d0**ab	!car données en LOG10
	CASE DEFAULT
	 chain=nom_abon
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1011)TRIM(chain) ; WRITE(2,1011)TRIM(chain)
1011	  FORMAT('STOP, the initial mixture : ',a,/,'is unknown',/,
	1 'known initial mixtures : meteorites_ag, meteorites_gs',/,
	2 'solaire_gn, enhan_cha, enhan_w, enhan_al, mixture')	  
	 CASE DEFAULT	  
	  WRITE(*,11)TRIM(chain) ; WRITE(2,11)TRIM(chain)
11	  FORMAT('ARRET, nom de mixture initiale inconnue: ',a,/,
	1 'mixtures connues : meteorites_ag, meteorites_gs',/,
	2 'solaire_gn, enhan_cha, enhan_w, enhan_al, mixture') 
	 END SELECT
	 STOP 	   
	END SELECT
	
! faut-il adapter la mixture initiale suivant un fichier type modif_mix?
! pour sécurité, on ne recherche ce fichier que si modif_chim
	INQUIRE(file=TRIM(chain),exist=ok) !; PRINT*,ok ; PAUSE'ok1'	 
	IF(.NOT.ok)THEN
	 chain='modif_mix'
	 INQUIRE(file=TRIM(chain),exist=ok) !; PRINT*,ok ; PAUSE'ok2'
	ENDIF
	IF(ok)THEN
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1012)TRIM(nom_abon),TRIM(chain)
	  WRITE(2,1012)TRIM(nom_abon),TRIM(chain)
1012	  FORMAT(/,'The initial mixture : ',a,/,
	1 'is modified according to the data of the file : ',a) 	 
	 CASE DEFAULT
	  WRITE(*,12)TRIM(nom_abon),TRIM(chain)
	  WRITE(2,12)TRIM(nom_abon),TRIM(chain)
12	  FORMAT(/,'La mixture initiale : ',a,/,
	1 'est adaptée suivant les données du fichier : ',a)
	 END SELECT
	 READ(3,nl_modif_mix) ; WRITE(*,nl_modif_mix)
	 WRITE(2,nl_modif_mix) ; CLOSE(unit=3)	  
	 add=0.d0 ; add(3)=add_Li ; add(4)=add_Be
	 add(5)=add_B ; add(6)=add_C ; add(7)=add_N ; add(8)=add_O
	 add(9)=add_F ; add(10)=add_Ne ; add(11)=add_Na ; add(12)=add_Mg
	 add(13)=add_Al ; add(14)=add_Si ; add(15)=add_P ; add(16)=add_S
	 add(17)=add_Cl ; add(18)=add_Ar ; add(19)=add_K ; add(20)=add_Ca
	 add(21)=add_Sc ; add(22)=add_Ti ; add(23)=add_V ; add(24)=add_Cr
	 add(25)=add_Mn ; add(26)=add_Fe ; add(27)=add_Co
	 add(28)=add_Ni ; add(1)=0.d0	 
	 ab=LOG10(ab/ab(1))+12.d0	!normalisation à H=1.d12
	 ab=ab+add		!modification des abondances en DeX
	 ab(3:28)=ab(3:28)+add_Z
	 ab=10.d0**ab	!car données en LOG10
	
! mesure de sécurité : le fichier modif_mix ne peut être pris en compte que si
! modif_chim = .TRUE.
	 IF(.NOT.modif_chim)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1010) ; WRITE(2,1010)
1010	   FORMAT('STOP, with a file modif_mix, modif_chim = .TRUE. is needed',/,
	1  'erase the file modif_mix, or modify the input file')
	  CASE DEFAULT
	   WRITE(*,10) ; WRITE(2,10)
10	   FORMAT('ARRET, avec un fichier modif_mix il faut modif_chim = .TRUE.',/,
	1  'supprimer le fichier modif_mix , ou modifier le fichier de données')
	  END SELECT
	  CALL sortie	 
	 ENDIF
	ENDIF	

! si garde_xish = .TRUE. les rapports Xi/H seront ceux de la
!	mixture, l'abondance d'hélium sera Y0 lue dans le fichier *.don.
!	X0 et Z0 seront déduits de la mixture
! si garde_xish = .FALSE. les rapports Xi/H seront modifiés
!	les rapports Xi/Z seront ceux de la
!	mixture, X0, Y0, Z0 seront ceux du fichier *.don
!	(exemple : la calibration solaire en Z/X)

	IF(garde_xish)THEN	
! ajustement de ab(2) pour obtenir Y=Y0		
	 nisnh=ab*m*1.d-12
	 ab(2)=y0/(1.d0-y0)/m(2)*(m(1)+SUM(nisnh(3:28)))*1.d12	 
!	 WRITE(*,2000)ab(2) ; PAUSE'ab(2)'

! détermination de X0, Y0, Z0
	 sniai=SUM(ab*m)
	 x0=ab(1)*m(1)/sniai ; y0=ab(2)*m(2)/sniai ; z0=1.d0-x0-y0
	 zsx=z0/x0 ; fesh=LOG10(zsx/zsx_sol)
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1014)x0,z0,zsx,fesh,y0
	  WRITE(2,1014)x0,z0,zsx,fesh,y0
1014	  FORMAT(/,'X0:',es10.3,', Z0:',es10.3,', Z0/X0:',es10.3,
	1 ', [Z/X]:',es10.3,/,'are different from the input data',/,
	2 'the ratios metal/H of the mixture are used, Y0=',es10.3,
	3 ' is the input value') 	  
	 CASE DEFAULT	 
	  WRITE(*,14)x0,z0,zsx,fesh,y0
	  WRITE(2,14)x0,z0,zsx,fesh,y0
14	  FORMAT(/,'X0:',es10.3,', Z0:',es10.3,', Z0/X0:',es10.3,
	1 ', [Z/X]:',es10.3,/,'diffèrent des données',/,
	2 'les rapports métaux/H de la mixture sont conservés, Y0=',es10.3,
	2 ' est celui des données')
	 END SELECT
	ELSE
	 SELECT CASE(langue)	  
	 CASE('english')
c	  WRITE(*,1018)x0,y0,z0,z0/x0
c	  WRITE(2,1018)x0,y0,z0,z0/x0
1018	  FORMAT(/,'One uses the input data X0:',es10.3,', Y0:',es10.3,
	1 ', Z0:',es10.3,/,'Z0/X0:',es10.3,
	2 ' and the ratios metal/Z of the mixture')	
	 CASE DEFAULT	 
c	  WRITE(*,18)x0,y0,z0,z0/x0
c	  WRITE(2,18)x0,y0,z0,z0/x0
18	  FORMAT(/,'On utilise les données X0:',es10.3,', Y0:',es10.3,
	1 ', Z0:',es10.3,/,'Z0/X0:',es10.3,
	2 ' et les rapports métaux/Z de la mixture')
	 END SELECT
	ENDIF
		
	SELECT CASE(langue)	  
	CASE('english')	
c	 WRITE(*,1020) ; WRITE(2,1020)
1020	 FORMAT(/,'mixture used')	
	CASE DEFAULT	 
c	 WRITE(*,20) ; WRITE(2,20)
20	 FORMAT(/,'mixture retenue')	
	END SELECT	 
c	WRITE(*,2000)LOG10(ab) ; WRITE(2,2000)LOG10(ab)

! abondances(en nombre) des métaux dans Z
	abon_rela=ab/SUM(ab(3:nelem_ini))
	SELECT CASE(langue)	  
	CASE('english')	
c	 WRITE(*,1019) ; WRITE(2,1019)
1019	 FORMAT(/,'relative abundances(by number) of metals within Z')
	CASE DEFAULT	 
c	 WRITE(*,19) ; WRITE(2,19)
19	 FORMAT(/,'abondances(en nombre) des métaux dans Z')	
	END SELECT	 
c	WRITE(*,2000)abon_rela(3:nelem_ini)
c	WRITE(2,2000)abon_rela(3:nelem_ini)

! produit abondance * masse pour calcul des abondances relatives en masse
!	et écriture des abondances des métaux dans Z	
	abon_m=ab*m
	abon_relam=abon_m/SUM(abon_m(3:nelem_ini))	
	SELECT CASE(langue)	  
	CASE('english')	
c	 WRITE(*,1025) ; WRITE(2,1025)
1025	 FORMAT(/,'relative abundances(by mass) of metals within Z')
	CASE DEFAULT	 
c	 WRITE(*,25) ; WRITE(2,25)
25	 FORMAT(/,'abondances(en masse) des métaux dans Z')	
	END SELECT	 
c	WRITE(*,2000)abon_relam(3:nelem_ini) ; WRITE(*,*)
c	WRITE(2,2000)abon_relam(3:nelem_ini) ; WRITE(2,*)
	
! définition des rapports isotopiques, recherche d'un fichier de rap.iso.
	chain=TRIM(nom_fich2)//'.rap_iso'	 
	INQUIRE(file=TRIM(chain),exist=ok)	 
	IF(.NOT.ok)THEN  
	 chain='rap_iso'
	 INQUIRE(file=TRIM(chain),exist=ok)
	ENDIF
	
! écritures
	IF(ok)THEN
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 PRINT* ; WRITE(2,*)	!rapports isotopiques du fichier rap_iso
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1016)TRIM(chain)
	  WRITE(2,1016)TRIM(chain)
1016	  FORMAT('Isotopic ratios from the file : ',a)
	 CASE DEFAULT
	  WRITE(*,16)TRIM(chain)
	  WRITE(2,16)TRIM(chain)
16	  FORMAT('Rapports isotopiques du fichier : ',a)
	 END SELECT
	 READ(3,nl_rap_iso) ; WRITE(*,nl_rap_iso) ; WRITE(2,nl_rap_iso) 
	 CLOSE(unit=3)
	 
! mesure de sécurité : le fichier rap_iso ne peut être pris en compte que si
! modif_chim = .TRUE.
	 IF(.NOT.modif_chim)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1027) ; WRITE(2,1027)
1027	   FORMAT('STOP, with a file modif_mix, modif_chim = .TRUE. is needed',/,
	1  'erase the file vent, or modify the input file')
	  CASE DEFAULT
	   WRITE(*,27) ; WRITE(2,27)
27	   FORMAT('ARRET, avec un fichier modif_mix il faut modif_chim = .TRUE.',/,
	1  'supprimer le fichier vent, ou modifier le fichier de données')
	  END SELECT
	  CALL sortie	 
	 ENDIF	 
	ELSE		!rapports isotopiques d'Anders & Grevesse
	 SELECT CASE(langue)
	 CASE('english')
c	  WRITE(*,1006) ; WRITE(2,1006)
1006	  FORMAT('Isotopic ratios according to Anders & Grevesse',/,
	1 'but He3/He4 by Gautier & Morel 97')	
	 CASE DEFAULT
c	  WRITE(*,6) ; WRITE(2,6)
6	  FORMAT('Rapports isotopiques suivant Anders & Grevesse',/,
	1 'sauf He3/He4 de Gautier & Morel 97')
	 END SELECT	
	 h2sh1=3.01d-5           !Gautier & Morel
	 he3she4=1.1d-4          !Gautier & Morel
	 he3she4z=4.185d-4       !Gautier & Morel sur la ZAMS H2-->He3
	 li6sli7=7.5d0/92.5d0    !table 3 de Anders & Grevesse
	 be7sbe9=1.d-25		 !valeur arbitraire non nulle
	 c13sc12=1.10d0/98.90d0      !table 3 de Anders & Grevesse
	 n15sn14=0.366d0/99.634d0    !table 3 de Anders & Grevesse
	 o17so16=0.038d0/99.762d0    !table 3 de Anders & Grevesse
	 o18so16=0.200d0/99.762d0    !table 3 de Anders & Grevesse	 
	 ne22sne20=6.79d0/92.99d0    !table 3 de Anders & Grevesse
	 mg25smg24=10.00d0/78.99d0   !table 3 de Anders & Grevesse
	 mg26smg25=11.01d0/10.00d0     !table 3 de Anders & Grevesse
	 be7sz=1.d-29    		!Be7 arbitraire	
	ENDIF
c	WRITE(*,21)h2sh1,he3she4,he3she4z
c	WRITE(2,21)h2sh1,he3she4,he3she4z	
21	FORMAT('H2/H1=',es10.3,', He3/He4(pms)=',es10.3,', He3/He4(Zams)=',
	1 es10.3)
c	WRITE(*,22)li6sli7,be7sbe9,c13sc12
c	WRITE(2,22)li6sli7,be7sbe9,c13sc12
22	FORMAT('Li6/Li7=',es10.3,', Be7/Be9=',es10.3, ', C13/C12=',es10.3)
c	WRITE(*,23)n15sn14,o17so16,o18so16,ne22sne20
c	WRITE(2,23)n15sn14,o17so16,o18so16,ne22sne20	
23	FORMAT('N15/N14=',es10.3,', O17/O16=',es10.3,', O18/O16=',es10.3,
	1 ', Ne22/Ne20=',es10.3)
c	WRITE(*,24)mg25smg24,mg26smg25,be7sz
c	WRITE(2,24)mg25smg24,mg26smg25,be7sz	
24	FORMAT('Mg25/Mg24=',es10.3,', Mg26/Mg25=',es10.3,', Be7/Z=',
	1 es10.3,/)	
		
	RETURN

	END SUBROUTINE abon_ini
