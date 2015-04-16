
c******************************************************************

	SUBROUTINE modif_mix

c	routine private du module mod_nuc
c	personnalisation d'une mixture

c	on recherche le fichier mon_modele.modif_mix
c	s'il existe, on y lit des modifications en DeX des abondances
c	initiales. On modifie les abondances des éléments lourds

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : langue, nom_abon, nom_fich2, x0, y0, z0,
	1 zsx_sol
	USE mod_kind

	IMPLICIT NONE
	
	REAL(kind=dp), DIMENSION(28) :: add	
	REAL(kind=dp) :: add_H,add_He,add_Li,add_Be,add_B,add_C,
	1 add_N,add_O,add_F,add_Ne,add_Na,add_Mg,add_Al,add_Si,add_P,
	2 add_S,add_Cl,add_Ar,add_K,add_Ca,add_Sc,add_Ti,add_V,add_Cr,
	3 add_Mn,add_Fe,add_Co,add_Ni,fesh, sniai, zsx
	
	NAMELIST/nl_modif_mix/add_H,add_He,add_Li,add_Be,add_B,add_C,
	1 add_N,add_O,add_F,add_Ne,add_Na,add_Mg,add_Al,add_Si,add_P,
	2 add_S,add_Cl,add_Ar,add_K,add_Ca,add_Sc,add_Ti,add_V,add_Cr,
	3 add_Mn,add_Fe,add_Co,add_Ni

c-------------------------------------------------------------------

c faut-il adapter la mixture initiale?
	
	IF(TRIM(nom_abon) == 'mixture')RETURN
	
c recherche du fichier mon_modele.modif_mix
	
	chain=TRIM(nom_fich2)//'.modif_mix'	 
	INQUIRE(file=TRIM(chain),exist=ok)	 
	IF(ok)THEN
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1012)TRIM(nom_abon),TRIM(chain)
	  WRITE(2,1012)TRIM(nom_abon),TRIM(chain)
1012	  FORMAT('The initial mixture : ',a,/,
	1 'is changed according to the data of the file : ',a) 	 
	 CASE DEFAULT
	  WRITE(*,12)TRIM(nom_abon),TRIM(chain)
	  WRITE(2,12)TRIM(nom_abon),TRIM(chain)
12	  FORMAT('La mixture initiale : ',a,/,
	1 'est adaptée suivant les données du fichier : ',a)
	 END SELECT
	 READ(3,nl_modif_mix) ; WRITE(*,nl_modif_mix)
	 WRITE(2,nl_modif_mix)
	 CLOSE(unit=3)	  
	 add(1)=add_H ; add(2)=add_He ; add(3)=add_Li ; add(4)=add_Be
	 add(5)=add_B ; add(6)=add_C ; add(7)=add_N ; add(8)=add_O
	 add(9)=add_F ; add(10)=add_Ne ; add(11)=add_Na ; add(12)=add_Mg
	 add(13)=add_Al ; add(14)=add_Si ; add(15)=add_P ; add(16)=add_S
	 add(17)=add_Cl ; add(18)=add_Ar ; add(19)=add_K ; add(20)=add_Ca
	 add(21)=add_Sc ; add(22)=add_Ti ; add(23)=add_V ; add(24)=add_Cr
	 add(25)=add_Mn ; add(26)=add_Fe ; add(27)=add_Co
	 add(28)=add_Ni

c	modification des abondances en DeX
	 
	 ab=ab+add	 
	 
	ELSE
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1013)TRIM(chain),TRIM(nom_abon)
	  WRITE(2,1013)TRIM(chain),TRIM(nom_abon)
1013	  FORMAT('The file : ',a,/,'is unknown, the mixture : ',a,/,
	1 'is kept') 	 
	 CASE DEFAULT
	  WRITE(*,13)TRIM(chain),TRIM(nom_abon)
	  WRITE(2,13)TRIM(chain),TRIM(nom_abon)
13	  FORMAT('Pas de fichier : ',a,/,'on conserve la mixture : ',a)
	 END SELECT
	ENDIF

	RETURN

	END SUBROUTINE modif_mix

