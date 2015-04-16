
c***********************************************************************

	SUBROUTINE ini_ctes

c routine PUBLIC des modules mod_donnees et mod_exploit

c initialisation de l'ensemble des principales constantes physiques

c constantes physiques selon  CRC Handbook of
c chemistry and physics 75th edition 1994 D.R. Lide et al. CRC press
c Boka Raton : Ann Arbor : London : Tokyo,  p.11.35

c les données de base sont celles de ctes_94
c les corrections nécessaires sont apportées pour les autres ensembles

c la routine print_ctes.f permet de les écrire

c Auteur: P.Morel, Département J.D. Cassini : O.C.A., CESAM2k

c-------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	CHARACTER (len=50) :: chain

c-------------------------------------------------------------------

c le fichier de données
	chain=TRIM(nom_fich2)//'.don'

c origine des principales données	
	source='Handbook chem. phys. 94 + NACRE'

c constantes de base début: ctes_94 c'est le défaut----------------------

c constantes de base
	amu=1.6605402d-24		!masse atom. unité, Avogadro=1/amu
	aradia=7.565912199839849d-15	!constante de la radiation
	clight=2.99792458d10		!célérité de la lumière
	echarg=4.8032068d-10		!charge de l'électron
	eve=1.60217733d-12		!électron volt
	gmsol=1.32712438d26		!produit G Msol, cte. de Newton
	hpl=6.6260755d-27		!planck	
	kbol=1.380658d-16		!boltzman
	granr=kbol/amu			!constante des gaz parfaits	
	ln10=LOG(10.d0)
	lbol0=3.055d33			!point 0 des Mbol	
	me=9.1093897d-28		!masse électron			
	pi=acos(-1.d0)
	secon6=365.25d0*86400.d6	!nb. de s. dans 1Myr (norme de 1967)
	sigma=aradia*clight/4.d0	!constante de Stefan
	zsx_sol=0.0245
	fesh_sol=7.5d0-12.00d0		![Fe/H] de GN 93
	ua=149.59787066d6*1.d5		!unité astronomique

c constantes astro		
	lsol=3.846d33			!luminosité solaire		
	msol=1.98919d33			!masse solaire	
	rsol=6.9599d10			!rayon solaire
	mterre=5.977d27			!masse terrestre

c la gravité dépend de la valeur retenue pour la masse solaire	
	g=gmsol/msol

c masses des isotopes en unités atomiques	
	an=1.008665d0		!masse atomique en amu du neutron
	ap=1.00727647d0		!masse atomique en amu du proton
	ah=1.007825d0		!masse atomique en amu de l'hydrogène
	ah2=2.0141018d0		!masse atomique en amu du deutérium
	ahe3=3.0160293d0	!masse atomique en amu de l'hélium 3
	ahe4=4.0026033d0	!masse atomique en amu de l'hélium 4
	ali6=6.015121d0		!masse atomique en amu du lithium 6	
	ali7=7.0160040d0	!masse atomique en amu du lithium 7	
	abe7=7.0169292d0	!masse atomique en amu du béryllium 7
	abe9=9.0121821d0	!masse atomique en amu du béryllium 9
	ab11=11.0093055d0       !masse atomique en amu du bore 11	
	ac12=12.d0		!masse atomique en amu du carbone 12
	ac13=13.0033548d0	!masse atomique en amu du carbone 13
	an13=13.0057386d0	!masse atomique en amu de l'azote 13
	an14=14.003074d0	!masse atomique en amu de l'azote 14
	an15=15.001089d0	!masse atomique en amu de l'azote 15
	ao16=15.9949146d0	!masse atomique en amu de l'oxygène 16
	ao17=16.9991315d0	!masse atomique en amu de l'oxygène 17
	ao18=17.9991604d0	!masse atomique en amu de l'oxygène 18
	afe56=55.847d0		!masse atomique en amu du fer 56	
	af18=18.0009377d0	!masse atomique en amu du fluor 18
	af19=18.9984032d0	!masse atomique en amu du fluor 19
	ane20=19.9924402d0	!masse atomique en amu du néon 20
	ane21=20.9938467d0	!masse atomique en amu du néon 21
	ane22=21.9913855d0	!masse atomique en amu du néon 22
	ana23=22.9897697d0	!masse atomique en amu du sodium 23
	amg23=22.9941249d0	!masse atomique en amu du magnésium 23
	amg24=23.9850419d0	!masse atomique en amu du magnésium 24
	amg25=24.985837d0	!masse atomique en amu du magnésium 25
	amg26=25.982593d0	!masse atomique en amu du magnésium 26
	asi28=28.0855d0		!masse atomique en amu du silicium 28
	aal27=26.9854d0		!masse atomique en amu de l'aluminium 27
	as32=31.972070d0	!masse atomique en amu du soufre 32
	ap31=30.973762d0	!masse atomique en amu du phosphore 31
	
c constantes de base fin ----------------------

c commentaires
	SELECT CASE(langue)
	CASE('english')	
	 WRITE(*,1001)nom_ctes ; WRITE(2,1001)nom_ctes
1001	 FORMAT('Use of the physical constants of the set ',a)
	CASE DEFAULT
	 WRITE(*,1)nom_ctes ; WRITE(2,1)nom_ctes
1	 FORMAT('Utilisation des constantes physiques de l''ensemble: ',a)
	END SELECT

c corrections suivant les ensembles, pas de corrections pour ctes_94
	SELECT CASE(nom_ctes)

c constantes de Toulouse (département de haute Garonne 31)
	CASE ('ctes_31')
	 lsol=3.815d33			!luminosité solaire		
	 rsol=6.9575d10			!rayon solaire

c constantes du solar model comparison project (JCD)
	CASE ('ctes_85')	
	 msol=1.9891d33			!masse solaire	
	 g=gmsol/msol			!gravité

c constantes avec valeurs entieres pour les masses atomiques des isotopes	 
	CASE ('ctes_94m')
	 an=1.d0		!masse atomique en amu du neutron
	 ap=1.d0		!masse atomique en amu du proton
	 ah=1.d0		!masse atomique en amu de l'hydrogène
	 ah2=2.d0		!masse atomique en amu du deutérium
	 ahe3=3.d0	!masse atomique en amu de l'hélium 3
	 ahe4=4.d0	!masse atomique en amu de l'hélium 4
	 ali6=6.d0		!masse atomique en amu du lithium 6	
	 ali7=7.d0	!masse atomique en amu du lithium 7	
	 abe7=7.d0	!masse atomique en amu du béryllium 7
	 abe9=9.d0	!masse atomique en amu du béryllium 9
	 ab11=11.d0       !masse atomique en amu du bore 11	
	 ac12=12.d0		!masse atomique en amu du carbone 12
	 ac13=13.d0	!masse atomique en amu du carbone 13
	 an13=13.d0	!masse atomique en amu de l'azote 13
	 an14=14.d0	!masse atomique en amu de l'azote 14
	 an15=15.d0	!masse atomique en amu de l'azote 15
	 ao16=16.d0	!masse atomique en amu de l'oxygène 16
	 ao17=17.d0	!masse atomique en amu de l'oxygène 17
	 ao18=18.d0	!masse atomique en amu de l'oxygène 18
	 afe56=56.d0		!masse atomique en amu du fer 56	
	 af18=18.d0	!masse atomique en amu du fluor 18
	 af19=19.d0	!masse atomique en amu du fluor 19
	 ane20=20.d0	!masse atomique en amu du néon 20
	 ane21=21.d0	!masse atomique en amu du néon 21
	 ane22=22.d0	!masse atomique en amu du néon 22
	 ana23=23.d0	!masse atomique en amu du sodium 23
	 amg23=23.d0	!masse atomique en amu du magnésium 23
	 amg24=24.d0	!masse atomique en amu du magnésium 24
	 amg25=25.d0	!masse atomique en amu du magnésium 25
	 amg26=26.d0	!masse atomique en amu du magnésium 26
	 asi28=28.d0		!masse atomique en amu du silicium 28
	 aal27=27.d0		!masse atomique en amu de l'aluminium 27
	 as32=32.d0	!masse atomique en amu du soufre 32
	 ap31=31.d0	!masse atomique en amu du phosphore 31	

c constantes de Basu & Antia
	CASE ('ctes_ba')
	 lsol=3.8418d33			!luminosité solaire		
	 msol=1.98892d33		!masse solaire	
	 g=gmsol/msol			!gravité

c constantes de GAIA
	CASE ('ctes_gaia')
	 lsol=3.8515d33			!luminosité solaire		
	 msol=1.9884d33			!masse solaire	
	 rsol=6.9599767d10		!rayon solaire	
	 g=gmsol/msol			!gravité
  
	CASE default
	 IF(TRIM(nom_ctes) /= 'ctes_94')THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1002)nom_ctes,chain ; WRITE(2,1002)nom_ctes,chain
1002	   FORMAT('STOP, unkown set of physical constants: ',a,/
	1  'choose among: ctes_85, ctes_94, ctes_94m, ctes_31, ctes_gaia',/,
	2  'and correct the input data file ',a)
	  CASE DEFAULT
	   WRITE(*,2)nom_ctes,chain ; WRITE(2,2)nom_ctes,chain
2	   FORMAT('ARRET, ensemble de constantes physiques inconnu: ',a,/,
	1  'Ensembles reconnus: ctes_85, ctes_94, ctes_94m, ctes_ba, ctes_31, ctes_gaia',/,
	2  'corriger le fichier de données ',a)
	  END SELECT
	  STOP
	 ENDIF
	END SELECT	

	RETURN

	END SUBROUTINE ini_ctes
