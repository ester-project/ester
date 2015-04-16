
c**************************************************************************

	SUBROUTINE planetoides(xchim,dxchim,jac,m_planet,mw_planet)
	
c routine public du module mod_nuc

c modifie les dérivées temporelles de la composition chimique
c calcule le moment cinétique par unité de masse apporté par les planétoîdes

c planétoides initialisé dans cesam, n'est appelé par coeff_rota3/4
c et les routines de perte de masse,
c que s'il y a chute de planétoïdes, alors l_planet=.TRUE.

c entrées
c	xchim : composition chimique

c entrées/sorties
c	dxchim : dérivée temporelle modifiée
c	jac : jacobien

c sorties
c	m_planet : flux de masse en Msol/My
c	mw_planet(2) : flux de moment cinétique	

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ihe4, i_ex, Krot, langue, modif_chim, msol,
	1 mterre, nchim, nom_elem, nom_fich2, nucleo, pi, secon6, rsol, ua
	USE mod_kind	
	USE mod_variables, ONLY : age, lim, mstar, m_zc, sortie, x_planet

	IMPLICIT NONE
	
	REAL (kind=dp), OPTIONAL, INTENT(in), DIMENSION(:) :: xchim		
	REAL (kind=dp), OPTIONAL, INTENT(inout), DIMENSION(:,:) :: jac
	REAL (kind=dp), OPTIONAL, INTENT(inout), DIMENSION(:) :: dxchim
	REAL (kind=dp), OPTIONAL, INTENT(out), DIMENSION(2) :: mw_planet	
	REAL (kind=dp), OPTIONAL, INTENT(out) :: m_planet
	
	REAL (kind=dp), DIMENSION(28) :: vt
	
	REAL (kind=dp), PARAMETER :: beta=3.d0, npas=10.d0		
	REAL (kind=dp), SAVE ::	flux0, mw_planet0=0.d0,
	1 n_planet, sig_gauss, t12	
	REAL (kind=dp) ::  flux, m_pla, r_giration, t_giration, vt_H, vt_He,
	1 vt_Li, vt_Be, vt_B, vt_C, vt_N, vt_O, vt_F, vt_Ne, vt_Na, vt_Mg,
	2 vt_Al, vt_Si, vt_P, vt_S, vt_Cl, vt_Ar, vt_K, vt_Ca, vt_Sc, vt_Ti,
	3 vt_V,vt_Cr, vt_Mn, vt_Fe, vt_Co, vt_Ni, xp, xpl, yp, ypl, zp, zpl

	INTEGER :: i, ili7=0
		
	LOGICAL, SAVE :: init=.TRUE.
	
	CHARACTER (len=20), SAVE :: profil
	CHARACTER (len=50) :: chain
	
	NAMELIST/nl_planet/ vt_H, vt_He, vt_Li, vt_Be, vt_B, vt_C,
	1 vt_N, vt_O, vt_F, vt_Ne, vt_Na, vt_Mg, vt_Al, vt_Si,
	2 vt_P, vt_S, vt_Cl, vt_Ar, vt_K, vt_Ca, vt_Sc, vt_Ti,
	3 vt_V,vt_Cr, vt_Mn, vt_Fe, vt_Co, vt_Ni, ypl, zpl, n_planet, 
	4 profil, age_deb, age_fin,r_giration,t_giration
		
c----------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(10es8.1)

	IF(init)THEN
	 init=.FALSE.
	 
c recherche d'une namelist de composition chimique des planétoïdes
	 chain=TRIM(nom_fich2)//'.planet'	 
	 INQUIRE(file=TRIM(chain),exist=l_planet)	 
	 IF(.NOT.l_planet)THEN		!l_planet 1
	  chain='planet'	 
	  INQUIRE(file=TRIM(chain),exist=l_planet)
	  IF(.NOT.l_planet)THEN		!l_planet 2
	   SELECT CASE(langue)
	   CASE('english')	
	    WRITE(*,1007) ; WRITE(2,1007)
1007	    FORMAT('No infall of planets')
	   CASE DEFAULT
	    WRITE(*,7) ; WRITE(2,7)
7	    FORMAT('Pas de chute de planétoïdes')
	   END SELECT
	   RETURN
	  ENDIF				!l_planet 2
	 ENDIF				!l_planet 1
	 	  
c lecture de la namelist de la composition chimique des planetoïdes	  
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 WRITE(*,*) ; WRITE(2,*)	
	 READ(3,nl_planet) ; WRITE(*,nl_planet) ; WRITE(2,nl_planet)
	 WRITE(*,*) ; WRITE(2,*)
	 
c normalisation	 
	 vt(1)=vt_H ; vt(2)=vt_He ; vt(3)=vt_Li ; vt(4)=vt_Be
	 vt(5)=vt_B ; vt(6)=vt_C ; vt(7)=vt_N ; vt(8)=vt_O
	 vt(9)=vt_F ; vt(10)=vt_Ne ; vt(11)=vt_Na ; vt(12)=vt_Mg
	 vt(13)=vt_Al ; vt(14)=vt_Si ; vt(15)=vt_P ; vt(16)=vt_S
	 vt(17)=vt_Cl ; vt(18)=vt_Ar ; vt(19)=vt_K ; vt(20)=vt_Ca
	 vt(21)=vt_Sc ; vt(22)=vt_Ti ; vt(23)=vt_V ; vt(24)=vt_Cr
	 vt(25)=vt_Mn ; vt(26)=vt_Fe ; vt(27)=vt_Co ; vt(28)=vt_Ni

c normalisation par masse  	 
	 vt=vt/SUM(vt)
	  
c affectation des abondances par masse des principaux isotopes
c ATTENTION l'algorithme n'est totalement indépendant
c de l'ordre et de la présence des éléments	  
	 ALLOCATE(x_planet(nchim)) ; x_planet=0.d0

c hydrogène, deutérium et hélium 2 et 4	  
	 IF(nom_elem(2) == ' H2 ')THEN
	  x_planet(1)=vt(1)*(1.d0-h2sh1) ; x_planet(2)=vt(1)*h2sh1
	  x_planet(ihe4-1)=vt(2)*he3she4 ; x_planet(ihe4)=vt(2)*(1.d0-he3she4)
	 ELSE
	  x_planet(1)=vt(1) ; xp=x_planet(1)
	  
c on suppose H2 dans He3
	  x_planet(ihe4-1)=vt(2)*he3she4z ; x_planet(ihe4)=vt(2)*(1.d0-he3she4z)
	 ENDIF

c autres éléments 
	 i=ihe4
	 B1: DO	   
	  i=i+1 ; IF(i > nchim)EXIT B1
	  IF(nom_elem(i) == 'Li7 ')THEN
	   x_planet(i)=vt(3) ; ili7=i ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Li6 ')THEN
	   IF(ili7 /= 0)THEN
	    x_planet(i)=vt(3)*li6sli7 ; x_planet(ili7)=vt(3)*(1.d0-li6sli7)
	    CYCLE B1
	   ELSE
	    x_planet(i)=vt(3) ; CYCLE B1
	   ENDIF  
	  ELSEIF(nom_elem(i) == 'Be9 ')THEN
	   x_planet(i)=vt(4) ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'B11 ')THEN
	   x_planet(i)=vt(5) ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'C12 ')THEN
	   x_planet(i)=vt(6)*(1.d0-c13sc12) ; x_planet(i+1)=vt(6)*c13sc12
	   i=i+1 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'N14 ')THEN
	   x_planet(i)=vt(7)*(1.d0-n15sn14) ; x_planet(i+1)=vt(7)*n15sn14
	   i=i+1 ; CYCLE B1  
	  ELSEIF(nom_elem(i) == 'O16 ')THEN
	   x_planet(i)=vt(8)*(1.d0-o17so16-o18so16)
	   x_planet(i+1)=vt(8)*o17so16
	   x_planet(i+2)=vt(8)*o18so16
	   i=i+1 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'O18 ')THEN
	   x_planet(i)=vt(8)*o18so16 ; CYCLE B1	    
	  ELSEIF(nom_elem(i) == 'Ne20')THEN
	   x_planet(i)=vt(10)*(1.d0-ne22sne20) ; x_planet(i+1)=vt(10)*ne22sne20
	   i=i+1 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Mg24')THEN	    
	   x_planet(i)=vt(12)*(1.d0-mg25smg24*(1.d0+mg26smg25)) ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Mg25')THEN	    
	   x_planet(i)=vt(12)*mg25smg24 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Mg26')THEN	    
	   x_planet(i)=vt(12)*mg25smg24*(1.d0+mg26smg25) ; CYCLE B1    
	  ELSEIF(nom_elem(i) == 'Fe56')THEN
	   x_planet(i)=vt(26) ; CYCLE B1	    
	  ENDIF	     
	 ENDDO B1	  
	  
c affectation de l'abondances par masse pour l'extra élément Ex
	 x_planet(i_ex)=1.d0-SUM(x_planet)
	 
c renormalisation	
	 x_planet=x_planet/SUM(x_planet)
		 
c X, Y, Z pour les planétoïdes
	 xp=SUM(x_planet(1:ihe4-2))
	 yp=SUM(x_planet(ihe4-1:ihe4))
	 zp=1.d0-xp-yp
	 
c Avec xp, yp, zp non nuls on modifie les abondances relatives pour obtenir
c ypl, zpl, xpl=1-ypl-zpl du fichier planet

	 IF(xp*yp*zp > 0.d0)THEN
	  xpl=1.d0-ypl-zpl
	  x_planet(1:ihe4-2)=x_planet(1:ihe4-2)*xpl/xp
	  xp=SUM(x_planet(1:ihe4-2))	 
	  x_planet(ihe4-1:ihe4)=x_planet(ihe4-1:ihe4)*ypl/yp
	  yp=SUM(x_planet(ihe4-1:ihe4))
	  x_planet(ihe4+1:nchim)=x_planet(ihe4+1:nchim)*zpl/zp
	  zp=1.d0-xp-yp
	 ENDIF	  
	  
c abondances par mole
	 x_planet(1:nchim)=x_planet(1:nchim)/nucleo(1:nchim)
	  
c flux de masse total en Msol/My
	 flux0=mterre*n_planet/msol/(age_fin-age_deb)	 

c cas de la gaussiènne
	 IF(profil == 'gauss')THEN
	  sig_gauss=(age_fin-age_deb)/beta
	  flux0=flux0*(age_fin-age_deb)/SQRT(2.d0*pi)/sig_gauss
	  t12=(age_deb+age_fin)/2.d0
	 ENDIF

c moment cinétique apporté par les planétoïdes
c un temps de giration négatif correspond à un apport rétrograde	
	 IF( (Krot == 3 .OR. Krot == 4) .AND. t_giration /= 0.d0)THEN

c le rayon de giration est en UA, le temps de giration est le temps en années
c nécessaire pour une rotation de pi/2
c le temps de giration est donné en années d'où 1.d6/secon6
c mw_planet0: moment cinétique par unité de masse et de temps apporté par les
c planétoïdes dont la vitesse angulaire est non nulle
c mw_planet0 est en cgs
	  mw_planet0=pi/2.d0*1.d6/secon6*(r_giration*ua)**2/t_giration
c	  WRITE(*,2000)mw_planet0,flux0 ; PAUSE'mw_planet0'

	 ENDIF
	 
c pas temporel maximal pendant l'accrétion des planétoïdes
	 dt_planet=ABS(age_fin-age_deb)/npas
	 
c	 WRITE(*,2000)mterre,n_planet,age,age_deb,age_fin,flux0
c	 WRITE(*,2000)x_planet
c	 PRINT*,profil ; PAUSE'initialisation'

c écritures	  
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1004)TRIM(chain),xp,yp,zp ; WRITE(2,1004)TRIM(chain),xp,yp,zp 
1004	  FORMAT('Infalls of planetoids, abundances from the file :',a,/,
	1 'X=',es10.3,', Y=',es10.3,', Z=',es10.3,/,'abundances per mole :')
	 CASE DEFAULT
	  WRITE(*,4)TRIM(chain),xp,yp,zp ; WRITE(2,4)TRIM(chain),xp,yp,zp
4	  FORMAT('Chutes de planétoïdes, abondances du fichier :',a,/,
	1 'X=',es10.3,', Y=',es10.3,' ,Z=',es10.3,/,'abondances par mole :')
	 END SELECT
	 WRITE(*,6)(nom_elem(i),i=1,MIN(nchim,8))
	 WRITE(2,6)(nom_elem(i),i=1,MIN(nchim,8))	  
	 WRITE(*,5)(x_planet(i),i=1,MIN(nchim,8))
	 WRITE(2,5)(x_planet(i),i=1,MIN(nchim,8))
5	 FORMAT(1x,12es10.3)
	 IF(nchim > 8)THEN	  	
	  WRITE(*,6)(nom_elem(i),i=9,MIN(nchim,16))
	  WRITE(2,6)(nom_elem(i),i=9,MIN(nchim,16))
6	  FORMAT(12(3x,a4,3x))
	  WRITE(*,5)(x_planet(i),i=9,MIN(nchim,16))
	  WRITE(2,5)(x_planet(i),i=9,MIN(nchim,16))
	 ENDIF
	 	 
c mesure de sécurité : les chutes de planétoïdes ne peuvent être pris en
c compte que si modif_chim = .TRUE.
	 IF(.NOT.modif_chim)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1009) ; WRITE(2,1009)
1009	   FORMAT('STOP, with a file planet, modif_chim = .TRUE. is needed',/,
	1  'erase the file planet, or modify the input file')
	  CASE DEFAULT
	   WRITE(*,9) ; WRITE(2,9)
9	   FORMAT('ARRET, avec un fichier planet il faut modif_chim = .TRUE.',/,
	1  'supprimer le fichier planet, ou modifier le fichier de données')
	  END SELECT
	  CALL sortie	 
	 ENDIF	 
	 
	 RETURN	 
	ENDIF		!init
	
c---------------fin des initialisations-------------------------------	
	
c pendant l'intervalle de chute
	IF(age >= age_deb .AND. age < age_fin)THEN
	
c masse de planétoïdes (Msol) par unité de temps (My) suivant le type de profil
c flux0: flux de masse total en Msol/My
	 SELECT CASE(profil)
	 CASE('rectangle')	 
	  m_pla=flux0
	 CASE('triangle')	 
	  m_pla=flux0*2.d0*(1.d0-ABS(2.d0*age-age_deb-age_fin)/
	1 (age_fin-age_deb))
	 CASE('parabole')
	  m_pla=flux0*6.d0*(age-age_deb)*(age_fin-age)/(age_fin-age_deb)**2
	 CASE('gauss')
	  m_pla=flux0*EXP(-((age-t12)/sig_gauss)**2/2.d0)
	 CASE DEFAULT	 
	  SELECT CASE(langue)
	  CASE('english')	
	  CASE DEFAULT
	   WRITE(*,1008)profil ; WRITE(2,1008)profil	 
1008	   FORMAT('ERROR, unknown shape :',a,/,
	1  'Shapes known : rectangle, triangle, parabole, gauss',/,
	2  'The shape rectangle is used')
	   WRITE(*,8)profil ; WRITE(2,8)profil	 
8	   FORMAT('ERREUR, profil inconnu :',a,/,
	1  'Profils connus : rectangle, triangle, parabole, gauss',/,
	2  'On utilise le profil rectangle')	  
	  END SELECT
	  m_pla=flux0	  
	 END SELECT

c masse (Msol) de l'accrétion pendant l'unité de temps (My)
	 IF(PRESENT(m_planet))m_planet=m_pla	 
	 
c accrétion par unité de masse de ZC, par Myr, mzc_ext masse (Msol) de la ZC ext
	 flux=m_pla/mzc_ext

c apport de moment cinétique par unité de temps et de masse	 
	 IF(PRESENT(mw_planet))THEN
	  mw_planet(1)=flux*mw_planet0/secon6	!apport direct cgs par giration
	  mw_planet(2)=SIGN(flux,mw_planet0)	!apport fictif par masse ajoutée

c	  PRINT*,'mw_planet0,flux0,m_pla,mzc_ext,flux,mw_planet'
c	  WRITE(*,2000)mw_planet0,flux0,m_pla,mzc_ext,flux,mw_planet


	 ENDIF	 

c corrections à la composition chimique
	 IF(PRESENT(xchim))THEN	 
	  DO i=1,nchim
	   dxchim(i)=dxchim(i)+flux*(x_planet(i)-xchim(i))
	   jac(i,i)=jac(i,i)-flux
	  ENDDO
c	  WRITE(*,2000)age
c	  WRITE(*,2000)dxchim(1:8) ;  WRITE(*,2000)xchim(1:8)	 
c	  PAUSE'planetoides2'	  	  
	 ENDIF
c	 WRITE(*,2000)m_planet,flux0,flux,age,mzc_ext ; PRINT*,profil
	 
	ELSE	!en dehors de l'intervalle de chute
	 IF(PRESENT(xchim))THEN
	  dxchim=0.d0 ; jac=0.d0
	 ENDIF	 
	 IF(PRESENT(m_planet))m_planet=0.d0	 
	 IF(PRESENT(mw_planet))mw_planet=0.d0	 
	ENDIF
	
c	PRINT*,'dxchim' ; WRITE(*,2000)dxchim
c	WRITE(*,2000)age,m_planet	
c	PAUSE'planet'
	
	RETURN
	
	END SUBROUTINE planetoides
		
