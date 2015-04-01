
c*******************************************************************

	SUBROUTINE ascii(nglob,tglob,totvar,tvar,var,glob,itot,
	1 logm,reverse,titre_ascii)
	
c	routine private du module mod_cesam
c	création d'un fichier de sortie ASCII personnalisé
c	les quantités à écrire et leur ordre sont définies dans le fichier
c	sortie_ascii
c	le fichier créé aura pour extension -ascii
c	Exemple: mon_modele-ascii

c	on entre dans var(P, T, L...) et glob(age,d2ro...) toutes les
c	quantités dont on extrait dans les tableaux eglob(nglob) et
c	evar(nvar,itot) celles à écrire suivant l'ordre définit dans
c	les tableaux d'indices tglob(nglob) et tvar(nvar)
c	ces tableaux sont ALLOCATABLES de façon à pouvoir facilement
c	en changer les dimensions et les contenus

c	glob: variables globales

c	glob(1)=mstar*msol
c	glob(2)=rtot*rsol
c	glob(3)=ltot*lsol
c	glob(4)=z0
c	glob(5)=x0
c	glob(6)=alpha
c	glob(7)=9./4.
c	glob(8)=1./162.
c	glob(9)=X dans ZC
c	glob(10)=Y dans ZC
c	glob(11)=d2p
c	glob(12)=d2ro
c	glob(13)=age
c	glob(14)=wrot vitesse de rotation globale
c	glob(15)=w_rot initial

c	nglob=15

c	var: variables

c	var(1,i): Rayon
c	var(2,i): Ln M/Mtot
c	var(3,i): Température
c	var(4,i): Pression
c	var(5,i): Densité
c	var(6,i): Gradient
c	var(7,i): Gradient
c	var(8,i): Luminosité
c	var(9,i): Opacité
c	var(10,i): Energie nuc+grav
c	var(11,i): Grand Gamma1
c	var(12,i): Gradient adiabatique
c	var(13,i): Delta
c	var(14,i): Cp
c	var(15,i): Mue^(-1)
c	var(16,i): Mu      
c	var(17,i): Vaissala
c	var(18,i): Omega
c	var(19,i): dln kappa/dln T
c	var(20,i): dln kappa/dln ro
c	var(21,i): d epsilon(nuc) / d ln T
c	var(22,i): d epsilon(nuc) / d ln ro
c	var(23,i): Ptot/Pgaz
c	var(24,i): Gradient radiatif       
c	var(25,i): d Gamma1 / d lnP (TY)
c	var(26,i): d Gamma1 / d lnT (PY)
c	var(27,i): d Gamma1 / dY	(PT)
c	var(28,i): dP / dro (TX)
c	var(29,i): dP / dT (roX)
c	var(30,i): dP / dX (Tro)
c	var(31,i): du / dro (TX)
c	var(32,i): du / dT (roX)
c	var(33,i): du / dX(Tro)
c	var(34,i): énergie interne
c	var(35,i): d^2P / dro^2 (TX)
c	var(36,i): d^2P / dro dT (X)	 
c	var(37,i): d^2P / dT^2(roX)		 	 
c	var(38,i): d^2U / dro^2 (TX)	 
c	var(39,i): d^2U / dro dT (X) 
c	var(40,i): d^2U / dT^2 (X) 	  
c	var(41,i): dK / dX
c	var(42,i): d^2K / dT^2	
c	var(43,i): d epsi / dX
c	var(44,i): dX / dR
c	var(45,i): J-B	  
c	var(46,i): Edding. facteur
c	var(ivar+j,i): xchim1g(j=1,nchim) Abondances / gramme

c	ivar=46 ! YLD ivar=48 (ajout eps_grav, eps_nuc) 

c	P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY: all_output, Krot, methode, nchim, nom_abon,
	1 nom_atm, nom_conv, nom_ctes, nom_diffm, nom_difft, nom_diffw,
	2 nom_elem, nom_etat, nom_fich2, nom_nuc, nom_nuc_cpl,
	3 nom_opa, nom_pertm, nom_pertw, nom_tdetau
	USE mod_kind
	USE mod_variables, ONLY : model_num
      
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: var	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: glob
	INTEGER, INTENT(in), DIMENSION(:) :: tglob, tvar	
	INTEGER, INTENT(in) :: itot, nglob, totvar
	LOGICAL, INTENT(in) :: logm, reverse
	CHARACTER (len=*), INTENT(in) :: titre_ascii	
		
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: evar
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: eglob
	
	INTEGER, DIMENSION(8) :: values
	INTEGER :: i
	
	CHARACTER (len=4) :: number	
	CHARACTER (len=5) :: zone	
	CHARACTER (len=8) :: date	
	CHARACTER (len=9), PARAMETER, DIMENSION(12) :: mois=(/
	1 'Janvier  ','Février  ','Mars     ','Avril    ',
	2 'Mai      ','Juin     ','Juillet  ','Aout     ',
	3 'Septembre','Octobre  ','Novembre ','Décembre ' /)
	CHARACTER (len=10) :: time	
	
	CHARACTER (len=80) :: fich
	CHARACTER (len=95) :: chain
		
c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	ALLOCATE(eglob(nglob),evar(totvar,itot))

c transposition dans les tableaux d'écriture
	eglob(1:nglob)=glob(tglob(1:nglob))
	evar(1:totvar,:)=var(tvar(1:totvar),:)
	
c la fraction de masse en réel et non en dex, 0 au centre
	IF(.NOT.logm)THEN
	 evar(2,:)=EXP(evar(2,:)) ; evar(2,1)=0.d0
	ENDIF	

c tabulation du centre vers la surface
	IF(reverse)THEN
	 DO i=1,totvar
	  evar(i,1:itot:1)=evar(i,itot:1:-1)
	 ENDDO
	ENDIF
		
c écritures
	IF(all_output)THEN
	 WRITE(number,10)model_num
10	 FORMAT(i4.4)
	 fich=TRIM(nom_fich2)//number//'-ascii'	 	
	ELSE
	 fich=TRIM(nom_fich2)//'-ascii'	
	ENDIF
	WRITE(*,1)ADJUSTL(fich) ; WRITE(2,1)ADJUSTL(fich)
1	FORMAT(/,'création du fichier ascii: ',a)
	OPEN(unit=30,form='formatted',status='unknown',file=TRIM(fich))
	
	chain=TRIM(titre_ascii)//' '//ADJUSTL(fich)
	WRITE(30,2)ADJUSTL(chain)
2	FORMAT(a95)
	CALL date_and_time(date,time,zone,values)	
	chain=TRIM(methode)//', '//date(7:8)//' '
	1 //TRIM(mois(values(2)))//' '//date(1:4)//' '//time(1:2)//
	2 'h'//time(3:4)
	WRITE(30,2)ADJUSTL(chain)
	chain='Physique utilisée: '//TRIM(nom_etat)//', '//TRIM(nom_opa)
	1 //', '//TRIM(nom_conv)//', '//TRIM(nom_nuc)
	2 //', '//TRIM(nom_nuc_cpl)	
	WRITE(30,2)ADJUSTL(chain)
	chain='   '//TRIM(nom_abon)//', '//TRIM(nom_atm)
	1 //', '//TRIM(nom_tdetau)//', '//TRIM(nom_pertm)
	2 //', '//TRIM(nom_pertw)//', '//TRIM(nom_diffm)
	3 //', '//TRIM(nom_difft)//', '//TRIM(nom_diffw)
	4 //', '//TRIM(nom_ctes)
	WRITE(30,2)ADJUSTL(chain)
	WRITE(30,90)nchim,(nom_elem(i),i=1,nchim)
90	FORMAT(i3,14(1x,a4))
	
	WRITE(30,137)itot,nglob,totvar,nchim,Krot
137	FORMAT(4i10)
	WRITE(30,138)eglob
138	FORMAT(5es19.12)
	DO i=1,itot
	 WRITE(30,138)evar(1:totvar,i)
	ENDDO
	CLOSE(unit=30)

	DEALLOCATE(eglob,evar)
	
	RETURN	
	
	END SUBROUTINE ascii
