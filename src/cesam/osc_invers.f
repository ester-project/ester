
c*******************************************************************

	SUBROUTINE osc_invers(var,glob,itot,ivar)
	
c	routine private du module mod_cesam
c	création du fichier pour les inversions
c	le fichier créé aura pour extension -inv.osc
c	Exemple: mon_modele-inv.osc

c	on entre dans var(P, T, L...) et glob(age,d2ro...) toutes les
c	quantités dont on extrait dans les tableaux eglob(nglob) et
c	evar(nvar,itot) celles à écrire suivant l'ordre défini dans
c	les tableaux d'indices tglob(nglob) et tvar(nvar)
c	ces tableaux sont ALLOCATABLES de façon à pouvoir facilement
c	en changer les dimensions et les contenus

c	glob: variables globales du fichier mon_modele-inv.osc
c		glob(1)=mstar*msol
c		glob(2)=rtot*rsol
c		glob(3)=ltot*lsol
c		glob(4)=z0
c		glob(5)=x0
c		glob(6)=alpha
c		glob(7)=X dans ZC
c		glob(8)=Y dans ZC
c		glob(9)=d2p
c		glob(10)=d2ro
c		glob(11)=age
c		glob(12)=wrot vitesse de rotation globale
c		glob(13)=w_rot initial
c		glob(14)=g constante de la gravitaion
c		glob(15)=msol
c		glob(16)=rsol
c		glob(17)=lsol

c	var: variables locales utilisées

c	nvar=25 pour inversion
c		var(1,i)=r*rsol
c	 	var(2,i)=log(m/mstar) -1.d38 au centre
c		var(3,i)=t
c		var(4,i)=Ptot
c		var(5,i)=ro
c		var(6,i)=gradient reel d ln T / d ln P
c		var(7,i)=l
c		var(8,i)=kap
c		var(9,i)=énergie thermo+gravifique
c		var(10,i)=grand Gamma1
c		var(11,i)=gradient adiabatique
c		var(12,i)=delta
c		var(13,i)=cp
c		var(14,i)=mu elec.
c		var(15,i)=vaissala, 0 au centre
c	 	var(16,i)=vitesse angulaire, radian/sec
c	 	var(17,i)=d ln kappa / d ln T
c	 	var(18,i)=d ln kappa / d ln ro
c	 	var(19,i)=d epsilon(nuc) / d ln T
c	 	var(20,i)=d epsilon(nuc) / d ln ro
c		var(21,i)=Ptot / Pgas ou grad_mu
c		var(22,i)=gradient radiatif
c	        var(23,i)=d Gamma1 / d log P
c	        var(24,i)=d Gamma1 / d log T
c	        var(25,i)=d Gamma1 / dY = d Gamma1 / dZ
c	  	var(25+j,i)=xchim(j)*nucleo(j), j=1,nchim

c	P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY: all_output, Krot, methode, nchim, nom_abon,
	1 nom_atm, nom_conv, nom_ctes, nom_diffm, nom_difft, nom_diffw,
	2 nom_elem, nom_etat, nom_fich2, nom_nuc, nom_nuc_cpl,
	3 nom_opa, nom_pertm, nom_pertw, nom_tdetau, nom_thw
	USE mod_kind
	USE mod_variables, ONLY : model_num
      
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: var	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: glob
	INTEGER, INTENT(in) :: itot, ivar
		
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: evar
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: eglob

	INTEGER, ALLOCATABLE, DIMENSION(:) :: tglob, tvar
		
	INTEGER, DIMENSION(8) :: values
	INTEGER :: i, j, nglob, nvar, totvar
	
	CHARACTER (len=4) :: number	
	CHARACTER (len=5) :: zone	
	CHARACTER (len=8) :: date	
	CHARACTER (len=9), PARAMETER, DIMENSION(12) :: mois=(/
	1 'Janvier  ','Février  ','Mars     ','Avril    ',
	2 'Mai      ','Juin     ','Juillet  ','Aout     ',
	3 'Septembre','Octobre  ','Novembre ','Décembre ' /)
	CHARACTER (len=10) :: time	
	CHARACTER (len=50) :: titre_osc
	CHARACTER (len=80) :: fich
	CHARACTER (len=95) :: chain	
	
c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	titre_osc='Fichier pour inversion:'
	nglob=17 ; nvar=25 ; totvar=nvar+nchim
	ALLOCATE(eglob(nglob),evar(totvar,itot),tglob(nglob),tvar(totvar))
	
	tglob=(/ 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 /)
	tvar=(/ (i,i=1,5),(i,i=7,15),(i,i=17,27),(i,i=ivar+1,ivar+nchim) /)	

c transposition dans les tableaux d'écriture
	eglob(1:nglob)=glob(tglob(1:nglob))
	evar(1:totvar,:)=var(tvar(1:totvar),:)
		
c écritures
	IF(all_output)THEN
	 WRITE(number,10)model_num
10	 FORMAT(i4.4)
	 fich=number//'-'//TRIM(nom_fich2)//'-inv.osc'	 	
	ELSE
	 fich=TRIM(nom_fich2)//'-inv.osc'		
	ENDIF
	WRITE(*,1)ADJUSTL(fich) ; WRITE(2,1)ADJUSTL(fich)
1	FORMAT(/,'création du fichier ascii: ',a)
	OPEN(unit=30,form='formatted',status='unknown',file=TRIM(fich))
		
	chain=TRIM(titre_osc)//' '//ADJUSTL(fich)
	WRITE(30,2)ADJUSTL(chain)
2	FORMAT(a95)
	CALL date_and_time(date,time,zone,values)	
	chain=TRIM(methode)//', '//date(7:8)//' '
	1 //TRIM(mois(values(2)))//' '//date(1:4)//' '//time(1:2)//
	2 'h'//time(3:4)
	WRITE(30,2)ADJUSTL(chain)
	chain='Physique utilisée: '//TRIM(nom_etat)//', '//TRIM(nom_opa)
	1 //', '//TRIM(nom_conv)//', '//TRIM(nom_nuc)//', '//TRIM(nom_ctes)
	2 //', '//TRIM(nom_nuc_cpl)//', '//TRIM(nom_thw)		
	WRITE(30,2)ADJUSTL(chain)
	chain='   '//TRIM(nom_abon)//', '//TRIM(nom_atm)
	1 //', '//TRIM(nom_tdetau)//', '//TRIM(nom_pertm)
	2 //', '//TRIM(nom_pertw)//', '//TRIM(nom_diffm)
	3 //', '//TRIM(nom_difft)//', '//TRIM(nom_diffw)
	WRITE(30,2)ADJUSTL(chain)	
	WRITE(30,90)nchim,(nom_elem(i),i=1,nchim)
90	FORMAT(i3,14(1x,a4))	
	WRITE(30,137)itot,nglob,nvar,nchim,Krot
137	FORMAT(5i10)
	WRITE(30,138)eglob
138	FORMAT(5es19.12)
	DO j=1,itot
	 WRITE(30,138)evar(1:totvar,j)
	ENDDO
	CLOSE(unit=30)

	DEALLOCATE(eglob,evar,tglob,tvar)
	
	RETURN	
	
	END SUBROUTINE osc_invers
