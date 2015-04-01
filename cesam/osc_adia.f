
c*******************************************************************

	SUBROUTINE osc_adia(var,glob,itot,ivar)
	
c	routine private du module mod_cesam
c	création du fichier pour le calcul des oscillations adiabatiques
c	le fichier créé aura pour extension -ad.osc
c	Exemple: mon_modele-ad.osc

c	on entre dans var(P, T, L...) et glob(age, d2ro...) toutes les
c	quantités dont on extrait dans les tableaux eglob(nglob) et
c	evar(nvar,itot) celles à écrire suivant l'ordre défini dans
c	les tableaux d'indices tglob(nglob) et tvar(nvar)
c	ces tableaux sont ALLOCATABLES de façon à pouvoir facilement
c	en changer les dimensions et les contenus

c	glob: variables globales du fichier mon_modele-ad.osc
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
c	        glob(14)=g constante de la gravitation
c		glob(15)=msol
c		glob(16)=rsol
c		glob(17)=lsol

c	var: variables locales utilisées

c	ivar=22 pour oscillations adiabatiques
c		var(1,i)=r*rsol
c	 	var(2,i)=log(m/mstar) -1.d38 au centre
c		var(3,i)=t
c		var(4,i)=Ptot
c		var(5,i)=ro
c		var(6,i)=gradient réel d ln T / d ln P
c		var(7,i)=l*lsol
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

c	  	var(22+j,i)=xchim(j)*nucleo(j), j=1,nchim

!       var(47,i)=nuclear energy ; var(48,i)=gravitational energy !YLD

c	P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY: all_output, Krot, methode, nchim, nom_abon,
	1 nom_atm, nom_conv, nom_ctes, nom_diffm, nom_difft, nom_diffw,
	2 nom_elem, nom_etat, nom_fich2, nom_nuc, nom_nuc_cpl,
	3    nom_opa, nom_pertm, nom_pertw, nom_tdetau, nom_thw, nom_output ! YLD add nom_output
!	3 nom_opa, nom_pertm, nom_pertw, nom_tdetau, nom_thw
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
	CHARACTER (len=100) :: titre_osc, titre_losc1, titre_losc2
	CHARACTER (len=80) :: fich, fich_losc
	CHARACTER (len=95) :: chain, chain1, chain2
	
c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	titre_osc='Fichier pour oscillations adiabatiques:'
	titre_losc1=' 1:r/R, 2:m/M, 3:t, 4:p, 5:rho, 6:L, 7: kappa, 8: eps_tot'
	titre_losc2='9:Xc 10:Yc, 11:nablaT, 12:nabla_ad; 13:nabla_rad, 14:Gamma1, 15: A'
!	nglob=17 ; nvar=22 ; totvar=nvar+nchim !YLD
    
      SELECT CASE(nom_output) !YLD
	  CASE ('osc_adia_plus')  !YLD
	     nglob=17 ; nvar=24 ; totvar=nvar+nchim  !YLD
      CASE default !YLD
	      nglob=17 ; nvar=22 ; totvar=nvar+nchim  !YLD
      END SELECT !YLD

	ALLOCATE(eglob(nglob),evar(totvar,itot),tglob(nglob),tvar(totvar))
	
	tglob=(/ 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 /)
!	tvar=(/ (i,i=1,5),(i,i=7,15),(i,i=17,24),(i,i=ivar+1,ivar+nchim) /)

      SELECT CASE(nom_output) !YLD
	  CASE ('osc_adia_plus')  !YLD
	     tvar=(/ (i,i=1,5),(i,i=7,15),(i,i=17,24),(i,i=47,48), (i,i=ivar+1,ivar+nchim) /) !YLD
      CASE default !YLD
	    tvar=(/ (i,i=1,5),(i,i=7,15), (i,i=17,24),(i,i=ivar+1,ivar+nchim) /) !YLD
      END SELECT !YLD

c transposition dans les tableaux d'écriture
	eglob(1:nglob)=glob(tglob(1:nglob))
	evar(1:totvar,:)=var(tvar(1:totvar),:)
		
c écritures
	IF(all_output)THEN
	 WRITE(number,10)model_num
10	 FORMAT(i4.4)
	 fich=number//'-'//TRIM(nom_fich2)//'-ad.osc'
	 fich_losc=number//'-'//TRIM(nom_fich2)//'.losc'
	ELSE
	 fich=TRIM(nom_fich2)//'-ad.osc'
	 fich_losc=TRIM(nom_fich2)//'.losc'
	ENDIF
	WRITE(*,1)ADJUSTL(fich) ; WRITE(2,1)ADJUSTL(fich)
1	FORMAT(/,'création du fichier ascii: ',a)
	OPEN(unit=30,form='formatted',status='unknown',file=TRIM(fich))	
	OPEN(unit=40,form='formatted',status='unknown',file=TRIM(fich_losc))
	chain=TRIM(titre_osc)//' '//ADJUSTL(fich)
	chain1='#'//TRIM(titre_losc1)
	chain2='#'//TRIM(titre_losc2)
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
        WRITE(40,2)ADJUSTL(chain1)
        WRITE(40,2)ADJUSTL(chain2)
   	DO j=itot,1,-1
         WRITE(40,14) evar(1,j)/glob(2),exp(evar(2,j)),evar(3,j),
	1 evar(4,j),evar(5,j),evar(7,j), evar(8,j),
	2 evar(9,j),evar(23,j),evar(25,j), evar(6,j),
	4 evar(11,j),evar(22,j),evar(10,j), evar(15,j) 
14	FORMAT(15(1x,es19.12))
	ENDDO
	CLOSE(unit=40)

	DEALLOCATE(eglob,evar,tglob,tvar)
	
	RETURN	
	
	END SUBROUTINE osc_adia
