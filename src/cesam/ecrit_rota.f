
c************************************************************************

	SUBROUTINE ecrit_rota(dt)

c routine public du module mod_evol
c écriture des variables de la rotation et des coefficients pour dessins
c dans mon_modele_coeff_rota.dat, utilisé par des2k_rot

c	'no_des', Kdes_rot=0
c	'end_evol', Kdes_rot=1
c	'all_mod', Kdes_rot=2
c	'all_iter', Kdes_rot=3
c	'end_mod', Kdes_rot=4	 

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c-----------------------------------------------------------------

	USE mod_donnees, ONLY : ihe4, ili7, Kdes_rot, Krot, methode,
	1 nom_fich2, nchim, nom_abon, nom_atm, nom_conv, nom_ctes, nom_des,
	2 nom_diffm, nom_difft, nom_diffw, nom_elem, nom_etat, nom_fich2,
	3 nom_nuc, nom_nuc_cpl, nom_opa, nom_pertm, nom_pertw, nom_tdetau,
	4 nom_thw, nom_thw, nrl, nrot, ord_rot
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss	
	USE mod_variables, ONLY : knotr, model_num, mrot, mrott, n_rot,
	1 rota, sortie	 

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: ecrit	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: frl
	REAL (kind=dp), DIMENSION(nrot,0:1) :: y	
	REAL (kind=dp), DIMENSION(ncoeff+nchim) :: coeff
	REAL (kind=dp), DIMENSION(nrot) :: bd, bs
	
	INTEGER, DIMENSION(8) :: values	
	INTEGER, SAVE :: l=1, nb_var
	INTEGER :: i
	
	LOGICAL, SAVE :: init=.TRUE.

	CHARACTER (len=4) :: number	
	CHARACTER (len=5) :: zone	
	CHARACTER (len=8) :: date	
	CHARACTER (len=9), PARAMETER, DIMENSION(12) :: mois=(/
	1 'Janvier  ','Février  ','Mars     ','Avril    ',
	2 'Mai      ','Juin     ','Juillet  ','Aout     ',
	3 'Septembre','Octobre  ','Novembre ','Décembre ' /)
	CHARACTER (len=10) :: time	
	CHARACTER (len=50) :: titre	
	CHARACTER (len=80) :: fich
	CHARACTER (len=95), SAVE  :: ord_var		
	CHARACTER (len=95) :: chain
 
c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

c initialisation
	IF(init)THEN
	 init=.FALSE.
	 nb_var=ncoeff+nchim
	 ALLOCATE(frl(nrl))	 
	 titre='Fichier pour la diffusion du moment cinétique:'
	 SELECT CASE(Krot)
	 CASE(3)	
	  ord_var='Variables : R, M, Omega, U, Theta, Psi, Lambda, Flux, Deff, Dh, Dv, f_eps, T, ro, grad_mu, Xchim'
	 CASE(4)
	  ord_var='Variables : R, M, Omega, U, Psi, Lambda, Phi, Flux, Deff, Dh, Dv, f_eps, T, ro, grad_mu, Xchim'
	 CASE DEFAULT
	  PRINT*,'ERREUR, dans ecrit_rota Krot=',Krot
	  PRINT*,'Krot=3 ou 4'
	  CALL sortie
	 END SELECT	 
	  
	ENDIF

c allocation
	ALLOCATE(ecrit(nb_var,n_rot))

c construction de la table des quantités à écrire
	DO i=2,n_rot
	
c variables et dérivées
	 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 mrot(i),l,bs,bd)
	 IF(no_croiss)PRINT*,'ecrit-rota, Pb. en 1'
	 y(:,0)=bs ; y(:,1)=bd
	
c les coefficients de rotation 
	 CALL coeff_rota(dt,mrot(i),y,frl,coeff=coeff)
	 ecrit(:,i)=coeff
	ENDDO
	
c écriture de U
c	WRITE(*,2000)ecrit(4,:) ; PAUSE'U'
	
c pour le centre	
	ecrit(:,1)=ecrit(:,2) ; ecrit(1:2,1)=0.d0
	
c écriture dans mon_modele_coeff_rota.dat
	IF(Kdes_rot == 2)THEN
	 WRITE(number,10)model_num
10	 FORMAT(i4.4)
	 fich=number//'-'//TRIM(nom_fich2)//'_coeff_rota.dat'
	ELSE
	 fich=TRIM(nom_fich2)//'_coeff_rota.dat'
	ENDIF

	chain=TRIM(titre)//' '//ADJUSTL(fich)
	OPEN(unit=30,form='formatted',status='unknown',file=TRIM(fich))
	WRITE(30,2)ADJUSTL(chain)
2	FORMAT(a)
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
	4 //', '//TRIM(nom_thw)
	WRITE(30,2)ADJUSTL(chain)
	WRITE(30,2)ADJUSTL(ord_var)
	WRITE(30,90)nchim,(nom_elem(i),i=1,nchim)
90	FORMAT(i3,14(1x,a4))	
	WRITE(30,137)n_rot,nb_var,ncoeff,nchim,Krot,ihe4,model_num,ili7
137	FORMAT(8i10)
	DO i=1,n_rot
	 WRITE(30,3)ecrit(:,i)
3	 FORMAT(5es19.12)
	ENDDO

	CLOSE(unit=30)

c s'il n'y a pas de dessin en ligne, dessin des variables de la rotation
c sur mon_modele.ps ou en ligne en adaptant le device dans plot_rota
	IF(nom_des == 'no_des')CALL plot_rota
	
	DEALLOCATE(ecrit)

	RETURN
	
	CONTAINS
	 INCLUDE 'plot_rota.f'	
	
	END SUBROUTINE ecrit_rota
