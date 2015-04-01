
c***************************************************************

	SUBROUTINE lit_hr(init,chaine,fin,erreur,log_l,log_r,log_teff,
	1 vrot,wrot)

c subroutine public du module mod_exploit	
c lecture du fichier mon_modele.HR	

c entrées
c	init=.TRUE. : initialisation
c	chaine : nom_générique du modèle

c sorties :
c	fin=.TRUE. : le fichier .HR a été lu en entier
c	erreur=.TRUE. : erreur de lecture 
c	log_l,log_r,log_teff : sens évident
c	vrot : vitesse couche externe en km/s
c	wrot : vitesse angulaire de la couche externe en radian/s

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c----------------------------------------------------------------------

	USE mod_donnees, ONLY : nchim, nom_elem
	USE mod_kind
	USE mod_variables, ONLY : age, lconv, lim, model_num, mstar, m_zc,
	1 r_zc, r_ov
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(out) :: log_l, log_r, log_teff, vrot, wrot	
	INTEGER :: i, Krot

	LOGICAL, INTENT(inout) :: init
	LOGICAL, INTENT(out) :: erreur, fin
	
	CHARACTER (len=4) :: nom_vwrot
	CHARACTER (len=*), INTENT(in) :: chaine
	
c---------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

c initialisation
	IF(init)THEN
	 OPEN(unit=25,form='FORMATTED',status='old',file=chaine)
	 READ(25,1,END=30,ERR=4)age,nchim,Krot,lim,model_num,lconv(1:lim)
1	 FORMAT(es22.15,3i3,i4,20l2)
	 REWIND(unit=25) ; init=.FALSE.
	 IF(ALLOCATED(compc))DEALLOCATE(compc,compe,nom_elem)
	 ALLOCATE(compc(nchim),compe(nchim),nom_elem(nchim))
	ENDIF
	READ(25,1,END=30,ERR=4)age,nchim,Krot,lim,model_num,(lconv(i),i=1,lim)
	READ(25,2,END=30,ERR=4)log_teff,log_l,log_r,mstar,
	1(m_zc(i),r_zc(i),r_ov(i),i=1,lim)
2	FORMAT(6es13.6)	
	READ(25,3,END=30,ERR=4)(nom_elem(i),compc(i),compe(i),i=1,nchim)
3	FORMAT(a4,2es12.5)
	IF(Krot >= 3)THEN
	 READ(25,3,END=30,ERR=4)nom_vwrot,vrot,wrot
	ELSE
	 vrot=0.d0 ; wrot=0.d0
	ENDIF
	fin=.FALSE. ; erreur=.FALSE.; RETURN	
30	fin=.TRUE.  ; erreur=.FALSE.; CLOSE(unit=25) ; RETURN
4	fin=.FALSE. ; erreur=.TRUE. ; PRINT*,'Dans lit_hr, erreur=',erreur
	WRITE(*,5)model_num,age,log_teff,log_l
5	FORMAT('model_num=',i5,', age=',es10.3,', log_teff=',es10.3,',log_l=',
	1 es10.3)	
	RETURN

c	WRITE(*,1)age,nchim,lim,(lconv(i),i=1,lim)	
c	WRITE(*,2)log_teff,log_l,log_r,mstar,
c	1 (m_zc(i),r_zc(i),r_ov(i),i=1,lim)
c	WRITE(*,3)(nom_elem(i),compc(i),compe(i),i=1,nchim)	
c	PAUSE
c	CLOSE(unit=25) ; RETURN
	
	END SUBROUTINE lit_hr
