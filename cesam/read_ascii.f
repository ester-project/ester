		
c*****************************************************************

	SUBROUTINE read_ascii(nom_fich,itot,nglob,nvar,abid)

c	routine public du module mod_exploit	
c	lecture des fichiers d'oscillation de cesam2k

c	le fichier débute par 4 lignes de commentaires pour faciliter son
c	identification

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrée:
c	nom_fich: nom du fichier d'oscillation sans l'extension .osc

c sorties:
c	itot: nombre de points
c	nglob: nombre de constantes global
c	nvar: nb. de variables du tableau var sans les éléments chimiques
c	abid(4): chaine de caractères a80

c	glob: variables globales
c		glob(1)=mstar*msol
c		glob(2)=rtot*rsol
c		glob(3)=ltot*lsol
c		glob(4)=z0
c		glob(5)=x0
c		glob(6)=alpha
c		glob(7)=9.d0/4.d0 paramètre de convection arbitraire
c		glob(8)=1./162. paramètre de convection arbitraire
c		glob(9)=X dans ZC
c		glob(10)=Y dans ZC
c		glob(11)=d2p
c		glob(12)=d2ro
c		glob(13)=age
c		glob(14)=wrot vitesse de rotation globale
c		glob(15)=w_rot vitesse de rotation initiale

c	var: variables

c	nvar=22 pour oscillations adiabatiques
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
c		var(21,i)=Ptot / Pgas
c		var(22,i)=gradient radiatif

c	nvar=25 pour inversion

c	        var(23,i)=d Gamma1 / d log P
c	        var(24,i)=d Gamma1 / d log T
c	        var(25,i)=d Gamma1 / dY = d Gamma1 / dZ

c	nvar=44 pour oscillations non adiabatiques

c		var(26,i)=dP / dro (TX)
c		var(27,i)=dP / dT (roX)
c		var(28,i)=dP / dX (Tro)
c		var(29,i)=du / dro (TX)
c		var(30,i)=du / dT (roX)
c		var(31,i)=du / dx(Tro)
c		var(32,i)=énergie interne
c		var(33,i)=d^2P / dro^2 (TX)
c		var(34,i)=d^2P / dro dT (X)	 
c		var(35,i)=d^2P / dT^2(roX)		 	 
c		var(36,i)=d^2U / dro^2 (TX)	 
c		var(37,i)=d^2U / dro dT (X) 
c		var(38,i)=d^2U / dT^2 (X) 	  
c		var(39,i)=dK / dX
c		var(40,i)=d^2K / dT^2	
c		var(41,i)=d epsi / dX
c		var(42,i)=dX / dR
c		var(43,i)=J-B	  
c		var(44,i)=Edding. facteur

c	vecteur de composition chimique

c	  	var(nvar+j,i)=xchim(j)*nucleo(j), j=1,nchim

c----------------------------------------------------------------

	USE mod_donnees, ONLY : Krot, nchim, nom_elem
	USE mod_kind
	USE mod_numerique, ONLY : pause

	IMPLICIT NONE

	CHARACTER (len=*), INTENT(in) :: nom_fich  

	INTEGER, INTENT(out) :: nglob, itot, nvar
	CHARACTER (len=*), INTENT(out), DIMENSION(:) :: abid
		
	INTEGER :: j

	LOGICAL :: ok
	
c------------------------------------------------------------------	

	INQUIRE(file=TRIM(nom_fich),exist=ok)
	IF(ok)THEN
	 OPEN(unit=30,form='formatted',status='old',file=TRIM(nom_fich))

	 READ(30,'(a)')abid(1)	!m 100X707a173_$.dat du 18/ 8/1993 etc...
c	 PRINT*,abid(1)
	 READ(30,'(a)')abid(2)	!fichier pour calcul des oscillations etc..
c	 PRINT*,abid(2)
	 READ(30,'(a)')abid(3)	!methode: CESAM 3.0.0.0 colloc. etc...
c	 PRINT*,abid(3)
	 READ(30,'(a)')abid(4)
c	 PRINT*,abid(4)
	 READ(30,89)nchim	!  9  H  He3 He4 etc...
89	 FORMAT(i3)
	 ALLOCATE(nom_elem(nchim))
	 BACKSPACE (unit=30)
	 READ(30,90)nchim,nom_elem(1:nchim)
90	 FORMAT(i3,14(1x,a4))
	 READ(30,137)itot,nglob,nvar,nchim,Krot	 
137	 FORMAT(5i10)
c	 PRINT*,itot,nglob,nvar,nchim,Krot ; PAUSE'itot,nglob,nvar,nchim,Krot'
	 ALLOCATE(glob(nglob),var(nvar+nchim,itot))	 
	 READ(30,138)glob	 	 
138	 FORMAT(5es19.12)
c	 PRINT*,nglob ; WRITE(*,138)glob ; PAUSE'glob'
	 DO j=1,itot
	  READ(30,138)var(:,j)
c	  PRINT*,nvar,j ; WRITE(*,138)var(:,j) ; PRINT* 
c	  WRITE(*,138)var(nvar+1:nvar+nchim,j) ; CALL pause('var') 
	 ENDDO
	 CLOSE(unit=30)
	ELSE
	 WRITE(*,1)nom_fich
1	 FORMAT('fichier inconnu: ',a80)
	 PRINT*,'ARRET' ; STOP
	ENDIF 

	RETURN

	END SUBROUTINE read_ascii		
