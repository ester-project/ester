	
c*****************************************************************	

	SUBROUTINE output(var,glob,itot,ivar)

c	routine private du module mod_cesam

c	routine générique de création des fichiers de sortie en ASCII
c	pour les fichiers d'oscillation coder:
c	nom_output='osc_adia', 'osc_invers', 'osc_nadia'
c	Pour personaliser le fichier ASCII de sortie coder:
c	nom_output='ascii'
c	créer un fichier de nom sortie_ascii dans lequel sera indiqué
c	les variables que l'on désire écrire
c	EXEMPLE:
c	4  5    ==> nglob=12 global et nvar=3 variables
c	1 11 5 2 ==> indices des 4 global : mstar*msol, d2p, x0, rtot*rsol
c	43 1 7 3 8 ==> indices des 5 variables :J-B, r*rsol, l, t, kap  
c	logm=.FALSE. ==> m/mstar en fraction et non en LOG10
c	reverse=.TRUE. ==> fichier construit du centre vers la surface
c	comp_chim=.TRUE. ==> insertion de la composition chimique

c	Fichier pour reprise ascii:

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c---------------------------------------------------------------

	USE mod_donnees, ONLY: nchim, nom_fich2, nom_output
	USE mod_kind
	 	 
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: var	 
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: glob
	INTEGER, INTENT(in) :: itot, ivar
	
	INTEGER, ALLOCATABLE, DIMENSION(:) :: tglob, tvar, tvart
	INTEGER :: i, nglob, nvar, totvar
	
	LOGICAL :: comp_chim, logm, ok, reverse

	CHARACTER (len=80) :: chain, titre_ascii
			 
c----------------------------------------------------------------

	SELECT CASE(nom_output)
!	CASE ('osc_adia','all_adia')
	case ('osc_adia','all_adia', 'osc_adia_plus') !YLD
	 CALL osc_adia(var,glob,itot,ivar)
	CASE ('osc_invers','all_invers')
	 CALL osc_invers(var,glob,itot,ivar)	 
	CASE ('osc_nadia','all_nadia')
	 CALL osc_nadia(var,glob,itot,ivar)
	CASE ('no_output')
	 RETURN
	CASE ('ascii','all_ascii')	
	 chain=TRIM(nom_fich2)//'.ascii'	 
	 INQUIRE(file=TRIM(chain),exist=ok)	 
	 IF(ok)THEN
	  OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 ELSE
	  WRITE(*,8)TRIM(chain) ; WRITE(2,8)TRIM(chain)
8	  FORMAT('fichier de codage de sortie en ASCII : ',a,/,
	1 'non trouvé, recherche et utilisation du fichier sortie_ascii')
	  chain='sortie_ascii'	
	  INQUIRE(file=TRIM(chain),exist=ok)
	  IF(ok)THEN
	   OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1  file=TRIM(chain))
	   WRITE(*,7)TRIM(chain) ; WRITE(2,7)TRIM(chain)
7	   FORMAT('indices des global et variables du fichier: ',a)
	   READ(3,*)nglob,nvar
	   ALLOCATE(tglob(nglob),tvar(nvar))
	   READ(3,*)tglob,tvar,logm,reverse,comp_chim
	   READ(3,'(a)')titre_ascii	  	  
	   CLOSE(unit=3)
	
c on ajoute les indices pour la composition chimique
	   totvar=nvar	
	   IF(comp_chim)THEN
	    totvar=totvar+nchim
	    ALLOCATE(tvart(totvar))
	    tvart(1:nvar)=tvar
	    tvart(nvar+1:totvar)=(/ (i,i=ivar+1,ivar+nchim) /)
	    DEALLOCATE(tvar) ; ALLOCATE(tvar(totvar)) ; tvar=tvart
	    DEALLOCATE(tvart)
	   ENDIF
	  
	   CALL ascii(nglob,tglob,totvar,tvar,var,glob,itot,
	1  logm,reverse,titre_ascii)	  
	   DEALLOCATE(tglob,tvar)	  
	  ELSE	  
	   WRITE(*,9) ; WRITE(2,9) ; nom_output='no_output' ; RETURN
9	   FORMAT('fichier sortie_ascii non trouvé, no sortie ASCII')
	  ENDIF
	 ENDIF
	CASE default
	 chain=TRIM(nom_fich2)//'.ascii'
	 WRITE(*,10)nom_output,chain ; WRITE(2,10)nom_output,chain
	 nom_output='no_output' ; RETURN
10	 FORMAT('pas de sortie ASCII, nom de fichier ASCII inconnu: ',
	1 a,/,'noms connus: osc_adia, osc_invers, osc_nadia, ascii',/,
	2 'ou encore :',a)
	END SELECT
	
	RETURN

	END SUBROUTINE output
