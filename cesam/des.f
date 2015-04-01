		
c********************************************************************	

	SUBROUTINE des(fin,dt,teff)
	
c	routine private du module mod_cesam	

c       routine générique du dessin des variables en cours d'évolution

c       Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c       CESAM2k
C
c       25 08 97 : mise en place des variables eulériennes
C
c entrées, significations évidentes
c       p, t, m, l, r, ro, grad_ad, grad_mj,alfa,delta,kap,cp,teff,age
c       n: nombre de couches
c       chim, mc, mct, nc, knotc: pour interpolation de la comp. chim.
c       lim: nombre de limites ZR/ZC
c       lconv: .TRUE. si debut de ZC
c       m_zc: masse aux limites ZC/ZR
c       mstar: masse totale au temps du dessin

c------------------------------------------------------------     

	USE mod_donnees, ONLY : device, dl, dh, dfesh_des, dl_des,
	1 dteff_des, fesh_des, h, langue,ld, l_des, logteff_max,
	2 logteff_min, logl_max, logl_min, nom_des, nom_fich2, teff_des,
	3 xleft, ybot, zoom_l, zoom_t

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt, teff
	LOGICAL, INTENT(in) :: fin

	CHARACTER (len=50) :: chain		
	LOGICAL, SAVE :: pass=.TRUE.
	LOGICAL :: ok1, ok2
	
	NAMELIST/nl_des/teff_des,dteff_des,zoom_t,l_des,dl_des,zoom_l,
	1 fesh_des,dfesh_des,logteff_max,logteff_min,logl_max,logl_min
	NAMELIST/nl_device/h,dh,ld,dl,xleft,ybot,device
	
c------------------------------------------------------------

	SELECT CASE(nom_des)	!select 1
	CASE('des_r' , 'des_m', 'zoom')
	 IF(pass)THEN	
	  chain='device'
	  INQUIRE(file=TRIM(chain),exist=ok1)	 
	  IF(ok1)THEN
	   OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1  file=TRIM(chain))
	   READ(3,nl_device) ; WRITE(*,nl_device)
	   CLOSE(unit=3)	   
	  ENDIF !ok1
	 ENDIF	!pass 	
	 SELECT CASE(nom_des)	!select 2
	 CASE ('des_r')		!select 2
	  CALL des_r(fin,dt,teff)
	 CASE ('des_m') 	!select 2
	  CALL des_m(fin,dt,teff)
	 CASE('zoom')		!select 2
	  IF(pass)THEN
	   chain=TRIM(nom_fich2)//'.zoom'	 
	   INQUIRE(file=TRIM(chain),exist=ok1)	 
	   IF(ok1)THEN !ok1
	    OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1   file=TRIM(chain))
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1001)TRIM(chain) ; WRITE(2,1001)TRIM(chain)
1001	     FORMAT(/,'data for the zoom read in the file : ',a,/)	
	    CASE DEFAULT	
	     WRITE(*,1)TRIM(chain) ; WRITE(2,1)TRIM(chain)
1	     FORMAT(/,'données pour le zoom lues dans le fichier : ',a,/)
	    END SELECT
	    READ(3,nl_des) ; WRITE(*,nl_des) ; WRITE(2,nl_des)
	    CLOSE(unit=3)	
	   ELSE !ok1
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1008)TRIM(chain) ; WRITE(2,1008)TRIM(chain)
1008	     FORMAT('The file zoom for the model : ',a,/,
	1    'is unknown, CESAM looks for a file named zoom')
	    CASE DEFAULT
	     WRITE(*,8)TRIM(chain) ; WRITE(2,8)TRIM(chain)
8	     FORMAT('le fichier zoom pour le modèle : ',a,/,
	1    'non trouvé, recherche et utilisation du fichier zoom')
	    END SELECT	
	    chain='zoom' ; INQUIRE(file=TRIM(chain),exist=ok2)
	    IF(ok2)THEN
	     OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1    file=TRIM(chain))
	     READ(3,nl_des) ; WRITE(*,nl_des) ; WRITE(2,nl_des)
	     CLOSE(unit=3)
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1001)TRIM(chain) ; WRITE(2,1001)TRIM(chain)	
	     CASE DEFAULT	
	      WRITE(*,1)TRIM(chain) ; WRITE(2,1)TRIM(chain)	
	     END SELECT		
	    ELSE	!ok2
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1002) ; WRITE(2,1002)
1002	      FORMAT('No zoom, the file named zoom is unknown')
	     CASE DEFAULT	    	  
	      WRITE(*,2) ; WRITE(2,2)
2	      FORMAT('pas de fichier zoom dans le directory : pas de zoom')
	     END SELECT
	    ENDIF	!ok2	   
	   ENDIF 	!ok1
	  ENDIF 	!pass
	  CALL des_m(fin,dt,teff)
	 END SELECT		!select 2
	CASE ('no_des')		!select 1
	 IF(pass)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1003) ; WRITE(2,1003)	    
1003	   FORMAT('No plot on line')	    
	  CASE DEFAULT
	   WRITE(*,3) ; WRITE(2,3)	    
3	   FORMAT('pas de dessin on line')	    
	  END SELECT	    	  
	 ENDIF			!pass
	CASE DEFAULT		!select 1
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1004)nom_des ; WRITE(2,1004)nom_des	    
1004	  FORMAT('STOP, unknown subroutine for plotting : ',a,/,
	1 'known subroutines : des_r, des_m, no_des, zoom')	    
	 CASE DEFAULT
	  WRITE(*,4)nom_des ; WRITE(2,4)nom_des	    
4	  FORMAT('ARRET, routine de dessin inconnue: ',a,/,
	1 'routines connues: des_r, des_m, no_des, zoom')	    
	 END SELECT	    	  
	 STOP
	END SELECT		!select 1
	pass=.FALSE.	

	RETURN

	END SUBROUTINE des
 
