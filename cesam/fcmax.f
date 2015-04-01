	
	SUBROUTINE fcmax(clogic)	
	
c subroutine subordonnée de resout du module mod_static
c clogic = .TRUE. s'il faut utiliser le nombre maximum de couches 

c ATTENTION ne pas mettre en 1-ère ligne c**********************
c c'est refusé par le compilateur

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------

	LOGICAL, INTENT(out) :: clogic
	LOGICAL :: logic

c---------------------------------------------------------------------
	
	logic=age+dt+0.1d0 >= agemax .OR. lt_stop .OR. lhe_stop
	1 .OR. lx_stop .OR. nb_modeles >= nb_max_modeles-1
	
	IF(precision == 'mx')THEN
	 clogic=.TRUE. ; RETURN
	ELSEIF(nom_output == 'no_output')THEN
	 clogic=.FALSE. ; RETURN
	ELSEIF(logic)THEN	
	 SELECT CASE(precision)
	 CASE('sa')
	  precix=5.d-6	!précision de l'intégration globale
	  clogic=.TRUE.	!nombre max de couches
	  loc_zc=precix !précision de la localisation des limites ZR/ZC
	 CASE('co')
	  clogic=.TRUE.	!nombre max de couches
	 CASE DEFAULT
	  clogic=nc_max	!nombre max de couches   
	 END SELECT  
	ELSE
	 clogic=.FALSE. 	  	  
	ENDIF
	
c	PRINT*,age+dt+1.d0,agemax,lt_stop,lhe_stop,lx_stop,nb_modeles,
c	1 nb_max_modeles-1,nom_output	
c	PRINT*,logic ; PRINT*,clogic ; IF(clogic)PAUSE'fcmax'	
  
	RETURN
  
	END SUBROUTINE fcmax
