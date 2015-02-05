
c***********************************************************************

	SUBROUTINE ini_ctes
	
c       routine générique d'initialisation des constantes physiques   
c       Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k
	
c------------------------------------------------------------------
		
	IMPLICIT NONE
	
c-----------------------------------------------------------------

	SELECT CASE(nom_ctes)
	CASE ('ctes_85')    
	 CALL ctes_85    
	CASE ('ctes_94')
	 CALL ctes_94    
	CASE default
	 PRINT*,'routine de constantes inconnue: ',nom_ctes
	 PRINT*,'routines connues: ctes_94, ctes_85'
	 STOP "arrêt car problème dans ini_ctes"
	END SELECT
	
	RETURN

	END SUBROUTINE ini_ctes
