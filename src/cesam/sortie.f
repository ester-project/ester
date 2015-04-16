
c********************************************************

	SUBROUTINE sortie

c routine public du module mod_variables
c en cas de difficulté fermeture des fichiers .HR, .lis pgend

c Auteur: P.Morel, Département Cassiopée, O.C.A.
c CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : baratine

	CLOSE(UNIT=2)		!FICHIER .LIS	
	CLOSE(UNIT=53)		!fichier .HR
	IF(.NOT.baratine)THEN
	 CLOSE(UNIT=101)
	 CLOSE(UNIT=102)
	 CLOSE(UNIT=103)
	ENDIF
	CALL pgend
	STOP
	
	END SUBROUTINE sortie
