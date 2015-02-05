
c********************************************************

	SUBROUTINE sortie

c	routine public du module mod_variables
c	en cas de difficulté fermeture des fichiers .HR, .lis pgend

c	Auteur: P.Morel, Département Cassiopée, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------

	CLOSE(unit=2)		!fichier .lis	
	CLOSE(unit=53)		!fichier .HR	
	CALL pgend
	STOP
	
	END SUBROUTINE sortie
