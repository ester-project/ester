
c*********************************************************************

	PROGRAM cesam2k_dbg
	
c	execution de cesam2k_dbg
c	Programme du sous directory CESAM2k/SOURCE

c	la routine cesam constituant le programme principal a été placée
c	dans le module mod_cesam. Cette disposition permet de mettre
c	cesam dans la bibliothèque et évite de recompiler le programme
c	principal lors de tests ou de mises au point

c	cesam2k-dbg est destiné aux modifications importantantes,
c	typiquement aux restructuration
	
c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c    	CESAM2k

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	USE mod_cesam, ONLY : cesam
	
	IMPLICIT NONE		
	
	CALL cesam
	
	STOP
	
	END PROGRAM cesam2k_dbg
