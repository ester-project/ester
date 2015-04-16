
c*******************************************************

	SUBROUTINE base_rota

c routine subordonnée de evol
c Formation de la base pour la rotation

c r_qs=1 ordre des équations différentielles, est un PARAMETER de mod_donnees
c ord_rot est défini dans initialise_rota et initialise rota4

c Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c------------------------------------------------------------------
	
	INTEGER, ALLOCATABLE, DIMENSION(:) :: mult

c------------------------------------------------------------------

2000	FORMAT(8es10.3)
		
c multiplicités
	ALLOCATE(mult(n_rot))
	
c discontinuité de la dérivée 1-ière et  discontinuité à chaque limite ZR/ZC
	ord_rot=m_rot+r_qs ; mult=m_rot 			
	mult(1)=ord_rot ; mult(n_rot)=ord_rot
	 
c vecteur nodal	
	knotr=SUM(mult) ; DEALLOCATE(mrott) ; ALLOCATE(mrott(knotr))
	CALL noeud(mrot,mrott,mult,n_rot,knotr) ; dim_rot=knotr-ord_rot
	DEALLOCATE(mult)
	
	RETURN
	
	END SUBROUTINE base_rota
