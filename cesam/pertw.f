	
c************************************************************************

         SUBROUTINE pertw(nu,omega,r,mw_dot)

c	routine générique pour le calcul de la perte/gain de moment
c	cinétique
c	il y a des appels différents suivant nom_pertw

c	routine private du module mod_evol

c	Auteur: P.Morel, Département  Cassiopée, O.C.A.
c	CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, nom_pertw
	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: nu, omega, r	
	REAL (kind=dp), INTENT(out) :: mw_dot
	
	LOGICAL, SAVE :: pass=.TRUE.
		
c---------------------------------------------------------------- 

	SELECT CASE(nom_pertw)
	CASE ('pertw_sch')
	 CALL pertw_sch(nu,omega,r,mw_dot)
	CASE ('pertw_loc')
	 CALL pertw_loc(nu,omega,r,mw_dot)
	CASE ('pertw_ptm')
	 CALL pertw_ptm(nu,r,mw_dot)
	CASE ('pertw_0')
	 mw_dot=0.d0
	 IF(pass)THEN
	  pass=.FALSE.
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1000) ; WRITE(2,1000)
1000	   FORMAT('Without angular momentum losses')
	  CASE DEFAULT	 
	   WRITE(*,1) ; WRITE(2,1)
1	   FORMAT('Sans perte de moment cinétique')
	  END SELECT
	 ENDIF	 
	CASE DEFAULT	
	 SELECT CASE(langue)
	 CASE('english')
	  PRINT*,'STOP, unknown routine of angular momentum lost : ',
	1  nom_pertw
	  PRINT*,'known routines: pertw_sch, pertw_loc, pertw_ptm, pertw_0'	 
	 CASE DEFAULT
	  PRINT*,'ARRET, routine de perte de moment cinétique inconnue: ',
	1   nom_pertw
	  PRINT*,'routines connues: pertw_sch, pertw_loc, pertw_ptm, pertw_0'
	 END SELECT
	 STOP
	END SELECT

	RETURN
	
	END SUBROUTINE pertw
