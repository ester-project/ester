     
c****************************************************************
	 
	SUBROUTINE diffw(dlnro,grand_k,nu,n2mu,n2t,r,ro,dfrot,frot,
	1 deff,dh,dv)

c	subroutine private du module mod_evol
	 	
c	subroutine générique de calcul des coefficients de
c	diffusion turbulente pour la rotation

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : langue, nom_diffw, nrot
	USE mod_kind	 

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(nrot) :: dfrot, frot 	
	REAL (kind=dp), INTENT(in) :: dlnro, grand_k, nu, n2mu, n2t, r, ro
	REAL (kind=dp), INTENT(out) :: deff, dh, dv
	 
c-------------------------------------------------------------------------
	 
	SELECT CASE(nom_diffw)
	CASE('diffw_mpz')
	 CALL diffw_mpz(dlnro,grand_k,nu,n2mu,n2t,r,ro,dfrot,frot,
	1 deff,dh,dv)
	CASE('diffw_0')
	 deff=0.d0 ; dh=0.d0 ; dv=0.d0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1002)nom_diffw ; WRITE(2,1002)nom_diffw
1002	  FORMAT('WARNING, the diff. routine used for the rotation : ',a,/,
	1  'has no interest')	 
	 CASE DEFAULT
	  WRITE(*,2)nom_diffw ; WRITE(2,2)nom_diffw
2	  FORMAT('ATTENTION, la routine de diff. pour la rotation : ',a,/,
	1 'n''a guère d''intérêt')
	 END SELECT
	CASE DEFAULT
	 PRINT*
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)nom_diffw ; WRITE(2,1001)nom_diffw
1001	  FORMAT('STOP, unknown subroutine of diff. for the rotation',a,/,
	1 'known subroutines : diffw_mpz')
	 CASE DEFAULT
	  WRITE(*,1)nom_diffw ; WRITE(2,1)nom_diffw
1	  FORMAT('ARRET, routine de diff. pour la rotation inconnue: ',a,/,
	1 'routines connues : diffw_0, diffw_mpz')
	 END SELECT
	 STOP
	END SELECT
	 	
	RETURN
	 
	END SUBROUTINE diffw
