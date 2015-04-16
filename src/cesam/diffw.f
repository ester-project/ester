
c****************************************************************

	SUBROUTINE diffw(dlnro,grand_k,nu,nu_cin,n2mu,n2t,r,ro,y,deff,dh,dv)

c subroutine public du module mod_evol

c subroutine générique de calcul des coefficients de diffusion pour la rotation

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : langue, nom_diffw
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(2,0:1) :: y
	REAL (kind=dp), INTENT(in) :: dlnro, grand_k, nu, nu_cin, n2mu, n2t,
	1 r, ro
	REAL (kind=dp), INTENT(out) :: deff, dh, dv

c-------------------------------------------------------------------------

	SELECT CASE(nom_diffw)
	CASE('diffw_mpz')
	 CALL diffw_mpz(dlnro,grand_k,nu,nu_cin,n2mu,n2t,r,ro,y,deff,dh,dv)
	CASE('diffw_p03')
	 CALL diffw_p03(dlnro,grand_k,nu,nu_cin,n2mu,n2t,r,ro,y,deff,dh,dv)
	CASE('diffw_toul')
	 CALL diffw_toul(r,y,deff,dh,dv)
	CASE('diffw_0')
	 deff=0.d0 ; dh=0.d0 ; dv=0.d0
	CASE('diffw_cte')
	 CALL diffw_cte(deff,dh,dv)
	CASE DEFAULT
	 PRINT*
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)nom_diffw ; WRITE(2,1001)nom_diffw
1001	  FORMAT('STOP, unknown subroutine of diff. for the rotation',a,/,
	1 'known subroutines : diffw_0, diffw_cte, diffw_mpz, diffw_p03')
	 CASE DEFAULT
	  WRITE(*,1)nom_diffw ; WRITE(2,1)nom_diffw
1	  FORMAT('ARRET, routine de diffusion. pour la rotation inconnue: ',a,/,
	1 'SSP connues: diffw_0, diffw_cte, diffw_mpz, diffw_p03, diffw_toul')
	 END SELECT
	 STOP
	END SELECT
	
	RETURN

	END SUBROUTINE diffw
