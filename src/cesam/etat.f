	
c***********************************************************************

	SUBROUTINE etat(p,t,xchim,deriv,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c subroutine générique pour le calcul de l'équation d'état
c routine public du module mod_etat	

c entrées:
c	  p : pression
c	  t : temperature
c	  xchim : composition chimique par gramme
c	  deriv=.FALSE. : évite le calcul  de certaines dérivées

c sorties:
c	  ro : densité et dérivées
c	  u : énergie interne et dérivées

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : nom_etat
	USE mod_kind
	
	IMPLICIT NONE
	  			
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim	
	REAL (kind=dp), INTENT(in) :: p, t
	LOGICAL, INTENT(in) :: deriv	  		
	REAL (kind=dp), INTENT(out) :: ro, drop, drot, drox, u, dup,
	1 dut, dux,delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
	2 gradad, dgradadp, dgradadt, dgradadx, alfa, beta, gamma1
	
c---------------------------------------------------------------------
	
	SELECT CASE(nom_etat)	  
	CASE ('etat_ceff')
	 CALL etat_ceff(p,t,xchim,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_eff')
	 CALL etat_eff(p,t,xchim,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_gong1')
	 CALL etat_gong1(p,t,xchim,deriv,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_gong2')
	 CALL etat_gong2(p,t,xchim,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_mhd')
	 CALL etat_mhd(p,t,xchim,deriv,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_opal')
	 CALL etat_opal(p,t,xchim,deriv,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_opalX')
	 CALL etat_opalX(p,t,xchim,deriv,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	CASE ('etat_opalZ')
	 CALL etat_opalZ(p,t,xchim,deriv,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)	
	CASE DEFAULT
	 PRINT*,'routine d''équation d''état inconnue: ',nom_etat
	 PRINT*,'routines connues: etat_ceff, etat_eff, etat_gong1'	   
	 PRINT*,'etat_gong2, etat_mhd, etat_opal, etat_opalX, etat_opalZ'
	 PRINT*,'arret' ; STOP
	END SELECT
	  
c les dérivées / X ne sont fiables que pour X > 1.d-2 
c	IF(xchim(1) > 1.d-2)THEN
c	 drox=0.d0 ; dux=0.d0 ; deltax=0.d0 ; dcpx=0.d0 ; dgradadx=0.d0
c	ENDIF
	  
	RETURN

	END SUBROUTINE etat
