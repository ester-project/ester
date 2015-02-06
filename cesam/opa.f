
c*********************************************************************

	 SUBROUTINE opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
	
c	 routine générique pour le calcul de l'opacité	
c	 il y a des appels différents suivant suivant nom_opa

c	 routine public du module mod_opa

c entrées :
c        xchim(1)=X : comp. chim. par gramme
c        t : température K
c        ro : densité cgs
  
c sorties :
c        kappa : opacité gr / cm2)
c        dkapdt : kappa / d t
c        dkapdr : kappa / d densité              
c        dkapdx : kappa / d xchim(1)
  
c        Z est obtenu par 1-X-Y

c	 Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	 CESAM2k
	
c--------------------------------------------------------------------

	 USE mod_donnees, ONLY : langue, nchim, nom_opa	
	 USE mod_kind

	 IMPLICIT NONE
   
         REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	 REAL (kind=dp), INTENT(in) :: t, ro
	 REAL (kind=dp), INTENT(out) ::	kap, dkapt, dkapro, dkapx
	 
         REAL (kind=dp), DIMENSION(nchim) :: xchi	 
	
c--------------------------------------------------------------------

	xchi=ABS(xchim(1:nchim))	
	SELECT CASE(nom_opa)
        CASE ('opa_gong')
         CALL opa_gong(t,ro,kap,dkapt,dkapro,dkapx)	 
        CASE ('opa_houdek9')
         CALL opa_houdek9(xchi,t,ro,kap,dkapt,dkapro,dkapx)	
c	CASE ('opa_houdek04')
c        CALL opa_houdek04(xchi,t,ro,kap,dkapt,dkapro,dkapx)	 
	CASE ('opa_int_zsx')
         CALL opa_int_zsx(xchi,t,ro,kap,dkapt,dkapro,dkapx)
        CASE ('opa_opalCO')
         CALL opa_opalCO(xchi,t,ro,kap,dkapt,dkapro,dkapx)
        CASE ('opa_opal2_cno')
         CALL opa_opal2(xchi,t,ro,kap,dkapt,dkapro,dkapx,.TRUE.)
        CASE ('opa_opal2_co')
         CALL opa_opal2(xchi,t,ro,kap,dkapt,dkapro,dkapx,.FALSE.)
	CASE ('opa_yveline')
         CALL opa_yveline(xchi,t,ro,kap,dkapt,dkapro,dkapx)
c        CASE ('opa_yveline04')
c         CALL opa_yveline04(xchi,t,ro,kap,dkapt,dkapro,dkapx)
c       CASE ('opa_yveline_lisse')
c        CALL opa_yveline_lisse(xchi,t,ro,kap,dkapt,dkapro,dkapx)
        CASE default
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1001)nom_opa ; WRITE(2,1001)nom_opa
1001	  FORMAT('STOP, unknown opacity : ',a,/,
	1 'known routines : opa_gong, opa_houdek9, opa_int_zsx,'
	2 'opa_opalCO, opa_opal2_co, opa_opal2_cno, opa_yveline,'
	3 'opa_yveline_lisse','opa_houdek04')	
	 CASE DEFAULT	  	 
	  WRITE(*,1)nom_opa ; WRITE(2,1)nom_opa
1	  FORMAT('ARRET, routine d''opacité inconnue : ',a,/,
	1 'routines connues: opa_gong, opa_houdek9, opa_int_zsx,'
	2 'opa_opalCO, opa_opal2_co, opa_opal2_cno, opa_yveline,'
	3 'opa_yveline_lisse','opa_houdek04')
	 END SELECT
	 STOP
	END SELECT
	
c opacités conductives

	 CALL kappa_cond(xchi,t,ro,kap,dkapt,dkapro,dkapx)
		
	RETURN
	
	END SUBROUTINE opa
