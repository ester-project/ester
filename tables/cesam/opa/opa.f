
c*********************************************************************

	 SUBROUTINE opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)

c routine public du module mod_opa
	
c routine générique pour le calcul de l'opacité	
c il y a des appels différents suivant suivant nom_opa

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

	USE mod_donnees, ONLY : langue, f_opa, nchim, nom_opa	
	USE mod_kind
	USE mod_numerique, ONLY : polyder
	
	IMPLICIT NONE
   
        REAL (kind=dp), INTENT(in), DIMENSION(nchim) :: xchim
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out) ::	kap, dkapt, dkapro, dkapx

	REAL (kind=dp), PARAMETER, DIMENSION(0:3) ::
	1 a=(/ 0.1298d0,-0.1856d0 ,0.1064d0,-0.0345d0 /)
        REAL (kind=dp), DIMENSION(nchim) :: xchi
	REAL (kind=dp), DIMENSION(0:3) :: p
	REAL (kind=dp), PARAMETER :: t_compton=6.d7
	REAL (kind=dp) :: f, zeta		
	 	 	 	
c--------------------------------------------------------------------

2000	FORMAT(8es10.3)	 

c normalisation de la composition chimique, retrait des valeurs négatives
	print*,'Je suis la 1 xchi=',xchim
	 xchi=ABS(xchim) ; WHERE(xchi < 1.d-30)xchi=0.d0 ; xchi=xchi/SUM(xchi)	 

	print*,'Je suis la xchi=',xchi

c Si t > t_compton les éléments chimiques sont totalement ionisés
c utilisation de l'opacité compton

	IF(t >= t_compton)THEN
	 CALL opa_compton(xchi,t,ro,kap,dkapt,dkapro,dkapx)
	 RETURN
	ENDIF
	
c opacités tabulées	
	SELECT CASE(nom_opa)
        CASE ('opa_gong')
         CALL opa_gong(t,ro,kap,dkapt,dkapro,dkapx)	 
c        CASE ('opa_houdek9')
c         CALL opa_houdek9(xchi,t,ro,kap,dkapt,dkapro,dkapx) 
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

c correction de Jorgen astro-ph 0811.100v1 6 Nov 2008
	CASE ('opa_yveline_jorgen')	
         CALL opa_yveline(xchi,t,ro,kap,dkapt,dkapro,dkapx)
	 IF(TRIM(f_opa(1)) == 'opa_yveline_ags.bin')THEN
	  IF(t < 3.2d7)THEN	
	   zeta=MAX(LOG10(t)-6.3d0,0.d0)
	   CALL polyder(a,3,1,zeta,p)
	   f=10.d0**p(0)
	   kap=kap*f ; dkapt=f*(dkapt+p(1)/t) 
	  ENDIF
	 ENDIF	
	 
        CASE ('opa_yveline_lisse')
         CALL opa_yveline_lisse(xchi,t,ro,kap,dkapt,dkapro,dkapx)
        CASE DEFAULT
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1001)nom_opa ; WRITE(2,1001)nom_opa
1001	  FORMAT('STOP, unknown opacity : ',a,/,
	1 'known routines : opa_gong, opa_houdek9, opa_int_zsx,'
	2 'opa_opalCO, opa_opal2_co, opa_opal2_cno, opa_yveline,'
	3 'opa_yveline_lisse, opa_yv_liss_corr')	
	 CASE DEFAULT	  	 
	  WRITE(*,1)nom_opa ; WRITE(2,1)nom_opa
1	  FORMAT('ARRET, routine d''opacité inconnue : ',a,/,
	1 'routines connues: opa_gong, opa_houdek9, opa_int_zsx,'
	2 'opa_opalCO, opa_opal2_co, opa_opal2_cno, opa_yveline,'
	3 ' opa_yveline_lisse, opa_yveline_jorgen')
	 END SELECT
	 STOP
	END SELECT
	
c opacités conductives
	IF(nom_opa /= 'opa_houdek9')CALL kappa_cond(xchi,t,ro,kap,dkapt,
	1 dkapro,dkapx)
		
	RETURN
	
	END SUBROUTINE opa
