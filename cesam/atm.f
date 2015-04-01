
c***********************************************************************

	SUBROUTINE atm(list,l_rac,r_rac,xchim,pt_rac,dptsdl,dptsdr,
	1 t_rac,dtsdl,dtsdr,m_rac,dmsdl,dmsdr,p_rac,dpsdl,dpsdr,t_eff)

c	routine générique pour la restitution de l'atmosphère
c	il y a des appels différents suivant nom_atm

c	routine public du module mod_atm

c entrées :
c	list=.true. : calcul réduit pour une liste
c	r_rac : rayon au raccord
c	l_rac : luminosité au raccord
c	xchim : composition chimique par gramme

c sorties :
c	  pt_rac : pression totale au raccord,
c	  dptsdl : dérivée / L de la pression totale au raccord,   
c	  dptsdr : dérivée / R de la pression totale au raccord,   
c	  t_rac : température au raccord,
c	  dtsdl : dérivée / L de la température au raccord,    
c	  dtsdr : dérivée / R de la température au raccord,    
c	  m_rac : masse au raccord,
c	  dmsdl : dérivée / L de la masse au raccord,   
c	  dmsdr : dérivée / R de la masse au raccord,       
c	  p_rac : pression gazeuse au raccord,
c	  dpsdl : dérivée / L de la pression gazeuse au raccord,   
c	  dpsdr : dérivée / R de la pression gazeuse au raccord,    
c	  t_eff  : température effective. 

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c--------------------------------------------------------------------

	USE mod_donnees, ONLY : baratine, nom_atm, nom_fich2
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: l_rac, r_rac
	LOGICAL, INTENT(in) :: list    
	REAL (kind=dp), INTENT(out) :: pt_rac, dptsdl, dptsdr, t_rac,
	1 dtsdl, dtsdr, m_rac, dmsdl, dmsdr, p_rac, dpsdl, dpsdr, t_eff

	LOGICAL, SAVE :: init=.TRUE.

c-------------------------------------------------------------

	IF(init)THEN
	 init=.FALSE.
	 
c baratine=.FALSE. permet de dérouter sur le fichier mon_modèle_102 les
c informations concernant le déroulement des calculs pour l'atmosphère	
	 IF(baratine)THEN
	  usl_atm=6
	 ELSE	  	 
	  usl_atm=102
	  OPEN(unit=102,form='formatted',status='unknown',!access='append',
	1 file=TRIM(nom_fich2)//'_atm')
	 ENDIF
	ENDIF

	SELECT CASE(nom_atm)
	CASE('lim_gong1')
	 CALL lim_gong1(l_rac,r_rac,xchim,pt_rac,dptsdl,dptsdr,
	1 t_rac,dtsdl,dtsdr,m_rac,dmsdl,dmsdr,p_rac,dpsdl,dpsdr,t_eff)  
	CASE('lim_tau1')
	 CALL lim_tau1(l_rac,r_rac,xchim,pt_rac,dptsdl,dptsdr,
	1  t_rac,dtsdl,dtsdr,m_rac,dmsdl,dmsdr,p_rac,dpsdl,dpsdr,t_eff)
	CASE('lim_atm')
	 CALL lim_atm(list,l_rac,r_rac,xchim,pt_rac,dptsdl,dptsdr,
	1 t_rac,dtsdl,dtsdr,m_rac,dmsdl,dmsdr,p_rac,dpsdl,dpsdr,t_eff)
	CASE DEFAULT
	 PRINT*,'routine de restitution d''atmosphère inconnue: ',nom_atm
	 PRINT*,'routines connues: lim_gong1, lim_tau1, lim_atm'
	 PRINT*,'arrêt' ; STOP
	END SELECT
	
	RETURN

	END SUBROUTINE atm
