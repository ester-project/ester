
c******************************************************************

	MODULE mod_opa

c	module regroupant les routines de CESAM2k
c	concernant les routines d'opacité et leurs routines propres
c	d'exploitation

c	le calcul de l'opacité est effectué par la routine générique opa
c	les opacités sont différenciées par leur nom: nom_opa
c	lu dans lit_nl mis dans mod_donnees

c fonctions private:
c	kappa_cond : opacités conductives
c	opa_gong : opacités simplifiees (Kramers ameliore)
c	opa_houdek9 : opacités de Houdek version 9, (OPAL+Alexander),
c	interpolation par rational B-spline
c	opa_int_zsx : opacités OPAL interpolation linéaires
c	opa_yveline : opacités OPAL+Alexander interp. et raccord d'Yveline
c	extension arbitraire X=1
c	opa_yveline_lisse : opacités OPAL+Alexander raccord d'Yveline,
c	interpolation linéaire et extension arbitraire X=1

c fonction public:
c	opa : routine générique d'opacité

c	Auteurs: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-------------------------------------------------------

	PRIVATE
	PUBLIC :: opa
	
	CONTAINS
	
c------------------------------------------------------------------	
	
	INCLUDE 'kappa_cond.f'
	INCLUDE 'opa.f'
	INCLUDE 'opa_gong.f'
	INCLUDE 'opa_houdek9.f'   !YLD
c	INCLUDE 'opa_houdek04.f'
	INCLUDE 'opa_int_zsx.f'
	INCLUDE 'opa_opalCO.f'
	INCLUDE 'opa_opal2.f'
	INCLUDE 'opa_yveline.f'
c	INCLUDE 'opa_yveline_lisse.f'			
	
	END MODULE mod_opa

	INCLUDE 'z14xcotrin21.f'
