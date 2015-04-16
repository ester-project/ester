
c*****************************************************************
	
	SUBROUTINE df_rotx(p,ro,t,alfa,delta,dfdp_tx,dfdt_px,dfdx_pt,
	1 drodx_pt,df_dro,df_dt,df_dx)
	
c	transformation des dérivées partielles P, T, X --> ro, T, X
c	pour une fonction f

c	routine public du module mod_etat

c entrées:
c	p: pression
c	ro: densité
c	t: température
c	alfa: dln ro/dln P
c	delta: -dln ro/dln T
c	dfdp_tx: df/dP(T,X)
c	dfdt_px: df/dT(P,X)
c	dfdx_pt: df/dX(P,T)
c	drodx_pt: dro/dX(P,T)

c sorties:
c	df_dro=df/dro(T,X) 
c	df_dt=df/dT(ro,X) 
c	df_dx=df/dX(ro,T) 

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM5

c---------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: alfa, delta,
	1 dfdp_tx, dfdt_px, dfdx_pt, drodx_pt, p, ro, t
	REAL (kind=dp), INTENT(out) :: df_dro, df_dt, df_dx	
	
	df_dro=p/ro/alfa*dfdp_tx			!df/dro(T,X)
	df_dt=dfdt_px+p*delta/t/alfa*dfdp_tx		!df/dT(ro,X)	
	df_dx=dfdx_pt-p/ro/alfa*drodx_pt*dfdp_tx	!df/dX(ro,T)
	
	RETURN
		
	END SUBROUTINE df_rotx
