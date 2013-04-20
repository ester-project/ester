	subroutine nuc_cesam_init
      
	USE mod_donnees, ONLY : lit_nl, nchim, nom_fich2, w_rot
	USE mod_kind
	USE mod_nuc, ONLY : nuc

	IMPLICIT NONE
	
	REAL (kind=dp), DIMENSION(0,0) :: jac
	REAL (kind=dp), DIMENSION(0) :: comp, dcomp, ex
	REAL (kind=dp) :: et, ero, hhe, be7e, b8e, n13e, o15e, f17e
	REAL (kind=dp), DIMENSION(5) :: epsilon	

	INTEGER :: fait

c lecture du fichier de donnees
	nom_fich2=ESTER_ROOT//'/tables/cesam_nuc/test'
	CALL lit_nl(w_rot)    

	fait=0
	CALL nuc(0d0,0d0,comp,dcomp,jac,.FALSE.,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	
	end
	
	
	
	
	subroutine nuc_cesam_init_abon(X,Z,comp)
      
	USE mod_donnees, ONLY : x0,y0,z0, nchim,nucleo
	USE mod_kind
	USE mod_nuc, ONLY : nuc

	IMPLICIT NONE
	
	REAL (kind=dp), DIMENSION(0,0) :: jac
	REAL (kind=dp), DIMENSION(0) :: dcomp, ex
	REAL (kind=dp), DIMENSION(nchim) :: comp
	REAL (kind=dp) :: et, ero, hhe, be7e, b8e, n13e, o15e, f17e, X, Z
	REAL (kind=dp), DIMENSION(5) :: epsilon	

	INTEGER :: i,fait
	
	x0=X;z0=Z;y0=1.-X-Z
	
	fait=1
	
	CALL nuc(0d0,0d0,comp,dcomp,jac,.FALSE.,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	
	comp=comp*nucleo
	
	end
	
	
	
	
	subroutine nuc_cesam_eps(t,ro,comp,epsilon,et,ero,ex)
	
	USE mod_donnees, ONLY : nucleo,nchim
	USE mod_kind
	USE mod_nuc, ONLY : nuc
	
	IMPLICIT NONE
	
	REAL (kind=dp), DIMENSION(0,0) :: jac
	REAL (kind=dp), DIMENSION(0) :: dcomp
	REAL (kind=dp), DIMENSION(nchim) :: comp,ex
	REAL (kind=dp) :: t,ro,et, ero, hhe, be7e, b8e, n13e, o15e, f17e
	REAL (kind=dp), DIMENSION(5) :: epsilon	
	INTEGER :: fait
	
	comp=comp/nucleo
	fait=3
	CALL nuc(t,ro,comp,dcomp,jac,.TRUE.,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	comp=comp*nucleo
	
	end
	
	
