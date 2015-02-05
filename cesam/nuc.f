
c***********************************************************************

	SUBROUTINE nuc(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	
c	subroutine générique de calcul des taux des réactions
c	thermonucléaires

c	routine public du module mod_nuc

c-------------------------------------------------------------------

	USE mod_donnees, ONLY : nom_nuc	
	USE mod_kind
	
	IMPLICIT NONE
	  
	REAL (kind=dp), INTENT(in):: t, ro		
	INTEGER, INTENT(in) :: fait
	LOGICAL, INTENT(in) :: deriv	  
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: comp
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: jac
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dcomp, ex, epsilon
	REAL (kind=dp), INTENT(out) :: et, ero, hhe, be7e, b8e, n13e,
	1 o15e, f17e
	
c--------------------------------------------------------------------

2000	FORMAT(8es10.3)
  
c------- Pertes neutriniques--------------------------------------- !YLD
c
c        CALL neutrinos()
c--------------------------------------------------------------------!YLD
	  
	SELECT CASE(nom_nuc)
	CASE ('iben')
	 CALL iben(t,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)	  
	CASE ('pp1')
	 CALL pp1(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE ('pp3')
	 CALL pp3(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)		  
	CASE ('ppcno9')
	 CALL ppcno9(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE ('ppcno9Fe')
	 CALL ppcno9Fe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)	
	CASE('ppcno10')
	 CALL ppcno10(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno10Fe')
	 CALL ppcno10Fe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno10K')
	 CALL ppcno10K(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno10BeBFe')
	 CALL ppcno10BeBFe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)		
	CASE('ppcno11')
	 CALL ppcno11(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno12')
	 CALL ppcno12(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno12Be')
	 CALL ppcno12Be(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno12BeBFe')
	 CALL ppcno12BeBFe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)		
	CASE('ppcno12Li')
	 CALL ppcno12Li(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno3a9')
	 CALL ppcno3a9(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno3ac10')
	 CALL ppcno3ac10(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)	
	CASE DEFAULT
	 PRINT*,'routine de réactions nucléaires inconnue: ',nom_nuc
	 PRINT*,'connues: iben, pp1, pp3, ppcno9, ppcno9Fe, ppcno10'
	 PRINT*,'ppcno10Fe, ppcno10K, ppcno10BeBFe, ppcno11, ppcno12'
	 PRINT*,'ppcno12Be, ppcno3a9, ppcno3ac10, ppcno12Li, ppcno12BeBFe'
	 PRINT*,'arrêt' ; STOP
	END SELECT
	
	RETURN
	
	END SUBROUTINE nuc
