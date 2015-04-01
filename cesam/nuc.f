
c***********************************************************************

	SUBROUTINE nuc(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	
c routine public du module mod_nuc
	
c subroutine générique de calcul des taux des réactions thermonucléaires

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c-------------------------------------------------------------------

	USE mod_donnees, ONLY : nom_nuc	
	USE mod_kind
	
	IMPLICIT NONE
	  
	REAL (kind=dp), INTENT(in):: ro, t		
	INTEGER, INTENT(in) :: fait
	LOGICAL, INTENT(in) :: deriv	  
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: comp
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: jac
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dcomp, epsilon, ex
	REAL (kind=dp), INTENT(out) ::  be7e, b8e, ero, et, f17e, hhe, n13e,
	1 o15e

c--------------------------------------------------------------------

2000	FORMAT(8es10.3)
  	  
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
	CASE('ppcno3a12Ne')
	 CALL ppcno3a12Ne(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)	
	CASE('ppcno3a9')
	 CALL ppcno3a9(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE('ppcno3aco')
	 CALL ppcno3aco(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	CASE DEFAULT
	 PRINT*,'routine de réactions nucléaires inconnue: ',nom_nuc
	 PRINT*,'connues: iben, pp1, pp3, ppcno9, ppcno9Fe, ppcno10'
	 PRINT*,'ppcno10Fe, ppcno10K, ppcno10BeBFe, ppcno11, ppcno12'
	 PRINT*,'ppcno12Be, ppcno3a9, ppcno12Li, ppcno12BeBFe'
	 PRINT*,'ppcno3a12Ne, ppcno3aco' 
	 PRINT*,'arrêt dans nuc' ; STOP
	END SELECT
	
	RETURN
	
	END SUBROUTINE nuc
