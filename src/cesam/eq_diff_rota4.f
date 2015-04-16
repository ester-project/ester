
c***********************************************************************

	SUBROUTINE eq_diff_rota4(dt,fait,nu,y,as,bs)

c routine private du module mod_evol

c formation des équations à résoudre dans resout_rota4
c pour la rotation par la méthode decollocation
c formalisme de Mathis & Zahn 2004, Krot=4

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c entrées
c	fait=0 : équations aux points de collocation
c	fait=1 : conditions au centre
c	fait=2 : conditions à la limite externe toujours convectif
c	nu: point de calcul en m**2/3
c	y: solution
c	dt : pas temporel

c sorties
c	as : coefficients de la spline/dérivée
c	bs : seconds membres

c Variables
c	y(1,*): Omega
c	y(2,*): U
c	y(3,*): Psi
c	y(4,*): Lambda
c	y(5,*): Tau
c	y(6,*): Upsilon
c	y(7,*): Phi
c	y(8,*): Pi

c---------------------------------------------------------------------

	USE mod_donnees, ONLY :  nchim, nom_rot, nrl, nrot
	USE mod_kind
	USE mod_variables, ONLY : mrot

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nrot,0:1) :: y
	REAL (kind=dp), INTENT(in) :: dt, nu
	INTEGER, INTENT(in) :: fait

	REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: as
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: bs

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: ys
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: bss, frls
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: coeff
	REAL (kind=dp), DIMENSION(nrl) :: frl
	REAL (kind=dp), PARAMETER :: dx=1.d-5, unpdx=1.d0+dx

	REAL (kind=dp) :: dstor, stor, stor0

	INTEGER :: i, id, j

	LOGICAL :: deriv=.FALSE., zc

c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c le nombre d'inconnues est nrot=8
c y(1:nrot,der): dérivée d'ordre der=0 pour la fonction, der=1 pour la dérivée

c mises à 0
	as=0.d0 ; bs=0.d0
	
c ZC?
	zc=lmix(nu)	

c les équations aux points de collocation, fait=0
	SELECT CASE(fait)
	CASE(0)

c les coefficients
	 CALL coeff_rota4(dt,nu,y,frl)
c dans ZC
	 IF(zc)THEN

c homogènisation de la vitesse angulaire par diffusion
	  bs(1)=frl(26)*y(1,0)-frl(30)-y(5,1)
	  as(1,1,0)=frl(26)
	  as(1,5,1)=-1.d0

c U fictif
	  bs(2)=frl(31)*y(1,0)*y(2,0)-y(5,0)
	  as(2,1,0)=frl(31)*y(2,0)	  	  
	  as(2,2,0)=frl(31)*y(1,0)
	  as(2,5,0)=-1.d0
	  
c Psi=0
	  bs(3)=y(3,0) ; as(3,3,0)=1.d0

c Lambda = 0
	  bs(4)=y(4,0) ; as(4,4,0)=1.d0
	  
c définition de Tau (Flux de moment cinétique)
	  bs(5)=frl(32)*y(1,1)-y(5,0)  
	  as(5,1,1)=frl(32)	  
	  as(5,5,0)=-1.d0

c Upsilon = 0
	  bs(6)=y(6,0) ; as(6,6,0)=1.d0
	  
c équation de Poisson, Phi
	  bs(7)=frl(29)*y(1,0)**2+frl(28)*y(1,0)*y(1,1)+frl(27)*y(7,0)-y(8,1)
	  as(7,1,0)=2.d0*frl(29)*y(1,0)+frl(28)*y(1,1)
	  as(7,1,1)=frl(28)*y(1,0)
	  as(7,7,0)=frl(27)
	  as(7,8,1)=-1.d0

c définition de Pi
	  bs(8)=y(7,0)+frl(25)*y(7,1)-y(8,0)
	  as(8,7,0)=1.0d0
	  as(8,7,1)=frl(25)
	  as(8,8,0)=-1.0d0

c dans ZR
	 ELSE

c diffusion du moment cinétique
	  bs(1)=frl(26)*y(1,0)-frl(30)-y(5,1)
	  as(1,1,0)=frl(26)
	  as(1,5,1)=-1.d0

c définition de U
	  bs(2)=frl(1)*y(1,0)**4+frl(2)*y(1,0)**3*y(1,1)
	1 +frl(3)*y(1,0)**2+frl(4)*y(1,0)*y(1,1)
	2 +frl(5)*y(1,0)**2*y(7,1)+frl(6)*y(1,0)**2*y(7,0)
	3 +frl(7)*y(1,0)*y(1,1)*y(7,1)+frl(8)*y(1,0)*y(7,0)*y(1,1)  
	4 +frl(9)*y(7,0)+frl(10)*y(7,1)+frl(12)*y(3,1)
	5 +frl(13)*y(3,0)+frl(14)*y(4,0)+frl(15)
	6 +frl(16)*y(2,0)-y(6,1)
 
	  as(2,1,0)=4.d0*frl(1)*y(1,0)**3
	1 +3.d0*frl(2)*y(1,0)**2*y(1,1)+2.d0*frl(3)*y(1,0)
	2 +frl(4)*y(1,1)+2.d0*frl(5)*y(1,0)*y(7,1)+2.d0*frl(6)*y(1,0)*y(7,0)
	3 +frl(7)*y(1,1)*y(7,1)+frl(8)*y(7,0)*y(1,1) 		
	  as(2,1,1)=frl(2)*y(1,0)**3+frl(4)*y(1,0)
	1 +frl(7)*y(1,0)*y(7,1)+frl(8)*y(1,0)*y(7,0)	
	  as(2,2,0)=frl(16) 	  	  
	  as(2,3,0)=frl(13)	  	  
	  as(2,3,1)=frl(12)
	  as(2,4,0)=frl(14)
	  as(2,6,1)=-1.d0
	  as(2,7,0)=frl(6)*y(1,0)**2+frl(8)*y(1,0)*y(1,1)+frl(9) 
	  as(2,7,1)=frl(5)*y(1,0)**2+frl(7)*y(1,0)*y(1,1)+frl(10) 

c définition de Psi
	  bs(3)=frl(19)*y(1,0)*y(1,1)+frl(20)*y(4,0)+frl(21)*y(3,0)
	  as(3,1,0)=frl(19)*y(1,1)
	  as(3,1,1)=frl(19)*y(1,0)
	  as(3,3,0)=frl(21)
	  as(3,4,0)=frl(20)
	  
c évolution de Lambda
	  bs(4)=frl(22)*y(2,0)+frl(23)*y(4,0)+frl(24)  
	  as(4,2,0)=frl(22)
	  as(4,4,0)=frl(23)	  

c définition de Tau (Flux de moment cinétique)
	  bs(5)=frl(31)*y(1,0)*y(2,0)+frl(32)*y(1,1)-y(5,0)
	  as(5,1,0)=frl(31)*y(2,0)	  
	  as(5,1,1)=frl(32)	  
	  as(5,2,0)=frl(31)*y(1,0)
	  as(5,5,0)=-1.d0

c définition de Upsilon
	  bs(6)=frl(17)*y(3,1)+frl(11)*y(4,0)+frl(18)*y(3,0)-y(6,0)
	  as(6,3,0)=frl(18)
	  as(6,3,1)=frl(17)
	  as(6,4,0)=frl(11)	  
	  as(6,6,0)=-1.d0
	  
c équation de Poisson, Phi
	  bs(7)=frl(29)*y(1,0)**2+frl(28)*y(1,0)*y(1,1)+frl(27)*y(7,0)-y(8,1)
	  as(7,1,0)=2.d0*frl(29)*y(1,0)+frl(28)*y(1,1)
	  as(7,1,1)=frl(28)*y(1,0)
	  as(7,7,0)=frl(27)
	  as(7,8,1)=-1.d0
	  
c définition de Pi
	  bs(8)=y(7,0)+frl(25)*y(7,1)-y(8,0)
	  as(8,7,0)=1.0d0 
	  as(8,7,1)=frl(25)
	  as(8,8,0)=-1.0d0

	 ENDIF

c~~~~~~~~~limites ZR/ZC~~~~~~~~~~~~~~~~~~~~~~~~~

c fait=1, centre, 5 conditions
	CASE(1)
	
c condition dOmega=0
	 bs(1)=y(1,1) ; as(1,1,1)=1.d0
	  
c condition U=0
	 bs(2)=y(2,0) ; as (2,2,0)=1.d0
	 
c condition Lambda=0
	 bs(3)=y(4,0) ; as(3,4,0)=1.d0
	 
c condition Upsilon=0
	 bs(4)=y(6,0) ; as(4,6,0)=1.d0
	 
c condition Phi=0 
	 bs(5)=y(7,0) ; as(5,7,0)=1.d0
	 
c condition Psi = 0
c	 bs(6)=y(3,0) ; as(6,3,0)=1.d0
	 
	 
c limite externe fait=2, 3 conditions
	CASE(2)
	
c les coefficients
c	 CALL coeff_rota4(dt,nu,y,frl)

c condition flux = Tau = 0	  
	 bs(1)=y(5,0) ; as(1,5,0)=1.d0
	  
c condition sur Pi=-2 Phi	 
	 bs(2)=2.d0*y(7,0)+y(8,0)
	 as(2,7,0)=2.d0	 
 	 as(2,8,0)=1.d0	 
	  	  
c condition Psi = 0
	 bs(3)=y(3,0) ; as(3,3,0)=1.d0

	END SELECT

c-----------------vérification des dérivées-------------------

c	deriv=.TRUE.
c	deriv=(nu >= 6.E-01 .AND. nu <= 6.03E-01)	!ZR
c	1 .OR. (nu >= 1.E-01 .AND. nu <= 1.03E-01)	!ZC
c	2 .OR. (fait > 0)				!limite
c	deriv=lmix(nu)
c	deriv=.NOT.lmix(nu)
c	deriv= fait == 0
c	deriv= fait == 1
c	deriv= fait == 2
c	deriv= fait > 0
c	deriv= nu > 0.8d0 .AND. nu < 0.9d0
c	deriv= nu > 0.141d0 .AND. nu < 0.142d0
c	deriv= deriv .OR. nu > 1.155d0 .AND. nu < 1.157d0
c	deriv= nu > mrot(convf(1)-1) .AND. nu < mrot(convf(1)+1)
c	deriv= deriv .OR. nu > mrot(convd(2)-1) .AND. nu < mrot(convd(2)+1)
c	deriv= nu > mrot(convf(1)-3) .AND. nu < mrot(convf(1)+5)
c	deriv= nu > mrot(convd(1)-1) .AND. nu < mrot(convd(1)+1)

	IF(deriv)THEN
	 IF(zc)THEN
	  PRINT*,'convectif'
	 ELSE
	  PRINT*,'radiatif'
	 ENDIF
	 PRINT*,'nu=',nu,', fait=',fait
	 PRINT*,'y(:,0)' ; WRITE(*,2000)y(:,0)
	 PRINT*,'y(:,1)' ; WRITE(*,2000)y(:,1)
	 PRINT*,'bs' ; WRITE(*,2000)bs
	 PRINT*,'frl' ; WRITE(*,2000)frl
	 ALLOCATE(bss(nrot),coeff(ncoeff+nchim),frls(nrl),ys(nrot,0:1))
	 CALL coeff_rota4(dt,nu,y,frl,coeff)
	 PRINT*,'coeff' ; WRITE(*,2000)coeff(1:ncoeff)
	 DO id=0,1
	  DO i=1,nrot
	   ys=y
	   stor0=ys(i,id) ; stor=stor0*unpdx ; IF(stor == 0.d0)stor=dx
	   dstor=stor-stor0 ; ys(i,id)=stor ; bss=0.d0
	   CALL coeff_rota4(dt,nu,ys,frls)
	   SELECT CASE(fait)
	   CASE(0)
	    IF(zc)THEN	    
	     bss(1)=frls(26)*ys(1,0)-frls(30)-ys(5,1)
	     bss(2)=frls(31)*ys(1,0)*ys(2,0)-ys(5,0)
	     bss(3)=ys(3,0)
	     bss(4)=ys(4,0)
	     bss(5)=frls(32)*ys(1,1)-ys(5,0)     
	     bss(6)=ys(6,0)
	     bss(7)=frls(29)*ys(1,0)**2+frls(28)*ys(1,0)*ys(1,1)
	1    +frls(27)*ys(7,0)-ys(8,1)
	     bss(8)=ys(7,0)+frls(25)*ys(7,1)-ys(8,0)
	     	     
	    ELSE
	     bss(1)=frls(26)*ys(1,0)-frls(30)-ys(5,1)	      
	     bss(2)=frls(1)*ys(1,0)**4+frls(2)*ys(1,0)**3*ys(1,1)
	1    +frls(3)*ys(1,0)**2+frls(4)*ys(1,0)*ys(1,1)
	2    +frls(5)*ys(1,0)**2*ys(7,1)+frls(6)*ys(1,0)**2*ys(7,0)
	3    +frls(7)*ys(1,0)*ys(1,1)*ys(7,1)+frls(8)*ys(1,0)*ys(7,0)*ys(1,1)
	4    +frls(9)*ys(7,0)+frls(10)*ys(7,1)+frls(12)*ys(3,1)
	5    +frls(13)*ys(3,0)+frls(14)*ys(4,0)+frls(15)
	6    +frls(16)*ys(2,0)-ys(6,1)	
	     bss(3)=frls(19)*ys(1,0)*ys(1,1)+frls(20)*ys(4,0)+frls(21)*ys(3,0)
	     bss(4)=frls(22)*ys(2,0)+frls(23)*ys(4,0)+frls(24)
	     bss(5)=frls(31)*ys(1,0)*ys(2,0)+frls(32)*ys(1,1)-ys(5,0)
	     bss(6)=frls(17)*ys(3,1)+frls(11)*ys(4,0)+frls(18)*ys(3,0)-ys(6,0)
	     bss(7)=frls(29)*ys(1,0)**2+frls(28)*ys(1,0)*ys(1,1)
	1    +frls(27)*ys(7,0)-ys(8,1)
	     bss(8)=ys(7,0)+frls(25)*ys(7,1)-ys(8,0)	
	    ENDIF
	   CASE(1)	!centre
	    bss(1)=ys(1,1)
	    bss(2)=ys(2,0)
	    bss(3)=ys(4,0)
	    bss(4)=ys(6,0)
	    bss(5)=ys(7,0)
	   CASE(2)	!extérieur
	    bss(1)=frls(31)*ys(1,0)*ys(2,0)-ys(5,0)
	    bss(2)=ys(3,0)
	    bss(3)=2.d0/3.d0*frls(25)*ys(7,1)-ys(8,0)
	   END SELECT

c écritures
	   IF(id == 0)THEN
	    PRINT*,'eq_diff_rota4 dérivées ',nom_rot(i)
	   ELSE
	    PRINT*,'eq_diff_rota4 dérivées d',nom_rot(i)
	   ENDIF
c	   WRITE(*,2000)bss(:) ; WRITE(*,2000)bs(:)
	   WRITE(*,2000)(as(j,i,id),(bss(j)-bs(j))/dstor,j=1,nrot)
	   PAUSE'eq_diff_rot test dérivées'
	  ENDDO
	 ENDDO
	 DEALLOCATE(bss,coeff,frls,ys) ; deriv=.FALSE.
	ENDIF

c------------------fin des tests --------------------------

	RETURN

	END SUBROUTINE eq_diff_rota4
