
c***********************************************************************
	
	SUBROUTINE eq_diff_rota3(dt,fait,nu,y,as,bs)	
	
c routine private du module mod_evol

c formation des équations à résoudre dans resout_rota
c pour la rotation par la méthode de collocation
c formalisme de Talon & Zahn 1997, Krot=3

c Auteur: P.Morel, Département Cassiopée, O.C.A.
c CESAM2k

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
c	y(3,*): Théta
c	y(4,*): Lambda
c	y(5,*): Psi
c	y(6,*): Tau
c	y(7,*): Upsilon

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : nom_rot, nrl, nrot
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
	REAL (kind=dp), DIMENSION(nrl) :: frl	
	REAL (kind=dp), PARAMETER :: dx=1.d-5, unpdx=1.d0+dx	
	REAL (kind=dp) :: dstor, stor, stor0	
	
	INTEGER :: i, id, j

	LOGICAL :: deriv=.FALSE., zc
		
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c le nombre d'inconnues est nrot=7
c y(1:nrot,der): dérivée d'ordre der=0, der=1 pour dérivée première

c mises à 0
	as=0.d0 ; bs=0.d0
	
c ZC?
	zc=lmix(nu)

c les équations aux points de collocation, fait=0
	SELECT CASE(fait)
	CASE(0)
	
c les coefficients
	 CALL coeff_rota3(dt,nu,y,frl)

c dans ZC
	 IF(zc)THEN
	 	 	  
c homogènisation de la vitesse angulaire par diffusion, Dv >> 1
	  bs(1)=frl(1)*y(1,0)-frl(2)-y(6,1)
	  as(1,1,0)=frl(1) ; as(1,6,1)=-1.d0
	   	 
c dans ZC, variation de moment cinétique --> U 
	  bs(2)=frl(3)*y(1,0)*y(2,0)-y(6,0)
	  as(2,1,0)=frl(3) ; as(2,2,0)=frl(3) ; as(2,6,0)=-1.d0
	 
c Théta=0
	  bs(3)=y(3,0) ; as(3,3,0)=1.d0

c Lambda=0
	  bs(4)=y(4,0) ; as(4,4,0)=1.d0
	  
c Psi = 0
	  bs(5)=y(5,0) ; as(5,5,0)=1.d0
	  
c définition de Tau (Flux de moment cinétique)	 
	  bs(6)=frl(4)*y(1,1)-y(6,0)
	  as(6,1,1)=frl(4) ; as(6,6,0)=-1.d0

c Upsilon=0	 
	  bs(7)=y(7,0) ; as(7,7,0)=1.d0
	 	 
	 ELSE
	 	 	 
c dans les ZR advection + diffusion du moment cinétique
	  bs(1)=frl(1)*y(1,0)-frl(2)-y(6,1)
	  as(1,1,0)=frl(1) ; as(1,6,1)=-1.d0
	 
c définition de U
	  bs(2)=(frl(5)*y(1,0)**2-frl(7))*y(1,0)**2+frl(10)+frl(6)*y(2,0)	
	1 +frl(11)*y(3,0)+frl(15)*y(4,0)+frl(12)*y(5,0)-frl(13)*y(5,1)-y(7,1)	
	  as(2,1,0)=4.d0*frl(5)*y(1,0)**3-2.d0*frl(7)*y(1,0)
	  as(2,2,0)=frl(6) ; as(2,3,0)=frl(11) ; as(2,4,0)=frl(15)
	  as(2,5,0)=frl(12) ; as(2,5,1)=-frl(13) ; as(2,7,1)=-1.d0
	  	 	  
c définition de Théta 
	  bs(3)=frl(19)*y(1,0)*y(1,1)-y(3,0)
	  as(3,1,0)=frl(19)*y(1,1) ; as(3,1,1)=frl(19)*y(1,0) ; as(3,3,0)=-1.d0
	  	 
c évolution de Lambda
	  bs(4)=frl(22)+frl(21)*y(2,0)-frl(20)*y(4,0)
	  as(4,2,0)=frl(21) ; as(4,4,0)=-frl(20)
	  
c définition de Psi
	  bs(5)=frl(24)*y(4,0)-frl(23)*y(3,0)-y(5,0)
	  as(5,3,0)=-frl(23)
	  as(5,4,0)=frl(24)
	  as(5,5,0)=-1.d0
	  
c définition de Tau (Flux de moment cinétique)	 
	  bs(6)=frl(3)*y(1,0)*y(2,0)+frl(4)*y(1,1)-y(6,0)
	  as(6,1,0)=frl(3)*y(2,0) ; as(6,1,1)=frl(4)	  	  
	  as(6,2,0)=frl(3)*y(1,0) ; as(6,6,0)=-1.d0 	 
	  
c définition de Upsilon	
	  bs(7)=frl(17)*y(3,0)+frl(18)*y(4,0)+frl(16)*y(5,0)+frl(14)*y(5,1)
	1 -y(7,0)
	  as(7,3,0)=frl(17) ; as(7,4,0)=frl(18)	; as(7,5,0)=frl(16)	
	  as(7,5,1)=frl(14) ; as(7,7,0)=-1.d0	  
	  	
	 ENDIF
	 
c~~~~~~~~~~Conditions limites~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c fait=0 : équations aux points de collocation
c fait=1 : conditions au centre
c fait=2 : conditions à la limite externe toujours convectif

c conditions limites au centre, fait=1, 4 limites
	CASE(1)
		 
c condition dOmega=0
	 bs(1)=y(1,1) ; as(1,1,1)=1.d0
	 	 
c condition U=0
	 bs(2)=y(2,0) ; as (2,2,0)=1.d0
	 
c condition Lambda=0
	 bs(3)=y(4,0) ; as(3,4,0)=1.d0	
	   	  
c condition Upsilon=0
	 bs(4)=y(7,0) ; as(4,7,0)=1.d0	  

c conditions limites externes toujours ZC, fait=2, 3 limites
	CASE(2)		

c condition U = 0
	 bs(1)=y(2,0) ; as(1,2,0)=1.d0

c condition Théta = 0
	 bs(2)=y(3,0) ; as(2,3,0)=1.d0
	 
c condition Psi = 0
	 bs(3)=y(5,0) ; as(3,5,0)=1.d0	 	 
	END SELECT	 
		 
c-----------------vérification des dérivées-------------------

c	deriv=.TRUE.
c	deriv=(nu >= 6.E-01 .AND. nu <= 6.03E-01)	!ZR
c	1 .OR. (nu >= 1.E-01 .AND. nu <= 1.03E-01)	!ZC
c	2 .OR. (fait > 0)				!limite
c	deriv=zc
c	deriv=.NOT.zc
c	deriv= fait > 0
c	deriv= nu > 0.8d0 .AND. nu < 0.9d0
c	deriv= nu > 0.141d0 .AND. nu < 0.142d0
c	deriv= deriv .OR. nu > 1.155d0 .AND. nu < 1.157d0
c	deriv= nu > mrot(convf(1)-1) .AND. nu < mrot(convf(1)+1)
c	deriv= deriv .OR. nu > mrot(convd(2)-1) .AND. nu < mrot(convd(2)+1)
c	deriv= nu > mrot(convf(1)-1) .AND. nu < mrot(convf(1)+5)

	IF(deriv)THEN
	 IF(zc)THEN 
	  PRINT*,'convectif'
	 ELSE
	  PRINT*,'radiatif'
	 ENDIF
	 PRINT*,'nu=',nu
	 PRINT*,'y(:,0)' ; WRITE(*,2000)y(:,0) 
	 PRINT*,'y(:,1)' ; WRITE(*,2000)y(:,1)
	 PRINT*,'frl' ; WRITE(*,2000)frl
	 PRINT*,'bs' ; WRITE(*,2000)bs
	 ALLOCATE(bss(nrot),frls(nrl),ys(nrot,0:1))	   
	 DO id=0,1
	  DO i=1,nrot
	   ys=y
	   stor0=ys(i,id) ; stor=stor0*unpdx ; IF(stor == 0.d0)stor=dx
	   dstor=stor-stor0 ; ys(i,id)=stor
	   CALL coeff_rota3(dt,nu,ys,frls)
	   IF(zc)THEN	 
	    bss(1)=frls(1)*ys(1,0)-frls(2)-ys(6,1)
	    bss(2)=frls(3)*ys(1,0)*ys(2,0)-ys(6,0)
	    bss(3)=ys(3,0)
	    bss(4)=ys(4,0)
	    bss(5)=ys(5,0)
	    bss(6)=frls(4)*ys(1,1)-ys(6,0)
	    bss(7)=ys(7,0)
	   ELSE
	    bss(1)=frls(1)*ys(1,0)-frls(2)-ys(6,1)
	    bss(2)=(frls(5)*ys(1,0)**2-frls(7))*ys(1,0)**2+frls(10)
	1   +frls(6)*ys(2,0)+frls(11)*ys(3,0)+frls(15)*ys(4,0)
	2   +frls(12)*ys(5,0)-frls(13)*ys(5,1)-ys(7,1)	
	    bss(3)=frls(19)*ys(1,0)*ys(1,1)-ys(3,0)
	    bss(4)=frls(22)+frls(21)*ys(2,0)-frls(20)*ys(4,0)
	    bss(5)=frls(24)*ys(4,0)-frls(23)*ys(3,0)-ys(5,0)
	    bss(6)=frls(3)*ys(1,0)*ys(2,0)+frls(4)*ys(1,1)-ys(6,0)
	    bss(7)=frls(17)*ys(3,0)+frls(18)*ys(4,0)+frls(16)*ys(5,0)
	1   +frls(14)*ys(5,1)-ys(7,0)	   
	   ENDIF	 
	   IF(id == 0)THEN
	    PRINT*,'eq_diff_rota3 dérivées ',nom_rot(i)
	   ELSE
	    PRINT*,'eq_diff_rota3 dérivées d',nom_rot(i)
	   ENDIF
c	   WRITE(*,2000)bss(:) ; WRITE(*,2000)bs(:)
	   WRITE(*,2000)(as(j,i,id),(bss(j)-bs(j))/dstor,j=1,nrot) 
	   PAUSE'eq_diff_rota3 test dérivées'
	  ENDDO
	 ENDDO
	 DEALLOCATE(bss,frls,ys) ; deriv=.FALSE.	   
	ENDIF
	 	
	RETURN

	END SUBROUTINE eq_diff_rota3
