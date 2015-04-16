
c***********************************************************************
	
	SUBROUTINE eq_diff_rota(fait,nu,y,ad,as,bd,bs)	
	
c	routine private du module mod_evol

c	formation des équations à résoudre dans resout_rota
c	pour la rotation par la méthode des éléments finis Galerkin

c	Auteur: P.Morel, Département Cassiopée, O.C.A.
c	CESAM2k

c entrées
c	nu: point de calcul en m**2/3
c	y: solution
c	dt: pas temporel

c sorties	
c	as, ad: coefficients de la spline/dérivée
c	bs, bd: seconds membres

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : m_rl, nom_rot, nrot
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
				
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,0:) :: y
	REAL (kind=dp), INTENT(in) :: nu
	INTEGER, INTENT(in) :: fait	
	REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: ad, as
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: bd, bs
	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: yd
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: bds, bss

	REAL (kind=dp), DIMENSION(nrl) :: dfrl, frl
	REAL (kind=dp), DIMENSION(1) :: dfm, fm
	REAL (kind=dp), PARAMETER :: dx=1.d-5, unpdx=1.d0+dx
	
	REAL (kind=dp) :: dstor, stor, stor0	
	
	INTEGER :: i, id, j, l=1
	
	LOGICAL :: deriv=.FALSE.
		
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	le nombre d'inconnues est nrot

c	y(1:nrot,der): dérivée d'ordre der=0, der=1
c	(der=0 pour la fonction, der=1 pour dérivée première)

c	bs(i) : coefficient de S dans < bs . S > = 0	
c	as(i,j,k) : dérivée bs(i) / y(j,k)
c	i=1,...,nv, j=1,...,nv, k=0,1
	
c	bd(i) : coefficient de dS dans < bd . dS > = 0
c	ad(i,j,k) : dérivée bd(i) / y(j,k)
c	i=1,...,nv, j=1,...,nv, k=0,1	

c	mises à 0
	ad=0.d0 ; as=0.d0 ; bd=0.d0 ; bs=0.d0

	SELECT CASE(fait)
	
	CASE(0)
c limites extérieures et intérieures		  
	 bs(1)=y(1,1) ; as(1,1,1)=1.d0 !dOmega/dnu=0
	 bs(2)=y(2,0) ; as(2,2,0)=1.d0 !U=0  	  	  	  
	 bs(3)=y(3,0) ; as(3,3,0)=1.d0 !theta=0
	 bs(4)=y(4,0) ; as(4,4,0)=1.d0 !Lambda=0
	 bs(5)=y(5,0) ; as(5,5,0)=1.d0 !psi=0 

c les équations au point courant	 
	CASE(1)

c les coefficients
	 CALL bsp1dn(nrl,coef_rl,mrl,mrlt,n_rl,m_rl,knotrl,
	1 .TRUE.,nu,l,frl,dfrl)
	 IF(no_croiss)PRINT*,'Pb. en 1 dans eq_diff_rota'
	
c diffusion de Omega
	 bs(1)=frl(1)*y(1,0)-frl(2)
	 as(1,1,0)=frl(1)
	 	
	 bd(1)=frl(3)*y(1,0)*y(2,0)+frl(4)*y(1,1)
	 ad(1,1,0)=frl(3)*y(2,0)
	 ad(1,2,0)=frl(3)*y(1,0)
	 ad(1,1,1)=frl(4)
	 
c est-on dans une zone mélangée
	 CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1  nu,l,fm,dfm)
	 IF(fm(1) < 0.d0)THEN	 

c d2U/dnu2=0
	  bd(2)=y(2,1)
	  ad(2,2,1)=1.d0
 
c theta=0		 
	  bs(3)=y(3,0)
	  as(3,3,0)=1.d0
c Lambda=0
	  bs(4)=y(4,0)
	  as(4,4,0)=1.d0
c psi=0	 
	  bs(5)=y(5,0)
	  as(5,5,0)=1.d0

	 ELSE	 
	 
c zone non mélangée	 

c définition de U
	 bs(2)=(frl(5)*y(1,0)**2-frl(7)-frl(8)*y(2,0)
	1  +frl(9)*y(3,0))*y(1,0)**2+frl(10)+frl(6)*y(2,0)	
	2  +frl(11)*y(3,0)+frl(15)*y(4,0)+frl(12)*y(5,0)-frl(13)*y(5,1)
	 as(2,1,0)=2.d0*(2.d0*frl(5)*y(1,0)**2-frl(7)-frl(8)*y(2,0)
	1 +frl(9)*y(3,0))*y(1,0) 
	 as(2,2,0)=-frl(8)*y(1,0)**2+frl(6)	
	 as(2,3,0)=frl(9)*y(1,0)**2+frl(11)
	 as(2,4,0)=frl(15)
	 as(2,5,0)=frl(12)
	 as(2,5,1)=-frl(13)

	 bd(2)=frl(17)*y(3,0)-frl(18)*y(4,0)-frl(16)*y(5,0)+frl(14)*y(5,1)
	 ad(2,3,0)=frl(17)
	 ad(2,4,0)=-frl(18)	
	 ad(2,5,0)=-frl(16)	
	 ad(2,5,1)=frl(14)
	  
c définition de theta 
	 bs(3)=y(3,0)-frl(19)*y(1,0)*y(1,1)
	 as(3,1,0)=-frl(19)*y(1,1)
	 as(3,1,1)=-frl(19)*y(1,0) 	 
	 as(3,3,0)=1.d0

c évolution de Lambda
	  bs(4)=frl(20)*y(4,0)-frl(22)-frl(21)*y(2,0) 
	  as(4,2,0)=-frl(21)
	  as(4,4,0)=frl(20)
	
c définition de psi
	  bs(5)=y(5,0)-frl(24)*y(4,0)+frl(23)*y(3,0)
	  as(5,3,0)=frl(23)
	  as(5,4,0)=-frl(24)
	  as(5,5,0)=1.d0
	  
	 ENDIF
	 
	END SELECT 
	 
c-----------------vérification des dérivées-------------------

c          deriv=.TRUE.
c	 deriv= nu > 0.5d0 .AND. nu < 0.503d0	   
	 IF(deriv)THEN
	  IF(frl(25) >= 1.d12)THEN 
	   PRINT*,'convectif'
	  ELSE
	   PRINT*,'radiatif'
	  ENDIF
	  PRINT*,'nu=',nu
	  PRINT*,'y(:,0)' ; WRITE(*,2000)y(:,0) 
	  PRINT*,'y(:,1)' ; WRITE(*,2000)y(:,1)
	  PRINT*,'frl' ; WRITE(*,2000)frl
	  PRINT*,'bs' ; WRITE(*,2000)bs
	  PRINT*,'bd' ; WRITE(*,2000)bd
	  ALLOCATE(bss(nrot),bds(nrot),yd(nrot,0:1))
	  yd=y	   
	  DO id=0,1
	   DO i=1,nrot
	    stor0=yd(i,id) ; stor=stor0*unpdx ; IF(stor == 0.d0)stor=dx
	    dstor=stor-stor0 ; yd(i,id)=stor
	    bss=0.d0 ; bds=0.d0
	    IF(frl(23) >= 1.d12)THEN 
	     bss(1)=frl(1)*yd(1,0)-frl(2) 	
	     bds(1)=frl(3)*yd(1,0)*yd(2,0)
	     bss(2)=yd(1,1) 
	     bss(3)=yd(3,0)
	     bss(4)=yd(4,0)
	     bss(5)=yd(5,0)
	    ELSE	    
	     bss(1)=frl(1)*yd(1,0)-frl(2)
	     bds(1)=frl(3)*yd(1,0)*yd(2,0)+frl(4)*yd(1,1)
	     bss(2)=(frl(5)*yd(1,0)**2-frl(7)-frl(8)*yd(2,0)
	1      +frl(9)*yd(3,0))*yd(1,0)**2+frl(10)+frl(6)*yd(2,0)
	2      +frl(11)*yd(3,0)+frl(15)*yd(4,0)+frl(12)*yd(5,0)
	1      -frl(13)*yd(5,1)	
	     bds(2)=frl(17)*yd(3,0)-frl(18)*yd(4,0)-frl(16)*yd(5,0)
	1      +frl(14)*yd(5,1)
	     bss(3)=yd(3,0)-frl(19)*yd(1,0)*yd(1,1)
	     bss(4)=frl(20)*yd(4,0)-frl(22)-frl(21)*yd(2,0)
	     bss(5)=yd(5,0)-frl(24)*yd(4,0)+frl(23)*yd(3,0)
	    ENDIF		 
	    IF(id == 0)THEN
	     PRINT*,'eq_diff_rot dérivées ',nom_rot(i)
	    ELSE
	     PRINT*,'eq_diff_rot dérivées d',nom_rot(i)
	    ENDIF
c	      WRITE(*,2000)bss(:) ; WRITE(*,2000)bs(:)
c	      WRITE(*,2000)bds(:) ; WRITE(*,2000)bd(:)
	    WRITE(*,2000)(as(j,i,id),(bss(j)-bs(j))/dstor,j=1,nrot)
	    WRITE(*,2000)(ad(j,i,id),(bds(j)-bd(j))/dstor,j=1,nrot)
	    yd(i,id)=stor0    
	    PAUSE'eq_diff_rot test dérivées'
	   ENDDO
	  ENDDO
	  DEALLOCATE(bds,bss,yd) ; deriv=.FALSE.	   
	 ENDIF
	 	  	 	 	      
c------------------fin des tests --------------------------
	 	
	RETURN

	END SUBROUTINE eq_diff_rota
