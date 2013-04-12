	
c***********************************************************************

	PROGRAM test_tabul_reac

c test pour la tabulation des reactions nucleaires
c permet la comparaison avec les taux tabules
c utilise le fichier de donnees test.don

c Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c CESAM2k

c--------------------------------------------------------------------

	USE mod_donnees, ONLY : lit_nl, nchim, nom_fich2, w_rot
	USE mod_kind
	USE mod_nuc, ONLY : knot_temp, m_temp, nreac, nuc, n_temp, taux_reac,
	1 temp, ttemp, t_sup
	USE mod_numerique, ONLY : bsp1dn
		
	IMPLICIT NONE

!	REAL (kind=dp), DIMENSION(0,0)	:: jac
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: jac
!	REAL (kind=dp), DIMENSION(0) :: comp, dcomp, ex	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: ri, dri, comp, dcomp, ex
	REAL (kind=dp), DIMENSION(5) :: epsilon	
	REAL (kind=dp) :: t, ro, et, ero, hhe, be7e, b8e, n13e, o15e, f17e
	
	INTEGER :: l=1, fait

c---------------------------------------------------------------
	
2000	FORMAT("temperature=",8es10.3)
2001	FORMAT("taux=",8es10.3)

c lecture du fichier de donnees
	nom_fich2='test' ; CALL lit_nl(w_rot)
	
c initialisations	
!	t= 1.589E+08 ; ro= 7.850E+03
	t= 1.5E+07 ; ro= 1.50E+02
        fait=0
        ALLOCATE(comp(0),dcomp(0),jac(0,0),ex(0))
	CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

        print*,'nchim = ',nchim
        DEALLOCATE(comp) ; ALLOCATE(comp(nchim)); fait=1
	CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)


        fait=3
        comp(1)=0.7
        comp(2)=0.937d-5
        comp(3)=0.07    
	CALL nuc(t,ro,comp,dcomp,jac,.FALSE.,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)	

        print*,'Les neutrinos hh',hhe

c interpolation en ln t		
	PRINT*,'ATTENTION: les facteurs 1/2!, ou 1/3!'//
     s  'sont inclus dans les taux'
	ALLOCATE(dri(nreac),ri(nreac))

	print*,'eps = ',epsilon
	print*,'composition = ',comp
	t=1.d6
	WRITE(*,2000)t
	t=LOG(t)	
	CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	WRITE(*,2001)EXP(ri)	 

	t=1.d7
         print*,' '
	WRITE(*,2000)t
	t=LOG(t)	
	CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)
	WRITE(*,2001)EXP(ri)	 

	t=15.d6
         print*,' '
	WRITE(*,2000)t
	t=LOG(t)	
	CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	WRITE(*,2001)EXP(ri)	 

	t=25.d6
	IF(t_sup >= t)THEN
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	 WRITE(*,2001)EXP(ri)
	ENDIF 
	 
	t=50.d6	
	IF(t_sup >= t)THEN
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)  
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)
	 WRITE(*,2001)EXP(ri)
	ENDIF
	 
	t=1.d8
	IF(t_sup >= t)THEN	  
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)
	 WRITE(*,2001)EXP(ri)
	ENDIF
		 
	t=1.5d8
	IF(t_sup >= t)THEN	
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	 WRITE(*,2001)EXP(ri)
	ENDIF
	
	t=3.d8
	IF(t_sup >= t)THEN	
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	 WRITE(*,2001)EXP(ri)
	ENDIF

	t=5.0d8
	IF(t_sup >= t)THEN	
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	 WRITE(*,2001)EXP(ri)
	ENDIF
	
	t=7.d8
	IF(t_sup >= t)THEN	
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	 WRITE(*,2001)EXP(ri)
	ENDIF

	t=9.d8
	IF(t_sup >= t)THEN	
         print*,' '
	 WRITE(*,2000)t
	 t=LOG(t)
	 CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,t,l,ri,dri)	    
	 WRITE(*,2001)EXP(ri)
	ENDIF
	
	STOP	 	 	
	
	END PROGRAM test_tabul_reac
