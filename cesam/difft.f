         	
c****************************************************************
	 
	SUBROUTINE difft(l,melange,p,t,r,ro,drox,kap,dkapr,dkapx,deff,gradad,
	1 gradrad,m,y,d,dd)

c subroutine private du module mod_evol	 	
c subroutine générique de calcul des coefficients de diffusion turbulente

c entrées :
c	melange=.TRUE.: on est dans une ZC
c	p : pression
c	t : température
c	ro, drox : densité et dérivée / X (par gramme)
c	kap, dkapx: opacité et dérivée / X (par gramme)
c	dkapr : dérivée / ro
c	deff : diffusivité turbulente due à la rotation
c	gradad, gradrad : gradients adiabatique et radiatif
c	m : masse
c	y : composition chimique / MOLE

c sorties :
c	d, dd : coefficients d_ij de d x_j / d m et dérivées / x_k par MOLE

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : d_conv, iLi7, nchim, nom_difft
	USE mod_kind
	USE mod_variables, ONLY : tot_conv

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nchim,0:1) :: y	
	REAL (kind=dp), INTENT(in) :: deff, dkapr, dkapx, drox, gradad, gradrad,
	1 kap, l, m, p, r, ro, t	
	LOGICAL, INTENT(in) :: melange			
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d

	REAL (kind=dp) :: nu_tac
		
	INTEGER :: i
	 
c-------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c mélange total
	IF(tot_conv)THEN
	 DO i=1,nchim
	  d(i,i)=d_conv
	 ENDDO

c mélange partiel
	ELSE

c diffusion turbulente		 
	 SELECT CASE(nom_difft)
	 CASE('difft_nu')
	  CALL difft_nu(melange,t,ro,drox,kap,dkapr,dkapx,deff,d,dd)

c ajout de la diffusion dans la tachocline selon Castro & al. A&A 463,755, 2007	  
	 CASE('difft_nut_ext')
	  CALL difft_nu(melange,t,ro,drox,kap,dkapr,dkapx,deff,d,dd)
	  nu_tac=tacho_ext(r)
	  DO i=1,nchim
	   d(i,i)=d(i,i)+nu_tac
	  ENDDO
	  
	 CASE('difft_nut_full')
	  CALL difft_nu(melange,t,ro,drox,kap,dkapr,dkapx,deff,d,dd)
	  nu_tac=tacho_full(r)
	  DO i=1,nchim
	   d(i,i)=d(i,i)+nu_tac
	  ENDDO

	 CASE('difft_gab')
	  CALL difft_gab(melange,t,deff,d)
	 CASE('difft_smc')
	  CALL difft_smc(melange,deff,l,m,p,r,t,y,d,dd)	 
	 CASE('difft_sun')
	  CALL difft_sun(melange,deff,gradad,gradrad,m,d)	 
	 CASE DEFAULT
	  PRINT*,'routine de diffusion turbulente inconnue: ',nom_difft
	  PRINT*,'routines connues: difft_gab, difft_nu, difft_smc, difft_sun'
	  PRINT*,'ARRÊT' ; STOP
	 END SELECT
	ENDIF
	
	RETURN
	
	CONTAINS
	
	 INCLUDE 'tacho_ext.f'
	 INCLUDE 'tacho_full.f'
	 
	END SUBROUTINE difft
