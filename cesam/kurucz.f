
c****************************************************************
C
      SUBROUTINE kurucz(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
     1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c     routine private du module mod_tdetau

c     loi t(tau) de C van Veer issue du modèle Kurucz 5777, 4.47, X=.7, Z=.02
c     mis sous la forme t=Teff(3/4(tau+q(tau)))**1/4, l/Hp=1.78

c     Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c     YL modifications de k5777.f pour lire n'importe quelle loi

c entrées :
c     tau : profondeur optique Rosseland
c     teff : température effective
c     grav : gravité

c sorties :
c     t : température
c     dtsd* : dérivées t/ tau, teff, grav
c     dtsd* : dérivées t/ tau, teff, grav
c     dro_** : dérivées ro_ext/ teff, grav
c     f_tau, df_tau, df_tau2 : f, d f / d tau, d2 f / d2 tau

c---------------------------------------------------------------------

      USE mod_kind
      USE mod_numerique, ONLY : bsp1ddn, no_croiss
      
      IMPLICIT NONE

      REAL (kind=dp), INTENT(in) :: tau, teff
      REAL (kind=dp), INTENT(out) :: df_tau, dro_grav, dro_teff, dtsdg,
     1 dtsdtau, dtsdteff, f_tau, d2f_tau, ro_ext, t

cYLD   modifications pour lecture d'une loi T(tau) dans un fichier
cYLD      INTEGER, PARAMETER :: mx=4, nx=143
      INTEGER, PARAMETER :: mx=4, nx=285   !YLD
      INTEGER :: nxloi                     !YLD
      REAL (kind=dp), SAVE, DIMENSION(nx) :: taus    !YLD
      REAL (kind=dp), SAVE, DIMENSION(nx) :: t1      !YLD


      REAL (kind=dp), SAVE, DIMENSION(nx+2*mx) :: taut
      REAL (kind=dp), SAVE, DIMENSION(1,nx) :: ts
      REAL (kind=dp), DIMENSION(mx,0:mx-1) :: fx
c     REAL (kind=dp), PARAMETER :: teffc=5777.d0, ro_ext0=3.55d-9        !YLD
      REAL (kind=dp), SAVE :: teffc, ro_ext0                       !YLD

      REAL (kind=dp), SAVE :: cte1
      INTEGER, SAVE :: knot, l=1
      LOGICAL, SAVE :: init=.TRUE.
      INTEGER :: i
c---------------------------------------------------------------------

      IF(init)THEN
cYLD    lecture loi T(tau)  
       OPEN(unit=30,form='formatted',status='old',file='Ttau')    !YLD
       READ(30,*); read(30,*) ro_ext0 ;	read(30,*) nxloi           !YLD
       READ(30,*) teffc ; READ(30,*) tau_min	                  !YLD
       DO i=1,nxloi ; READ(30,*) t1(i),taus(i) ; ENDDO            !YLD
       CLOSE(UNIT=30)	                                          !YLD
cYLD    fin lecture
       init=.FALSE. ; cte1=(3.d0/4.d0)**0.25d0
cYLD       ts=RESHAPE(t1,(/1,nx/)) ; ts=(ts/teffc)**4*4.d0/3.d0
       ts=RESHAPE(t1,(/1,nxloi/)) ; ts=(ts/teffc)**4*4.d0/3.d0    !YLD
       CALL bsp1ddn(1,ts,taus,taut,nxloi,mx,knot,.FALSE.,taus(1),l,fx)   !YLD
        IF(no_croiss)THEN
        PRINT*,'Arrêt 1 dans kurucz' ; STOP
       ENDIF
cYLD       rad=.FALSE. ; tau_min=1.d-4              !YLD
       rad=.FALSE.
       WRITE(*,1)tau_min,ro_ext0 ; WRITE(2,1)tau_min,ro_ext0
1      FORMAT(/,1x,'loi t(tau,teff,grav), k5777, non purement radiative',/,
     1  'tau_min=',1pd10.3,' ro_ext=',1pd10.3,/)       
      ENDIF

      CALL bsp1ddn(1,ts,taus,taut,nxloi,mx,knot,.TRUE.,tau,l,fx)    !YLD
      f_tau=fx(1,0) ; df_tau=fx(1,1) ; d2f_tau=fx(1,2)
      t=teff*cte1*f_tau**0.25d0 ; dtsdtau=t/4.d0/f_tau*df_tau
      dtsdteff=t/teff ; dtsdg=0.d0 ; ro_ext=ro_ext0 ; dro_grav=0.d0
      dro_teff=0.d0

	RETURN
	
	END SUBROUTINE kurucz
