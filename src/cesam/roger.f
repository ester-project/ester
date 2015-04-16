
c*********************************************************************

      SUBROUTINE roger(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
C
c     routine private du module mod_atm
C    
c     interface pour les lois T(tau) de R. cayrel avec [Fe/H]=0.0 
c     interface pour les lois T(tau) de R. cayrel avec [Fe/H]=0.0 
c     interface pour les lois T(tau) de R. cayrel avec [Fe/H]=0.0 
c     interface pour les lois T(tau) de R. cayrel avec [Fe/H]=0.0 
c     interface pour les lois T(tau) de R. cayrel avec [Fe/H]=0.0 
C
C     Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k, B.Pichon
C
c entrées :
c     tau : profondeur optique Rosseland
c     teff : température effective
c     grav : gravité
C
c sorties :
c     t : température
c     dtsd* : dérivées t/ tau, teff, grav
c     dtsd* : dérivées t/ tau, teff, grav
c     dro_** : dérivées ro_ext/ teff, grav
c     f_tau, df_tau, df_tau2 : f, d f / d tau, d2 f / d2 tau
C
c---------------------------------------------------------------------------
C
      USE mod_donnees, ONLY : langue, nom_chemin, nom_tdetau
      USE mod_kind
C
      IMPLICIT NONE
C
      REAL (kind=dp), INTENT(in) :: tau, teff, grav
      REAL (kind=dp), INTENT(out) :: t, dtsdtau, dtsdteff, dtsdg,
     1 ro_ext, dro_grav, dro_teff, f_tau, df_tau, d2f_tau
C
      REAL (kind=dp) :: logg, tef
      CHARACTER (len=4) :: fesh
      CHARACTER (len=80) :: nom
      LOGICAL, SAVE :: init=.TRUE.

c-----------------------------------------------------------------

      IF(init)THEN
         init=.FALSE.
         SELECT CASE ( nom_tdetau )
         CASE ( "roger+00" )
            nom=TRIM(nom_chemin)//'fesh+00.data' ; fesh='+00 '
         CASE ( "roger+02" ) 
            nom=TRIM(nom_chemin)//'fesh+02.data' ; fesh='+02 ' 
         CASE ( "roger-05" ) 
            nom=TRIM(nom_chemin)//'fesh-05.data' ; fesh='-05 '
         CASE ( "roger-10a" ) 
            nom=TRIM(nom_chemin)//'fesh-10a.data' ; fesh='-10a'
         CASE DEFAULT
            STOP " ROGER : PB - 0 "
         END SELECT
         rad=.FALSE.
         SELECT CASE(langue)	 
         CASE('english')
            WRITE(*,1001)fesh
            WRITE(2,1001)fesh
         CASE DEFAULT	  
            WRITE(*,1)fesh
            WRITE(2,1)fesh
         END SELECT	 
	 
	 
      ENDIF

1     FORMAT(/,'loi T(tau,teff,grav), roger***, non purement radiative',/,
     1 '[Fe/H]=',a4,'. Interpolations en Teff et gravité',/,
     2 'table étendue artificiellement a Teff = 3500K',/,
     3 '3500K < Teff < 7000K, 1 < Log10(g) < 5, tau_min=1.d-4',/)

1001  FORMAT(/,'T(tau,teff,grav) law, roger**, no purely radiative',/,
     1 '[Fe/H]=',a4,'. Interpolations in Teff and Log g',/,
     2 'artificial extend of the data to Teff = 3500K',/,
     3 '3500K < Teff < 7000K, 1 < Log10(g) < 5, tau_min=1.d-4',/)	  

      logg=LOG10(grav)        ! peut être changé dans trho
      tef=teff                ! peut être changé dans trho
      CALL trho(nom,tau,tef,logg,t,dtsdtau,dtsdteff,
     1 dtsdg,ro_ext,dro_teff,dro_grav,f_tau,df_tau,d2f_tau,"roger")
     
	RETURN

      END SUBROUTINE roger
