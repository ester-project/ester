      SUBROUTINE marcs(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
C
c     routine private du module mod_atm
C    
c     interface pour les lois T(tau) de MARCS
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
      INTEGER :: ind
      CHARACTER (len=5) :: fesh
      CHARACTER (len=80) :: nom
      LOGICAL, SAVE :: init=.TRUE. , new
C
c-----------------------------------------------------------------
C
      IF(init)THEN
         init=.FALSE.
         new = ( INDEX(nom_tdetau,"NEW") /= 0 )
C
C                                  MARCS_Z+0.00.data
         nom = TRIM(nom_chemin) // TRIM(nom_tdetau) // ".data"
         rad=.FALSE.
	 ind=INDEX(nom_tdetau,"Z")+1
	 fesh(1:5)=nom_tdetau(ind:ind+4)
         If ( new ) Then                   ! table étendue
            SELECT CASE(langue)	 
            CASE('english')
               WRITE(*,1002)fesh
               WRITE(2,1002)fesh
            CASE DEFAULT	  
               WRITE(*,2)fesh
               WRITE(2,2)fesh
            END SELECT
         Else                             ! table non étendue
            SELECT CASE(langue)	 
            CASE('english')
               WRITE(*,1001)fesh
               WRITE(2,1001)fesh
            CASE DEFAULT	  
               WRITE(*,1)fesh
               WRITE(2,1)fesh
            END SELECT
         End If
C
      ENDIF
C
      logg=LOG10(grav)        ! peut être changé dans trho
      tef=teff                ! peut être changé dans trho
      If ( new ) Then
         CALL trho_4000(nom,tau,tef,logg,t,dtsdtau,dtsdteff,
     1   dtsdg,ro_ext,dro_teff,dro_grav,f_tau,df_tau,d2f_tau,"MARCS")
      Else
         CALL trho(nom,tau,tef,logg,t,dtsdtau,dtsdteff,
     1   dtsdg,ro_ext,dro_teff,dro_grav,f_tau,df_tau,d2f_tau,"MARCS")
      End If
C
1     FORMAT(/,'loi T(tau,teff,grav), MARCS***, non purement radiative',/,
     1 '[Fe/H]=',a5,'. interpolations en Teff et gravité',/,
     2 'provenance de MARCS, extrapolation pour Teff < 4000 K ',/,
     3 '3000K <= Teff <= 8000K, 3 <= Log10(g) <= 5, tau_min=1.d-5',/)
C
1001  FORMAT(/,'T(tau,teff,grav) law, MARCS**, no purely radiative',/,
     1 '[Fe/H]=',a5,'. interpolations in Teff and Log g',/,
     2 'from MARCS, extrapoled of data for Teff > 3500 K',/,
     3 '3000K <= Teff <= 8000K, 3 <= Log10(g) <= 5, tau_min=1.d-5',/)	  
C
2     FORMAT(/,'loi T(tau,teff,grav), MARCS***, non purement radiative',/,
     1 '[Fe/H]=',a5', pour le moment . interpolations en Teff et gravité',/,
     2 'provenance de MARCS, avec des atmosphères de type p et/ou s',/,
     3 '2500K <= Teff <= 8000K, 0 <= Log10(g) <= 5, tau_min=1.d-5',/)
C
1002  FORMAT(/,'T(tau,teff,grav) law, MARCS**, no purely radiative',/,
     1 '[Fe/H]=',a5,' for now. interpolations in Teff and Log g',/,
     2 'from MARCS, with data of s and p type',/,
     3 '2500K <= Teff <= 8000K, 0 <= Log10(g) <= 5, tau_min=1.d-5',/)	  
C
      END SUBROUTINE marcs
