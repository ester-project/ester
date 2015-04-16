C
c*****************************************************************
C
      MODULE mod_atm
C
C     Module regroupant les routines de restitution de l'atmosphère
C     La signification des variables est décrite au paragraphe F8 de la notice de CESAM2k
C
c     Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k avec Bernard PICHON, derniere modification 2008 OCT 
C
c-----------------------------------------------------------------
C
      USE mod_kind
      
      IMPLICIT NONE

      REAL(dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:):: bp_atm, bp_atm_t 
      REAL(dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) :: dlpp_atm,m_atm, 
     1 p_atm, pt_atm, r_atm, tau, t_atm, x_atm, x_atm_t, xt_atm, xt_atm_t
      REAL(dp), SAVE, PRIVATE :: delfim, delfip, ltaue, ltauf, tau_min
C
      INTEGER, SAVE, PRIVATE :: Ipgt, ne_atm, n23_atm, usl_atm
      LOGICAL, SAVE, PRIVATE :: rad
c
      CHARACTER (len=3), PARAMETER, PRIVATE, DIMENSION(8) :: nom_atm=
     1 (/ ' Pt', ' T ', ' R ', ' R*', ' M ', 'to*', ' to', ' Pg' /)
c
      PRIVATE
      PUBLIC :: atm, tdetau, thermo_atm
C
      CONTAINS
C
c----------------------------------------------------------------
C
      INCLUDE 'atm.f'
      INCLUDE 'coll_atm.f'
      INCLUDE 'edding.f'
      INCLUDE 'eq_atm.f'
      INCLUDE 'hopf.f'
      INCLUDE 'hsra.f'
      INCLUDE 'k5750.f'
      INCLUDE 'k5777.f'
      INCLUDE 'lim_atm.f'
      INCLUDE 'lim_gong1.f'
      INCLUDE 'lim_tau1.f'
      INCLUDE 'marcs.f'          ! ajout Bernard PICHON 2008
      INCLUDE 'roger.f'
      INCLUDE 'taueff.f'
      INCLUDE 'tdetau.f'
      INCLUDE 'thermo_atm.f'
      INCLUDE 'trho.f'
      INCLUDE 'trho_4000.f'      ! ajout Bernard PICHON 2008 OCT
C
      END MODULE mod_atm
