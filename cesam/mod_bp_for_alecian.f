      Module MOD_BP_FOR_ALECIAN
      USE MOD_DONNEES, only: msol, lsol, rsol, g, nomch => nom_chemin, 
     +                       zi, pi, sigma, ihe4, i_ex
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: alecian2        !!!!!!! , id_ppxi
C
      Integer, Parameter :: DP = KIND(1.d0)    
C
C grav devrait etre egal a g_grav d'Alecian
C         et donc plus besoin pour Geaoges de le calculer
C    Dans ce cas, dans   modul_g_rad2    mtot ne sert que lors de l'initialisation
C                                        MAIS     mass ne sert plus
C
C-Cesam      Real(DP), Public :: grav, teff, mass
C
      Real(DP) :: g_grav , mass    ! Cesam , Rmq : mass ne sert pas ..... 
C
C LES ANCIENS COMMONS DE GEORGES
C
	integer, parameter                       :: pnchim = 13 , pzi=26
C
      integer, save  :: igrad = 1 , tgrad = 2      ! on l'impose en interne

	integer, parameter                      :: nion_ppxi=30
	real(DP), dimension(0:nion_ppxi,pnchim) :: phistar,psistar
	real(DP), dimension(0:nion_ppxi,pnchim) :: xistar,alphstar
	real(DP), dimension(0:nion_ppxi,pnchim) :: acstar,betstar
C
	integer, dimension(pnchim)                     :: id_ppxi
	integer, dimension(pnchim,pnchim)              :: isotops
C
	character(Len=2), dimension(pnchim)           :: el_svp     ! ne sert pas ?????
	real(DP), dimension(pnchim)                   :: C_sol
C
C variables common
C tableaux ions
C
	integer, dimension(0:pzi,pnchim)        :: niv_flag,niv_nb,niv_z
	integer, dimension(0:pzi,pnchim), SAVE  :: rar_flag
	real(DP), dimension(0:pzi,pnchim)       :: el_pot
C
C tableaux ions,niveaux: tableau(niveau,deg_ionisation,numero_element_cesam)
C
	integer, dimension(30,0:pzi,pnchim)  :: niv_g
	real(DP), dimension(30,0:pzi,pnchim) :: niv_e,niv_q

      CONTAINS
C===================================================================
      SUBROUTINE alecian2 (lum,ray,t,nel,grav,ychim,ioni,g_rad,dg_rad)

c     routine private du module mod_evol
c     Version 0.0 (2004/09/22)
c     cette subroutine calcule l'accélération radiative et ses dérivées
c     par rapport aux l'abondances. Elle doit être appelée pour chaque
c     couche, et à chaque fois que la diffusion est calculée.
c
c     Adaptation a CESAM2k par B. Pichon, OCA
c
c entrées:
c     ray: rayon
c     t : température
c     nel : nombre électrons par volume,
c     ychim : composition chimique/mole
c     ioni : taux d'ionisation
c     grav : gravité
c
c sorties
c     g_rad(i) : accélérations radiatives la gravité est grav+g_rad sur i
c     dg_rad(i,j) : dérivée de g_rad(i) / ion j
c
c--------------------------------------------------------------
c
      USE mod_donnees, ONLY: nchim, nucleo
      USE mod_variables, ONLY : mtot => mstar , rtot => rstar

      IMPLICIT NONE

      Real(DP), INTENT(in) :: lum, ray, t, nel, grav
      Real(DP), INTENT(in), DIMENSION(:) :: ychim
      Real(DP), INTENT(in), DIMENSION(0:,:) :: ioni
      Real(DP), INTENT(out), DIMENSION(:,:) :: dg_rad
      Real(DP), INTENT(out), DIMENSION(:) :: g_rad

      Logical, Save :: init_acc_rad2 = .False.

      If ( pnchim < nchim ) Then
         Write(*,*) " MAUVAISE VALEUR DE pnchim = ", pnchim
         STOP " NOS SOMMES DANS acc_rad2 " 
      End If

c      If ( .NOT. f_rad ) Return
      g_rad = 0.d0       ! securite
      dg_rad = 0.d0       ! securite     
      If ( .NOT. init_acc_rad2 ) Then
         CALL modul_ppxi2(mtot*msol,nchim,nucleo)
         init_acc_rad2 = .True.
      Else
         If ( grav >= 0.d0 ) Return    
         call modul_g_rad2_bp(nchim,nucleo,lum*lsol,ray*rsol,t,nel,
     +                        grav,ychim(:nchim),ioni,g_rad,dg_rad)
C
C Il samble qu'il faut changer de signe (convention de signe de g_rad
C   differente entre P.Morel et G.Alecian
C
         g_rad  = - g_rad
         dg_rad = - dg_rad
C
C-BP
C-BP ON DECIDE DE PRENDRE : 
C-BP
      dg_rad = 0.
C-BP
      g_rad(1)      = 0.    ! H1
      g_rad(ihe4-1) = 0.    ! He3
      g_rad(ihe4)   = 0.    ! He4
CCC      g_rad(i_ex)   = 0.    ! Ex
C-BP

C
C         Write(99,*) " lum,ray,t,grav == ",REAL(lum), REAL(ray), REAL(t), REAL(grav)
C         Write(99,*) " YCHIM == ",REAL(ychim)
C         Write(99,*) " G_RAD == ",REAL(g_rad) 
C         Write(99,"(1x,75('-'),1x)")
      End If
C 
c-------------------------------------------------------------------- 
      END SUBROUTINE alecian2
!=======================================================================
! subroutines et function a mettre dans CESAM, a partir d'ici !!!
!=======================================================================
      Include "from_alecian.f"
!=======================================================================
!=======================================================================
      End Module MOD_BP_FOR_ALECIAN
