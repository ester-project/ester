
c******************************************************************

	MODULE mod_evol

c	module contenant les routines concernant l'initialisation
c	et l'�volution temporelle de la composition chimique avec
c	et sans diffusion � l'exclusion des routines relatives au nucl�aire

c fonctions private:
c	acc_rad : calcul des acc�l�rations radiatives
c	collision : calcul des taux de collision selon Paquette et al.
c	coulomb : logarithme de Coulomb
c	diffm : routine g�n�rique de diffusion microscopique
c	diffm_br : calcul des coefficients de diffusion microscopique selon
c	Burgers
c	diffm_mp : calcul des coefficients de diffusion microscopique selon
c	Michaud & Proffit
c	difft :  routine g�n�rique de diffusion turbulente
c	difft_nu : calcul des coefficients de diff. turbulente + radiative
c	diffus : routine d'int�gration des �quations de diffusion
c	eq_diff_chim : routine de calcul des coefficients des �quations de
c	diffusion
c	rk_imps : routine d'int�gration des �quations d'�volution de la
c	composition chimique sans diffusion	

c fonctions public:
c	evol : gestion de l'�volution temporelle de la composition chimique

c	Auteur: P.Morel, D�partement J.D. Cassini, O.C.A.,
c	CESAM2k

c NOTATIONS (h�las incoh�rentes) pour les d�veloppements en B-splines
c	n_ch : nombre VARIABLE de points �l�ment de mod_variables
c	nch : nombre FIXE de fonctions �l�ment de mod_donnees
c	m_ch : ordre FIXE des splines �l�ment de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES �l�ment de mod_variables

c eps : taux des �carts � droite pour les discontinuit�s
c mroct : base continue pour la rotation

c-------------------------------------------------------------------------

	USE mod_kind
	USE mod_nuc, ONLY : nuc
	
	REAL (kind=dp),	PARAMETER :: eps=1.d-3 , un_eps=1.d0-eps
	INTEGER, PARAMETER :: m_mix=2, nlim_rot=2, nmix=1, nrl=25	
	
	REAL (kind=dp), PRIVATE, SAVE, ALLOCATABLE, DIMENSION(:,:) :: coef_rl,
	1 fmix
	REAL (kind=dp), PRIVATE, SAVE, ALLOCATABLE, DIMENSION(:) :: xlim_rot	
	REAL (kind=dp), PRIVATE, SAVE, ALLOCATABLE, DIMENSION(:) :: mrl, mrlt,
	1 xcin, xcint, x_mix, x_mixt
	REAL (kind=dp), SAVE, PRIVATE :: mlim_vent=HUGE(1.d0),
	1 mlim_w=HUGE(1.d0)
	REAL (kind=dp), SAVE, PUBLIC :: nu_min
	
	INTEGER, SAVE, PRIVATE :: knotrl, knot_mix, n_rl, n_mix=0
	
	LOGICAL, SAVE, PUBLIC :: ini_rot=.FALSE.

	PRIVATE
	PUBLIC :: evol, v2_cgm

	CONTAINS

c********************************************************************

!	INCLUDE 'alecian1.f'
!	INCLUDE 'coeff_rota.f'	
!	INCLUDE 'collision.f'
!	INCLUDE 'coulomb.f'
!	INCLUDE 'diffm.f'
!	INCLUDE 'diffm_br.f'
!	INCLUDE 'diffm_mp.f'	
!	INCLUDE 'difft.f'
!	INCLUDE 'difft_nu.f'
!	INCLUDE 'difft_ventura.f'
!	INCLUDE 'diffus.f'	
!	INCLUDE 'diffw.f'
!	INCLUDE 'diffw_mpz.f'
!	INCLUDE 'eq_diff_chim.f'
!	INCLUDE 'eq_diff_rota.f'		
!	INCLUDE 'evol.f'
!	INCLUDE 'f_rad.f'	
!	INCLUDE 'pertw.f'
!	INCLUDE 'pertw_loc.f'
!	INCLUDE 'pertw_ptm.f'
!	INCLUDE 'pertw_sch.f'	
!	INCLUDE 'resout_chim.f'
!	INCLUDE 'resout_rota.f'		
	INCLUDE 'rk_imps.f'
!	INCLUDE 'v2_cgm.f'

	END MODULE mod_evol
