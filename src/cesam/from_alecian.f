
c*******************************************************************

c	subroutine PRIVATE du module mod_bp_for_alecian

	subroutine modul_ppxi2(mtot,nchim,nucleo)

!  Version 2.2 (11/10/2004)
!  Cette subroutine effectue la lecture des options, des donnees et
!  prepare les tableaux necessaires au calcul des
!  accelerations (modul_grad2).

!  Cette subroutine fait partie du module grcesam2.
!  Auteur:
!     Georges ALECIAN
!     DAEC - UMR 8631, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

	Implicit None

	Real(DP), Intent(In) :: mtot
      Integer, Intent(In) :: nchim
      Real(DP), Dimension(nchim), Intent(In) :: nucleo

!** declarations Morel mod. Pichon (BP)
	integer                                  :: i,j,k
	
!** declarations Alecian
! mes variables locales
	integer, parameter                      :: nsvp=300,nlist=100
	integer                                 :: iunit_tdb,iunit_svp,izn,ios
	integer                                 :: icompte,ic,icf,numion
	integer, dimension(nsvp)                :: np,iz,numa
	integer, dimension(3,nsvp)              :: n_cross
	integer, dimension(nlist)               :: mass_svp,numions
	integer, dimension(nsvp,nlist)          :: nps,izs
	character(1)                            :: b
	character(2), dimension(nsvp)           :: el
	character(4), dimension(nsvp,nlist)     :: listis
	character(4), dimension(nsvp)           :: listi
c	character(64), dimension(nlist)         :: listf
	CHARACTER (len=15), DIMENSION(nlist) :: listf
	real, dimension(nsvp,nlist)             :: phis,psis,xis,alphs,acs,bets
	real, dimension(nsvp)                   :: phi,psi,xi,alph,ac,bet,absol
	real(KIND=KIND(0D0))                    :: mass_st_sol
	real(KIND=KIND(0D0)), dimension(nlist)  :: X,F1,F2,F3,F4,F5,F6

! igrad  (si =0 le g_rad des elements non-identifies sera zero)
!        (si =1 le g_rad des elements non-identifies sera -gravite)
! tgrad  (si =0 le g_rad de tous les elements identifies sera zero)
!        (si =1 le g_rad de tous les elements identifies sera -gravite)
!        (si =2 le g_rad de tous les elements identifies sera calcule)
! np: nb de protons du noyau
! iz: degre d'ionization (0 pour neutre)
! izn: (non utilisee)
! listi: nom de l'ion (C1 = carbone neutre)
! numion: nb d'ions dans la base
!**************************
!
! unites de lecture
	iunit_tdb=70
	iunit_svp=71

! initialisation des tableaux
	id_ppxi  =  0
	isotops  =  0
	phistar  =  0.
	psistar  =  0.
	xistar   =  0.
	alphstar = -0.5
	acstar   =  1.
	betstar  =  0.

	if(igrad.eq.0) then
		print*,'Le g_rad des elements non-identifies sera zero (igrad=0)'
		write(10,*)'Le g_rad des elements non-identifies sera zero (igrad=0)'
	else if(igrad.eq.1) then
		print*,'Le g_rad des elements non-identifies sera -gravite (igrad=1)'
		write(10,*)'g_rad des elements non-identifies sera -gravite (igrad=1)'
	end if
	if(tgrad.eq.0) then
		print*,'Le g_rad des elements identifies sera zero (tgrad=0)'
		write(10,*) 'g_rad des elements identifies sera zero (tgrad=0)'
	else if(tgrad.eq.1) then
		print*,'Le g_rad des elements identifies sera -gravite (tgrad=1)'
		write(10,*) 'g_rad des elements identifies sera -gravite (tgrad=1)'
	else if(tgrad.eq.2) then
		print*,'g_rad des elements identifies sera calcule (tgrad=2)'
		write(10,*) 'g_rad des elements identifies sera calcule (tgrad=2)'
	end if

!***	

! Tableau des isotopes de CESAM, matrice (isotops) du genre:
!         H1  He3  He4  C12  C13  N14  N15  O16  O17  Fe56
!  H1     1    0    0    0    0    0    0    0    0    0
! He3     0    1    1    0    0    0    0    0    0    0
! He4     0    1    1    0    0    0    0    0    0    0
! C12     0    0    0    1    1    0    0    0    0    0
! C13     0    0    0    1    1    0    0    0    0    0
! N14     0    0    0    0    0    1    1    0    0    0
! N15     0    0    0    0    0    1    1    0    0    0
! O16     0    0    0    0    0    0    0    1    1    0
! O17     0    0    0    0    0    0    0    1    1    0
!Fe56     0    0    0    0    0    0    0    0    0    1

	j=1
	do while (j.le.nchim)
		icompte=0
		do k=1,nchim
			if(k.lt.j) then
				isotops(j,k)=0
			else if(k.eq.j) then
				isotops(j,k)=1
			else
				if(zi(j).eq.zi(k)) then
					isotops(j,k)=1
					icompte=icompte+1
				else
					isotops(j,k)=0
				end if
			end if
		end do
		ic=1
		do while (ic.le.icompte)
			j=j+1
			do k=1,nchim
				isotops(j,k)=isotops(j-1,k)
			end do
			ic=ic+1
		end do
		j=j+1
	end do

! lecture de la liste des fichiers contenant les tables phi, psi, etc.
	OPEN (UNIT = iunit_svp,
     +      FILE = TRIM(nomch)//'datai_gr2/liste_svp.data',
     +      STATUS = 'old')
	icf=0
	
	read(iunit_svp,'(a)') b
c	PRINT*,b	
	read(iunit_svp,'(a)') b
c	PRINT*,b
		
	do i=1,100
c		PRINT*,i,icf
		read(iunit_svp,'(a)',end=90) listf(i)
		icf=icf+1
c		PRINT*,i,icf,listf(i)
	end do
90	continue
	close(iunit_svp)

! determination des masses disponibles dans la base
	do i=1,icf
		read(listf(i),'(8x,i3)',iostat=ios) mass_svp(i)
		if(ios/=0) then
			print*,'erreur identification de mass_svp => stop'
			stop
		end if
	end do

	mass_st_sol=(mtot/msol)*100
	if((mass_st_sol.lt.mass_svp(1)).or.
     +	(mass_st_sol.gt.mass_svp(icf))) then
		print*
		print*,'Attention mtot hors limites dans modul_ppxi2:'
		print*,'mtot*100/msol=',mass_st_sol
		print*,'min*100=',mass_svp(1)
		print*,'max*100=',mass_svp(icf)
		print*
	end if

! lecture de la base des phi et psi
	
	do i=1,icf
	
c	PRINT*,i,nomch,listf(i)
c	PRINT*,TRIM(nomch)//'datai_gr2/'//listf(i)
	
		OPEN (UNIT = iunit_svp,
     +      FILE = TRIM(nomch)//'datai_gr2/'//TRIM(listf(i)),
     +      STATUS = 'old')

		read(iunit_svp,'(a1)') b
		numions(i)=0
		do j=1,1000
		   read(iunit_svp,'(3i4,1x,a4,6e12.3)',end=100) 
     +	        nps(j,i),izs(j,i),izn,listis(j,i),
     +	        phis(j,i),psis(j,i),xis(j,i),
     +	        alphs(j,i),acs(j,i),bets(j,i)

		   numions(i)=numions(i)+1

		end do
100		continue
		close(unit=iunit_svp)
	end do
10	format(i4,a1,i4,a1,e10.3,a1,e10.3,a1,i4,a1,a4)

! control des fichiers lus
	do i=2,icf
		if(numions(i-1).ne.numions(i)) then
			print*, 'Anomalie des numions => STOP'
			stop
		end if
		do j=1,numions(i)
			if(listis(j,i-1).ne.listis(j,i)) then
				print*, 'Anomalie des listis => STOP'
				stop
			end if
		end do
	end do	

! preparation
	numion=numions(1)
	n_cross= 0                              !tableau(3,nsvp)
	np(1:numion)     = nps(1:numion,1)      ! nb protons du noyau
	iz(1:numion)     = izs(1:numion,1)      ! charge ion
	listi(1:numion)  = listis(1:numion,1)
	
! interpolations
	do i=1,icf
		X(i)=float(mass_svp(i))
	end do
	do j=1,numion
		F1(1:icf) = phis(j,1:icf)
		F2(1:icf) = psis(j,1:icf)
		F3(1:icf) = xis(j,1:icf)
		F4(1:icf) = alphs(j,1:icf)
		F5(1:icf) = acs(j,1:icf)
		F6(1:icf) = bets(j,1:icf)

		phi(j)  = FT(mass_st_sol,100,X,F1)
		psi(j)  = FT(mass_st_sol,100,X,F2)
		xi(j)   = FT(mass_st_sol,100,X,F3)
		alph(j) = FT(mass_st_sol,100,X,F4)
		ac(j)   = FT(mass_st_sol,100,X,F5)
		bet(j)  = FT(mass_st_sol,100,X,F6)
	end do

! Cross-identification ppxi/cesam. On pose id_ppxi(j)=1 si au moins un ion OK
! On suppose que les proprietes atomiques (pour les g_rad) sont 
! les memes entre isotopes d'un meme element.
	do j=1,nchim
		do k=1,numion
			if(zi(j).eq.np(k)) then
				n_cross(1,k)       = iz(k)    !correspondance svp/cesam
				n_cross(2,k)       = j        !correspondance svp/cesam
				n_cross(3,k)       = n_cross(3,k) +1      ! nb.isotopes
				id_ppxi(j)         = 1
				phistar(iz(k),j)   = phi(k)
				psistar(iz(k),j)   = psi(k)
				xistar(iz(k),j)    = xi(k)
				alphstar(iz(k),j)  = alph(k)
				acstar(iz(k),j)    = ac(k)
				betstar(iz(k),j)   = bet(k)
			end if
		end do
	end do


! lecture des niveaux atomiques pour le calcul du g_cont (photoionisation)

	call modul_niv(listi,numion,n_cross,np,iunit_svp,nchim)
	
! lecture de la liste des abondances solaires.
	OPEN (UNIT = iunit_svp,
     +      FILE = TRIM(nomch)//'datai_gr2/solabnum.dat',
     +      STATUS = 'old')
	ic=0
	do k=1,nsvp
		read(iunit_svp,*,end=101) el(k),numa(k),absol(k)		
		ic=ic+1
	end do
101	continue
	close (iunit_svp)
! concentration par rapport a H pour le calcul de khi dans modul_grad_2
	print*, 'ic= ',ic
	do j=1,nchim
		do i=1,ic
			if(zi(j).eq.numa(i)) then
				el_svp(j) = el(i)
				C_sol(j)  = absol(i)
				exit
			end if
		end do
	end do
	do j=2,nchim
		C_sol(j)=C_sol(j)/C_sol(1)
	end do
	C_sol(1)=1.

	end subroutine modul_ppxi2
!=======================================================================
	subroutine modul_niv(listi,numion,n_cross,np,iunit_svp,nchim)
!  Version 1.1 (11/10/2004)
!  Cette subroutine lit les niveaux atomiques necessaires au calcul du g_cont.
!  Les donnees se trouvent dans le repertoire tbniv/, les niveaux atomiques
!  dans les fichiers tbion... sont des niveaux regroupes approximatifs, adaptes
!  pour le calcul des accelerations radiatives.

!  Cette subroutine fait partie du module grcesam2.
!  Version 2.1 (14/06/2004)
!  Auteur:
!     Georges ALECIAN
!     LUTH - UMR 8102, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20


	implicit none

! input
	integer, parameter                      :: nsvp=300
	integer                                 :: numion,iunit_svp,nchim
	integer, dimension(nsvp)                :: np
	integer, dimension(3,nsvp)              :: n_cross
	character(4), dimension(nsvp)           :: listi

! mes variables locales
	character(319)                          :: nfe
	integer                                 :: i,j,k,ii,jj,ios
	integer                                 :: dim
	character(4)                            :: nom
	integer, dimension(nsvp)                :: jz=0,ns=0,elflag
	integer, dimension(30,nsvp)             :: saq
	real, dimension(30,nsvp)                :: sae,sag
	real, dimension(nsvp)                   :: x=0.

	elflag=0         ! initialisation tableau

! On charge toutes les donnees atomiques disponibles
	do k=1,numion
		nfe=TRIM(nomch)//
     +	'datai_gr2/tbniv/tbion'//TRIM(ADJUSTL(listi(k)))//'.niv'
			
		OPEN (UNIT = iunit_svp,
     +		FILE = nfe,
     +            IOSTAT=ios,
     +		STATUS = 'old')

		if(ios/=0) CYCLE	

		elflag(k)=1      ! le fichier de l'ion existe

		read(iunit_svp,1075) nom                       ! nom de l'ion
		read(iunit_svp,1081) x(k)                      ! pot. ioni. (eV)
		read(iunit_svp,1073) jz(k)                     ! charge
		read(iunit_svp,1073) ns(k)                     ! nombre de niveaux
		read(iunit_svp,1045) dim                       ! inutilise
		read(iunit_svp,1059) (sae(i,k),i=1,ns(k))      ! energie (eV)
		read(iunit_svp,1057) (sag(i,k),i=1,ns(k))      ! poids stat.
		read(iunit_svp,1073) (saq(i,k),i=1,ns(k))      ! nb quantique pple
		close (iunit_svp)

	end do

! On remplit les tableaux utiles pour cesam aux memes positions que phistar,...
	do k=1,numion
		if(n_cross(2,k).ne.0) then    ! si l'element est calcule pour cesam
			ii = n_cross(1,k)
			jj = n_cross(2,k)
			do j=jj-n_cross(3,k)+1 , jj     ! n_cross(3,k): nb d'isotopes
				niv_flag(ii,j) = elflag(k)
				niv_nb(ii,j)   = ns(k)
				niv_z(ii,j)    = jz(k)
				el_pot(ii,j)   = x(k)
				do i=1,ns(k)
					niv_e(i,ii,j) = sae(i,k)
					niv_q(i,ii,j) = saq(i,k)
					niv_g(i,ii,j) = sag(i,k)
				end do
				call gr_gazrare(np(k),jz(k),rar_flag(ii,j))
			end do
		end if
	end do

1075	format(1x,a4)
1081	format(16f9.2)
1073	format(20i4)
1045	format(16i5)
1057	format(12f6.0)
1059	format(8f8.3)

	end subroutine modul_niv
!=======================================================================
      subroutine  modul_g_rad2_bp(nchim,nucleo,lum,ray,t,nel,g_grav,
     +                            ychim,ioni,g_rad,dg_rad)

! Changement de nom "modul_g_rad2"  en "modul_g_rad2_bp" (BP)
!     du fait de legers changements dans la liste d'arguments
!
!  Version 2.1 (14/06/2004)
!  Cette subroutine calcule l'acceleration radiative et sa derivee
!  par rapport a l'abondance. Elle doit etre appelee pour chaque couche,
!  et a chaque fois que la diffusion est calculee.

!  Cette subroutine fait partie du module grcesam2.
!  Version 2.1 (14/06/2004)
!  Auteur:
!     Georges ALECIAN
!     LUTH - UMR 8102, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

	implicit none

!** declarations Morel mod. Pichon (BP)
	integer, Intent(In)       :: nchim
	Real(DP), Intent(In)      :: lum,ray,t,nel,g_grav

	Real(DP), dimension(nchim), Intent(In)        :: nucleo, ychim
	Real(DP), dimension(0:,:), Intent(In)        :: ioni
	
!** declarations Alecian
!	OUTPUTS
      Real(DP), DIMENSION(nchim+1), INTENT(out) :: g_rad
      Real(DP), DIMENSION(nchim+1,nchim+1), INTENT(out) :: dg_rad

! pour le g_cont
	real(KIND=KIND(0D0)), dimension(0:pzi) :: fparti,theta
	real(KIND=KIND(0D0)), dimension(0:pzi,pnchim) :: popi

! mes variables locales
	integer                  :: iunit,iso,i,j,k,m,dim_1
	real,parameter           :: pondrare = 1.5
	real(KIND=KIND(0D0))     :: pondt
	real(KIND=KIND(0D0))	 :: gr_kj,dgr_kj,q,b,Ck_s,CX
	real(KIND=KIND(0D0))	 :: gc_kj,dgc_kj,bco

	real(KIND=KIND(0D0)), dimension(nchim)          :: N_chim,khi
	real(KIND=KIND(0D0)), dimension(0:pzi,nchim)    :: N_ion,C_ion
	real(KIND=KIND(0D0)), dimension(0:pzi,nchim)    :: pond_ion
	real(KIND=KIND(0D0)), dimension(30,0:pzi)        :: pniv

!**************************

! securites
! (BP) - debut
!
      If ( SIZE(ioni,Dim=2) /= nchim ) Then
         Write(*,*) " MAUVAISE DIMENSION (dim=2) DE 'ioni' "
         STOP " NOUS SOMMES DANS 'modul_g_rad2_bp' " 
      End If
      dim_1 = SIZE(ioni,Dim=1)

	popi=0.          ! tableau
      pond_ion = 0.    ! tableau   NEW NEW NEW NEW NEW NEW NEW 
!
! (BP) - fin
!
	g_rad  = 0.
	dg_rad = 0.

	if(ychim(1)*t*ray.lt.1.d-50) return     ! bizarre (BP)
	khi    = 1. !tableau
	gc_kj  = 0.
	dgc_kj = 0.

! module du vecteur gravite
C
! On considere que les isotopes ont leurs transitions atomiques
! presque aux memes frequences (bf comme bb).
! Cela implique que l'effet de saturation doit etre presque le meme
! pour tous les isotopes d'un meme element. Cette saturation est reglee par
! la concentration C_ion(k,j) (en pratique: abondance en nb par rapport a H).
! On veut donc que C_ion(k,j) soit le meme pour tous les isotopes d'un element
! et soit donne par la concentration totale pour l'element.
! Infos: Somme[ychim(j)*nucleo(j)]=1

	do j=1,nchim
		N_chim(j)=0.
		do iso=1,nchim
			N_chim(j)=N_chim(j)+ychim(iso)*isotops(j,iso)
		end do
	end do

	do j=1,nchim
! (BP) - debut
         If ( NINT(zi(j)) > dim_1 ) Then
            Write(*,*) " MAUVAISE DIMENSION (dim=1) DE 'ioni' "
            STOP " NOUS SOMMES DANS 'modul_g_rad2_bp' " 
         End If
! (BP) - fin
		do k=0,nint(zi(j))
			C_ion(k,j)=ioni(k,j)*N_chim(j)/N_chim(1)
		end do
		khi(j)=(N_chim(j)/N_chim(1)) / C_sol(j)
	end do

! coefficients de ponderation pour le g_rad total

	do j=1,nchim
		pondt=0.
		do k=0,nint(zi(j))
			pondt = pondt + ioni(k,j)*
     +		       (1.+rar_flag(k,j)*(pondrare-1.))
		end do
		do k=0,nint(zi(j))
			pond_ion(k,j) = ioni(k,j)*
     +		                (1.+rar_flag(k,j)*(pondrare-1.))/pondt
		end do
	end do

! Substitution par les populations relatives calculees par gr_popion.
! (Recalcul de: C_ion, pond_ion.)
! Cette sequence peut-etre enlevee si on prefere utiliser les populations ioni(k,j).
! Dans ce cas, activer la sequence suivante a la place:
!	do j=1,nchim
!		if(id_ppxi(j).eq.1) then
!			do k=0,nint(zi(j))
!				popi(k,j) = ioni(k,j)
!			end do
!		end if
!	end do
	do j=1,nchim
		if(id_ppxi(j).eq.1) then
			call gr_parti2(t,fparti,j)
			call gr_popion(t,nel,fparti,popi,j)
			do k=0,nint(zi(j))
				C_ion(k,j)=popi(k,j)*N_chim(j)/N_chim(1)
			end do
			pondt=0.
			do k=0,nint(zi(j))
				pondt = pondt + popi(k,j)*
     +		       (1. + rar_flag(k,j)*(pondrare-1.) )
			end do
			do k=0,nint(zi(j))
				pond_ion(k,j) = popi(k,j)*
     +		             (1. + rar_flag(k,j)*(pondrare-1.) )/pondt
			end do
		end if
	end do
! Fin substitution

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! constantes:
!                pi**2 * kbol**3 * echarg**2
!   5.57E-05 = ---------------------------------
!             2.* clight**4 * hpl**2 * me * amu
!
!
!               me * mp* clight
!   9.83E-23 = -----------------
!                    2. * pi
!
!   CX : fraction de masse de l'hydrogene
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	if(tgrad.eq.0) then
		do j=1,nchim
			 if(id_ppxi(j).eq.1) then
			 	g_rad(j) = 0.
			 else
				if(igrad.eq.0) then
				 	g_rad(j) = 0.
				else if(igrad.eq.1) then
				 	g_rad(j) = g_grav *
     +				(1-pond_ion(nint(zi(j)),j))    ! correction noyau nu,
! le noyau nu ne peut pas avoir g_rad = -g, on suppose que la photoionisation est 0
				end if
			 end if
		end do
	else if(tgrad.eq.1) then
		do j=1,nchim
			 if(id_ppxi(j).eq.1) then
			 	g_rad(j) = g_grav
			 else
				if(igrad.eq.0) then
				 	g_rad(j) = 0.
				 else if(igrad.eq.1) then
				 	g_rad(j) = g_grav *
     +				(1-pond_ion(nint(zi(j)),j))    ! correction noyau nu,
! le noyau nu ne peut pas avoir g_rad = -g, on suppose que la photoionisation est 0
				 end if
			 end if
		end do
	else if(tgrad.eq.2) then
		do j=1,nchim

			if(id_ppxi(j).eq.1) then
				call gr_bolm2(t,pniv,j)
				call gr_parti2(t,fparti,j)
				call gr_theta(t,theta,j)
			end if

			CX  = N_chim(1) * nucleo(1)
!
! MODIFICATION (BP) pour pouvoir utiliser le flux en tout point de l'etoile
!
            q   = 5.57E-05 * lum / ( 4.d0 * pi * sigma ) / t / ray**2 / nucleo(j)
			b   = 9.83E-23 * nel * (1./sqrt(t)) / CX

			do k=0,nint(zi(j))
			   if(id_ppxi(j).eq.1) then
			      if(psistar(k,j).gt.1.E-30) then

				   Ck_s   = b * psistar(k,j)**2

				   if(niv_z(k,j).gt.0) then
					   bco = 7.16e-26 * lum / ( 4.d0 * pi * sigma ) / ray**2 * nel
     +	                      /(t**1.5)/(nucleo(j)*niv_z(k,j)*niv_z(k,j))
					   gc_kj  = bco * theta(k)
				   else
					   gc_kj = 0.
				   end if

				   gr_kj  = q * phistar(k,j) 
     +					  * (1. + xistar(k,j)*C_ion(k,j))
     +					  * (1. + C_ion(k,j)/Ck_s)**alphstar(k,j)
     +				    + acstar(k,j) * gc_kj
     +					  * (khi(j)/(khi(j)+1.))**betstar(k,j)

				   dgr_kj = popi(k,j)*gr_kj* (
     +					xistar(k,j)/(1. + xistar(k,j)*C_ion(k,j))
     +					+ alphstar(k,j)/(C_ion(k,j) + Ck_s) )
				   dgc_kj = betstar(k,j)*gc_kj*
     +					(khi(j)**betstar(k,j))
     +					/((khi(j)+1.)**(betstar(k,j)+1))

				 else
			 	   gr_kj  = 0.
			 	   dgr_kj = 0.
				 end if
			   else
				if(igrad.eq.0) then
			 	   gr_kj = 0.
			 	else if(igrad.eq.1) then
			 	   gr_kj = g_grav
			 	end if
			 	dgr_kj = 0.
			   end if
			   g_rad(j)    = g_rad(j) + pond_ion(k,j)*gr_kj
			   dg_rad(j,j) = dg_rad(j,j) + pond_ion(k,j)* (
     +				     dgr_kj/N_chim(1) +
     +				     dgc_kj*N_chim(j) )
			end do
		end do
	end if

	end subroutine modul_g_rad2_bp
!=======================================================================
	subroutine gr_gazrare(nelem,jz,rarflag)
! identifie une configuration de gaz rare (G.Alecian, 26/06/2001)
	implicit none
	integer, dimension (6)      :: zgrar
	integer                     :: nelem,i,j,k,jz,rarflag
	data zgrar/2,10,18,36,54,86/
	do j=1,6
		if(zgrar(j).eq.(nelem-jz)) then
			rarflag=1
			exit
		end if
	end do
	end subroutine gr_gazrare
!=======================================================================
	subroutine gr_theta(t,theta,j)

! Version 1.0 (06/05/2004)
! Determination de la fonction Theta pour le g_cont de tous les ions de j,
! a une temperature donnee

!  Cette subroutine fait partie du module grcesam2.
!  Version 2.0 (7/01/2004)
!  Auteur:
!     Georges ALECIAN
!     LUTH - UMR 8102, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

	implicit none
	integer                     :: j
	real(DP)	                   :: t

! output du subroutine
	real(KIND=KIND(0D0)), dimension(0:pzi) :: theta

! variables locales
	integer                             :: k,m,ks
	real, parameter                     :: secu=100.
	real(KIND=KIND(0D0))                :: tke,uk,Qk
	real(KIND=KIND(0D0)), dimension(0:pzi)        :: fparti
	real(KIND=KIND(0D0)), dimension(30,0:pzi)     :: pniv

!=========== initialisation tableaux
	theta=0.
!===========
	
	call gr_bolm2(t,pniv,j)
	call gr_parti2(t,fparti,j)

	tke=t * 8.625e-5     ! temperature en eV

	do k=1,nint(zi(j)) ! on boucle sur les ions de j (sauf le premier)

		ks = niv_nb(k,j)
		if(ks.ge.1) then
			do m=1,ks   ! on boucle sur les niveaux de chaque ion k
				uk  = (el_pot(k-1,j)-niv_e(m,k-1,j))/tke
				if(uk.lt.1.e-2)  CYCLE
				if(uk.gt.20.)    CYCLE
				
				Qk=(uk**3)*(uk/(1.-exp(-uk)) - exp(uk)*log(1.-exp(-uk)))
				theta(k) = theta(k) + 
     +	              niv_q(m,k-1,j)*(niv_g(m,k-1,j)/niv_g(1,k-1,j)) * Qk
	
			end do
			theta(k)=pniv(1,k-1)*(fparti(k-1)/fparti(k))*theta(k)
		end if
	end do

	end subroutine gr_theta
!=======================================================================
	subroutine gr_bolm2(t,pniv,j)

! Version 1.0 (16/02/2004)
! Determination des populations des niveaux de tous les ions de j,
! a une temperature donnee, et en fraction: somme(pniv)=1 pour chq ion

!  Cette subroutine fait partie du module grcesam2.
!  Version 2.0 (7/01/2004)
!  Auteur:
!     Georges ALECIAN
!     LUTH - UMR 8102, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

	implicit none

	integer                                    :: j
	real(KIND=KIND(0D0))	                   :: t

! output du subroutine
	real(KIND=KIND(0D0)), dimension(30,0:pzi)   :: pniv

! variables locales
	integer                                     :: k,m,ks
	real(KIND=KIND(0D0)), dimension(30)         :: p,aa
	real(KIND=KIND(0D0))                        :: g1,gm,e,tke,so

	pniv=0.            ! pniv(niveau,degre_ionis) de l'element j

	tke=t*8.625e-5     ! temperature en eV

	do k=0,nint(zi(j)) ! on boucle sur les ions de j
		ks = niv_nb(k,j)
		if(ks.eq.1) then
			pniv(1,k) = 1.
		else
			so=1.
			p(1)=1.
			g1 = REAL(niv_g(1,k,j),KIND(p))
			do m=2,ks   ! on boucle sur les niveaux de chaque ion k
				gm   = REAL(niv_g(m,k,j),KIND(p))
				e    = niv_e(m,k,j)
				p(m) = (gm/g1)*exp(-e/tke)
				so   = so + p(m)
			end do
			aa(1)     = 1./so
			pniv(1,k) = aa(1)
			do m=2,ks
				aa(m)     = p(m)*aa(1)
				pniv(m,k) = aa(m)
			end do
		end if
	end do
      end subroutine gr_bolm2
!=======================================================================
	subroutine gr_parti2(t,fparti,j)

! Version 1.0 (06/05/2004)
! Determination des fonctions de partition de tous les ions de j,
! a une temperature donnee. Methode simple, suffisante pour le calcul des
! accelerations, a n'utiliser qu'avec les niveaux d'energies lus avec modul_niv.

!  Cette subroutine fait partie du module grcesam2.
!  Version 2.0 (7/01/2004)
!  Auteur:
!     Georges ALECIAN
!     LUTH - UMR 8102, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

	implicit none
	integer                                    :: j
	real(KIND=KIND(0D0))	                   :: t

! output du subroutine
	real(KIND=KIND(0D0)), dimension(0:pzi) :: fparti

! variables locales
	integer                             :: k,m,ks
	real, parameter                     :: secu=100.
	real(KIND=KIND(0D0))                :: tke

!=========== initialisation tableaux
	fparti=0.
!===========
	
	tke=t * 8.625e-5     ! temperature en eV

	do k=0,nint(zi(j)) ! on boucle sur les ions de j

		ks = niv_nb(k,j)
		
		do m=1,ks   ! on boucle sur les niveaux de chaque ion k

			if(niv_e(m,k,j)/tke.gt.secu) CYCLE   ! m suivant
			fparti(k) = fparti(k) + 
     +			      niv_g(m,k,j) * exp(-niv_e(m,k,j)/tke)

		end do

	end do
	
	end subroutine gr_parti2
!=======================================================================
	subroutine gr_popion(t,nel,fparti,popi,j)

! Version 1.0 (11/05/2004)
! Determination des populations de tous les ions de j,
! a une temperature donnee.

!  Cette subroutine fait temporairement partie du module grcesam2.
!  Version 2.0 (7/01/2004)
!  Auteur:
!     Georges ALECIAN
!     LUTH - UMR 8102, Observatoire de Meudon
!     F-92195 MEUDON CEDEX, FRANCE
!     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

	implicit none
	integer                                    :: j
	real(KIND=KIND(0D0))	                   :: t,nel

! input du subroutine
	real(KIND=KIND(0D0)), dimension(0:pzi) :: fparti
! output du subroutine
	real(KIND=KIND(0D0)), dimension(0:pzi,pnchim) :: popi

! variables locales
	integer                                :: k,ks,na0,na
	real, parameter                        :: secu=100.
	real(KIND=KIND(0D0))                   :: tke,r,d,an,c
	real(KIND=KIND(0D0)), dimension(0:pzi) :: ax,p

!=========== initialisations
	r   = 1.
	d   = 1.
	na0 = 0
	na  = nint(zi(j))
	an  = 1.          ! car on se contente des populations relatives
!===========
	
	tke = t * 8.625e-5     ! temperature en eV
	c = (1.415e+26 * 1.38e-16 * t)**1.5

	do k=0,nint(zi(j))-1  ! on boucle sur les ions de j

		ks = niv_nb(k,j)
		if(ks.ge.1) then		
			if(el_pot(k,j)/tke.gt.secu) EXIT
			ax(k) = (2./nel)*(fparti(k+1)/fparti(k))*c*exp(-el_pot(k,j)/tke)
			na=k
		else
			na0=k+1    ! premier ion de la base de donnees atomiques
		endif
	end do
	
	do k=na0,na
		ks = niv_nb(k,j)
		if(ks.ge.1) then		
			d=d*ax(k)
			r=r+d
		end if
	end do
	do k=0,na0
		popi(k,j)=0.
	end do
	p(na0)=1./r
	popi(na0,j)=p(na0)*an
	do k=na0+1,na+1
		p(k)=p(k-1)*ax(k-1)
		popi(k,j)=p(k)*an
	end do
	
	end subroutine gr_popion
!=======================================================================
	FUNCTION FT(T,N,X,F)
!
!	INTERPOLATION D'UNE LISTE A PAS NON CONSTANT.
!	Provenance: exablacop.for (extraction 15/9/98).
!     Auteur de FT inconnu (origine avant 1980 probablement).
!     Quelques modifications f77 et tests par G.Alecian
!	X est le tableau de la variable discrete 
!	F est le tableau de la fonction F(X) 
!	T est la valeur de X a laquelle on veut interpoler (resultat FT=F(T))
!	N est la dimension du tableau

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER :: N, I, J
	DIMENSION X(N),F(N)

	DO J=1,N
		I=J
		IF(T-X(J)) 3,2,1
2		FT=F(J)
		GOTO 1000
1		CONTINUE
	END DO
3	CONTINUE
	IF(I.EQ.1) I=2
	IF(I.GE.N) I=N-1
	T0=T-X(I-1)
	T1=T-X(I)
	T2=T-X(I+1)
	U0=X(I+1)-X(I-1)
	U1=X(I+1)-X(I)
	U2=X(I)-X(I-1)
	if((U0.eq.0.).or.(U1.eq.0.).or.(U2.eq.0.)) then
		print*,'I,T,U0,U1,U2: ', I, T,U0,U1,U2
		print*,'Probleme: modele hors limite => STOP'
		stop
	end if
	A=T0/U0
	B=T1/U1
	C=T2/U2
	D=T0/U1
	E=T1/U0
	FT=F(I+1)*A*B - F(I)*D*C + F(I-1)*E*C

1000	CONTINUE

	END FUNCTION FT
