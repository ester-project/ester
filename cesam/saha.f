
c*********************************************************	

	SUBROUTINE saha(xchim,t,ro,ioni,z_bar,nel,eta)	

c	routine public du module mod_etat
	
c	calcul des taux d'ionisation et charges moyennes
c	potentiels d'ionisation du Handbook of Chemistry and Physics
c	75-ieme edition 1995 p. 10-205 

c	on écrit la formule de saha sous la forme de Mihalas stellar
c	atmosphere, formule 5-17, on résout par Newton-Raphson 
c	ro/amu \sum Z_bar(ne) xchim - ne = 0

c	on utilise une adaptation de la correction de Clayton (2-254)
c	pour éviter la recombinaison  des ions à haute température
c	(voir notice de CESAM)

c	on regroupe les isotopes d'une même espèce
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	CESAM2k

c entrées:
c	xchim: composition chimique par mole
c	t: température
c	ro: densité

c sorties:
c	ioni(0:zi(i),i): taux d'ionisation pour chaque élément i
c	z_bar: charges moyennes 
c	nel: nombre d'électrons libres par volume
c	eta : dégénérescence

c-------------------------------------------------------------------------

	USE mod_donnees, ONLY : amu, echarg, eve, hpl, kbol, me, nchim,
	1 nom_elem, pi, zi
	USE mod_kind	
	USE mod_numerique, ONLY : bsp1dn, ferdir => fermi_dirac, no_croiss

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out), DIMENSION(0:,:) :: ioni	  	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: z_bar		  
	REAL (kind=dp), INTENT(out) :: nel, eta	  
	
	INTEGER, PARAMETER :: pf12=101,	!nb. de points de tabulation de F12
	1 mdg=4	!ordre des splines
	
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: ki, u, phi,
	1 gsg
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: dphi, ioni_esp,
	1 dioni_esp
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: zi_esp
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: dz_esp, z_esp, x_esp
	REAL (kind=dp), SAVE, DIMENSION(1,pf12) :: dgce
	REAL (kind=dp), SAVE, DIMENSION(pf12+mdg) :: f12t	
	REAL (kind=dp), SAVE, DIMENSION(pf12) :: f12
	REAL (kind=dp), DIMENSION(5) :: degen
	REAL (kind=dp), SAVE, DIMENSION(1) :: bid, dbid
	REAL (kind=dp), PARAMETER :: epsi=1.d-3, den_max=1.d15, p=4.d0
	REAL (kind=dp), SAVE :: cte1, cte2, phi_max, fin		
	REAL (kind=dp) :: a, cla, da, dden, dep, den, deta, df12l,
	1 dnel, f12l, kisrd, nelp, q, x, zeta
		
	INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: iso	
	INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: niso
	INTEGER, PARAMETER :: tourmax=30 
	INTEGER, SAVE :: knot, nesp, nzi	
	INTEGER i, j, k, tour, l
	
	LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: no_iden	
	LOGICAL, SAVE :: init=.TRUE., no_conv=.TRUE.

c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

c	cubique pour lissage de la correction des poids statistiques
c	entre x=0 et x=q
c	zeta(0,q)=0, zeta'(0,q)=0, zeta(q,q)=1, zeta'(q,q)=0,
c	zeta(x<0,q)=0, zeta(x>q,q)=1
	
	zeta(x,q)=(x/q)**2*(-2.d0*x/q+3.d0)

	IF(init)THEN
	 init=.FALSE. ; nzi=NINT(maxval(zi))
	 ALLOCATE(no_iden(nchim),niso(nchim),iso(nchim,nchim),zi_esp(nzi))
	
	 no_iden=.TRUE. ; niso=0
	 
c	 constantes diverses	   		 
c	 epsi : précision sur nel
c	 tourmax : limitation du nombre d'itérations sur tour
c	 den_max : max pour le dénominateur
c	 p : correction de recombinaison des ions

	 cte1=SQRT(hpl/2.d0/me/kbol*hpl)**3/4.d0/pi	!pour calcul F1/2
	 cte2=echarg**3*SQRT(4.d0*pi/kbol/amu)/kbol	 
	 phi_max=LOG(den_max)	!max pour ki(j,i)/t+ksi
	 	 
c	 identification des espèces
c	 on recherche les d'espèces chimiques distinctes à partir
c	 des charges, pour les espèces reconnues on construit les
c	 tableaux des potentiels
c	 d'ionisation et des poids statistiques

c	 pour ajouter une espèce chimique il suffit d'insérer les tables
c	 des potentiels d'ionisation et des poids statistiques de cette
c	 espèce chimique

	 nesp=0		!nombre d'espèces distinctes	 
	 DO i=1,nchim	
	  IF(no_iden(i))THEN	!l'espèce n'est pas encore identifiée
	   nesp=nesp+1		!nb. d'espèces d'isotopes  
	   no_iden(i)=.FALSE.	!l'espèce est désormais identifiée
	   zi_esp(nesp)=zi(i)	!charge de l'espèce
	   niso(nesp)=1		!nb. isotope de l'espèce nesp	   	   
	   iso(nesp,1)=i	!iden. isotope i de l'espèce nesp
	   DO j=i+1,nchim	!recherche des isotopes	de i
	    IF(no_iden(j) .AND. NINT(zi(i)) == NINT(zi(j)))THEN !j est un isotope de i
	     no_iden(j)=.FALSE.		!j est identifie d'espèce nesp=i
	     niso(nesp)=niso(nesp)+1	!nb. isotope de l'espèce nesp
	     iso(nesp,niso(nesp))=j	!iden. isotope j de l'espèce nesp
	    ENDIF
	   ENDDO	!j
	  ENDIF
	 ENDDO		!i
c	 PRINT*,'nb espèces, nb isotopes par espèce'
c	 PRINT*,nesp,(niso(k),k=1,nesp)
c	 PRINT*,'identification des isotopes'
c	 DO k=1,nesp
c	  PRINT*,(iso(k,l),l=1,niso(k))
c	 ENDDO
c	 WRITE(*,2000)(zi_esp(k),k=1,nesp)
c	 PAUSE'iso'
	
	 DEALLOCATE(no_iden)
	 ALLOCATE(ki(nzi,nesp),u(0:nzi,nesp),phi(nzi,nesp),gsg(nzi,nesp))
	
c	 formation des tableaux des potentiels d'ionisation, poids statistiques

	 DO i=1,nesp
	  IF(NINT(zi_esp(i)) == 1)THEN 			!H
	   ki(1,i)=13.595d0	!potentiel d'ionisation 1
	   u(0,i)=2.d0	!poids statistique du fondamental
	   u(1,i)=1.d0	!poids statistique du niveau 1
	   
	  ELSEIF(NINT(zi_esp(i)) == 2)THEN		!He
	   ki(1,i)=24.580d0 !potentiel d'ionisation 1
	   ki(2,i)=54.403d0 !potentiel d'ionisation 2	   
	   u(0,i)=1.d0	!poids statistique du fondamental
	   u(1,i)=2.d0	!poids statistique du niveau 1
	   u(2,i)=1.d0	!poids statistique du niveau 2
	   
	  ELSEIF(NINT(zi_esp(i)) == 3)THEN		!Li
	   ki(1,i)=5.39172d0 !potentiel d'ionisation 1
	   ki(2,i)=75.64018d0 !potentiel d'ionisation 2
	   ki(3,i)=122.45429d0 !potentiel d'ionisation 3	   
	   u(0,i)=2.d0	!poids statistique du fondamental
	   u(1,i)=1.d0	!poids statistique du niveau 1
	   u(2,i)=2.d0	!poids statistique du niveau 2
	   u(3,i)=2.d0	!poids statistique du niveau 3
	   
	  ELSEIF(NINT(zi_esp(i)) == 4)THEN		!Be
	   ki(1,i)=9.32263d0 !potentiel d'ionisation 1
	   ki(2,i)=18.21116d0 !potentiel d'ionisation 2
	   ki(3,i)=153.89661d0 !potentiel d'ionisation 3
	   ki(4,i)=217.71865d0 !potentiel d'ionisation 4
	   u(0,i)=1.d0	!poids statistique du fondamental
	   u(1,i)=2.d0	!poids statistique du niveau 1
	   u(2,i)=1.d0	!poids statistique du niveau 2
	   u(3,i)=2.d0	!poids statistique du niveau 3
	   u(4,i)=2.d0	!poids statistique du niveau 4

	  ELSEIF(NINT(zi_esp(i)) == 5)THEN		!B	
	   ki(1,i)=8.29803d0 !potentiel d'ionisation 1
	   ki(2,i)=25.15484d0 !potentiel d'ionisation 2
	   ki(3,i)=37.93064d0 !potentiel d'ionisation 3
	   ki(4,i)=259.37521d0 !potentiel d'ionisation 4
	   ki(5,i)=340.22580d0 !potentiel d'ionisation 5
	   u(0,i)=6.d0	!poids statistique du fondamental
	   u(1,i)=1.d0	!poids statistique du niveau 1 
	   u(2,i)=2.d0	!poids statistique du niveau 2
	   u(3,i)=1.d0	!poids statistique du niveau 3
	   u(4,i)=2.d0	!poids statistique du niveau 4
	   u(5,i)=2.d0	!poids statistique du niveau 5

	  ELSEIF(NINT(zi_esp(i)) == 6)THEN		!C	
	   ki(1,i)=11.26030d0 !potentiel d'ionisation 1
	   ki(2,i)=24.38332d0 !potentiel d'ionisation 2
	   ki(3,i)=47.8878d0 !potentiel d'ionisation 3
	   ki(4,i)=64.4939d0 !potentiel d'ionisation 4
	   ki(5,i)=391.087d0 !potentiel d'ionisation 5
	   ki(6,i)=489.99334d0 !potentiel d'ionisation 6
	   u(0,i)=9.d0	!poids statistique du fondamental
	   u(1,i)=6.d0	!poids statistique du niveau 1 
	   u(2,i)=1.d0	!poids statistique du niveau 2
	   u(3,i)=2.d0	!poids statistique du niveau 3
	   u(4,i)=1.d0	!poids statistique du niveau 4
	   u(5,i)=2.d0	!poids statistique du niveau 5
	   u(6,i)=2.d0	!poids statistique du niveau 6

	  ELSEIF(NINT(zi_esp(i)) == 7)THEN		!N
	   ki(1,i)=14.53414d0 !potentiel d'ionisation 1
	   ki(2,i)=29.296013d0 !potentiel d'ionisation 2
	   ki(3,i)=47.44924d0 !potentiel d'ionisation 3
	   ki(4,i)=77.4735d0 !potentiel d'ionisation 4
	   ki(5,i)=97.8902d0 !potentiel d'ionisation 5
	   ki(6,i)=552.0718d0 !potentiel d'ionisation 6
	   ki(7,i)=667.046d0 !potentiel d'ionisation 7
	   u(0,i)=4.d0	!poids statistique du fondamental
	   u(1,i)=9.d0	!poids statistique du niveau 1
	   u(2,i)=6.d0	!poids statistique du niveau 2
	   u(3,i)=1.d0	!poids statistique du niveau 3
	   u(4,i)=2.d0	!poids statistique du niveau 4
	   u(5,i)=1.d0	!poids statistique du niveau 5
	   u(6,i)=2.d0	!poids statistique du niveau 6
	   u(7,i)=2.d0	!poids statistique du niveau  7  	   
	   	   	   	   
	  ELSEIF(NINT(zi_esp(i)) == 8)THEN		!O
	   ki(1,i)=13.61806d0 !potentiel d'ionisation 1
	   ki(2,i)=35.11730d0 !potentiel d'ionisation 2
	   ki(3,i)=54.9355d0 !potentiel d'ionisation 3
	   ki(4,i)=77.41353d0 !potentiel d'ionisation 4
	   ki(5,i)=113.8990d0 !potentiel d'ionisation 5
	   ki(6,i)=138.1197d0 !potentiel d'ionisation 6
	   ki(7,i)=739.29d0 !potentiel d'ionisation 7
	   ki(8,i)=871.4101d0 !potentiel d'ionisation 8
	   u(0,i)=9.d0	!poids statistique du fondamental
	   u(1,i)=4.d0	!poids statistique du niveau 1
	   u(2,i)=9.d0	!poids statistique du niveau 2
	   u(3,i)=6.d0	!poids statistique du niveau 3  
	   u(4,i)=1.d0	!poids statistique du niveau 4	   
	   u(5,i)=2.d0	!poids statistique du niveau 5	   
	   u(6,i)=1.d0	!poids statistique du niveau 6
	   u(7,i)=2.d0	!poids statistique du niveau 7	   
	   u(8,i)=2.d0	!poids statistique du niveau 8

	  ELSEIF(NINT(zi_esp(i)) == 10)THEN	!Ne
	   ki(1,i)=21.56454d0 !potentiel d'ionisation 1
	   ki(2,i)=40.96328d0 !potentiel d'ionisation 2
	   ki(3,i)=63.45d0 !potentiel d'ionisation 3
	   ki(4,i)=97.12d0 !potentiel d'ionisation 4
	   ki(5,i)=126.21d0 !potentiel d'ionisation 5
	   ki(6,i)=157.93d0 !potentiel d'ionisation 6
	   ki(7,i)=207.2759d0 !potentiel d'ionisation 7
	   ki(8,i)=239.0989d0 !potentiel d'ionisation 8
	   ki(9,i)=1195.8286d0 !potentiel d'ionisation 9
	   ki(10,i)=1362.1995d0 !potentiel d'ionisation 10
	   u(0,i)=1.d0	!poids statistique du fondamental 10 électrons
	   u(1,i)=6.d0	!poids statistique du niveau 1 9 électrons
	   u(2,i)=9.d0	!poids statistique du niveau 2 8 électrons
	   u(3,i)=4.d0	!poids statistique du niveau 3 7 électrons
	   u(4,i)=9.d0	!poids statistique du niveau 4 6 électrons
	   u(5,i)=6.d0	!poids statistique du niveau 5 5 électrons
	   u(6,i)=1.d0	!poids statistique du niveau 6 4 électrons
	   u(7,i)=2.d0	!poids statistique du niveau 7 3 électrons
	   u(8,i)=1.d0	!poids statistique du niveau 8 2 électrons
	   u(9,i)=3.d0	!poids statistique du niveau 9 1 électrons
	   u(10,i)=3.d0	!poids statistique du niveau 10 0 électrons
	   
	  ELSEIF(NINT(zi_esp(i)) == 11)THEN		!Na
	   ki(1,i)=5.13908d0	!potentiel d'ionisation 1
	   ki(2,i)=47.2864d0	!potentiel d'ionisation 2
	   ki(3,i)=71.62d0	!potentiel d'ionisation 3 
	   ki(4,i)=98.91d0	!potentiel d'ionisation 4 
	   ki(5,i)=138.40d0	!potentiel d'ionisation 5
	   ki(6,i)=172.18d0	!potentiel d'ionisation 6
	   ki(7,i)=208.50d0	!potentiel d'ionisation 7
	   ki(8,i)=264.25d0	!potentiel d'ionisation 8
	   ki(9,i)=299.864d0	!potentiel d'ionisation 9
	   ki(10,i)=1465.121d0	!potentiel d'ionisation 10
	   ki(11,i)=1648.702d0	!potentiel d'ionisation 11
	   u(0,i)=2.d0	!poids statistique du fondamental
	   u(1,i)=1.d0	!poids statistique du niveau 1
	   u(2,i)=6.d0	!poids statistique du niveau 2
	   u(3,i)=9.d0	!poids statistique du niveau 3
	   u(4,i)=4.d0	!poids statistique du niveau 4
	   u(5,i)=9.d0	!poids statistique du niveau 5
	   u(6,i)=6.d0	!poids statistique du niveau 6
	   u(7,i)=1.d0	!poids statistique du niveau 7
	   u(8,i)=2.d0	!poids statistique du niveau 8
	   u(9,i)=1.d0	!poids statistique du niveau 9
	   u(10,i)=2.d0	!poids statistique du niveau 10
	   u(11,i)=2.d0	!poids statistique du niveau 11

	  ELSEIF(NINT(zi_esp(i)) == 12)THEN		!Mg
	   ki(1,i)=7.64624d0	!potentiel d'ionisation 1
	   ki(2,i)=15.03528d0	!potentiel d'ionisation 2
	   ki(3,i)=80.1437d0	!potentiel d'ionisation 3 
	   ki(4,i)=109.2655d0	!potentiel d'ionisation 4 
	   ki(5,i)=141.27d0	!potentiel d'ionisation 5
	   ki(6,i)=186.75d0	!potentiel d'ionisation 6
	   ki(7,i)=225.02d0	!potentiel d'ionisation 7
	   ki(8,i)=265.96d0	!potentiel d'ionisation 8
	   ki(9,i)=328.06d0	!potentiel d'ionisation 9
	   ki(10,i)=367.50d0	!potentiel d'ionisation 10
	   ki(11,i)=1761.805d0	!potentiel d'ionisation 11
	   ki(12,i)=1962.6650d0	!potentiel d'ionisation 12
	   u(0,i)=1.d0	!poids statistique du fondamental
	   u(1,i)=2.d0	!poids statistique du niveau 1
	   u(2,i)=1.d0	!poids statistique du niveau 2
	   u(3,i)=6.d0	!poids statistique du niveau 3
	   u(4,i)=9.d0	!poids statistique du niveau 4
	   u(5,i)=4.d0	!poids statistique du niveau 5
	   u(6,i)=9.d0	!poids statistique du niveau 6
	   u(7,i)=6.d0	!poids statistique du niveau 7
	   u(8,i)=1.d0	!poids statistique du niveau 8
	   u(9,i)=2.d0	!poids statistique du niveau 9
	   u(10,i)=1.d0	!poids statistique du niveau 10
	   u(11,i)=2.d0	!poids statistique du niveau 11
	   u(12,i)=2.d0	!poids statistique du niveau 12
   
	  ELSEIF(NINT(zi_esp(i)) == 13)THEN		!Al
	   ki(1,i)=5.98577d0	!potentiel d'ionisation 1
	   ki(2,i)=18.82856d0	!potentiel d'ionisation 2
	   ki(3,i)=28.44765d0	!potentiel d'ionisation 3 
	   ki(4,i)=119.992d0	!potentiel d'ionisation 4 
	   ki(5,i)=153.825d0	!potentiel d'ionisation 5
	   ki(6,i)=190.49d0	!potentiel d'ionisation 6
	   ki(7,i)=241.76d0	!potentiel d'ionisation 7
	   ki(8,i)=284.66d0	!potentiel d'ionisation 8
	   ki(9,i)=330.13d0	!potentiel d'ionisation 9
	   ki(10,i)=398.75d0	!potentiel d'ionisation 10
	   ki(11,i)=442.00d0	!potentiel d'ionisation 11
	   ki(12,i)=2085.98d0	!potentiel d'ionisation 12
	   ki(13,i)=2304.1410d0	!potentiel d'ionisation 13	   
	   u(0,i)=6.d0	!poids statistique du fondamental
	   u(1,i)=1.d0	!poids statistique du niveau 1
	   u(2,i)=2.d0	!poids statistique du niveau 2
	   u(3,i)=1.d0	!poids statistique du niveau 3
	   u(4,i)=6.d0	!poids statistique du niveau 4
	   u(5,i)=9.d0	!poids statistique du niveau 5
	   u(6,i)=4.d0	!poids statistique du niveau 6
	   u(7,i)=9.d0	!poids statistique du niveau 7
	   u(8,i)=6.d0	!poids statistique du niveau 8
	   u(9,i)=1.d0	!poids statistique du niveau 9
	   u(10,i)=2.d0	!poids statistique du niveau 10
	   u(11,i)=1.d0	!poids statistique du niveau 11
	   u(12,i)=2.d0	!poids statistique du niveau 12
	   u(13,i)=2.d0	!poids statistique du niveau 13
	   
	  ELSEIF(NINT(zi_esp(i)) == 14)THEN		!Si
	   ki(1,i)=8.15169d0	!potentiel d'ionisation 1
	   ki(2,i)=16.34585d0	!potentiel d'ionisation 2
	   ki(3,i)=33.49302d0	!potentiel d'ionisation 3 
	   ki(4,i)=45.14181d0	!potentiel d'ionisation 4 
	   ki(5,i)=166.767d0	!potentiel d'ionisation 5
	   ki(6,i)=205.27d0	!potentiel d'ionisation 6
	   ki(7,i)=246.5d0	!potentiel d'ionisation 7
	   ki(8,i)=303.54d0	!potentiel d'ionisation 8
	   ki(9,i)=351.12d0	!potentiel d'ionisation 9
	   ki(10,i)=401.37d0	!potentiel d'ionisation 10
	   ki(11,i)=476.36d0	!potentiel d'ionisation 11
	   ki(12,i)=523.42d0	!potentiel d'ionisation 12
	   ki(13,i)=2437.63d0	!potentiel d'ionisation 13
	   ki(14,i)=2673.182d0	!potentiel d'ionisation 14	   	   
	   u(0,i)=9.d0	!poids statistique du fondamental
	   u(1,i)=6.d0	!poids statistique du niveau 1
	   u(2,i)=1.d0	!poids statistique du niveau 2
	   u(3,i)=2.d0	!poids statistique du niveau 3
	   u(4,i)=1.d0	!poids statistique du niveau 4
	   u(5,i)=6.d0	!poids statistique du niveau 5
	   u(6,i)=9.d0	!poids statistique du niveau 6
	   u(7,i)=4.d0	!poids statistique du niveau 7
	   u(8,i)=9.d0	!poids statistique du niveau 8
	   u(9,i)=6.d0	!poids statistique du niveau 9
	   u(10,i)=1.d0	!poids statistique du niveau 10
	   u(11,i)=2.d0	!poids statistique du niveau 11
	   u(12,i)=1.d0	!poids statistique du niveau 12
	   u(13,i)=2.d0	!poids statistique du niveau 13
	   u(14,i)=2.d0	!poids statistique du niveau 14	
	      
	  ELSEIF(NINT(zi_esp(i)) == 17)THEN		!Cl
	   ki(1,i)=12.96764	 !potentiel d'ionisation 1
	   ki(2,i)=23.814	 !potentiel d'ionisation 2
	   ki(3,i)=39.61	 !potentiel d'ionisation 3
	   ki(4,i)=53.4652	 !potentiel d'ionisation 4
	   ki(5,i)=67.8		 !potentiel d'ionisation 5
	   ki(6,i)=97.03	 !potentiel d'ionisation 6
	   ki(7,i)=114.1958	 !potentiel d'ionisation 7
	   ki(8,i)=348.28	 !potentiel d'ionisation 8
	   ki(9,i)=400.06	 !potentiel d'ionisation 9
	   ki(10,i)=455.63	 !potentiel d'ionisation 10
	   ki(11,i)=529.28	 !potentiel d'ionisation 11
	   ki(12,i)=591.99	 !potentiel d'ionisation 12
	   ki(13,i)=656.71	 !potentiel d'ionisation 13
	   ki(14,i)=749.76	 !potentiel d'ionisation 14
	   ki(15,i)=809.40	 !potentiel d'ionisation 15
	   ki(16,i)=3658.521	 !potentiel d'ionisation 16
	   ki(17,i)=3946.2960	 !potentiel d'ionisation 17
	   u(0,i)=6	!poids statistique du fondamental
	   u(1,i)=9	!poids statistique du niveau 1
	   u(2,i)=4	!poids statistique du niveau 2
	   u(3,i)=9	!poids statistique du niveau 3
	   u(4,i)=6	!poids statistique du niveau 4
	   u(5,i)=1	!poids statistique du niveau 5
	   u(6,i)=2	!poids statistique du niveau 6
	   u(7,i)=1	!poids statistique du niveau 7
	   u(8,i)=6	!poids statistique du niveau 8
	   u(9,i)=9	!poids statistique du niveau 9
	   u(10,i)=4	!poids statistique du niveau 10
	   u(11,i)=9	!poids statistique du niveau 11
	   u(12,i)=6	!poids statistique du niveau 12
	   u(13,i)=1	!poids statistique du niveau 13
	   u(14,i)=2	!poids statistique du niveau 14
	   u(15,i)=1	!poids statistique du niveau 15
	   u(16,i)=2	!poids statistique du niveau 16
	   u(17,i)=2	!poids statistique du niveau 17

	  ELSEIF(NINT(zi_esp(i)) == 20)THEN	!Ca
	   ki(1,i)=6.11316d0 !potentiel d'ionisation 1
	   ki(2,i)=11.87172d0 !potentiel d'ionisation 2
	   ki(3,i)=50.9131d0 !potentiel d'ionisation 3
	   ki(4,i)=67.27d0 !potentiel d'ionisation 4
	   ki(5,i)=84.50d0 !potentiel d'ionisation 5
	   ki(6,i)=108.78d0 !potentiel d'ionisation 6
	   ki(7,i)=127.2d0 !potentiel d'ionisation 7
	   ki(8,i)=147.24d0 !potentiel d'ionisation 8
	   ki(9,i)=188.54d0 !potentiel d'ionisation 9
	   ki(10,i)=211.275d0 !potentiel d'ionisation 10
	   ki(11,i)=591.9d0 !potentiel d'ionisation 11
	   ki(12,i)=657.2d0 !potentiel d'ionisation 12
	   ki(13,i)=726.6d0 !potentiel d'ionisation 13
	   ki(14,i)=817.6d0 !potentiel d'ionisation 14
	   ki(15,i)=894.5d0 !potentiel d'ionisation 15
	   ki(16,i)=974.d0 !potentiel d'ionisation 16
	   ki(17,i)=1087.d0 !potentiel d'ionisation 17
	   ki(18,i)=1157.8d0 !potentiel d'ionisation 18
	   ki(19,i)=5128.8d0 !potentiel d'ionisation 19
	   ki(20,i)=5469.864d0!potentiel d'ionisation 20	    
	   u(0,i)=1.d0	!poids statistique du fondamental 20 électrons
	   u(1,i)=2.d0	!poids statistique du niveau 1 19 électrons
	   u(2,i)=1.d0	!poids statistique du niveau 2 18 électrons
	   u(3,i)=6.d0	!poids statistique du niveau 3 17 électrons
	   u(4,i)=9.d0	!poids statistique du niveau 4 16 électrons
	   u(5,i)=4.d0	!poids statistique du niveau 5 15 électrons
	   u(6,i)=9.d0	!poids statistique du niveau 6 14 électrons
	   u(7,i)=6.d0	!poids statistique du niveau 7 13 électrons
	   u(8,i)=1.d0	!poids statistique du niveau 8 12 électrons
	   u(9,i)=2.d0	!poids statistique du niveau 9 11 électrons
	   u(10,i)=1.d0	!poids statistique du niveau 10 10 électrons
	   u(11,i)=6.d0	!poids statistique du niveau 11 9 électrons
	   u(12,i)=9.d0	!poids statistique du niveau 12 8 électrons
	   u(13,i)=4.d0	!poids statistique du niveau 13 7 électrons
	   u(14,i)=9.d0	!poids statistique du niveau 14 6 électrons
	   u(15,i)=6.d0	!poids statistique du niveau 15 5 électrons
	   u(16,i)=1.d0	!poids statistique du niveau 16 4 électrons
	   u(17,i)=2.d0	!poids statistique du niveau 17 3 électrons
	   u(18,i)=1.d0	!poids statistique du niveau 18 2 électrons
	   u(19,i)=2.d0	!poids statistique du niveau 19 1 électrons
	   u(20,i)=2.d0	!poids statistique du niveau 20 0 électrons
	    
	  ELSEIF(NINT(zi_esp(i)) == 26)THEN		!Fe
	   ki(1,i)=7.9024d0	!potentiel d'ionisation 1
	   ki(2,i)=16.1878d0	!potentiel d'ionisation 2
	   ki(3,i)=30.652d0	!potentiel d'ionisation 3
	   ki(4,i)=54.8d0	!potentiel d'ionisation 4
	   ki(5,i)=75.0d0	!potentiel d'ionisation 5
	   ki(6,i)=99.1d0	!potentiel d'ionisation 6
	   ki(7,i)=124.98d0	!potentiel d'ionisation 7
	   ki(8,i)=151.06d0	!potentiel d'ionisation 8
	   ki(9,i)=233.6d0	!potentiel d'ionisation 9
	   ki(10,i)=262.1d0	!potentiel d'ionisation 10
	   ki(11,i)=290.2d0	!potentiel d'ionisation 11
	   ki(12,i)=330.8d0	!potentiel d'ionisation 12
	   ki(13,i)=361.0d0	!potentiel d'ionisation 13
	   ki(14,i)=392.2d0	!potentiel d'ionisation 14
	   ki(15,i)=457.d0	!potentiel d'ionisation 15
	   ki(16,i)=489.256d0	!potentiel d'ionisation 16
	   ki(17,i)=1266.d0	!potentiel d'ionisation 17
	   ki(18,i)=1358.d0	!potentiel d'ionisation 18 
	   ki(19,i)=1456.d0	!potentiel d'ionisation 19
	   ki(20,i)=1582.d0	!potentiel d'ionisation 20
	   ki(21,i)=1689.d0	!potentiel d'ionisation 21
	   ki(22,i)=1799.d0	!potentiel d'ionisation 22
	   ki(23,i)=1950.d0	!potentiel d'ionisation 23
	   ki(24,i)=2023.d0	!potentiel d'ionisation 24
	   ki(25,i)=8828.d0	!potentiel d'ionisation 25
	   ki(26,i)=9277.69d0	!potentiel d'ionisation 26	    
	   u(0,i)=25.d0	!poids statistique du fondamental 26 électrons
	   u(1,i)=30.d0	!poids statistique du niveau 1 25 électrons
	   u(2,i)=25.d0	!poids statistique du niveau 2 24 électrons
	   u(3,i)=6.d0	!poids statistique du niveau 3 23 électrons
	   u(4,i)=25.d0	!poids statistique du niveau 4 22 électrons
	   u(5,i)=12.d0	!poids statistique du niveau 5 21 électrons
	   u(6,i)=21.d0	!poids statistique du niveau 6 20 électrons
	   u(7,i)=10.d0	!poids statistique du niveau 7 19 électrons
	   u(8,i)=1.d0	!poids statistique du niveau 8 18 électrons
	   u(9,i)=6.d0	!poids statistique du niveau 9  17 électrons
	   u(10,i)=9.d0	!poids statistique du niveau 10 16 électrons
	   u(11,i)=4.d0	!poids statistique du niveau 11 15 électrons
	   u(12,i)=9.d0	!poids statistique du niveau 12 14 électrons
	   u(13,i)=6.d0	!poids statistique du niveau 13 13 électrons
	   u(14,i)=1.d0	!poids statistique du niveau 14 12 électrons
	   u(15,i)=2.d0	!poids statistique du niveau 15 11 électrons
	   u(16,i)=1.d0	!poids statistique du niveau 16 10 électrons 
	   u(17,i)=6.d0	!poids statistique du niveau 17 9 électrons
	   u(18,i)=9.d0	!poids statistique du niveau 18 8 électrons 
	   u(19,i)=4.d0	!poids statistique du niveau 19 7 électrons
	   u(20,i)=9.d0	!poids statistique du niveau 20 6 électrons
	   u(21,i)=6.d0	!poids statistique du niveau 21 5 électrons
	   u(22,i)=1.d0	!poids statistique du niveau 22 4 électrons
	   u(23,i)=2.d0	!poids statistique du niveau 23 3 électrons
	   u(24,i)=1.d0	!poids statistique du niveau 24 2 électrons
	   u(25,i)=2.d0	!poids statistique du niveau 25 1 electron
	   u(26,i)=1.d0	!poids statistique du niveau 26 0 electron
	   	    
c	  ELSEIF(NINT(zi_esp(i)) == )THEN		!isotopes de
c	   ki(,i)= !potentiel d'ionisation 
c	   ki(,i)= !potentiel d'ionisation 
c	   u(0,i)=	!poids statistique du fondamental
c	   u(,i)=	!poids statistique du niveau 
c	   u(,i)=	!poids statistique du niveau 

	  ELSE
	   PRINT*,'dans saha, erreur élément non prévu nom_elem=',
	1 (nom_elem(iso(i,j)),j=1,niso(i))
	   PRINT*,'charge:',NINT(zi_esp(i))	; PRINT*,'arret' ; STOP
	  ENDIF
	 ENDDO
	
c	 les ki/k en cgs

	 DO i=1,nesp
	  DO j=1,NINT(zi_esp(i))
	   ki(j,i)=ki(j,i)*eve/kbol
	   gsg(j,i)=u(j-1,i)/u(j,i)	!rapport des poids statistiques
	  ENDDO
	 ENDDO
 
	 DEALLOCATE(u)
	 
c	 tabulation des fonctions F12 de Fermi-Dirac de dep à fin
c	 signe opposé de celui de Clayton 

	 dep=-30.d0 ; fin=25.d0 ; a=(fin-dep)/DBLE(pf12-1)
	 DO i=1,pf12
	  dgce(1,i)=dep+a*(i-1) ; CALL ferdir(dgce(1,i),degen)
	  f12(i)=degen(2)
	 ENDDO
	 	  	 
	 CALL bsp1dn(1,dgce,f12,f12t,pf12,mdg,knot,.FALSE.,f12(2),l,
	1 bid,dbid)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 1 dans saha' ; STOP
         ENDIF
	
c	 eta=-1.3d0	!test de la tabulation ; eta=3.5d0 ; eta=15.d0 
c	 CALL ferdir(eta,degen) ; WRITE(*,2000)eta,degen(2) ; f12l=degen(2)
c	 CALL bsp1dn(1,dgce,f12,f12t,pf12,mdg,knot,.TRUE.,f12l,l,bid,dbid)
c	 WRITE(*,2000)f12l,bid(1) ; PAUSE'f12l'	
	 	
	ENDIF

	ALLOCATE(x_esp(nesp),z_esp(nesp),dz_esp(nesp),
	1 ioni_esp(0:nzi,nesp),dioni_esp(0:nzi,nesp),dphi(nzi,nesp))
		
c	determination des abondances/mole de chaque espèce

	DO k=1,nesp	!pour chaque espèce
	 x_esp(k)=0.d0
	 DO l=1,niso(k)	!pour chaque isotope d'espèce k
	  x_esp(k)=x_esp(k)+xchim(iso(k,l))
	 ENDDO		!iso(k,l): indice de l'isotope l d'espèce k
	ENDDO
c	PRINT*,'abondances' ; WRITE(*,2000)xchim
c	PRINT*,'abondance par espèce' ; WRITE(*,2000)x_esp ; PAUSE'x_esp'
	
c	initialisation du nombre d'électrons libres

	IF(t < 5.d3)THEN	!x_esp(1) est H, x_esp(2) est He
	 nel=x_esp(1)*1.d-10	!H est neutre
	 z_esp=zi_esp(1:nesp)*1.d-5	!init. à 1.d-5 ionisation tot.
	ELSEIF(t < 1.d4)THEN
	 nel=x_esp(1)/1000.d0	!H est 1/1000 ionisé
	 z_esp=zi_esp(1:nesp)*1.d-1	!init. à 10% ionisation tot.
	ELSEIF(t < 5.d4)THEN
	 nel=x_esp(1)/10.d0	!H est 1/10 ionisé
	 z_esp=zi_esp(1:nesp)*3.d-1	!init. à 30% ionisation totale
	ELSEIF(t < 1.d5)THEN
	 nel=x_esp(1)	!H est ionisé
	 IF(nesp > 1)nel=nel+x_esp(2) !He est ionisé 1 fois
	 z_esp=zi_esp(1:nesp)*6.d-1	!init. à 60% ionisation totale
	ELSE
	 nel=x_esp(1)!H est ionisé
	 IF(nesp > 1)nel=nel+2.d0*x_esp(2)	!He est ionisé 2 fois
	 z_esp=zi_esp(1:nesp)	!initialisation à 100% ionisation totale
	ENDIF
	nel=nel*ro/amu
c	PRINT*,'nel initial' ; WRITE(*,2000)nel ; PRINT*,'convergence'
c	PRINT*,'t,ro' ; WRITE(*,2000)t,ro ; PAUSE't,ro'

	tour=0 ; nelp=nel ; dnel=1.d3	!pour le WHILE
	
	DO WHILE(dnel > epsi .AND. tour < tourmax)
	 tour=tour+1
	 
c	 dégénérescence

	 df12l=cte1/SQRT(t)**3 ; f12l=nel*df12l
	 IF(f12l < f12(1))THEN
	  eta=dgce(1,1) ; deta=0.d0
	 ELSEIF(f12l > f12(pf12))THEN
	  eta=dgce(1,pf12) ; deta=0.d0	  
	 ELSE	 
	  CALL bsp1dn(1,dgce,f12,f12t,pf12,mdg,knot,.TRUE.,f12l,l,bid,dbid)
	  eta=bid(1) ; deta=dbid(1)*df12l
	 ENDIF
	 
c	 calcul des ki/Rd (Clayton 2-237, 2-254)
	 kisrd=sum(z_esp*(z_esp+1.d0)*x_esp) ; kisrd=cte2*SQRT(ro*kisrd/t)
c	 PRINT*,kisrd

c	 calcul des phi

	 DO i=1,nesp
c	  PRINT*,nom_elem(i)
	  DO j=1,NINT(zi_esp(i))	    
	   cla=ki(j,i)/kisrd/DBLE(j)-1.d0	!corr. de Clayon Xi/Xid-1
	   IF(cla >= p)THEN			!cubique entre 0 et p
	    cla=1.d0
	   ELSEIF(cla <= 0.d0)THEN
	    cla=0.d0
	   ELSE
	    cla=zeta(cla,p)			!la cubique
	   ENDIF

	   a=ki(j,i)/t+eta ; da=deta
	   
	   IF(a < phi_max)THEN
	    phi(j,i)=gsg(j,i)*exp(a)*cla ; dphi(j,i)=phi(j,i)*da 
	   ELSE
	    phi(j,i)=den_max ; dphi(j,i)=0.d0
	   ENDIF
	  ENDDO
	 ENDDO

c	 le dénominateur ne doit pas depasser den_max (1.d15)

	 ioni_esp=0.d0 ; dioni_esp=0.d0	
	 DO i=1,nesp
	  den=0.d0 ; dden=0.d0
	  j=0
	  DO WHILE(den < den_max .AND. j < NINT(zi_esp(i)))
	   j=j+1
	   dden=dden*phi(j,i)+(den+1.d0)*dphi(j,i); den=(den+1.d0)*phi(j,i)
	  ENDDO
	  den=den+1.d0

c	  le dernier niveau j atteint assure un den < den_max
c	  on met tous les ioni des autres niveaux a 0
	  
	  ioni_esp(j,i)=1.d0/den		!le niveau j
	  dioni_esp(j,i)=-ioni_esp(j,i)/den*dden
	  DO k=j-1,0,-1
	   ioni_esp(k,i)=ioni_esp(k+1,i)*phi(k+1,i)
	  ENDDO	  
	 ENDDO
	 
c	 calcul des charges moyennes

	 z_esp=0.d0 ; dz_esp=0.d0
	 DO i=1,nesp
	  DO j=1,NINT(zi_esp(i))
	   z_esp(i)=z_esp(i)+ioni_esp(j,i)*j	!nb e libres du à élément i
	   dz_esp(i)=dz_esp(i)+dioni_esp(j,i)*j
	  ENDDO
	 ENDDO
	 
c	 nombre d'électrons libres
	 
	 nel=sum(z_esp*x_esp)*ro/amu ; dnel=sum(dz_esp*x_esp)*ro/amu-1.d0
	 dnel=(nel-nelp)/dnel		!second terme de Newton-Raphson
	 nel=nelp-dnel			!nouveau nel
	 dnel=ABS(dnel/nel)		!précision relative
c	 IF(t < 4.d3)THEN
c	  WRITE(*,2000)t,ro,dnel,nel,nelp ; PAUSE'saha nel'
c	 ENDIF
	 nelp=nel
	ENDDO
	
	IF(tour >= tourmax .AND. no_conv)THEN
1	 FORMAT('pas de convergence dans saha au moins 1 fois',/,
	1't=',es10.3,', ro=',es10.3,', dnel=',es10.3)
	 WRITE(*,1)t,ro,dnel ; 	WRITE(2,1)t,ro,dnel ; no_conv=.FALSE.
	ENDIF
	
c	rétablissement des taux pour chaque isotope

	DO k=1,nesp	!pour chaque espèce
	 DO l=1,niso(k)	!pour chaque isotope d'espèce k
	  i=iso(k,l)	!indice de l'isotope l d'espèce k
	  z_bar(i)=z_esp(k)
	  DO j=0,NINT(zi(i))	!pour tous les ions de l'isotope i
	   ioni(j,i)=ioni_esp(j,k)  !memes taux pour les ions d'espèce k
	  ENDDO
	 ENDDO
	ENDDO
	
	DEALLOCATE(x_esp,z_esp,dz_esp,ioni_esp,dioni_esp,dphi)

	RETURN
			 
	END SUBROUTINE saha
