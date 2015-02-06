
c************************************************************************ 

      SUBROUTINE alecian1(lum,ray,t,kap,dkapx,nel,ychim,ioni,grav,
     1 g_rad,dg_rad)

c     routine private du module mod_evol
c     Version 1.0 (24/11/1999)
c     cette subroutine calcule l'accélération radiative et ses dérivées
c     par rapport aux l'abondances. Elle doit être appelée pour chaque
c     couche, et à chaque fois que la diffusion est calculée.

c     Auteur:
c     Georges ALECIAN
c     DAEC - UMR 8631, Observatoire de Meudon
c     F-92195 MEUDON CEDEX, FRANCE
c     Tel: 01 45 07 74 20, + 33 1 45 07 74 20

c     Adaptation a CESAM2k par P. Morel, OCA

c entrées:
c     lum : luminosité locale
c     ray: rayon
c     t : température
c     kap : opacité
c     dkapx : d kap / d X    / kap (ie. kap pres)
c     nel : nombre électrons par volume,
c     ychim : composition chimique/mole
c     ioni : taux d'ionisation
c     grav : gravité

c sorties
c     g_rad(i) : accélérations radiatives la gravité est grav+g_rad sur i
c     dg_rad(i,j) : dérivée de g_rad(i) / ion j

c--------------------------------------------------------------

      USE mod_donnees, ONLY : amu, ap, aradia, clight, echarg,
     1 hpl, kbol, lsol, me, nchim, nom_elem, nucleo, pi, rsol, zi
      USE mod_kind
      USE mod_numerique, ONLY : bsp1dn, no_croiss
      USE mod_variables, ONLY : mstar

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(0:,:) :: ioni
      REAL (kind=dp), INTENT(in), DIMENSION(:) :: ychim
      REAL (kind=dp), INTENT(in) :: grav, dkapx, kap, lum, nel, ray, t
      REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: dg_rad
      REAL (kind=dp), INTENT(out), DIMENSION(:) :: g_rad

      INTEGER, PARAMETER :: pnspin=158, pnlist=8, m=4

      REAL (kind=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: d_ciony
      REAL (kind=dp), SAVE, DIMENSION(:,:), ALLOCATABLE :: phipsi2
      REAL (kind=dp), SAVE, DIMENSION(:,:), ALLOCATABLE :: c_ion,d_chimy
      REAL (kind=dp), SAVE, DIMENSION(:), ALLOCATABLE :: c_chim, df,
     1 d_grkj, f, m_spin, m_spint
      REAL (kind=dp), SAVE :: cte1, cte2, mstarp=-1.d0
      REAL (kind=dp) :: psi2, phi, gr_kj, q, b, q0, ci_s

      INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: isotops
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nps,izs
      INTEGER, DIMENSION(:), ALLOCATABLE :: numions
      INTEGER, SAVE :: knot, mz, nbf, nb
      INTEGER :: i, j, k, nf

      LOGICAL, DIMENSION(:), SAVE, ALLOCATABLE :: iden
      LOGICAL, SAVE :: init=.TRUE., lg_rad, lgtot0, no_iden

      CHARACTER (len=4), DIMENSION(:,:), ALLOCATABLE :: listis
      CHARACTER (len=4), DIMENSION(:), ALLOCATABLE :: tampon
      CHARACTER (len=80) :: nom_chemin = "/data1/sdeheuve/local/src/cesam2k_v1.1.8_ESTA/SUN_STAR_DATA/"
      CHARACTER (len=80), DIMENSION(:), ALLOCATABLE :: list

c-------------------------------------------------------------------- 

c     lecture des options, des donnees et
c     preparation des tableaux necessaires au calcul des
c     accélérations (modul_grad)

c     nps: nb de protons du noyau
c     izs: degre d'ionisation (0 pour neutre)
c     listis: nom des ions (c1 = carbone neutre)
c     numion: nb d'ions dans la base

2000  FORMAT(8es10.3)

      IF(init)THEN
       init=.FALSE.
       WRITE(*,20) ; WRITE(2,20)
20     FORMAT('Détermination des accélérations radiatives',/,
     1 'utilisation des algorithmes et tables de G.Alécian')

c      si iden(i) (pour les ions identifiés) 
c          si lg_rad  on calcule g_rad
c          sinon si lgtot0    g_rad=gravite ie. gtot=0
c            sinon    g_rad=0 ie. gtot=gravite
c      sinon  (pour les ions non identifiés)
c          si lgtot0  g_rad=gravite ie. gtot=0
c          sinon  g_rad=0 ie. gtot=gravite

       lg_rad=.TRUE. ; lgtot0=.FALSE.

c      nombre MAX. d'ions

       mz=0
       DO i=1,nchim
        mz=MAX(mz,NINT(zi(i)))
       ENDDO
       mz=mz+1 ; nbf=2*nchim*mz !nombre de fonctions phi, psi2

c      Tableau des isotopes de cESAM, matrice (isotops) du genre:

c             H1   He3  He4  c12  c13  N14  N15  O16  O17  Ex
c      H1     1    0    0    0    0    0    0    0    0    0
c      He3    0    1    1    0    0    0    0    0    0    0
c      He4    0    1    1    0    0    0    0    0    0    0
c      c12    0    0    0    1    1    0    0    0    0    0
c      c13    0    0    0    1    1    0    0    0    0    0
c      N14    0    0    0    0    0    1    1    0    0    0
c      N15    0    0    0    0    0    1    1    0    0    0
c      O16    0    0    0    0    0    0    0    1    1    0
c      O17    0    0    0    0    0    0    0    1    1    0
c      Ex     0    0    0    0    0    0    0    0    0    1

       ALLOCATE(isotops(nchim,nchim))
       isotops=0
       DO j=1,nchim
        DO k=1,nchim
         IF(NINT(zi(j)) == NINT(zi(k)))isotops(j,k)=1
        ENDDO
       ENDDO

c      liste des fichiers contenant les tables phi & psi2

       ALLOCATE(list(pnlist),m_spin(pnlist),m_spint(pnlist+m))
       list(1)='phi_psi2_100.data' ; m_spin(1)=1.d0
       list(2)='phi_psi2_110.data' ; m_spin(2)=1.1d0
       list(3)='phi_psi2_120.data' ; m_spin(3)=1.2d0
       list(4)='phi_psi2_130.data' ; m_spin(4)=1.3d0
       list(5)='phi_psi2_140.data' ; m_spin(5)=1.4d0
       list(6)='phi_psi2_150.data' ; m_spin(6)=1.5d0
       list(7)='phi_psi2_175.data' ; m_spin(7)=1.75d0
       list(8)='phi_psi2_200.data' ; m_spin(8)=2.d0

       IF(.NOT.(mstar >= m_spin(1) .AND. mstar <= m_spin(pnlist)))THEN
        PRINT*
        PRINT*,'mstar hors des tables disponibles => STOP'
        WRITE(*,21)mstar,m_spin(1),m_spin(pnlist)
        WRITE(2,21)mstar,m_spin(1),m_spin(pnlist)   
21      FORMAT(1x,'mstar=',es10.3,' hors de [',es10.3,',',es10.3,']')
        STOP
       ENDIF

c      lecture de la base des phi et psi2=psi**2
c      suivant nouveau FORMAT par le programme réécrit_phi_psi2
c      pour chaque ion on interpole en fonction de la masse de l'*
c      les quantites phi et phi2
c      il y a moyen de n'interpoler que les données relatives
c      aux éléments utilisés et non pas toutes ainsi qu'il est fait

c      nps(j,i) est la charge de l'élément, repétée pour tous les isotopes
c      et tous les fichiers list
c      izs(j,i) est le taux d'ionisation, k le nombre de charges de l'ion
c      k=nps(j,i)-izs(j,i), il y a charge-1=nps(j,i)-1 donnees par élément

       ALLOCATE(nps(pnspin,pnlist),izs(pnspin,pnlist),listis(pnspin,pnlist),
     1 phipsi2(nbf,pnlist),df(nbf),f(nbf),numions(pnlist))

       phipsi2=0.d0
       DO i=1,pnlist
c        PRINT*,TRIM(nom_chemin)//list(i)
        OPEN(UNIT=70,FILE=TRIM(nom_chemin)//list(i),STATUS='old')
c       PRINT*,list(i)
        READ(70,*); READ(70,*)
        DO j=1,pnspin
         read(70,11,end=100)nps(j,i),izs(j,i),phi,psi2,k,listis(j,i)
c        WRITE(*,11)nps(j,i),izs(j,i),phi,psi2,k,listis(j,i)
11       FORMAT(2i4,2es10.3,i4,4x,a4)      
         DO k=1,nchim
          IF(nps(j,i) == NINT(zi(k)))THEN
           nf=2*(nchim*izs(j,i)+k-1)  !phipsi2(2,nchim,pzi+1,pnlist)
           phipsi2(nf+1,i)=phi ; phipsi2(nf+2,i)=psi2
          ENDIF
         ENDDO
        ENDDO
100     close(unit=70)
        numions(i)=j-1
       ENDDO
c      CALL pause('nmions')

c      interpolation des phi et psi2 en mstar

       CALL bsp1dn(nbf,phipsi2,m_spin,m_spint,pnlist,m,knot,.FALSE.,mstar,
     1  i,f,df)
        IF(no_croiss)PRINT*,'Pb. dans alecian1'
       mstarp=mstar

c      cross-identification on identifie les ions
c      On supposera que les propriétés atomiques (pour les g_rad) sont 
c      les mêmes entre isotopes d'un même élément.

c      dans les tables, nps est le nombre de protons ie. la charge
c      on identifie l'élément par sa charge
c      iden(i)=.TRUE. l'ion i est dans la table de Georges
c      c2 est le c ionise izs=1 fois, c3 est le c ionise izs=2 fois, etc..
c      izs est le degre d'ionisation 

       ALLOCATE(iden(nchim),tampon(nchim))

       iden=.FALSE.
c      PRINT*,'numion',numions(1)
c      WRITE(*,"(10(1x,i3))")(nps(k,1),k=1,numions(1))
       DO j=1,nchim
        DO k=1,numions(1)
         iden(j)=iden(j) .or. nint(zi(j)) == nps(k,1)
        ENDDO
c       PRINT*,j,' ',iden(j),' ',NINT(zi(j))
       ENDDO
c      CALL pause('numoins')

c      si on ne veut pas calculer d'accélération radiative pour l'élément
c      d'indice i imposer iden(i)=.FALSE.
c      l'élément est alors considéré comme non identifié

c      iden(nchim)=.FALSE.    !pour l'élément Ex

c      DO j=1,nchim
c       IF(iden(j))THEN
c        PRINT*,nom_elem(j),NINT(zi(j))
c        WRITE(*,2000)(f(2*(nchim*k+j-1)+2),k=0,NINT(zi(j)))
c       ENDIF
c      ENDDO
c      CALL pause('iden') 

c      informations concernant l'identification des ions

       i=0        !ions identifiés
       DO j=1,nchim
        IF(iden(j))THEN
         i=i+1
         tampon(i)=nom_elem(j)
        ENDIF
       ENDDO
       IF(i >= 1)THEN 
        WRITE(*,1)(tampon(j),j=1,i) ; WRITE(2,1)(tampon(j),j=1,i)
1       FORMAT(/,'pour les isotopes : ',/,30(a4,' '))
        IF(lg_rad)THEN
         WRITE(*,2) ; WRITE(2,2)
2        FORMAT('calcul des accélérations radiatives')
        ELSE
         IF(lgtot0)THEN
          WRITE(*,3) ; WRITE(2,3)
3         FORMAT('accélérations radiatives egales a -gravite')    
         ELSE
          WRITE(*,4) ; WRITE(2,4)
4         FORMAT('accélérations radiatives nulles')
         ENDIF
        ENDIF
       ELSE
        WRITE(*,6) ; WRITE(2,6)
6       FORMAT('aucun isotopes n''a été identifié')
        no_iden=.TRUE.
       ENDIF

       i=0        !ions non identifiés
       DO j=1,nchim
        IF(.NOT.iden(j))THEN
         i=i+1 ; tampon(i)=nom_elem(j)
        ENDIF
       ENDDO
       IF(i >= 1)THEN
        WRITE(*,5)(tampon(j),j=1,i) ; WRITE(2,5)(tampon(j),j=1,i)
5       FORMAT(/,'pour les isotopes : ',/,30(a4,' '))
        IF(lgtot0)THEN
         WRITE(*,3) ; WRITE(2,3)
        ELSE
         WRITE(*,4) ; WRITE(2,4)
        ENDIF
       ELSE
        WRITE(*,7) ; WRITE(2,7)
7       FORMAT('tous les isotopes ont été identifiés')
       ENDIF
       WRITE(*,*) ; WRITE(2,*)

c      controle des fichiers lus

       DO i=2,pnlist
        IF(numions(i-1) /= numions(i)) THEN
        WRITE(*,30)i-1,i,listis(j-1,i),listis(j,i)
        WRITE(2,30)i-1,i,listis(j-1,i),listis(j,i)      
30      FORMAT(/,"Anomalie dans le nombre d'ions entre les tables",
     1  i3,' et',i3,2(1x,2a4),' => STOP')
         STOP
        ENDIF
        DO j=1,numions(i)
         IF(listis(j,i-1) /= listis(j,i)) THEN
          WRITE(*,31)j-1,j,i-1,i,listis(j-1,i),listis(j,i)
          WRITE(2,31)j-1,j,i-1,i,listis(j-1,i),listis(j,i)       
31        FORMAT(1x,'Anomalie dans les noms des ions',i4,' et',i4,
     1    ' entre les tables',i3,' et',i3,2(1x,2a4),' => STOP')
          STOP
         ENDIF
         IF(izs(j,i-1) /= izs(j,i)) THEN
          WRITE(*,32)j-1,j,i-1,i,listis(j-1,i),listis(j,i)
          WRITE(2,32)j-1,j,i-1,i,listis(j-1,i),listis(j,i)       
32        FORMAT(1x,'Anomalie dans les taux izs des ions',i4,' et',i4,
     1    ' entre les tables',i3,' et',i3,2(1x,2a4),' => STOP')
          STOP
         ENDIF
         IF(nps(j,i-1) /= nps(j,i)) THEN
          WRITE(*,33)j-1,j,i-1,i,listis(j-1,i),listis(j,i)
          WRITE(2,33)j-1,j,i-1,i,listis(j-1,i),listis(j,i)       
33        FORMAT(1x,'Anomalie dans les charges nps des ions',i4,' et',i4,
     1   ' entre les tables',i3,' et',i3,2(1x,2a4),' => STOP')
          STOP
         ENDIF
        ENDDO
       ENDDO

c      CALL pause('anomalie')

c      les deux constantes de Georges      

       cte1=pi**2*kbol**3/hpl**2*echarg**2/(2.d0*clight**4*me*amu)
       cte2=me*ap*clight*amu/2.d0/echarg**2

c              pi**2 * kbol**3 * echarg**2
c      cte1=5.57E-05 = ---------------------------------
c              2.* clight**4 * hpl**2 * me * amu

c      exactement : cte1=5.585427439162097E-005

c              me * mp* clight
c      cte2=9.83E-23 = ---------------
c              2. * echarg**2

c      exactement : cte2=9.899551080572617E-023

c      PRINT*,'cte1,cte2' ; PRINT*,cte1,cte2 ; CALL pause('cte1,cte2')

c      cte1=5.57d-5       !valeurs de Georges pour tests
c      cte2=9.83d-23

c      on divise cte1 par 4 pi sigma = pi aradia clight
c      pour eviter les calculs de Teff^4 Rtot^2 Rsol^2 / lsol

       cte1=cte1/pi/aradia/clight/rsol**2*lsol ; nb=nchim+1

       DEALLOCATE(list,tampon,listis,numions,nps,izs)

       i=MAXVAL(zi)
       ALLOCATE(c_chim(nchim),c_ion(0:i,nchim),d_chimy(nchim,nchim),
     1 d_ciony(0:i,nchim,nb),d_grkj(nb))

      ENDIF   !initialisations

c-------  calcul des accélérations radiatives    ------------------------

c     initialisations, retour si grav=-Gm/r^2 > 0
c     ou si on ne tient pas compte des accélérations radiatives

      g_rad=0.d0 ; dg_rad=0.d0 ; IF(grav >= 0.d0)RETURN

c     tests sur l'eventualite du calcul des accélérations radiatives

      IF(no_iden)THEN
       IF(lgtot0)g_rad=-grav  !les acc. rad. equilibrent la gravite
       RETURN 
      ENDIF

c     si mstar a change de 2% reactualisation des interpolations en mstar

      IF(ABS(1.d0-mstar/mstarp) > 2.d-2)THEN
       CALL bsp1dn(nbf,phipsi2,m_spin,m_spint,pnlist,m,knot,.TRUE.,mstar,
     1 i,f,df)
       mstarp=mstar
      ENDIF

      if_lg_rad: IF(lg_rad)THEN

c      On considère que les isotopes ont leurs raies presque aux mêmes
c      fréquences. Ce qui implique que l'effet de saturation doit être
c      presque le même pour tous les isotopes d'un même élément.
c      Cette saturation est réglée par la concentration c_ion(k,j)
c      (en pratique: abondance par rapport a H). On veut donc que
c      c_ion(k,j) soit le même pour tous les isotopes d'un élément
c      et soit donné par la concentration totale pour l'élément.

c      d_chimy(i,j)   d c_chim(i) / d ychim(j)
c      d_ciony(i,j,k) d c_ion(i,j) /d ychim(k)
c      d_grkj(j)  d gr_kj / d ychim(j)

       d_chimy=isotops ; d_ciony=0.d0 ; c_chim=0.d0

c      calcul des abondances sommées sur tous les isotopes

       DO i=1,nchim
        c_chim(:)=c_chim(:)+ychim(i)*isotops(:,i)
       ENDDO
c      WRITE(*,2000)(c_chim(j),j=1,nchim)

c      calcul des abondances des ions suivant les taux d'ionisation

c      PRINT*,'c_ion acc'
       DO i=1,nchim
        DO k=0,NINT(zi(i))
         c_ion(k,:)=ioni(k,:)*c_chim(:)
         DO j=1,nchim
          d_ciony(k,:,j)=ioni(k,:)*d_chimy(:,j) !d c_ion / d ychim
         ENDDO    !j
        ENDDO     !k   
c       WRITE(*,2000)(c_ion(k,i),k=0,NINT(zi(i)))
       ENDDO      !i

c      calcul des accélérations radiatives et des dérivées
c      pour les ions identifiés        

       q0=cte1*lum/t/ray**2 ; b=cte2*nel/sqrt(t)*kap/nucleo(1)
       b_i: DO i=1,nchim
        IF(.NOT.iden(i))THEN  !pour les éléments non identifiés
         IF(lgtot0)g_rad(i)=0.d0 ; CYCLE b_i
        ENDIF
        q=q0/nucleo(i)
c       PRINT*,nom_elem(i),i ; WRITE(*,2000)q
        b_k: DO k=0,NINT(zi(i))
         nf=2*(nchim*k+i-1)
         IF(f(nf+2) < 1.E-30)CYCLE b_k
         ci_s=b*f(nf+2)       !psi**2 avec le carré interpolé
         gr_kj=q*f(nf+1)/sqrt(1.d0+c_ion(k,i)/ci_s)   !phistar interpolé
         d_grkj(:)=-gr_kj/2.d0/(1.d0+c_ion(k,i)/ci_s)*d_ciony(k,i,:)/ci_s
         d_grkj(1)=d_grkj(1)+gr_kj/2.d0/(1.d0+c_ion(k,i)/ci_s)*
     1   c_ion(k,i)/ci_s*dkapx
         g_rad(i)=g_rad(i)+ioni(k,i)*gr_kj
c        PRINT*,i,k ; WRITE(*,2000)g_rad(i),gr_kj,ioni(k,i)
         dg_rad(i,:)=dg_rad(i,:)+ioni(k,i)*d_grkj(:)!d g_rad / d ychim
        ENDDO b_k
       ENDDO b_i
      ENDIF if_lg_rad

      RETURN

      END SUBROUTINE alecian1
