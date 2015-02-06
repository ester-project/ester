
C*****************************************************************

      SUBROUTINE diffm_br(t,r,lum,m,ro,drox,kap,dkapx,w,
     1 gradrad,dgradradx,xchim,d,dd,v,dv)

c     routine private du module diffm_br_mod

c     inversion des équations de Burgers avec accélérations radiatives

c     P est recalculée à partir de ro, T, mu
c     on tient compte approximativement des ionisations

             
c     Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k

c entrées
c     p,t,r,lum,ltot,m : pression, temperature, lum. locale et totale, masse
c     ro,drox,kap,dkapx : densite, dérivée/ X, opacite, dérivée / X 
c     gradrad,dgradradx,gradad,dgradadx : gradient rad., ad. et dérivées / X
c     xchim : comp. chim. par mole
c     les dérivées / X, issues de thermo, sont par gramme 

c sorties
c     d, dd : coefficients d_ij de d x_j / d m et dérivées / x_k
c     v, dv : coefficients v_i de x_i et dérivées / x_k

c variables non utilisees pour cette routine de dIFfusion
c     gradad,dgradadx

c     Dimensions dans le programme appelant
c     d(nchim,nchim), dd(nchim,nchim,nchim), v(nchim), dv(nchim,nchim)

c     les tableaux d, dd, v, dv doivent etre initialisés à 0 avant l'appel

c     convention de notation :
c     equation de diffusion dXi/dt=dFi/dm + nuclear, i=1,nchim
c     Fi=4 pi r**2 ro (4 pi r**2 ro D.dX/dm - Vi Xi)

c     d=D=(di1, di2,... din) avec Dij coefficient de d Xj / dm
c     dans le produit scalaire D.dX/dm=sum d_ij d Xj / dm

c     pour ligne d'indice i
c     v(i) coefficient de x_i,
c     dv(i,k)=dérivée v_i / x_k
c     d(i,j)=coefficient d_ij de d x_j / dm
c     dd(i,j,k)= dérivée de d_ij / x_k

c     on calcule les dérivées / x_k, k=1,nb, y compris pour les électrons x_e,
c     puis on considere que x_e=x_e(x_i), i=1,nchim

c-------------------------------------------------------------------

      USE mod_donnees, ONLY : amu, echarg, g, granr, kbol, langue,
     1 me, msol, nchim, nucleo, rsol, zi
      USE mod_etat, ONLY : saha
      USE mod_kind		
      USE mod_numerique, ONLY : matinv, pause

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim   
      REAL (kind=dp), INTENT(in) :: dkapx, drox, gradrad, dgradradx, kap,
     1 m, lum, r, ro, t, w
      REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd   
      REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d, dv
      REAL (kind=dp), INTENT(inout), DIMENSION(:) :: v    

      REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: da, db,
     1 dbqij, dcd, dcf, dfij, dkij, dpij, dqij, dom11, dzij, dzpij, dzsij        
      REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  a, b, bqij,
     1 cd, cf, ddiag, dgamma, dg_rad, dmu_nuir_gtsro, dnuir, fij, ioni, kij,
     2 mij, mijk, mismj, mi1, mi2, mi3, mt, om11, pij, qij,
     3 zij, zpij, zsij
      REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: dgt, dgtsro,
     1 diag, dmu, dmu_gtsro_grad, gamma, g_rad, mi, mu_nuir_gtsro, nui,
     2 nuir, unpzb, vt, xi, zb      
      REAL (kind=dp), SAVE :: cgrav, dro_n0_ex, n0, n02k, ro_n0_e
      REAL (kind=dp) :: bid, dgrad_x, dmu_px, drox1, eta, gt, gtsro,
     1 grav, mu, mu_gtsro_grad, mu_p, nel

      INTEGER, SAVE :: nb, nl
      INTEGER :: i, j, k, l

      LOGICAL, SAVE :: init=.TRUE.
      LOGICAL :: inversible           

c*******************pour tests de dérivation*******************   
c     REAL (kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: da0
c     REAL (kind=dp), DIMENSION(:,:), ALLOCATABLE :: a0   
c     REAL (kind=dp), parameter :: dx=1.d-5, unpdx=1.d0+dx
c     REAL (kind=dp) :: stor, stor0, dstor
c     INTEGER :: der
c*******************pour tests de dérivation*******************

c----------------------------------------------------------------------   

2000  FORMAT(8es10.3)
2001  FORMAT(10es8.1)

c     initialisations

      IF(init)THEN
       init=.FALSE.
       SELECT CASE(langue)
       CASE('english')
        WRITE(*,1010) ; WRITE(2,1010)      
1010    FORMAT(/,'diff. micro. of Burgers, with approached ionisation',/)       
       CASE DEFAULT
        WRITE(*,10) ; WRITE(2,10)      
10      FORMAT(/,'diff. micro. par eq. Burgers, dans ZR',/,
     1  'On tient compte approximativement de l''ionisation',/)
       END SELECT
       
       cgrav=-msol*g/rsol**2  !pour le calcul de la gravite
       n0=1.d0/amu        !n0 Avogadro     
       n02k=16.d0/3.d0*n0**2  !pour les K_ij
       nb=nchim+1 !nombre d'elements pour Burgers, nchim+1 (électron)
       nl=nchim+nb+2  !nombre de lignes de A, les wi, ri et E

       ALLOCATE(mij(nb,nb),mijk(nb,nb),mi1(nb,nb),mi2(nb,nb),mi3(nb,nb),
     1 mismj(nb,nb),nui(nb),mi(nb))

       nui(1:nchim)=nucleo(1:nchim)   !masses atomiques
       mi(1:nchim)=nui(1:nchim)*amu   !masses
       nui(nb)=me/amu ; mi(nb)=me !pour l'électron
              
c      masses reduites

       DO j=1,nb
        mi1(:,j)=mi(:)/(mi(:)+mi(j)) ; mij(:,j)=mi(:)*mi(j)/(mi(:)+mi(j))
        mismj(:,j)=-mi(:)/mi(j)   !- pour \bar q
       ENDDO
       DO j=1,nb   
        mi3(:,j)=mi1(:,j)*mi1(j,:) 
       ENDDO         
       mi2=mi1**2 ; mijk=mij*n02k !coefficient des k_ij

c      appel fictif à f_rad, avec grav=+100, pour init. et écritures

c	PRINT*,nint(maxval(zi)),nchim

       ALLOCATE(g_rad(nb),dg_rad(nb,nb),ioni(0:NINT(MAXVAL(zi)),nchim)) 
      
       CALL f_rad(lum,r,t,kap,dkapx/kap,nel,xchim,ioni,100.d0,
     1 g_rad,dg_rad)
           
       ALLOCATE(om11(nb,nb), dom11(nb,nb,nb), cd(nb,nb), dcd(nb,nb,nb),
     1 kij(nb,nb), dkij(nb,nb,nb), pij(nb,nb), dpij(nb,nb,nb),
     2 fij(nb,nb), dfij(nb,nb,nb),
     3 zij(nb,nb), dzij(nb,nb,nb), zpij(nb,nb), dzpij(nb,nb,nb),
     4 zsij(nb,nb), dzsij(nb,nb,nb), cf(nb,nb), dcf(nb,nb,nb),
     5 qij(nb,nb), dqij(nb,nb,nb), ddiag(nb,nb),
     6 bqij(nb,nb), dbqij(nb,nb,nb), mt(nl,nb),
     7 dnuir(nb,nb), dmu_nuir_gtsro(nb,nb), dgamma(nl,nb),
     8 a(nl,nl), da(nl,nl,nb), b(nl,nl), db(nl,nl,nb),
     9 gamma(nl), vt(nl), dmu(nb), dgt(nb), xi(nb),
     1 dgtsro(nb), diag(nb), nuir(nb), mu_nuir_gtsro(nb),
     2 zb(nb), unpzb(nb), dmu_gtsro_grad(nb))
      ENDIF   !fin des initialisations

c     transfert des abondances, xi = X_i / A_i, xi(nb)=0 pour X_e,    

      xi(1:nchim)=xchim(1:nchim)

c     appel à Saha, calcul des taux d'ionisation, des charges

      CALL saha(xi,t,ro,ioni,zb,nel,eta)
c      WRITE(*,2000)(zb(i),i=1,nchim),nel  
c      CALL pause('après saha')

c*******************pour tests de dérivation*******************
c     ALLOCATE(da0(nl,nl,nb),a0(nl,nl))
c     DO der=0,nchim
c     IF(der /= 0)THEN
c      stor0=xi(der)
c      stor=stor0*unpdx
c      dstor=stor-stor0
c      xi(der)=stor
c     ENDIF
c*******************pour tests de dérivation*******************   

c     transformation xi(nb): nombre d'e par volume ==> par mole

      unpzb(1:nchim)=1.d0+zb(1:nchim) !1+Zi
      xi(nb)=0.d0
      DO i=1,nchim        !xi(nb) : nombre d'électrons par mole
       xi(nb)=xi(nb)+zb(i)*xi(i)   
      ENDDO
      zb(nb)=-1.d0 ; unpzb(nb)=0.d0   !pour les électrons

c     quantites qui dependent des x_i
c     il y a incohérence entre mu, ro et P

      mu=sum(xi*unpzb)    !poids moleculaire moyen  on a unpzb(nb)=0
      mu=1.d0/mu ; dmu=-mu**2*unpzb
      mu_p=granr*t/ro         !R T / ro
      drox1=drox*nui(1)/ro        !dro / x_1 par mole à ro pres   
      dmu_px=-drox1           !a mu_p pres
      dgrad_x=dgradradx/gradrad*nui(1)

c     la gravite

      IF(r*m > 0.d0)THEN
       grav=cgrav*ABS(m)/r**2         !gravite
      ELSE
       grav=0.d0
      ENDIF		!ADAPTER
      grav=grav+2.d0/3.d0*r*w**2    !acc. centrifuge    

c     accelerations radiatives

      CALL f_rad(lum,r,t,kap,dkapx/kap,nel,xi,ioni,grav,
     1 g_rad,dg_rad)
c     WRITE(*,2001)g_rad
c     DO i=1,nb
c      WRITE(*,2001)dg_rad(i,:)
c     ENDDO
c     PAUSE'acc rad'      

c     gt=grav + somme ponderee des accelerations radiatives
c     g_rad(nb)=dg_rad(nb,k)=0 pour les ``électrons'' 

      gt=grav+sum(xi*nui*g_rad)
      dgt=nui*g_rad
      DO j=1,nb
       dgt(j)=dgt(j)+sum(xi*nui*dg_rad(:,j))
      ENDDO
      gtsro=gt/ro ; dgtsro=dgt/ro ; dgtsro(1)=dgtsro(1)-gtsro*drox1
      mu_gtsro_grad=mu*gtsro*gradrad
      dmu_gtsro_grad=mu_gtsro_grad*(dmu/mu+dgtsro/gtsro)
      dmu_gtsro_grad(1)=dmu_gtsro_grad(1)+mu_gtsro_grad*(dgrad_x-drox1)

      nuir=nui*(1.d0-g_rad/gt)
      DO j=1,nb
       dnuir(:,j)=-nui/gt*(dg_rad(:,j)-g_rad*dgt/gt)
      ENDDO
      
      mu_nuir_gtsro=(mu-nuir)*gtsro
      DO i=1,nb
       DO j=1,nb
        dmu_nuir_gtsro(i,j)=dmu(j)*gtsro
     1  -dnuir(i,j)*gtsro+(mu-nuir(i))*dgtsro(j)
       ENDDO
      ENDDO

c     WRITE(*,2000)mu,gt
c     WRITE(*,2000)nuir
c     WRITE(*,2000)mu-nuir
c     WRITE(*,2000)mu_nuir_gtsro
c     PAUSE'mu,nuir,mu-nui'
c     PRINT*,'dmu_nuir_gtsro' 
c     DO i=1,nb
c      WRITE(*,2000)(dmu_nuir_gtsro(i,j),j=1,nb)
c     ENDDO 
c     PAUSE'dérivées'

c     quantites diverses

      ro_n0_e=echarg*n0/ro
      dro_n0_ex=-ro_n0_e*drox1

c     CALL pause('quantites diverses')

c-------------test local de dérivation----------------------------------

c      pour un scalaire --------debut-------------------
c       IF(der == 0)THEN
c        a0(1,1)=mu_gtsro_grad
c        DO i=1,nchim
c         da0(1,1,i)=dmu_gtsro_grad(i)
c        ENDDO
c       ELSE
c        v(1)=mu_gtsro_grad
c        PRINT*,'dérivée inconnue ',nom_elem(der)
c        WRITE(*,2000)(v(1)-a0(1,1))/dstor,da0(1,1,der)
c        xi(der)=stor0
c        PAUSE
c       ENDIF
c      pour un scalaire --------fin------------------   

c      pour un vecteur --------debut-------------------
c      IF(der == 0)THEN    
c       DO i=1,nchim
c        a0(1,i)=mu_nuir_gtsro(i) 
c        DO j=1,nchim
c         da0(i,j,1)=dmu_nuir_gtsro(i,j)      
c        ENDDO       
c       ENDDO
c      ELSE
c       DO i=1,nchim
c        v(i)=mu_nuir_gtsro(i)
c       ENDDO
c       PRINT*,'dérivée inconnue ',nom_elem(der)
c       WRITE(*,2000)((v(i)-a0(1,i))/dstor,da0(i,der,1),i=1,nchim)
c       PAUSE'test local de dérivation'
c       xi(der)=stor0
c      ENDIF   
c      ENDDO  !boucle sur les xi   
c      pour un vecteur --------fin-------------------

c-------------fin du test local de dérivation----------------------------------   

c     integrales de collision

      CALL collision(nb,zb,mij,xi,ro,drox1,t,om11,zij,zpij,zsij,
     1 dom11,dzij,dzpij,dzsij)

      DO i=1,nb   !coefficients pour d et f
       DO j=1,nb
        cd(i,j)=2.d0/5.d0*xi(j)*(3.d0*mi2(i,j)+zpij(i,j)*mi2(j,i)
     1  +4.d0/5.d0*zsij(i,j)*mi3(i,j))
        cf(i,j)=-2.d0/5.d0*xi(j)*mi3(i,j)*(3.d0+zpij(i,j)-4.d0/5.d0*zsij(i,j))
        DO k=1,nb
         dcd(i,j,k)=2.d0/5.d0*xi(j)*(dzpij(i,j,k)*mi2(j,i)
     1   +4.d0/5.d0*dzsij(i,j,k)*mi3(i,j))
         dcf(i,j,k)=-2.d0/5.d0*mi3(i,j)*xi(j)*(dzpij(i,j,k)
     1   -4.d0/5.d0*dzsij(i,j,k))
        ENDDO
        dcd(i,j,j)=dcd(i,j,j)+cd(i,j)/xi(j)
        dcf(i,j,j)=dcf(i,j,j)+cf(i,j)/xi(j)     
       ENDDO
      ENDDO
       
c     kij, dérivées

      DO i=1,nb
       DO j=i,nb
        kij(i,j)=mijk(i,j)*om11(i,j)    
        IF(j /= i)kij(j,i)=kij(i,j)
        DO k=1,nb
         dkij(i,j,k)=mijk(i,j)*dom11(i,j,k)
         IF(j /= i)dkij(j,i,k)=dkij(i,j,k)
        ENDDO
       ENDDO
      ENDDO
       
      DO i=1,nb
       diag(i)=4.d0/25.d0*kij(i,i)*zsij(i,i)*xi(i)
       DO k=1,nb
        ddiag(i,k)=diag(i)*(dkij(i,i,k)/kij(i,i)+dzsij(i,i,k)/zsij(i,i))
       ENDDO
       ddiag(i,i)=ddiag(i,i)+diag(i)/xi(i)    !pour k=i 
       DO j=1,nb
        pij(i,j)=kij(i,j)*xi(j)
        DO k=1,nb
         dpij(i,j,k)=dkij(i,j,k)*xi(j)
        ENDDO
        dpij(i,j,j)=dpij(i,j,j)+kij(i,j)

        qij(i,j)=pij(i,j)*zij(i,j)*mi1(j,i)
        DO k=1,nb
         dqij(i,j,k)=qij(i,j)*(dpij(i,j,k)/pij(i,j)+dzij(i,j,k)/zij(i,j))
        ENDDO

        bqij(i,j)=qij(i,j)*mismj(i,j)
        DO k=1,nb
         dbqij(i,j,k)=dqij(i,j,k)*mismj(i,j)
        ENDDO

        IF(i /= j)THEN
         diag(i)=diag(i)+kij(i,j)*cd(i,j)
         DO k=1,nb
          ddiag(i,k)=ddiag(i,k)+dkij(i,j,k)*cd(i,j)+kij(i,j)*dcd(i,j,k)
         ENDDO
        ENDIF
        fij(i,j)=kij(i,j)*cf(i,j)
        DO k=1,nb
         dfij(i,j,k)=dkij(i,j,k)*cf(i,j)+kij(i,j)*dcf(i,j,k)
        ENDDO
       ENDDO
      ENDDO
      
c     FORMATion du systeme

c     mise à 0

      a=0.d0 ; da=0.d0

c     la sous matrice Ae

      DO i=1,nchim
       a(i,nl)=ro_n0_e*zb(i)      !derniere colonne : matrice Ae
       da(i,nl,1)=dro_n0_ex*zb(i) !dérivée k=1, ie /X  
      ENDDO
          
c     la sous matrice Ac

      DO j=1,nchim
       a(nl-1,j)=nui(j)*xi(j) !avant derniere ligne : matrice Ac
       da(nl-1,j,j)=nui(j)    !dérivée k=j
      ENDDO
      DO j=1,nb
       a(nl,j)=zb(j)*xi(j)    !derniere ligne : matrice Ac
       da(nl,j,j)=zb(j)   !dérivée k=j         
      ENDDO
      
c     la sous matrice Aw

      DO i=1,nchim    !partie superieure
       DO j=1,nb
        IF(j == i)THEN        !on est sur la diagonale
         DO l=1,nb    !somme pour les indices l non "diagonaux"
          IF(l /= j)THEN
           a(i,j)=a(i,j)-pij(i,l)
           DO k=1,nb
            da(i,j,k)=da(i,j,k)-dpij(i,l,k)
           ENDDO
          ENDIF
         ENDDO    !l
        ELSE
         a(i,j)=pij(i,j)
         DO k=1,nb
          da(i,j,k)=dpij(i,j,k)
         ENDDO
        ENDIF
       ENDDO  !j
      ENDDO   !i
      DO i=nb,nchim+nb    !partie inferieure
       DO j=1,nb
        IF(j == i-nchim)THEN  !on est sur la "diagonale"
         DO l=1,nb        !somme pour les indices l non "diagonaux"
          IF(l /= j)THEN
           a(i,j)=a(i,j)-qij(i-nchim,l)
           DO k=1,nb
            da(i,j,k)=da(i,j,k)-dqij(i-nchim,l,k)
           ENDDO
          ENDIF
         ENDDO
        ELSE
         a(i,j)=qij(i-nchim,j)
         DO k=1,nb
          da(i,j,k)=dqij(i-nchim,j,k)    
         ENDDO
        ENDIF
       ENDDO  !j
      ENDDO   !i

c     la sous matrice Ar

      DO i=1,nchim    !partie superieure
       DO j=1,nb
        IF(j == i)THEN        !la "diagonale"
         a(i,nb+j)=-a(i+nchim,j)  !anti symetrie
         DO k=1,nb
          da(i,nb+j,k)=-da(i+nchim,j,k)   !idem
         ENDDO
        ELSE
         a(i,nb+j)=bqij(i,j)
         DO k=1,nb
          da(i,nb+j,k)=dbqij(i,j,k)
         ENDDO
        ENDIF
       ENDDO  !j
      ENDDO   !i
      DO i=nb,nchim+nb        !partie inferieure
       DO j=1,nb
        IF(j == i-nchim)THEN  !la "diagonale"
         a(i,nb+j)=diag(j)
         DO k=1,nb
          da(i,nb+j,k)=ddiag(j,k)
         ENDDO      
        ELSE    
         a(i,nb+j)=fij(i-nchim,j)
         DO k=1,nb       
          da(i,nb+j,k)=dfij(i-nchim,j,k)
         ENDDO       
        ENDIF
       ENDDO  !j
      ENDDO   !i

c     PRINT*,'debut A avant inversion'
c     DO i=1,nl
c      WRITE(*,2000)(a(i,j),j=1,nl)
c     ENDDO
c     PAUSE'fin A avant inversion'

c     inverse A**(-1)

      CALL matinv(a,nl,inversible)
      IF(.not.inversible)PAUSE'matrice non inversible'

c     CALL svd_cond(a,nl,nl,bid)  
c     CALL svd_inv(a,nl)      
c     CALL minv(a,nl) 

c     PRINT*,'A apres inversion'
c     DO i=1,nl
c      WRITE(*,2000)a(i,:)
c     ENDDO
c     CALL pause('A apres inversion')
           
c     vecteur gamma

      gamma=0.d0 ; dgamma=0.d0
      gamma(1:nchim)=mu_nuir_gtsro(1:nchim)
      DO i=1,nchim        !pour les nchim premières lignes
       dgamma(i,1:nb)=dmu_nuir_gtsro(i,1:nb)
      ENDDO
      gamma(nb:nl-2)=-mu_gtsro_grad   !pour les nb lignes suivantes
      DO i=nb,nl-2        !pour les nb lignes suivantes       
       dgamma(i,1:nb)=-dmu_gtsro_grad(1:nb)
      ENDDO

c     on extrait vt = A^(-1) gamma

      vt=matmul(a,gamma)

c     pour V on prend les nchim premières composantes

      v(1:nchim)=vt(1:nchim)

c     dérivées des nchim composantes de V

      DO k=1,nb
       DO i=1,nl          !d gamma / dx_k - dA/dx_k  V
        DO j=1,nl
         dgamma(i,k)=dgamma(i,k)-da(i,j,k)*vt(j)
        ENDDO
       ENDDO
       DO i=1,nchim   !produit A^(-1) (d gamma / dx_k - dA/dx_k  V)
        bid=0.d0  !il suffit des nchim premières lignes
        DO j=1,nl
         bid=bid+a(i,j)*dgamma(j,k)
        ENDDO
        IF(k <= nchim)THEN
         dv(i,k)=bid
        ELSE          !k=nb ie pour les électrons
         DO j=1,nchim     !d v_i /d x_j + d v_i / d x_e d x_e/ dx_j
          dv(i,j)=dv(i,j)+bid*zb(j)
         ENDDO
        ENDIF
       ENDDO
      ENDDO

c-------------test local de dérivation----------------------------------

c       IF(der == 0)THEN   
c        DO i=1,nchim
c         a0(1,i)=v(i) 
c         DO j=1,nchim
c          da0(i,j,1)=dv(i,j)     
c         ENDDO
c        ENDDO
c       ELSE    
c        PRINT*,'dérivée inconnue ',nom_elem(der)
c        WRITE(*,2000)((v(i)-a0(1,i))/dstor,da0(i,der,1),i=1,nchim)
c        xi(der)=stor0       
c        PAUSE'test local de dérivation'     
c       ENDIF
c      ENDDO  !boucle sur les xi 

c-------------fin du test local de dérivation----------------------------------   

c     matrice B :  R T / ro ( Diag (1/x_i) - mu M )

      b=0.d0 ; mt=0.d0 ; db=0.d0

      DO i=1,nchim
       b(i,i)=1.d0/xi(i)          !Diag (1/x_i)
       db(i,i,i)=-b(i,i)/xi(i)    !der/ x_i
       DO j=1,nchim
        b(i,j)=mu_p*(b(i,j)-mu*unpzb(j))  !- mu M
        DO k=1,nchim
         db(i,j,k)=mu_p*(db(i,j,k)-dmu(k)*unpzb(j))
        ENDDO
        db(i,j,1)=db(i,j,1)+dmu_px*b(i,j) !k=1=X  
       ENDDO
      ENDDO

c     PRINT*,'B'
c     DO i=1,nchim
c      WRITE(*,2001)(b(j,i),j=1,nl)
c     ENDDO
c     PAUSE

c--------------debut test local de dérivation sur B-----------------------

c     IF(der == 0)THEN
c      DO i=1,nchim
c       DO j=1,nchim
c        a0(i,j)=b(i,j)
c        DO k=1,nchim
c         da0(i,j,k)=db(i,j,k)
c        ENDDO
c       ENDDO  
c      ENDDO
c     ELSE
c      DO i=1,nchim
c       PRINT*,'B ligne',i,' / element ',nom_elem(der)
c      WRITE(*,2000)((b(i,j)-a0(i,j))/dstor,
c     1   da0(i,j,der),j=1,nchim)
c       PAUSE'test de dérivation sur B'
c      ENDDO
c      xi(der)=stor0
c     ENDIF
c     ENDDO

c--------------fin test local de dérivation sur B-----------------------

c     calcul de mt = A^(-1) B

      DO i=1,nl
       DO j=1,nchim           !nchim colonnes suffisent
        mt(i,j)=0.d0
        DO l=1,nl
         mt(i,j)=mt(i,j)+a(i,l)*b(l,j)
        ENDDO
       ENDDO
      ENDDO

c     on extrait G, nchim premières lignes et colonnes de G = mt =  A^(-1) B
c     mis dans D

      d(1:nchim,1:nchim)=mt(1:nchim,1:nchim)

c     dérivées des nchim X nchim premières composantes de G 

      DO k=1,nb
       DO i=1,nl          !d B / dx_k - dA/dx_k  G
        DO j=1,nchim          !nchim colonnes suffisent
         DO l=1,nl
          db(i,j,k)=db(i,j,k)-da(i,l,k)*mt(l,j)
         ENDDO
        ENDDO
       ENDDO
       DO i=1,nchim   !produit A^(-1) (d B / dx_k - dA/dx_k  G)
        DO j=1,nchim  !il suffit des nchim premières lignes/colonnes
         vt(j)=0.d0     
         DO l=1,nl
          vt(j)=vt(j)+a(i,l)*db(l,j,k)
         ENDDO
        ENDDO
        IF(k <= nchim)THEN
         DO j=1,nchim
          dd(i,j,k)=vt(j)
         ENDDO
        ELSE          !k=nb ie pour les électrons
         DO j=1,nchim
          DO l=1,nchim      
           dd(i,j,l)=dd(i,j,l)+vt(j)*zb(l)
          ENDDO
         ENDDO          
        ENDIF
       ENDDO
      ENDDO

c     d*xi(i) : notation de CESAM
c     pour d, dd nchim X nchim  d, dd

      DO i=1,nchim
       DO j=1,nchim
        d(i,j)=-d(i,j)*xi(i)
        DO k=1,nchim
         dd(i,j,k)=-dd(i,j,k)*xi(i)
        ENDDO
        dd(i,j,i)=dd(i,j,i)+d(i,j)/xi(i)      !pour k=i     
       ENDDO
      ENDDO

c****************************pour test de dérivation********************
c     WRITE(*,2000)v(1:nchim)
c     WRITE(*,*)
c     DO i=1,nchim    
c      WRITE(*,2000)(dv(i,j),j=1,nchim)
c     ENDDO
c     WRITE(*,*)  
c     DO i=1,nchim    
c      WRITE(*,2000)(d(i,j),j=1,nchim)
c     ENDDO
c     PAUSE

c      IF(der == 0)THEN
c       DO i=1,nchim
c        a0(1,i)=v(i)
c        DO j=1,nchim
c         da0(i,j,1)=dv(i,j)
c        ENDDO
c       ENDDO
c      ELSE
c       PRINT*,'dérivée / ',nom_elem(der)
c       WRITE(*,2000)((v(i)-a0(1,i))/dstor,da0(i,der,1),i=1,nchim)
c       xi(der)=stor0
c       PAUSE'test de dérivation sur V'
c      ENDIF
c     ENDDO

c     IF(der == 0)THEN
c      DO i=1,nchim
c       DO j=1,nchim
c        a0(i,j)=d(i,j)
c        DO k=1,nchim
c         da0(i,j,k)=dd(i,j,k)
c        ENDDO
c       ENDDO
c      ENDDO
c     ELSE
c      DO i=1,nchim
c       PRINT*,'D ligne',i,' / element ',nom_elem(der)
c       WRITE(*,2000)((d(i,j)-a0(i,j))/dstor,da0(i,j,der),j=1,nchim)
c       PAUSE'test de dérivation sur D'
c      ENDDO
c      xi(der)=stor0
c     ENDIF
c     ENDDO

c****************************pour test de dérivation********************

	RETURN

	END SUBROUTINE diffm_br
