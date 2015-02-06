
C******************************************************************

      SUBROUTINE etat_gong2(pp,tt,xchim,ro,drop,drot,drox,u,dup,dut,dux,
     1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c     routine public du module mod_etat

c     équation d'état pour GONG étape 2
c     d'après la note de J.Christensen-Dalsgard 1 / 3 / 88

c     modifs:
c     19 11 99 : élimination en sortie de nh1, nhe1, nhe2, lamb
c     20 08 01 : F95

c     Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k

c entree :
c     p : pression
c     t : température
c     xchim : composition chimique en grammes
c     deriv=.false. : évite le calcul de certaines dérivées (ineffectif)

c sortie :
c     ro : densité et dérivees
c     u : énergie interne et dérivées

c----------------------------------------------------------------------

      USE mod_donnees, ONLY : ah, amu, ahe4, eve, kbol, hpl, pi, me,
     1 echarg, aradia, granr, z0
      USE mod_kind
      USE mod_numerique, ONLY : gauss_band

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim   
      REAL (kind=dp), INTENT(in) :: pp, tt    
      REAL (kind=dp), INTENT(out) :: ro, drop, drot, drox, u,dup, dut, dux,
     1 delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
     2 gradad, dgradadp, dgradadt, dgradadx, alfa, beta, gamma1

      REAL (kind=dp), DIMENSION(2,2) :: aa
      REAL (kind=dp), DIMENSION(1,2) :: bb

      REAL (kind=dp), PARAMETER :: kih=13.595d0, kihe=24.580d0,
     1 kihe1=54.403d0, dx=1.d-8, unpdx=1.d0+dx,
     2 zisai=0.5d0, zih=8.d0, cf=15.d0, unsai=0.0625d0

      REAL (kind=dp), SAVE :: ahmu, ahemu, kihk, kihek, kihe1k, kihah,
     1 kihehe, kihe1he, a03, cte1, cte14, cte15,
     2 aradias3, dmunx, cte16, ah2ahe

      REAL (kind=dp), SAVE :: zizaiz, zihz, zai, tpre=0.d0, ropre=0.d0,
     1 nepre=1.d0, cte2

      REAL (kind=dp) ::  x, mum1, corr, stor, stor0, dstor, ro0, u0,
     1 nh1, nhe1, nhe2, y, cpp, cte19, cte20, t, p,
     2 cte21, cte22, cte23, cte24, cte25, cte26,
     3 drott, drotp, dutt ,dutp, drotx, dutx,
     4 dmux, dgx, dnex, dpx, drox0, dux0, fifi,
     6 ysahe, kihkt, kihekt, kt, kihe1kt,
     7 psi, dnh1ne, dnhene, muem1, bid, t32,
     8 deltp, deltpne, deltpro, munm1, ne00, zbar, zbar2, zbar3, dmu,
     5 xsah, dmuem1ne, dnero, dune,
     6 kihx, kihey, kihe1y,
     7 dgne, dgro, dgt, dnet, dpro, dpt, dpne, duro, drop0, drot0,
     8 dup0, dut0, ne0, nel, phi1, phi2, phi3, nemne, dphi1t,
     9 dphi2t, dphi3t, npart

      INTEGER, DIMENSION(2) :: indpc  
      INTEGER :: ntour, kder

      LOGICAL, SAVE :: init=.true.
      LOGICAL :: inversible

c----------------------------------------------------------------------

2000  FORMAT(8es10.3)

      IF(init)THEN    !initialisations
       init=.false.

c      Z proportion d'elements lourds doit avoir ete
c      initialise par appel (fictif) a opa
       ahmu=ah*amu ; ahemu=ahe4*amu ; kihk=kih*eve/kbol
       kihek=kihe*eve/kbol ; kihe1k=kihe1*eve/kbol
       kihah=kih*eve/ah/amu ; kihehe=kihe*eve/ahe4/amu
       kihe1he=(kihe+kihe1)*eve/ahe4/amu
       a03=(hpl/4.d0/pi/pi/me*hpl/echarg/echarg)**3
       cte1=3.d0/2.d0*kbol ; cte14=cf*a03*20.d0*kih*eve
       cte15=cf*a03 ; aradias3=aradia/3.d0
       dmunx=1.d0/ah-1.d0/ahe4 ; ah2ahe=(1.d0/ah-2.d0/ahe4)/amu
       cte16=20.d0*kih*eve ; zizaiz=zisai*z0/amu  ![ Zi /Ai ]h Z /mu
       zihz=zih*z0
       cte2=SQRT(hpl/me*hpl/kbol/2.d0/pi)**3/EXP(1.d0)**2/2.d0 !f/4=cte2*nel/T p.8
       WRITE(2,1) ; WRITE(*,1)
1      FORMAT(1x,'équation d''état de GONG2:',/,
     1 'équation de Saha évitant la recombinaison,',/,
     2 'pas de pression de radiation, pas de dégénérescence')
      ENDIF

c      PRINT*,nchim ; WRITE(*,2000)(xchim(kder),kder=1,nchim)

      t=tt ; p=pp ; x=xchim(1) ; y=1.d0-x-z0 ; xsah=x/ahmu ; ysahe=y/ahemu

      kihx=kihah*x        !pour 2.9
      kihey=kihehe*y ; kihe1y=kihe1he*y
c      PRINT*'kihx,kihey,kihe1y' ; WRITE(*,2000)kihx,kihey,kihe1y

      zai=unsai*z0/amu    ![ Zi / Ai ]h / mu
      munm1=xsah+ysahe+zai            !2.7/amu
      ne00=xsah+2.d0*ysahe+zizaiz     !2.18/amu
      zbar=x+2.d0*y+zihz ; zbar2=zbar**2 ; zbar3=zbar*zbar2
      cte19=cte15/zbar3 ; cte20=cte16*zbar2 ; cte21=-cte14/zbar
      cte22=cte19/2.d0 ; cte26=cte20*cte19 ; cte23=cte26/zbar
      cte24=cte23*2.d0 ; cte25=cte21/2.d0

c      solution de l'equation de SAHA par la methode de Mihalas
c      Stellar Atmosphere chapitre 5 formules 5-17 et 5-20 et Newton Raphson
c      2 variables : ro et nel

c      PRINT*'p,t' ; WRITE(*,2000)p,t

      IF(ABS(t-tpre)/t > .2)THEN  !initialisation de ro et nel si T a
       IF(t >= 1.d5)THEN  !varie de plus de 20% depuis l'appel precedent
        ro=p/kbol/t/(ne00+munm1) ; nel=ro*ne00            !ionisation totale
       ELSE
        ro=p/kbol/t/munm1 ; nel=ro*zizaiz !ionisation partielle
       ENDIF
      ELSE
       ro=ropre ; nel=nepre
      ENDIF
c      PRINT*'ro,nel provisoires' ; WRITE(*,2000)ro,nel,zizaiz

19    kder=1  !indice pour dérivées 1(p,t), 2dt
20    kihkt=kihk/t        !ki / k /T
      kihekt=kihek/t ; kihe1kt=kihe1k/t ; kt=kbol*t
      t32=SQRT(t)**3 ; bid=cte19*(kt+cte20) ; ntour=0
c      PRINT*'kder / ro,nel',kder ; WRITE(*,2000)ro,nel
21    ntour=ntour+1
      IF(ntour > 30)THEN
       PRINT*,'pas de conv. dans saha pour etat_gong2/P,T,epsi,nel,Tp'
       WRITE(*,2000)p,t,corr,nel,tpre ; PRINT*,'x,y,_local'
       WRITE(*,2000)x,y,z0
      ENDIF
      
      IF(nel <= 0.d0)THEN
       PRINT*,'nel <= 0.d0, ntour=',ntour ; PRINT*,'nel,p,t,xchim(1), ro ini, nel ini'
       WRITE(*,2000)nel,p,t,xchim(1),p/kbol/t/munm1,ro*zizaiz ; STOP
      ENDIF

      psi=2.d0+LOG(cte2*nel/t32)      !2.19
      dmu=cte19*(kt+cte20)*nel/kt     !2.13
      phi1=psi+kihkt-dmu ; phi1=EXP(phi1)*2.d0    ! 1 / 2.22
      nh1=1.d0/(1.d0+phi1) ; phi2=psi+kihekt-dmu ; phi2=EXP(phi2)/2.d0    ! 1 / 2.23
      phi3=psi+kihe1kt-dmu ; phi3=EXP(phi3)*2.d0              ! 1 / 2.24
      nhe2=1.d0/(phi3*(phi2+1.d0)+1.d0) ; nhe1=phi3*nhe2
      dnh1ne=-nh1**2*phi1 ; dnhene=-nhe2**2*phi3*(phi2*(phi3+4.d0)+1.d0)
      muem1=nh1*xsah+(nhe1+2.d0*nhe2)*ysahe+zizaiz    !a 1/amu pres
      mum1=munm1+muem1 ; dmuem1ne=(xsah*dnh1ne+dnhene*ysahe)*(1.d0-dmu)/nel
      ne0=ne00*ro ; nemne=(ne0-nel)*(ne0+nel) ; deltp=bid*nemne/2.d0
      deltpne=-bid*nel ; deltpro= bid*ne00*ne0
      bb(1,1)=kt*ro*mum1+deltp-p  !2.4  1/amu dans mum1
      bb(1,2)=ro*muem1-nel        !2.8
      aa(1,1)=kt*ro*dmuem1ne+deltpne  !jacobien    d 1 /dne
      aa(1,2)=kt*mum1+deltpro         !d 1 /dro
      aa(2,1)=ro*dmuem1ne-1.d0        !d 2 /dne
      aa(2,2)=muem1                   !d 2 /dro

c      solution

      indpc=1 ; CALL gauss_band(aa,bb,indpc,2,2,2,1,inversible)
      IF(.not.inversible)THEN
       PRINT*,'dans etat_gong2, matrice non inversible, ARRET' ; STOP
      ENDIF

      corr=1.d0
27    IF(corr*bb(1,1) > 0.6d0*nel .or. corr*bb(1,2) > 0.6d0*ro)THEN
       corr=corr/2.d0 ; GOTO 27
      ENDIF
      nel=nel-bb(1,1)*corr
      ro=ro-bb(1,2)*corr
      corr=max(ABS(bb(1,1))/nel,ABS(bb(1,2)/ro))

c      PRINT*'ntour kder/ ro,nel,bb(1,1),bb(1,2),corr',ntour,kder
c      WRITE(*,2000)ro,nel,bb(1,1),bb(1,2),corr
      IF(corr > 1.d-13)GOTO 21

c      PRINT*  'nel,ne0,ne0-nel,ro,nh1,nhe1,nhe2'
c      WRITE(*,2000)nel,ne0,ne0-nel,ro,nh1,nhe1,nhe2

c      calcul de delta, cp, gradad

      dmux=3.d0*dmu/zbar-cte24*nel/kt ; dgne=dmuem1ne*ro  ; dgro=muem1
      dphi1t=-phi1*(cte21*nel/kt+1.5d0+kihkt)/t
      dphi2t=-phi2*(cte21*nel/kt+1.5d0+kihekt)/t
      dphi3t=-phi3*(cte21*nel/kt+1.5d0+kihe1kt)/t
      dgt=-ro*(nh1**2*dphi1t*xsah+nhe2**2*(dphi2t*phi3*(phi3+2.d0)+
     1 dphi3t*(2.d0*phi2+1.d0))*ysahe)
      dgx=ro*(nh1/ahmu-(nhe1+2.d0*nhe2)/ahemu+dmux*(nh1**2*phi1*xsah+
     1 nhe2**2*phi3*(phi2*(phi3+4.d0)+1.)*ysahe))
      dnero=dgro/(1.d0-dgne) ; dnet =dgt /(1.d0-dgne) ; dnex =dgx /(1.d0-dgne)

      npart=munm1*ro+nel  !nombre de particules libres
      dpro=kt*munm1+deltpro ; dpt=kbol*npart+cte22*nemne*kbol
      dpne=kt-bid*nel
      dpx=granr*ro*t*dmunx+bid*ne0*ah2ahe*ro+3.d0*deltp/zbar-nemne*cte23

      dpro=dpro+dpne*dnero ; dpt= dpt +dpne*dnet ; dpx= dpx +dpne*dnex

      drop=1.d0/dpro ; drot=-drop*dpt ; drox=-drop*dpx

      fifi=nh1**2*phi1*kihx+
     1 nhe2**2*(phi3*(phi3*phi2-1.d0)*kihey+phi3*(2.d0*phi2+1.d0)*kihe1y)
      u=(1.5d0*kt*npart-cte25*nemne)/ro+kihx*nh1+nhe1*kihey+nhe2*kihe1y
      duro=-(1.5d0*kt*nel+cte25*(ne0**2+nel**2))/ro/ro
      dut=cte1*npart/ro-(nh1**2*dphi1t*kihx+nhe2**2*((dphi2t*phi3**2-
     1 dphi3t)*kihey+(dphi2t*phi3+dphi3t*(phi2+1.d0))*kihe1y))
      dune=(1.5d0*kt+cte21*nel)/ro-(1.d0-dmu)/nel*fifi
      dux=1.5d0*granr*t*dmunx-cte25/zbar*nemne/ro+cte26*ah2ahe*ne0+
     1 nh1*kihah-nhe1*kihehe-nhe2*kihe1he+dmux*fifi

c      PRINT*'drox,dux,dune,dnex' ; WRITE(*,2000)drox,dux,dune,dnex

      duro=duro+dune*dnero ; dut= dut +dune*dnet ; dux= dux +dune*dnex

      dup=duro*drop ; dut=dut +duro*drot ; dux=dux +duro*drox
      
      GOTO(30,31),kder
30    tpre=t ; ropre=ro ; nepre=nel

      ro0=ro ; drop0=drop ; drot0=drot ; drox0=drox
      u0=u ; dup0=dup ; dut0=dut ; dux0=dux
      stor0=t ; stor=stor0*unpdx ; dstor=stor-stor0 ; t=stor ; kder=2
      GOTO 20

31    t=stor0     !dérivée/t
      drott=(drot-drot0)/dstor ; drotp=(drop-drop0)/dstor

      drotx=(drox-drox0)/dstor ; dutt =(dut -dut0 )/dstor
      dutp =(dup -dup0 )/dstor ; dutx =(dux -dux0 )/dstor

c      PRINT*'d.,d.-d.0'
c      WRITE(*,2000)drot,drot-drot0,drop,drop-drop0,dut,dut-dut0,dup,dup-dup0

      ro=ro0 ; drop=drop0 ; drot=drot0 ; drox=drox0
      u=u0 ; dup=dup0 ; dut=dut0 ; dux=dux0 ; 

      delta=-t/ro*drot ; deltap=delta*(-drop/ro+drotp/drot)
      deltat=delta*(-drot/ro+drott/drot+1.d0/t)
      deltax=delta*(-drox/ro+drotx/drot)

      cpp=p/ro/t*delta ; cp=dut+cpp   
      dcpp=dutp+cpp*(-drop/ro+deltap/delta+1.d0/p)
      dcpt=dutt+cpp*(-drot/ro+deltat/delta-1.d0/t)
      dcpx=dutx+cpp*(-drox/ro+deltax/delta)

      gradad=p/ro*delta/cp/t      !gradient adiabatique
      dgradadp=gradad*(-drop/ro+deltap/delta-dcpp/cp+1.d0/p)
      dgradadt=gradad*(-drot/ro+deltat/delta-dcpt/cp-1.d0/t)
      dgradadx=gradad*(-drox/ro+deltax/delta-dcpx/cp)
      
      alfa=p/ro*drop ; beta=1.d0-aradias3*t**4/p
      gamma1=1.d0/(alfa-delta*gradad)

      END SUBROUTINE etat_gong2
      
