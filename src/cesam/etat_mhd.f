
c*************************************************************

      SUBROUTINE etat_mhd(pp,tt,xchim,deriv,
     1 ro,drop,drot,drox,u,dup,dut,dux,
     2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c     routine public du module mod_etat

c     équation d'état MHD de W. Dappen

c     interface MHD/CESAM
c     auteurs G.Berthomieu, P.Morel
c     version 25 02 92

c    modifs
c    06 06 92 W. DAPPEN nouvelles tables, Pression de radiation
c    19 11 99 : suppression de nh1, nhe1, nhe2, lamb

c     adaptation CESAM2k P.Morel

c     les noms des fichiers d'equation d'etat peuvent etre lus
c     dans le fichier *don puis eventuellement modifies dans cette routine

c entree :
c     p : pression
c     t : temperature
c     xchim : composition chimique
c     deriv=.FALSE. : evite le calcul  de certaines derivees

c sortie :
c     ro : densite et derivees
c     u : energie interne et derivees
c     delta, cp, gradad : notations thermodynamiques evidentes, et derivees
c     alpha, beta, gamma1 : notations thermodynamiques evidentes

c----------------------------------------------------------------

      USE mod_donnees, ONLY : aradia, f_eos, ln10, nom_chemin
      USE mod_kind

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
      REAL (kind=dp), INTENT(in) :: pp, tt
      LOGICAL, INTENT(in) :: deriv 
      REAL (kind=dp), INTENT(out) :: ro,drop,drot,drox,u,dup,dut,dux,
     1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1

      REAL (kind=dp), PARAMETER :: dx=1.d-5, unpdx=1.d0+dx
      REAL (kind=dp), SAVE :: aradias3

      REAL (kind=dp) :: prad,pcrit,duro,pg,rl0,pg10,t10,xh,gradad0,
     1 delta0,cp0,pg10s,t10s,stor,stor0,dstor,p,t

      INTEGER, DIMENSION(8) :: ir
      INTEGER :: ier, i, ntab, ir0

      LOGICAL, SAVE :: init=.TRUE.

      INTEGER, PARAMETER :: ivarx=25
      REAL (kind=dp), DIMENSION(ivarx) :: varmhd
      COMMON/mhdout/varmhd

      CHARACTER (len=80), DIMENSION(8) :: eostb

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccc CALLs mhd tables ccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

2000  FORMAT(8es10.3)

      p=pp ; t=tt

      IF(init)THEN ! calculs effectues au premier passage
       init=.FALSE. ; aradias3=aradia/3.d0

c      identification des fichiers de donnees

       f_eos(1)='mhd1.bin' ; f_eos(2)='mhd2.bin'
       f_eos(3)='mhd3.bin' ; f_eos(4)='mhd4.bin'
       f_eos(5)='mhd5.bin' ; f_eos(6)='mhd6.bin'
       f_eos(7)='mhd7.bin' ; f_eos(8)='mhd8.bin'
       eostb=TRIM(nom_chemin)//f_eos

c ...... set mhd equation of state ................

c======================================================================


c uses 3 ZAMS-type tables AND optionally 5 center-type tables.
c======================================================================

c      WRITE(*,*)'Equation d''etat MHD, entrer le nombre de tables'
c      WRITE(*,*) 'number of tables (3 [ZAMS] or 8 [ZAMS + center])?'
c      READ(*,*) ntab
       ntab=8
       IF(ntab /= 3 .AND. ntab /= 8)THEN
        WRITE(*,*)'illegal number of tables. error.' ; STOP
       ENDIF
       WRITE(*,*)'Utilisation de l''equation d''etat MHD'
       WRITE(2,*)'Utilisation de l''equation d''etat MHD'
       IF(ntab == 3)THEN
        WRITE(*,*)'utilisation des 3 tables de ZAMS'
        WRITE(2,*)'utilisation des 3 tables de ZAMS'
       ELSE
        WRITE(*,*)'utilisation de 3 tables de ZAMS et de 5 tables de centre'
        WRITE(2,*)'utilisation de 3 tables de ZAMS et de 5 tables de centre'
       ENDIF

       ir(1:8)=0 ; ir0=10
       DO i=1,ntab
        ir(i) = ir0 + i
        OPEN(unit=ir(i),status='old',file=eostb(i),form='unformatted')
       ENDDO
       CALL mhdst(ir(1),ir(2),ir(3),ir(4),ir(5),ir(6),ir(7),ir(8))
      ENDIF        !fin des operations du premier passage

c       type *,'p,t',p,t
c     Modif W. Dappen ----- pour eviter des appels de p<p_rad
c     pg=p-7.564d-15*t**4/3.
c     ---------------- la valeur de la constante crmhd est identique a
c     ---------------- celle employe dans le code MHD qui a produit les
c     ---------------- tables. Le 1.d-6 est assez arbitraire.

      ier=0 ; prad=aradias3*t**4 ; pcrit=prad+1.d-6*prad ; pg=p-prad
      pg=MAX(pg,pcrit) ; pg10=LOG10(pg) ; t10=LOG10(t) ; xh=xchim(1)
      CALL mhdpx(pg10,t10,xh,rl0,ier)
      IF(ier == 1 .OR. varmhd(4) == 0.d0)THEN
       PRINT*,'hors des tables MHD, appel a EFF'
       CALL etat_eff(p,t,xchim,
     1 ro,drop,drot,drox,u,dup,dut,dux,
     2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
       RETURN
      ENDIF

      ro=10.d0**rl0 ; drop=ro/varmhd(4) ; drot=-drop*varmhd(5)/t
      drox=-drop*varmhd(21)*ln10 ; drop=drop/p
      u=10.d0**varmhd(3)/ro ; duro=varmhd(6)-1.d0
      dut=varmhd(7) ; dux=varmhd(22)*ln10 ; dup=duro*drop*u/ro
      dut=(dut/t+duro/ro*drot)*u ; dux=(dux+duro/ro*drox)*u
      cp=10.**varmhd(9) ; gradad=varmhd(8) ; delta=varmhd(5)/varmhd(4)
      alfa=p/ro*drop ; beta=1.d0-prad/p ; gamma1=1./(alfa-delta*gradad)

c     nh1 = varmhd(14)    !taux d'ionisation H, He4
c     nhe1 = varmhd(15)
c     nhe2 = varmhd(16)
c     lamb = varmhd(18)   !degenerescence

      IF(deriv)THEN
       pg10s=pg10 ; t10s=t10

c      derivees / P

       stor0=p ; stor=stor0*unpdx
       IF(stor < dx)stor=dx
       dstor=stor-stor0 ; p=stor ; pcrit=prad+1.d-6*prad
       pg=p-prad ; pg=MAX(pg,pcrit) ; pg10=LOG10(pg)
       CALL mhdpx(pg10,t10,xh,rl0,ier) ; p=stor0 ; pg10=pg10s
       gradad0=varmhd(8) ; delta0=varmhd(5)/varmhd(4)
       cp0=10.**varmhd(9) ; dgradadp=(gradad0-gradad)/dstor
       deltap=(delta0-delta)/dstor ; dcpp=(cp0-cp)/dstor

c      derivees / T

       stor0=t ; stor=stor0*unpdx
       IF(stor < dx)stor=dx
       dstor=stor-stor0 ; t=stor ; prad=aradias3*t**4
       pcrit=prad+1.d-6*prad ; pg = p - prad ; pg=MAX(pg,pcrit)
       pg10=LOG10(pg) ; t10=LOG10(t) ; CALL mhdpx(pg10,t10,xh,rl0,ier)
       t=stor0 ; pg10=pg10s ; t10=t10s
       gradad0=varmhd(8) ; delta0=varmhd(5)/varmhd(4)
       cp0=10.**varmhd(9) ; dgradadt=(gradad0-gradad)/dstor
       deltat=(delta0-delta)/dstor ; dcpt=(cp0-cp)/dstor

c      derivees / X

       stor0=xh ; stor=stor0*unpdx
       IF(stor < dx)stor=dx
       dstor=stor-stor0 ; xh=stor ; CALL mhdpx(pg10,t10,xh,rl0,ier)
       gradad0=varmhd(8) ; delta0=varmhd(5)/varmhd(4)
       cp0=10.**varmhd(9) ; dgradadx=(gradad0-gradad)/dstor
       deltax=(delta0-delta)/dstor ; dcpx=(cp0-cp)/dstor
      ENDIF

      RETURN
      
      END SUBROUTINE etat_mhd

c****************************************************************************

      SUBROUTINE mhdpx(pgl,tl,xc,rl,ierr)

c computes derivative with respect to X also in ZAMS-table situation.
c changed SUBROUTINE is mhdpx1, CALLing quintd instead of quint.

c=====================================================================
c=====       MHD interpolation package      ==========================
c=====          pressure as argument        ==========================
c=====           interpolation in X         ==========================
c=====                Z fixed               ==========================
c=====================================================================

c     input   pgl = log10 (gas pressure [dyn/cm2])
c             tl  = log10 (temperature [K])
c             xc  = hydrogen abundance (by mass)

c     output  rl  = log10 (density [g/cm3])        : ARGUMENT
c             other (see separate instructions)    : COMMON/MHDOUT/...

c     error   ierr = 1 signals pgl,tl outside the domain of tables
c             (otherwise ierr = 0).

c           s/r mhdstx must be called in main.


	USE mod_kind
	
      IMPLICIT NONE

      INTEGER, PARAMETER :: ivarc=25

      REAL (kind=dp), DIMENSION(ivarc) :: varmhd
      COMMON/mhdout/varmhd

      INTEGER irange,ierr

      REAL (kind=dp) :: pgl,tl,xc,rl

      CALL mhdpx1(pgl,tl,xc,irange)

      ierr=0
      IF(irange == 0)THEN
       WRITE(*,*) 'out of table range. RETURN'
       ierr=1
       RETURN
      ENDIF

      rl=varmhd(1)

      END SUBROUTINE mhdpx

C****************************************************************************

      SUBROUTINE mhdpx1(pgl,tl,x,irange)
      
c     mhdst must be CALLed in main.

c     interpolation in tables with different x AND fixed z

	USE mod_kind

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ivarx=25,ndimt=8
      
      INTEGER ix(4),irange,i,itbl,iv,ixmin,ii,id,np,l,inter

      REAL (kind=dp) :: varout(ndimt,ivarx),xc(ndimt),dvdx(ivarx),vari1(ivarx,4),xx(4),
     1 pgl,tl,x,x0,h,y0,y1,y2,povere,dledlp

c     quantities for interpolation in x
c     tl < tlim1:         lower part of zams tables
c     tlim1 < tl < tlim2: upper part of zams tables
c     tl > tlim2:         centre tables
c     tmini,tmaxi:        temperature interval covered by the tables

      REAL (kind=dp) :: tlim1,tlim2,tmini,tmaxi
      COMMON/tttt/tlim1,tlim2,tmini,tmaxi

      REAL (kind=dp) :: varfin(ivarx) !output
      COMMON/mhdout/varfin

1000  FORMAT('results from mhdpx2, itbl,x = ',i6,es15.7/)
1001  FORMAT(12(/1x,1p5e15.6))
5001  FORMAT(1x,'******* warning: extrapolation in X (quint) ',
     1 'pgl,tl,x= ',3f12.6)
5011  FORMAT(1x,'******* warning: extrapolation in X (lir) ',
     1 'pgl,tl,x = ',3f12.6)
9001  FORMAT(' error in mhdpx1. tl out of range. tl,tmini,tmaxi=',
     1 1x,3es13.6)

c             READ from appropriate tables
c             AND fill arrays varout AND xc

c------------ can we do it?

      irange=1
      IF(tl < tmini .OR. tl > tmaxi)THEN
       PRINT 9001,tl,tmini,tmaxi ; irange=0
       RETURN
      ENDIF

      IF(tl < tlim1)THEN  !lower zams tables
       DO i=1,3
        itbl=-i ; CALL mhdpx2(pgl,tl,itbl,varout,xc,ndimt,irange)
c       PRINT 1000,itbl,xc(i) ; PRINT 1001,(varout(i,jj),jj=1,ivarx)
       ENDDO

      ELSEIF(tl < tlim2)THEN  !upper zams tables
       DO i=1,3
       itbl=i ; CALL mhdpx2(pgl,tl,itbl,varout,xc,ndimt,irange)
c       PRINT 1000,itbl,xc(i) ; PRINT 1001,(varout(i,jj),jj=1,ivarx)
       ENDDO

      ELSE    !centre tables
       DO i=1,5   !offset (+3) to access center tables
        itbl=i+3 ; CALL mhdpx2(pgl,tl,itbl,varout,xc,ndimt,irange)
c       PRINT 1000,itbl,xc(itbl) ; PRINT 1001,(varout(itbl,jj),jj=1,ivarx)
       ENDDO
      ENDIF

c===================== interpolation in x ========================

      IF(tl <= tlim2)THEN  !quadratic newton (equidistant xc's)
       x0=xc(1) ; h=xc(2)-xc(1)
       IF(ABS(xc(3)-xc(2)-h) > 1.e-4)THEN
        WRITE(*,*) 'non-equidistant zams tables. error stop'
        WRITE(*,*) 'xc(1-3)= ',xc(1),xc(2),xc(3) ; STOP
       ENDIF

       IF((xc(3) > xc(1) .AND. (x>xc(3) .OR. x<xc(1))) .OR.
     1  ( xc(3) < xc(1) .AND. (x > xc(1) .OR. x < xc(3))))
     2  WRITE(*,5001)pgl,tl,x
       DO iv=1,ivarx
        y0 = varout(1,iv) ; y1 = varout(2,iv) ; y2 = varout(3,iv)
        CALL quintd(x,x0,h,y0,y1,y2,varfin(iv),dvdx(iv))
       ENDDO
      
       varfin(21)=-dvdx(1)*varfin(4) ; povere=10.d0**(varfin(2)-varfin(3))
       dledlp=(povere*(1.d0-varfin(5))+1.d0)/varfin(4)
       varfin(22)=dvdx(3)+varfin(21)*dledlp

        varfin(23)=dvdx(8)-dvdx(1)*varfin(10)
       varfin(24)=dvdx(9)-dvdx(1)*varfin(12) ; varfin(25)=-5555533333.d0
      ELSE    !cubic lagrangian (arbitrarily spaced xc's)

c      use 4 highest xc(i)'s, if possible.
c      distinguish case where the xc(i), i=4..8,
c      increase or decrease. 

       IF(xc(4) < xc(8))THEN

        ixmin = 5
        IF(x < xc(5))ixmin=4
       ELSE
        ixmin = 4
        IF(x < xc(7)) ixmin=5
       ENDIF

       DO i=1,4
        ix(i)=ixmin+i-1
       ENDDO

       DO i =1,4
        xx(i) = xc(ix(i))
c       WRITE(*,*)' xx(i) for lir: i,xx(i),x ',i,xx(i),x
        DO iv=1,ivarx
         vari1(iv,i) = varout(ix(i),iv)
        ENDDO
       ENDDO
       ii=ivarx ; id=ivarx ; np=4 ; l=1 ; inter=1
       CALL lir(x,xx,varfin,vari1,ii,id,np,l,inter)
       IF(inter == 0)WRITE(*,5011) pgl,tl,x
      ENDIF
C
      END SUBROUTINE mhdpx1
C
C****************************************************************************
C
      SUBROUTINE mhdpx2(pgl,tl,itbl,varout,xcomp,ndimt,irange)
      
	USE mod_kind
	
      IMPLICIT NONE

      INTEGER, PARAMETER :: ivarc=20, ivarx=25, nchem0=6,
     1 nt1m=39, nt2m=63, ntxm=22, nr1m=139, nr2m=61, nrxm=61

      INTEGER :: itbl,ndimt,irange

      REAL (kind=dp) :: pgl,tl,varout(ndimt,ivarx),xcomp(ndimt),varo(ivarx),
     1 vvar1(ivarx,4),vvar2(ivarx,4),vy(ivarx)

c.    zams tables (labelled by a,b,c)

      INTEGER :: nt1,nr1
      REAL (kind=dp) :: tdvr1a(nt1m,nr1m,ivarc),tlog1(nt1m),drh1
      COMMON/tab1a/tdvr1a,tlog1,nt1,nr1,drh1

      INTEGER :: nt2,nr2
      REAL (kind=dp) :: tdvr2a(nt2m,nr2m,ivarc),tlog2(nt2m),drh2
      COMMON/tab2a/tdvr2a,tlog2,nt2,nr2,drh2

      REAL (kind=dp) :: tdvr1b(nt1m,nr1m,ivarc)
      COMMON/tab1b/tdvr1b

      REAL (kind=dp) :: tdvr2b(nt2m,nr2m,ivarc)
      COMMON/tab2b/tdvr2b

      REAL (kind=dp) :: tdvr1c(nt1m,nr1m,ivarc)
      COMMON/tab1c/tdvr1c

      REAL (kind=dp) :: tdvr2c(nt2m,nr2m,ivarc)
      COMMON/tab2c/tdvr2c

      REAL (kind=dp) :: atwta(nchem0),abuna(nchem0),abfrca(nchem0),gasma
      COMMON/chea/atwta,abuna,abfrca,gasma

      REAL (kind=dp) :: atwtb(nchem0),abunb(nchem0),abfrcb(nchem0),gasmb
      COMMON/cheb/atwtb,abunb,abfrcb,gasmb

      REAL (kind=dp) :: atwtc(nchem0),abunc(nchem0),abfrcc(nchem0),gasmc
      COMMON/chec/atwtc,abunc,abfrcc,gasmc

c     centre tables (labelled by 1,2,3,4,5)

      INTEGER :: ntx,nrx
      REAL (kind=dp) :: tdvrx1(ntxm,nrxm,ivarx),tlogx(ntxm),drhx
      COMMON/tabx1/tdvrx1,tlogx,ntx,nrx,drhx

      REAL (kind=dp) :: tdvrx2(ntxm,nrxm,ivarx)
      COMMON/tabx2/tdvrx2

      REAL (kind=dp) :: tdvrx3(ntxm,nrxm,ivarx)
      COMMON/tabx3/tdvrx3

      REAL (kind=dp) :: tdvrx4(ntxm,nrxm,ivarx)
      COMMON/tabx4/tdvrx4
      
      REAL (kind=dp) :: tdvrx5(ntxm,nrxm,ivarx)
      COMMON/tabx5/tdvrx5

      REAL (kind=dp) :: atwt1(nchem0),abun1(nchem0),abfrc1(nchem0),gasm1
      COMMON/che1/atwt1,abun1,abfrc1,gasm1

      REAL (kind=dp) :: atwt2(nchem0),abun2(nchem0),abfrc2(nchem0),gasm2
      COMMON/che2/atwt2,abun2,abfrc2,gasm2

      REAL (kind=dp) ::   atwt3(nchem0),abun3(nchem0),abfrc3(nchem0),gasm3
      COMMON/che3/atwt3,abun3,abfrc3,gasm3

      REAL (kind=dp) ::   atwt4(nchem0),abun4(nchem0),abfrc4(nchem0),gasm4
      COMMON/che4/atwt4,abun4,abfrc4,gasm4

      REAL (kind=dp) ::  atwt5(nchem0),abun5(nchem0),abfrc5(nchem0),gasm5
      COMMON/che5/atwt5,abun5,abfrc5,gasm5

c=================================================================

c     nomenclature for accessing the tables
c     =====================================

c     zams tables

c     itbl = -1   : tdvr1a
c     itbl =  1   : tdvr2a
c     itbl = -2   : tdvr1b
c     itbl =  2   : tdvr2b
c     itbl = -3   : tdvr1c
c     itbl =  3   : tdvr2c

c     contral (variable x) tables

c     itbl =  4   : tdvrx1
c     itbl =  5   : tdvrx2
c     itbl =  6   : tdvrx3
c     itbl =  7   : tdvrx4
c     itbl =  8   : tdvrx5

c     note: the same output variables varout(itbl,.) is used
c           for itbl=1 AND -1, 2 AND -2, 3 AND -3, respectively.

c=================================================================

      irange=1
      varo(:ivarx)=0.d0       !!CALL zero(varo,ivarx) PM

c------------------------------------------------------------------
c-------------- main selection of tables --------------------------
c------------------------------------------------------------------

c               zams tables
c               ***********

      IF(itbl == -1)THEN
       CALL intpt(pgl,tl,tdvr1a,nt1m,nr1m,ivarc,tlog1,nt1,nr1,
     1  vvar1,vvar2,vy,varo,irange)
       varout(1,1:ivarc)=varo(1:ivarc) ; xcomp(1)=abfrca(1)
      ELSEIF(itbl == 1)THEN
       CALL intpt(pgl,tl,tdvr2a,nt2m,nr2m,ivarc,tlog2,nt2,nr2,
     .  vvar1,vvar2,vy,varo,irange)
       varout(1,1:ivarc) = varo(1:ivarc) ; xcomp(1)    = abfrca(1)
      ELSEIF(itbl==-2)THEN
       CALL intpt(pgl,tl,tdvr1b,nt1m,nr1m,ivarc,tlog1,nt1,nr1,
     1  vvar1,vvar2,vy,varo,irange)
       varout(2,1:ivarc)=varo(1:ivarc) ; xcomp(2)=abfrcb(1)
      ELSEIF(itbl == 2)THEN
       CALL intpt(pgl,tl,tdvr2b,nt2m,nr2m,ivarc,tlog2,nt2,nr2,
     1  vvar1,vvar2,vy,varo,irange)
       varout(2,1:ivarc)=varo(1:ivarc) ; xcomp(2)=abfrcb(1)
      ELSEIF(itbl == -3)THEN
       CALL intpt(pgl,tl,tdvr1c,nt1m,nr1m,ivarc,tlog1,nt1,nr1,
     1  vvar1,vvar2,vy,varo,irange)
       varout(3,1:ivarc) = varo(1:ivarc) ; xcomp(3)    = abfrcc(1)

      ELSEIF(itbl== 3)THEN
       CALL intpt (pgl,tl,tdvr2c,nt2m,nr2m,ivarc,tlog2,nt2,nr2,
     1  vvar1,vvar2,vy,varo,irange)
       varout(3,1:ivarc)=varo(1:ivarc) ; xcomp(3)=abfrcc(1)

c               center tables
c               *************

      ELSEIF(itbl == 4)THEN
       CALL intpt(pgl,tl,tdvrx1,ntxm,nrxm,ivarx,tlogx,ntx,nrx,
     1  vvar1,vvar2,vy,varo,irange)
       varout(4,1:ivarx)=varo(1:ivarx) ; xcomp(4)=abfrc1(1)
      ELSEIF(itbl== 5)THEN
       CALL intpt(pgl,tl,tdvrx2,ntxm,nrxm,ivarx,tlogx,ntx,nrx,
     1  vvar1,vvar2,vy,varo,irange)
       varout(5,1:ivarx)=varo(1:ivarx) ; xcomp(5)=abfrc2(1)
      ELSEIF(itbl == 6)THEN
       CALL intpt(pgl,tl,tdvrx3,ntxm,nrxm,ivarx,tlogx,ntx,nrx,
     1  vvar1,vvar2,vy,varo,irange)
       varout(6,1:ivarx)=varo(1:ivarx) ; xcomp(6)=abfrc3(1)
      ELSEIF(itbl == 7)THEN
       CALL intpt (pgl,tl,tdvrx4,ntxm,nrxm,ivarx,tlogx,ntx,nrx,
     1  vvar1,vvar2,vy,varo,irange)
       varout(7,1:ivarx)=varo(1:ivarx) ; xcomp(7)    = abfrc4(1)
      ELSEIF(itbl == 8)THEN
       CALL intpt(pgl,tl,tdvrx5,ntxm,nrxm,ivarx,tlogx,ntx,nrx,
     1  vvar1,vvar2,vy,varo,irange)
       varout(8,1:ivarx)=varo(1:ivarx) ; xcomp(8)=abfrc5(1)
      ENDIF

c-------------- end selection of tables --------------------------
C
      END SUBROUTINE mhdpx2
C
C****************************************************************************
C
      SUBROUTINE intpt(pgl,tl,tdvar,ntm,nrm,ivar,tlog,nt,nr,
     1 var1,var2,y,varout,irange)

	USE mod_kind

      IMPLICIT NONE

      INTEGER :: irm(4),ir(4,4),it(4),ntm,nrm,ivar,nt,nr,irange,n,i,indt,itmax,
     1  m,irmax,j,ilir  ,indr,iv,ii,id,np,l,inter

      REAL (kind=dp) :: tdvar(ntm,nrm,ivar),tlog(ntm),var1(ivar,4),pgl,tl,
     1   var2(ivar,4),varout(ivar),y(ivar),xx(4),pmini,pmaxi

      DO n=1,nt
       IF(tlog(n)>=tl)GOTO 101
       it(1)=n
      ENDDO

101   IF(it(1)>=2)it(1)=it(1)-1
      itmax=nt-3
      IF(it(1)>itmax)it(1)=itmax
      DO i=2,4
       it(i)=it(1)+i-1

      ENDDO

      DO i=1,4
       irm(i)=1 ; indt=it(i) ; pmini=tdvar(indt,1,20) ; pmaxi=tdvar(indt,nr,20)
       IF(pgl>pmaxi)THEN
        PRINT 9002,tl,tlog(indt),pgl,pmaxi ; irange=0 ; RETURN
       ENDIF
       IF(pgl<pmini)PRINT 9003,tl,pgl,pmini
       DO m=1,nr
        IF(tdvar(indt,m,20)>=pgl) GOTO 201
        irm(i)=m
       ENDDO
201    IF(irm(i)>=2)irm(i)=irm(i)-1
       irmax=nr-3
       IF(irm(i)>irmax)irm(i)=irmax
      ENDDO
      DO i=1,4
       DO j=1,4
        ir(j,i)=irm(i)+j-1
       ENDDO

      ENDDO

      DO ilir=1,4
       indt=it(ilir)
       DO i=1,4
        indr=ir(i,ilir) ; xx(i)=tdvar(indt,indr,20)
       ENDDO
       DO i=1,4
        indr=ir(i,ilir)
        DO iv=1,ivar
         var1(iv,i)=tdvar(indt,indr,iv)
        ENDDO

       ENDDO
       ii=ivar ; id=ivar ; np=4 ; l=1 ; inter=1
       CALL lir(pgl,xx,y,var1,ii,id,np,l,inter)
       DO iv=1,ivar
        var2(iv,ilir)=y(iv)

       ENDDO

      ENDDO
      DO i=1,4
       indt=it(i) ; xx(i)=tlog(indt)

      ENDDO
      ii=ivar ; id=ivar ; np=4 ; l=1 ; inter=1
      CALL lir(tl,xx,varout,var2,ii,id,np,l,inter)

9002  FORMAT(' error in intpt. pgl out of range.',
     1 ' tl,tlog(it(i)),pgl,pglmax: ',/,4es13.5)
9003  FORMAT(' ********** warning: extrapolation in s/r intpt.',
     1 ' tl,pgl,pglmin = ',3f10.6)

      END SUBROUTINE intpt
C
C****************************************************************************
C
      SUBROUTINE mhdst(ir1,ir2,ir3,ir4,ir5,ir6,ir7,ir8)

	USE mod_kind

      IMPLICIT NONE

      INTEGER, PARAMETER :: ivarc=20, ivarx=25, nchem0=6,
     1 nt1m=39, nt2m=63, ntxm=22, nr1m=139, nr2m=61, nrxm=61

      INTEGER :: ir1,ir2,ir3,ir4,ir5,ir6,ir7,ir8,nchem,ntdum,idx,nrdum

      REAL (kind=dp) :: tddum(nt1m,nr1m,ivarc),tldum(nt1m),tlogd(nt2m),tlogu(nt2m),
     1 tddif0(nt2m,nr2m,ivarc),tddifd(nt2m,nr2m,ivarc),
     2 tddifu(nt2m,nr2m,ivarc),atwd(nchem0),abud(nchem0),
     3 abfrcd(nchem0),atwu(nchem0),abuu(nchem0),abfrcu(nchem0),drdum      

c....... tl<tlim1:       lower part of zams tables
c....... tlim1<tl<tlim2: upper part of zams tables
c....... tl>tlim2:       variable x tables

      REAL (kind=dp) :: tlim1,tlim2,tmini,tmaxi
      COMMON/tttt/tlim1,tlim2,tmini,tmaxi

c..zams tables (labelled by a,b,c)

      INTEGER :: nt1,nr1
      REAL (kind=dp) :: tdvr1a(nt1m,nr1m,ivarc),tlog1(nt1m),drh1
      COMMON/tab1a/tdvr1a,tlog1,nt1,nr1,drh1

      INTEGER :: nt2,nr2
      REAL (kind=dp) :: tdvr2a(nt2m,nr2m,ivarc),tlog2(nt2m),drh2
      COMMON/tab2a/tdvr2a,tlog2,nt2,nr2,drh2

      REAL (kind=dp) :: tdvr1b(nt1m,nr1m,ivarc)
      COMMON/tab1b/tdvr1b

      REAL (kind=dp) :: tdvr2b(nt2m,nr2m,ivarc)
      COMMON/tab2b/tdvr2b

      REAL (kind=dp) :: tdvr1c(nt1m,nr1m,ivarc)
      COMMON/tab1c/tdvr1c

      REAL (kind=dp) :: tdvr2c(nt2m,nr2m,ivarc)
      COMMON/tab2c/tdvr2c

      REAL (kind=dp) :: atwta(nchem0),abuna(nchem0),abfrca(nchem0),gasma
      COMMON/chea/atwta,abuna,abfrca,gasma

      REAL (kind=dp) :: atwtb(nchem0),abunb(nchem0),abfrcb(nchem0),gasmb    
      COMMON/cheb/atwtb,abunb,abfrcb,gasmb

      REAL (kind=dp) :: atwtc(nchem0),abunc(nchem0),abfrcc(nchem0),gasmc
      COMMON/chec/atwtc,abunc,abfrcc,gasmc

c     centre tables (labelled by 1,2,3,4,5)

      INTEGER :: ntx,nrx
      REAL (kind=dp) :: tdvrx1(ntxm,nrxm,ivarx),tlogx(ntxm),drhx
      COMMON/tabx1/tdvrx1,tlogx,ntx,nrx,drhx

      REAL (kind=dp) :: tdvrx2(ntxm,nrxm,ivarx)
      COMMON/tabx2/tdvrx2

      REAL (kind=dp) :: tdvrx3(ntxm,nrxm,ivarx)
      COMMON/tabx3/tdvrx3

      REAL (kind=dp) :: tdvrx4(ntxm,nrxm,ivarx)
      COMMON/tabx4/tdvrx4

      REAL (kind=dp) :: tdvrx5(ntxm,nrxm,ivarx) 
      COMMON/tabx5/tdvrx5

      REAL (kind=dp) :: atwt1(nchem0),abun1(nchem0),abfrc1(nchem0),gasm1
      COMMON/che1/atwt1,abun1,abfrc1,gasm1

      REAL (kind=dp) :: atwt2(nchem0),abun2(nchem0),abfrc2(nchem0),gasm2
      COMMON/che2/atwt2,abun2,abfrc2,gasm2

      REAL (kind=dp) :: atwt3(nchem0),abun3(nchem0),abfrc3(nchem0),gasm3
      COMMON/che3/atwt3,abun3,abfrc3,gasm3

      REAL (kind=dp) :: atwt4(nchem0),abun4(nchem0),abfrc4(nchem0),gasm4
      COMMON/che4/atwt4,abun4,abfrc4,gasm4
      
      REAL (kind=dp) :: atwt5(nchem0),abun5(nchem0),abfrc5(nchem0),gasm5 
      COMMON/che5/atwt5,abun5,abfrc5,gasm5

c------------ define, with unused statements, storage for variables
c------------ that would otherwise only appear as formal PARAMETERs
c------------ in this SUBROUTINE. vicious bugs can be the result if
c------------ this storage were not provided (remember: fortran is not
c------------ a recursive language.)

      nchem=0 ; drdum=0.d0 ; ntdum=0 ; nrdum=0

c======================= READ zams tables ===========================

      IF(ir1 > 0)THEN
       idx=0
       CALL mhdst1(ir1,idx,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0,
     1  nt1,nr1,nt2,nr2,tlog1,tlog2,tdvr1a,tdvr2a,
     2  drh1,drh2,nchem,atwta,abuna,abfrca,gasma,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF

      IF(ir2 > 0)THEN
       idx=0
       CALL mhdst1(ir2,idx,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0,
     1  nt1,nr1,nt2,nr2,tlog1,tlog2,tdvr1b,tdvr2b,
     2  drh1,drh2,nchem,atwtb,abunb,abfrcb,gasmb,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF

      IF(ir3 > 0)THEN
       idx=0
       CALL mhdst1(ir3,idx,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0,
     1  nt1,nr1,nt2,nr2,tlog1,tlog2,tdvr1c,tdvr2c,
     2  drh1,drh2,nchem,atwtc,abunc,abfrcc,gasmc,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF

c======================= READ centre tables =====================

      IF(ir4 > 0)THEN
       idx=1
       CALL mhdst1(ir4,idx,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0,
     1  ntdum,nrdum,ntx,nrx,tldum,tlogx,tddum,tdvrx1,
     2  drdum,drhx,nchem,atwt1,abun1,abfrc1,gasm1,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4 abud,abuu,abfrcd,abfrcu)
      ENDIF

      IF(ir5 > 0)THEN
       idx=1
       CALL mhdst1(ir5,idx,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0,
     1  ntdum,nrdum,ntx,nrx,tldum,tlogx,tddum,tdvrx2,
     2  drdum,drhx,nchem,atwt2,abun2,abfrc2,gasm2,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF

      IF(ir6>0)THEN
       idx=1
       CALL mhdst1(ir6,idx,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0,
     1  ntdum,nrdum,ntx,nrx,tldum,tlogx,tddum,tdvrx3,
     2  drdum,drhx,nchem,atwt3,abun3,abfrc3,gasm3,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF

      IF(ir7 > 0)THEN
       idx=1
       CALL mhdst1(ir7,idx,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0,
     1  ntdum,nrdum,ntx,nrx,tldum,tlogx,tddum,tdvrx4,
     2  drdum,drhx,nchem,atwt4,abun4,abfrc4,gasm4,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF
c
      IF(ir8 > 0)THEN
       idx=1
       CALL mhdst1(ir8,idx,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0,
     1  ntdum,nrdum,ntx,nrx,tldum,tlogx,tddum,tdvrx5,
     2  drdum,drhx,nchem,atwt5,abun5,abfrc5,gasm5,
     3  tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4  abud,abuu,abfrcd,abfrcu)
      ENDIF

c======================= temperature limits =====================

      tmini=tlog1(1) ; tlim1=tlog1(nt1)
      IF(ir4 <= 0)THEN
       tlim2=tlog2(nt2) ; tmaxi=tlim2
      ELSE
       tlim2=tlogx(1) ; tmaxi=tlogx(ntx)
      ENDIF

      WRITE(*,*)' before leaving mhdst: compositions of 8 tables '

c      PRINT 8002,(atwta(ic),abuna(ic),abfrca(ic),ic=1,nchem) ; PRINT 8003,gasma
c      PRINT 8002,(atwtb(ic),abunb(ic),abfrcb(ic),ic=1,nchem) ; PRINT 8003,gasmb
c      PRINT 8002,(atwtc(ic),abunc(ic),abfrcc(ic),ic=1,nchem) ; PRINT 8003,gasmc

c      IF(ir4 <= 0) GOTO 500

c      PRINT 8002,(atwt1(ic),abun1(ic),abfrc1(ic),ic=1,nchem) ; PRINT 8003,gasm1
c      PRINT 8002,(atwt2(ic),abun2(ic),abfrc2(ic),ic=1,nchem) ; PRINT 8003,gasm2
c      PRINT 8002,(atwt3(ic),abun3(ic),abfrc3(ic),ic=1,nchem) ; PRINT 8003,gasm3
c      PRINT 8002,(atwt4(ic),abun4(ic),abfrc4(ic),ic=1,nchem) ; PRINT 8003,gasm4
c      PRINT 8002,(atwt5(ic),abun5(ic),abfrc5(ic),ic=1,nchem) ; PRINT 8003,gasm5

 500  WRITE(*,*)' tlim1,tlim2 = ',tlim1,tlim2
8002  FORMAT('      at. weight     number ',
     1 'abundance  mass fraction',(/1x,1p3g16.7))
8003  FORMAT(/' mean molecular weight = ',f12.7//)

      END SUBROUTINE mhdst
C
C****************************************************************************
C
      SUBROUTINE mhdst1(ir,idx,nt1m,nr1m,ivar1,nt2m,nr2m,ivar2,nchem0,
     1 nt1,nr1,nt2,nr2,tlog1,tlog2,tdvar1,tdvar2,
     2 drh1,drh2,nchem,atwt,abun,abfrcs,gasmu,
     3 tlogd,tlogu,tddif0,tddifd,tddifu,atwd,atwu,
     4 abud,abuu,abfrcd,abfrcu)

	USE mod_kind

      IMPLICIT NONE
      
      INTEGER :: ir,idx,nt1m,nr1m,ivar1,nt2m,nr2m,ivar2,nchem0,
     1 nt1,nr1,nt2,nr2,nchem,ifiles,iidx,ivarr,idxr,irescr,iresco,
     2 ic,jt,n,m,iv    
      
      REAL (kind=dp) :: tlog1(nt1m),tdvar1(nt1m,nr1m,ivar1),ddx,
     1 tlog2(nt2m),tdvar2(nt2m,nr2m,ivar2),
     2 atwt(nchem0),abun(nchem0),abfrcs(nchem0),
     3 tlogd(nt2m),tlogu(nt2m),tddif0(nt2m,nr2m,ivar1),
     4 tddifd(nt2m,nr2m,ivar1),tddifu(nt2m,nr2m,ivar1),
     5 atwd(nchem0),abud(nchem0),abfrcd(nchem0),gasmu,
     6 atwu(nchem0),abuu(nchem0),abfrcu(nchem0),drh1,drh2,
     7 gasxd,gasxu,cxlim
      
c---------------------------------------------------------------

      IF(idx == 0)THEN
       ifiles=1
      ELSE
       ifiles=3
      ENDIF

      DO iidx=1,ifiles
c      READ(ir,98,end=1000) ivarr,idxr,irescr,ddx
       READ(ir,   end=1000) ivarr,idxr,irescr,ddx
       IF(ivar1<ivarr)THEN
        PRINT 9006,ivar1,ivarr ; STOP
       ENDIF
       IF(idx /= idxr)THEN
        PRINT 9007,idx,idxr ; STOP
       ENDIF
       IF(iidx == 1)THEN
        iresco= 0
       ELSEIF(iidx == 2)THEN
        iresco=-1
       ELSEIF(iidx == 3)THEN
        iresco= 1
       ENDIF

       IF(iresco/=irescr)THEN
        PRINT 9008,iresco,irescr ; STOP
       ENDIF

       WRITE(*,*) 'ivarr,idxr,iresco,ddx ',ivarr,idxr,iresco,ddx

       IF(iidx==1)CALL rabu(ir,nchem0,nchem,atwt,abun,abfrcs,gasmu)
       IF(iidx==2)CALL rabu(ir,nchem0,nchem,atwd,abud,abfrcd,gasxd)
       IF(iidx==3)CALL rabu(ir,nchem0,nchem,atwu,abuu,abfrcu,gasxu)

c      READ(ir,1001) nt1,nt2,drh1,drh2
       READ(ir)nt1,nt2,drh1,drh2

       IF(idx==1 .AND. nt1/=0)THEN
        PRINT 9010,idx,nt1,nt2 ; STOP
       ENDIF
       WRITE(*,*) 'nt1,nt2,drh1,drh2', nt1,nt2,drh1,drh2
       IF(nt1 > 0)CALL rtab(ir,nt1m,nr1m,ivar1,nt1,nr1,tlog1,tdvar1)
       IF(idx==1)THEN
        IF(iidx==1) CALL rtab(ir,nt2m,nr2m,ivar1,nt2,nr2,tlog2,tddif0)
        IF(iidx==2) CALL rtab(ir,nt2m,nr2m,ivar1,nt2,nr2,tlogd,tddifd)
        IF(iidx==3) CALL rtab(ir,nt2m,nr2m,ivar1,nt2,nr2,tlogu,tddifu)
       ELSE IF(idx==0)THEN
        CALL rtab(ir,nt2m,nr2m,ivar2,nt2,nr2,tlog2,tdvar2)
       ENDIF
      ENDDO
      IF(idx == 0)GOTO 450

c======================================================================
c........ if idx=1: check tables for correct composition construction
c........           AND perform numerical derivatives w.r.t. x
c======================================================================
c
      cxlim=0.05*ABS(ddx)
c
      IF(ABS(abfrcs(1)-abfrcd(1)-ddx)>cxlim .OR.
     1   ABS(abfrcs(1)-abfrcu(1)+ddx)>cxlim .OR.
     2   ABS(abfrcs(2)-abfrcd(2)+ddx)>cxlim .OR.
     3   ABS(abfrcs(2)-abfrcu(2)-ddx)>cxlim )GOTO 500
c
      DO ic=3,nchem
       IF( ABS(abfrcs(ic)-abfrcu(ic))>cxlim )GOTO 500
       IF( ABS(abfrcs(ic)-abfrcd(ic))>cxlim )GOTO 500
      ENDDO
c
      DO jt=1,nt2
       IF(tlog2(jt) /= tlogd(jt))GOTO 600
       IF(tlog2(jt) /= tlogu(jt))GOTO 600
      ENDDO
      PRINT 8001


c======================================================================
c....... numerical derivatives w.r.t. x
c======================================================================

      DO n=1,nt2
       DO m=1,nr2
        DO iv=1,ivar1
         tdvar2(n,m,iv)=tddif0(n,m,iv)

        ENDDO


c>>>>>>>>>> extended set of variables (tdvar2(n,m,ivar1+1...ivar2))
c>>>>>>>>>> for t-rho regions with inhomogeneous composition.
c>>>>>>>>>> in the comments,r AND t denote log10(rho) AND log10(t).
c---------- dlog10(p)/dx,dlog10(u)/dx,ddelad/dx,dlog10(cp)/dx

        tdvar2(n,m,21)=(tddifu(n,m,2)-tddifd(n,m,2))/(2.d0*ddx)
        tdvar2(n,m,22)=(tddifu(n,m,3)-tddifd(n,m,3))/(2.d0*ddx)
        tdvar2(n,m,23)=(tddifu(n,m,8)-tddifd(n,m,8))/(2.d0*ddx)
        tdvar2(n,m,24)=(tddifu(n,m,9)-tddifd(n,m,9))/(2.d0*ddx)
        tdvar2(n,m,25)=8888844444.d0    !space-holder variable (like var(20))
       ENDDO
      ENDDO


c======================================================================
c                   normal exit
c======================================================================

 450  PRINT 8002,(atwt(ic),abun(ic),abfrcs(ic),ic=1,nchem) ; PRINT 8003,gasmu

      RETURN

c======================================================================
c                   error exit AND error messages
c======================================================================

 500  PRINT 9800
      DO ic=1,nchem
       PRINT 9810,ic,abfrcs(ic) ; PRINT 9810,ic,abfrcd(ic)
       PRINT 9810,ic,abfrcu(ic) ; PRINT 9820
      ENDDO
      STOP
 600  PRINT 9850
      DO jt=1,nt2
       PRINT 9860,jt,tlogd(jt),tlog2(jt),tlogu(jt)
      ENDDO
      STOP
 1000 PRINT 9900,ir,idx
      STOP
c
c======================================================================
c
 98   FORMAT(1x,3i5,f13.5)
 99   FORMAT(1x,i5,(/1x,3e15.7))
1001  FORMAT(2i5,2f10.6)

8001  FORMAT('correct table construction for x-derivatives.',
     1 ' centroid composition is:'//)
8002  FORMAT('      at. weight     number ',
     1 'abundance  mass fraction',(/1x,1p3g16.7))
8003  FORMAT(/' mean molecular weight = ',f12.7//)
9006  FORMAT(' error in mhdst1. ivarr READ from table is',
     1 ' bigger than the value used in the commons.',
     2 ' ivar,ivarr= ',/1x,2i8)
9007  FORMAT(' error in mhdst1. idxr READ from table is incorrect',
     1 ' idx,idxr= ',/1x,2i8)
9008  FORMAT(' error in mhdst1. irescr READ from table is incorrect',
     1 ' iresco,irescr= ',/1x,2i8)
9010  FORMAT(' error in mhdst1. nt1 AND idx are inconsistent',
     1 ' idx,nt1,nt2 ',/1x,3i8)
9800  FORMAT(' error in table construction for x-derivatives',
     1 ' central, lower, upper table: n(element),abfrcs(n)'//)
9810  FORMAT(1x,i5,f15.9)
9820  FORMAT(/)
9850  FORMAT(' error in table construction for x-derivatives:',
     1 ' temperatures wrong: j,tlow(j),tcent(j),tupp(j)'//)
9860  FORMAT(1x,i5,3f15.9)
9900  FORMAT(' eof reached in input file. error stop. ir,idx = ',2i5)

      END SUBROUTINE mhdst1
C
C****************************************************************************
C
      SUBROUTINE rtab(ir,ntm,nrm,ivar,nt,nr,tl,tdvar)
            
c to allow derivatives with respect to X also in the ZAMS-tables situation,
c the SUBROUTINEs mhdpx1 AND mhdrx1 have been modified. they CALL an 
c enhanced quadratic-interpolation SUBROUTINE, quintd, instead of the
c usual quint. the file here contains both quint AND quintd.

	USE mod_kind

      IMPLICIT NONE
c
c....... nt is input; nr,tl,tdvar are output,n,nrr

      INTEGER :: ir,ntm,nrm,ivar,nt,nr,n,nrr,j

      REAL (kind=dp) :: tl(ntm),tdvar(ntm,nrm,ivar)

      DO n=1,nt
       READ(ir)nrr,tl(n)
       IF(n == 1)nr=nrr
       IF(nt > ntm .OR.  nr > nrm)THEN
        PRINT 9001,nt,nr,ntm,nrm ; STOP
       ENDIF
       IF(n > 1 .AND. nrr /= nr)THEN
        PRINT 9011,n,nr,nrr ; STOP
       ENDIF
       DO j=1,nr
        READ(ir)tdvar(n,j,1:ivar)
       ENDDO

	ENDDO
c
      RETURN


1001  FORMAT(i5,e16.8)
1002  FORMAT(5e16.8)
1011  FORMAT(i5,es16.8)
1012  FORMAT(5es16.8)
9001  FORMAT(' error in rtab. too small PARAMETERs: ',
     1 /' nt,nr,ntm,nrm = ',4i5)
9011  FORMAT(' error in rtab. wrong density construction: ',
     1 /' n,nr,nrr = ',3i5)


      END SUBROUTINE rtab
      
c****************************************************************      
      
      SUBROUTINE rabu(ir,nchem0,nchem,atwt,abun,abfrcs,gasmu)

	USE mod_kind
      
      IMPLICIT NONE
      
      INTEGER :: ir,nchem0,nchem,ic

      REAL (kind=dp) ::   atwt(nchem0),abun(nchem0),abfrcs(nchem0),gasmu

c....... nchem,atwt,abun,abfrcs are output

      READ(ir) nchem,(atwt(ic),abun(ic),abfrcs(ic),ic=1,nchem),gasmu
      IF(nchem0<nchem)THEN
       PRINT 9009,nchem0,nchem ; STOP
      ENDIF

      RETURN


 99   FORMAT(1x,i5,(/1x,3es15.7))
9009  FORMAT(' error in rabu. nchem READ from table is',
     1 ' bigger than the value used in the commons.',
     2 ' nchem0,nchem= ',/1x,2i8)


      END SUBROUTINE rabu
      
c********************************************************

      SUBROUTINE quintd(x,x0,h,y0,y1,y2,y,dydx)

	USE mod_kind
      
      IMPLICIT NONE
      
      REAL(kind=dp) :: x0,h,y0,y1,y2,y,dydx,d1,d2,rh,t,x,dydt
      
c...... quadratic interpolation for equidistant points
c...... y0=y(x0),y1=y(x1),y2=y(x2); h=x1-x0=x2-x1;
c...... computes y=y(x) AND dydx. 

      d1=y1-y0 ; d2=y2-2.d0*y1+y0 ; rh=1.d0/h ; t=(x-x0)*rh
      y=y0+t*d1+0.5d0*t*(t-1.d0)*d2
      dydt=d1+(t-0.5d0)*d2 ;  dydx=dydt*rh

      RETURN
      
      END SUBROUTINE quintd
      
c***************************************************************      
      
      SUBROUTINE lir(z,zi,y,yi,ii,id,nt,l,inter)
      
	USE mod_kind
	
      IMPLICIT NONE
      
      INTEGER :: ii,id,nt,l,inter,n=-1,il,ir,ir1,ird,iir,j,i,m,k
          
      REAL (kind=dp) :: z,zi(1),y(1),yi(1),a(4),diff,y1,y2,y3,y4,

     1 z4,z1,z2,z3,z12,z34,yy   

c                interpolation/extrapolation routine

c     for a such that z=zi(a),  sets y(i)=yi(i,a), i=1,ii

c     zi(n),yi(i,n) must be supplied for n=1,nt AND i=1,ii
c     id is first DIMENSION of yi

c     inter is set to 1 for interpolation AND 0 for extrapolation

c     if l<=1, scan to find the zi(n) which immediately bound z
c               starts at n=1
c     if l>1, scan starts from value of n from previous CALL of lir

c     lir use cubic interpolation/extrapolation unless nt<4
c     lir1 use linear interpolation/extrapolation

c     note
c     ****
c     most of the computation is performed in single precision

      il=0 ; ir=1

c     check nt AND reset il if necessary


      IF(nt < 2)GOTO 101
      IF(nt < 4) il=1

c     addressing constants


      inter=1 ; ir1=ir-1 ; ird=ir*id ; iir=(ii-1)*ir+1
      j=(nt-1)*ir+1 ; diff=zi(j)-zi(1)

c     set index for start of search


      n=(n-2)*ir+1 ; IF(l <= 1 .OR. n < 1)n=1

c     determine position of z within zi

    2 IF(n > j)GOTO 8
      IF(diff)4,102,3
    3 IF(zi(n)-z)5,6,9
    4 IF(zi(n)-z)9,6,5
    5 n=n+ir
      GOTO 2

c     set y when z lies on a mesh point


    6 j=(n-1)*id
      DO i=1,iir
       y(i)=yi(i+j) ; IF(y(i) == 0.d0)y(i+ir1)=0.d0
      ENDDO
      GOTO 30

c     control when z does not lie on a mesh point


    8 inter=0
    9 IF(n <= 1)inter=0
      IF(il == 1)GOTO 20

c     cubic interpolation/extrapolation
c     pivotal point (m) AND point (k) closest to z


   10 m=n ; k=3
      IF(n > 1+ir)GOTO 11
      m=1+ir+ir ; k=n
   11 IF(n < j)GOTO 12
      m=j-ir ; k=4

c     weighting factors


   12 y1=zi(m-ir*2) ; y2=zi(m-ir) ; y3=zi(m) ; y4=zi(m+ir)
      z1=z-y1 ; z2=z-y2 ; z3=z-y3 ; z4=z-y4

   13 z12=z1*z2 ; z34=z3*z4

   14 a(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      a(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      a(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      a(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))

c     correct a(k)

   15 diff=a(1)+a(2)+a(3)+a(4) ; a(k)=(1.d0+a(k))-diff

c     compute y

   16 m=(m-1)/ir-3 ; m=m*ird
      DO i=1,iir
        k=i+m ; yy=0.d0
        DO j=1,4
           k=k+ird ; diff=yi(k) ; yy=yy+a(j)*diff
        ENDDO
        y(i)=yy ; IF(y(i) == 0.d0)y(i+ir1)=0.d0
	ENDDO
      GOTO 30

c     linear interpolation/extrapolation


   20 IF(n == 1)n=1+ir
      IF(n > j)n=j
      z1=zi(n) ; y1=(z1-z)/(z1-zi(n-ir)) ; y2=1.0-y1 ; j=(n-1)*id ; m=j-ird
      DO i=1,iir,ir
       y(i)=y1*yi(i+m)+y2*yi(i+j)
       IF(y(i) == 0.d0)y(i+ir1)=0.d0

      ENDDO

c     reset n

   30 n=(n+ir-1)/ir

      RETURN

c     diagnostics

  101 WRITE(*,1001)nt ; RETURN
 1001 FORMAT(/,10('*'),5x,'there are fewer than two DATA points in',
     1 ' lir     nt =',i4,5x,10('*')/)  
  102 WRITE(*,1002)zi(1),nt,zi(j) ; RETURN
 1002 FORMAT(/,10('*'),5x,'extreme values of independent variable',
     1 ' equal in lir',5x,10('*')/16x,'zi(1) =',1pe13.5,',   ',
     2 'zi(',i4,') =',1pe13.5/)
     
      END SUBROUTINE lir
